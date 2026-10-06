// Exercises computeTracklets's kernel mechanics (allocation, concurrent
// reuse, allocation-failure recovery) against a small hand-built
// timeframe. The candidate-search algorithm itself (window projection, ROF
// navigation, vertex policy) is validated far more thoroughly against the
// real production pipeline in testComputeLayerTrackletsOrchestration.cxx
// (CPU vs GPU, multiple geometries, edge cases); this file stays
// deliberately minimal and uses a single trivial bin (nRowBins = nColBins =
// 1) so index-table arithmetic cannot introduce ambiguity, keeping the
// expected candidate set obvious.

#include <boost/test/unit_test.hpp>
#include <array>
#include <future>
#include <mutex>
#include <new>
#include <utility>
#include <vector>
#include <stdexcept>
#include "ITSMFTTrackingGPU/TimeFrameGPU.h"

using namespace o2::itsmft::tracking;

namespace
{
// One source ROF (3 clusters), one target ROF (2 clusters), single global
// bin. source[0] and source[2] are geometrically identical and should each
// find target[0]; source[1] is marked used and must be skipped entirely;
// target[1] is marked used and must never be matched. Used flags are
// indexed by cluster id, which is the cluster's position here.
struct Fixture {
  int layer{0};
  int toLayer{1};
  std::vector<GlobalMeasurement> sources, targets;
  std::vector<uint8_t> sourceUsed{0, 1, 0};
  std::vector<uint8_t> targetUsed{0, 1};
  std::vector<int> sourceROFBoundaries{0, 3};
  std::vector<int> targetROFBoundaries{0, 2};
  std::vector<uint8_t> sourceROFEnabled{1};
  std::vector<uint8_t> targetROFEnabled{1};
  std::vector<RuntimeROFTableEntry> overlapRanges{RuntimeROFTableEntry{0, 1}};
  // tableSize = nRowBins * nColBins + 1 = 1*1 + 1 = 2; both target clusters
  // are in the single bin, so the device-built table is {0, 2} (used-flag
  // filtering happens in evaluateTracklet).
  std::vector<int> targetClusterBins{0, 0};
  std::vector<uint8_t> targetTableBuilt{1};
  o2::itsmft::IndexTableUtilsCore indexTableUtils;
  Fixture()
  {
    const std::array<float, 1> colHalfExtent{1000.f};
    indexTableUtils.setIndexTableParams(o2::itsmft::IndexTableCoordType::PhiR, 1, 1, 0.f,
                                        o2::constants::math::TwoPI,
                                        gsl::span<const float>{colHalfExtent});
    for (uint32_t i = 0; i < 3; ++i) {
      sources.push_back(cluster(3.f, 0.3f, i));
    }
    for (uint32_t i = 0; i < 2; ++i) {
      targets.push_back(cluster(4.f, 0.4f, i));
    }
  }
  static GlobalMeasurement cluster(float r, float z, uint32_t id)
  {
    GlobalMeasurement measurement{};
    measurement.x = r;
    measurement.y = 0.f;
    measurement.z = z;
    measurement.radius = r;
    measurement.phi = 0.f;
    measurement.clusterId = id;
    return measurement;
  }

  // Uploads both layers as TrackerTraitsGPU does for a real timeframe.
  void load(gpu::TimeFrameGPU& frame) const
  {
    const size_t nLayers = std::max(layer, toLayer) + 1;
    frame.initialise(nLayers);
    std::vector<std::span<const uint8_t>> used(nLayers), enabled(nLayers);
    for (size_t index = 0; index < nLayers; ++index) {
      gpu::LayerHostData host;
      if (static_cast<int>(index) == layer) {
        host.clusters = sources;
        host.rofClusters = sourceROFBoundaries;
        used[index] = sourceUsed;
        enabled[index] = sourceROFEnabled;
      } else if (static_cast<int>(index) == toLayer) {
        host.clusters = targets;
        host.rofClusters = targetROFBoundaries;
        host.clusterBins = targetClusterBins;
        host.tableBuilt = targetTableBuilt;
        host.tableSize = 2;
        used[index] = targetUsed;
        enabled[index] = targetROFEnabled;
      }
      frame.loadLayer(static_cast<int>(index), host);
    }
    frame.loadUsedClusters(used);
    frame.loadROFEnabled(enabled);
    frame.loadVertices({});
    frame.finishLoading();
  }

  TrackletSearch makeInput() const
  {
    TrackletSearch input{};
    input.overlapRanges = overlapRanges.data();
    input.nOverlapRanges = overlapRanges.size();
    input.edgeCache = TrackletProjectionCache{layer, toLayer, 3.f, 4.f, -1000.f, 1000.f, 0.f, 0.f, 1.e-4f, 0.f, 1.f};
    input.indexTableUtils = indexTableUtils;
    input.nSigmaCut = 100.f;
    input.beamPositionVariance = 1.e-4f;
    input.useDiamond = true;
    input.diamondBase = Vertex{};
    input.fromLayerTiming = o2::itsmft::tracking::ROFTimingLayer{1, 40, 0, 0, 1000};
    input.toLayerTiming = o2::itsmft::tracking::ROFTimingLayer{1, 40, 0, 0, 1000};
    input.selectUPC = false;
    input.iVertex = -1;
    input.kind = SurfaceKind::Cylinder;
    return input;
  }
};

// A loaded timeframe.
struct Setup {
  int device{gpu::currentDevice()};
  gpu::TimeFrameGPU frame{device};
  explicit Setup(const Fixture& fixture, size_t nEdges = 1)
  {
    fixture.load(frame);
    frame.initialiseTracklets(nEdges);
  }
  // Runs one edge with the given buffer capacity; returns the number found.
  int search(const TrackletSearch& input, int capacity, int edge = 0)
  {
    frame.tracklets(edge).reserve(capacity);
    return gpu::computeTracklets(frame, input, edge);
  }
  // As the traits do: grows the buffer until everything fits, then downloads.
  std::vector<Tracklet> run(const TrackletSearch& input, int edge = 0)
  {
    int capacity = 0;
    for (int found = search(input, capacity, edge); found > capacity; found = search(input, capacity, edge)) {
      capacity = found;
    }
    std::vector<Tracklet> tracklets;
    frame.download(frame.tracklets(edge).items, tracklets);
    frame.download(frame.tracklets(edge).lookup, lookup);
    return tracklets;
  }
  std::vector<int> lookup;
};

// A device memory pool behaving as GPUReconstruction::AllocateDirectMemory:
// plain allocations grow from its start until clear(), stack allocations
// from its end until the pop of the tag pushed before them.
class PoolAllocator final : public ExternalAllocator
{
 public:
  explicit PoolAllocator(size_t bytes) : mEnd{bytes} { mPool.reserve(bytes); }
  void* allocateDevice(size_t bytes, bool stack) final
  {
    std::lock_guard lock{mMutex};
    bytes = (bytes + 63) / 64 * 64;
    if (bytes > mEnd - mBegin) {
      throw std::bad_alloc{};
    }
    ++(stack ? nStack : nPlain);
    if (stack) {
      mEnd -= bytes;
      return static_cast<char*>(mPool.data) + mEnd;
    }
    mBegin += bytes;
    return static_cast<char*>(mPool.data) + mBegin - bytes;
  }
  void pushTagOnStack(uint64_t tag) final { mStack.emplace_back(tag, mEnd); }
  void popTagOffStack(uint64_t tag) final
  {
    BOOST_REQUIRE(!mStack.empty());
    BOOST_CHECK_EQUAL(mStack.back().first, tag);
    mEnd = mStack.back().second;
    mStack.pop_back();
  }
  void clear()
  {
    BOOST_CHECK(mStack.empty());
    mBegin = 0;
    mEnd = mPool.capacity;
  }
  bool contains(const void* p) const { return p >= mPool.data && p < static_cast<const char*>(mPool.data) + mPool.capacity; }
  size_t plainBytes() const { return mBegin; }
  size_t stackBytes() const { return mPool.capacity - mEnd; }
  int nPlain{0}, nStack{0};

 private:
  gpu::DeviceBuffer mPool;
  size_t mBegin{0}, mEnd;
  std::vector<std::pair<uint64_t, size_t>> mStack;
  std::mutex mMutex;
};

// Sources 0 and 2 both reach target 0, in (first, second) cluster order.
void checkExpectedPair(const std::vector<Tracklet>& results)
{
  BOOST_REQUIRE_EQUAL(results.size(), 2u);
  const float expectedTanLambda = (0.3f - 0.4f) / (3.f - 4.f);
  const float expectedPhi = o2::gpu::CAMath::ATan2(0.f, -1.f);
  for (size_t i = 0; i < results.size(); ++i) {
    BOOST_CHECK_EQUAL(results[i].firstClusterIndex, 2 * static_cast<int>(i));
    BOOST_CHECK_EQUAL(results[i].secondClusterIndex, 0);
    BOOST_CHECK_SMALL(results[i].tanLambda - expectedTanLambda, 2.e-6f);
    BOOST_CHECK_SMALL(results[i].phi - expectedPhi, 2.e-6f);
  }
}
} // namespace

BOOST_AUTO_TEST_CASE(DeviceTrackletSearchBasicAndEdgeCases)
{
  const Fixture fixture;
  Setup setup{fixture};
  checkExpectedPair(setup.run(fixture.makeInput()));
  // Lookup table: first tracklet of each source cluster, plus the total.
  const std::vector<int> expectedLookup{0, 1, 1, 2};
  BOOST_CHECK_EQUAL_COLLECTIONS(setup.lookup.begin(), setup.lookup.end(), expectedLookup.begin(), expectedLookup.end());

  // Too small a buffer: everything is counted, nothing else is done.
  BOOST_CHECK_EQUAL(setup.search(fixture.makeInput(), 1), 2);
  BOOST_CHECK_EQUAL(setup.frame.tracklets(0).items.size, 0u);
  checkExpectedPair(setup.run(fixture.makeInput()));

  // Without a loaded timeframe there is nothing to search.
  gpu::TimeFrameGPU unloaded{setup.device};
  unloaded.initialiseTracklets(1);
  BOOST_CHECK_THROW(gpu::computeTracklets(unloaded, fixture.makeInput(), 0), std::logic_error);

  // Invalid navigation-table dimensions must throw, not read out of bounds.
  auto badOverlap = fixture.makeInput();
  badOverlap.nOverlapRanges = 0;
  BOOST_CHECK_THROW(setup.search(badOverlap, 2), std::out_of_range);
  auto badLayer = fixture.makeInput();
  badLayer.edgeCache.toLayer = 5;
  BOOST_CHECK_THROW(setup.search(badLayer, 2), std::out_of_range);
  gpu::TimeFrameGPU invalid{setup.device};
  invalid.initialise(1);
  gpu::LayerHostData badBins{fixture.targets, fixture.targetROFBoundaries, std::span<const int>{fixture.targetClusterBins}.first(1),
                             fixture.targetTableBuilt, 2};
  BOOST_CHECK_THROW(invalid.loadLayer(0, badBins), std::invalid_argument);
  checkExpectedPair(setup.run(fixture.makeInput()));
}

// A target ROF whose index table the host never built has an all-zero table,
// so it yields no candidates even though it is enabled and has clusters.
BOOST_AUTO_TEST_CASE(UnbuiltTargetIndexTableYieldsNoCandidates)
{
  Fixture fixture;
  fixture.targetTableBuilt = {0};
  Setup setup{fixture};
  BOOST_CHECK(setup.run(fixture.makeInput()).empty());
}

BOOST_AUTO_TEST_CASE(PersistentWorkspaceConcurrentTrackletCalls)
{
  const Fixture fixture;
  Setup setup{fixture, 2};
  // Edges run concurrently on independent stream/buffer slots (see
  // TrackingKernels.cu), each into its own buffers. Warm up first so every
  // slot pays its one-time allocation cost before the "before" snapshot.
  const auto search = [&](int edge) {
    std::vector<Tracklet> tracklets;
    if (setup.search(fixture.makeInput(), 8, edge) == 2) {
      setup.frame.download(setup.frame.tracklets(edge).items, tracklets);
    }
    return tracklets;
  };
  {
    auto warmFirst = std::async(std::launch::async, search, 0);
    auto warmSecond = std::async(std::launch::async, search, 1);
    checkExpectedPair(warmFirst.get());
    checkExpectedPair(warmSecond.get());
  }
  const auto before = setup.frame.statistics();
  auto first = std::async(std::launch::async, search, 0);
  auto second = std::async(std::launch::async, search, 1);
  checkExpectedPair(first.get());
  checkExpectedPair(second.get());
  BOOST_CHECK_EQUAL(setup.frame.statistics().allocations, before.allocations);
}

// A new timeframe must be loaded before use and its data replaces the
// previous one's.
BOOST_AUTO_TEST_CASE(NewTimeframeReplacesResidentData)
{
  Fixture first;
  Setup setup{first};
  checkExpectedPair(setup.run(first.makeInput()));

  setup.frame.initialise(0); // a new timeframe starts
  setup.frame.initialiseTracklets(1);
  BOOST_CHECK_THROW(setup.run(first.makeInput()), std::logic_error);

  // Target z moved from 0.4 to 0.42 changes the expected tanLambda from
  // 0.1 to 0.12 = (0.42 - 0.3) / hypot(4 - 3, 0 - 0).
  Fixture second;
  for (auto& target : second.targets) {
    target.z = 0.42f;
  }
  second.load(setup.frame);
  setup.frame.initialiseTracklets(1);
  const auto output = setup.run(second.makeInput());
  BOOST_REQUIRE_EQUAL(output.size(), 2u);
  for (const auto& r : output) {
    BOOST_CHECK_SMALL(r.tanLambda - 0.12f, 2.e-6f);
  }
}

// Used flags reloaded after accepted tracks (as TrackerTraitsGPU does)
// drop newly used clusters from the next search.
BOOST_AUTO_TEST_CASE(ReloadedUsedFlagsApply)
{
  Fixture fixture;
  Setup setup{fixture};
  checkExpectedPair(setup.run(fixture.makeInput()));
  fixture.sourceUsed[2] = 1;
  const std::array<std::span<const uint8_t>, 2> used{fixture.sourceUsed, fixture.targetUsed};
  setup.frame.loadUsedClusters(used);
  const auto output = setup.run(fixture.makeInput());
  BOOST_REQUIRE_EQUAL(output.size(), 1u);
  BOOST_CHECK_EQUAL(output[0].firstClusterIndex, 0);
}

// With a framework allocator the frame takes all device memory from it:
// the timeframe's data as plain allocations, stage results and scratch on
// the memory stack of the pass, and everything again after each release.
BOOST_AUTO_TEST_CASE(FrameworkAllocatorScopesDeviceMemory)
{
  const Fixture fixture;
  Setup setup{fixture}; // its runtime-allocated buffers are replaced
  checkExpectedPair(setup.run(fixture.makeInput()));
  PoolAllocator pool{size_t(1) << 24};
  auto& frame = setup.frame;
  frame.setFrameworkAllocator(&pool);
  BOOST_CHECK(frame.hasFrameworkAllocator());
  size_t timeframeBytes = 0;
  for (int timeframe = 0; timeframe < 3; ++timeframe) {
    pool.clear();
    fixture.load(frame);
    BOOST_CHECK_GT(pool.plainBytes(), 0u);
    BOOST_CHECK_EQUAL(pool.stackBytes(), 0u);
    BOOST_CHECK(pool.contains(frame.view().layers[fixture.layer].clusters));
    BOOST_CHECK(pool.contains(frame.deviceView()));
    if (timeframe) {
      BOOST_CHECK_EQUAL(pool.plainBytes(), timeframeBytes);
    }
    timeframeBytes = pool.plainBytes();
    for (int iteration = 0; iteration < 2; ++iteration) {
      const int nPlain = pool.nPlain;
      frame.pushMemoryStack(iteration);
      frame.initialiseTracklets(1);
      checkExpectedPair(setup.run(fixture.makeInput()));
      BOOST_CHECK(pool.contains(frame.tracklets(0).items.data()));
      BOOST_CHECK(pool.contains(frame.tracklets(0).lookup.data()));
      BOOST_CHECK_GT(pool.stackBytes(), 0u);
      frame.popMemoryStack(iteration);
      BOOST_CHECK_EQUAL(pool.stackBytes(), 0u);
      BOOST_CHECK_EQUAL(pool.nPlain, nPlain);
    }
  }
  // Back on the device runtime nothing refers to the pool any more.
  frame.setFrameworkAllocator(nullptr);
  const int nAllocations = pool.nPlain + pool.nStack;
  fixture.load(frame);
  frame.initialiseTracklets(1);
  checkExpectedPair(setup.run(fixture.makeInput()));
  BOOST_CHECK(!pool.contains(frame.tracklets(0).items.data()));
  BOOST_CHECK_EQUAL(pool.nPlain + pool.nStack, nAllocations);
}

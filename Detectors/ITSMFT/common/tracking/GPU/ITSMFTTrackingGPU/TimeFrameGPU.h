// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file TimeFrameGPU.h
/// \brief Device-resident TimeFrame and the CA stage handlers, as the legacy ITS TimeFrameGPU and TrackingKernels
///

#ifndef ALICEO2_ITSMFT_TRACKING_TIMEFRAMEGPU_H_
#define ALICEO2_ITSMFT_TRACKING_TIMEFRAMEGPU_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <span>
#include <stdexcept>
#include <vector>

#include "ITSMFTTracking/ExternalAllocator.h"
#include "ITSMFTTracking/GenericTrack.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingKernels.h"

// The CUDA and HIP libraries are built from the same sources and may be
// loaded in one process, so each has its own namespace (headers are not
// hipified: the HIP library defines ITSMFT_TRACKING_HIP for its users).
#if defined(__HIPCC__) || defined(ITSMFT_TRACKING_HIP)
#define ITSMFT_TRACKING_GPU_BACKEND hip
#else
#define ITSMFT_TRACKING_GPU_BACKEND cuda
#endif

namespace o2::itsmft::tracking::gpu
{
inline namespace ITSMFT_TRACKING_GPU_BACKEND
{

int currentDevice();

// Where a TimeFrameGPU's buffers come from. By default the device runtime
// (cudaMalloc): buffers only grow and are kept across timeframes. With a
// framework allocator, as in the legacy tracker, they come from the
// framework's pool, which frees nothing individually: the timeframe's data
// lives until the framework clears it, stage results and scratch on the
// memory stack of one traversal pass. Buffers notice lazily that their
// memory was released (the scope's epoch moved on) and start again from
// nothing.
class DeviceMemory
{
 public:
  enum class Scope : uint8_t { TimeFrame,
                               Pass };
  // May be changed between timeframes; null selects the device runtime.
  void setFrameworkAllocator(ExternalAllocator* allocator)
  {
    mExternal = allocator;
    release(Scope::TimeFrame);
    release(Scope::Pass);
  }
  bool hasFrameworkAllocator() const noexcept { return mExternal != nullptr; }
  ExternalAllocator* getFrameworkAllocator() const noexcept { return mExternal; }
  // The scope's buffers are gone.
  void release(Scope scope) noexcept { ++mEpochs[static_cast<int>(scope)]; }
  uint32_t epoch(Scope scope) const noexcept { return mEpochs[static_cast<int>(scope)]; }

 private:
  ExternalAllocator* mExternal{};
  std::array<uint32_t, 2> mEpochs{};
};

// Device memory that only grows. Unbound buffers use the device runtime.
struct DeviceBuffer {
  void* data{};
  size_t capacity{}; // bytes
  DeviceBuffer() = default;
  DeviceBuffer(const DeviceMemory* memory, DeviceMemory::Scope scope) : mMemory{memory}, mScope{scope} {}
  DeviceBuffer(const DeviceBuffer&) = delete;
  DeviceBuffer& operator=(const DeviceBuffer&) = delete;
  DeviceBuffer(DeviceBuffer&& other) noexcept;
  DeviceBuffer& operator=(DeviceBuffer&& other) noexcept;
  ~DeviceBuffer();
  // Returns true if the buffer had to grow (its content is then lost).
  bool reserve(size_t bytes);
  // The capacity that is still allocated.
  size_t usable() const noexcept { return released() ? 0 : capacity; }
  // An empty buffer allocating as this one.
  DeviceBuffer sibling() const { return {mMemory, mScope}; }

 private:
  bool released() const noexcept { return mMemory && mEpoch != mMemory->epoch(mScope); }
  void free() noexcept;
  const DeviceMemory* mMemory{};
  DeviceMemory::Scope mScope{DeviceMemory::Scope::TimeFrame};
  uint32_t mEpoch{};
  bool mOwned{false}; // allocated by the device runtime
};

template <typename T>
struct DeviceVector {
  DeviceBuffer buffer;
  size_t size{};
  DeviceVector() = default;
  DeviceVector(const DeviceMemory* memory, DeviceMemory::Scope scope) : buffer{memory, scope} {}
  T* data() const noexcept { return static_cast<T*>(buffer.data); }
  T* reserve(size_t n)
  {
    buffer.reserve(n * sizeof(T));
    return data();
  }
};

// One edge's or path's results: the items, their slab capacity (see
// runOnSlab) and their lookup table, the first item of every source element
// plus the total.
template <typename T>
struct LinkedResults {
  DeviceVector<T> items;
  size_t capacity{};
  DeviceVector<int> lookup;
  LinkedResults() = default;
  LinkedResults(const DeviceMemory* memory, DeviceMemory::Scope scope) : items{memory, scope}, lookup{memory, scope} {}
  void reserve(size_t n)
  {
    items.reserve(n);
    capacity = n;
  }
};

// One layer as the tracker keeps it on the host. Index tables are built on
// the device from each cluster's bin and which ROFs have a table (see
// TimeFrame::getClusterBins).
struct LayerHostData {
  std::span<const GlobalMeasurement> clusters; // sorted by ROF, then index-table bin
  std::span<const int> rofClusters;            // first cluster of every ROF, plus the total
  std::span<const int> clusterBins;
  std::span<const uint8_t> tableBuilt;
  size_t tableSize{}; // entries per ROF: bins + 1
  std::span<const SurfaceMeasurement> surfaceMeasurements; // by cluster id
};

struct WorkspaceStatistics {
  size_t allocations{}, reservedBytes{}, uploadedBytes{};
};
struct Workspace; // the handlers' streams and scratch

// Everything the CA stages keep on the device, as the legacy TimeFrameGPU:
// the timeframe, uploaded once like its load*Device calls (the used flags,
// which change with every accepted track, and the per-iteration ROF masks
// are reloaded explicitly), and every stage's results. view() mirrors on the
// host the FrameView the kernels read from deviceView(). Loading is not
// thread-safe; handlers on distinct edges, paths or targets may run
// concurrently.
//
// With a framework allocator (setFrameworkAllocator, as the legacy
// TimeFrameGPU) all device memory comes from the framework: initialise()
// assumes the framework cleared the previous timeframe's, and stage results
// and scratch must be used between pushMemoryStack and popMemoryStack.
class TimeFrameGPU
{
 public:
  explicit TimeFrameGPU(int device = currentDevice());
  ~TimeFrameGPU();
  TimeFrameGPU(const TimeFrameGPU&) = delete;
  TimeFrameGPU& operator=(const TimeFrameGPU&) = delete;
  WorkspaceStatistics statistics() const;
  Workspace& workspace() const noexcept { return *mWorkspace; }

  void setFrameworkAllocator(ExternalAllocator* allocator) { mMemory.setFrameworkAllocator(allocator); }
  bool hasFrameworkAllocator() const noexcept { return mMemory.hasFrameworkAllocator(); }
  // Marks and releases the device memory of one traversal pass of an
  // iteration: no-ops without a framework allocator, whose buffers are kept.
  void pushMemoryStack(int iteration);
  void popMemoryStack(int iteration);

  // A new timeframe: every layer must be loaded again before use.
  void initialise(size_t nLayers);
  bool isLoaded() const noexcept { return mLoaded; }
  void finishLoading();
  // Uploads a layer on its own stream (see waitLayer).
  void loadLayer(int layer, const LayerHostData& host);
  void loadUsedClusters(std::span<const std::span<const uint8_t>> used);
  void loadROFEnabled(std::span<const std::span<const uint8_t>> enabled);
  void loadVertices(std::span<const Vertex> vertices);
  // As the tracker loads the frame: every measurement surface, with index
  // tables for the first nTableLayers.
  void loadTimeFrame(TimeFrame& frame, size_t nTableLayers);
  void loadUsedClusters(TimeFrame& frame);
  void loadROFEnabled(const TimeFrame& frame);
  // Makes a cudaStream_t/hipStream_t wait for a layer's upload.
  void waitLayer(void* stream, int layer) const;
  const FrameView& view() const noexcept { return mView; }
  const FrameView* deviceView() const noexcept { return mDeviceView.data(); }

  // Tracklets of every edge (legacy createTrackletsBuffers), sorted by
  // (first, second) cluster without duplicates, indexed by source cluster.
  void initialiseTracklets(size_t nEdges) { resetResults(mTracklets, nEdges); }
  LinkedResults<Tracklet>& tracklets(int edge) { return mTracklets.at(edge); }
  // Cells of every path (legacy createCellsBuffers), sorted by (first,
  // second) tracklet, indexed by first-edge tracklet. A new cell stage makes
  // the road graph stale.
  void initialiseCells(size_t nPaths)
  {
    resetResults(mCells, nPaths);
    ++mCellGeneration;
  }
  LinkedResults<Triplet>& cells(int path) { return mCells.at(path); }
  uint64_t getCellGeneration() const noexcept { return mCellGeneration; }
  // Neighbours of every target path (legacy createNeighboursDevice), sorted
  // by (target cell, source path, source cell), indexed by target cell.
  void initialiseNeighbours(size_t nPaths) { resetResults(mNeighbours, nPaths); }
  LinkedResults<CellNeighbour>& neighbours(int path) { return mNeighbours.at(path); }
  size_t getNPaths() const noexcept { return mCells.size(); }

  // The road stage (legacy processNeighboursHandler): the graph of the
  // current cells, and the results of the last extension.
  DeviceVector<RoadGraphPath> roadGraph{&mMemory, DeviceMemory::Scope::Pass};
  uint64_t roadGraphGeneration{0};
  std::array<DeviceVector<RoadResult>, 2> roads{{{&mMemory, DeviceMemory::Scope::Pass}, {&mMemory, DeviceMemory::Scope::Pass}}};
  int lastRoads{0};
  std::mutex roadMutex; // the road and refit handlers are serialised
  // The track seeds of one road start level, appended start path by start
  // path; their refit results, and the accepted tracks sorted for acceptance.
  DeviceVector<TrackSeed> trackSeeds{&mMemory, DeviceMemory::Scope::Pass};
  DeviceVector<RefitResult> refitResults{&mMemory, DeviceMemory::Scope::Pass};
  DeviceVector<TrackingCandidate> tracks{&mMemory, DeviceMemory::Scope::Pass};

  template <typename T, typename Container>
  void download(const DeviceVector<T>& from, Container& to) const
  {
    to.resize(from.size);
    copy(to.data(), from.data(), from.size * sizeof(T), false);
  }
  template <typename T>
  void upload(DeviceVector<T>& to, std::span<const T> from)
  {
    to.reserve(from.size());
    to.size = from.size();
    copy(to.data(), from.data(), from.size_bytes(), true);
  }

 private:
  struct Layer {
    explicit Layer(const DeviceMemory* memory)
      : clusters{memory, DeviceMemory::Scope::TimeFrame}, rofClusters{clusters.sibling()}, bins{clusters.sibling()}, built{clusters.sibling()}, indexTables{clusters.sibling()}, measurements{clusters.sibling()} {}
    DeviceBuffer clusters, rofClusters, bins, built, indexTables, measurements;
    void* stream{};
    void* ready{};
  };
  template <typename T>
  void resetResults(std::vector<LinkedResults<T>>& results, size_t n)
  {
    while (results.size() < n) {
      results.emplace_back(&mMemory, DeviceMemory::Scope::Pass);
    }
    for (auto& result : results) {
      result.items.size = result.lookup.size = result.capacity = 0;
    }
  }
  void copy(void* to, const void* from, size_t bytes, bool toDevice) const;
  void uploadAsync(void* stream, void* to, const void* from, size_t bytes);
  // Concatenates per-layer flags into one buffer, pointed to by the view.
  void loadLayers(DeviceBuffer& buffer, std::span<const std::span<const uint8_t>> layers, const uint8_t* LayerView::*flags);
  void publishView();

  int mDevice;
  bool mLoaded{false};
  DeviceMemory mMemory;
  std::vector<Layer> mLayers;
  DeviceBuffer mUsed{&mMemory, DeviceMemory::Scope::TimeFrame}, mEnabled{mUsed.sibling()}, mVertices{mUsed.sibling()};
  FrameView mView;
  DeviceVector<FrameView> mDeviceView{&mMemory, DeviceMemory::Scope::TimeFrame};
  std::vector<LinkedResults<Tracklet>> mTracklets;
  std::vector<LinkedResults<Triplet>> mCells;
  std::vector<LinkedResults<CellNeighbour>> mNeighbours;
  uint64_t mCellGeneration{0};
  void* stream() const; // of loads and copies, created on first use
  mutable std::once_flag mStreamCreated;
  mutable void* mStream{};
  std::unique_ptr<Workspace> mWorkspace; // last: its streams drain first
};

// The handlers, as the legacy TrackingKernels handlers. Those that find a
// stage's results fill the edge's or path's buffer, reserved beforehand,
// then sort them and build their lookup table; they return the number found,
// and if that exceeds the capacity nothing else was done and the caller
// retries with a larger buffer (runOnSlab). They are thread-safe across
// edges, paths, or targets whose sources are final.

// Search pointers are host pointers: the ROF ranges are uploaded.
int computeTracklets(TimeFrameGPU& frame, const TrackletSearch& search, int edge);
// Search tracklet pointers are ignored: they are the edges' resident ones.
int computeCells(TimeFrameGPU& frame, const CellSearch& search, int firstEdge, int secondEdge, int path);
// Every source path's neighbours among the target path's cells, raising the
// target cells' levels; search cell pointers are ignored.
int computeNeighbours(TimeFrameGPU& frame, int targetPath, std::span<const NeighbourSearch> sources);

// One road extension of startPath's cells, or of the previous extension's
// results, built on the device from the graph bound by bindRoadGraph (call
// after the neighbour stage) as TrackerTraits::buildRoadJobs does on the
// host. The results stay on the device (frame.roads[frame.lastRoads]).
struct RoadExtension {
  int startPath{-1};
  bool reusePreviousSeeds{false};
  int currentLevel{};
  int nLayers{};
  float bz{}, maxChi2{};
  std::span<const SurfaceDescriptor> surfaces; // indexed by layer
};
void bindRoadGraph(TimeFrameGPU& frame);
size_t extendRoads(TimeFrameGPU& frame, const RoadExtension& extension);
// Appends the last extension's results that pass the selector, in order, to
// the track seeds; returns their total.
size_t selectRoadSeeds(TimeFrameGPU& frame, const RoadSeedSelector& selector);
// Refits every track seed, then keeps the accepted ones sorted as the host
// sorts tracks (longer, then lower chi2, ties in seed order), as the legacy
// computeTrackSeedHandler; returns their number.
size_t refitTrackSeeds(TimeFrameGPU& frame, std::span<const SurfaceDescriptor> surfaces, const RefitParameters& parameters,
                       std::span<const float> minPt);

// Kernel harness for tests: evaluates host-built road or refit jobs with the
// device kernels and downloads the results. Road jobs reusing the previous
// extension's seeds read them by job.source, of nSources.
struct RoadInput {
  std::span<const RoadJob> jobs;
  std::span<const RoadTarget> targets;
  float bz, maxChi2;
  bool reusePreviousSeeds{false};
  size_t nSources{};
};
void computeRoads(TimeFrameGPU& frame, const RoadInput& input, bounded_vector<RoadResult>& output);
struct RefitInput {
  std::span<const RefitJob> jobs;
  std::span<const SurfaceDescriptor> surfaces;
  RefitParameters parameters;
};
void refitTracks(TimeFrameGPU& frame, const RefitInput& input, bounded_vector<RefitResult>& output);

} // namespace ITSMFT_TRACKING_GPU_BACKEND
} // namespace o2::itsmft::tracking::gpu

#endif

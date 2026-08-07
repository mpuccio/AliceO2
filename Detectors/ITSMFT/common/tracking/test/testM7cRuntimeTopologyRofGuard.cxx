// Gate 4 M7c sparse-topology and runtime-ROF ownership coverage.
//
// The common production scan below has no frozen-legacy exclusions. Frozen
// ITStracking sources are outside this common-tracking scan and remain
// unchanged.

#define BOOST_TEST_MODULE ITSMFT M7c runtime topology and ROF guard
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/ROFViews.h"
#include "ITSMFTTracking/SurfaceGraph.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
#include "ITStracking/ROFLookupTables.h"

namespace fs = std::filesystem;

using namespace o2::itsmft::tracking;

namespace
{

fs::path trackingRoot()
{
  return fs::path{__FILE__}.parent_path().parent_path();
}

std::string readFile(const fs::path& path)
{
  std::ifstream input{path};
  if (!input.good()) {
    return {};
  }
  return {std::istreambuf_iterator<char>{input}, {}};
}

bool isAdapterEdge(const fs::path& root, const fs::path& path)
{
  (void)root;
  (void)path;
  return false;
}

void scanCommonProductionSources(const fs::path& root)
{
  static constexpr std::array<std::string_view, 6> forbidden{
    "TrackingTopology<",
    "ROFOverlapTable<",
    "ROFVertexLookupTable<",
    "ROFMaskTable<",
    "checkSupportedNLayers",
    "IndexTableUtils<"};

  BOOST_REQUIRE_MESSAGE(fs::is_directory(root / "include"), "cannot inspect " << root.string());
  for (const auto& directory : {root / "include", root / "src"}) {
    for (const auto& entry : fs::recursive_directory_iterator(directory)) {
      if (!entry.is_regular_file()) {
        continue;
      }
      const auto extension = entry.path().extension().string();
      if (extension != ".h" && extension != ".hpp" && extension != ".c" && extension != ".cc" && extension != ".cxx") {
        continue;
      }
      const auto source = readFile(entry.path());
      BOOST_REQUIRE_MESSAGE(!source.empty(), "cannot inspect " << entry.path().string());
      if (isAdapterEdge(root, entry.path())) {
        BOOST_TEST_MESSAGE("adapter-edge exception: " << entry.path().lexically_relative(root).generic_string());
        continue;
      }
      for (const auto token : forbidden) {
        BOOST_CHECK_MESSAGE(source.find(token) == std::string::npos,
                            "forbidden M7c compatibility object in common production source " << entry.path().string() << ": " << token);
      }
    }
  }
  BOOST_CHECK_MESSAGE(!fs::exists(root / "include/ITSMFTTracking/TrackingTopology.h"),
                      "the deleted layer-indexed common TrackingTopology header was recreated");
}

template <int NLayers>
void checkRuntimeROFViewMatchesFrozenBuilder()
{
  o2::its::LayerTiming timing{};
  timing.mNROFsTF = 3;
  timing.mROFLength = 40;
  timing.mROFDelay = 7;
  timing.mROFBias = 11;
  timing.mROFAddTimeErr = 5;

  o2::its::ROFOverlapTable<NLayers> overlap;
  o2::its::ROFVertexLookupTable<NLayers> vertexLookup;
  for (int layer = 0; layer < NLayers; ++layer) {
    overlap.defineLayer(layer, timing);
    vertexLookup.defineLayer(layer, timing);
  }
  overlap.init();

  std::array<o2::its::Vertex, 2> vertices{};
  vertices[0].setTimeStamp(o2::its::TimeEstBC{16, 3});
  vertices[1].setTimeStamp(o2::its::TimeEstBC{95, 4});
  vertexLookup.init(vertices.data(), vertices.size());

  o2::its::ROFMaskTable<NLayers> mask{overlap};
  mask.resetMask();
  for (int layer = 0; layer < NLayers; ++layer) {
    mask.setROFEnabled(layer, 0, 1);
    if ((layer & 1) == 0) {
      mask.setROFEnabled(layer, 2, 1);
    }
  }

  const auto expectedOverlap = overlap.getView();
  const auto expectedVertexLookup = vertexLookup.getView();
  const auto expectedMask = mask.getView();
  const RuntimeROFViews runtime{expectedOverlap, expectedVertexLookup, expectedMask, {}};

  BOOST_CHECK_EQUAL(runtime.overlap.mLayerCount, NLayers);
  BOOST_CHECK_EQUAL(runtime.vertexLookup.mLayerCount, NLayers);
  BOOST_CHECK_EQUAL(runtime.mask.mLayerCount, NLayers);
  for (int layer = 0; layer < NLayers; ++layer) {
    const auto& expectedTiming = expectedOverlap.getLayer(layer);
    const auto& actualTiming = runtime.overlap.getLayer(layer);
    BOOST_CHECK_EQUAL(actualTiming.mNROFsTF, expectedTiming.mNROFsTF);
    BOOST_CHECK_EQUAL(actualTiming.mROFLength, expectedTiming.mROFLength);
    BOOST_CHECK_EQUAL(actualTiming.mROFDelay, expectedTiming.mROFDelay);
    BOOST_CHECK_EQUAL(actualTiming.mROFBias, expectedTiming.mROFBias);
    BOOST_CHECK_EQUAL(actualTiming.mROFAddTimeErr, expectedTiming.mROFAddTimeErr);
    BOOST_CHECK_EQUAL(actualTiming.getROFStartInBC(0), expectedTiming.getROFStartInBC(0));
    BOOST_CHECK_EQUAL(actualTiming.getROFEndInBC(2), expectedTiming.getROFEndInBC(2));
    BOOST_CHECK(actualTiming.getROFTimeBounds(0, true) == expectedTiming.getROFTimeBounds(0, true));
    BOOST_CHECK(actualTiming.getROFTimeBounds(2, true) == expectedTiming.getROFTimeBounds(2, true));
    BOOST_CHECK_EQUAL(runtime.mask.isROFEnabled(layer, 0), expectedMask.isROFEnabled(layer, 0));
    BOOST_CHECK_EQUAL(runtime.mask.isROFEnabled(layer, 1), expectedMask.isROFEnabled(layer, 1));
    BOOST_CHECK_EQUAL(runtime.mask.isROFEnabled(layer, 2), expectedMask.isROFEnabled(layer, 2));
    for (int rof = 0; rof < 3; ++rof) {
      const auto expectedVertices = expectedVertexLookup.getVertices(layer, rof);
      const auto actualVertices = runtime.vertexLookup.getVertices(layer, rof);
      BOOST_CHECK_EQUAL(actualVertices.getFirstEntry(), expectedVertices.getFirstEntry());
      BOOST_CHECK_EQUAL(actualVertices.getEntries(), expectedVertices.getEntries());
      for (const auto& vertex : vertices) {
        BOOST_CHECK_EQUAL(runtime.vertexLookup.isVertexCompatible(layer, rof, vertex),
                          expectedVertexLookup.isVertexCompatible(layer, rof, vertex));
      }
    }
  }

  // Includes empty, first, and last ROF pairs and exercises the non-zero
  // diamond timing envelope through both the overlap and vertex views.
  BOOST_CHECK_GT(runtime.overlap.getLayer(0).mROFAddTimeErr, 0u);
  for (int from = 0; from < NLayers; ++from) {
    for (int to = 0; to < NLayers; ++to) {
      for (int rof = 0; rof < 3; ++rof) {
        BOOST_CHECK_EQUAL(runtime.overlap.doROFsOverlap(from, rof, to, rof),
                          expectedOverlap.doROFsOverlap(from, rof, to, rof));
        if (runtime.overlap.doROFsOverlap(from, rof, to, rof)) {
          BOOST_CHECK(runtime.overlap.getTimeStamp(from, rof, to, rof) ==
                      expectedOverlap.getTimeStamp(from, rof, to, rof));
        }
      }
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(CommonProductionUsesOnlySparseTopologyAndRuntimeROFViews)
{
  scanCommonProductionSources(trackingRoot());
}

BOOST_AUTO_TEST_CASE(SparseTopologyViewRetainsExplicitNonIdentityOrder)
{
  SurfaceGraph topology{8};
  const auto first = topology.addTransition(SurfaceTransition{SurfaceId{5}, SurfaceId{2}, {}, 0});
  const auto second = topology.addTransition(SurfaceTransition{SurfaceId{2}, SurfaceId{7}, {}, 0});
  BOOST_REQUIRE(first.isValid());
  BOOST_REQUIRE(second.isValid());
  const auto cell = topology.addCell(first, second);
  BOOST_REQUIRE(cell.isValid());
  BOOST_REQUIRE(topology.finalize());

  const auto view = topology.getView();
  BOOST_REQUIRE_EQUAL(view.nTransitions, 2u);
  BOOST_REQUIRE_EQUAL(view.nCells, 1u);
  BOOST_CHECK(view.getTransition(first).from == SurfaceId{5});
  BOOST_CHECK(view.getTransition(first).to == SurfaceId{2});
  BOOST_CHECK(view.getTransition(second).from == SurfaceId{2});
  BOOST_CHECK(view.getTransition(second).to == SurfaceId{7});
  BOOST_CHECK(view.getCell(cell).hitSurfaces.has(SurfaceId{5}));
  BOOST_CHECK(view.getCell(cell).hitSurfaces.has(SurfaceId{2}));
  BOOST_CHECK(view.getCell(cell).hitSurfaces.has(SurfaceId{7}));
}

BOOST_AUTO_TEST_CASE(RuntimeROFViewsMatchITSAndMFTFixedBuilders)
{
  checkRuntimeROFViewMatchesFrozenBuilder<o2::itsmft::tracking::ITSNLayers>();
  checkRuntimeROFViewMatchesFrozenBuilder<o2::itsmft::tracking::MFTNLayers>();
}

BOOST_AUTO_TEST_CASE(ScratchResetInvalidatesRuntimeROFContext)
{
  o2::its::LayerTiming timing{};
  timing.mNROFsTF = 1;
  timing.mROFLength = 40;
  o2::its::ROFOverlapTable<o2::itsmft::tracking::ITSNLayers> overlap;
  o2::its::ROFVertexLookupTable<o2::itsmft::tracking::ITSNLayers> vertexLookup;
  for (int layer = 0; layer < o2::itsmft::tracking::ITSNLayers; ++layer) {
    overlap.defineLayer(layer, timing);
    vertexLookup.defineLayer(layer, timing);
  }
  overlap.init();
  vertexLookup.init();
  o2::its::ROFMaskTable<o2::itsmft::tracking::ITSNLayers> mask{overlap};
  mask.setROFsEnabled(0, 0, 1, 1);

  SurfaceTrackingScratch scratch;
  scratch.setROFViews(RuntimeROFViews{overlap.getView(), vertexLookup.getView(), mask.getView(), mask.getView()});
  scratch.useUPCMask();
  BOOST_CHECK_EQUAL(scratch.getROFOverlapView().mLayerCount, o2::itsmft::tracking::ITSNLayers);
  BOOST_CHECK(scratch.getROFMaskView().mFlatMask != nullptr);

  scratch.reset();
  const auto& views = scratch.getROFViews();
  BOOST_CHECK_EQUAL(views.overlap.mLayerCount, 0);
  BOOST_CHECK_EQUAL(views.vertexLookup.mLayerCount, 0);
  BOOST_CHECK_EQUAL(views.mask.mLayerCount, 0);
  BOOST_CHECK_EQUAL(views.upcMask.mLayerCount, 0);
  BOOST_CHECK(views.overlap.mLayers == nullptr);
  BOOST_CHECK(views.vertexLookup.mLayers == nullptr);
  BOOST_CHECK(views.mask.mFlatMask == nullptr);
  BOOST_CHECK(views.upcMask.mFlatMask == nullptr);
}

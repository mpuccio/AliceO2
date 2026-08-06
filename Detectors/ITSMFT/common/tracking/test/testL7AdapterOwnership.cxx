// Gate 4 L7: detector timing tables and publication compatibility belong to
// workflow/adapter contexts; the common core receives borrowed event views.

#define BOOST_TEST_MODULE ITSMFT L7 adapter ownership
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/ClockTimingPublicationView.h"
#include "ITSMFTTracking/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/ROFViews.h"
#include "ITSMFTTracking/SurfaceTiming.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITStracking/ROFLookupTables.h"

namespace fs = std::filesystem;
using o2::its::LayerTiming;
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

template <typename Table>
void defineUniformTiming(Table& table, const std::array<LayerTiming, ITSNLayers>& timing)
{
  for (int layer = 0; layer < ITSNLayers; ++layer) {
    table.defineLayer(layer, timing[layer]);
  }
  table.init();
}

BOOST_AUTO_TEST_CASE(generic_core_has_no_detector_timing_or_publication_ownership)
{
  const auto root = trackingRoot();
  static constexpr std::array<std::string_view, 6> coreFiles{
    "include/ITSMFTTracking/TimeFrame.h",
    "include/ITSMFTTracking/MultiSourceTimeFrameLoader.h",
    "include/ITSMFTTracking/Tracker.h",
    "include/ITSMFTTracking/TrackerTraits.h",
    "include/ITSMFTTracking/TrackingInterface.h",
    "include/ITSMFTTracking/SurfacePlanTrackingParticipant.h",
  };
  static constexpr std::array<std::string_view, 10> forbiddenOwnership{
    "ROFOverlapTable<",
    "ROFVertexLookupTable<",
    "ROFMaskTable<",
    "mPublicationClock",
    "mPublicationValid",
    "mMFTROFrameLengthInBC",
    "ITSSharedClusterCompatibilityOwner",
    "MFTPublicationCompatibilityOwner",
    "mMultiplicityMask",
    "mUPCMask",
  };

  for (const auto relative : coreFiles) {
    const auto path = root / relative;
    const auto source = readFile(path);
    BOOST_REQUIRE_MESSAGE(!source.empty(), "cannot inspect " << path.string());
    for (const auto token : forbiddenOwnership) {
      BOOST_CHECK_MESSAGE(source.find(token) == std::string::npos,
                          relative << " retains detector/workflow ownership " << token);
    }
  }

  const auto itsWorkflow = readFile(root.parent_path().parent_path() / "ITS/workflow-ca/include/ITSCAWorkflow/CATrackerSpec.h");
  const auto mftWorkflow = readFile(root.parent_path().parent_path() / "MFT/workflow/include/MFTWorkflow/CATrackerSpec.h");
  BOOST_REQUIRE(!itsWorkflow.empty());
  BOOST_REQUIRE(!mftWorkflow.empty());
  BOOST_CHECK(itsWorkflow.find("ROFOverlapTable<") != std::string::npos);
  BOOST_CHECK(itsWorkflow.find("mPublicationClock") != std::string::npos);
  BOOST_CHECK(mftWorkflow.find("ROFOverlapTable<") != std::string::npos);
  BOOST_CHECK(mftWorkflow.find("mPublicationClock") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(workflow_owned_runtime_views_preserve_timing_and_mask_semantics)
{
  std::array<LayerTiming, ITSNLayers> timing{};
  for (auto& layer : timing) {
    layer = LayerTiming{.mNROFsTF = 4, .mROFLength = 40, .mROFDelay = 3, .mROFBias = 5, .mROFAddTimeErr = 2};
  }

  o2::its::ROFOverlapTable<ITSNLayers> overlap;
  o2::its::ROFVertexLookupTable<ITSNLayers> vertex;
  defineUniformTiming(overlap, timing);
  defineUniformTiming(vertex, timing);
  o2::its::ROFMaskTable<ITSNLayers> mask{overlap};
  mask.resetMask();
  mask.setROFsEnabled(0, 1, 2, 1);

  RuntimeROFViews views{overlap.getView(), vertex.getView(), mask.getView(), {}};
  BOOST_REQUIRE_EQUAL(views.overlap.mLayerCount, ITSNLayers);
  BOOST_CHECK_EQUAL(views.overlap.getClockLayer().getROFStartInBC(0), 8);
  BOOST_CHECK_EQUAL(views.overlap.getClockLayer().getROFStartInBC(3), 128);
  BOOST_CHECK(!views.mask.isROFEnabled(0, 0));
  BOOST_CHECK(views.mask.isROFEnabled(0, 1));
  BOOST_CHECK(views.mask.isROFEnabled(0, 2));
  BOOST_CHECK(!views.mask.isROFEnabled(0, 3));

  ClockTimingPublicationView publicationClock{views.overlap.getClockLayer()};
  CommonTrackTimestamp timestamp{8, 48};
  const auto output = publicationClock.makeOutputTimestamp(timestamp);
  BOOST_REQUIRE(output.has_value());
  BOOST_CHECK_EQUAL(publicationClock.getROF(*output), 0);
  BOOST_CHECK_EQUAL(publicationClock.getROFCount(), 4);

  const auto borrowedFlatTable = views.overlap.mFlatTable;
  views = {};
  BOOST_CHECK_EQUAL(views.overlap.mLayerCount, 0);
  BOOST_CHECK_EQUAL(views.overlap.mFlatTable, nullptr);
  BOOST_CHECK(borrowedFlatTable != nullptr);
}

BOOST_AUTO_TEST_CASE(publication_sidecar_is_adapter_local_and_resettable)
{
  ITSSharedClusterCompatibility sidecar;
  sidecar.clear();
  BOOST_CHECK(sidecar.entries().empty());
  BOOST_CHECK(!sidecar.isSealed());

  MFTPublicationCompatibility mftSidecar;
  mftSidecar.clear();
  BOOST_CHECK(mftSidecar.entries().empty());
}

} // namespace

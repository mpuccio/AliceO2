// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

// Gate 4 M7d: the shared Tracker/TrackerTraits orchestration is one
// non-templated runtime-plan core. ITS/MFT templates are permitted only at
// the adapter edge; this test also executes the same core with four sparse,
// non-identity ordered surfaces.

#define BOOST_TEST_MODULE ITSMFT M7d non - templated Tracker guard
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IndexTableConfiguration.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/detail/TimeFrameScratch.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITStracking/Constants.h"

#include "TraversalTestSupport.h"

using namespace o2::itsmft::tracking;

namespace
{

template <typename T>
concept HasParametersByKindMember = requires(T value) { value.parametersByKind; };

template <typename T>
concept HasParametersByKindAccessor = requires(const T& value) { value.getTrackingParametersByKind(); };

static_assert(!HasParametersByKindMember<TrackerIterationConfiguration>);
static_assert(!HasParametersByKindAccessor<TimeFrame>);

std::string readFile(const std::filesystem::path& path)
{
  std::ifstream input{path};
  BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path.string());
  return {std::istreambuf_iterator<char>{input}, {}};
}

std::string withoutCommentsAndPreprocessor(std::string_view source)
{
  std::string result;
  result.reserve(source.size());
  bool blockComment = false;
  std::istringstream lines{std::string{source}};
  std::string line;
  while (std::getline(lines, line)) {
    if (!line.empty() && line.front() == '#') {
      continue;
    }
    std::string clean;
    for (std::size_t i = 0; i < line.size();) {
      if (blockComment) {
        if (i + 1 < line.size() && line[i] == '*' && line[i + 1] == '/') {
          blockComment = false;
          i += 2;
        } else {
          ++i;
        }
      } else if (i + 1 < line.size() && line[i] == '/' && line[i + 1] == '*') {
        blockComment = true;
        i += 2;
      } else if (i + 1 < line.size() && line[i] == '/' && line[i + 1] == '/') {
        break;
      } else {
        clean.push_back(line[i++]);
      }
    }
    result += clean;
    result.push_back('\n');
  }
  return result;
}

std::filesystem::path trackingRoot()
{
  return std::filesystem::path{__FILE__}.parent_path().parent_path();
}

bool noopSeedRefit(const TrackSeed&,
                   const o2::itsmft::TrackingParameters&,
                   float,
                   TimeFrameScratch&,
                   gsl::span<const gsl::span<const GlobalMeasurement>>,
                   gsl::span<const gsl::span<const SurfaceMeasurement>>,
                   SurfaceCatalogView,
                   gsl::span<const LayerId>,
                   TrackingCandidate&)
{
  return false;
}

BOOST_AUTO_TEST_CASE(CommonProductionHasOneNonTemplatedTrackerCore)
{
  const auto root = trackingRoot();
  const std::regex forbidden{R"(\bTrackerTraits\s*<|\bTracker\s*<|\bCATracker\b|\bTrackerITS\b|\bTrackerMFT\b|\bparametersByKind\b|\bparametersForKind\b|\bmTrkParamsByKind\b|array\s*<\s*(?:o2::itsmft::)?TrackingParameters\s*,\s*2\s*>)"};
  for (const auto& tree : {root / "include", root / "src"}) {
    for (const auto& entry : std::filesystem::recursive_directory_iterator(tree)) {
      if (!entry.is_regular_file()) {
        continue;
      }
      const auto extension = entry.path().extension().string();
      if (extension != ".h" && extension != ".hpp" && extension != ".c" && extension != ".cc" && extension != ".cxx") {
        continue;
      }
      const auto code = withoutCommentsAndPreprocessor(readFile(entry.path()));
      BOOST_CHECK_MESSAGE(!std::regex_search(code, forbidden),
                          "retired templated/compatibility Tracker spelling remains in " << entry.path().string());
    }
  }
}

BOOST_AUTO_TEST_CASE(GenericCoreHeadersHaveNoApplicationDependencies)
{
  const auto root = trackingRoot();
  const std::array<std::filesystem::path, 2> headers{
    root / "include/ITSMFTTracking/Tracker.h",
    root / "include/ITSMFTTracking/TrackerTraits.h"};
  const std::array<std::string_view, 7> forbidden{
    "DetID", "DPL", "Workflow", "Writer", "ITS", "MFT", "SurfaceKindPair"};
  for (const auto& header : headers) {
    const auto code = withoutCommentsAndPreprocessor(readFile(header));
    for (const auto token : forbidden) {
      BOOST_CHECK_MESSAGE(code.find(token) == std::string::npos,
                          "generic core header " << header.filename().string() << " names application dependency " << token);
    }
  }
}

BOOST_AUTO_TEST_CASE(NonSevenOrTenPlanExecutesTheNonTemplatedCore)
{
  constexpr std::array<uint16_t, 4> orderedIds{1, 4, 7, 10};
  std::vector<SurfaceDescriptor> catalog;
  catalog.reserve(12);
  for (uint16_t id = 0; id < 12; ++id) {
    SurfaceDescriptor surface{LayerId{id}, id, static_cast<uint8_t>(o2::detectors::DetID::MFT), SurfaceKind::Disk,
                              0, static_cast<float>(id + 1), 0.f, 100.f};
    const auto position = std::find(orderedIds.begin(), orderedIds.end(), id);
    if (position != orderedIds.end()) {
      const auto layer = static_cast<std::size_t>(position - orderedIds.begin());
      surface.material.xOverX0 = kNominalMFTLayerX0[layer];
      surface.material.arealDensityGPerCm2 = kNominalMFTLayerX0[layer] * o2::its::constants::Radl * o2::its::constants::Rho;
    }
    catalog.push_back(surface);
  }

  o2::itsmft::TrackingParameters parameters;
  o2::itsmft::resetDetectorDefaults(parameters, o2::detectors::DetID::MFT);
  parameters.NLayers = 4;
  parameters.StartLayerMask = 0x0f;
  parameters.HoleLayerMask = 0;
  parameters.MaxHoles = 0;
  for (std::size_t position = 0; position < orderedIds.size(); ++position) {
    parameters.LayerxX0[position] = catalog[orderedIds[position]].material.xOverX0;
  }
  parameters.LayerxX0.resize(orderedIds.size());
  const std::vector<LayerId> ordered{LayerId{1}, LayerId{4}, LayerId{7}, LayerId{10}};
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  Tracker tracker{noopSeedRefit};
  TrackerInitialization configuration;
  configuration.catalog = {catalog.data(), static_cast<uint32_t>(catalog.size())};
  configuration.memoryPool = pool;
  TrackerIterationConfiguration iteration;
  iteration.layout = makeSurfaceLayoutChain(ordered, parameters.MaxHoles,
                                            positionalSurfaceMask(parameters.HoleLayerMask, ordered, ordered.size()),
                                            positionalSurfaceMask(parameters.StartLayerMask, ordered, ordered.size()));
  iteration.parameters = parameters;
  configuration.iterations.push_back(std::move(iteration));
  const auto configured = tracker.initialize(frame, configuration);
  BOOST_REQUIRE(configured.ok());

  const auto& layout = frame.getLayout(0);
  auto& workspace = frame.getWorkspace().getTraversalWorkspace(0);
  TrackerTestAccess::preparePlan(tracker, workspace, layout);
  const auto topology = workspace.getTopologyView();
  TrackerTraits traits;
  std::shared_ptr<tbb::task_arena> arena;
  traits.setNThreads(1, arena);
  BOOST_REQUIRE_EQUAL(frame.getWorkspace().getNOwnedSurfaces(), 4u);
  BOOST_CHECK(topology.orderedSurfaces[0] == LayerId{1});
  BOOST_CHECK(topology.orderedSurfaces[2] == LayerId{7});

  TrackSeed seed;
  SurfaceMask activePositions;
  for (uint16_t position = 0; position < ordered.size(); ++position) {
    activePositions.set(LayerId{position});
    seed.setCluster(position, static_cast<int>(100 + position));
  }
  seed.setSurfaceMask(activePositions);
  BOOST_CHECK_EQUAL(seed.getActiveSurfaceCount(), 4);
  BOOST_CHECK_EQUAL(seed.getCluster(0), 100);
  BOOST_CHECK_EQUAL(seed.getCluster(3), 103);
  // Static layout configuration retains the sparse, non-identity order. The
  // Tracker derives the pass-local workspace topology from it at traversal time.
  BOOST_CHECK_EQUAL(topology.nOrderedSurfaces, ordered.size());
}

} // namespace

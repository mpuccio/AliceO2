// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

// M6f (Detectors/ITSMFT/common/tracking/doc/design/0002-m6-generic-workspace-migration.md
// Sec 9, 10): focused coverage for the sole detector-neutral
// TimeFrameScratch after the temporary workspace/binding bridge was
// retired. This file keeps container, dependency-boundary, and deleted-file
// assertions; production traffic and replay coverage live in the participant,
// workflow, and validation tests.

#define BOOST_TEST_MODULE ITSMFT TimeFrameScratch
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <memory>
#include <regex>
#include <string>
#include <vector>

#include "ITSMFTTracking/detail/TimeFrameScratch.h"
#include "ITSMFTTracking/TraversalTopology.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITStracking/BoundedAllocator.h"

namespace fs = std::filesystem;

namespace
{
using namespace o2::itsmft::tracking;
using o2::itsmft::TrackingParameters;

std::vector<LayerId> ordered(uint16_t first, uint16_t count)
{
  std::vector<LayerId> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(LayerId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

SurfaceDescriptor surfaceWithOwner(uint16_t id, SurfaceKind kind, uint8_t detectorId)
{
  return SurfaceDescriptor{id, detectorId, kind};
}

// A synthetic, detector-neutral Cylinder chain of `n` surfaces. Deliberately
// uses a synthetic detectorId (250), never
// o2::detectors::DetID::ITS/MFT, and never touches TimeFrameScratch
// itself with a detector identity of any kind -- only the three plain
// counts exposed directly by the graph.
struct SyntheticChain {
  std::vector<SurfaceDescriptor> surfaces;
  DetectorLayout layout;
  TraversalTopology topology;

  explicit SyntheticChain(uint16_t n)
    : surfaces{makeSurfaces(n)}, layout{buildLayout(surfaces)}, topology{buildTopology(layout)}
  {
  }

  std::size_t nOwnedSurfaces() const { return layout.size(); }
  std::size_t nEdges() const { return topology.edges.size(); }
  std::size_t nCells() const { return topology.paths.size(); }

 private:
  static std::vector<SurfaceDescriptor> makeSurfaces(uint16_t n)
  {
    std::vector<SurfaceDescriptor> result;
    result.reserve(n);
    for (uint16_t id = 0; id < n; ++id) {
      result.push_back(surfaceWithOwner(id, SurfaceKind::Cylinder, 250));
    }
    return result;
  }
  static DetectorLayout buildLayout(const std::vector<SurfaceDescriptor>& surfaces)
  {
    const auto orderedSurfaces = ordered(0, static_cast<uint16_t>(surfaces.size()));
    return DetectorLayout{gsl::span<const SurfaceDescriptor>{surfaces.data(), surfaces.size()},
                          makeDetectorLayout()};
  }

  static TraversalTopology buildTopology(const DetectorLayout& layout)
  {
    TrackingParameters parameters;
    parameters.NLayers = static_cast<int>(layout.size());
    parameters.StartLayerMask = LayerMask{1u << (parameters.NLayers - 1)};
    auto built = deriveTraversalTopology(layout, parameters);
    BOOST_REQUIRE(built.ok());
    return std::move(*built.topology);
  }
};

std::shared_ptr<o2::its::BoundedMemoryResource> makePool()
{
  return std::make_shared<o2::its::BoundedMemoryResource>();
}

} // namespace

BOOST_AUTO_TEST_CASE(ConfiguresStorageWithDistinctRuntimeCountsWithoutDetectorOrLayerCountAssumption)
{
  // Two synthetic plans of deliberately different owned-surface/edge/
  // cell cardinality -- 4 surfaces (3 edges, 2 cells) and 6 surfaces (5
  // edges, 4 cells) -- adopted in turn by the *same* scratch instance,
  // proving sizing tracks whatever plan was last adopted rather than any
  // baked-in constant (no NLayers, no ITS/MFT-shaped assumption anywhere in
  // this test or in TimeFrameScratch itself).
  SyntheticChain small{4};
  SyntheticChain large{6};
  BOOST_CHECK_EQUAL(small.nOwnedSurfaces(), 4u);
  BOOST_CHECK_EQUAL(small.nEdges(), 3u);
  BOOST_CHECK_EQUAL(small.nCells(), 2u);
  BOOST_CHECK_EQUAL(large.nOwnedSurfaces(), 6u);
  BOOST_CHECK_EQUAL(large.nEdges(), 5u);
  BOOST_CHECK_EQUAL(large.nCells(), 4u);
  BOOST_REQUIRE_NE(small.nOwnedSurfaces(), large.nOwnedSurfaces());

  TimeFrameScratch scratch;
  scratch.setMemoryPool(makePool());

  scratch.configureStorage(small.nEdges(), small.nCells());
  BOOST_CHECK_EQUAL(scratch.getNEdges(), 3u);
  BOOST_CHECK_EQUAL(scratch.getNCells(), 2u);
  BOOST_CHECK_EQUAL(scratch.mTracklets.size(), 3u);
  BOOST_CHECK_EQUAL(scratch.mCells.size(), 2u);

  // Reconfiguring differently-shaped storage on the same instance must fully
  // re-size every container -- no residue from the previous plan.
  scratch.configureStorage(large.nEdges(), large.nCells());
  BOOST_CHECK_EQUAL(scratch.getNEdges(), 5u);
  BOOST_CHECK_EQUAL(scratch.getNCells(), 4u);
  BOOST_CHECK_EQUAL(scratch.mTracklets.size(), 5u);
  BOOST_CHECK_EQUAL(scratch.mCells.size(), 4u);
}

BOOST_AUTO_TEST_CASE(PerSurfaceAndEdgeCellContainersHaveExpectedRuntimeSizes)
{
  SyntheticChain chain{5};
  TimeFrameScratch scratch;
  scratch.setMemoryPool(makePool());
  scratch.configureStorage(chain.nEdges(), chain.nCells());

  const auto nTr = chain.nEdges();
  const auto nCe = chain.nCells();

  // Group B: sparse edge/cell counts.
  BOOST_CHECK_EQUAL(scratch.mTracklets.size(), nTr);
  BOOST_CHECK_EQUAL(scratch.mTrackletsLookupTable.size(), nTr);
  BOOST_CHECK_EQUAL(scratch.mTrackletLabels.size(), nTr);
  BOOST_CHECK_EQUAL(scratch.mEdgePhiCuts.size(), nTr);
  BOOST_CHECK_EQUAL(scratch.mEdgeMSAngles.size(), nTr);
  BOOST_CHECK_EQUAL(scratch.mCells.size(), nCe);
  BOOST_CHECK_EQUAL(scratch.mCellsLookupTable.size(), nCe);
  BOOST_CHECK_EQUAL(scratch.mCellsNeighbours.size(), nCe);
  BOOST_CHECK_EQUAL(scratch.mCellsNeighboursTopology.size(), nCe);
  BOOST_CHECK_EQUAL(scratch.mCellsNeighboursLUT.size(), nCe);
  BOOST_CHECK_EQUAL(scratch.mCellLabels.size(), nCe);
}

BOOST_AUTO_TEST_CASE(ResetClearsWorkingStateWithoutMutatingAPopulatedTimeFrameOrTheAdoptedPlanSize)
{
  SyntheticChain chain{4};
  TimeFrameScratch scratch;
  scratch.setMemoryPool(makePool());
  scratch.configureStorage(chain.nEdges(), chain.nCells());

  // Populate a handful of containers with observable content.
  scratch.mTracklets[0].emplace_back();
  scratch.mCells[0].emplace_back();
  scratch.getEdgePhiCuts()[0] = 7.f;

  // An unrelated, populated TimeFrame: reset() takes no TimeFrame parameter
  // at all, so it structurally cannot reach it -- this proves that
  // structural guarantee holds at runtime too.
  TimeFrame frame;
  frame.setBz(5.f);
  frame.setBeamPosition(1.f, 2.f, 0.1f);

  scratch.reset();

  // Vector-of-bounded_vector containers: outer element
  // count -- the configured storage size -- survives reset(); only each element's
  // *contents* are cleared. Mirrors TimeFrameScratch::reset()
  // exactly (it never shrinks its own NLayers-wide outer arrays either).
  BOOST_CHECK_EQUAL(scratch.mTracklets.size(), chain.nEdges());
  BOOST_CHECK(scratch.mTracklets[0].empty());
  BOOST_CHECK_EQUAL(scratch.mCells.size(), chain.nCells());
  BOOST_CHECK(scratch.mCells[0].empty());
  BOOST_CHECK_EQUAL(scratch.getNEdges(), chain.nEdges());
  BOOST_CHECK_EQUAL(scratch.getNCells(), chain.nCells());

  // Flat bounded_vector containers are fully cleared to empty.
  BOOST_CHECK(scratch.getEdgePhiCuts().empty());

  BOOST_CHECK_EQUAL(frame.getBz(), 5.f);
  BOOST_CHECK_EQUAL(frame.getBeamX(), 1.f);
}

namespace
{
void checkNoForbiddenToken(const fs::path& path, const std::string& token, const std::string& description)
{
  std::ifstream input{path};
  BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path.string());
  const std::regex tokenRegex{"\\b" + token + "\\b"};
  std::string line;
  size_t lineNumber = 0;
  while (std::getline(input, line)) {
    ++lineNumber;
    BOOST_CHECK_MESSAGE(!std::regex_search(line, tokenRegex),
                        path.filename().string() << ":" << lineNumber << " mentions " << token << " (" << description << "): " << line);
  }
}
} // namespace

// Scratch is deliberately unaware of TimeFrame, detector configuration and
// workflow/output types. Tracker translates those domains into plain sizes.
BOOST_AUTO_TEST_CASE(NewHeadersPullNoITSWorkflowOutputOrSurfaceKindDependency)
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  const fs::path trackingRoot = fs::path(testDirectory) / "..";

  const std::vector<fs::path> newFiles = {
    trackingRoot / "include/ITSMFTTracking/detail/TimeFrameScratch.h",
    trackingRoot / "src/TimeFrameScratch.cxx"};
  for (const auto& path : newFiles) {
    BOOST_REQUIRE_MESSAGE(fs::is_regular_file(path), "cannot find " << path.string());
  }

  const std::vector<std::pair<std::string, std::string>> forbidden = {
    {"SurfaceKind", "detail/-confined dispatch-key type"},
    {"SurfaceKind", "Barrel/Forward state-family selector"},
    {"DPL", "DPL workflow machinery"},
    {"Workflow", "workflow-layer naming"},
    {"workflow", "workflow-layer naming"},
  };
  for (const auto& path : newFiles) {
    for (const auto& [token, description] : forbidden) {
      checkNoForbiddenToken(path, token, description);
    }
  }

  const std::vector<std::string> forbiddenIncludes = {
    "DetectorTraversalBinding.h",
    "IOUtils.h",
    "LegacyTrackerScratch.h",
    "TimeFrame.h",
    "SurfacePlanTrackingParticipant.h"};
  for (const auto& path : newFiles) {
    std::ifstream input{path};
    BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path.string());
    std::string line;
    size_t lineNumber = 0;
    while (std::getline(input, line)) {
      ++lineNumber;
      for (const auto& forbiddenInclude : forbiddenIncludes) {
        BOOST_CHECK_MESSAGE(line.find(forbiddenInclude) == std::string::npos,
                            path.filename().string() << ":" << lineNumber << " must not include " << forbiddenInclude << ": " << line);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(RetiredLegacyWorkspaceFilesAreAbsentAfterMigration)
{
  // M6f removes the temporary workspace/binding bridge once every production
  // participant uses the plan-driven model. Keep the filesystem assertion
  // here alongside the source guard so an accidental reintroduction is
  // caught by both the focused scratch test and the migration guard.
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  const fs::path trackingRoot = fs::path(testDirectory) / "..";
  BOOST_CHECK(!fs::exists(trackingRoot / "include/ITSMFTTracking/LegacyTrackerScratch.h"));
  BOOST_CHECK(!fs::exists(trackingRoot / "include/ITSMFTTracking/detail/DetectorTraversalBinding.h"));
  BOOST_CHECK(!fs::exists(trackingRoot / "src/LegacyTrackerScratch.cxx"));
  BOOST_CHECK(fs::is_regular_file(trackingRoot / "include/ITSMFTTracking/detail/TimeFrameScratch.h"));
  BOOST_CHECK(!fs::exists(trackingRoot / "include/ITSMFTTracking/detail/SurfacePlanBinding.h"));
}

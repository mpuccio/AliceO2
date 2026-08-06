// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

// M6f (Detectors/ITSMFT/common/tracking/doc/design/0002-m6-generic-workspace-migration.md
// Sec 9, 10): focused coverage for the sole detector-neutral
// SurfaceTrackingScratch after the temporary workspace/binding bridge was
// retired. This file keeps container, dependency-boundary, and deleted-file
// assertions; production traffic and replay coverage live in the participant,
// workflow, and validation tests.

#define BOOST_TEST_MODULE ITSMFT SurfaceTrackingScratch
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

#include "ITSMFTTracking/DetectorLayoutBuilder.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "ITStracking/BoundedAllocator.h"

namespace fs = std::filesystem;

namespace
{
using namespace o2::itsmft::tracking;

std::vector<SurfaceId> ordered(uint16_t first, uint16_t count)
{
  std::vector<SurfaceId> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

SurfaceDescriptor surfaceWithOwner(uint16_t id, SurfaceKind kind, uint8_t detectorId)
{
  return SurfaceDescriptor{SurfaceId{id}, id, detectorId, kind};
}

// A synthetic, detector-neutral Cylinder chain of `n` surfaces -- same
// construction technique as testSurfacePlanBinding.cxx's own synthetic
// fixtures (this file deliberately duplicates it locally rather than
// sharing a new header, matching that file's own stated per-file-local
// fixture convention). Deliberately uses a synthetic detectorId (250), never
// o2::detectors::DetID::ITS/MFT, and never touches SurfaceTrackingScratch
// itself with a detector identity of any kind -- only the three plain
// counts a SurfacePlanBinding already exposes.
struct SyntheticChain {
  std::vector<SurfaceDescriptor> surfaces;
  DetectorLayout layout;
  DetectorLayoutView view;
  SurfacePlanBinding::BuildResult binding;

  explicit SyntheticChain(uint16_t n)
    : surfaces{makeSurfaces(n)}, layout{build(surfaces)}, view{layout.getView(surfaces, allCylinder(n), SurfaceMask{})}, binding{buildBinding(view, n)}
  {
    BOOST_REQUIRE(binding.ok());
  }

  std::size_t nOwnedSurfaces() const { return static_cast<std::size_t>(binding.binding->getOwnedSurfaces().count()); }
  std::size_t nTransitions() const { return binding.binding->getGlobalTransitions().size(); }
  std::size_t nCells() const { return binding.binding->getGlobalCells().size(); }

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
  static SurfaceMask allCylinder(uint16_t n)
  {
    SurfaceMask mask;
    for (uint16_t id = 0; id < n; ++id) {
      mask.set(SurfaceId{id});
    }
    return mask;
  }
  static DetectorLayout build(const std::vector<SurfaceDescriptor>& surfaces)
  {
    DetectorLayoutBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}};
    SurfaceMask seed;
    seed.set(SurfaceId{static_cast<uint16_t>(surfaces.size() - 1)});
    auto built = builder.addSubgraph({ordered(0, static_cast<uint16_t>(surfaces.size())), 0, SurfaceMask{}, seed}).build();
    BOOST_REQUIRE(built.ok());
    return std::move(*built.layout);
  }
  static SurfacePlanBinding::BuildResult buildBinding(const DetectorLayoutView& view, uint16_t n)
  {
    SurfaceMask owned;
    for (uint16_t id = 0; id < n; ++id) {
      owned.set(SurfaceId{id});
    }
    return SurfacePlanBinding::build(view, ClusterSourceId{0}, owned, ordered(0, n), SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  }
};

std::shared_ptr<o2::its::BoundedMemoryResource> makePool()
{
  return std::make_shared<o2::its::BoundedMemoryResource>();
}

} // namespace

BOOST_AUTO_TEST_CASE(AdoptsPlansWithDistinctRuntimeCountsWithoutDetectorOrLayerCountAssumption)
{
  // Two synthetic plans of deliberately different owned-surface/transition/
  // cell cardinality -- 4 surfaces (3 transitions, 2 cells, per
  // testSurfacePlanBinding.cxx's own equivalent fixture) and 6 surfaces (5
  // transitions, 4 cells) -- adopted in turn by the *same* scratch instance,
  // proving sizing tracks whatever plan was last adopted rather than any
  // baked-in constant (no NLayers, no ITS/MFT-shaped assumption anywhere in
  // this test or in SurfaceTrackingScratch itself).
  SyntheticChain small{4};
  SyntheticChain large{6};
  BOOST_CHECK_EQUAL(small.nOwnedSurfaces(), 4u);
  BOOST_CHECK_EQUAL(small.nTransitions(), 3u);
  BOOST_CHECK_EQUAL(small.nCells(), 2u);
  BOOST_CHECK_EQUAL(large.nOwnedSurfaces(), 6u);
  BOOST_CHECK_EQUAL(large.nTransitions(), 5u);
  BOOST_CHECK_EQUAL(large.nCells(), 4u);
  BOOST_REQUIRE_NE(small.nOwnedSurfaces(), large.nOwnedSurfaces());

  SurfaceTrackingScratch scratch;
  scratch.setMemoryPool(makePool());

  scratch.adoptPlan(small.nOwnedSurfaces(), small.nTransitions(), small.nCells());
  BOOST_CHECK_EQUAL(scratch.getNOwnedSurfaces(), 4u);
  BOOST_CHECK_EQUAL(scratch.getNTransitions(), 3u);
  BOOST_CHECK_EQUAL(scratch.getNCells(), 2u);
  BOOST_CHECK_EQUAL(scratch.mClusters.size(), 4u);
  BOOST_CHECK_EQUAL(scratch.mTracklets.size(), 3u);
  BOOST_CHECK_EQUAL(scratch.mCells.size(), 2u);

  // Re-adopting a differently-shaped plan on the same instance must fully
  // re-size every container -- no residue from the previous plan.
  scratch.adoptPlan(large.nOwnedSurfaces(), large.nTransitions(), large.nCells());
  BOOST_CHECK_EQUAL(scratch.getNOwnedSurfaces(), 6u);
  BOOST_CHECK_EQUAL(scratch.getNTransitions(), 5u);
  BOOST_CHECK_EQUAL(scratch.getNCells(), 4u);
  BOOST_CHECK_EQUAL(scratch.mClusters.size(), 6u);
  BOOST_CHECK_EQUAL(scratch.mTracklets.size(), 5u);
  BOOST_CHECK_EQUAL(scratch.mCells.size(), 4u);
}

BOOST_AUTO_TEST_CASE(PerSurfaceAndTransitionCellContainersHaveExpectedRuntimeSizes)
{
  SyntheticChain chain{5};
  SurfaceTrackingScratch scratch;
  scratch.setMemoryPool(makePool());
  scratch.adoptPlan(chain.nOwnedSurfaces(), chain.nTransitions(), chain.nCells());

  const auto nSurf = chain.nOwnedSurfaces();
  const auto nTr = chain.nTransitions();
  const auto nCe = chain.nCells();

  // Group A: one slot per owned surface.
  BOOST_CHECK_EQUAL(scratch.mClusters.size(), nSurf);
  BOOST_CHECK_EQUAL(scratch.mUnsortedClusters.size(), nSurf);
  BOOST_CHECK_EQUAL(scratch.mTrackingFrameInfo.size(), nSurf);
  BOOST_CHECK_EQUAL(scratch.mClusterExternalIndices.size(), nSurf);
  BOOST_CHECK_EQUAL(scratch.mClusterSize.size(), nSurf);
  BOOST_CHECK_EQUAL(scratch.mROFramesClusters.size(), nSurf);
  BOOST_CHECK_EQUAL(scratch.mClusterLabels.size(), nSurf);
  BOOST_CHECK_EQUAL(scratch.mIndexTables.size(), nSurf);
  BOOST_CHECK_EQUAL(scratch.mUsedClusters.size(), nSurf);
  BOOST_CHECK_EQUAL(scratch.mNClustersPerROF.size(), nSurf);
  BOOST_CHECK_EQUAL(scratch.mMinR.size(), nSurf);
  BOOST_CHECK_EQUAL(scratch.mMaxR.size(), nSurf);
  BOOST_CHECK_EQUAL(scratch.mBogusClusters.size(), nSurf);
  BOOST_CHECK_EQUAL(scratch.mPositionResolution.size(), nSurf);

  // Group B: sparse transition/cell counts.
  BOOST_CHECK_EQUAL(scratch.mTracklets.size(), nTr);
  BOOST_CHECK_EQUAL(scratch.mTrackletsLookupTable.size(), nTr);
  BOOST_CHECK_EQUAL(scratch.mTrackletLabels.size(), nTr);
  BOOST_CHECK_EQUAL(scratch.mTransitionPhiCuts.size(), nTr);
  BOOST_CHECK_EQUAL(scratch.mTransitionMSAngles.size(), nTr);
  BOOST_CHECK_EQUAL(scratch.mCells.size(), nCe);
  BOOST_CHECK_EQUAL(scratch.mCellsLookupTable.size(), nCe);
  BOOST_CHECK_EQUAL(scratch.mCellsNeighbours.size(), nCe);
  BOOST_CHECK_EQUAL(scratch.mCellsNeighboursTopology.size(), nCe);
  BOOST_CHECK_EQUAL(scratch.mCellsNeighboursLUT.size(), nCe);
  BOOST_CHECK_EQUAL(scratch.mCellLabels.size(), nCe);

  // Group D: never plan-sized, unaffected by adoptPlan().
  BOOST_CHECK_EQUAL(scratch.mLines.size(), 0u);
  BOOST_CHECK_EQUAL(scratch.mTrackletClusters.size(), 0u);
}

BOOST_AUTO_TEST_CASE(ResetClearsWorkingStateWithoutMutatingAPopulatedTimeFrameOrTheAdoptedPlanSize)
{
  SyntheticChain chain{4};
  SurfaceTrackingScratch scratch;
  scratch.setMemoryPool(makePool());
  scratch.adoptPlan(chain.nOwnedSurfaces(), chain.nTransitions(), chain.nCells());

  // Populate a handful of containers with observable content.
  scratch.mClusters[0].emplace_back(1.f, 2.f, 3.f, 0);
  scratch.mCells[0].emplace_back();
  scratch.mBogusClusters[0] = 7;
  BOOST_REQUIRE_EQUAL(scratch.mBogusClusters.size(), chain.nOwnedSurfaces());

  // An unrelated, populated TimeFrame: reset() takes no TimeFrame parameter
  // at all, so it structurally cannot reach it -- this proves that
  // structural guarantee holds at runtime too.
  TimeFrame frame;
  frame.setBz(5.f);
  frame.setBeamPosition(1.f, 2.f, 0.1f);

  scratch.reset();

  // Vector-of-bounded_vector (Group A/B outer) containers: outer element
  // count -- the adopted plan size -- survives reset(); only each element's
  // *contents* are cleared. Mirrors SurfaceTrackingScratch::reset()
  // exactly (it never shrinks its own NLayers-wide outer arrays either).
  BOOST_CHECK_EQUAL(scratch.mClusters.size(), chain.nOwnedSurfaces());
  BOOST_CHECK(scratch.mClusters[0].empty());
  BOOST_CHECK_EQUAL(scratch.mCells.size(), chain.nCells());
  BOOST_CHECK(scratch.mCells[0].empty());
  BOOST_CHECK_EQUAL(scratch.getNOwnedSurfaces(), chain.nOwnedSurfaces());
  BOOST_CHECK_EQUAL(scratch.getNTransitions(), chain.nTransitions());
  BOOST_CHECK_EQUAL(scratch.getNCells(), chain.nCells());

  // Flat bounded_vector (Group A) containers: fully cleared to empty, not
  // preserved at the adopted plan size -- mirrors reset()'s
  // deepVectorClear(mBogusClusters) exactly.
  BOOST_CHECK_EQUAL(scratch.mBogusClusters.size(), 0u);
  BOOST_CHECK_EQUAL(scratch.mPositionResolution.size(), 0u);

  BOOST_CHECK_EQUAL(frame.getBz(), 5.f);
  BOOST_CHECK_EQUAL(frame.getBeamX(), 1.f);
}

BOOST_AUTO_TEST_CASE(AllocatorsMatchDetectsSharedVersusDistinctPools)
{
  SyntheticChain chain{4};
  auto poolA = makePool();
  auto poolB = makePool();

  SurfaceTrackingScratch live;
  live.setMemoryPool(poolA);
  live.adoptPlan(chain.nOwnedSurfaces(), chain.nTransitions(), chain.nCells());

  SurfaceTrackingScratch stagedSamePool;
  stagedSamePool.setMemoryPool(poolA);
  stagedSamePool.adoptPlan(chain.nOwnedSurfaces(), chain.nTransitions(), chain.nCells());
  BOOST_CHECK(live.allocatorsMatch(stagedSamePool));

  SurfaceTrackingScratch stagedDifferentPool;
  stagedDifferentPool.setMemoryPool(poolB);
  stagedDifferentPool.adoptPlan(chain.nOwnedSurfaces(), chain.nTransitions(), chain.nCells());
  BOOST_CHECK(!live.allocatorsMatch(stagedDifferentPool));
}

BOOST_AUTO_TEST_CASE(SwapExchangesContentAndPreservesLiveAllocatorIdentity)
{
  SyntheticChain chain{4};
  auto poolA = makePool();

  SurfaceTrackingScratch live;
  live.setMemoryPool(poolA);
  live.adoptPlan(chain.nOwnedSurfaces(), chain.nTransitions(), chain.nCells());
  live.mBogusClusters[0] = 1;

  SurfaceTrackingScratch staged;
  staged.setMemoryPool(poolA); // same resource -- allocatorsMatch() precondition for swap().
  staged.adoptPlan(chain.nOwnedSurfaces(), chain.nTransitions(), chain.nCells());
  staged.mBogusClusters[0] = 42;
  staged.mClusters[0].emplace_back(9.f, 8.f, 7.f, 0);

  BOOST_REQUIRE(live.allocatorsMatch(staged));
  live.swap(staged);

  BOOST_CHECK_EQUAL(live.mBogusClusters[0], 42);
  BOOST_CHECK_EQUAL(live.mClusters[0].size(), 1u);
  BOOST_CHECK_EQUAL(staged.mBogusClusters[0], 1);
  BOOST_CHECK(staged.mClusters[0].empty());

  // Allocator identity is owner-bound, never staged data: both sides still
  // report the exact same pool object after swap().
  BOOST_CHECK(live.getMemoryPool() == poolA);
  BOOST_CHECK(staged.getMemoryPool() == poolA);
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

// M6d revision: SurfaceTrackingScratch is now wired into production for
// MFT specifically (SurfacePlanTrackingParticipantMFT), so it legitimately
// knows a handful of things M6c's own additive-only version could not yet:
// o2::detectors::DetID::MFT (loadNormalizedSource()'s own detector
// preflight), the literal "MFT" token (MFTNLayers,
// o2::mft::constants::mft::LayersNumber -- the M6c design note's own
// explicitly-flagged narrow exception for the auxiliary NLayers-templated
// types this milestone hardcodes, see the header's file-level doc), and
// TimeFrame.h (initialise()/loadNormalizedSource()/getPrimaryVertices()/
// updateROFVertexLookupTable() all cooperate with TimeFrame directly now --
// M6c's "never touches TimeFrame" framing applied only to that milestone's
// unwired scope).
//
// M6e2 revision: this scratch type became shared by ITS too (the combined
// and standalone ITS common-CA participants now both back onto
// SurfaceTrackingScratch, not just MFT's), so o2::detectors::DetID::ITS is
// now a legitimate mention too (loadNormalizedSource()'s preflight accepts
// both), removed from the forbidden list below. What remains genuinely
// forbidden: workflow/DPL/output-layer naming, and the detail/-confined
// TransitionPolicyTag/StateFamily policy-key types -- SurfaceTrackingScratch
// itself must still never reintroduce those, and no detector-specific
// *switch* (as opposed to a DetID::ITS/DetID::MFT preflight/comparison) may
// appear here.
BOOST_AUTO_TEST_CASE(NewHeadersPullNoITSWorkflowOutputOrTransitionPolicyTagDependency)
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  const fs::path trackingRoot = fs::path(testDirectory) / "..";

  const std::vector<fs::path> newFiles = {
    trackingRoot / "include/ITSMFTTracking/SurfaceTrackingScratch.h",
    trackingRoot / "src/SurfaceTrackingScratch.cxx"};
  for (const auto& path : newFiles) {
    BOOST_REQUIRE_MESSAGE(fs::is_regular_file(path), "cannot find " << path.string());
  }

  const std::vector<std::pair<std::string, std::string>> forbidden = {
    {"TransitionPolicyTag", "detail/-confined policy-key type"},
    {"StateFamily", "Barrel/Forward state-family selector"},
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
    "SurfacePlanBinding.h",
    "DetectorTraversalBinding.h",
    "MultiSourceTimeFrameLoader.h",
    "LegacyTrackerScratch.h",
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
  BOOST_CHECK(fs::is_regular_file(trackingRoot / "include/ITSMFTTracking/SurfaceTrackingScratch.h"));
  BOOST_CHECK(fs::is_regular_file(trackingRoot / "include/ITSMFTTracking/detail/SurfacePlanBinding.h"));
}

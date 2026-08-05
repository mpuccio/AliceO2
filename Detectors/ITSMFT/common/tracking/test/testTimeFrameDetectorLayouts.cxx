// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 B2 Slice 2 rewrote this file. Its earlier content tested
// TimeFrame::ensureDetectorLayouts()/invalidateDetectorLayouts() and the
// runtime-provider-backed catalog build/validate/currency/geometry-epoch
// machinery -- all of it removed in this slice: TimeFrame owns no
// catalog/layout/plan/epoch at all any more, and buildDetectorLayoutSet()
// (DetectorLayoutSet.h) performs no runtime catalog validation (the static
// catalog it borrows in production is already proven valid at compile time,
// see SurfaceSpec.h). Every test that only existed to exercise that removed
// machinery is gone; the tests that exercise TrackerTraits::
// initialiseTimeFrame()'s own validation logic (legacy parity, material
// compatibility, mixed-policy/invalid-schedule rejection, iteration bounds)
// are kept, migrated to build a local DetectorLayoutSet via
// buildDetectorLayoutSet() and pass it to initialiseTimeFrame() as its
// explicit plan parameter, exactly as ITSMFTTrackingInterface does in
// production. One new focused test (buildDetectorLayoutSetRejects...) covers
// buildDetectorLayoutSet()'s own two failure modes, which had no other
// coverage after the deleted tests.

#define BOOST_TEST_MODULE ITSMFT TimeFrame detector layouts
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <TGeoGlobalMagField.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

#include "Field/MagneticField.h"
#include "ITSMFTTracking/DetectorLayoutSet.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/NominalSurfaceMaterialDefaults.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITStracking/Constants.h"

using namespace o2::itsmft::tracking;
using o2::itsmft::TrackingParameters;

struct TraversalPropagatorFieldFixture {
  TraversalPropagatorFieldFixture()
  {
    if (!TGeoGlobalMagField::Instance()->GetField()) {
      TGeoGlobalMagField::Instance()->SetField(o2::field::MagneticField::createNominalField(5, true));
      TGeoGlobalMagField::Instance()->Lock();
    }
  }
};

BOOST_GLOBAL_FIXTURE(TraversalPropagatorFieldFixture);

namespace
{

/// DetectorLayout no longer owns a surface copy (Slice 3, shared ownership):
/// test fixtures that build one in isolation keep the surfaces alongside it.
struct BuiltLayout {
  DetectorLayout layout;
  std::vector<SurfaceDescriptor> surfaces;
};

static_assert(std::is_nothrow_move_constructible_v<DetectorLayoutSet>);

// Nominal material matching this detector's default TrackingParameters::LayerxX0
// (o2::itsmft::resetDetectorDefaults()), so fixtures that reach
// TrackerTraits::initialiseTimeFrame() with unperturbed parameters satisfy the
// LegacyMaterialMismatch compatibility check by construction. Indexed by
// global surface id, which equals the legacy layer index for every identity-
// ordered fixture in this file.
float nominalXOverX0(o2::detectors::DetID::ID detector, uint16_t surfaceIndex)
{
  if (detector == o2::detectors::DetID::MFT) {
    return kNominalMFTLayerX0[surfaceIndex % MFTNLayers];
  }
  return kNominalITSLayerX0[surfaceIndex % ITSNLayers];
}

std::vector<SurfaceDescriptor> catalog(size_t count, SurfaceKind kind = SurfaceKind::Cylinder,
                                       o2::detectors::DetID::ID detector = o2::detectors::DetID::ITS)
{
  std::vector<SurfaceDescriptor> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(SurfaceDescriptor{SurfaceId{i}, i, static_cast<uint8_t>(detector), kind, 0, static_cast<float>(i + 1), 0.f, 100.f});
    const float xOverX0 = nominalXOverX0(detector, i);
    result.back().material.xOverX0 = xOverX0;
    result.back().material.arealDensityGPerCm2 = xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho;
  }
  return result;
}

std::vector<SurfaceId> order(size_t count)
{
  std::vector<SurfaceId> result;
  for (uint16_t i = 0; i < count; ++i) {
    result.emplace_back(i);
  }
  return result;
}

// Owns both the catalog and the DetectorLayoutSet built from it: the plan
// borrows a SurfaceCatalogView into `catalog` (Gate 4 B2 Slice 2), so both
// must be moved together and the vector's underlying buffer address (which
// std::vector's move never changes) stays valid for the plan's lifetime.
struct BuiltPlan {
  std::vector<SurfaceDescriptor> catalog;
  DetectorLayoutSet plan;
};

BuiltPlan buildPlan(std::vector<SurfaceDescriptor> surfaces, gsl::span<const SurfaceId> ordered,
                    TransitionPolicyTag tag, gsl::span<const TrackingParameters> params)
{
  const SurfaceCatalogView view{surfaces.data(), static_cast<uint32_t>(surfaces.size())};
  auto result = buildDetectorLayoutSet(view, ordered, params);
  BOOST_REQUIRE(result.ok());
  return BuiltPlan{std::move(surfaces), std::move(*result.layout)};
}

// Wraps an already-built DetectorLayout (e.g. a deliberately cyclic or
// mixed-policy one that buildDetectorLayoutSet() itself would never produce)
// into a one-iteration DetectorLayoutSet, so initialiseTimeFrame()'s own
// fail-closed checks can be exercised against it directly -- no TimeFrame-
// subclass injection needed, since the plan is an explicit parameter now.
BuiltPlan wrapLayout(BuiltLayout built)
{
  DetectorLayoutConfigurationKey key;
  key.orderedSurfaces = order(built.surfaces.size());
  std::vector<DetectorLayout> layouts;
  layouts.push_back(std::move(built.layout));
  const SurfaceCatalogView view{built.surfaces.data(), static_cast<uint32_t>(built.surfaces.size())};
  return BuiltPlan{std::move(built.surfaces), DetectorLayoutSet{std::move(key), view, std::move(layouts)}};
}

BuiltLayout cyclicDiskLayout()
{
  SparseTrackingTopology topology{10};
  BOOST_REQUIRE(topology.addTransition(SurfaceTransition{SurfaceId{0}, SurfaceId{1}, {}, 0}).isValid());
  BOOST_REQUIRE(topology.addTransition(SurfaceTransition{SurfaceId{1}, SurfaceId{2}, {}, 0}).isValid());
  BOOST_REQUIRE(topology.addTransition(SurfaceTransition{SurfaceId{2}, SurfaceId{0}, {}, 0}).isValid());
  BOOST_REQUIRE(topology.finalize());
  auto surfaces = catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  return BuiltLayout{DetectorLayout{surfaces, std::move(topology)}, std::move(surfaces)};
}

BuiltLayout mixedDisconnectedLayout()
{
  SparseTrackingTopology topology{10};
  for (uint16_t id = 0; id < 4; ++id) {
    BOOST_REQUIRE(topology.addTransition(SurfaceTransition{SurfaceId{id}, SurfaceId{static_cast<uint16_t>(id + 1)}, {}, 0}).isValid());
  }
  for (uint16_t id = 5; id < 9; ++id) {
    BOOST_REQUIRE(topology.addTransition(SurfaceTransition{SurfaceId{id}, SurfaceId{static_cast<uint16_t>(id + 1)}, {}, 0}).isValid());
  }
  BOOST_REQUIRE(topology.finalize());
  auto surfaces = catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  for (uint16_t id = 0; id < 5; ++id) {
    surfaces[id].kind = SurfaceKind::Cylinder;
  }
  return BuiltLayout{DetectorLayout{surfaces, std::move(topology)}, std::move(surfaces)};
}

TrackingParameters parameters(int activeCount, int maxHoles = 0, uint16_t holes = 0, uint16_t starts = 0xffff)
{
  TrackingParameters result;
  result.NLayers = activeCount;
  result.MaxHoles = maxHoles;
  result.HoleLayerMask = holes;
  result.StartLayerMask = starts;
  return result;
}

template <int NLayers>
void prepareTraversalFrame(TimeFrame& frame,
                           SurfaceTrackingScratch& scratch,
                           TrackerTraits<NLayers>& traits,
                           const std::shared_ptr<BoundedMemoryResource>& pool,
                           const std::vector<TrackingParameters>& params)
{
  frame.setMemoryPool(pool);
  scratch.setMemoryPool(pool);
  scratch.adoptPlan(NLayers, 0, 0);
  for (auto& rofOffsets : scratch.mROFramesClusters) {
    rofOffsets.resize(1, 0);
  }
  scratch.initTrackerTopologies<NLayers>(params);
  traits.setMemoryPool(pool);
  traits.adoptScratch(&scratch);
  traits.adoptFrame(&frame);
  traits.updateTrackingParameters(params);
}

std::vector<TrackingParameters> mftTraversalParameters()
{
  std::vector<TrackingParameters> params(1);
  o2::itsmft::resetDetectorDefaults(params.front(), o2::detectors::DetID::MFT);
  return params;
}

std::vector<TrackingParameters> itsTraversalParameters()
{
  std::vector<TrackingParameters> params(1);
  o2::itsmft::resetDetectorDefaults(params.front(), o2::detectors::DetID::ITS);
  return params;
}

template <int NLayers>
void checkLegacyParity(SurfaceKind kind, TransitionPolicyTag policyTag, uint16_t startMask)
{
  const auto detector = kind == SurfaceKind::Disk ? o2::detectors::DetID::MFT : o2::detectors::DetID::ITS;
  auto ordered = order(NLayers);
  std::vector<TrackingParameters> params{parameters(NLayers, 1, uint16_t{1} << (NLayers / 2), startMask)};
  auto built = buildPlan(catalog(NLayers, kind, detector), ordered, policyTag, params);

  TrackingTopology<NLayers> legacy;
  legacy.init(NLayers, params[0].MaxHoles, params[0].HoleLayerMask);
  const auto legacyView = legacy.getView();
  const auto layoutView = built.plan.getLayoutView(0);
  const auto sparse = layoutView.topology;
  BOOST_REQUIRE_EQUAL(sparse.nTransitions, legacyView.nTransitions);
  BOOST_REQUIRE_EQUAL(sparse.nCells, legacyView.nCells);
  for (uint32_t i = 0; i < sparse.nTransitions; ++i) {
    BOOST_CHECK_EQUAL(sparse.transitions[i].from.value(), legacyView.transitions[i].fromLayer);
    BOOST_CHECK_EQUAL(sparse.transitions[i].to.value(), legacyView.transitions[i].toLayer);
  }
  // Independent legacy oracle for which CellTopologyIds the *former*
  // findRoadsForPolicy predicate (hitLayerMask.last() + StartLayerMask.has())
  // would have selected as road starts, recomputed here from the frozen
  // legacy TrackingTopology<NLayers> view -- not by calling into any
  // TransitionPolicyGrouping/TrackerTraits production code.
  std::vector<CellTopologyId> legacyOracleStarts;
  for (uint32_t i = 0; i < sparse.nCells; ++i) {
    const auto& sparseCell = sparse.cells[i];
    const auto& legacyCell = legacyView.cells[i];
    BOOST_CHECK_EQUAL(sparseCell.firstTransition.value(), legacyCell.firstTransition);
    BOOST_CHECK_EQUAL(sparseCell.secondTransition.value(), legacyCell.secondTransition);
    BOOST_CHECK_EQUAL(sparseCell.hitSurfaces.value(), static_cast<uint16_t>(legacyCell.hitLayerMask));
    const auto legacyStartLayer = legacyCell.hitLayerMask.last();
    const auto sparseStartSurface = SurfaceId{static_cast<uint16_t>(sparseCell.hitSurfaces.last())};
    BOOST_CHECK_EQUAL(sparse.seedingSurfaces.has(sparseStartSurface), params[0].StartLayerMask.has(legacyStartLayer));
    if (params[0].StartLayerMask.has(legacyStartLayer)) {
      legacyOracleStarts.push_back(CellTopologyId{static_cast<uint16_t>(i)});
    }
  }

  // Item 7: exact identity-layout parity between roadStartCellsForTag() and
  // the former StartLayerMask/hitLayerMask.last() selection, for both
  // TrackerTraits<7> (ITS-like) and TrackerTraits<10> (MFT-like) call sites
  // of this helper.
  TransitionPolicyGrouping grouping{layoutView};
  BOOST_REQUIRE(grouping.valid());
  const auto starts = grouping.roadStartCellsForTag(policyTag);
  BOOST_CHECK(std::is_sorted(starts.begin(), starts.end()));
  BOOST_REQUIRE_EQUAL(starts.size(), legacyOracleStarts.size());
  for (size_t i = 0; i < starts.size(); ++i) {
    BOOST_CHECK(starts[i] == legacyOracleStarts[i]);
  }
}
} // namespace

// buildDetectorLayoutSet()'s own two failure modes (InvalidActiveCount,
// LayoutBuilderFailure) -- the only coverage this specific function had came
// from tests that also exercised the now-deleted runtime-provider machinery;
// this replaces that lost coverage narrowly.
BOOST_AUTO_TEST_CASE(buildDetectorLayoutSetRejectsInvalidActiveCountAndLayoutBuilderFailure)
{
  {
    // NLayers exceeds orderedSurfaces.size().
    const auto surfaces = catalog(7);
    const SurfaceCatalogView view{surfaces.data(), static_cast<uint32_t>(surfaces.size())};
    const auto ordered = order(7);
    const std::vector<TrackingParameters> params{parameters(8)};
    const auto result = buildDetectorLayoutSet(view, ordered, params);
    BOOST_CHECK(result.error == DetectorLayoutSetBuildError::InvalidActiveCount);
    BOOST_CHECK_EQUAL(result.failedIteration, 0u);
    BOOST_CHECK(!result.ok());
  }
  {
    // A structurally invalid subgraph (negative maxHoles) surfaces as
    // DetectorLayoutBuilder's own LayoutBuilderFailure, propagated through
    // buildDetectorLayoutSet() unchanged.
    const auto surfaces = catalog(7);
    const SurfaceCatalogView view{surfaces.data(), static_cast<uint32_t>(surfaces.size())};
    const auto ordered = order(7);
    const std::vector<TrackingParameters> params{parameters(7), parameters(6), parameters(5, -1)};
    const auto result = buildDetectorLayoutSet(view, ordered, params);
    BOOST_CHECK(result.error == DetectorLayoutSetBuildError::LayoutBuilderFailure);
    BOOST_CHECK_EQUAL(result.failedIteration, 2u);
    BOOST_CHECK(result.layoutBuildError == DetectorLayoutBuildError::NegativeMaxHoles);
  }
}

BOOST_AUTO_TEST_CASE(catalog_identity_active_count_and_mask_mapping)
{
  const std::vector<SurfaceId> ordered{SurfaceId{3}, SurfaceId{0}, SurfaceId{6}, SurfaceId{2}, SurfaceId{5}, SurfaceId{1}, SurfaceId{4}};
  std::vector<TrackingParameters> params{
    parameters(7, 1, uint16_t{1} << 1, (uint16_t{1} << 0) | (uint16_t{1} << 4)),
    parameters(4, 1, uint16_t{1} << 1, (uint16_t{1} << 0) | (uint16_t{1} << 4) | (uint16_t{1} << 6))};
  auto built = buildPlan(catalog(7), ordered, TransitionPolicyTag::CylinderCylinder, params);
  const auto* full = built.plan.getLayout(0);
  const auto* reduced = built.plan.getLayout(1);
  BOOST_REQUIRE(full != nullptr && reduced != nullptr);
  // Slice 3: the surface catalog is now owned exactly once (by this test's
  // BuiltPlan, borrowed by DetectorLayoutSet as a SurfaceCatalogView) --
  // there is no longer a separate per-DetectorLayout copy to compare
  // element-by-element. Both iterations trivially observe the same single
  // catalog by construction.
  const auto sharedCatalog = built.plan.getSurfaceCatalog();
  BOOST_REQUIRE_EQUAL(sharedCatalog.nSurfaces, 7u);
  BOOST_CHECK_LT(reduced->getTopology().getTransitions().size(), full->getTopology().getTransitions().size());
  BOOST_CHECK(full->getTopology().getView().seedingSurfaces.has(SurfaceId{3}));
  BOOST_CHECK(full->getTopology().getView().seedingSurfaces.has(SurfaceId{5}));
  BOOST_CHECK(!full->getTopology().getView().seedingSurfaces.has(SurfaceId{0}));
  BOOST_CHECK(reduced->getTopology().getView().seedingSurfaces.has(SurfaceId{3}));
  BOOST_CHECK_EQUAL(reduced->getTopology().getView().seedingSurfaces.count(), 1);
  const auto& reducedTransitions = reduced->getTopology().getTransitions();
  BOOST_REQUIRE(!reducedTransitions.empty());
  BOOST_CHECK(reducedTransitions.front().from == SurfaceId{3});
  BOOST_CHECK(reducedTransitions.front().to == SurfaceId{0});
  const auto skipped = std::find_if(reducedTransitions.begin(), reducedTransitions.end(), [](const auto& transition) {
    return transition.from == SurfaceId{3} && transition.to == SurfaceId{6};
  });
  BOOST_REQUIRE(skipped != reducedTransitions.end());
  BOOST_CHECK(skipped->skippedSurfaces.has(SurfaceId{0}));

  // Selected-iteration isolation (item 6): iterations 0 and 1 use different
  // StartLayerMask-derived seeding masks over the same catalog/ordering
  // (`ordered`: position 0 = SurfaceId 3, position 4 = SurfaceId 5).
  // Position 0 is a seeding surface in *both* iterations' masks but can
  // never be a cell's transition endpoint (transitions only run from an
  // earlier to a later position) -- "a seeded surface with no terminating
  // cell" in both layouts simultaneously. Iteration 0's other seeding
  // surface (5, position 4) is reachable within the full 7-active-surface
  // layout; iteration 1's activeCount=4 mask only ever tests position 0
  // (positions 4 and 6 of its starts bitset are out of range and ignored by
  // positionalSurfaceMask), so iteration 1's own layout must select no road
  // starts at all -- proving each iteration's roadStartCellsForTag reflects
  // only its own layout's seedingSurfaces.
  const auto fullView = built.plan.getLayoutView(0);
  const auto reducedView = built.plan.getLayoutView(1);
  TransitionPolicyGrouping fullGrouping{fullView};
  TransitionPolicyGrouping reducedGrouping{reducedView};
  BOOST_REQUIRE(fullGrouping.valid());
  BOOST_REQUIRE(reducedGrouping.valid());
  const auto fullStarts = fullGrouping.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder);
  const auto reducedStarts = reducedGrouping.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder);
  BOOST_CHECK(std::is_sorted(fullStarts.begin(), fullStarts.end()));
  BOOST_CHECK(reducedStarts.empty());
  BOOST_REQUIRE_GT(fullStarts.size(), 0u);
  for (const auto id : fullStarts) {
    const auto endpoint = fullView.topology.getTransition(fullView.topology.getCell(id).secondTransition).to;
    BOOST_CHECK(endpoint == SurfaceId{5});
  }
}

// Gate 4 B2 Slice 2: initialiseTimeFrame() takes the plan as an explicit
// `const DetectorLayoutSet&` parameter, so "missing"/"stale" plan states are
// no longer constructible at all -- only IterationOutOfRange, from the
// removed traversal_initialisation_classifies_missing_and_stale_layouts
// test, still applies and is kept here under its own name.
BOOST_AUTO_TEST_CASE(traversal_initialisation_rejects_iteration_beyond_configured_layout_set)
{
  auto twoIterations = mftTraversalParameters();
  twoIterations.push_back(twoIterations.front());
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame shortLayoutFrame;
  SurfaceTrackingScratch shortLayoutScratch;
  TrackerTraits<10> shortLayoutTraits;
  prepareTraversalFrame(shortLayoutFrame, shortLayoutScratch, shortLayoutTraits, pool, twoIterations);
  std::vector<TrackingParameters> oneLayout{twoIterations.front()};
  auto built = buildPlan(catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), TransitionPolicyTag::DiskDisk, oneLayout);
  try {
    shortLayoutTraits.initialiseTimeFrame(1, built.plan);
    BOOST_FAIL("iteration beyond the configured layout set must throw");
  } catch (const TraversalException& error) {
    BOOST_CHECK_EQUAL(error.getIteration(), 1);
    BOOST_CHECK(error.getReason() == TraversalFailureReason::IterationOutOfRange);
  }
  BOOST_CHECK(!shortLayoutTraits.hasTraversalCache());
}

BOOST_AUTO_TEST_CASE(traversal_cache_groups_once_across_repeated_neighbour_and_road_calls)
{
  // M4b removed TrackerTraits::getPolicyBindingCount(StateFamily) (a
  // test-only introspection seam using StateFamily as a policy-binding-key
  // proxy); getTraversalGroupingCount() staying at 1 across every repeated
  // call below remains the public-API evidence that initialiseTimeFrame()'s
  // traversal grouping/policy-parameter binding is not redone per call.
  auto params = mftTraversalParameters();
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  TrackerTraits<10> traits;
  prepareTraversalFrame(frame, scratch, traits, pool, params);
  auto built = buildPlan(catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), TransitionPolicyTag::DiskDisk, params);

  std::shared_ptr<tbb::task_arena> arena;
  traits.setNThreads(1, arena);
  traits.initialiseTimeFrame(0, built.plan);
  BOOST_REQUIRE(traits.hasTraversalCache());
  BOOST_CHECK_EQUAL(traits.getTraversalGroupingCount(), 1);

  traits.findCellsNeighbours(0);
  traits.findCellsNeighbours(0);
  traits.findRoads(0);
  traits.findRoads(0);
  BOOST_CHECK_EQUAL(traits.getTraversalGroupingCount(), 1);

  std::vector<TrackingParameters> itsParams{parameters(7, 0, 0, 0x7f)};
  TimeFrame itsFrame;
  SurfaceTrackingScratch itsScratch;
  TrackerTraits<7> itsTraits;
  prepareTraversalFrame(itsFrame, itsScratch, itsTraits, pool, itsParams);
  auto itsBuilt = buildPlan(catalog(7, SurfaceKind::Cylinder, o2::detectors::DetID::ITS), order(7), TransitionPolicyTag::CylinderCylinder, itsParams);
  itsTraits.setNThreads(1, arena);
  itsTraits.initialiseTimeFrame(0, itsBuilt.plan);
  itsTraits.findRoads(0);
  itsTraits.findRoads(0);
  BOOST_CHECK_EQUAL(itsTraits.getTraversalGroupingCount(), 1);
}

BOOST_AUTO_TEST_CASE(traversal_empty_road_start_span_is_valid_and_produces_no_tracks)
{
  // StartLayerMask=0 -> an empty seeding mask -> an empty roadStartCellsForTag
  // span (Architecture.md Sec 10, item 7: "empty road-start span is valid").
  // Unlike testCATrackerFailureContract.cxx's
  // ValidEmptyInputCompletesWithoutErrorAndProducesNoTracks (full StartLayerMask,
  // no cluster data), this exercises the *topologically empty* road-start case
  // through the same initialiseTimeFrame()/findRoads() pair.
  auto params = mftTraversalParameters();
  params[0].StartLayerMask = 0;
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  TrackerTraits<10> traits;
  prepareTraversalFrame(frame, scratch, traits, pool, params);
  auto built = buildPlan(catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), TransitionPolicyTag::DiskDisk, params);

  std::shared_ptr<tbb::task_arena> arena;
  traits.setNThreads(1, arena);
  BOOST_CHECK_NO_THROW(traits.initialiseTimeFrame(0, built.plan));
  BOOST_REQUIRE(traits.hasTraversalCache());
  BOOST_CHECK_NO_THROW(traits.findRoads(0));
  BOOST_CHECK_EQUAL(frame.getCommonTracks().size(), 0u);
}

BOOST_AUTO_TEST_CASE(traversal_legacy_cell_container_size_mismatch_fails_before_indexing)
{
  // Item 4/7: findRoads() indexes mScratch->getCells() with sparse
  // CellTopologyId values; a desync between that legacy container and the
  // cached sparse layout must fail with LegacyIndexMismatch rather than
  // index out of bounds. Reached here through
  // LegacyTrackerScratch::getCells(), the existing public, non-const
  // production scratch accessor -- no new mutation API is
  // added for this test.
  auto params = mftTraversalParameters();
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  TrackerTraits<10> traits;
  prepareTraversalFrame(frame, scratch, traits, pool, params);
  auto built = buildPlan(catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), TransitionPolicyTag::DiskDisk, params);

  std::shared_ptr<tbb::task_arena> arena;
  traits.setNThreads(1, arena);
  traits.initialiseTimeFrame(0, built.plan);
  BOOST_REQUIRE(traits.hasTraversalCache());
  BOOST_REQUIRE(!scratch.getCells().empty());
  scratch.getCells().pop_back();

  try {
    traits.findRoads(0);
    BOOST_FAIL("legacy cell-container size mismatch must throw before indexing");
  } catch (const TraversalException& error) {
    BOOST_CHECK(error.getReason() == TraversalFailureReason::LegacyIndexMismatch);
  }
}

BOOST_AUTO_TEST_CASE(traversal_preflight_rejects_legacy_mismatch_state_mismatch_and_bad_parameters)
{
  auto checkFailure = [](std::vector<TrackingParameters> params,
                         std::vector<SurfaceDescriptor> surfaces,
                         std::vector<SurfaceId> ordered,
                         TransitionPolicyTag tag,
                         TraversalFailureReason expected) {
    auto pool = std::make_shared<BoundedMemoryResource>();
    TimeFrame frame;
    SurfaceTrackingScratch scratch;
    TrackerTraits<10> traits;
    prepareTraversalFrame(frame, scratch, traits, pool, params);
    auto built = buildPlan(std::move(surfaces), ordered, tag, params);
    try {
      traits.initialiseTimeFrame(0, built.plan);
      BOOST_FAIL("invalid traversal preflight must throw");
    } catch (const TraversalException& error) {
      BOOST_CHECK(error.getReason() == expected);
    }
    BOOST_CHECK(!traits.hasTraversalCache());
  };

  auto params = mftTraversalParameters();
  checkFailure(params, catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT),
               {SurfaceId{3}, SurfaceId{0}, SurfaceId{6}, SurfaceId{2}, SurfaceId{9}, SurfaceId{5}, SurfaceId{1}, SurfaceId{8}, SurfaceId{4}, SurfaceId{7}},
               TransitionPolicyTag::DiskDisk, TraversalFailureReason::LegacyIndexMismatch);
  checkFailure(params, catalog(10, SurfaceKind::Cylinder, o2::detectors::DetID::MFT), order(10),
               TransitionPolicyTag::CylinderCylinder, TraversalFailureReason::StateFamilyMismatch);
  params[0].MaxChi2ClusterAttachment = -1.f;
  checkFailure(params, catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), TransitionPolicyTag::DiskDisk,
               TraversalFailureReason::InvalidPolicyParameters);

  // Perturbing the temporary legacy LayerxX0 away from the catalog's
  // authoritative material now fails the material-compatibility check
  // (LegacyMaterialMismatch) before ever reaching attachHitConfig's own
  // finite/non-negative validation -- the catalog's material is unperturbed
  // and finite, so the mismatch is purely numeric.
  params = mftTraversalParameters();
  params[0].LayerxX0[3] = std::numeric_limits<float>::quiet_NaN();
  checkFailure(params, catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), TransitionPolicyTag::DiskDisk,
               TraversalFailureReason::LegacyMaterialMismatch);

  params = mftTraversalParameters();
  params[0].LayerxX0[7] = std::numeric_limits<float>::infinity();
  checkFailure(params, catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), TransitionPolicyTag::DiskDisk,
               TraversalFailureReason::LegacyMaterialMismatch);

  params = mftTraversalParameters();
  params[0].CorrType = static_cast<o2::base::PropagatorF::MatCorrType>(99);
  checkFailure(params, catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), TransitionPolicyTag::DiskDisk,
               TraversalFailureReason::InvalidPolicyParameters);
}

BOOST_AUTO_TEST_CASE(every_iteration_resolves_identical_authoritative_material)
{
  auto params = mftTraversalParameters();
  params.push_back(params.front());
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  TrackerTraits<10> traits;
  prepareTraversalFrame(frame, scratch, traits, pool, params);
  auto built = buildPlan(catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), TransitionPolicyTag::DiskDisk, params);

  traits.initialiseTimeFrame(0, built.plan);
  const auto firstIterationMaterial = traits.getLayerMaterial();
  std::vector<NominalSurfaceMaterial> firstIteration(firstIterationMaterial.begin(), firstIterationMaterial.end());

  traits.initialiseTimeFrame(1, built.plan);
  const auto secondIteration = traits.getLayerMaterial();
  BOOST_REQUIRE_EQUAL(secondIteration.size(), firstIteration.size());
  for (size_t layer = 0; layer < firstIteration.size(); ++layer) {
    BOOST_CHECK_EQUAL(secondIteration[layer].xOverX0, firstIteration[layer].xOverX0);
    BOOST_CHECK_EQUAL(secondIteration[layer].arealDensityGPerCm2, firstIteration[layer].arealDensityGPerCm2);
  }
}

BOOST_AUTO_TEST_CASE(rejected_initialisation_does_not_mutate_surface_descriptor_material)
{
  auto params = mftTraversalParameters();
  params[0].LayerxX0[4] = std::numeric_limits<float>::quiet_NaN();
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  TrackerTraits<10> traits;
  prepareTraversalFrame(frame, scratch, traits, pool, params);
  auto built = buildPlan(catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), TransitionPolicyTag::DiskDisk, params);
  const NominalSurfaceMaterial materialBefore = built.plan.getSurfaceCatalog().getSurface(SurfaceId{4}).material;

  try {
    traits.initialiseTimeFrame(0, built.plan);
    BOOST_FAIL("perturbed LayerxX0 must throw");
  } catch (const TraversalException& error) {
    BOOST_CHECK(error.getReason() == TraversalFailureReason::LegacyMaterialMismatch);
  }

  const auto& materialAfter = built.plan.getSurfaceCatalog().getSurface(SurfaceId{4}).material;
  BOOST_CHECK_EQUAL(materialAfter.xOverX0, materialBefore.xOverX0);
  BOOST_CHECK_EQUAL(materialAfter.arealDensityGPerCm2, materialBefore.arealDensityGPerCm2);
}

BOOST_AUTO_TEST_CASE(non_monotonic_ordered_surfaces_maps_correctly_then_resets_on_later_failure)
{
  // Distinct-per-surface material (not the uniform MFT nominal default) so an
  // identity-assuming mapping bug (mLayerMaterial[legacyLayer] read from
  // catalog[legacyLayer] instead of catalog[orderedSurfaces[legacyLayer]])
  // would be observable as a numeric difference, not masked by every layer
  // happening to carry the same value.
  auto surfaces = catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  for (uint16_t id = 0; id < surfaces.size(); ++id) {
    const float xOverX0 = 0.001f * static_cast<float>(id + 1);
    surfaces[id].material.xOverX0 = xOverX0;
    surfaces[id].material.arealDensityGPerCm2 = xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho;
  }
  const std::vector<SurfaceId> nonMonotonicOrder{
    SurfaceId{3}, SurfaceId{0}, SurfaceId{6}, SurfaceId{2}, SurfaceId{9},
    SurfaceId{5}, SurfaceId{1}, SurfaceId{8}, SurfaceId{4}, SurfaceId{7}};

  auto params = mftTraversalParameters();
  for (size_t legacyLayer = 0; legacyLayer < nonMonotonicOrder.size(); ++legacyLayer) {
    params[0].LayerxX0[legacyLayer] = surfaces[nonMonotonicOrder[legacyLayer].value()].material.xOverX0;
  }

  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  TrackerTraits<10> traits;
  prepareTraversalFrame(frame, scratch, traits, pool, params);
  auto built = buildPlan(std::move(surfaces), nonMonotonicOrder, TransitionPolicyTag::DiskDisk, params);

  // The non-identity mapping is structurally incompatible with the separate
  // legacy-topology-parity check (validateLegacyParity), so the overall call
  // still fails -- with LegacyIndexMismatch, not a material reason. Getting
  // LegacyIndexMismatch here (rather than LegacyMaterialMismatch) is itself
  // the proof that the material-compatibility check used the correct
  // (non-identity) orderedSurfaces mapping: params.LayerxX0 above was built
  // from surfaces[nonMonotonicOrder[legacyLayer]], so an identity-mapping bug
  // (reading surfaces[legacyLayer] instead) would disagree with it at every
  // permuted position and fail earlier with LegacyMaterialMismatch instead.
  try {
    traits.initialiseTimeFrame(0, built.plan);
    BOOST_FAIL("non-monotonic ordering must fail legacy topology parity");
  } catch (const TraversalException& error) {
    BOOST_CHECK(error.getReason() == TraversalFailureReason::LegacyIndexMismatch);
  }

  // Material validation passing is not the same as initialisation succeeding:
  // this later (unrelated) failure must still leave mLayerMaterial exactly at
  // its resetTraversalCache() state, never partially populated from the
  // staged-but-never-committed resolution above. hasTraversalCache() is the
  // existing single source of truth that the call did not succeed, and
  // mAttachHitConfig (material-dependent, committed in the same final block
  // as mTraversalGrouping) is therefore equally not in effect.
  BOOST_CHECK(!traits.hasTraversalCache());
  const auto resolvedMaterial = traits.getLayerMaterial();
  BOOST_REQUIRE_EQUAL(resolvedMaterial.size(), nonMonotonicOrder.size());
  for (const auto& material : resolvedMaterial) {
    BOOST_CHECK_EQUAL(material.xOverX0, 0.f);
    BOOST_CHECK_EQUAL(material.arealDensityGPerCm2, 0.f);
  }
}

BOOST_AUTO_TEST_CASE(traversal_preflight_reports_invalid_schedule_and_mixed_policy_layout)
{
  auto checkInstalledLayout = [](BuiltLayout layout, TraversalFailureReason expected) {
    auto params = mftTraversalParameters();
    auto pool = std::make_shared<BoundedMemoryResource>();
    TimeFrame frame;
    SurfaceTrackingScratch scratch;
    TrackerTraits<10> traits;
    prepareTraversalFrame(frame, scratch, traits, pool, params);
    auto built = wrapLayout(std::move(layout));
    try {
      traits.initialiseTimeFrame(0, built.plan);
      BOOST_FAIL("invalid installed layout must throw");
    } catch (const TraversalException& error) {
      BOOST_CHECK(error.getReason() == expected);
    }
    BOOST_CHECK(!traits.hasTraversalCache());
  };

  checkInstalledLayout(cyclicDiskLayout(), TraversalFailureReason::InvalidTraversalSchedule);
  checkInstalledLayout(mixedDisconnectedLayout(), TraversalFailureReason::MixedPolicyLayout);
}

BOOST_AUTO_TEST_CASE(its_legacy_topology_and_road_start_parity)
{
  checkLegacyParity<7>(SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder,
                       (uint16_t{1} << 6) | (uint16_t{1} << 3));
}

BOOST_AUTO_TEST_CASE(mft_legacy_topology_and_road_start_parity)
{
  checkLegacyParity<10>(SurfaceKind::Disk, TransitionPolicyTag::DiskDisk,
                        (uint16_t{1} << 9) | (uint16_t{1} << 5));
}

// Gate 4 M5b: TrackerTraits<NLayers>'s dispatchActivePolicy()/computeLayerTracklets()
// et al. no longer gate which TransitionPolicyTraits<Tag> orchestration body
// is even compiled on NLayers -- both are now compiled for every NLayers (see
// TrackerTraits.cxx's dispatchActivePolicy() doc). This is a compile-time
// instantiation-footprint change only: the runtime StateFamilyMismatch gate
// in initialiseTimeFrame() is unchanged, so an ITS-shaped (NLayers=7)
// instantiation fed a Disk-kind catalog must still fail exactly like an
// MFT-shaped one fed a Cylinder-kind catalog does above
// (traversal_preflight_rejects_legacy_mismatch_state_mismatch_and_bad_parameters).
// This is the missing symmetric direction of that existing coverage.
BOOST_AUTO_TEST_CASE(its_family_rejects_disk_kind_catalog_after_m5b_dispatch_removal)
{
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  TrackerTraits<7> traits;
  auto params = itsTraversalParameters();
  prepareTraversalFrame(frame, scratch, traits, pool, params);
  auto built = buildPlan(catalog(7, SurfaceKind::Disk, o2::detectors::DetID::ITS), order(7),
                         TransitionPolicyTag::DiskDisk, params);
  try {
    traits.initialiseTimeFrame(0, built.plan);
    BOOST_FAIL("invalid traversal preflight must throw");
  } catch (const TraversalException& error) {
    BOOST_CHECK(error.getReason() == TraversalFailureReason::StateFamilyMismatch);
  }
  BOOST_CHECK(!traits.hasTraversalCache());
}

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 B2 Slice 2 rewrote this file. Its earlier content tested
// TimeFrame::ensureSurfaceGraphs()/invalidateSurfaceGraphs() and the
// runtime-provider-backed catalog build/validate/currency/geometry-epoch
// machinery -- all of it removed in this slice: TimeFrame owns no
// catalog/layout/plan/epoch at all any more, and buildSurfaceGraphs()
// (std::vector<SurfaceGraph>.h) performs no runtime catalog validation (the static
// catalog it borrows in production is already proven valid at compile time,
// see SurfaceSpec.h). Every test that only existed to exercise that removed
// machinery is gone; the tests that exercise TrackerTraits::
// initialiseTimeFrame()'s own validation logic (legacy parity, material
// compatibility, mixed-kind/invalid-schedule rejection, iteration bounds)
// are kept, migrated to build a local std::vector<SurfaceGraph> via
// buildSurfaceGraphs() and pass it to initialiseTimeFrame() as its
// explicit plan parameter, exactly as the standalone workflow does in
// production. One new focused test (buildSurfaceGraphsRejects...) covers
// buildSurfaceGraphs()'s own two failure modes, which had no other
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
#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
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

bool noopSeedRefit(const TrackSeed&,
                   const TrackingParameters&,
                   float,
                   SurfaceTrackingScratch&,
                   gsl::span<const gsl::span<const GlobalMeasurement>>,
                   gsl::span<const gsl::span<const SurfaceMeasurement>>,
                   SurfaceCatalogView,
                   gsl::span<const SurfaceId>,
                   TrackingCandidate&)
{
  return false;
}

/// SurfaceGraph no longer owns a surface copy (Slice 3, shared ownership):
/// test fixtures that build one in isolation keep the surfaces alongside it.
struct BuiltLayout {
  SurfaceGraph layout;
  std::vector<SurfaceDescriptor> surfaces;
};

static_assert(std::is_nothrow_move_constructible_v<std::vector<SurfaceGraph>>);

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

// Owns both the catalog and the std::vector<SurfaceGraph> built from it: the plan
// keeps the catalog and graph vector together as one fixture so callers can
// exercise graph construction and runtime-plan installation as one unit.
struct BuiltPlan {
  std::vector<SurfaceDescriptor> catalog;
  std::vector<SurfaceGraph> plan;
  std::vector<std::unique_ptr<SurfacePlanBinding>> bindings;
};

void adoptPlanBinding(BuiltPlan& built, TrackerTraits& traits, int iteration)
{
  const auto& graph = built.plan[iteration];
  const auto ordered = graph.getOrderedSurfaces();
  SurfaceMask owned;
  for (const auto surface : ordered) {
    owned.set(surface);
  }
  const auto kind = ordered.empty() ? SurfaceKind::Cylinder : graph.getView().getSurface(ordered.front()).kind;
  auto result = SurfacePlanBinding::build(graph.getView(), owned, ordered);
  BOOST_REQUIRE(result.ok());
  traits.adoptSurfacePlanBinding(result.binding.get());
  built.bindings.push_back(std::move(result.binding));
}

BuiltPlan buildPlan(std::vector<SurfaceDescriptor> surfaces, gsl::span<const SurfaceId> ordered,
                    SurfaceKind tag, gsl::span<const TrackingParameters> params)
{
  const SurfaceCatalogView view{surfaces.data(), static_cast<uint32_t>(surfaces.size())};
  auto result = buildSurfaceGraphs(view, ordered, params);
  BOOST_REQUIRE(result.ok());
  return BuiltPlan{std::move(surfaces), std::move(result.graphs), {}};
}

// Wraps an already-built SurfaceGraph (e.g. a deliberately cyclic or
// mixed-kind one that buildSurfaceGraphs() itself would never produce)
// into a one-iteration std::vector<SurfaceGraph>, so initialiseTimeFrame()'s own
// fail-closed checks can be exercised against it directly -- no TimeFrame-
// subclass injection needed, since the plan is an explicit parameter now.
BuiltPlan wrapLayout(BuiltLayout built)
{
  std::vector<SurfaceGraph> layouts;
  layouts.push_back(std::move(built.layout));
  return BuiltPlan{std::move(built.surfaces), std::move(layouts), {}};
}

BuiltLayout cyclicDiskLayout()
{
  SurfaceGraph topology{10};
  BOOST_REQUIRE(topology.addTransition(SurfaceTransition{SurfaceId{0}, SurfaceId{1}, {}, 0}).isValid());
  BOOST_REQUIRE(topology.addTransition(SurfaceTransition{SurfaceId{1}, SurfaceId{2}, {}, 0}).isValid());
  BOOST_REQUIRE(topology.addTransition(SurfaceTransition{SurfaceId{2}, SurfaceId{0}, {}, 0}).isValid());
  BOOST_REQUIRE(topology.finalize());
  auto surfaces = catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  return BuiltLayout{SurfaceGraph{surfaces, std::move(topology)}, std::move(surfaces)};
}

BuiltLayout mixedDisconnectedLayout()
{
  SurfaceGraph topology{10};
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
  return BuiltLayout{SurfaceGraph{surfaces, std::move(topology)}, std::move(surfaces)};
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

void prepareTraversalFrame(TimeFrame& frame,
                           SurfaceTrackingScratch& scratch,
                           TrackerTraits& traits,
                           const std::shared_ptr<BoundedMemoryResource>& pool,
                           const std::vector<TrackingParameters>& params)
{
  frame.setMemoryPool(pool);
  scratch.setMemoryPool(pool);
  // The fixtures using this helper build the identity, no-hole plan below.
  // Seed the workspace with that plan's sparse runtime extents rather than
  // leaving the transition/cell vectors at the old zero placeholder.
  const int activeCount = params.empty() ? 0 : params.front().NLayers;
  scratch.adoptPlan(activeCount, activeCount > 0 ? activeCount - 1 : 0,
                    activeCount > 1 ? activeCount - 2 : 0);
  for (auto& rofOffsets : scratch.mROFramesClusters) {
    rofOffsets.resize(1, 0);
  }
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

} // namespace

// buildSurfaceGraphs()'s own two failure modes (GraphRejected,
// GraphRejected) -- the only coverage this specific function had came
// from tests that also exercised the now-deleted runtime-provider machinery;
// this replaces that lost coverage narrowly.
BOOST_AUTO_TEST_CASE(buildSurfaceGraphsRejectsGraphRejectedAndGraphRejected)
{
  {
    // NLayers exceeds orderedSurfaces.size().
    const auto surfaces = catalog(7);
    const SurfaceCatalogView view{surfaces.data(), static_cast<uint32_t>(surfaces.size())};
    const auto ordered = order(7);
    const std::vector<TrackingParameters> params{parameters(8)};
    const auto result = buildSurfaceGraphs(view, ordered, params);
    BOOST_CHECK(result.error == SurfaceGraphBuildError::GraphRejected);
    BOOST_CHECK_EQUAL(result.failedIteration, 0u);
    BOOST_CHECK(!result.ok());
  }
  {
    // A structurally invalid subgraph (negative maxHoles) surfaces as
    // SurfaceGraphBuilder's own GraphRejected, propagated through
    // buildSurfaceGraphs() unchanged.
    const auto surfaces = catalog(7);
    const SurfaceCatalogView view{surfaces.data(), static_cast<uint32_t>(surfaces.size())};
    const auto ordered = order(7);
    const std::vector<TrackingParameters> params{parameters(7), parameters(6), parameters(5, -1)};
    const auto result = buildSurfaceGraphs(view, ordered, params);
    BOOST_CHECK(result.error == SurfaceGraphBuildError::GraphRejected);
    BOOST_CHECK_EQUAL(result.failedIteration, 2u);
    BOOST_CHECK(result.detail == SurfaceGraphBuildError::NegativeMaxHoles);
  }
}

BOOST_AUTO_TEST_CASE(catalog_identity_active_count_and_mask_mapping)
{
  const std::vector<SurfaceId> ordered{SurfaceId{3}, SurfaceId{0}, SurfaceId{6}, SurfaceId{2}, SurfaceId{5}, SurfaceId{1}, SurfaceId{4}};
  std::vector<TrackingParameters> params{
    parameters(7, 1, uint16_t{1} << 1, (uint16_t{1} << 0) | (uint16_t{1} << 4)),
    parameters(4, 1, uint16_t{1} << 1, (uint16_t{1} << 0) | (uint16_t{1} << 4) | (uint16_t{1} << 6))};
  auto built = buildPlan(catalog(7), ordered, SurfaceKind::Cylinder, params);
  BOOST_REQUIRE_EQUAL(built.plan.size(), 2u);
  const auto& full = built.plan[0];
  const auto& reduced = built.plan[1];
  // Slice 3: the surface catalog is now owned exactly once (by this test's
  // BuiltPlan, borrowed by std::vector<SurfaceGraph> as a SurfaceCatalogView) --
  // there is no longer a separate per-SurfaceGraph copy to compare
  // element-by-element. Both iterations trivially observe the same single
  // catalog by construction.
  const auto sharedCatalog = full.getSurfaceCatalog();
  BOOST_REQUIRE_EQUAL(sharedCatalog.nSurfaces, 7u);
  BOOST_CHECK_LT(reduced.getTransitions().size(), full.getTransitions().size());
  BOOST_CHECK(full.getView().seedingSurfaces.has(SurfaceId{3}));
  BOOST_CHECK(full.getView().seedingSurfaces.has(SurfaceId{5}));
  BOOST_CHECK(!full.getView().seedingSurfaces.has(SurfaceId{0}));
  BOOST_CHECK(reduced.getView().seedingSurfaces.has(SurfaceId{3}));
  BOOST_CHECK_EQUAL(reduced.getView().seedingSurfaces.count(), 1);
  const auto& reducedTransitions = reduced.getTransitions();
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
  const auto fullView = full.getView();
  const auto reducedView = reduced.getView();
  SurfaceMask fullOwned;
  for (const auto surface : full.getOrderedSurfaces()) {
    fullOwned.set(surface);
  }
  SurfaceMask reducedOwned;
  for (const auto surface : reduced.getOrderedSurfaces()) {
    reducedOwned.set(surface);
  }
  const auto fullBinding = SurfacePlanBinding::build(fullView, fullOwned, full.getOrderedSurfaces());
  const auto reducedBinding = SurfacePlanBinding::build(reducedView, reducedOwned, reduced.getOrderedSurfaces());
  BOOST_REQUIRE(fullBinding.ok());
  BOOST_REQUIRE(reducedBinding.ok());
  const auto fullStarts = fullBinding.binding->getGlobalRoadStartCells();
  const auto reducedStarts = reducedBinding.binding->getGlobalRoadStartCells();
  BOOST_CHECK(std::is_sorted(fullStarts.begin(), fullStarts.end()));
  BOOST_CHECK(reducedStarts.empty());
  BOOST_REQUIRE_GT(fullStarts.size(), 0u);
  for (const auto id : fullStarts) {
    const auto endpoint = fullView.getTransition(fullView.getCell(id).secondTransition).to;
    BOOST_CHECK(endpoint == SurfaceId{5});
  }
}

// Gate 4 B2 Slice 2: initialiseTimeFrame() takes the plan as an explicit
// `const std::vector<SurfaceGraph>&` parameter, so "missing"/"stale" plan states are
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
  TrackerTraits shortLayoutTraits;
  prepareTraversalFrame(shortLayoutFrame, shortLayoutScratch, shortLayoutTraits, pool, twoIterations);
  std::vector<TrackingParameters> oneLayout{twoIterations.front()};
  auto built = buildPlan(catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), SurfaceKind::Disk, oneLayout);
  try {
    shortLayoutTraits.initialiseTimeFrame(1, built.plan);
    BOOST_FAIL("iteration beyond the configured layout set must throw");
  } catch (const TraversalException& error) {
    BOOST_CHECK_EQUAL(error.getIteration(), 1);
    BOOST_CHECK(error.getReason() == TraversalFailureReason::IterationOutOfRange);
  }
  BOOST_CHECK(!shortLayoutTraits.hasTraversalCache());
}

BOOST_AUTO_TEST_CASE(repeated_neighbour_and_road_stages_preserve_traversal_cache)
{
  auto params = mftTraversalParameters();
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  TrackerTraits traits;
  prepareTraversalFrame(frame, scratch, traits, pool, params);
  auto built = buildPlan(catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), SurfaceKind::Disk, params);

  std::shared_ptr<tbb::task_arena> arena;
  traits.setNThreads(1, arena);
  adoptPlanBinding(built, traits, 0);
  traits.initialiseTimeFrame(0, built.plan);
  BOOST_REQUIRE(traits.hasTraversalCache());

  traits.findCellsNeighbours(0);
  traits.findCellsNeighbours(0);
  traits.findRoads(0, noopSeedRefit);
  traits.findRoads(0, noopSeedRefit);

  std::vector<TrackingParameters> itsParams{parameters(7, 0, 0, 0x7f)};
  TimeFrame itsFrame;
  SurfaceTrackingScratch itsScratch;
  TrackerTraits itsTraits;
  prepareTraversalFrame(itsFrame, itsScratch, itsTraits, pool, itsParams);
  auto itsBuilt = buildPlan(catalog(7, SurfaceKind::Cylinder, o2::detectors::DetID::ITS), order(7), SurfaceKind::Cylinder, itsParams);
  itsTraits.setNThreads(1, arena);
  adoptPlanBinding(itsBuilt, itsTraits, 0);
  itsTraits.initialiseTimeFrame(0, itsBuilt.plan);
  itsTraits.findRoads(0, noopSeedRefit);
  itsTraits.findRoads(0, noopSeedRefit);
}

BOOST_AUTO_TEST_CASE(traversal_empty_road_start_span_is_valid_and_produces_no_tracks)
{
  // StartLayerMask=0 -> an empty seeding mask -> an empty roadStartCellsForTag
  // span (Architecture.md Sec 10, item 7: "empty road-start span is valid").
  // Unlike testTrackerFailureContract.cxx's
  // ValidEmptyInputCompletesWithoutErrorAndProducesNoTracks (full StartLayerMask,
  // no cluster data), this exercises the *topologically empty* road-start case
  // through the same initialiseTimeFrame()/findRoads() pair.
  auto params = mftTraversalParameters();
  params[0].StartLayerMask = 0;
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  TrackerTraits traits;
  prepareTraversalFrame(frame, scratch, traits, pool, params);
  auto built = buildPlan(catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), SurfaceKind::Disk, params);

  std::shared_ptr<tbb::task_arena> arena;
  traits.setNThreads(1, arena);
  adoptPlanBinding(built, traits, 0);
  BOOST_CHECK_NO_THROW(traits.initialiseTimeFrame(0, built.plan));
  BOOST_REQUIRE(traits.hasTraversalCache());
  BOOST_CHECK_NO_THROW(traits.findRoads(0, noopSeedRefit));
  BOOST_CHECK_EQUAL(frame.getGenericTracks().size(), 0u);
}

BOOST_AUTO_TEST_CASE(traversal_legacy_cell_container_size_mismatch_fails_before_indexing)
{
  // Item 4/7: findRoads() indexes mScratch->getCells() with sparse
  // CellTopologyId values; a desync between that legacy container and the
  // cached sparse layout must fail with SparseTopologyMismatch rather than
  // index out of bounds. Reached here through
  // SurfaceTrackingScratch::getCells(), the existing public, non-const
  // production scratch accessor -- no new mutation API is
  // added for this test.
  auto params = mftTraversalParameters();
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  TrackerTraits traits;
  prepareTraversalFrame(frame, scratch, traits, pool, params);
  auto built = buildPlan(catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), SurfaceKind::Disk, params);

  std::shared_ptr<tbb::task_arena> arena;
  traits.setNThreads(1, arena);
  adoptPlanBinding(built, traits, 0);
  traits.initialiseTimeFrame(0, built.plan);
  BOOST_REQUIRE(traits.hasTraversalCache());
  BOOST_REQUIRE(!scratch.getCells().empty());
  scratch.getCells().pop_back();

  try {
    traits.findRoads(0, noopSeedRefit);
    BOOST_FAIL("legacy cell-container size mismatch must throw before indexing");
  } catch (const TraversalException& error) {
    BOOST_CHECK(error.getReason() == TraversalFailureReason::SparseTopologyMismatch);
  }
}

BOOST_AUTO_TEST_CASE(traversal_preflight_rejects_bad_parameters_but_not_detector_shape)
{
  auto checkFailure = [](std::vector<TrackingParameters> params,
                         std::vector<SurfaceDescriptor> surfaces,
                         std::vector<SurfaceId> ordered,
                         SurfaceKind tag,
                         TraversalFailureReason expected) {
    auto pool = std::make_shared<BoundedMemoryResource>();
    TimeFrame frame;
    SurfaceTrackingScratch scratch;
    TrackerTraits traits;
    prepareTraversalFrame(frame, scratch, traits, pool, params);
    auto built = buildPlan(std::move(surfaces), ordered, tag, params);
    adoptPlanBinding(built, traits, 0);
    try {
      traits.initialiseTimeFrame(0, built.plan);
      BOOST_FAIL("invalid traversal preflight must throw");
    } catch (const TraversalException& error) {
      BOOST_CHECK(error.getReason() == expected);
    }
    BOOST_CHECK(!traits.hasTraversalCache());
  };

  auto params = mftTraversalParameters();
  // The runtime-plan core no longer infers a state family from detector
  // identity or layer count. A valid cylinder surface kind is therefore accepted
  // even when this adapter-shaped fixture carries the MFT detector label.
  {
    auto pool = std::make_shared<BoundedMemoryResource>();
    TimeFrame frame;
    SurfaceTrackingScratch scratch;
    TrackerTraits traits;
    prepareTraversalFrame(frame, scratch, traits, pool, params);
    auto built = buildPlan(catalog(10, SurfaceKind::Cylinder, o2::detectors::DetID::MFT), order(10),
                           SurfaceKind::Cylinder, params);
    adoptPlanBinding(built, traits, 0);
    BOOST_CHECK_NO_THROW(traits.initialiseTimeFrame(0, built.plan));
    BOOST_CHECK(traits.hasTraversalCache());
  }
  params[0].MaxChi2ClusterAttachment = -1.f;
  checkFailure(params, catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), SurfaceKind::Disk,
               TraversalFailureReason::InvalidSurfaceParameters);

  // Perturbing the temporary legacy LayerxX0 away from the catalog's
  // authoritative material now fails the material-compatibility check
  // (LegacyMaterialMismatch) before ever reaching attachHitConfig's own
  // finite/non-negative validation -- the catalog's material is unperturbed
  // and finite, so the mismatch is purely numeric.
  params = mftTraversalParameters();
  params[0].LayerxX0[3] = std::numeric_limits<float>::quiet_NaN();
  checkFailure(params, catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), SurfaceKind::Disk,
               TraversalFailureReason::LegacyMaterialMismatch);

  params = mftTraversalParameters();
  params[0].LayerxX0[7] = std::numeric_limits<float>::infinity();
  checkFailure(params, catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), SurfaceKind::Disk,
               TraversalFailureReason::LegacyMaterialMismatch);

  params = mftTraversalParameters();
  params[0].CorrType = static_cast<o2::base::PropagatorF::MatCorrType>(99);
  checkFailure(params, catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), SurfaceKind::Disk,
               TraversalFailureReason::InvalidSurfaceParameters);
}

BOOST_AUTO_TEST_CASE(every_iteration_resolves_identical_authoritative_material)
{
  auto params = mftTraversalParameters();
  params.push_back(params.front());
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  TrackerTraits traits;
  prepareTraversalFrame(frame, scratch, traits, pool, params);
  auto built = buildPlan(catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), SurfaceKind::Disk, params);

  adoptPlanBinding(built, traits, 0);
  traits.initialiseTimeFrame(0, built.plan);
  const auto firstIterationMaterial = traits.getLayerMaterial();
  std::vector<NominalSurfaceMaterial> firstIteration(firstIterationMaterial.begin(), firstIterationMaterial.end());

  adoptPlanBinding(built, traits, 1);
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
  TrackerTraits traits;
  prepareTraversalFrame(frame, scratch, traits, pool, params);
  auto built = buildPlan(catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), SurfaceKind::Disk, params);
  const NominalSurfaceMaterial materialBefore = built.plan.front().getSurfaceCatalog().getSurface(SurfaceId{4}).material;

  try {
    adoptPlanBinding(built, traits, 0);
    traits.initialiseTimeFrame(0, built.plan);
    BOOST_FAIL("perturbed LayerxX0 must throw");
  } catch (const TraversalException& error) {
    BOOST_CHECK(error.getReason() == TraversalFailureReason::LegacyMaterialMismatch);
  }

  const auto& materialAfter = built.plan.front().getSurfaceCatalog().getSurface(SurfaceId{4}).material;
  BOOST_CHECK_EQUAL(materialAfter.xOverX0, materialBefore.xOverX0);
  BOOST_CHECK_EQUAL(materialAfter.arealDensityGPerCm2, materialBefore.arealDensityGPerCm2);
}

BOOST_AUTO_TEST_CASE(non_monotonic_ordered_surfaces_maps_material_and_traversal_slots)
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
  TrackerTraits traits;
  prepareTraversalFrame(frame, scratch, traits, pool, params);
  auto built = buildPlan(std::move(surfaces), nonMonotonicOrder, SurfaceKind::Disk, params);

  adoptPlanBinding(built, traits, 0);
  BOOST_CHECK_NO_THROW(traits.initialiseTimeFrame(0, built.plan));
  BOOST_REQUIRE(traits.hasTraversalCache());
  const auto resolvedMaterial = traits.getLayerMaterial();
  BOOST_REQUIRE_EQUAL(resolvedMaterial.size(), nonMonotonicOrder.size());
  for (size_t slot = 0; slot < resolvedMaterial.size(); ++slot) {
    const auto& material = resolvedMaterial[slot];
    const auto expected = built.plan.front().getSurfaceCatalog().getSurface(nonMonotonicOrder[slot]).material;
    BOOST_CHECK_EQUAL(material.xOverX0, expected.xOverX0);
    BOOST_CHECK_EQUAL(material.arealDensityGPerCm2, expected.arealDensityGPerCm2);
  }
}

BOOST_AUTO_TEST_CASE(traversal_preflight_rejects_cycles_and_accepts_disconnected_surface_kinds)
{
  auto bind = [](BuiltLayout layout) {
    const auto view = layout.layout.getView();
    const auto ordered = layout.layout.getOrderedSurfaces();
    SurfaceMask owned;
    for (const auto surface : ordered) {
      owned.set(surface);
    }
    return SurfacePlanBinding::build(view, owned, ordered);
  };

  BOOST_CHECK(!bind(cyclicDiskLayout()).ok());

  const auto mixed = bind(mixedDisconnectedLayout());
  BOOST_REQUIRE(mixed.ok());
  BOOST_CHECK_EQUAL(mixed.binding->getOwnedSurfaces().count(), 10);
  BOOST_CHECK_EQUAL(mixed.binding->getGlobalTransitions().size(), 8u);
  BOOST_CHECK_EQUAL(mixed.binding->getGlobalCells().size(), 0u);
}

// The runtime-plan core accepts a valid surface kind independently of the
// detector-labelled adapter that supplied the descriptors.
BOOST_AUTO_TEST_CASE(runtime_plan_accepts_disk_kind_without_detector_family_dispatch)
{
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  TrackerTraits traits;
  auto params = itsTraversalParameters();
  prepareTraversalFrame(frame, scratch, traits, pool, params);
  auto built = buildPlan(catalog(7, SurfaceKind::Disk, o2::detectors::DetID::ITS), order(7),
                         SurfaceKind::Disk, params);
  adoptPlanBinding(built, traits, 0);
  BOOST_CHECK_NO_THROW(traits.initialiseTimeFrame(0, built.plan));
  BOOST_CHECK(traits.hasTraversalCache());
}

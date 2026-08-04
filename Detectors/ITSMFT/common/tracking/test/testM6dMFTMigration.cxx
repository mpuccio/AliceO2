// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// M6d (doc/design/0002-m6-generic-workspace-migration.md Sec 9) focused
// tests for the MFT common-CA participant's migration from
// LegacyTrackerScratch<MFTNLayers>/DetectorTraversalBinding to
// SurfaceTrackingScratch/SurfacePlanBinding. Exercises the real production
// wiring (ITSMFTLegacyParticipantSet, LegacyCATrackingParticipantMFT), not a
// synthetic seam-only fixture -- testSurfacePlanBinding.cxx already covers
// SurfacePlanBinding's own generic slot-mapping equivalence;
// testCombinedTrackingComposition.cxx already covers the full load/track/
// publish composition. This file proves specifically:
//  - the MFT participant's concrete scratch/binding types actually changed
//    (compile-time type proof, not just behavior);
//  - ITS's own types did not (byte-for-byte, same instantiation);
//  - MFT load failure remains atomic across the mixed-storage-type pair;
//  - resetting only MFT's scratch never mutates TimeFrame or ITS's own
//    (differently-typed) scratch;
//  - the production MFT SurfacePlanBinding construction (real combined
//    catalog, real ClusterSourceId{1}/SurfaceKind::Disk/
//    TransitionPolicyTag::DiskDisk parameters) resolves to the same
//    transition/cell slot counts and owned-surface indices the old
//    DetectorTraversalBinding construction at the same parameters would
//    have;
//  - the MFT publication/sidecar export path still works with the new
//    scratch type backing tracking;
//  - ITS and MFT scratch storage stay isolated during the coexistence
//    phase (different concrete types, never cross-referenced).

#define BOOST_TEST_MODULE ITSMFT M6dMFTMigration
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <type_traits>
#include <vector>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/DetectorLayoutBuilder.h"
#include "ITSMFTTracking/ITSMFTLegacyParticipantSet.h"
#include "ITSMFTTracking/LegacyCATrackingParticipant.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/detail/DetectorTraversalBinding.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

// --- 1/2: compile-time type proof -------------------------------------------

static_assert(std::is_same_v<LegacyCATrackingParticipantMFT::ScratchN, SurfaceTrackingScratch>);
static_assert(std::is_same_v<LegacyCATrackingParticipantITS::ScratchN, LegacyTrackerScratchITS>);
static_assert(std::is_same_v<LegacyCATrackingParticipantITS::ScratchN, LegacyTrackerScratch<ITSNLayers>>);
static_assert(!std::is_same_v<LegacyCATrackingParticipantMFT::ScratchN, LegacyTrackerScratch<MFTNLayers>>);
// getScratch() itself, not just the nested alias -- proves the actual member
// type, not merely a type-level artifact nothing reads.
static_assert(std::is_same_v<decltype(std::declval<LegacyCATrackingParticipantMFT&>().getScratch()), SurfaceTrackingScratch&>);
static_assert(std::is_same_v<decltype(std::declval<LegacyCATrackingParticipantITS&>().getScratch()), LegacyTrackerScratchITS&>);
// The binding parameter type LegacyCATrackingParticipantMFT actually adopts.
static_assert(std::is_invocable_v<decltype(&LegacyCATrackingParticipantMFT::adoptDetectorTraversalBinding), LegacyCATrackingParticipantMFT&, std::unique_ptr<SurfacePlanBinding>>);
static_assert(std::is_invocable_v<decltype(&LegacyCATrackingParticipantITS::adoptDetectorTraversalBinding), LegacyCATrackingParticipantITS&, std::unique_ptr<DetectorTraversalBinding>>);
// ITSMFTLegacyParticipantSet's own public readback types.
static_assert(std::is_same_v<decltype(std::declval<const ITSMFTLegacyParticipantSet&>().getMFTScratch()), const SurfaceTrackingScratch&>);
static_assert(std::is_same_v<decltype(std::declval<const ITSMFTLegacyParticipantSet&>().getITSScratch()), const LegacyTrackerScratchITS&>);

BOOST_AUTO_TEST_CASE(CompileTimeTypeProofsHoldAtRuntimeToo)
{
  // The static_asserts above are the real proof (compile-time, so a
  // regression fails the build, not just this test); this case exists only
  // so the module reports it.
  BOOST_CHECK(true);
}

// --- shared fixtures ---------------------------------------------------------

namespace
{

TrackingParameters makeItsParams()
{
  TrackingParameters p;
  resetDetectorDefaults(p, o2::detectors::DetID::ITS);
  return p;
}

TrackingParameters makeMftParams()
{
  TrackingParameters p;
  resetDetectorDefaults(p, o2::detectors::DetID::MFT);
  return p;
}

ITSMFTLegacyParticipantSet makeSet()
{
  return ITSMFTLegacyParticipantSet{std::vector<TrackingParameters>{makeItsParams()},
                                    std::vector<TrackingParameters>{makeMftParams()}};
}

std::vector<SurfaceId> orderedRange(uint16_t first, uint16_t count)
{
  std::vector<SurfaceId> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

// A minimal decoder never actually invoked (every source below carries zero
// clusters): only its address must be a valid, dereferenceable
// ClusterDecoder for loadSources()'s own preflight to accept the source.
class StubDecoder final : public ClusterDecoder
{
 public:
  o2::itsmft::ioutils::SurfaceMeasurementDecodeResult decode(
    const o2::itsmft::CompClusterExt&, BoundedPatternCursor&, const o2::itsmft::TopologyDictionary*,
    gsl::span<const SurfaceId>, ClusterSourceId, uint32_t, uint32_t, bool) const override
  {
    return {};
  }
};

StubDecoder& stubDecoder()
{
  static StubDecoder decoder;
  return decoder;
}

// A structurally valid, empty-cluster source for `det`/`id`, with a
// correctly-sized layerToSurface (`nLayers` entries starting at
// `surfaceOffset`) unless `corruptLayerToSurfaceSize` overrides it -- the one
// deliberate malformation the atomic-failure test below needs.
// `layerToSurfaceStorage` is caller-owned (the returned ClusterSourceInput's
// `layerToSurface` span is non-owning and must not outlive it) -- every call
// site needs its own distinct storage, never a shared/static buffer that a
// second call would silently overwrite out from under the first.
ClusterSourceInput makeEmptySource(ClusterSourceId id, o2::detectors::DetID::ID det, uint16_t surfaceOffset, uint16_t nLayers,
                                   std::vector<SurfaceId>& layerToSurfaceStorage, int corruptLayerToSurfaceSize = -1)
{
  ClusterSourceInput input{};
  input.id = id;
  input.detector = det;
  input.timing = ROFTimingConfig{40, 0, 0, 0};
  input.decoder = &stubDecoder();
  const auto count = corruptLayerToSurfaceSize >= 0 ? static_cast<uint16_t>(corruptLayerToSurfaceSize) : nLayers;
  layerToSurfaceStorage = orderedRange(surfaceOffset, count);
  input.layerToSurface = layerToSurfaceStorage;
  return input;
}

} // namespace

// --- 3: atomic MFT load failure ----------------------------------------------

BOOST_AUTO_TEST_CASE(AtomicMFTLoadFailureLeavesSharedTimeFrameAndBothParticipantScratchesUntouched)
{
  auto participants = makeSet();
  TimeFrame frame;
  participants.adoptFrame(frame);

  // ITS's own source is structurally valid (correctly-sized layerToSurface);
  // MFT's is deliberately malformed (9 entries instead of MFTNLayers=10) --
  // SurfaceTrackingScratch::loadNormalizedSource()'s own preflight
  // (orderedSurfaces.size() != mNOwnedSurfaces) must reject it before any
  // commit.
  std::vector<SurfaceId> itsLayerToSurfaceStorage;
  std::vector<SurfaceId> mftLayerToSurfaceStorage;
  const auto itsSource = makeEmptySource(ClusterSourceId{0}, o2::detectors::DetID::ITS, 0, ITSNLayers, itsLayerToSurfaceStorage);
  const auto mftSource = makeEmptySource(ClusterSourceId{1}, o2::detectors::DetID::MFT, ITSNLayers, MFTNLayers, mftLayerToSurfaceStorage, MFTNLayers - 1);

  BOOST_REQUIRE(!participants.validateSources(itsSource, mftSource).has_value());
  const auto bindings = participants.loadBindings(itsSource, mftSource);
  const auto result = MultiSourceTimeFrameLoader::loadEvent(
    frame, gsl::span<const MultiSourceTimeFrameLoader::AtomicLoadBinding>{bindings}, participants.catalogView(), o2::InteractionRecord{50, 5});

  BOOST_REQUIRE(!result.ok());
  // Both LegacyTrackerScratch<NLayers>::loadNormalizedSource() and
  // SurfaceTrackingScratch::loadNormalizedSource() re-stage every source
  // under a fixed ClusterSourceId{0} (each stage() call rewrites `one.id`
  // before delegating, mirroring the single-source dense-ID contract), so
  // the failing stage's reported source is always 0, not the caller-facing
  // MFT binding position (1).
  BOOST_CHECK(result.source == ClusterSourceId{0});
  BOOST_CHECK(result.error == MultiSourceLoadError::InvalidLayerMapping);

  // Nothing committed anywhere: not the shared normalized frame, not ITS's
  // own scratch (its own stage() would have succeeded individually, exactly
  // the same "earlier participant's target left untouched" property
  // testMultiSourceTimeFrameLoader.cxx proves generically -- this is that
  // same property proven across the real, differently-typed ITS/MFT pair).
  BOOST_CHECK_EQUAL(frame.getNormalizedFrame().getTotalMeasurements(), 0u);
  BOOST_CHECK_EQUAL(participants.getITSScratch().getTotalClusters(), 0);
  BOOST_CHECK_EQUAL(participants.getMFTScratch().getTotalClusters(), 0);
}

// --- 4: MFT-only reset isolation ---------------------------------------------

BOOST_AUTO_TEST_CASE(MFTOnlyResetDoesNotMutateTimeFrameOrITSScratch)
{
  auto participants = makeSet();
  // Real production code (CombinedCATrackerSpec.cxx) always calls this
  // before any scratch container is touched -- needed here too since this
  // test writes into scratch containers directly rather than through a real
  // load()/track() call.
  participants.setMemoryPool(std::make_shared<o2::its::BoundedMemoryResource>());
  TimeFrame frame;
  participants.adoptFrame(frame);
  frame.setBz(5.f);
  frame.setBeamPosition(1.f, 2.f, 0.1f);

  // Populate observable working state directly on each scratch (no real
  // load/track needed for this isolation proof). Only getScratch()'s const
  // overload is reachable through ITSMFTLegacyParticipantSet's own public
  // API; const_cast here is test-only instrumentation, mirroring how
  // testCombinedTrackingComposition.cxx's own composer wrapper reaches into
  // scratch state for assertions.
  auto& itsParticipantScratch = const_cast<LegacyTrackerScratchITS&>(participants.getITSScratch());
  auto& mftParticipantScratch = const_cast<SurfaceTrackingScratch&>(participants.getMFTScratch());
  itsParticipantScratch.getUnsortedClusters()[0].emplace_back(1.f, 2.f, 3.f, 0);
  BOOST_REQUIRE(!mftParticipantScratch.getUnsortedClusters().empty());
  mftParticipantScratch.getUnsortedClusters()[0].emplace_back(4.f, 5.f, 6.f, 0);
  const auto mftPlanSurfaces = mftParticipantScratch.getNOwnedSurfaces();

  // eventReset() on the MFT TrackingParticipant interface only -- never
  // touches `frame` (TrackingParticipant.h's own contract) and, by
  // construction (ITSMFTLegacyParticipantSet composes two independent
  // participant objects of different concrete scratch types), cannot reach
  // ITS's own scratch either.
  TrackingParticipant& mftParticipant = *participants.schedule()[1];
  BOOST_REQUIRE(mftParticipant.id() == ParticipantId{1});
  mftParticipant.eventReset(frame);

  BOOST_CHECK(mftParticipantScratch.getUnsortedClusters()[0].empty());
  BOOST_CHECK_EQUAL(mftParticipantScratch.getNOwnedSurfaces(), mftPlanSurfaces); // reset() never un-adopts the plan

  // ITS's own scratch and the shared TimeFrame are both completely
  // unaffected.
  BOOST_CHECK_EQUAL(itsParticipantScratch.getUnsortedClusters()[0].size(), 1u);
  BOOST_CHECK_EQUAL(frame.getBz(), 5.f);
  BOOST_CHECK_EQUAL(frame.getBeamX(), 1.f);
}

// --- 5: production MFT SurfacePlanBinding matches the old construction -----

namespace
{
SurfaceMask surfaceRangeMaskForTest(uint16_t first, uint16_t count)
{
  SurfaceMask result;
  for (uint16_t i = 0; i < count; ++i) {
    result.set(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

// Byte-for-byte the same combined-layout construction
// ITSMFTLegacyParticipantSet.cxx's own buildCombinedLayout() performs
// (duplicated locally per this test directory's own established per-file
// fixture convention -- testSurfacePlanBinding.cxx's CombinedLayout does the
// same for its own synthetic case).
DetectorLayout buildProductionCombinedLayoutForTest()
{
  const auto itsSurfaces = orderedRange(0, ITSNLayers);
  const auto mftSurfaces = orderedRange(ITSNLayers, MFTNLayers);
  const auto itsParams = makeItsParams();
  const auto mftParams = makeMftParams();

  DetectorLayoutSubgraph itsSubgraph;
  itsSubgraph.orderedSurfaces.assign(itsSurfaces.begin(), itsSurfaces.end());
  itsSubgraph.maxHoles = itsParams.MaxHoles;
  itsSubgraph.holeSurfaces = positionalSurfaceMask(itsParams.HoleLayerMask, itsSurfaces, static_cast<uint32_t>(itsSurfaces.size()));
  itsSubgraph.seedingSurfaces = positionalSurfaceMask(itsParams.StartLayerMask, itsSurfaces, static_cast<uint32_t>(itsSurfaces.size()));

  DetectorLayoutSubgraph mftSubgraph;
  mftSubgraph.orderedSurfaces.assign(mftSurfaces.begin(), mftSurfaces.end());
  mftSubgraph.maxHoles = mftParams.MaxHoles;
  mftSubgraph.holeSurfaces = positionalSurfaceMask(mftParams.HoleLayerMask, mftSurfaces, static_cast<uint32_t>(mftSurfaces.size()));
  mftSubgraph.seedingSurfaces = positionalSurfaceMask(mftParams.StartLayerMask, mftSurfaces, static_cast<uint32_t>(mftSurfaces.size()));

  const SurfaceCatalogView catalog{kITSMFTCombinedStaticSurfaceCatalog.data(), static_cast<uint32_t>(kITSMFTCombinedStaticSurfaceCatalog.size())};
  DetectorLayoutBuilder builder{catalog};
  builder.addSubgraph(std::move(itsSubgraph));
  builder.addSubgraph(std::move(mftSubgraph));
  auto built = builder.build();
  BOOST_REQUIRE(built.ok());
  return std::move(*built.layout);
}
} // namespace

BOOST_AUTO_TEST_CASE(ProductionMFTSurfacePlanBindingMatchesDetectorTraversalBindingAtRealParameters)
{
  // Exactly the parameters ITSMFTLegacyParticipantSet.cxx's own constructor
  // uses for MFT: ClusterSourceId{1}, surfaceRangeMask(ITSNLayers,
  // MFTNLayers), the MFT ordered-surface range, SurfaceKind::Disk,
  // TransitionPolicyTag::DiskDisk.
  const auto layout = buildProductionCombinedLayoutForTest();
  const auto masks = computeSurfaceKindMasks(kITSMFTCombinedStaticSurfaceCatalog);
  const auto view = layout.getView(kITSMFTCombinedStaticSurfaceCatalog, masks.first, masks.second);
  const auto mftSurfaces = orderedRange(ITSNLayers, MFTNLayers);
  const auto mftMask = surfaceRangeMaskForTest(ITSNLayers, MFTNLayers);

  const auto oldBinding = DetectorTraversalBinding::build(view, o2::detectors::DetID::MFT, ClusterSourceId{1}, mftMask, mftSurfaces);
  const auto newBinding = SurfacePlanBinding::build(view, ClusterSourceId{1}, mftMask, mftSurfaces, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk);
  BOOST_REQUIRE(oldBinding.ok());
  BOOST_REQUIRE(newBinding.ok());

  // Owned-surface count feeds SurfaceTrackingScratch::adoptPlan()'s Group A
  // sizing; transition/cell counts feed its Group B sizing -- these three
  // runtime numbers are exactly what production wiring
  // (LegacyCATrackingParticipant<...>::adoptDetectorLayoutSet()) passes
  // through, so agreement here is the concrete "resolves to the same
  // compact slots" proof at real parameters, not just the generic synthetic
  // fixture testSurfacePlanBinding.cxx already covers.
  BOOST_CHECK_EQUAL(newBinding.binding->getOwnedSurfaces().count(), static_cast<int>(mftSurfaces.size()));
  BOOST_CHECK_EQUAL(newBinding.binding->getGlobalTransitions().size(), oldBinding.binding->getGlobalTransitions().size());
  BOOST_CHECK_EQUAL(newBinding.binding->getGlobalCells().size(), oldBinding.binding->getGlobalCells().size());

  for (uint16_t s = ITSNLayers; s < ITSNLayers + MFTNLayers; ++s) {
    const auto oldSlot = oldBinding.binding->getLegacyLayer(SurfaceId{s});
    const auto newSlot = newBinding.binding->getOwnedSurfaceIndex(SurfaceId{s});
    BOOST_REQUIRE_EQUAL(oldSlot.has_value(), newSlot.has_value());
    if (oldSlot.has_value()) {
      BOOST_CHECK_EQUAL(*oldSlot, *newSlot);
    }
  }
  for (size_t t = 0; t < oldBinding.binding->getGlobalTransitions().size(); ++t) {
    BOOST_CHECK(oldBinding.binding->getGlobalTransitions()[t] == newBinding.binding->getGlobalTransitions()[t]);
  }
  for (size_t c = 0; c < oldBinding.binding->getGlobalCells().size(); ++c) {
    BOOST_CHECK(oldBinding.binding->getGlobalCells()[c] == newBinding.binding->getGlobalCells()[c]);
  }

  // The real ITSMFTLegacyParticipantSet construction (production code path)
  // adopts a plan whose sizing agrees with this direct comparison.
  auto participants = makeSet();
  BOOST_CHECK_EQUAL(participants.getMFTScratch().getNOwnedSurfaces(), static_cast<size_t>(mftSurfaces.size()));
  BOOST_CHECK_EQUAL(participants.getMFTScratch().getNTransitions(), newBinding.binding->getGlobalTransitions().size());
  BOOST_CHECK_EQUAL(participants.getMFTScratch().getNCells(), newBinding.binding->getGlobalCells().size());
}

// --- 6: MFT sidecar/publication export remain valid --------------------------

BOOST_AUTO_TEST_CASE(MFTSidecarAndPublicationExportRemainValidAfterMigration)
{
  auto participants = makeSet();
  TimeFrame frame;
  participants.adoptFrame(frame);

  std::vector<SurfaceId> itsLayerToSurfaceStorage;
  std::vector<SurfaceId> mftLayerToSurfaceStorage;
  const auto itsSource = makeEmptySource(ClusterSourceId{0}, o2::detectors::DetID::ITS, 0, ITSNLayers, itsLayerToSurfaceStorage);
  const auto mftSource = makeEmptySource(ClusterSourceId{1}, o2::detectors::DetID::MFT, ITSNLayers, MFTNLayers, mftLayerToSurfaceStorage);
  participants.configureRofTables(itsSource, mftSource);

  // The MFT compatibility sidecar itself is still reachable and clearable --
  // ownership/lifetime unaffected by the scratch/binding type swap (M6d
  // requirement 4: "MFT compatibility sidecar ownership" preserved).
  const auto& sidecarBefore = participants.getMFTPublicationCompatibility();
  (void)sidecarBefore;
  participants.clearPublicationSidecars();

  // markPublicationValid() reads mMFTParticipant.getScratch().getROFOverlapTableView().getClockLayer()
  // -- the exact call site this migration's own SurfaceTrackingScratch
  // extension (getROFOverlapTableView()) exists to keep working.
  participants.markPublicationValid();
  const auto mftExport = participants.getMFTPublicationExport();
  BOOST_REQUIRE(mftExport.has_value());
  BOOST_CHECK(mftExport->detector == o2::detectors::DetID::MFT);
  BOOST_CHECK(mftExport->source == ClusterSourceId{1});
  BOOST_CHECK_EQUAL(mftExport->orderedSurfaces.size(), static_cast<size_t>(MFTNLayers));

  participants.invalidatePublication();
  BOOST_CHECK(!participants.getMFTPublicationExport().has_value());
}

// --- 7: ITS/MFT scratch storage isolation during coexistence ---------------

BOOST_AUTO_TEST_CASE(CombinedExecutionKeepsITSAndMFTScratchStorageIsolated)
{
  // Different concrete C++ types (LegacyTrackerScratch<ITSNLayers> vs
  // SurfaceTrackingScratch) structurally cannot alias or be confused for one
  // another -- this test proves the runtime corollary: mutating one leaves
  // the other's every observable count exactly where it started, for both
  // directions.
  auto participants = makeSet();
  participants.setMemoryPool(std::make_shared<o2::its::BoundedMemoryResource>());
  TimeFrame frame;
  participants.adoptFrame(frame);

  auto& itsScratch = const_cast<LegacyTrackerScratchITS&>(participants.getITSScratch());
  auto& mftScratch = const_cast<SurfaceTrackingScratch&>(participants.getMFTScratch());

  const auto mftClustersBefore = mftScratch.getTotalClusters();
  const auto mftOwnedBefore = mftScratch.getNOwnedSurfaces();
  itsScratch.getUnsortedClusters()[0].emplace_back(1.f, 1.f, 1.f, 0);
  BOOST_CHECK_EQUAL(mftScratch.getTotalClusters(), mftClustersBefore);
  BOOST_CHECK_EQUAL(mftScratch.getNOwnedSurfaces(), mftOwnedBefore);

  const auto itsClustersBefore = itsScratch.getTotalClusters();
  mftScratch.getUnsortedClusters()[0].emplace_back(2.f, 2.f, 2.f, 0);
  BOOST_CHECK_EQUAL(itsScratch.getTotalClusters(), itsClustersBefore);

  // Resetting one never touches the other or the shared TimeFrame (already
  // proven in detail above for the MFT direction; this closes the loop for
  // the ITS direction too, at the composed-set level).
  frame.setBz(9.f);
  TrackingParticipant& itsParticipant = *participants.schedule()[0];
  itsParticipant.eventReset(frame);
  BOOST_CHECK_EQUAL(mftScratch.getTotalClusters(), 1); // the mutation just above survives
  BOOST_CHECK_EQUAL(frame.getBz(), 9.f);
}

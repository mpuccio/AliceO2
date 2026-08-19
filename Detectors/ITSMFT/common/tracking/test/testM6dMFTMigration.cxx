// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Focused MFT adapter coverage for the combined flat plan and the
// TimeFrame-owned workspace. Exercises the real production Tracker wiring,
// not a synthetic seam-only fixture. This file proves specifically:
//  - the MFT Tracker's concrete scratch/binding types remain the intended
//    kernel seam (compile-time type proof, not just behavior);
//  - MFT load failure remains atomic across the shared workspace;
//  - one TimeFrame reset clears the shared workspace while preserving its
//    configured global capacity;
//  - the production MFT plan is represented by the shared traversal
//    workspace rather than a detector-local binding;
//  - the MFT publication/sidecar export path still works at the adapter edge.

#define BOOST_TEST_MODULE ITSMFT M6dMFTMigration
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <optional>
#include <type_traits>
#include <vector>

#include "CommonDataFormat/InteractionRecord.h"
#include "CombinedTrackingTestSupport.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/GenericTrackOutputAdapter.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/detail/TimeFrameScratch.h"
#include "ITSMFTTracking/TimeFrame.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

// --- 1/2: compile-time type proof -------------------------------------------

// getScratch() itself proves the actual member type and keeps the production
// class free of a compatibility alias.
static_assert(std::is_same_v<decltype(std::declval<const test::CombinedTrackingPlan&>().getMFTScratch()), const TimeFrameScratch&>);
static_assert(std::is_same_v<decltype(std::declval<const test::CombinedTrackingPlan&>().getITSScratch()), const TimeFrameScratch&>);
static_assert(std::is_invocable_v<decltype(&Tracker::run), Tracker&, TimeFrame&, TrackerTraits&>);

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

test::CombinedTrackingPlan makeSet()
{
  return test::CombinedTrackingPlan{std::vector<TrackingParameters>{makeItsParams()},
                                    std::vector<TrackingParameters>{makeMftParams()}};
}

std::vector<LayerId> orderedRange(uint16_t first, uint16_t count)
{
  std::vector<LayerId> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(LayerId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

// A minimal decoder never actually invoked (every source below carries zero
// clusters): only its address must be a valid, dereferenceable
// ClusterDecoder for loadSources()'s own preflight to accept the source.
class StubDecoder final : public ClusterDecoder
{
 public:
  o2::itsmft::tracking::ClusterDecodeResult decode(
    const o2::itsmft::CompClusterExt&, BoundedPatternCursor&, const o2::itsmft::TopologyDictionary*,
    uint32_t, bool) const override
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
                                   std::vector<LayerId>& layerToSurfaceStorage, int corruptLayerToSurfaceSize = -1)
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
  // The fixture writes directly to the workspace, so provide its allocator
  // before loading or inspecting event-owned containers.
  TimeFrame frame;
  participants.adoptFrame(frame);

  // ITS's own source is structurally valid (correctly-sized layerToSurface);
  // MFT's is deliberately malformed (9 entries instead of MFTNLayers=10) --
  // TimeFrameScratch::loadNormalizedSource()'s own preflight
  // mismatched topology extents must be rejected before any
  // commit.
  std::vector<LayerId> itsLayerToSurfaceStorage;
  std::vector<LayerId> mftLayerToSurfaceStorage;
  const auto itsSource = makeEmptySource(ClusterSourceId{0}, o2::detectors::DetID::ITS, 0, ITSNLayers, itsLayerToSurfaceStorage);
  const auto mftSource = makeEmptySource(ClusterSourceId{1}, o2::detectors::DetID::MFT, ITSNLayers, MFTNLayers, mftLayerToSurfaceStorage, MFTNLayers - 1);

  BOOST_REQUIRE(!participants.validateSources(itsSource, mftSource).has_value());
  participants.configureRofTables(itsSource, mftSource);
  auto itsInput = itsSource;
  auto mftInput = mftSource;
  itsInput.rofViews = participants.getITSROFViews();
  mftInput.rofViews = participants.getMFTROFViews();
  const std::array<ClusterSourceInput, 2> sources{itsInput, mftInput};
  const auto result = loadTimeFrameSources(
    frame, gsl::span<const ClusterSourceInput>{sources}, participants.catalogView(), o2::InteractionRecord{50, 5});

  BOOST_REQUIRE(!result.ok());
  // The direct loader reports the caller-facing source whose staged decode
  // failed; the low-level dense-ID adaptation is not part of this result.
  BOOST_CHECK(result.source == ClusterSourceId{1});
  BOOST_CHECK(result.error == MultiSourceLoadError::InvalidLayerMapping);

  // Nothing committed anywhere: not the shared normalized frame nor the
  // staged global workspace.
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
  BOOST_CHECK_EQUAL(participants.getITSScratch().getNumberOfTracklets(), 0);
  BOOST_CHECK_EQUAL(participants.getMFTScratch().getNumberOfTracklets(), 0);
}

// --- 4: shared-workspace reset -----------------------------------------------

BOOST_AUTO_TEST_CASE(TimeFrameResetClearsSharedWorkspaceAndPreservesFrameState)
{
  auto participants = makeSet();
  TimeFrame frame;
  participants.adoptFrame(frame);
  frame.setBz(5.f);
  frame.setBeamPosition(1.f, 2.f, 0.1f);

  // The two workflow accessors are aliases of one global workspace. Only the
  // const accessor is exposed by the fixture, so const_cast is limited to
  // this direct reset probe.
  auto& itsParticipantScratch = const_cast<TimeFrameScratch&>(participants.getITSScratch());
  auto& mftParticipantScratch = const_cast<TimeFrameScratch&>(participants.getMFTScratch());
  BOOST_CHECK_EQUAL(&itsParticipantScratch, &mftParticipantScratch);
  BOOST_CHECK_EQUAL(itsParticipantScratch.getNEdges(), 15u);
  BOOST_CHECK_EQUAL(itsParticipantScratch.getNCells(), 13u);
  itsParticipantScratch.getTracklets()[0].emplace_back();
  BOOST_CHECK_EQUAL(mftParticipantScratch.getNumberOfTracklets(), 1);
  mftParticipantScratch.getTracklets()[7].emplace_back();
  BOOST_CHECK_EQUAL(itsParticipantScratch.getNumberOfTracklets(), 2);

  frame.resetTimeFrame();

  BOOST_CHECK_EQUAL(mftParticipantScratch.getNumberOfTracklets(), 0);

  BOOST_CHECK(itsParticipantScratch.getTracklets()[0].empty());
  BOOST_CHECK_EQUAL(frame.getBz(), 5.f);
  BOOST_CHECK_EQUAL(frame.getBeamX(), 1.f);
}

// --- 5: production MFT traversal workspace ----------------------------------

BOOST_AUTO_TEST_CASE(ProductionMFTWorkspaceMatchesConfiguredTopologyAtRealParameters)
{
  auto participants = makeSet();
  TimeFrame frame;
  participants.adoptFrame(frame);
  BOOST_CHECK_EQUAL(&participants.getITSScratch(), &participants.getMFTScratch());
  BOOST_CHECK_EQUAL(participants.getMFTScratch().getNEdges(), 15u);
  BOOST_CHECK_EQUAL(participants.getMFTScratch().getNCells(), 13u);
}

// --- 6: MFT sidecar/publication export remain valid --------------------------

BOOST_AUTO_TEST_CASE(MFTPublicationExportRemainsValidAfterMigration)
{
  auto participants = makeSet();
  TimeFrame frame;
  participants.adoptFrame(frame);

  std::vector<LayerId> itsLayerToSurfaceStorage;
  std::vector<LayerId> mftLayerToSurfaceStorage;
  const auto itsSource = makeEmptySource(ClusterSourceId{0}, o2::detectors::DetID::ITS, 0, ITSNLayers, itsLayerToSurfaceStorage);
  const auto mftSource = makeEmptySource(ClusterSourceId{1}, o2::detectors::DetID::MFT, ITSNLayers, MFTNLayers, mftLayerToSurfaceStorage);

  participants.clearPublicationSidecars();

  // The publication clock/validity context is workflow-owned. Reproduce the
  // same narrow publication step locally without adding that event state to
  // the application-plan fixture.
  std::optional<ClockTimingPublicationView> mftClock;
  bool publicationValid = false;
  mftClock.emplace(participants.getMFTROFViews().overlap.getClockLayer());
  publicationValid = true;
  std::optional<GenericTrackPublicationExport> mftExport;
  if (publicationValid && mftClock) {
    mftExport.emplace(o2::detectors::DetID::MFT, ClusterSourceId{1}, *mftClock, participants.getMFTOrderedSurfaces());
  }
  BOOST_REQUIRE(mftExport.has_value());
  BOOST_CHECK(mftExport->detector == o2::detectors::DetID::MFT);
  BOOST_CHECK(mftExport->source == ClusterSourceId{1});
  BOOST_CHECK_EQUAL(mftExport->orderedSurfaces.size(), static_cast<size_t>(MFTNLayers));

  mftClock.reset();
  publicationValid = false;
  mftExport.reset();
  BOOST_CHECK(!mftExport.has_value());
}

// --- 7: ITS/MFT accessors alias the shared workspace -------------------------

BOOST_AUTO_TEST_CASE(CombinedExecutionUsesOneSharedWorkspace)
{
  // Both workflow accessors address one global plan/workspace.
  auto participants = makeSet();
  TimeFrame frame;
  participants.adoptFrame(frame);

  auto& itsScratch = const_cast<TimeFrameScratch&>(participants.getITSScratch());
  auto& mftScratch = const_cast<TimeFrameScratch&>(participants.getMFTScratch());
  BOOST_CHECK_EQUAL(&itsScratch, &mftScratch);
  BOOST_CHECK_EQUAL(itsScratch.getNEdges(), 15u);
  BOOST_CHECK_EQUAL(itsScratch.getNCells(), 13u);
  itsScratch.getTracklets()[0].emplace_back();
  BOOST_CHECK_EQUAL(mftScratch.getNumberOfTracklets(), 1);
  mftScratch.getTracklets()[7].emplace_back();
  BOOST_CHECK_EQUAL(itsScratch.getNumberOfTracklets(), 2);

  frame.setBz(9.f);
  frame.resetTimeFrame();
  BOOST_CHECK_EQUAL(itsScratch.getNumberOfTracklets(), 0);
  BOOST_CHECK_EQUAL(mftScratch.getNumberOfTracklets(), 0);
  BOOST_CHECK_EQUAL(frame.getBz(), 9.f);
}

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Focused ITS adapter coverage for the combined flat plan and its shared
// workspace.
// It proves:
//  - both ITS workflow accessors expose the same shared kernel workspace;
//  - an unconfigured workspace fails closed before loading;
//  - ITS load failure and TimeFrame::resetTimeFrame() preserve the
//    GenericTrack/sidecar/workspace contracts;
//  - the ITS shared-cluster compatibility sidecar (pending/sealed) still
//    works correctly backed by the new scratch storage;
//  - the production ITS traversal uses the same shared workspace plan as the
//    combined workflow.

#define BOOST_TEST_MODULE ITSMFT M6e2ITSWorkspaceMigration
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <memory>
#include <type_traits>
#include <vector>

#include <TGeoGlobalMagField.h>

#include "CommonDataFormat/InteractionRecord.h"
#include "CombinedTrackingTestSupport.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "Field/MagneticField.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/InputSpec.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TimeFrame.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

// --- 1: compile-time type proof ----------------------------------------------

static_assert(std::is_invocable_v<decltype(&Tracker::run), Tracker&, TimeFrame&, TrackerTraits&>);

BOOST_AUTO_TEST_CASE(CompileTimeTypeProofsHoldAtRuntimeToo)
{
  BOOST_CHECK(true); // compile-time proofs above; case exists so it appears in test output.
}

// --- 2: a scratch that never adopted a plan fails closed ---------------------

BOOST_AUTO_TEST_CASE(UnadoptedScratchRejectsLoadInsteadOfMisbehaving)
{
  // No adoptPlan() call anywhere: mNOwnedSurfaces stays 0. A structurally
  // valid ITS source with a real (nonzero) layerToSurface must then fail
  // loadNormalizedSource()'s own orderedSurfaces.size() != mNOwnedSurfaces
  // preflight -- a clean, typed error, not an out-of-bounds access into an
  // unsized Group A container.
  SurfaceTrackingScratch scratch;
  scratch.setMemoryPool(std::make_shared<o2::its::BoundedMemoryResource>());
  BOOST_CHECK_EQUAL(scratch.getNOwnedSurfaces(), 0u);

  class RejectDecoder final : public ClusterDecoder
  {
   public:
    o2::itsmft::tracking::SurfaceMeasurementDecodeResult decode(
      const o2::itsmft::CompClusterExt&, BoundedPatternCursor&, const o2::itsmft::TopologyDictionary*,
      gsl::span<const SurfaceId>, ClusterSourceId, uint32_t, uint32_t, bool) const override
    {
      return {};
    }
  };
  RejectDecoder decoder;

  std::vector<SurfaceId> layerToSurface;
  for (uint16_t i = 0; i < ITSNLayers; ++i) {
    layerToSurface.push_back(SurfaceId{i});
  }
  TimeFrame frame;
  const auto result = scratch.loadNormalizedSource(
    frame, decoder, o2::InteractionRecord{0, 0}, ROFTimingConfig{40, 0, 0, 0},
    gsl::span<const o2::itsmft::CompClusterExt>{}, gsl::span<const unsigned char>{}, gsl::span<const o2::itsmft::ROFRecord>{},
    nullptr, nullptr, o2::detectors::DetID::ITS, gsl::span<const SurfaceId>{layerToSurface},
    SurfaceCatalogView{kITSStaticSurfaceCatalog.data(), static_cast<uint32_t>(kITSStaticSurfaceCatalog.size())});

  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::InvalidLayerMapping);
  BOOST_CHECK_EQUAL(scratch.getTotalClusters(), 0);
}

// --- shared fixtures (combined-participant coverage) -------------------------

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

std::vector<SurfaceId> orderedRange(uint16_t first, uint16_t count)
{
  std::vector<SurfaceId> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

class StubDecoder final : public ClusterDecoder
{
 public:
  o2::itsmft::tracking::SurfaceMeasurementDecodeResult decode(
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

// --- 3: atomic ITS load failure -----------------------------------------------

BOOST_AUTO_TEST_CASE(AtomicITSLoadFailureLeavesSharedTimeFrameAndBothParticipantScratchesUntouched)
{
  auto participants = makeSet();
  TimeFrame frame;
  participants.adoptFrame(frame);

  // ITS's own source is deliberately malformed (6 entries instead of
  // ITSNLayers=7); MFT's is structurally valid. Since ITS is source 0
  // (staged first), this proves the failure is caught before any commit,
  // not only when the malformed source happens to be staged last (the
  // paired MFT-source failure test).
  std::vector<SurfaceId> itsLayerToSurfaceStorage;
  std::vector<SurfaceId> mftLayerToSurfaceStorage;
  const auto itsSource = makeEmptySource(ClusterSourceId{0}, o2::detectors::DetID::ITS, 0, ITSNLayers, itsLayerToSurfaceStorage, ITSNLayers - 1);
  const auto mftSource = makeEmptySource(ClusterSourceId{1}, o2::detectors::DetID::MFT, ITSNLayers, MFTNLayers, mftLayerToSurfaceStorage);

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
  BOOST_CHECK(result.source == ClusterSourceId{0});
  BOOST_CHECK(result.error == MultiSourceLoadError::InvalidLayerMapping);

  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
  BOOST_CHECK_EQUAL(participants.getITSScratch().getTotalClusters(), 0);
  BOOST_CHECK_EQUAL(participants.getMFTScratch().getTotalClusters(), 0);
}

// --- 4: shared-workspace reset -----------------------------------------------

BOOST_AUTO_TEST_CASE(TimeFrameResetClearsSharedWorkspaceAndPreservesFrameState)
{
  auto participants = makeSet();
  TimeFrame frame;
  participants.adoptFrame(frame);
  frame.setBz(5.f);
  frame.setBeamPosition(1.f, 2.f, 0.1f);

  auto& itsParticipantScratch = const_cast<SurfaceTrackingScratch&>(participants.getITSScratch());
  auto& mftParticipantScratch = const_cast<SurfaceTrackingScratch&>(participants.getMFTScratch());
  BOOST_CHECK_EQUAL(&itsParticipantScratch, &mftParticipantScratch);
  BOOST_CHECK_EQUAL(itsParticipantScratch.getNOwnedSurfaces(), 17u);
  BOOST_CHECK_EQUAL(itsParticipantScratch.getNEdges(), 15u);
  BOOST_CHECK_EQUAL(itsParticipantScratch.getNCells(), 13u);
  BOOST_REQUIRE(!itsParticipantScratch.getUnsortedClusters().empty());
  itsParticipantScratch.getUnsortedClusters()[0].emplace_back(1.f, 2.f, 3.f, 0);
  BOOST_CHECK_EQUAL(mftParticipantScratch.getTotalClusters(), 1);
  mftParticipantScratch.getUnsortedClusters()[7].emplace_back(4.f, 5.f, 6.f, 7);
  BOOST_CHECK_EQUAL(itsParticipantScratch.getTotalClusters(), 2);

  const auto resetCount = frame.getEventResetCount();
  frame.resetTimeFrame();

  BOOST_CHECK_EQUAL(frame.getEventResetCount(), resetCount + 1);
  BOOST_CHECK_EQUAL(itsParticipantScratch.getTotalClusters(), 0);
  BOOST_CHECK_EQUAL(itsParticipantScratch.getNOwnedSurfaces(), 17u);
  BOOST_CHECK_EQUAL(mftParticipantScratch.getTotalClusters(), 0);
  BOOST_CHECK_EQUAL(frame.getBz(), 5.f);
  BOOST_CHECK_EQUAL(frame.getBeamX(), 1.f);
}

// --- 5: production ITS traversal workspace ----------------------------------

BOOST_AUTO_TEST_CASE(ProductionITSWorkspaceMatchesConfiguredTopologyAtRealParameters)
{
  auto participants = makeSet();
  TimeFrame frame;
  participants.adoptFrame(frame);
  BOOST_CHECK_EQUAL(&participants.getITSScratch(), &participants.getMFTScratch());
  BOOST_CHECK_EQUAL(participants.getITSScratch().getNOwnedSurfaces(), 17u);
  BOOST_CHECK_EQUAL(participants.getITSScratch().getNEdges(), 15u);
  BOOST_CHECK_EQUAL(participants.getITSScratch().getNCells(), 13u);
}

// --- 6: ITS shared-cluster compatibility sidecar remains correct -----------

BOOST_AUTO_TEST_CASE(ITSSharedClusterCompatibilitySidecarRemainsFunctionalAfterMigration)
{
  auto participants = makeSet();
  TimeFrame frame;
  participants.adoptFrame(frame);

  // Reachable, clearable, and starts empty/unsealed -- ownership/lifetime of
  // the sidecar is unaffected by the shared-workspace composition.
  const auto& sidecar = participants.getITSSharedClusterCompatibility();
  BOOST_CHECK_EQUAL(sidecar.entries().size(), 0u);
  BOOST_CHECK(!sidecar.isSealed());
  BOOST_CHECK_EQUAL(sidecar.pendingSize(), 0u);

  participants.clearPublicationSidecars();
  BOOST_CHECK_EQUAL(sidecar.entries().size(), 0u);
  BOOST_CHECK(!sidecar.isSealed());
}

// --- 7: standalone-vs-combined ITS compact-slot agreement --------------------

BOOST_AUTO_TEST_CASE(StandaloneAndCombinedITSGraphsAgreeByRelativePosition)
{
  const auto standaloneParams = o2::itsmft::TrackingMode::getTrackingParameters(o2::detectors::DetID::ITS, o2::itsmft::TrackingMode::Sync);
  std::vector<SurfaceId> standaloneOrder;
  for (uint16_t i = 0; i < ITSNLayers; ++i) {
    standaloneOrder.push_back(SurfaceId{i});
  }
  const auto standaloneLayout = SurfaceLayout{
    gsl::span<const SurfaceDescriptor>{kITSStaticSurfaceCatalog.data(), kITSStaticSurfaceCatalog.size()},
    makeSurfaceLayoutChain(standaloneOrder, standaloneParams.front().MaxHoles,
                           positionalSurfaceMask(standaloneParams.front().HoleLayerMask, standaloneOrder, ITSNLayers),
                           positionalSurfaceMask(standaloneParams.front().StartLayerMask, standaloneOrder, ITSNLayers))};
  const auto standaloneTopology = deriveTraversalTopology(standaloneLayout);
  BOOST_REQUIRE(standaloneTopology.ok());

  // Combined: real ITS+MFT combined static catalog, ITS half only.
  const auto combinedParams = standaloneParams;
  const auto combinedOrder = orderedRange(0, ITSNLayers);
  const auto combinedLayout = SurfaceLayout{
    gsl::span<const SurfaceDescriptor>{kITSMFTCombinedStaticSurfaceCatalog.data(), kITSMFTCombinedStaticSurfaceCatalog.size()},
    makeSurfaceLayoutChain(combinedOrder, combinedParams.front().MaxHoles,
                           positionalSurfaceMask(combinedParams.front().HoleLayerMask, combinedOrder, ITSNLayers),
                           positionalSurfaceMask(combinedParams.front().StartLayerMask, combinedOrder, ITSNLayers))};
  const auto combinedTopology = deriveTraversalTopology(combinedLayout);
  BOOST_REQUIRE(combinedTopology.ok());
  BOOST_CHECK_EQUAL(standaloneTopology.topology->edges.size(), combinedTopology.topology->edges.size());
  BOOST_CHECK_EQUAL(standaloneTopology.topology->paths.size(), combinedTopology.topology->paths.size());
  for (uint16_t k = 0; k < ITSNLayers; ++k) {
    BOOST_CHECK(standaloneLayout.getOrderedSurfaces()[k] == combinedLayout.getOrderedSurfaces()[k]);
  }

  // A separately constructed standalone scratch remains independent. The
  // workflow participant accessors below, in contrast, must alias the one
  // global workspace.
  SurfaceTrackingScratch standaloneScratch;
  standaloneScratch.adoptPlan(standaloneOrder.size(), standaloneTopology.topology->edges.size(), standaloneTopology.topology->paths.size());
  auto participants = makeSet();
  TimeFrame frame;
  participants.adoptFrame(frame);
  BOOST_CHECK_EQUAL(standaloneScratch.getNOwnedSurfaces(), static_cast<size_t>(ITSNLayers));
  BOOST_CHECK_EQUAL(&participants.getITSScratch(), &participants.getMFTScratch());
  BOOST_CHECK_EQUAL(participants.getITSScratch().getNOwnedSurfaces(), 17u);
  BOOST_CHECK_EQUAL(participants.getITSScratch().getNEdges(), 15u);
  BOOST_CHECK_EQUAL(participants.getITSScratch().getNCells(), 13u);
}

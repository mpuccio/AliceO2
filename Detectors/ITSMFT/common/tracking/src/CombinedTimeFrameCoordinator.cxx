// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/CombinedTimeFrameCoordinator.h"

#include <stdexcept>
#include <utility>

#include "Framework/Logger.h"
#include "ITSMFTTracking/DetectorLayoutBuilder.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/TimeFrameLoadFailure.h"

namespace o2::itsmft::tracking
{

namespace
{

SurfaceCatalogView combinedCatalogView()
{
  return SurfaceCatalogView{kITSMFTCombinedStaticSurfaceCatalog.data(),
                            static_cast<uint32_t>(kITSMFTCombinedStaticSurfaceCatalog.size())};
}

std::vector<SurfaceId> orderedSurfaceRange(uint16_t first, uint16_t count)
{
  std::vector<SurfaceId> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

SurfaceMask surfaceRangeMask(uint16_t first, uint16_t count)
{
  SurfaceMask result;
  for (uint16_t i = 0; i < count; ++i) {
    result.set(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

// Builds the one shared ITS+MFT DetectorLayout -- both detectors' own
// disjoint subgraphs in a single DetectorLayoutBuilder call, exactly as
// testCombinedStaticSurfaceCatalogTopology.cxx /
// testDetectorTraversalBindingOrchestration.cxx already prove -- from each
// detector's own (single-iteration) TrackingParameters. Called exactly once
// by the constructor: this is the coordinator's one authoritative combined-
// topology construction; ownDetectorPlan() below only ever copies its
// result, never rebuilds it.
DetectorLayoutBuildResult buildCombinedLayout(gsl::span<const SurfaceId> itsSurfaces, const TrackingParameters& itsParams,
                                              gsl::span<const SurfaceId> mftSurfaces, const TrackingParameters& mftParams)
{
  DetectorLayoutSubgraph itsSubgraph;
  itsSubgraph.orderedSurfaces.assign(itsSurfaces.begin(), itsSurfaces.end());
  itsSubgraph.maxHoles = itsParams.MaxHoles;
  itsSubgraph.holeSurfaces = positionalSurfaceMask(itsParams.HoleLayerMask, itsSurfaces, static_cast<uint32_t>(itsSurfaces.size()));
  itsSubgraph.seedingSurfaces = positionalSurfaceMask(itsParams.StartLayerMask, itsSurfaces, static_cast<uint32_t>(itsSurfaces.size()));
  itsSubgraph.policyTag = TransitionPolicyTag::CylinderCylinder;

  DetectorLayoutSubgraph mftSubgraph;
  mftSubgraph.orderedSurfaces.assign(mftSurfaces.begin(), mftSurfaces.end());
  mftSubgraph.maxHoles = mftParams.MaxHoles;
  mftSubgraph.holeSurfaces = positionalSurfaceMask(mftParams.HoleLayerMask, mftSurfaces, static_cast<uint32_t>(mftSurfaces.size()));
  mftSubgraph.seedingSurfaces = positionalSurfaceMask(mftParams.StartLayerMask, mftSurfaces, static_cast<uint32_t>(mftSurfaces.size()));
  mftSubgraph.policyTag = TransitionPolicyTag::DiskDisk;

  DetectorLayoutBuilder builder{combinedCatalogView()};
  builder.addSubgraph(std::move(itsSubgraph));
  builder.addSubgraph(std::move(mftSubgraph));
  return builder.build();
}

// Wraps a *copy* of the coordinator's one authoritative combined
// ITS+MFT DetectorLayout (built exactly once, by buildCombinedLayout(),
// in the constructor -- never rebuilt here) into the DetectorLayoutSet
// `NLayers` will actually adopt: its key's orderedSurfaces is that
// detector's own global-id span, so TrackerTraits<NLayers>::
// initialiseTimeFrame()'s legacy-layer/material binding (step 2.5,
// TrackerTraits.cxx) resolves against the right detector, while the layout
// content itself carries both detectors' topology for the adopted
// DetectorTraversalBinding to scope.
//
// `layouts.push_back(authoritative)` is a plain copy-construction (DetectorLayout/
// SparseTrackingTopology are ordinary copyable value types -- no owning
// resource, no deleted/user-provided special members) of the one already-
// built object, not a second independent construction: DetectorLayoutSet
// has no reference/shared-ownership constructor (it always takes its
// `std::vector<DetectorLayout>` by value), so a passive copy is the only way
// two DetectorLayoutSets can each hold this coordinator's one authoritative
// topology. Both copies are therefore guaranteed byte-identical, including
// every global TransitionId/CellTopologyId's content, by construction --
// there is no code path that could make them diverge.
template <int NLayers>
DetectorLayoutSet ownDetectorPlan(const DetectorLayout& authoritative, gsl::span<const SurfaceId> ownSurfaces,
                                  const TrackingParameters& ownParams, TransitionPolicyTag ownPolicy)
{
  DetectorLayoutConfigurationKey key;
  key.orderedSurfaces.assign(ownSurfaces.begin(), ownSurfaces.end());
  key.policyTag = ownPolicy;
  key.iterations.push_back(DetectorLayoutIterationConfiguration{
    static_cast<uint32_t>(NLayers), ownParams.MaxHoles, ownParams.HoleLayerMask, ownParams.StartLayerMask});
  std::vector<DetectorLayout> layouts;
  layouts.push_back(authoritative);
  return DetectorLayoutSet{std::move(key), combinedCatalogView(), std::move(layouts)};
}

} // namespace

CombinedTimeFrameCoordinator::CombinedTimeFrameCoordinator(std::vector<o2::itsmft::TrackingParameters> itsParams,
                                                           std::vector<o2::itsmft::TrackingParameters> mftParams)
  : mITSParticipant(ParticipantId{0}, itsParams), mMFTParticipant(ParticipantId{1}, mftParams)
{
  if (itsParams.size() != 1 || mftParams.size() != 1) {
    throw std::invalid_argument("CombinedTimeFrameCoordinator requires exactly one TrackingParameters iteration per detector");
  }

  const auto itsSurfaces = orderedSurfaceRange(0, ITSNLayers);
  const auto mftSurfaces = orderedSurfaceRange(ITSNLayers, MFTNLayers);

  // The coordinator's one authoritative combined ITS+MFT topology
  // construction. mITSPlan/mMFTPlan below each get a passive copy of
  // `combinedLayout` (see ownDetectorPlan()'s own doc) -- never a second,
  // independent buildCombinedLayout() call.
  auto combinedBuild = buildCombinedLayout(itsSurfaces, itsParams[0], mftSurfaces, mftParams[0]);
  if (!combinedBuild.ok()) {
    throw std::runtime_error("CombinedTimeFrameCoordinator: failed to build the shared ITS+MFT DetectorLayout");
  }
  const DetectorLayout& combinedLayout = *combinedBuild.layout;

  mITSPlan.emplace(ownDetectorPlan<ITSNLayers>(combinedLayout, itsSurfaces, itsParams[0], TransitionPolicyTag::CylinderCylinder));
  mMFTPlan.emplace(ownDetectorPlan<MFTNLayers>(combinedLayout, mftSurfaces, mftParams[0], TransitionPolicyTag::DiskDisk));

  auto itsBindingResult = DetectorTraversalBinding::build(mITSPlan->getLayoutView(0), o2::detectors::DetID::ITS, ClusterSourceId{0},
                                                          surfaceRangeMask(0, ITSNLayers), itsSurfaces);
  if (!itsBindingResult.ok()) {
    throw std::runtime_error("CombinedTimeFrameCoordinator: failed to build the ITS DetectorTraversalBinding");
  }
  auto mftBindingResult = DetectorTraversalBinding::build(mMFTPlan->getLayoutView(0), o2::detectors::DetID::MFT, ClusterSourceId{1},
                                                          surfaceRangeMask(ITSNLayers, MFTNLayers), mftSurfaces);
  if (!mftBindingResult.ok()) {
    throw std::runtime_error("CombinedTimeFrameCoordinator: failed to build the MFT DetectorTraversalBinding");
  }

  mITSParticipant.adoptDetectorTraversalBinding(std::move(itsBindingResult.binding));
  mMFTParticipant.adoptDetectorTraversalBinding(std::move(mftBindingResult.binding));
  mITSParticipant.adoptDetectorLayoutSet(*mITSPlan);
  mMFTParticipant.adoptDetectorLayoutSet(*mMFTPlan);

  mSchedule = {&mITSParticipant, &mMFTParticipant};
}

void CombinedTimeFrameCoordinator::adoptFrame(TimeFrame& frame)
{
  mFrame = &frame;
  mITSParticipant.adoptFrame(frame);
  mMFTParticipant.adoptFrame(frame);
}

void CombinedTimeFrameCoordinator::setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
{
  if (mFrame != nullptr) {
    mFrame->setMemoryPool(pool);
  }
  mITSParticipant.setMemoryPool(pool);
  mMFTParticipant.setMemoryPool(pool);
}

void CombinedTimeFrameCoordinator::setBz(float bz)
{
  if (mFrame != nullptr) {
    mFrame->setBz(bz);
  }
  mITSParticipant.setBz(bz);
  mMFTParticipant.setBz(bz);
}

void CombinedTimeFrameCoordinator::setNThreads(int n)
{
  mITSParticipant.setNThreads(n);
  mMFTParticipant.setNThreads(n);
}

void CombinedTimeFrameCoordinator::resetCombinedEvent() noexcept
{
  // TrackingEngine::resetEvent() resets every scheduled participant's own
  // scratch/sidecar (eventReset(), schedule order) and then wipes the
  // shared TimeFrame exactly once -- the same sequencing
  // MultiSourceTimeFrameLoader::resetITSAndMFTEvent() applied directly
  // before this slice.
  mEngine.resetEvent(*mFrame, schedule());
  mITSClock.reset();
  mMFTClock.reset();
  mPublicationValid = false;
}

CombinedTimeFrameCoordinator::CombinedTrackingResult CombinedTimeFrameCoordinator::process(
  const ClusterSourceInput& itsSource, const ClusterSourceInput& mftSource, const o2::InteractionRecord& origin)
{
  if (mFrame == nullptr) {
    throw std::logic_error("CombinedTimeFrameCoordinator::process() called before adoptFrame()");
  }

  mPublicationValid = false;
  mITSClock.reset();
  mMFTClock.reset();

  // A prior successful process() call leaves both sidecars sealed
  // (ITSSharedClusterCompatibility::sealFromMarkedTracks(),
  // MFTPublicationCompatibility's per-TF entries): unlike
  // TimeFrame::getCommonTracks()/getTrackClusterIndices(), which
  // commitNormalizedFrame() clears atomically on the *next* load, neither
  // sidecar is cleared by a successful process() itself. Left sealed, the
  // very next TF's first accepted ITS track would fail
  // AcceptedTrackShadowPublisher<ITSNLayers>::publish()'s already-sealed
  // guard and fatal inside TrackerTraits<NLayers>::acceptTracks() ("CommonTrack
  // shadow construction failed"). Clearing both unconditionally at the top
  // of every process() call -- success or failure alike -- keeps every TF
  // starting from the same fresh state regardless of how the previous one
  // ended.
  mITSParticipant.clearPublicationSidecar();
  mMFTParticipant.clearPublicationSidecar();

  const auto loadResult = MultiSourceTimeFrameLoader::loadITSAndMFT(*mFrame, mITSParticipant.getScratch(), mMFTParticipant.getScratch(),
                                                                    itsSource, mftSource, combinedCatalogView(), origin);
  if (!loadResult.ok()) {
    // Reuse isRecoverableLoadError() (TimeFrameLoadFailure.h) rather than a
    // parallel taxonomy, then gate it by the *owning* detector's own
    // DropTFUponFailure -- ClusterSourceId{0}/{1} is loadITSAndMFT()'s own
    // fixed ITS/MFT position contract (MultiSourceTimeFrameLoader.h). An
    // unrecognized/missing source, a structural MultiSourceLoadError, or a
    // recoverable one whose owning detector has DropTFUponFailure=false is
    // always Structural -- never silently reclassified as dropped.
    //
    // This is a *load* failure: the event was never atomically committed,
    // so TrackingEngine::executeEvent() must never be called on it --
    // resetEvent() alone applies the same all-participant/shared-frame
    // reset contract without ever reaching track().
    const bool errorIsRecoverable = isRecoverableLoadError(loadResult.error, loadResult.timingDetail);
    const bool isITS = loadResult.source == ClusterSourceId{0};
    const bool isMFT = loadResult.source == ClusterSourceId{1};
    const bool sourceRecognized = isITS || isMFT;
    const bool dropAllowed = isITS ? mITSParticipant.getDropTFUponFailure() : isMFT ? mMFTParticipant.getDropTFUponFailure()
                                                                                    : false;
    if (!sourceRecognized) {
      LOGP(error, "Combined TF load failure reports an unrecognized source id {}; treating as structural", loadResult.source.value());
    }
    const auto outcome = errorIsRecoverable && sourceRecognized && dropAllowed
                           ? CombinedOutcome::RecoverableDropped
                           : CombinedOutcome::Structural;
    LOGP(error, "Combined TF loading failed (source={}, error={}, recoverable={}, dropAllowed={}): outcome={}",
         loadResult.source.value(), static_cast<int>(loadResult.error), errorIsRecoverable, dropAllowed,
         outcome == CombinedOutcome::RecoverableDropped ? "RecoverableDropped" : "Structural");
    resetCombinedEvent();
    return {outcome, 0, 0};
  }

  mITSParticipant.configureRofTables(itsSource.timing, static_cast<uint32_t>(itsSource.rofs.size()));
  mMFTParticipant.configureRofTables(mftSource.timing, static_cast<uint32_t>(mftSource.rofs.size()));

  // The load has committed: executeEvent() may now run. It executes the
  // explicit [ITS, MFT] schedule's track() calls in that exact order into
  // the shared TimeFrame (accepted CommonTracks therefore append
  // ITS-then-MFT), and on any non-Success outcome or exception already
  // applies the same whole-event reset resetCombinedEvent() above applies
  // for a load failure -- see TrackingEngine::executeEvent()'s own doc.
  const auto eventResult = mEngine.executeEvent(*mFrame, schedule());
  if (eventResult.outcome != ParticipantOutcome::Success) {
    LOGP(error, "Combined TF tracking failed via the delegated engine (outcome={})",
         eventResult.outcome == ParticipantOutcome::RecoverableDropped ? "RecoverableDropped" : "Structural");
    // executeEvent() already reset every participant and wiped the shared
    // TimeFrame; only the coordinator's own publication/timing bridge
    // state remains to invalidate.
    mITSClock.reset();
    mMFTClock.reset();
    mPublicationValid = false;
    const auto outcome = eventResult.outcome == ParticipantOutcome::RecoverableDropped
                           ? CombinedOutcome::RecoverableDropped
                           : CombinedOutcome::Structural;
    return {outcome, 0, 0};
  }

  mITSClock.emplace(mITSParticipant.getScratch().getROFOverlapTableView().getClockLayer());
  mMFTClock.emplace(mMFTParticipant.getScratch().getROFOverlapTableView().getClockLayer());
  mPublicationValid = true;

  CombinedTrackingResult result;
  result.outcome = CombinedOutcome::Success;
  result.nITSTracks = mITSParticipant.getScratch().getNumberOfTracks();
  result.nMFTTracks = mMFTParticipant.getScratch().getNumberOfTracks();
  return result;
}

std::optional<CommonTrackPublicationExport> CombinedTimeFrameCoordinator::getITSPublicationExport() const
{
  if (!mPublicationValid || !mITSClock) {
    return std::nullopt;
  }
  return CommonTrackPublicationExport{o2::detectors::DetID::ITS, ClusterSourceId{0}, *mITSClock,
                                      gsl::span<const SurfaceId>{mITSPlan->getConfigurationKey().orderedSurfaces}};
}

std::optional<CommonTrackPublicationExport> CombinedTimeFrameCoordinator::getMFTPublicationExport() const
{
  if (!mPublicationValid || !mMFTClock) {
    return std::nullopt;
  }
  return CommonTrackPublicationExport{o2::detectors::DetID::MFT, ClusterSourceId{1}, *mMFTClock,
                                      gsl::span<const SurfaceId>{mMFTPlan->getConfigurationKey().orderedSurfaces}};
}

} // namespace o2::itsmft::tracking

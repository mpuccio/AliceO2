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
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/TimeFrameLoadFailure.h"

namespace o2::itsmft::tracking
{

CombinedTimeFrameCoordinator::CombinedTimeFrameCoordinator(std::vector<o2::itsmft::TrackingParameters> itsParams,
                                                           std::vector<o2::itsmft::TrackingParameters> mftParams)
  : mParticipants(std::move(itsParams), std::move(mftParams))
{
}

void CombinedTimeFrameCoordinator::adoptFrame(TimeFrame& frame)
{
  mFrame = &frame;
  mParticipants.adoptFrame(frame);
}

void CombinedTimeFrameCoordinator::setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
{
  if (mFrame != nullptr) {
    mFrame->setMemoryPool(pool);
  }
  mParticipants.setMemoryPool(pool);
}

void CombinedTimeFrameCoordinator::setBz(float bz)
{
  if (mFrame != nullptr) {
    mFrame->setBz(bz);
  }
  mParticipants.setBz(bz);
}

void CombinedTimeFrameCoordinator::setNThreads(int n)
{
  mParticipants.setNThreads(n);
}

void CombinedTimeFrameCoordinator::resetCombinedEvent() noexcept
{
  // TrackingEngine::resetEvent() resets every scheduled participant's own
  // scratch/sidecar (eventReset(), schedule order) and then wipes the
  // shared TimeFrame exactly once -- the same sequencing
  // MultiSourceTimeFrameLoader::resetITSAndMFTEvent() applied directly
  // before this slice.
  mEngine.resetEvent(*mFrame, mParticipants.schedule());
  mParticipants.invalidatePublication();
}

CombinedTimeFrameCoordinator::CombinedTrackingResult CombinedTimeFrameCoordinator::process(
  const ClusterSourceInput& itsSource, const ClusterSourceInput& mftSource, const o2::InteractionRecord& origin)
{
  if (mFrame == nullptr) {
    throw std::logic_error("CombinedTimeFrameCoordinator::process() called before adoptFrame()");
  }

  mParticipants.invalidatePublication();

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
  mParticipants.clearPublicationSidecars();

  // M2b/M2c: the fixed ITS=0/MFT=1 source contract lives only in
  // mParticipants.validateSources() now, not in the generic loadEvent()
  // transaction below and not in this coordinator itself. A mismatch is
  // synthesized as the exact same LoadSourcesResult loadITSAndMFT()'s
  // equivalent guard used to return, so the generic error-handling block
  // right below classifies it identically.
  LoadSourcesResult loadResult;
  if (const auto rejected = mParticipants.validateSources(itsSource, mftSource)) {
    loadResult = *rejected;
  } else {
    const auto bindings = mParticipants.loadBindings(itsSource, mftSource);
    loadResult = MultiSourceTimeFrameLoader::loadEvent(*mFrame, gsl::span<const MultiSourceTimeFrameLoader::AtomicLoadBinding>{bindings},
                                                       mParticipants.catalogView(), origin);
  }
  if (!loadResult.ok()) {
    // Reuse isRecoverableLoadError() (TimeFrameLoadFailure.h) rather than a
    // parallel taxonomy, then gate it by the *owning* detector's own
    // DropTFUponFailure -- mParticipants.dropTFUponFailureFor() carries the
    // fixed ITS/MFT source-position contract this coordinator's own guard
    // used to enforce directly.
    //
    // This is a *load* failure: the event was never atomically committed,
    // so TrackingEngine::executeEvent() must never be called on it --
    // resetEvent() alone applies the same all-participant/shared-frame
    // reset contract without ever reaching track().
    const bool errorIsRecoverable = isRecoverableLoadError(loadResult.error, loadResult.timingDetail);
    const auto dropAllowed = mParticipants.dropTFUponFailureFor(loadResult.source);
    const bool sourceRecognized = dropAllowed.has_value();
    if (!sourceRecognized) {
      LOGP(error, "Combined TF load failure reports an unrecognized source id {}; treating as structural", loadResult.source.value());
    }
    const auto outcome = errorIsRecoverable && sourceRecognized && dropAllowed.value_or(false)
                           ? CombinedOutcome::RecoverableDropped
                           : CombinedOutcome::Structural;
    LOGP(error, "Combined TF loading failed (source={}, error={}, recoverable={}, dropAllowed={}): outcome={}",
         loadResult.source.value(), static_cast<int>(loadResult.error), errorIsRecoverable, dropAllowed.value_or(false),
         outcome == CombinedOutcome::RecoverableDropped ? "RecoverableDropped" : "Structural");
    resetCombinedEvent();
    return {outcome, 0, 0};
  }

  mParticipants.configureRofTables(itsSource, mftSource);

  // The load has committed: executeEvent() may now run. It executes the
  // explicit [ITS, MFT] schedule's track() calls in that exact order into
  // the shared TimeFrame (accepted CommonTracks therefore append
  // ITS-then-MFT), and on any non-Success outcome or exception already
  // applies the same whole-event reset resetCombinedEvent() above applies
  // for a load failure -- see TrackingEngine::executeEvent()'s own doc.
  const auto eventResult = mEngine.executeEvent(*mFrame, mParticipants.schedule());
  if (eventResult.outcome != ParticipantOutcome::Success) {
    LOGP(error, "Combined TF tracking failed via the delegated engine (outcome={})",
         eventResult.outcome == ParticipantOutcome::RecoverableDropped ? "RecoverableDropped" : "Structural");
    // executeEvent() already reset every participant and wiped the shared
    // TimeFrame; only the publication/timing bridge state remains to
    // invalidate.
    mParticipants.invalidatePublication();
    const auto outcome = eventResult.outcome == ParticipantOutcome::RecoverableDropped
                           ? CombinedOutcome::RecoverableDropped
                           : CombinedOutcome::Structural;
    return {outcome, 0, 0};
  }

  mParticipants.markPublicationValid();

  CombinedTrackingResult result;
  result.outcome = CombinedOutcome::Success;
  result.nITSTracks = mParticipants.getITSScratch().getNumberOfTracks();
  result.nMFTTracks = mParticipants.getMFTScratch().getNumberOfTracks();
  return result;
}

} // namespace o2::itsmft::tracking

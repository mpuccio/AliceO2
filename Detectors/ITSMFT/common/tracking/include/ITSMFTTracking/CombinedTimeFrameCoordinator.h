// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 C3/C4, E0/M2a/M2b/M2c (GenericTrackingEngineMigration.md; ADR
// 0007): a host-only combined disconnected-tracking coordinator, reachable
// in production from the opt-in Gate 4 C4 combined DPL workflow
// (Detectors/ITSMFT/common/workflow-combined-ca/). Executes ITS and MFT
// common CA against one shared TimeFrame.
//
// M2a wired the accepted TrackingEngine (TrackingEngine.h) into this
// coordinator's tracking phase: process() delegates tracking execution to
// TrackingEngine::executeEvent() over an explicit ordered schedule instead
// of hand-unrolling two clustersToTracks() calls. M2b replaced the atomic
// load path with MultiSourceTimeFrameLoader::loadEvent() (the participant-
// count-generic atomic transaction) over an explicit AtomicLoadBinding list.
//
// M2c moves every remaining ITS/MFT construction/configuration fact out of
// this coordinator and into one explicit application-layer factory/set,
// ITSMFTLegacyParticipantSet (ITSMFTLegacyParticipantSet.h): the two
// LegacyCATrackingParticipant instances (each already owning its own
// Tracker<NLayers>/TrackerTraits<NLayers>/LegacyTrackerScratch<NLayers>/
// DetectorTraversalBinding/detector compatibility sidecar), the static
// combined 17-surface ITS+MFT plan/binding configuration
// (StaticDetectorCatalogs.h; global ids ITS 0..6, MFT 7..16), the explicit
// [ITS, MFT] tracking schedule, the fixed ITS=0/MFT=1 atomic load-binding
// contract, and the temporary detector-specific publication/timing bridge
// (ClockTimingPublicationView, publication validity) all now live inside
// that one set. This coordinator holds only the set and the TrackingEngine,
// and delegates: construction to the set's constructor, atomic loading
// through the set's validateSources()/loadBindings() plus
// MultiSourceTimeFrameLoader::loadEvent(), and tracking/reset through
// TrackingEngine over the set's schedule().
//
// TrackingEngine::executeEvent() is only ever called once the atomic load
// has already committed. Any non-success -- a load failure, a non-Success
// tracking outcome, or an exception -- is a whole combined-TF failure: on a
// load failure, TrackingEngine::resetEvent() is called directly
// (executeEvent() is never reached on a partially loaded event); on a
// tracking-phase failure, executeEvent() itself already applies that same
// reset internally before returning. Both paths leave both legs' scratches,
// the shared TimeFrame's CommonTracks/references, and both detector
// compatibility sidecars empty. The failure is classified
// RecoverableDropped or Structural (see CombinedOutcome's own doc); the
// coordinator never fatals or throws past process()'s own boundary for a
// classified failure -- that decision belongs to the caller (e.g.
// skip-and-continue vs. fatal in the DPL workflow).
//
// This coordinator is now only a temporary compatibility wrapper over
// TrackingEngine plus the one ITSMFTLegacyParticipantSet (ADR 0007 decision
// 4): it contains no direct static catalog construction, DetectorLayoutSet,
// DetectorTraversalBinding, Tracker, TrackerTraits, LegacyTrackerScratch,
// source 0/1, ITS/MFT layer-count, or timing/publication-bridge members of
// its own. M3 deletes this coordinator once the C4 workflow constructs the
// engine and the participant set directly.

#ifndef ALICEO2_ITSMFT_TRACKING_COMBINEDTIMEFRAMECOORDINATOR_H_
#define ALICEO2_ITSMFT_TRACKING_COMBINEDTIMEFRAMECOORDINATOR_H_

#ifndef GPUCA_GPUCODE

#include <memory>
#include <optional>
#include <vector>

#include "CommonDataFormat/InteractionRecord.h"
#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/ITSMFTLegacyParticipantSet.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingEngine.h"
// Reused only for its detector-neutral CommonTrackPublicationExport type
// (source, detector, copied ClockTimingPublicationView, ordered-surface
// span); this coordinator never calls ITSMFTTrackingInterface's own
// source-0-only getCommonTrackPublicationExport() getter, since that always
// hardcodes ClusterSourceId{0}.
#include "ITSMFTTracking/TrackingInterface.h"

namespace o2::itsmft::tracking
{

class CombinedTimeFrameCoordinator
{
 public:
  // Reuses TrackingOutcome's own vocabulary (CATracker.h) rather than a
  // parallel taxonomy: RecoverableDropped means a per-TF data failure this
  // detector's own DropTFUponFailure opted to drop; every other non-success
  // -- an unrecognized load source, a structural MultiSourceLoadError, a
  // recoverable load error whose owning detector has DropTFUponFailure=false,
  // or any thrown exception (TraversalException, unclassified, or a
  // resource-exhaustion exception that already failed its own
  // DropTFUponFailure gate inside clustersToTracks()) -- is Structural.
  enum class CombinedOutcome : uint8_t {
    Success,
    RecoverableDropped,
    Structural
  };

  struct CombinedTrackingResult {
    CombinedOutcome outcome{CombinedOutcome::Structural};
    size_t nITSTracks{0};
    size_t nMFTTracks{0};
  };

  // `itsParams`/`mftParams` must each carry exactly one iteration -- the
  // only shape this slice's single shared combined layout can represent
  // (matching the only tracking-mode both detectors' common CA currently
  // supports end to end, Sync; see D009/D010, AgentCoordination.md). Throws
  // std::invalid_argument otherwise.
  CombinedTimeFrameCoordinator(std::vector<o2::itsmft::TrackingParameters> itsParams,
                               std::vector<o2::itsmft::TrackingParameters> mftParams);

  // Every participant's Tracker/TrackerTraits binds to a sibling member's
  // address at construction time; relocating this object would silently
  // dangle every one of those bound pointers. Never copyable or movable --
  // construct once, use by reference.
  CombinedTimeFrameCoordinator(const CombinedTimeFrameCoordinator&) = delete;
  CombinedTimeFrameCoordinator& operator=(const CombinedTimeFrameCoordinator&) = delete;
  CombinedTimeFrameCoordinator(CombinedTimeFrameCoordinator&&) = delete;
  CombinedTimeFrameCoordinator& operator=(CombinedTimeFrameCoordinator&&) = delete;

  // Binds the one shared TimeFrame both legs adopt. `frame` must outlive
  // every subsequent process() call.
  void adoptFrame(TimeFrame& frame);
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool);
  void setBz(float bz);
  void setNThreads(int n);

  // Phase 2+3 for both detectors against one shared TimeFrame: load via
  // MultiSourceTimeFrameLoader::loadEvent() (M2b's participant-count-generic
  // atomic transaction), over an explicit [ITS, MFT] AtomicLoadBinding list
  // built by mParticipants.loadBindings(); mParticipants.validateSources(),
  // not loadEvent() itself, enforces source 0 == ITS, source 1 == MFT --
  // then, once that load has committed, TrackingEngine::executeEvent() runs
  // the explicit [ITS, MFT] schedule: ITS's track() followed by MFT's,
  // serially, into the shared TimeFrame, so accepted CommonTracks append
  // ITS-then-MFT.
  //
  // Any non-success -- a load failure (LoadSourcesResult::ok() == false), a
  // non-Success ParticipantOutcome from either leg, or any exception from
  // either leg's track() -- is a whole combined-TF failure: on a load
  // failure this calls TrackingEngine::resetEvent() directly (executeEvent()
  // is never called on a partially loaded event); on a tracking-phase
  // failure executeEvent() itself already resets both legs and wipes the
  // shared TimeFrame before returning. Either way this returns
  // {RecoverableDropped or Structural, 0, 0} (see CombinedOutcome's own doc
  // for the exact classification). adoptFrame() must have been called
  // first.
  CombinedTrackingResult process(const ClusterSourceInput& itsSource,
                                 const ClusterSourceInput& mftSource,
                                 const o2::InteractionRecord& origin);

  const LegacyTrackerScratchITS& getITSScratch() const noexcept { return mParticipants.getITSScratch(); }
  const LegacyTrackerScratchMFT& getMFTScratch() const noexcept { return mParticipants.getMFTScratch(); }
  const TimeFrame* getFrame() const noexcept { return mFrame; }

  // Detector-local compatibility sidecars, populated fresh by a successful
  // process() call and cleared by resetCombinedEvent() on any non-success.
  // A DPL caller needs these to drive stageITSCommonTrackOutput()/
  // stageMFTCommonTrackOutput() (CommonTrackOutputAdapter.h), which take
  // them by const reference.
  const ITSSharedClusterCompatibility& getITSSharedClusterCompatibility() const noexcept { return mParticipants.getITSSharedClusterCompatibility(); }
  const MFTPublicationCompatibility& getMFTPublicationCompatibility() const noexcept { return mParticipants.getMFTPublicationCompatibility(); }

  // Each detector's ordered global-SurfaceId span (ITS: SurfaceId{0..6},
  // MFT: SurfaceId{7..16} -- the same fixed combined-catalog offset
  // ITSMFTLegacyParticipantSet's constructor uses). Unlike
  // getITSPublicationExport()/getMFTPublicationExport(), these are valid
  // immediately after construction and never invalidated by
  // process()/resetCombinedEvent(): the participant set's plans are
  // unconditionally engaged by its constructor and never rebuilt. A DPL
  // caller uses these to build each detector's own
  // ClusterSourceInput::layerToSurface before the first process() call,
  // instead of re-deriving the ITS/MFT global-id offset by hand.
  gsl::span<const SurfaceId> getITSOrderedSurfaces() const noexcept { return mParticipants.getITSOrderedSurfaces(); }
  gsl::span<const SurfaceId> getMFTOrderedSurfaces() const noexcept { return mParticipants.getMFTOrderedSurfaces(); }

  // Exposed for testing only: the DetectorLayoutView each detector's
  // DetectorTraversalBinding was built from. Both are passive copies of the
  // one combined ITS+MFT DetectorLayout ITSMFTLegacyParticipantSet builds
  // exactly once in its constructor (see buildCombinedLayout()/
  // ownDetectorPlan() in ITSMFTLegacyParticipantSet.cxx) -- a test can use
  // these to assert both copies carry byte-identical topology content,
  // including global TransitionId/CellTopologyId identity, and that neither
  // ever diverges from the other.
  DetectorLayoutView getITSLayoutView() const noexcept { return mParticipants.getITSLayoutView(); }
  DetectorLayoutView getMFTLayoutView() const noexcept { return mParticipants.getMFTLayoutView(); }

  // Immutable per-detector publication exports, populated only by a
  // successful process() call and invalidated by any subsequent reset
  // (failure or exception). ITS is always ClusterSourceId{0}, MFT always
  // ClusterSourceId{1} -- fixed positions, never the caller's choice, the
  // same contract ITSMFTLegacyParticipantSet::validateSources() enforces
  // before this coordinator ever calls
  // MultiSourceTimeFrameLoader::loadEvent().
  std::optional<CommonTrackPublicationExport> getITSPublicationExport() const { return mParticipants.getITSPublicationExport(); }
  std::optional<CommonTrackPublicationExport> getMFTPublicationExport() const { return mParticipants.getMFTPublicationExport(); }

 private:
  void resetCombinedEvent() noexcept;

  TimeFrame* mFrame = nullptr;

  // The one explicit ITS/MFT application-layer participant set (M2c): owns
  // every static plan/binding/schedule/publication-bridge fact this
  // coordinator used to hold directly (see this header's own file-level
  // doc). Declared before mEngine only for readability; mEngine holds no
  // pointer into it.
  ITSMFTLegacyParticipantSet mParticipants;

  TrackingEngine mEngine;
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_COMBINEDTIMEFRAMECOORDINATOR_H_ */

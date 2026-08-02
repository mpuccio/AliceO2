// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 C3/C4, E0/M2a (GenericTrackingEngineMigration.md; ADR 0007): a
// host-only combined disconnected-tracking coordinator, reachable in
// production from the opt-in Gate 4 C4 combined DPL workflow
// (Detectors/ITSMFT/common/workflow-combined-ca/). Executes ITS and MFT
// common CA against one shared TimeFrame, using the static combined
// 17-surface ITS+MFT catalog (StaticDetectorCatalogs.h; global ids ITS
// 0..6, MFT 7..16), one shared DetectorLayout built from that catalog's two
// disjoint per-detector subgraphs, and one DetectorTraversalBinding per
// detector (D010, AgentCoordination.md) scoping each leg's Tracker/
// TrackerTraits to its own owned surfaces/transitions/cells/road starts.
//
// M2a wires the accepted TrackingEngine (TrackingEngine.h) into this
// coordinator's tracking phase: the two detector legs are now
// LegacyCATrackingParticipant<ITSNLayers>/<MFTNLayers>
// (LegacyCATrackingParticipant.h), and process() delegates their execution
// to TrackingEngine::executeEvent() over an explicit [ITS, MFT] schedule
// instead of hand-unrolling two clustersToTracks() calls. The atomic
// two-source loading boundary is untouched: MultiSourceTimeFrameLoader::
// loadITSAndMFT() remains the sole load path, called directly by this
// coordinator exactly as before, and TrackingEngine::executeEvent() is only
// ever called once that load has already committed. Any non-success --
// a load failure, a non-Success tracking outcome, or an exception -- is a
// whole combined-TF failure: on a load failure, TrackingEngine::
// resetEvent() is called directly (executeEvent() is never reached on a
// partially loaded event); on a tracking-phase failure, executeEvent()
// itself already applies that same reset internally before returning. Both
// paths leave both legs' scratches, the shared TimeFrame's CommonTracks/
// references, and both detector compatibility sidecars empty. The failure
// is classified RecoverableDropped or Structural (see CombinedOutcome's own
// doc); the coordinator never fatals or throws past process()'s own
// boundary for a classified failure -- that decision belongs to the caller
// (e.g. skip-and-continue vs. fatal in the DPL workflow).
//
// This coordinator no longer directly owns Tracker<NLayers>,
// TrackerTraits<NLayers>, or LegacyTrackerScratch<NLayers>: those, and each
// leg's own DetectorTraversalBinding/TrackingParameters/detector
// compatibility sidecar, live inside its two LegacyCATrackingParticipant
// members instead. This coordinator still owns the TrackingEngine, the two
// participants, the static combined plan data (DetectorLayoutSet) both
// legs' DetectorTraversalBindings and Trackers reference, the explicit
// [ITS, MFT] schedule, and the source-specific publication/timing bridge
// state (ClockTimingPublicationView, publication validity) that builds each
// leg's richer CommonTrackPublicationExport. It remains the temporary
// two-detector adapter (ADR 0007 decision 4) and still retains the two-
// source loadITSAndMFT() call: participant-count-agnostic loading is out of
// scope for this slice (M2's later generalization).

#ifndef ALICEO2_ITSMFT_TRACKING_COMBINEDTIMEFRAMECOORDINATOR_H_
#define ALICEO2_ITSMFT_TRACKING_COMBINEDTIMEFRAMECOORDINATOR_H_

#ifndef GPUCA_GPUCODE

#include <array>
#include <memory>
#include <optional>
#include <vector>

#include "CommonDataFormat/InteractionRecord.h"
#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/DetectorLayoutSet.h"
#include "ITSMFTTracking/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/LegacyCATrackingParticipant.h"
#include "ITSMFTTracking/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingEngine.h"
#include "ITSMFTTracking/TrackingParticipant.h"
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
  // MultiSourceTimeFrameLoader::loadITSAndMFT() (source 0 == ITS, source 1
  // == MFT, per that loader's own fixed-position contract) -- unchanged,
  // the sole atomic load path -- then, once that load has committed,
  // TrackingEngine::executeEvent() runs the explicit [ITS, MFT] schedule:
  // ITS's track() followed by MFT's, serially, into the shared TimeFrame,
  // so accepted CommonTracks append ITS-then-MFT.
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

  const LegacyTrackerScratchITS& getITSScratch() const noexcept { return mITSParticipant.getScratch(); }
  const LegacyTrackerScratchMFT& getMFTScratch() const noexcept { return mMFTParticipant.getScratch(); }
  const TimeFrame* getFrame() const noexcept { return mFrame; }

  // Detector-local compatibility sidecars, populated fresh by a successful
  // process() call and cleared by resetCombinedEvent() on any non-success.
  // A DPL caller needs these to drive stageITSCommonTrackOutput()/
  // stageMFTCommonTrackOutput() (CommonTrackOutputAdapter.h), which take
  // them by const reference.
  const ITSSharedClusterCompatibility& getITSSharedClusterCompatibility() const noexcept { return *mITSParticipant.getITSSharedClusterCompatibility(); }
  const MFTPublicationCompatibility& getMFTPublicationCompatibility() const noexcept { return *mMFTParticipant.getMFTPublicationCompatibility(); }

  // Each detector's ordered global-SurfaceId span (ITS: SurfaceId{0..6},
  // MFT: SurfaceId{7..16} -- the same fixed combined-catalog offset
  // buildCombinedLayout() uses in the constructor, see
  // CombinedTimeFrameCoordinator.cxx). Unlike getITSPublicationExport()/
  // getMFTPublicationExport(), these are valid immediately after
  // construction and never invalidated by process()/resetCombinedEvent():
  // mITSPlan/mMFTPlan are unconditionally engaged by the constructor and
  // never rebuilt. A DPL caller uses these to build each detector's own
  // ClusterSourceInput::layerToSurface before the first process() call,
  // instead of re-deriving the ITS/MFT global-id offset by hand.
  gsl::span<const SurfaceId> getITSOrderedSurfaces() const noexcept { return mITSPlan->getConfigurationKey().orderedSurfaces; }
  gsl::span<const SurfaceId> getMFTOrderedSurfaces() const noexcept { return mMFTPlan->getConfigurationKey().orderedSurfaces; }

  // Exposed for testing only: the DetectorLayoutView each detector's
  // DetectorTraversalBinding was built from. Both are passive copies of the
  // one combined ITS+MFT DetectorLayout this coordinator builds exactly once
  // in its constructor (see buildCombinedLayout()/ownDetectorPlan() in
  // CombinedTimeFrameCoordinator.cxx) -- a test can use these to assert both
  // copies carry byte-identical topology content, including global
  // TransitionId/CellTopologyId identity, and that neither ever diverges
  // from the other.
  DetectorLayoutView getITSLayoutView() const noexcept { return mITSPlan->getLayoutView(0); }
  DetectorLayoutView getMFTLayoutView() const noexcept { return mMFTPlan->getLayoutView(0); }

  // Immutable per-detector publication exports, populated only by a
  // successful process() call and invalidated by any subsequent reset
  // (failure or exception). ITS is always ClusterSourceId{0}, MFT always
  // ClusterSourceId{1} -- fixed positions, never the caller's choice, same
  // contract as MultiSourceTimeFrameLoader::loadITSAndMFT().
  std::optional<CommonTrackPublicationExport> getITSPublicationExport() const;
  std::optional<CommonTrackPublicationExport> getMFTPublicationExport() const;

 private:
  void resetCombinedEvent() noexcept;
  gsl::span<TrackingParticipant* const> schedule() noexcept { return mSchedule; }

  std::optional<DetectorLayoutSet> mITSPlan;
  std::optional<DetectorLayoutSet> mMFTPlan;

  TimeFrame* mFrame = nullptr;

  // The two temporary legacy detector legs -- own everything ADR 0007
  // still classifies as temporary per-detector state (see this header's
  // own file-level doc). Declared before mEngine/mSchedule only for
  // readability; mEngine holds no pointer into either.
  LegacyCATrackingParticipantITS mITSParticipant;
  LegacyCATrackingParticipantMFT mMFTParticipant;

  TrackingEngine mEngine;
  // Explicit ITS-then-MFT execution order (ADR 0007 decision 6), supplied
  // as plan data to every TrackingEngine::executeEvent()/resetEvent() call
  // -- never re-derived from participant declaration order or any static
  // catalog concatenation. Populated once in the constructor from
  // &mITSParticipant/&mMFTParticipant, both of which are stable for this
  // non-relocatable object's lifetime.
  std::array<TrackingParticipant*, 2> mSchedule{};

  std::optional<ClockTimingPublicationView> mITSClock;
  std::optional<ClockTimingPublicationView> mMFTClock;
  bool mPublicationValid = false;
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_COMBINEDTIMEFRAMECOORDINATOR_H_ */

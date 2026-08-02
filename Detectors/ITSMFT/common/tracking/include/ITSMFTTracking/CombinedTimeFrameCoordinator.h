// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 C3: a plain host-only combined disconnected-tracking coordinator.
// Test/host-core slice only -- not reachable from any workflow. Executes ITS
// and MFT common CA against one shared TimeFrame, using the static combined
// 17-surface ITS+MFT catalog (StaticDetectorCatalogs.h; global ids ITS
// 0..6, MFT 7..16), one shared DetectorLayout built from that catalog's two
// disjoint per-detector subgraphs, and one DetectorTraversalBinding per
// detector (D010, AgentCoordination.md) scoping each Tracker/TrackerTraits
// to its own owned surfaces/transitions/cells/road starts. Loading goes
// through MultiSourceTimeFrameLoader::loadITSAndMFT() only; tracking runs
// serially, ITS then MFT, into the same shared TimeFrame, so accepted
// CommonTracks append in that same deterministic order (see
// AcceptedTrackShadowPublisher.h). Any load failure, non-Success tracking
// outcome, or exception from either leg is a whole combined-TF failure:
// MultiSourceTimeFrameLoader::resetITSAndMFTEvent() is called exactly once,
// and both detector compatibility sidecars are cleared, leaving both
// scratches, CommonTracks, references, and sidecars empty.

#ifndef ALICEO2_ITSMFT_TRACKING_COMBINEDTIMEFRAMECOORDINATOR_H_
#define ALICEO2_ITSMFT_TRACKING_COMBINEDTIMEFRAMECOORDINATOR_H_

#ifndef GPUCA_GPUCODE

#include <memory>
#include <optional>
#include <vector>

#include <oneapi/tbb/task_arena.h>

#include "CommonDataFormat/InteractionRecord.h"
#include "ITSMFTTracking/CATracker.h"
#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/DetectorLayoutSet.h"
#include "ITSMFTTracking/DetectorTraversalBinding.h"
#include "ITSMFTTracking/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/LegacyTrackerScratch.h"
#include "ITSMFTTracking/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
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
  enum class CombinedOutcome : uint8_t {
    Success,
    Failure
  };

  struct CombinedTrackingResult {
    CombinedOutcome outcome{CombinedOutcome::Failure};
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

  // Every Tracker/TrackerTraits member below binds to a sibling member's
  // address at construction time (adoptScratch()/adoptDetectorLayoutSet()/
  // adoptITSSharedClusterCompatibility()/adoptMFTPublicationCompatibility());
  // relocating this object would silently dangle every one of those bound
  // pointers. Never copyable or movable -- construct once, use by reference.
  CombinedTimeFrameCoordinator(const CombinedTimeFrameCoordinator&) = delete;
  CombinedTimeFrameCoordinator& operator=(const CombinedTimeFrameCoordinator&) = delete;
  CombinedTimeFrameCoordinator(CombinedTimeFrameCoordinator&&) = delete;
  CombinedTimeFrameCoordinator& operator=(CombinedTimeFrameCoordinator&&) = delete;

  // Binds the one shared TimeFrame both trackers adopt. `frame` must
  // outlive every subsequent process() call.
  void adoptFrame(TimeFrame& frame);
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool);
  void setBz(float bz);
  void setNThreads(int n);

  // Phase 2+3 for both detectors against one shared TimeFrame: load via
  // MultiSourceTimeFrameLoader::loadITSAndMFT() (source 0 == ITS, source 1
  // == MFT, per that loader's own fixed-position contract), then run ITS's
  // Tracker<ITSNLayers>::clustersToTracks() followed by MFT's
  // Tracker<MFTNLayers>::clustersToTracks() -- serially, ITS first -- into
  // the shared TimeFrame, so accepted CommonTracks append ITS-then-MFT.
  //
  // Any load failure (LoadSourcesResult::ok() == false), any non-Success
  // TrackingOutcome from either tracker, or any exception from either
  // tracker is a whole combined-TF failure: resetITSAndMFTEvent() runs
  // exactly once, both compatibility sidecars are cleared, and this
  // returns {Failure, 0, 0}. adoptFrame() must have been called first.
  CombinedTrackingResult process(const ClusterSourceInput& itsSource,
                                 const ClusterSourceInput& mftSource,
                                 const o2::InteractionRecord& origin);

  const LegacyTrackerScratchITS& getITSScratch() const noexcept { return mITSScratch; }
  const LegacyTrackerScratchMFT& getMFTScratch() const noexcept { return mMFTScratch; }
  const TimeFrame* getFrame() const noexcept { return mFrame; }

  // Immutable per-detector publication exports, populated only by a
  // successful process() call and invalidated by any subsequent reset
  // (failure or exception). ITS is always ClusterSourceId{0}, MFT always
  // ClusterSourceId{1} -- fixed positions, never the caller's choice, same
  // contract as MultiSourceTimeFrameLoader::loadITSAndMFT().
  std::optional<CommonTrackPublicationExport> getITSPublicationExport() const;
  std::optional<CommonTrackPublicationExport> getMFTPublicationExport() const;

 private:
  void resetCombinedEvent() noexcept;

  std::vector<o2::itsmft::TrackingParameters> mITSParams;
  std::vector<o2::itsmft::TrackingParameters> mMFTParams;

  std::optional<DetectorLayoutSet> mITSPlan;
  std::optional<DetectorLayoutSet> mMFTPlan;
  std::unique_ptr<DetectorTraversalBinding> mITSBinding;
  std::unique_ptr<DetectorTraversalBinding> mMFTBinding;

  TimeFrame* mFrame = nullptr;
  LegacyTrackerScratchITS mITSScratch;
  LegacyTrackerScratchMFT mMFTScratch;
  ITSSharedClusterCompatibility mITSCompatibility;
  MFTPublicationCompatibility mMFTCompatibility;

  TrackerTraits<ITSNLayers> mITSTraits;
  TrackerTraits<o2::mft::constants::mft::LayersNumber> mMFTTraits;
  Tracker<ITSNLayers> mITSTracker;
  Tracker<o2::mft::constants::mft::LayersNumber> mMFTTracker;

  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  std::shared_ptr<tbb::task_arena> mITSArena;
  std::shared_ptr<tbb::task_arena> mMFTArena;

  std::optional<ClockTimingPublicationView> mITSClock;
  std::optional<ClockTimingPublicationView> mMFTClock;
  bool mPublicationValid = false;
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_COMBINEDTIMEFRAMECOORDINATOR_H_ */

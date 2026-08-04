// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// M2c (GenericTrackingEngineMigration.md; ADR 0007): the one explicitly
// ITS/MFT-named application-layer factory/set. It owns the two concrete
// LegacyCATrackingParticipant instances (each already owning its own
// Tracker<NLayers>/TrackerTraits<NLayers>/LegacyTrackerScratch<NLayers>/
// DetectorTraversalBinding/detector compatibility sidecar), the static
// combined ITS+MFT plan/binding configuration (StaticDetectorCatalogs.h;
// global ids ITS 0..6, MFT 7..16), the explicit [ITS, MFT] tracking
// schedule, the fixed ITS=0/MFT=1 atomic load-binding contract, and the
// temporary detector-specific publication/timing bridge
// (ClockTimingPublicationView, publication validity).
//
// This is the sole owner of the current ITS/MFT source/layout facts: its
// caller -- the combined DPL task (CombinedCATrackerSpec.cxx) directly,
// as of M3 -- never needs to know source 0/1, ITS/MFT layer counts, or any
// static catalog/binding/tracker/scratch type to delegate construction,
// atomic loading, tracking, and publication export to this set instead.
// None of that knowledge is pushed down into
// TrackingEngine.h/TrackingParticipant.h/TimeFrame.h/CommonTrack.h/
// MultiSourceTimeFrameLoader.h, which stay exactly as detector-neutral as
// ADR 0007 requires.

#ifndef ALICEO2_ITSMFT_TRACKING_ITSMFTLEGACYPARTICIPANTSET_H_
#define ALICEO2_ITSMFT_TRACKING_ITSMFTLEGACYPARTICIPANTSET_H_

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
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingParticipant.h"
// Reused only for its detector-neutral CommonTrackPublicationExport type
// (source, detector, copied ClockTimingPublicationView, ordered-surface
// span); this set never calls ITSMFTTrackingInterface's own source-0-only
// getCommonTrackPublicationExport() getter, since that always hardcodes
// ClusterSourceId{0}.
#include "ITSMFTTracking/TrackingInterface.h"

namespace o2::itsmft::tracking
{

// The one explicit ITS/MFT-named application-layer factory/set (M2c). Owns
// every current ITS/MFT construction/configuration fact; exposes only the
// minimal operations its caller (the combined DPL task, M3) needs: the
// tracking schedule, atomic load bindings, a reset-compatible participant
// list, and ITS/MFT adapter publication exports.
class ITSMFTLegacyParticipantSet
{
 public:
  // `itsParams`/`mftParams` must each carry exactly one iteration -- the
  // only shape this set's single shared combined layout can represent
  // (matching the only tracking-mode both detectors' common CA currently
  // supports end to end, Sync; see D009/D010, AgentCoordination.md). Throws
  // std::invalid_argument otherwise.
  ITSMFTLegacyParticipantSet(std::vector<o2::itsmft::TrackingParameters> itsParams,
                             std::vector<o2::itsmft::TrackingParameters> mftParams);

  // Every participant's Tracker/TrackerTraits binds to a sibling member's
  // address at construction time; relocating this object would silently
  // dangle every one of those bound pointers. Never copyable or movable --
  // construct once, use by reference (the same non-relocatable contract
  // this set's own predecessor owner documented before M3).
  ITSMFTLegacyParticipantSet(const ITSMFTLegacyParticipantSet&) = delete;
  ITSMFTLegacyParticipantSet& operator=(const ITSMFTLegacyParticipantSet&) = delete;
  ITSMFTLegacyParticipantSet(ITSMFTLegacyParticipantSet&&) = delete;
  ITSMFTLegacyParticipantSet& operator=(ITSMFTLegacyParticipantSet&&) = delete;

  // --- Coordinator setup, forwarded to both participants (and, for
  // setMemoryPool()/setBz(), the shared TimeFrame is set by the coordinator
  // itself since this set never owns it). ---
  void adoptFrame(TimeFrame& frame);
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool);
  void setBz(float bz);
  void setNThreads(int n);

  // --- Tracking schedule (ADR 0007 decision 6): the explicit ITS-then-MFT
  // execution order, supplied as plan data to every TrackingEngine::
  // executeEvent()/resetEvent() call -- never re-derived from participant
  // declaration order or any static catalog concatenation. ---
  gsl::span<TrackingParticipant* const> schedule() noexcept { return mSchedule; }

  // --- Atomic load bindings (M2b/M2c): the fixed ITS=0/MFT=1 source
  // contract lives only here now. validateSources() returns the synthesized
  // LoadSourcesResult a mismatch produces (matching loadITSAndMFT()'s own
  // equivalent guard) so the coordinator's generic error-handling can
  // classify it identically; nullopt means the sources satisfy the
  // contract and loadBindings() may be used to build the actual
  // MultiSourceTimeFrameLoader::loadEvent() transaction. ---
  std::optional<LoadSourcesResult> validateSources(const ClusterSourceInput& itsSource,
                                                   const ClusterSourceInput& mftSource) const noexcept;
  std::array<MultiSourceTimeFrameLoader::AtomicLoadBinding, 2> loadBindings(const ClusterSourceInput& itsSource,
                                                                            const ClusterSourceInput& mftSource) noexcept;
  SurfaceCatalogView catalogView() const noexcept;

  // The owning detector's own DropTFUponFailure for a load failure's
  // `source`; nullopt if `source` is neither the fixed ITS nor MFT source
  // id -- the coordinator's own "unrecognized source" branch.
  std::optional<bool> dropTFUponFailureFor(ClusterSourceId source) const noexcept;

  void configureRofTables(const ClusterSourceInput& itsSource, const ClusterSourceInput& mftSource);

  // --- Reset-compatible participant list ---
  // Clears just both participants' publication sidecars -- the unconditional
  // top-of-trackFrame() step the combined DPL task's own composition
  // documents (CombinedCATrackerSpec.cxx).
  void clearPublicationSidecars() noexcept;
  // Invalidates the publication/timing bridge only (never touches either
  // participant's scratch or the shared TimeFrame) -- called on every
  // process() entry and on any non-success outcome.
  void invalidatePublication() noexcept;
  // Captures both participants' clocks and marks the publication/timing
  // bridge valid -- called only after a fully successful tracking phase.
  void markPublicationValid() noexcept;

  // --- ITS/MFT adapter publication exports: immutable per-detector exports,
  // valid only between a markPublicationValid() call and the next
  // invalidatePublication(). ITS is always ClusterSourceId{0}, MFT always
  // ClusterSourceId{1} -- fixed positions, never the caller's choice, the
  // same contract validateSources() enforces. ---
  std::optional<CommonTrackPublicationExport> getITSPublicationExport() const;
  std::optional<CommonTrackPublicationExport> getMFTPublicationExport() const;

  // --- Caller readback: forwards this API to the combined DPL task
  // without that task needing to name any of these concrete types itself. ---
  const LegacyTrackerScratchITS& getITSScratch() const noexcept { return mITSParticipant.getScratch(); }
  // M6d: MFT's own participant now owns SurfaceTrackingScratch, not
  // LegacyTrackerScratch<MFTNLayers> -- see LegacyCATrackingParticipant.h's
  // LegacyCATrackingParticipantMFT alias.
  const SurfaceTrackingScratch& getMFTScratch() const noexcept { return mMFTParticipant.getScratch(); }
  const ITSSharedClusterCompatibility& getITSSharedClusterCompatibility() const noexcept { return *mITSParticipant.getITSSharedClusterCompatibility(); }
  const MFTPublicationCompatibility& getMFTPublicationCompatibility() const noexcept { return *mMFTParticipant.getMFTPublicationCompatibility(); }
  gsl::span<const SurfaceId> getITSOrderedSurfaces() const noexcept { return mITSPlan->getConfigurationKey().orderedSurfaces; }
  gsl::span<const SurfaceId> getMFTOrderedSurfaces() const noexcept { return mMFTPlan->getConfigurationKey().orderedSurfaces; }
  DetectorLayoutView getITSLayoutView() const noexcept { return mITSPlan->getLayoutView(0); }
  DetectorLayoutView getMFTLayoutView() const noexcept { return mMFTPlan->getLayoutView(0); }

 private:
  std::optional<DetectorLayoutSet> mITSPlan;
  std::optional<DetectorLayoutSet> mMFTPlan;

  LegacyCATrackingParticipantITS mITSParticipant;
  LegacyCATrackingParticipantMFT mMFTParticipant;

  std::array<TrackingParticipant*, 2> mSchedule{};

  std::optional<ClockTimingPublicationView> mITSClock;
  std::optional<ClockTimingPublicationView> mMFTClock;
  bool mPublicationValid = false;
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE
#endif /* ALICEO2_ITSMFT_TRACKING_ITSMFTLEGACYPARTICIPANTSET_H_ */

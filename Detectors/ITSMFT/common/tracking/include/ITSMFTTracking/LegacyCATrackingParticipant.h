// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// E0/M2a (GenericTrackingEngineMigration.md; ADR 0007): the temporary
// concrete TrackingParticipant implementation CombinedTimeFrameCoordinator
// uses for each of its two legacy detector legs (ITS: NLayers=7, MFT:
// NLayers=10). Owns everything ADR 0007 still classifies as temporary
// per-detector legacy state -- LegacyTrackerScratch<NLayers>,
// TrackerTraits<NLayers>, Tracker<NLayers>, this leg's own
// DetectorTraversalBinding, its TrackingParameters, and its detector
// compatibility sidecar (ITSSharedClusterCompatibility for ITS,
// MFTPublicationCompatibility for MFT) -- entirely behind the generic
// TrackingParticipant boundary (TrackingParticipant.h): none of that state,
// or ITS/MFT/NLayers themselves, is nameable from TrackingEngine.h/
// TrackingParticipant.h's own public headers, which this class includes but
// never edits.
//
// This class never loads input itself, matching TrackingParticipant.h's own
// contract: the coordinator drives MultiSourceTimeFrameLoader::
// loadITSAndMFT() directly, unchanged, and only calls
// TrackingEngine::executeEvent()/resetEvent() once that atomic load has
// already committed -- or, on a load failure, without ever calling
// executeEvent() at all. track()/eventReset() below both assume that
// precondition; see TrackingParticipant.h for the exact contract.

#ifndef ALICEO2_ITSMFT_TRACKING_LEGACYCATRACKINGPARTICIPANT_H_
#define ALICEO2_ITSMFT_TRACKING_LEGACYCATRACKINGPARTICIPANT_H_

#ifndef GPUCA_GPUCODE

#include <memory>
#include <optional>
#include <vector>

#include <oneapi/tbb/task_arena.h>

#include "ITSMFTTracking/CATracker.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/DetectorLayoutSet.h"
#include "ITSMFTTracking/DetectorTraversalBinding.h"
#include "ITSMFTTracking/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/LegacyTrackerScratch.h"
#include "ITSMFTTracking/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/ParticipantId.h"
#include "ITSMFTTracking/TrackerTraits.h"
// Reused only for the ITSSharedClusterCompatibilityOwner<NLayers>/
// MFTPublicationCompatibilityOwner<NLayers> template-specialized mixins
// (empty for the "wrong" NLayers, one real sidecar member for the "right"
// one) -- the exact same per-NLayers sidecar-ownership mechanism
// ITSMFTTrackingInterface<NLayers> already uses, reused rather than
// reinvented.
#include "ITSMFTTracking/TrackingInterface.h"
#include "ITSMFTTracking/TrackingParticipant.h"

namespace o2::itsmft::tracking
{

template <int NLayers>
class LegacyCATrackingParticipant final : public TrackingParticipant,
                                          private ITSSharedClusterCompatibilityOwner<NLayers>,
                                          private MFTPublicationCompatibilityOwner<NLayers>
{
 public:
  static_assert(NLayers == ITSNLayers || NLayers == MFTNLayers,
                "LegacyCATrackingParticipant supports ITS (7) and MFT (10) layer counts only");
  static constexpr o2::detectors::DetID::ID DetId = detIdFromNLayers<NLayers>();

  using ScratchN = LegacyTrackerScratch<NLayers>;

  LegacyCATrackingParticipant(ParticipantId id, std::vector<TrackingParameters> params);

  // Tracker<NLayers>/TrackerTraits<NLayers> bind to sibling addresses at
  // construction (adoptScratch()/adoptITSSharedClusterCompatibility()/
  // adoptMFTPublicationCompatibility()) and to the addresses
  // adoptDetectorTraversalBinding()/adoptDetectorLayoutSet() bind
  // immediately afterward; relocating this object would silently dangle
  // every one of those bound pointers -- the same non-relocatable contract
  // CombinedTimeFrameCoordinator itself already documents.
  LegacyCATrackingParticipant(const LegacyCATrackingParticipant&) = delete;
  LegacyCATrackingParticipant& operator=(const LegacyCATrackingParticipant&) = delete;
  LegacyCATrackingParticipant(LegacyCATrackingParticipant&&) = delete;
  LegacyCATrackingParticipant& operator=(LegacyCATrackingParticipant&&) = delete;

  // --- Coordinator-only setup, called once before the first track() call.
  // Not part of the generic TrackingParticipant interface: the coordinator
  // holds this concrete type directly (never only a TrackingParticipant*)
  // specifically to reach these. ---
  void adoptDetectorTraversalBinding(std::unique_ptr<DetectorTraversalBinding> binding);
  void adoptDetectorLayoutSet(const DetectorLayoutSet& plan);
  void adoptFrame(TimeFrame& frame);
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool);
  void setBz(float bz);
  void setNThreads(int n);
  void configureRofTables(const ROFTimingConfig& timing, uint32_t nROFsTF);

  // Clears just this participant's publication sidecar -- the same
  // unconditional top-of-process() step CombinedTimeFrameCoordinator::
  // process() itself already documented before this slice: a prior
  // successful run leaves the sidecar sealed, and the very next TF's first
  // accepted track would otherwise fail the already-sealed guard.
  // Deliberately narrower than eventReset(): it never touches the scratch.
  void clearPublicationSidecar() noexcept;

  // This leg's own DropTFUponFailure, needed by the coordinator to classify
  // an atomic *load* failure (a decision made before track() is ever
  // reachable, so it cannot come from a ParticipantTrackingResult).
  bool getDropTFUponFailure() const noexcept { return mParams[0].DropTFUponFailure; }

  // --- TrackingParticipant (ADR 0007 decision 5) ---
  ParticipantId id() const noexcept override { return mId; }
  gsl::span<const SurfaceId> ownedSurfaces() const noexcept override;
  ParticipantTrackingResult track(TimeFrame& frame) override;
  void eventReset(TimeFrame& frame) noexcept override;
  std::optional<ParticipantPublicationExport> publicationExport() const override;

  // --- Coordinator-only readback, concrete type only: never exposed
  // through TrackingParticipant's own public contract. ---
  ScratchN& getScratch() noexcept { return mScratch; }
  const ScratchN& getScratch() const noexcept { return mScratch; }
  // Returns nullptr for the "wrong" NLayers instantiation (mirrors
  // ITSMFTTrackingInterface<NLayers>'s own getITSSharedClusterCompatibility
  // ()/getMFTPublicationCompatibility() exactly) -- the coordinator only
  // ever calls the accessor matching this leg's own DetId.
  const ITSSharedClusterCompatibility* getITSSharedClusterCompatibility() const noexcept;
  const MFTPublicationCompatibility* getMFTPublicationCompatibility() const noexcept;

 private:
  void clearCompatibility() noexcept;

  ParticipantId mId;
  std::vector<TrackingParameters> mParams;
  ScratchN mScratch;
  std::unique_ptr<DetectorTraversalBinding> mBinding;
  // Non-owning: the coordinator's own static combined plan data
  // (DetectorLayoutSet) outlives every participant -- see
  // CombinedTimeFrameCoordinator.h's own ownership doc.
  const DetectorLayoutSet* mPlan = nullptr;
  TrackerTraits<NLayers> mTraits;
  Tracker<NLayers> mTracker;
  std::shared_ptr<tbb::task_arena> mArena;
  // Engaged only between a successful track() and the next eventReset()
  // (or construction) -- gates publicationExport(), matching
  // TrackingParticipant.h's own "engaged only after a successful track()"
  // contract.
  bool mTracked = false;
};

using LegacyCATrackingParticipantITS = LegacyCATrackingParticipant<ITSNLayers>;
using LegacyCATrackingParticipantMFT = LegacyCATrackingParticipant<MFTNLayers>;

extern template class LegacyCATrackingParticipant<ITSNLayers>;
extern template class LegacyCATrackingParticipant<MFTNLayers>;

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE
#endif /* ALICEO2_ITSMFT_TRACKING_LEGACYCATRACKINGPARTICIPANT_H_ */

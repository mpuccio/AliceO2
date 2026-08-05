// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// E0/M2a-M2c (GenericTrackingEngineMigration.md; ADR 0007): this is the
// concrete plan-driven TrackingParticipant used by the ITS/MFT application
// adapter. It owns one SurfaceTrackingScratch, one SurfacePlanBinding, the
// layer-count-specific Tracker/TrackerTraits pair, the tracking parameters,
// and its detector compatibility sidecar entirely behind the generic
// TrackingParticipant boundary.
//
// This class never loads input itself, matching TrackingParticipant.h's own
// contract: it only exposes its own loadTarget() (the loader's surface-backed
// target, bound to its own scratch) for the owning workflow application to
// build an AtomicLoadBinding from; the actual
// MultiSourceTimeFrameLoader::loadEvent() call, and the subsequent
// TrackingEngine::executeEvent()/resetEvent() call once that atomic load has
// committed (or resetEvent() alone, without ever calling executeEvent(), on
// a load failure), both happen in the combined DPL task. track()/eventReset()
// below assume that precondition; see TrackingParticipant.h for the exact
// contract.

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEPLANTRACKINGPARTICIPANT_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEPLANTRACKINGPARTICIPANT_H_

#ifndef GPUCA_GPUCODE

#include <memory>
#include <optional>
#include <vector>

#include <oneapi/tbb/task_arena.h>

#include "ITSMFTTracking/DetectorPublicationAdapter.h"
#include "ITSMFTTracking/CATracker.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/DetectorLayoutSet.h"
#include "ITSMFTTracking/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/ParticipantId.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TrackingOperationAdapter.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
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
class SurfacePlanTrackingParticipant final : public TrackingParticipant,
                                             private TrackingOperationAdapter,
                                             private ITSSharedClusterCompatibilityOwner<NLayers>,
                                             private MFTPublicationCompatibilityOwner<NLayers>
{
 public:
  static_assert(NLayers == ITSNLayers || NLayers == MFTNLayers,
                "SurfacePlanTrackingParticipant supports ITS (7) and MFT (10) layer counts only");
  static constexpr o2::detectors::DetID::ID DetId = detIdFromNLayers<NLayers>();

  SurfacePlanTrackingParticipant(ParticipantId id, std::vector<TrackingParameters> params);

  // Tracker/TrackerTraits bind to sibling addresses at
  // construction (adoptScratch()/adoptITSSharedClusterCompatibility()/
  // adoptMFTPublicationCompatibility()) and to the addresses
  // adoptSurfacePlanBinding()/adoptDetectorLayoutSet() bind
  // immediately afterward; relocating this object would silently dangle
  // every one of those bound pointers -- the non-relocatable contract
  // required by the owning workflow application.
  SurfacePlanTrackingParticipant(const SurfacePlanTrackingParticipant&) = delete;
  SurfacePlanTrackingParticipant& operator=(const SurfacePlanTrackingParticipant&) = delete;
  SurfacePlanTrackingParticipant(SurfacePlanTrackingParticipant&&) = delete;
  SurfacePlanTrackingParticipant& operator=(SurfacePlanTrackingParticipant&&) = delete;

  // --- Owning-application setup, called once before the first track() call.
  // Not part of the generic TrackingParticipant interface: the owning
  // workflow holds this concrete type directly (never only a
  // TrackingParticipant*) specifically to reach these. ---
  void adoptSurfacePlanBinding(std::unique_ptr<SurfacePlanBinding> binding);
  void adoptDetectorLayoutSet(const DetectorLayoutSet& plan);
  void adoptFrame(TimeFrame& frame);
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool);
  void setBz(float bz);
  void setNThreads(int n);
  void configureRofTables(const ROFTimingConfig& timing, uint32_t nROFsTF);

  // This leg's own load target for MultiSourceTimeFrameLoader::loadEvent()
  // (M2b): bound once, at construction, to mScratch -- see mLoadTarget's
  // own doc below. Not part of the generic TrackingParticipant interface;
  // the owning workflow (or any other adapter driving the generic atomic
  // loading transaction) reaches this only through the concrete participant
  // type.
  MultiSourceTimeFrameLoader::LoadTarget& loadTarget() noexcept { return mLoadTarget; }

  // Clears just this participant's publication sidecar. The workflow calls
  // this unconditionally at the top of each trackFrame(): a prior successful
  // run leaves the sidecar sealed, and the next TF's first accepted track
  // would otherwise fail the already-sealed guard. Deliberately narrower than
  // eventReset(): it never touches the scratch.
  void clearPublicationSidecar() noexcept;

  // This leg's own DropTFUponFailure, needed by the owning workflow to
  // classify an atomic *load* failure (a decision made before track() is ever
  // reachable, so it cannot come from a ParticipantTrackingResult).
  bool getDropTFUponFailure() const noexcept { return mParams[0].DropTFUponFailure; }

  // --- TrackingParticipant (ADR 0007 decision 5) ---
  ParticipantId id() const noexcept override { return mId; }
  gsl::span<const SurfaceId> ownedSurfaces() const noexcept override;
  ParticipantTrackingResult track(TimeFrame& frame) override;
  void eventReset(TimeFrame& frame) noexcept override;
  std::optional<ParticipantPublicationExport> publicationExport() const override;

  // --- Owning-application readback, concrete type only: never exposed
  // through TrackingParticipant's own public contract. ---
  SurfaceTrackingScratch& getScratch() noexcept { return mScratch; }
  const SurfaceTrackingScratch& getScratch() const noexcept { return mScratch; }
  // Returns nullptr for the "wrong" NLayers instantiation (mirrors
  // ITSMFTTrackingInterface<NLayers>'s own getITSSharedClusterCompatibility
  // ()/getMFTPublicationCompatibility() exactly) -- the owning
  // the owning workflow only ever calls the accessor matching this leg's own
  // DetId.
  const ITSSharedClusterCompatibility* getITSSharedClusterCompatibility() const noexcept;
  const MFTPublicationCompatibility* getMFTPublicationCompatibility() const noexcept;

 private:
  void clearCompatibility() noexcept;

  bool refitSeed(const TrackSeed& seed,
                 const TrackingParameters& params,
                 float bz,
                 SurfaceTrackingScratch& scratch,
                 gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                 SurfaceCatalogView surfaceCatalog,
                 ClusterSourceId expectedSource,
                 TrackingCandidate& candidate) override;
  bool completeAccepted(gsl::span<const TrackingCandidate> candidates,
                        const TrackingParameters& params,
                        const SurfaceTrackingScratch& scratch,
                        bool final) override;
  void resetAdapterState() noexcept override { clearCompatibility(); }

  ParticipantId mId;
  std::vector<TrackingParameters> mParams;
  SurfaceTrackingScratch mScratch;
  std::unique_ptr<SurfacePlanBinding> mBinding;
  // Non-owning: the owning workflow's static combined plan data
  // (DetectorLayoutSet) outlives every participant.
  const DetectorLayoutSet* mPlan = nullptr;
  TrackerTraits mTraits;
  Tracker mTracker;
  DetectorPublicationAdapter<NLayers> mDetectorPublicationAdapter;
  // Adapter-edge ownership for this leg's fixed-capacity timing/mask tables.
  // SurfaceTrackingScratch receives only non-owning runtime views.
  o2::its::ROFOverlapTable<NLayers> mROFOverlapTable;
  o2::its::ROFVertexLookupTable<NLayers> mROFVertexLookupTable;
  o2::its::ROFMaskTable<NLayers> mMultiplicityMask;
  o2::its::ROFMaskTable<NLayers> mUPCMask;
  // Bound to mScratch above at construction (M2b) -- see
  // MultiSourceTimeFrameLoader::LoadTargetImplSurface's own doc. Declared
  // after mScratch so it is constructed after (and therefore only ever
  // binds an already-existing) mScratch.
  MultiSourceTimeFrameLoader::LoadTargetImplSurface mLoadTarget;
  std::shared_ptr<tbb::task_arena> mArena;
  // Engaged only between a successful track() and the next eventReset()
  // (or construction) -- gates publicationExport(), matching
  // TrackingParticipant.h's own "engaged only after a successful track()"
  // contract.
  bool mTracked = false;
};

using SurfacePlanTrackingParticipantITS = SurfacePlanTrackingParticipant<ITSNLayers>;
using SurfacePlanTrackingParticipantMFT = SurfacePlanTrackingParticipant<MFTNLayers>;

extern template class SurfacePlanTrackingParticipant<ITSNLayers>;
extern template class SurfacePlanTrackingParticipant<MFTNLayers>;

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE
#endif /* ALICEO2_ITSMFT_TRACKING_SURFACEPLANTRACKINGPARTICIPANT_H_ */

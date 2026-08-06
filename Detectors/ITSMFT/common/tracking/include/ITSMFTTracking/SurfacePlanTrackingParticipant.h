// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Concrete plan-driven TrackingParticipant used by the ITS/MFT application
// adapter. Static graph configuration and generic event workspace are borrowed
// from its reusable TimeFrame; this object retains only adapter sidecars.
//
// This class never loads input itself. The workflow owns the direct loading
// composition and invokes tracking only after a successful frame commit.

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEPLANTRACKINGPARTICIPANT_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEPLANTRACKINGPARTICIPANT_H_

#ifndef GPUCA_GPUCODE

#include <memory>
#include <optional>
#include <vector>

#include <oneapi/tbb/task_arena.h>

#include "ITSMFTTracking/DetectorPublicationAdapter.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/SurfaceGraph.h"
#include "ITSMFTTracking/ParticipantId.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TrackingOperationAdapter.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "ITSMFTTracking/TrackingParticipant.h"

namespace o2::itsmft::tracking
{

template <int NLayers>
class SurfacePlanTrackingParticipant final : public TrackingParticipant,
                                             private TrackingOperationAdapter
{
 public:
  static_assert(NLayers == ITSNLayers || NLayers == MFTNLayers,
                "SurfacePlanTrackingParticipant supports ITS (7) and MFT (10) layer counts only");
  static constexpr o2::detectors::DetID::ID DetId = detIdFromNLayers<NLayers>();

  SurfacePlanTrackingParticipant(ParticipantId id, ClusterSourceId source);

  // Tracker/TrackerTraits bind to sibling addresses at construction;
  // relocating this object would silently dangle those pointers.
  SurfacePlanTrackingParticipant(const SurfacePlanTrackingParticipant&) = delete;
  SurfacePlanTrackingParticipant& operator=(const SurfacePlanTrackingParticipant&) = delete;
  SurfacePlanTrackingParticipant(SurfacePlanTrackingParticipant&&) = delete;
  SurfacePlanTrackingParticipant& operator=(SurfacePlanTrackingParticipant&&) = delete;

  // --- Owning-application setup, called once before the first track() call.
  // Not part of the generic TrackingParticipant interface: the owning
  // workflow holds this concrete type directly (never only a
  // TrackingParticipant*) specifically to reach these. ---
  TrackerInitializationResult initialize(TimeFrame& frame, const TrackerInitialization& configuration);
  void adoptConfiguredFrame(TimeFrame& frame);
  void setBz(float bz);
  void setNThreads(int n);
  void bindPublicationAdapter(DetectorPublicationAdapter<NLayers>& adapter) noexcept { mDetectorPublicationAdapter = &adapter; }
  void setROFViews(RuntimeROFViews views) noexcept;
  void clearROFViews() noexcept;

  // Clears just this participant's publication sidecar. The workflow calls
  // this unconditionally at the top of each trackFrame(): a prior successful
  // run leaves the sidecar sealed, and the next TF's first accepted track
  // would otherwise fail the already-sealed guard. Deliberately narrower than
  // eventReset(): it never touches the scratch.
  void clearPublicationSidecar() noexcept;

  // This leg's own DropTFUponFailure, needed by the owning workflow to
  // classify an atomic *load* failure (a decision made before track() is ever
  // reachable, so it cannot come from a ParticipantTrackingResult).
  bool getDropTFUponFailure() const noexcept
  {
    const auto parameters = mFrame != nullptr && mFrame->isConfigured() ? mFrame->getTrackingParameters(mSource) : gsl::span<const TrackingParameters>{};
    return !parameters.empty() && parameters[0].DropTFUponFailure;
  }

  // --- TrackingParticipant (ADR 0007 decision 5) ---
  ParticipantId id() const noexcept override { return mId; }
  gsl::span<const SurfaceId> ownedSurfaces() const noexcept override;
  ParticipantTrackingResult track(TimeFrame& frame) override;
  void eventReset(TimeFrame& frame) noexcept override;
  std::optional<ParticipantPublicationExport> publicationExport() const override;

  // --- Owning-application readback, concrete type only: never exposed
  // through TrackingParticipant's own public contract. ---
  SurfaceTrackingScratch& getScratch() { return mFrame->getWorkspace(mSource); }
  const SurfaceTrackingScratch& getScratch() const { return mFrame->getWorkspace(mSource); }
  RuntimeROFViews getROFViews() const noexcept { return mFrame == nullptr ? RuntimeROFViews{} : getScratch().getROFViews(); }
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
  ClusterSourceId mSource;
  TimeFrame* mFrame = nullptr;
  TrackerTraits mTraits;
  Tracker mTracker;
  DetectorPublicationAdapter<NLayers>* mDetectorPublicationAdapter = nullptr;
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

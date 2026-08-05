// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKINGOPERATIONADAPTER_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKINGOPERATIONADAPTER_H_

#ifndef GPUCA_GPUCODE

#include <cstdint>

#include <gsl/span>

#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/SurfaceCatalogView.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TimeFrame.h"

namespace o2::itsmft::tracking
{

// Detector-neutral result of one successful seed refit. The candidate keeps
// only values needed by the shared acceptance/ordering flow and the narrow
// application adapter call. Typed output tracks and publication sidecars
// remain outside the common tracker.
struct TrackingCandidate {
  TrackSeed seed;
  SurfaceKinematicState innerState{};
  SurfaceKinematicState outerState{};
  o2::its::TimeStamp timestamp{};
  float chi2{0.f};
  float phi{0.f};
  float eta{0.f};
  double charge{0.};
  bool sharedClusters{false};

  int getNumberOfClusters() const noexcept { return seed.getActiveSurfaceCount(); }
  int getClusterIndex(int position) const noexcept { return seed.getCluster(position); }
  int getFirstClusterLayer() const noexcept { return seed.getSurfaceMask().first(); }
};

// Narrow application-operation seam used while typed ITS/MFT refit and
// publication hooks complete their separately gated M7e migration. It is an
// operation adapter, not a coordinator: the common core owns no detector
// identity, typed output, timing tables, or workflow state, and this object
// owns no event lifecycle.
class TrackingOperationAdapter
{
 public:
  virtual ~TrackingOperationAdapter() = default;

  virtual bool refitSeed(const TrackSeed& seed,
                         const TrackingParameters& params,
                         float bz,
                         SurfaceTrackingScratch& scratch,
                         gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                         SurfaceCatalogView surfaceCatalog,
                         ClusterSourceId expectedSource,
                         TrackingCandidate& candidate) = 0;

  virtual bool publishAccepted(TimeFrame& frame,
                               const TrackingCandidate& candidate,
                               gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                               SurfaceCatalogView surfaceCatalog) = 0;

  virtual bool haveSamePolarity(const TrackingCandidate& first,
                                const TrackingCandidate& second) const noexcept = 0;

  virtual bool sealAccepted(gsl::span<const TrackingCandidate> candidates) = 0;

  virtual void clearPublicationState() noexcept = 0;
};

} // namespace o2::itsmft::tracking

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_TRACKINGOPERATIONADAPTER_H_

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKINGOPERATIONADAPTER_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKINGOPERATIONADAPTER_H_

#ifndef GPUCA_GPUCODE

#include <cstdint>
#include <limits>

#include <gsl/span>

#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/CommonTrack.h"
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
  CommonTrack track{};
  float phi{0.f};
  float eta{0.f};
  double charge{0.};
  uint32_t commonTrackIndex{std::numeric_limits<uint32_t>::max()};

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

  // Called at each serial accepted-result boundary. The adapter may consume
  // the generic candidates and operation-local policy/scratch views to stage
  // detector compatibility state; only the final call may publish sidecars.
  // No typed track or compatibility flag is supplied by the core.
  virtual bool completeAccepted(gsl::span<const TrackingCandidate> candidates,
                                const TrackingParameters& params,
                                const SurfaceTrackingScratch& scratch,
                                bool final) = 0;

  // Lifecycle reset is deliberately an adapter-edge operation. It clears
  // only detector compatibility state; TimeFrame and scratch reset remain
  // owned by the generic tracker.
  virtual void resetAdapterState() noexcept = 0;
};

} // namespace o2::itsmft::tracking

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_TRACKINGOPERATIONADAPTER_H_

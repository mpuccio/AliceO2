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
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
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

// The only operation the generic road stage requires from an application is a
// call-scoped refit. Publication and reset are workflow operations after the
// generic tracker transaction has completed.
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
};

} // namespace o2::itsmft::tracking

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_TRACKINGOPERATIONADAPTER_H_

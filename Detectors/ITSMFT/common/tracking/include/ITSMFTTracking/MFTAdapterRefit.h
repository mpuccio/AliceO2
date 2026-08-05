// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_MFTADAPTERREFIT_H_
#define ALICEO2_ITSMFT_TRACKING_MFTADAPTERREFIT_H_

#ifndef GPUCA_GPUCODE

#include "ITSMFTTracking/MFTCATrack.h"
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/SurfaceKinematicStateLegacyAdapters.h"

namespace o2::itsmft::tracking
{

// Adapter-only typed export retained for the normalized-refit compatibility
// fixture and the MFT application edge. The generic tracker calls the
// SurfaceKinematicState overload in MFTFwdTrackHelpers.h instead.
bool refitTrackFwd(const TrackSeed& seed,
                   MFTCATrack& track,
                   const SurfaceTrackingScratch& tf,
                   const TrackingParameters& params,
                   float bz,
                   gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                   SurfaceCatalogView surfaceCatalog,
                   ClusterSourceId expectedSource);

} // namespace o2::itsmft::tracking

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_MFTADAPTERREFIT_H_

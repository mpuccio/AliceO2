// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#include "ITSMFTTracking/MFTAdapterRefit.h"

#include "ITSMFTTracking/SurfaceKinematicStateLegacyAdapters.h"
#include "ITStracking/Constants.h"
#include "MFTTracking/Constants.h"

namespace o2::itsmft::tracking
{

bool refitTrackFwd(const TrackSeed& seed,
                   MFTCATrack& track,
                   const SurfaceTrackingScratch& tf,
                   const TrackingParameters& params,
                   float bz,
                   gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                   SurfaceCatalogView surfaceCatalog,
                   ClusterSourceId expectedSource)
{
  SurfaceKinematicState paramIn{};
  SurfaceKinematicState paramOut{};
  float chi2 = 0.f;
  if (!refitTrackFwd(seed, tf, params, bz, layerMeasurements, surfaceCatalog, expectedSource, paramIn, paramOut, chi2)) {
    return false;
  }

  o2::track::TrackParCovFwd inFwd{};
  o2::track::TrackParCovFwd outFwd{};
  if (!legacy::exportLegacyForwardTrackParCov(paramIn, inFwd) ||
      !legacy::exportLegacyForwardTrackParCov(paramOut, outFwd)) {
    return false;
  }
  inFwd.setTrackChi2(chi2);
  auto& mftTr = track.getTrack();
  static_cast<o2::track::TrackParCovFwd&>(mftTr) = inFwd;
  mftTr.setOutParam(outFwd);
  mftTr.setCA(true);
  mftTr.setNumberOfPoints(seed.getHitLayerMask().count());
  track.setPattern(0);
  for (int layer = 0; layer < o2::mft::constants::mft::LayersNumber; ++layer) {
    if (!seed.getHitLayerMask().has(layer)) {
      track.setClusterIndex(layer, o2::its::constants::UnusedIndex);
      continue;
    }
    const int clIdx = seed.getCluster(layer);
    track.setClusterIndex(layer, clIdx);
    track.setClusterSize(layer, tf.getClusterSize(layer, clIdx));
  }
  return true;
}

} // namespace o2::itsmft::tracking

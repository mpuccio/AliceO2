// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file MFTFwdTrackHelpers.cxx
/// \brief MFT CA final forward refit using the shared native driver.
///

#include "ITSMFTTracking/detail/MFTFwdTrackHelpers.h"

#include "Framework/Logger.h"
#include "ITSMFTTracking/NativeRefitDriver.h"
#include "ITStracking/Constants.h"

namespace o2::itsmft::tracking
{

namespace
{
constexpr int kMFTLayers = o2::mft::constants::mft::LayersNumber;
}

bool refitTrackFwd(const TrackSeed& seed,
                   const SurfaceTrackingScratch& tf,
                   const TrackingParameters& params,
                   float bz,
                   gsl::span<const gsl::span<const GlobalMeasurement>> layerGlobals,
                   gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                   SurfaceCatalogView surfaceCatalog,
                   ClusterSourceId expectedSource,
                   SurfaceKinematicState& paramIn,
                   SurfaceKinematicState& paramOut,
                   float& chi2)
{
  const auto hitMask = seed.getHitLayerMask();

  for (int layer = 0; layer < static_cast<int>(layerMeasurements.size()); ++layer) {
    if (!hitMask.has(layer)) {
      continue;
    }
    const int clIdx = seed.getCluster(layer);
    if (clIdx == o2::its::constants::UnusedIndex) {
      continue;
    }
    if (clIdx < 0 || clIdx >= static_cast<int>(layerMeasurements[layer].size())) {
      LOGP(warn, "MFT CA forward refit: invalid cluster index {} on layer {}", clIdx, layer);
      return false;
    }
    const auto& measurement = layerMeasurements[layer][clIdx];
    const auto& global = layerGlobals[layer][clIdx];
    if (!surfaceCatalog.hasSurface(global.surface) || global.cluster.source != expectedSource ||
        global.cluster.index != static_cast<uint32_t>(tf.getClusterExternalIndex(layer, clIdx))) {
      LOGP(warn, "MFT CA forward refit: invalid global measurement identity on layer {} clIdx {}", layer, clIdx);
      return false;
    }
    if (!std::isfinite(measurement.frame.u) || !std::isfinite(measurement.frame.v) || !std::isfinite(measurement.frame.q) ||
        !std::isfinite(measurement.covariance.uu) || !std::isfinite(measurement.covariance.vv) ||
        measurement.covariance.uu < 0.f || measurement.covariance.vv < 0.f) {
      LOGP(warn, "MFT CA forward refit: invalid normalized measurement on layer {} clIdx {}", layer, clIdx);
      return false;
    }
  }

  paramIn = {};
  paramOut = {};
  chi2 = 0.f;
  OperationFailureReason reason{};
  if (!fitTrackSeedLegs(seed, layerGlobals, layerMeasurements, surfaceCatalog, bz,
                        params.ShiftRefToCluster, params.MaxChi2ClusterAttachment, params.MaxChi2NDF,
                        params.RepeatRefitOut, gsl::span<const float>(params.MinPt),
                        paramIn, paramOut, chi2, reason)) {
    LOGP(warn, "MFT CA forward refit: fitTrackSeedLegs failed, reason={}", static_cast<int>(reason));
    return false;
  }

  return true;
}

} // namespace o2::itsmft::tracking

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
/// \brief MFT CA final forward refit -- M5d: shared native driver, no frozen
/// legacy o2::mft::TrackFitter/TrackLTF Kalman engine.
///

#include "ITSMFTTracking/MFTFwdTrackHelpers.h"

#include <cmath>

#include "Framework/Logger.h"
#include "ITSMFTTracking/NativeRefitDriver.h"
#include "ITStracking/Constants.h"

namespace o2::itsmft::tracking
{

namespace
{
constexpr int kMFTLayers = o2::mft::constants::mft::LayersNumber;
}

// M5d: replaces the frozen o2::mft::TrackFitter<TrackLTF>/TrackLTFL Kalman
// engine (MFTTracking/TrackFitter.h) with fitTrackSeedLegs (NativeRefitDriver.h),
// the same shared, descriptor-driven driver the barrel/ITS branch
// (the ITS adapter's native refit operation) now uses -- the intentional, approved
// physics departure recorded in doc/decisions/0008-native-refit-activation.md.
// The leg structure this milestone activates (inward/outward/optional-repeat,
// Section "Required migration") is not a port of TrackLTF's own two-direction
// linear-track-finder algorithm; it is the same three-leg Kalman sequencing
// the barrel branch already uses, applied here to a Forward-family
// SurfaceKinematicState. Numerical output is expected to differ from the
// retired engine; see the design note for characterization evidence.
bool refitTrackFwd(const TrackSeed& seed,
                   const SurfaceTrackingScratch& tf,
                   const TrackingParameters& params,
                   float bz,
                   gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                   SurfaceCatalogView surfaceCatalog,
                   ClusterSourceId expectedSource,
                   SurfaceKinematicState& paramIn,
                   SurfaceKinematicState& paramOut,
                   float& chi2)
{
  const auto hitMask = seed.getHitLayerMask();

  // Defensive re-check of the ClusterRef identity contract that
  // TrackerTraits::initialiseTimeFrame() already established for every entry
  // of mLayerMeasurements (NormalizedMeasurementMismatch) -- unchanged from
  // the pre-M5d implementation, see MFTFwdTrackHelpers.h's own doc.
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
    const int extIdx = tf.getClusterExternalIndex(layer, clIdx);
    if (!measurement.cluster.isValid() || measurement.cluster.source != expectedSource ||
        extIdx < 0 || static_cast<uint32_t>(extIdx) != measurement.cluster.index) {
      LOGP(warn, "MFT CA forward refit: normalized measurement identity mismatch on layer {} clIdx {}", layer, clIdx);
      return false;
    }
    if (!std::isfinite(measurement.global.x) || !std::isfinite(measurement.global.y) || !std::isfinite(measurement.global.z) ||
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
  if (!fitTrackSeedLegs<kMFTLayers>(seed, layerMeasurements, surfaceCatalog, bz,
                                    params.ShiftRefToCluster, params.MaxChi2ClusterAttachment, params.MaxChi2NDF,
                                    params.RepeatRefitOut, gsl::span<const float>(params.MinPt),
                                    paramIn, paramOut, chi2, reason)) {
    LOGP(warn, "MFT CA forward refit: fitTrackSeedLegs failed, reason={}", static_cast<int>(reason));
    return false;
  }

  if (params.TrackletMinAbsX > 0.f) {
    if (std::abs(paramOut.parameters[0]) < params.TrackletMinAbsX) {
      return false;
    }
    for (int layer = 0; layer < static_cast<int>(layerMeasurements.size()); ++layer) {
      if (!hitMask.has(layer)) {
        continue;
      }
      const int clIdx = seed.getCluster(layer);
      if (clIdx != o2::its::constants::UnusedIndex &&
          std::abs(layerMeasurements[layer][clIdx].global.x) < params.TrackletMinAbsX) {
        return false;
      }
    }
  }

  return true;
}

} // namespace o2::itsmft::tracking

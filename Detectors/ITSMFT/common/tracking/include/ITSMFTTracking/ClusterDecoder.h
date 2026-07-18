// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_CLUSTERDECODER_H_
#define ALICEO2_ITSMFT_TRACKING_CLUSTERDECODER_H_

#include <cstdint>

#include <gsl/gsl>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/SurfaceId.h"

namespace o2::itsmft::tracking
{

// Host loading-boundary polymorphism only (Architecture.md 7.1): decoder
// implementations may call into detector geometry, but this interface, its
// implementations, and its result type must never enter device views or CA
// loops.
class ClusterDecoder
{
 public:
  virtual ~ClusterDecoder() = default;

  // Called once per source before decoding its first cluster (e.g. to fill a
  // detector geometry matrix cache). No-op by default.
  virtual void prepare() const {}

  virtual o2::itsmft::ioutils::SurfaceMeasurementDecodeResult decode(
    const o2::itsmft::CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const o2::itsmft::TopologyDictionary* dict,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool applySysErrors) const = 0;
};

// Production adapter: decodes through the integrated detector geometry
// singleton, exactly like the single-surface loadClusterSurfaceMeasurement
// path (same geometry lookup, pattern consumption, covariance calculation
// and systematic-error application, performed exactly once), plus mapping
// the decoded detector-local layer to a global SurfaceId.
template <o2::detectors::DetID::ID DetId>
class GeometryClusterDecoder final : public ClusterDecoder
{
 public:
  void prepare() const override { o2::itsmft::ioutils::fillMatrixCache(DetId); }

  o2::itsmft::ioutils::SurfaceMeasurementDecodeResult decode(
    const o2::itsmft::CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const o2::itsmft::TopologyDictionary* dict,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool applySysErrors) const override
  {
    // Check the required dictionary before asking GeometryTGeo::Instance()
    // for anything. In particular, this keeps a missing dictionary typed
    // even when detector geometry has not been loaded yet.
    if (dict == nullptr) {
      o2::itsmft::ioutils::SurfaceMeasurementDecodeResult result;
      result.error = ClusterDecodeError::MissingDictionary;
      return result;
    }
    return o2::itsmft::ioutils::loadClusterSurfaceMeasurement<DetId>(
      cluster, patterns, dict, layerToSurface, source, externalIndex, sourceROF, applySysErrors);
  }
};

using ITSGeometryClusterDecoder = GeometryClusterDecoder<o2::detectors::DetID::ITS>;
using MFTGeometryClusterDecoder = GeometryClusterDecoder<o2::detectors::DetID::MFT>;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_CLUSTERDECODER_H_ */

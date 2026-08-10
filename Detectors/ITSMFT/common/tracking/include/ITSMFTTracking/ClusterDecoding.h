// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_CLUSTERDECODING_H_
#define ALICEO2_ITSMFT_TRACKING_CLUSTERDECODING_H_

#include <cstddef>
#include <cstdint>
#include <cmath>

#include <gsl/gsl>

#include "DataFormatsITSMFT/ClusterPattern.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/GlobalMeasurement.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"

namespace o2::itsmft::tracking
{

// Typed failures at the host compact-cluster decoding boundary; the loader
// adds source, ROF, and external-cluster context when it maps them.
enum class ClusterDecodeError : uint8_t {
  None,
  MissingDictionary,
  TruncatedExplicitPattern,
  MalformedExplicitPattern,
  InvalidPatternId,
  InvalidSensor,
  InvalidLayer,
  InvalidLayerMapping,
  GeometryUnavailable,
  OtherMalformedInput
};

// Host-only cursor for the source-local explicit-pattern byte stream. It owns
// no storage and checks the complete encoded pattern before using the
// unbounded ClusterPattern iterator.
class BoundedPatternCursor
{
 public:
  explicit BoundedPatternCursor(gsl::span<const unsigned char> bytes) noexcept : mBytes(bytes) {}

  size_t consumed() const noexcept { return mPosition; }
  size_t remaining() const noexcept { return mBytes.size() - mPosition; }
  bool empty() const noexcept { return remaining() == 0; }

  ClusterDecodeError acquirePattern(o2::itsmft::ClusterPattern& pattern) noexcept
  {
    const auto available = remaining();
    if (available < 2) {
      return ClusterDecodeError::TruncatedExplicitPattern;
    }

    const auto rowSpan = mBytes[mPosition];
    const auto columnSpan = mBytes[mPosition + 1];
    if (rowSpan == 0 || columnSpan == 0 ||
        rowSpan > o2::itsmft::ClusterPattern::MaxRowSpan ||
        columnSpan > o2::itsmft::ClusterPattern::MaxColSpan) {
      return ClusterDecodeError::MalformedExplicitPattern;
    }

    const size_t nBits = static_cast<size_t>(rowSpan) * columnSpan;
    const size_t payloadBytes = (nBits + 7) / 8;
    const size_t encodedBytes = 2 + payloadBytes;
    if (available < encodedBytes) {
      return ClusterDecodeError::TruncatedExplicitPattern;
    }

    auto iterator = mBytes.begin() + mPosition;
    o2::itsmft::ClusterPattern decoded{iterator};
    if (decoded.getNPixels() == 0) {
      return ClusterDecodeError::MalformedExplicitPattern;
    }
    pattern = decoded;
    mPosition += encodedBytes;
    return ClusterDecodeError::None;
  }

 private:
  gsl::span<const unsigned char> mBytes{};
  size_t mPosition{0};
};

// Host-side facts produced by compact-cluster and geometry decoding.
struct DecodedCluster {
  GlobalPoint3F global{};
  // ITS geometry supplies its cylindrical tracking frame here. Disk
  // projection uses global coordinates directly.
  SurfaceFramePoint cylinderFrame{};
  // ALPIDE local row/column covariance. The detector projection determines
  // which normalized axes these values describe.
  SurfaceCovariance2F rowColumnCovariance{};
  ClusterShape shape{};
  uint32_t sensor{0};
  int layer{0};
};

// Decode result with the detector-local layer and its mapped SurfaceId. The
// measurement and kind are valid only when layerMapped is true; other fields
// are valid only on success, except layer on InvalidLayerMapping.
struct SurfaceMeasurementDecodeResult {
  GlobalMeasurement global{};
  SurfaceMeasurement measurement{};
  SurfaceKind kind{SurfaceKind::Cylinder};
  int layer{-1};
  bool layerMapped{false};
  ClusterDecodeError error{ClusterDecodeError::None};

  bool ok() const noexcept { return error == ClusterDecodeError::None; }
};

// Project decoded ITS facts into the accepted cylindrical convention.
inline GlobalMeasurement makeCylinderGlobalMeasurement(const DecodedCluster& decoded,
                                                       DetectorSensorId sensor,
                                                       SurfaceId surface,
                                                       ClusterRef cluster,
                                                       uint32_t sourceROF)
{
  const float sine = std::sin(decoded.cylinderFrame.frameAngle);
  const float cosine = std::cos(decoded.cylinderFrame.frameAngle);
  const auto& covariance = decoded.rowColumnCovariance;
  return GlobalMeasurement{
    decoded.global,
    std::hypot(decoded.global.x, decoded.global.y),
    {sine * sine * covariance.uu,
     -sine * cosine * covariance.uu,
     -sine * covariance.uv,
     cosine * cosine * covariance.uu,
     cosine * covariance.uv,
     covariance.vv},
    sensor,
    cluster,
    decoded.shape,
    sourceROF,
    surface};
}

// Project decoded MFT facts into z-normal, global-x/global-y disk coordinates.
// ALPIDE row is established as global x and column as global y by the MFT
// geometry decoder. No legacy TrackingFrameInfo participates in this mapping.
inline GlobalMeasurement makeDiskGlobalMeasurement(const DecodedCluster& decoded,
                                                   DetectorSensorId sensor,
                                                   SurfaceId surface,
                                                   ClusterRef cluster,
                                                   uint32_t sourceROF)
{
  return GlobalMeasurement{
    decoded.global,
    std::hypot(decoded.global.x, decoded.global.y),
    {decoded.rowColumnCovariance.uu, decoded.rowColumnCovariance.uv, 0.f,
     decoded.rowColumnCovariance.vv, 0.f, 0.f},
    sensor,
    cluster,
    decoded.shape,
    sourceROF,
    surface};
}

inline SurfaceMeasurement makeCylinderSurfaceMeasurement(const DecodedCluster& decoded)
{
  return {decoded.cylinderFrame, decoded.rowColumnCovariance};
}

inline SurfaceMeasurement makeDiskSurfaceMeasurement(const DecodedCluster& decoded)
{
  return {{decoded.global.z, decoded.global.x, decoded.global.y, 0.f},
          decoded.rowColumnCovariance};
}

inline SurfaceMeasurementDecodeResult makeCylinderMeasurementDecodeResult(
  const DecodedCluster& decoded, DetectorSensorId sensor, SurfaceId surface,
  ClusterRef cluster, uint32_t sourceROF)
{
  return {makeCylinderGlobalMeasurement(decoded, sensor, surface, cluster, sourceROF),
          makeCylinderSurfaceMeasurement(decoded), SurfaceKind::Cylinder,
          decoded.layer, true, ClusterDecodeError::None};
}

inline SurfaceMeasurementDecodeResult makeDiskMeasurementDecodeResult(
  const DecodedCluster& decoded, DetectorSensorId sensor, SurfaceId surface,
  ClusterRef cluster, uint32_t sourceROF)
{
  return {makeDiskGlobalMeasurement(decoded, sensor, surface, cluster, sourceROF),
          makeDiskSurfaceMeasurement(decoded), SurfaceKind::Disk,
          decoded.layer, true, ClusterDecodeError::None};
}

} // namespace o2::itsmft::tracking

namespace o2::itsmft::ioutils
{
void fillMatrixCache(o2::detectors::DetID::ID detId);

template <o2::detectors::DetID::ID DetId>
o2::itsmft::tracking::SurfaceMeasurementDecodeResult loadClusterSurfaceMeasurement(
  const o2::itsmft::CompClusterExt& c,
  o2::itsmft::tracking::BoundedPatternCursor& patterns,
  const o2::itsmft::TopologyDictionary* dict,
  gsl::span<const o2::itsmft::tracking::SurfaceId> layerToSurface,
  o2::itsmft::tracking::ClusterSourceId source,
  uint32_t externalClusterIndex,
  uint32_t sourceROF,
  bool applySysErrors = true);
} // namespace o2::itsmft::ioutils

namespace o2::itsmft::tracking
{

// Host-only loading boundary. Decoder implementations may call detector
// geometry, but this interface and its result never enter device views or CA
// loops.
class ClusterDecoder
{
 public:
  virtual ~ClusterDecoder() = default;

  // Called once per source before its first cluster; no-op by default.
  virtual void prepare() const {}

  virtual SurfaceMeasurementDecodeResult decode(
    const o2::itsmft::CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const o2::itsmft::TopologyDictionary* dict,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool applySysErrors) const = 0;
};

// Geometry-backed decoder. It performs the established single-pass geometry,
// pattern, covariance, and systematic-error operations, then maps the decoded
// detector layer to a global SurfaceId.
template <o2::detectors::DetID::ID DetId>
class GeometryClusterDecoder final : public ClusterDecoder
{
 public:
  void prepare() const override { o2::itsmft::ioutils::fillMatrixCache(DetId); }

  SurfaceMeasurementDecodeResult decode(
    const o2::itsmft::CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const o2::itsmft::TopologyDictionary* dict,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool applySysErrors) const override
  {
    // Preserve a typed dictionary failure before any geometry access.
    if (dict == nullptr) {
      SurfaceMeasurementDecodeResult result;
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

#endif /* ALICEO2_ITSMFT_TRACKING_CLUSTERDECODING_H_ */

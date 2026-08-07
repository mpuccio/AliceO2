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

#include <gsl/gsl>

#include "DataFormatsITSMFT/ClusterPattern.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"

namespace o2::itsmft::tracking
{

// Typed failures at the host compact-cluster decoding boundary. These values
// are intentionally independent of MultiSourceLoadError: a loader maps them
// while adding source/ROF/external-cluster context.
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

// Host-only cursor for the source-local explicit-pattern byte stream. The
// cursor owns no storage and always retains both the current position and the
// end through a span plus an offset. It validates the ClusterPattern encoding
// (row byte, column byte, then ceil(row*column/8) bitmap bytes) before calling
// ClusterPattern's iterator constructor, which is itself unbounded.
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

// Host-side facts produced by compact-cluster and geometry decoding. This is
// deliberately independent of detector output and legacy tracking types.
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

// Result of decoding a cluster when the target SurfaceId is not known ahead
// of decoding: `layer` is the detector-local layer discovered by the same
// geometry decode used to build `measurement`, and `layerMapped` reports
// whether the caller-supplied detector-layer-to-SurfaceId table covered that
// layer. `measurement`/`kind` are only meaningful when `layerMapped` is true.
// `kind` is the geometry kind (cylinder/disk) the decoder actually produced;
// it lets a caller validate the target surface's kind explicitly, without
// inferring detector geometry from surface count. `error` is the typed host
// decode failure; other fields are meaningful only when ok() is true (except
// `layer`, which may identify an InvalidLayerMapping failure).
struct SurfaceMeasurementDecodeResult {
  SurfaceMeasurement measurement{};
  SurfaceKind kind{SurfaceKind::Cylinder};
  int layer{-1};
  bool layerMapped{false};
  ClusterDecodeError error{ClusterDecodeError::None};

  bool ok() const noexcept { return error == ClusterDecodeError::None; }
};

// Project decoded ITS facts into the accepted cylindrical convention.
inline SurfaceMeasurement makeCylinderSurfaceMeasurement(const DecodedCluster& decoded,
                                                         DetectorSensorId sensor,
                                                         SurfaceId surface,
                                                         ClusterRef cluster,
                                                         uint32_t sourceROF)
{
  return SurfaceMeasurement{
    decoded.global,
    decoded.cylinderFrame,
    decoded.rowColumnCovariance,
    sensor,
    cluster,
    decoded.shape,
    sourceROF,
    surface};
}

// Project decoded MFT facts into z-normal, global-x/global-y disk coordinates.
// ALPIDE row is established as global x and column as global y by the MFT
// geometry decoder. No legacy TrackingFrameInfo participates in this mapping.
inline SurfaceMeasurement makeDiskSurfaceMeasurement(const DecodedCluster& decoded,
                                                     DetectorSensorId sensor,
                                                     SurfaceId surface,
                                                     ClusterRef cluster,
                                                     uint32_t sourceROF)
{
  return SurfaceMeasurement{
    decoded.global,
    {decoded.global.z, decoded.global.x, decoded.global.y, 0.f},
    decoded.rowColumnCovariance,
    sensor,
    cluster,
    decoded.shape,
    sourceROF,
    surface};
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
    // Check the required dictionary before asking GeometryTGeo::Instance()
    // for anything. In particular, this keeps a missing dictionary typed
    // even when detector geometry has not been loaded yet.
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

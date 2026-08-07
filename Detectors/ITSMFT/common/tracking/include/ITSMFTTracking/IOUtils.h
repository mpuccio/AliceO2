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
/// \file IOUtils.h
/// \brief Shared cluster I/O utilities for ITS and MFT (based on ITStracking/IOUtils.h)
///

#ifndef ALICEO2_ITSMFT_TRACKING_IOUTILS_H_
#define ALICEO2_ITSMFT_TRACKING_IOUTILS_H_

#include <array>
#include <vector>

#include <gsl/gsl>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "ReconstructionDataFormats/BaseCluster.h"
#include "DataFormatsITSMFT/ClusterPattern.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "MathUtils/Cartesian.h"

namespace o2::itsmft::ioutils
{

namespace detail
{
constexpr bool isSensorInGeometry(int sensor, int geometrySize) noexcept
{
  return sensor >= 0 && sensor < geometrySize;
}

constexpr bool isLayerInDetector(int layer, int detectorLayers) noexcept
{
  return layer >= 0 && layer < detectorLayers;
}

/// Whether a cluster-decoding systematic-error correction should be applied
/// for `DetId`. ITS common-CA is an explicit no-op here: the dedicated
/// ITSCommonCATrackerParam configuration does not support this feature, so
/// common ITS normalized decoding must not consult TrackerParamRef<ITS>::get()
/// (the frozen legacy ITSCATrackerParam's sysErrY2/sysErrZ2 -- see
/// TrackingConfigParam.h's doc comment on why the two namespaces are kept
/// distinct). MFT is unchanged: it keeps reading its own live
/// TrackerParamConfig<MFT> ("MFTCATrackerParam") sysErr2Row/sysErr2Col.
/// A free function (not a class member) purely so isolation tests can call it
/// directly without a geometry/dictionary fixture; not intended as broad
/// public API.
template <o2::detectors::DetID::ID DetId>
bool shouldApplySysErrors()
{
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    return false;
  } else {
    const auto& conf = o2::itsmft::tracking::TrackerParamRef<DetId>::get();
    for (int il = 0; il < o2::itsmft::tracking::TrackerParamRef<DetId>::nLayers(); il++) {
      if (conf.sysErr2Row[il] > 0.f || conf.sysErr2Col[il] > 0.f) {
        return true;
      }
    }
    return false;
  }
}

/// Adds the configured systematic-error correction to `sigma2Row`/`sigma2Col`
/// for `DetId`. ITS is an explicit no-op: see shouldApplySysErrors<ITS>()
/// above -- kept callable (rather than removed) so callers stay
/// detector-generic, but it never reads TrackerParamRef<ITS>::get().
template <o2::detectors::DetID::ID DetId>
void addSysErrors(int layerId, float& sigma2Row, float& sigma2Col)
{
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    (void)layerId;
    (void)sigma2Row;
    (void)sigma2Col;
  } else {
    const auto& conf = o2::itsmft::tracking::TrackerParamRef<DetId>::get();
    sigma2Row += conf.sysErr2Row[layerId];
    sigma2Col += conf.sysErr2Col[layerId];
  }
}
} // namespace detail

constexpr float DefClusErrorRow = o2::itsmft::SegmentationAlpide::PitchRow * 0.5f;
constexpr float DefClusErrorCol = o2::itsmft::SegmentationAlpide::PitchCol * 0.5f;
constexpr float DefClusError2Row = DefClusErrorRow * DefClusErrorRow;
constexpr float DefClusError2Col = DefClusErrorCol * DefClusErrorCol;

void fillMatrixCache(o2::detectors::DetID::ID detId);
int getClusterLayer(o2::detectors::DetID::ID detId, const CompClusterExt& cluster);

/// Decode a compact cluster through the detector geometry singleton directly
/// into the normalized surface representation.
template <o2::detectors::DetID::ID DetId>
o2::itsmft::tracking::SurfaceMeasurement loadClusterSurfaceMeasurement(
  const CompClusterExt& c,
  gsl::span<const unsigned char>::iterator& pattIt,
  const TopologyDictionary* dict,
  o2::itsmft::tracking::ClusterSourceId source,
  uint32_t externalClusterIndex,
  o2::itsmft::tracking::SurfaceId surface,
  uint32_t sourceROF,
  bool applySysErrors = true);

/// Result of decoding a cluster when the target SurfaceId is not known ahead
/// of decoding: `layer` is the detector-local layer discovered by the same
/// geometry decode used to build `measurement`, and `layerMapped` reports
/// whether the caller-supplied detector-layer-to-SurfaceId table covered that
/// layer. `measurement`/`kind` are only meaningful when `layerMapped` is
/// true. `kind` is the geometry kind (cylinder/disk) the decoder actually
/// produced; it lets a caller validate the target surface's kind explicitly,
/// without inferring detector geometry from surface count.
/// `error` is the typed host decode failure; other fields are meaningful only
/// when ok() is true (except `layer`, which may identify an InvalidLayerMapping
/// failure).
struct SurfaceMeasurementDecodeResult {
  o2::itsmft::tracking::SurfaceMeasurement measurement{};
  o2::itsmft::tracking::SurfaceKind kind{o2::itsmft::tracking::SurfaceKind::Cylinder};
  int layer{-1};
  bool layerMapped{false};
  o2::itsmft::tracking::ClusterDecodeError error{o2::itsmft::tracking::ClusterDecodeError::None};

  bool ok() const noexcept { return error == o2::itsmft::tracking::ClusterDecodeError::None; }
};

/// Decode a compact cluster through the detector geometry singleton exactly
/// once, then map the discovered detector-local layer to a global SurfaceId
/// via `layerToSurface`. This avoids the duplicate geometry/layer lookup that
/// a caller would otherwise need in order to know the SurfaceId before
/// decoding (see the single-surface overload above, preserved for callers
/// that already know the target surface).
template <o2::detectors::DetID::ID DetId>
SurfaceMeasurementDecodeResult loadClusterSurfaceMeasurement(
  const CompClusterExt& c,
  o2::itsmft::tracking::BoundedPatternCursor& patterns,
  const TopologyDictionary* dict,
  gsl::span<const o2::itsmft::tracking::SurfaceId> layerToSurface,
  o2::itsmft::tracking::ClusterSourceId source,
  uint32_t externalClusterIndex,
  uint32_t sourceROF,
  bool applySysErrors = true);

/// Convert compact clusters to 3D spacepoints.
/// \tparam DetId o2::detectors::DetID::ITS or DetID::MFT
template <o2::detectors::DetID::ID DetId>
void convertCompactClusters(gsl::span<const CompClusterExt> clusters,
                            gsl::span<const unsigned char>::iterator& pattIt,
                            std::vector<o2::BaseCluster<float>>& output,
                            const TopologyDictionary* dict);

template <class iterator, typename T>
o2::math_utils::Point3D<T> extractClusterData(const CompClusterExt& c, iterator& iter, const TopologyDictionary* dict, T& sig2Row, T& sig2Col, unsigned int* clusterSize = nullptr, o2::itsmft::tracking::ClusterShape* clusterShape = nullptr)
{
  auto pattID = c.getPatternID();
  sig2Row = DefClusError2Row;
  sig2Col = DefClusError2Col; // Dummy COG errors (about half pixel size)
  const auto setShape = [clusterSize, clusterShape](const ClusterPattern& patt, unsigned int nPixels) {
    if (clusterSize != nullptr) {
      *clusterSize = nPixels;
    }
    if (clusterShape != nullptr) {
      *clusterShape = o2::itsmft::tracking::ClusterShape{
        nPixels, static_cast<uint16_t>(patt.getRowSpan()), static_cast<uint16_t>(patt.getColumnSpan())};
    }
  };
  if (pattID != CompCluster::InvalidPatternID) {
    sig2Row = dict->getErr2X(pattID);
    sig2Col = dict->getErr2Z(pattID);
    if (!dict->isGroup(pattID)) {
      setShape(dict->getPattern(pattID), dict->getNpixels(pattID));
      return dict->getClusterCoordinates<T>(c);
    }
    ClusterPattern patt(iter);
    setShape(patt, patt.getNPixels());
    return dict->getClusterCoordinates<T>(c, patt);
  }
  ClusterPattern patt(iter);
  setShape(patt, patt.getNPixels());
  return dict->getClusterCoordinates<T>(c, patt, false);
}

template <typename T>
struct ClusterDataDecodeResult {
  o2::math_utils::Point3D<T> coordinates{};
  T sig2Row{DefClusError2Row};
  T sig2Col{DefClusError2Col};
  o2::itsmft::tracking::ClusterShape shape{};
  o2::itsmft::tracking::ClusterDecodeError error{o2::itsmft::tracking::ClusterDecodeError::None};

  bool ok() const noexcept { return error == o2::itsmft::tracking::ClusterDecodeError::None; }
};

// Bounded counterpart used by normalized loading. Legacy iterator-based
// helpers remain unchanged for existing workflows. Explicit and grouped
// patterns are acquired only after BoundedPatternCursor has proved that the
// complete encoded pattern is present.
template <typename T = float>
ClusterDataDecodeResult<T> extractClusterDataBounded(
  const CompClusterExt& c,
  o2::itsmft::tracking::BoundedPatternCursor& patterns,
  const TopologyDictionary* dict)
{
  ClusterDataDecodeResult<T> result;
  if (dict == nullptr) {
    result.error = o2::itsmft::tracking::ClusterDecodeError::MissingDictionary;
    return result;
  }

  const auto pattID = c.getPatternID();
  if (pattID != CompCluster::InvalidPatternID) {
    if (pattID >= dict->getSize()) {
      result.error = o2::itsmft::tracking::ClusterDecodeError::InvalidPatternId;
      return result;
    }
    result.sig2Row = dict->getErr2X(pattID);
    result.sig2Col = dict->getErr2Z(pattID);
    if (!dict->isGroup(pattID)) {
      const auto& pattern = dict->getPattern(pattID);
      result.shape = o2::itsmft::tracking::ClusterShape{
        static_cast<uint32_t>(dict->getNpixels(pattID)),
        static_cast<uint16_t>(pattern.getRowSpan()),
        static_cast<uint16_t>(pattern.getColumnSpan())};
      result.coordinates = dict->getClusterCoordinates<T>(c);
      return result;
    }
  }

  ClusterPattern pattern;
  result.error = patterns.acquirePattern(pattern);
  if (!result.ok()) {
    return result;
  }
  result.shape = o2::itsmft::tracking::ClusterShape{
    static_cast<uint32_t>(pattern.getNPixels()),
    static_cast<uint16_t>(pattern.getRowSpan()),
    static_cast<uint16_t>(pattern.getColumnSpan())};
  result.coordinates = dict->getClusterCoordinates<T>(c, pattern, pattID != CompCluster::InvalidPatternID);
  return result;
}

// same method returning coordinates as an array (suitable for the TGeoMatrix)
template <class iterator, typename T>
std::array<T, 3> extractClusterDataA(const CompClusterExt& c, iterator& iter, const TopologyDictionary* dict, T& sig2Row, T& sig2Col)
{
  auto pattID = c.getPatternID();
  sig2Row = DefClusError2Row;
  sig2Col = DefClusError2Col; // Dummy COG errors (about half pixel size)
  if (pattID != CompCluster::InvalidPatternID) {
    sig2Row = dict->getErr2X(pattID);
    sig2Col = dict->getErr2Z(pattID);
    if (!dict->isGroup(pattID)) {
      return dict->getClusterCoordinatesA<T>(c);
    }
    ClusterPattern patt(iter);
    return dict->getClusterCoordinatesA<T>(c, patt);
  }
  ClusterPattern patt(iter);
  return dict->getClusterCoordinatesA<T>(c, patt, false);
}

} // namespace o2::itsmft::ioutils

#endif /* ALICEO2_ITSMFT_TRACKING_IOUTILS_H_ */

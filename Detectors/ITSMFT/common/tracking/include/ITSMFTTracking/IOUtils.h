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

/// Whether to apply cluster-decoding systematic errors for `DetId`. ITS
/// common-CA is an explicit no-op; MFT reads its live tracker configuration.
/// The free function keeps this detector-generic boundary independently
/// callable.
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

/// Adds configured systematic-error corrections to `sigma2Row`/`sigma2Col`;
/// the ITS specialization is an explicit no-op.
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

/// Decode once through detector geometry, then map the discovered local layer
/// to a global SurfaceId through `layerToSurface`.
template <o2::detectors::DetID::ID DetId>
o2::itsmft::tracking::SurfaceMeasurementDecodeResult loadClusterSurfaceMeasurement(
  const CompClusterExt& c,
  o2::itsmft::tracking::BoundedPatternCursor& patterns,
  const TopologyDictionary* dict,
  gsl::span<const o2::itsmft::tracking::SurfaceId> layerToSurface,
  o2::itsmft::tracking::ClusterSourceId source,
  uint32_t externalClusterIndex,
  uint32_t sourceROF,
  bool applySysErrors);

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

// Bounded normalized-loading counterpart. Pattern bytes are acquired only
// after BoundedPatternCursor has proved that the complete encoding is present.
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

// Array-valued coordinate counterpart for TGeoMatrix callers.
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

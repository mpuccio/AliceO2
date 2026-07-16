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
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "MathUtils/Cartesian.h"

namespace o2::its
{
struct TrackingFrameInfo;
}

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
} // namespace detail

constexpr float DefClusErrorRow = o2::itsmft::SegmentationAlpide::PitchRow * 0.5f;
constexpr float DefClusErrorCol = o2::itsmft::SegmentationAlpide::PitchCol * 0.5f;
constexpr float DefClusError2Row = DefClusErrorRow * DefClusErrorRow;
constexpr float DefClusError2Col = DefClusErrorCol * DefClusErrorCol;

void fillMatrixCache(o2::detectors::DetID::ID detId);
int getClusterLayer(o2::detectors::DetID::ID detId, const CompClusterExt& cluster);

/// Decode a compact cluster into layer, size, and a TrackingFrameInfo (global + local frame).
template <o2::detectors::DetID::ID DetId>
void loadClusterTrackingFrameInfo(const CompClusterExt& c,
                                  gsl::span<const unsigned char>::iterator& pattIt,
                                  const TopologyDictionary* dict,
                                  int& layer,
                                  unsigned int& clusterSize,
                                  o2::its::TrackingFrameInfo& tfInfo,
                                  bool applySysErrors = true);

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
/// layer. `measurement.surface` is only meaningful when `layerMapped` is true.
struct SurfaceMeasurementDecodeResult {
  o2::itsmft::tracking::SurfaceMeasurement measurement{};
  int layer{-1};
  bool layerMapped{false};
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
  gsl::span<const unsigned char>::iterator& pattIt,
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

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
/// \file IOUtils.cxx
/// \brief Shared cluster I/O utilities for ITS and MFT (based on ITStracking/IOUtils.cxx)
///

#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/DecodedCluster.h"
#include "ITSMFTTracking/SurfaceMeasurementAdapters.h"
#include "ITSMFTTracking/TrackingFrameInfoAdapters.h"
#include "ITStracking/Cluster.h"

#include "Framework/Logger.h"
#include "ITSBase/GeometryTGeo.h"
#include "MFTBase/GeometryTGeo.h"
#include "MathUtils/Utils.h"

namespace
{
constexpr int PrimaryVertexLayerId{-1};
constexpr int EventLabelsSeparator{-1};

template <o2::detectors::DetID::ID DetId>
bool shouldApplySysErrors()
{
  const auto& conf = o2::itsmft::tracking::TrackerParamRef<DetId>::get();
  for (int il = 0; il < o2::itsmft::tracking::TrackerParamRef<DetId>::nLayers(); il++) {
    if constexpr (DetId == o2::detectors::DetID::MFT) {
      if (conf.sysErr2Row[il] > 0.f || conf.sysErr2Col[il] > 0.f) {
        return true;
      }
    } else {
      if (conf.sysErrY2[il] > 0.f || conf.sysErrZ2[il] > 0.f) {
        return true;
      }
    }
  }
  return false;
}

template <o2::detectors::DetID::ID DetId>
void addSysErrors(int layerId, float& sigma2Row, float& sigma2Col)
{
  const auto& conf = o2::itsmft::tracking::TrackerParamRef<DetId>::get();
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    sigma2Row += conf.sysErr2Row[layerId];
    sigma2Col += conf.sysErr2Col[layerId];
  } else {
    sigma2Row += conf.sysErrY2[layerId];
    sigma2Col += conf.sysErrZ2[layerId];
  }
}

template <o2::detectors::DetID::ID DetId, typename GeomT>
o2::itsmft::tracking::DecodedCluster decodeCluster(GeomT* geom,
                                                   const o2::itsmft::CompClusterExt& c,
                                                   gsl::span<const unsigned char>::iterator& pattIt,
                                                   const o2::itsmft::TopologyDictionary* dict,
                                                   bool applySysErrors)
{
  const auto sensorID = c.getSensorID();
  if (!o2::itsmft::ioutils::detail::isSensorInGeometry(sensorID, geom->getSize())) {
    LOGP(fatal, "Cluster sensorID {} is out of geometry range [0, {})", sensorID, geom->getSize());
  }
  const int layer = geom->getLayer(sensorID);
  if (!o2::itsmft::ioutils::detail::isLayerInDetector(layer, o2::itsmft::tracking::TrackerParamRef<DetId>::nLayers())) {
    LOGP(fatal, "Cluster sensorID {} maps to invalid layer {} (expected [0, {}))",
         sensorID, layer, o2::itsmft::tracking::TrackerParamRef<DetId>::nLayers());
  }
  const auto pattID = c.getPatternID();
  if (pattID != o2::itsmft::CompCluster::InvalidPatternID && dict != nullptr &&
      (pattID < 0 || pattID >= dict->getSize())) {
    LOGP(fatal, "Cluster patternID {} is out of dictionary range [0, {})", pattID, dict->getSize());
  }

  float sigma2Row{0.f};
  float sigma2Col{0.f};
  o2::itsmft::tracking::ClusterShape shape{};
  const auto locXYZ = o2::itsmft::ioutils::extractClusterData(c, pattIt, dict, sigma2Row, sigma2Col, nullptr, &shape);

  if (applySysErrors && shouldApplySysErrors<DetId>()) {
    addSysErrors<DetId>(layer, sigma2Row, sigma2Col);
  }

  if constexpr (DetId == o2::detectors::DetID::ITS) {
    const auto trkXYZ = geom->getMatrixT2L(sensorID) ^ locXYZ;
    const auto gloXYZ = geom->getMatrixL2G(sensorID) * locXYZ;
    return o2::itsmft::tracking::DecodedCluster{
      {gloXYZ.x(), gloXYZ.y(), gloXYZ.z()},
      {trkXYZ.x(), trkXYZ.y(), trkXYZ.z(), geom->getSensorRefAlpha(sensorID)},
      {sigma2Row, 0.f, sigma2Col},
      shape,
      static_cast<uint32_t>(sensorID),
      layer};
  } else {
    if (!geom->getCacheL2G().isFilled() || geom->getCacheL2G().getSize() <= sensorID) {
      LOGP(fatal, "MFT L2G matrix cache unavailable for sensorID {} (filled={}, size={})",
           sensorID, geom->getCacheL2G().isFilled(), geom->getCacheL2G().getSize());
    }
    const auto gloXYZ = geom->getMatrixL2G(sensorID) * locXYZ;
    // ALPIDE row (local X) -> global X, column (local Z) -> global Y
    return o2::itsmft::tracking::DecodedCluster{
      {gloXYZ.x(), gloXYZ.y(), gloXYZ.z()},
      {},
      {sigma2Row, 0.f, sigma2Col},
      shape,
      static_cast<uint32_t>(sensorID),
      layer};
  }
}

struct DecodedClusterResult {
  o2::itsmft::tracking::DecodedCluster decoded{};
  o2::itsmft::tracking::ClusterDecodeError error{o2::itsmft::tracking::ClusterDecodeError::None};

  bool ok() const noexcept { return error == o2::itsmft::tracking::ClusterDecodeError::None; }
};

// Safe normalized-loading boundary. The iterator-based decodeCluster above
// remains for unchanged legacy APIs; this overload performs every
// malformed-input check before geometry or ClusterPattern data are used.
template <o2::detectors::DetID::ID DetId, typename GeomT>
DecodedClusterResult decodeClusterBounded(GeomT* geom,
                                          const o2::itsmft::CompClusterExt& c,
                                          o2::itsmft::tracking::BoundedPatternCursor& patterns,
                                          const o2::itsmft::TopologyDictionary* dict,
                                          bool applySysErrors)
{
  using o2::itsmft::tracking::ClusterDecodeError;
  DecodedClusterResult result;
  if (dict == nullptr) {
    result.error = ClusterDecodeError::MissingDictionary;
    return result;
  }
  if (geom == nullptr) {
    result.error = ClusterDecodeError::GeometryUnavailable;
    return result;
  }

  const auto sensorID = c.getSensorID();
  if (!o2::itsmft::ioutils::detail::isSensorInGeometry(sensorID, geom->getSize())) {
    result.error = ClusterDecodeError::InvalidSensor;
    return result;
  }
  const int layer = geom->getLayer(sensorID);
  if (!o2::itsmft::ioutils::detail::isLayerInDetector(layer, o2::itsmft::tracking::TrackerParamRef<DetId>::nLayers())) {
    result.error = ClusterDecodeError::InvalidLayer;
    return result;
  }

  const auto clusterData = o2::itsmft::ioutils::extractClusterDataBounded(c, patterns, dict);
  if (!clusterData.ok()) {
    result.error = clusterData.error;
    return result;
  }
  float sigma2Row = clusterData.sig2Row;
  float sigma2Col = clusterData.sig2Col;
  if (applySysErrors && shouldApplySysErrors<DetId>()) {
    addSysErrors<DetId>(layer, sigma2Row, sigma2Col);
  }

  if constexpr (DetId == o2::detectors::DetID::ITS) {
    const auto trkXYZ = geom->getMatrixT2L(sensorID) ^ clusterData.coordinates;
    const auto gloXYZ = geom->getMatrixL2G(sensorID) * clusterData.coordinates;
    result.decoded = o2::itsmft::tracking::DecodedCluster{
      {gloXYZ.x(), gloXYZ.y(), gloXYZ.z()},
      {trkXYZ.x(), trkXYZ.y(), trkXYZ.z(), geom->getSensorRefAlpha(sensorID)},
      {sigma2Row, 0.f, sigma2Col},
      clusterData.shape,
      static_cast<uint32_t>(sensorID),
      layer};
  } else {
    if (!geom->getCacheL2G().isFilled() || geom->getCacheL2G().getSize() <= sensorID) {
      result.error = ClusterDecodeError::GeometryUnavailable;
      return result;
    }
    const auto gloXYZ = geom->getMatrixL2G(sensorID) * clusterData.coordinates;
    result.decoded = o2::itsmft::tracking::DecodedCluster{
      {gloXYZ.x(), gloXYZ.y(), gloXYZ.z()},
      {},
      {sigma2Row, 0.f, sigma2Col},
      clusterData.shape,
      static_cast<uint32_t>(sensorID),
      layer};
  }
  return result;
}

template <o2::detectors::DetID::ID DetId, typename GeomT>
void fillOutputClusters(GeomT* geom,
                        gsl::span<const o2::itsmft::CompClusterExt> clusters,
                        gsl::span<const unsigned char>::iterator& pattIt,
                        std::vector<o2::BaseCluster<float>>& output,
                        const o2::itsmft::TopologyDictionary* dict,
                        bool applyMisalignment)
{
  for (const auto& c : clusters) {
    float sigma2Row{0.f}, sigma2Col{0.f};
    const float sigmaRowCol{0.f}; // row-column covariance (unused)
    const auto locXYZ = o2::itsmft::ioutils::extractClusterData(c, pattIt, dict, sigma2Row, sigma2Col);
    if (applyMisalignment) {
      addSysErrors<DetId>(geom->getLayer(c.getSensorID()), sigma2Row, sigma2Col);
    }
    o2::math_utils::Point3D<float> outXYZ{};
    if constexpr (DetId == o2::detectors::DetID::ITS) {
      outXYZ = geom->getMatrixT2L(c.getSensorID()) ^ locXYZ;
    } else {
      outXYZ = geom->getMatrixL2G(c.getSensorID()) * locXYZ;
    }
    auto& cl3d = output.emplace_back(c.getSensorID(), outXYZ);
    cl3d.setErrors(sigma2Row, sigma2Col, sigmaRowCol);
  }
}

} // namespace

namespace o2::itsmft::ioutils
{

void fillMatrixCache(o2::detectors::DetID::ID detId)
{
  const auto mask = o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G);
  if (detId == o2::detectors::DetID::ITS) {
    o2::its::GeometryTGeo::Instance()->fillMatrixCache(mask);
  } else if (detId == o2::detectors::DetID::MFT) {
    o2::mft::GeometryTGeo::Instance()->fillMatrixCache(mask);
  } else {
    LOGP(fatal, "Unsupported detector id {} in fillMatrixCache", static_cast<int>(detId));
  }
}

int getClusterLayer(o2::detectors::DetID::ID detId, const CompClusterExt& cluster)
{
  if (detId == o2::detectors::DetID::ITS) {
    return o2::its::GeometryTGeo::Instance()->getLayer(cluster.getSensorID());
  }
  if (detId == o2::detectors::DetID::MFT) {
    return o2::mft::GeometryTGeo::Instance()->getLayer(cluster.getSensorID());
  }
  LOGP(fatal, "Unsupported detector id {} in getClusterLayer", static_cast<int>(detId));
  return -1;
}

template <o2::detectors::DetID::ID DetId>
void loadClusterTrackingFrameInfo(const CompClusterExt& c,
                                  gsl::span<const unsigned char>::iterator& pattIt,
                                  const TopologyDictionary* dict,
                                  int& layer,
                                  unsigned int& clusterSize,
                                  o2::its::TrackingFrameInfo& tfInfo,
                                  bool applySysErrors)
{
  o2::itsmft::tracking::DecodedCluster decoded;
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    decoded = decodeCluster<DetId>(o2::its::GeometryTGeo::Instance(), c, pattIt, dict, applySysErrors);
  } else {
    decoded = decodeCluster<DetId>(o2::mft::GeometryTGeo::Instance(), c, pattIt, dict, applySysErrors);
  }
  layer = decoded.layer;
  clusterSize = decoded.shape.nPixels;
  tfInfo = o2::itsmft::tracking::makeTrackingFrameInfo<DetId>(decoded);
}

template void loadClusterTrackingFrameInfo<o2::detectors::DetID::ITS>(const CompClusterExt& c,
                                                                      gsl::span<const unsigned char>::iterator& pattIt,
                                                                      const TopologyDictionary* dict,
                                                                      int& layer,
                                                                      unsigned int& clusterSize,
                                                                      o2::its::TrackingFrameInfo& tfInfo,
                                                                      bool applySysErrors);

template void loadClusterTrackingFrameInfo<o2::detectors::DetID::MFT>(const CompClusterExt& c,
                                                                      gsl::span<const unsigned char>::iterator& pattIt,
                                                                      const TopologyDictionary* dict,
                                                                      int& layer,
                                                                      unsigned int& clusterSize,
                                                                      o2::its::TrackingFrameInfo& tfInfo,
                                                                      bool applySysErrors);

template <o2::detectors::DetID::ID DetId>
o2::itsmft::tracking::SurfaceMeasurement loadClusterSurfaceMeasurement(
  const CompClusterExt& c,
  gsl::span<const unsigned char>::iterator& pattIt,
  const TopologyDictionary* dict,
  o2::itsmft::tracking::ClusterSourceId source,
  uint32_t externalClusterIndex,
  o2::itsmft::tracking::SurfaceId surface,
  uint32_t sourceROF,
  bool applySysErrors)
{
  o2::itsmft::tracking::DecodedCluster decoded;
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    decoded = decodeCluster<DetId>(o2::its::GeometryTGeo::Instance(), c, pattIt, dict, applySysErrors);
    return o2::itsmft::tracking::makeCylinderSurfaceMeasurement(
      decoded, {DetId, decoded.sensor}, surface, {source, externalClusterIndex}, sourceROF);
  } else {
    decoded = decodeCluster<DetId>(o2::mft::GeometryTGeo::Instance(), c, pattIt, dict, applySysErrors);
    return o2::itsmft::tracking::makeDiskSurfaceMeasurement(
      decoded, {DetId, decoded.sensor}, surface, {source, externalClusterIndex}, sourceROF);
  }
}

template o2::itsmft::tracking::SurfaceMeasurement loadClusterSurfaceMeasurement<o2::detectors::DetID::ITS>(
  const CompClusterExt&, gsl::span<const unsigned char>::iterator&, const TopologyDictionary*,
  o2::itsmft::tracking::ClusterSourceId, uint32_t, o2::itsmft::tracking::SurfaceId, uint32_t, bool);

template o2::itsmft::tracking::SurfaceMeasurement loadClusterSurfaceMeasurement<o2::detectors::DetID::MFT>(
  const CompClusterExt&, gsl::span<const unsigned char>::iterator&, const TopologyDictionary*,
  o2::itsmft::tracking::ClusterSourceId, uint32_t, o2::itsmft::tracking::SurfaceId, uint32_t, bool);

template <o2::detectors::DetID::ID DetId>
SurfaceMeasurementDecodeResult loadClusterSurfaceMeasurement(
  const CompClusterExt& c,
  o2::itsmft::tracking::BoundedPatternCursor& patterns,
  const TopologyDictionary* dict,
  gsl::span<const o2::itsmft::tracking::SurfaceId> layerToSurface,
  o2::itsmft::tracking::ClusterSourceId source,
  uint32_t externalClusterIndex,
  uint32_t sourceROF,
  bool applySysErrors)
{
  DecodedClusterResult decodedResult;
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    decodedResult = decodeClusterBounded<DetId>(o2::its::GeometryTGeo::Instance(), c, patterns, dict, applySysErrors);
  } else {
    decodedResult = decodeClusterBounded<DetId>(o2::mft::GeometryTGeo::Instance(), c, patterns, dict, applySysErrors);
  }

  SurfaceMeasurementDecodeResult result;
  if (!decodedResult.ok()) {
    result.error = decodedResult.error;
    return result;
  }
  const auto& decoded = decodedResult.decoded;
  result.layer = decoded.layer;
  if (decoded.layer < 0 || static_cast<size_t>(decoded.layer) >= layerToSurface.size()) {
    result.error = o2::itsmft::tracking::ClusterDecodeError::InvalidLayerMapping;
    return result;
  }
  result.layerMapped = true;

  const auto surface = layerToSurface[decoded.layer];
  const o2::itsmft::tracking::DetectorSensorId sensor{DetId, decoded.sensor};
  const o2::itsmft::tracking::ClusterRef cluster{source, externalClusterIndex};
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    result.kind = o2::itsmft::tracking::SurfaceKind::Cylinder;
    result.measurement = o2::itsmft::tracking::makeCylinderSurfaceMeasurement(decoded, sensor, surface, cluster, sourceROF);
  } else {
    result.kind = o2::itsmft::tracking::SurfaceKind::Disk;
    result.measurement = o2::itsmft::tracking::makeDiskSurfaceMeasurement(decoded, sensor, surface, cluster, sourceROF);
  }
  return result;
}

template SurfaceMeasurementDecodeResult loadClusterSurfaceMeasurement<o2::detectors::DetID::ITS>(
  const CompClusterExt&, o2::itsmft::tracking::BoundedPatternCursor&, const TopologyDictionary*,
  gsl::span<const o2::itsmft::tracking::SurfaceId>, o2::itsmft::tracking::ClusterSourceId, uint32_t, uint32_t, bool);

template SurfaceMeasurementDecodeResult loadClusterSurfaceMeasurement<o2::detectors::DetID::MFT>(
  const CompClusterExt&, o2::itsmft::tracking::BoundedPatternCursor&, const TopologyDictionary*,
  gsl::span<const o2::itsmft::tracking::SurfaceId>, o2::itsmft::tracking::ClusterSourceId, uint32_t, uint32_t, bool);

template <o2::detectors::DetID::ID DetId>
void convertCompactClusters(gsl::span<const CompClusterExt> clusters,
                            gsl::span<const unsigned char>::iterator& pattIt,
                            std::vector<o2::BaseCluster<float>>& output,
                            const TopologyDictionary* dict)
{
  const auto mask = o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G);

  const bool applyMisalignment = shouldApplySysErrors<DetId>();

  if constexpr (DetId == o2::detectors::DetID::ITS) {
    auto* geom = o2::its::GeometryTGeo::Instance();
    geom->fillMatrixCache(mask);
    fillOutputClusters<DetId>(geom, clusters, pattIt, output, dict, applyMisalignment);
  } else {
    auto* geom = o2::mft::GeometryTGeo::Instance();
    geom->fillMatrixCache(mask);
    fillOutputClusters<DetId>(geom, clusters, pattIt, output, dict, applyMisalignment);
  }
}

template void convertCompactClusters<o2::detectors::DetID::ITS>(gsl::span<const CompClusterExt> clusters,
                                                                gsl::span<const unsigned char>::iterator& pattIt,
                                                                std::vector<o2::BaseCluster<float>>& output,
                                                                const TopologyDictionary* dict);

template void convertCompactClusters<o2::detectors::DetID::MFT>(gsl::span<const CompClusterExt> clusters,
                                                                gsl::span<const unsigned char>::iterator& pattIt,
                                                                std::vector<o2::BaseCluster<float>>& output,
                                                                const TopologyDictionary* dict);

} // namespace o2::itsmft::ioutils

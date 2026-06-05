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
#include "ITSMFTTracking/TrackingConfigParam.h"
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
  const auto& conf = o2::itsmft::TrackerParamConfig<DetId>::Instance();
  for (int il = 0; il < o2::itsmft::TrackerParamConfig<DetId>::getNLayers(); il++) {
    if (conf.sysErr2Row[il] > 0.f || conf.sysErr2Col[il] > 0.f) {
      return true;
    }
  }
  return false;
}

template <o2::detectors::DetID::ID DetId, typename GeomT>
void loadClusterTrackingFrameInfoImpl(GeomT* geom,
                                      const o2::itsmft::CompClusterExt& c,
                                      gsl::span<const unsigned char>::iterator& pattIt,
                                      const o2::itsmft::TopologyDictionary* dict,
                                      int& layer,
                                      unsigned int& clusterSize,
                                      o2::its::TrackingFrameInfo& tfInfo,
                                      bool applySysErrors)
{
  const auto sensorID = c.getSensorID();
  layer = geom->getLayer(sensorID);
  clusterSize = o2::itsmft::ioutils::extractClusterSize(c, pattIt, dict);

  float sigma2Row{0.f};
  float sigma2Col{0.f};
  const auto locXYZ = o2::itsmft::ioutils::extractClusterData(c, pattIt, dict, sigma2Row, sigma2Col);
  if (applySysErrors && shouldApplySysErrors<DetId>()) {
    const auto layerId = geom->getLayer(sensorID);
    const auto& conf = o2::itsmft::TrackerParamConfig<DetId>::Instance();
    sigma2Row += conf.sysErr2Row[layerId];
    sigma2Col += conf.sysErr2Col[layerId];
  }

  if constexpr (DetId == o2::detectors::DetID::ITS) {
    const auto trkXYZ = geom->getMatrixT2L(sensorID) ^ locXYZ;
    const auto gloXYZ = geom->getMatrixL2G(sensorID) * locXYZ;
    tfInfo = o2::its::TrackingFrameInfo{
      gloXYZ.x(), gloXYZ.y(), gloXYZ.z(), trkXYZ.x(), geom->getSensorRefAlpha(sensorID),
      std::array<float, 2>{trkXYZ.y(), trkXYZ.z()},
      std::array<float, 3>{sigma2Row, 0.f, sigma2Col}};
  } else {
    const auto gloXYZ = geom->getMatrixL2G(sensorID) * locXYZ;
    // ALPIDE row (local X) -> global X, column (local Z) -> global Y
    tfInfo = o2::its::TrackingFrameInfo{
      gloXYZ.x(), gloXYZ.y(), gloXYZ.z(), gloXYZ.x(), 0.f,
      std::array<float, 2>{gloXYZ.y(), gloXYZ.z()},
      std::array<float, 3>{sigma2Row, 0.f, sigma2Col}};
  }
}

template <o2::detectors::DetID::ID DetId, typename GeomT>
void fillOutputClusters(GeomT* geom,
                        gsl::span<const o2::itsmft::CompClusterExt> clusters,
                        gsl::span<const unsigned char>::iterator& pattIt,
                        std::vector<o2::BaseCluster<float>>& output,
                        const o2::itsmft::TopologyDictionary* dict,
                        const o2::itsmft::TrackerParamConfig<DetId>& conf,
                        bool applyMisalignment)
{
  for (const auto& c : clusters) {
    float sigma2Row{0.f}, sigma2Col{0.f};
    const float sigmaRowCol{0.f}; // row-column covariance (unused)
    auto locXYZ = o2::itsmft::ioutils::extractClusterData(c, pattIt, dict, sigma2Row, sigma2Col);
    if (applyMisalignment) {
      const auto layerId = geom->getLayer(c.getSensorID());
      sigma2Row += conf.sysErr2Row[layerId];
      sigma2Col += conf.sysErr2Col[layerId];
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
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    loadClusterTrackingFrameInfoImpl<DetId>(o2::its::GeometryTGeo::Instance(), c, pattIt, dict, layer, clusterSize, tfInfo, applySysErrors);
  } else {
    loadClusterTrackingFrameInfoImpl<DetId>(o2::mft::GeometryTGeo::Instance(), c, pattIt, dict, layer, clusterSize, tfInfo, applySysErrors);
  }
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
void convertCompactClusters(gsl::span<const CompClusterExt> clusters,
                            gsl::span<const unsigned char>::iterator& pattIt,
                            std::vector<o2::BaseCluster<float>>& output,
                            const TopologyDictionary* dict)
{
  const auto mask = o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G);

  const auto& conf = TrackerParamConfig<DetId>::Instance();
  bool applyMisalignment = false;
  for (int il = 0; il < TrackerParamConfig<DetId>::getNLayers(); il++) {
    if (conf.sysErr2Row[il] > 0.f || conf.sysErr2Col[il] > 0.f) {
      applyMisalignment = true;
      break;
    }
  }

  if constexpr (DetId == o2::detectors::DetID::ITS) {
    auto* geom = o2::its::GeometryTGeo::Instance();
    geom->fillMatrixCache(mask);
    fillOutputClusters<DetId>(geom, clusters, pattIt, output, dict, conf, applyMisalignment);
  } else {
    auto* geom = o2::mft::GeometryTGeo::Instance();
    geom->fillMatrixCache(mask);
    fillOutputClusters<DetId>(geom, clusters, pattIt, output, dict, conf, applyMisalignment);
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

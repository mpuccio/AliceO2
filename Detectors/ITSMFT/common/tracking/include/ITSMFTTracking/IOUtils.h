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
#include "MathUtils/Cartesian.h"

namespace o2::its
{
struct TrackingFrameInfo;
}

namespace o2::itsmft::ioutils
{

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

/// Convert compact clusters to 3D spacepoints.
/// \tparam DetId o2::detectors::DetID::ITS or DetID::MFT
template <o2::detectors::DetID::ID DetId>
void convertCompactClusters(gsl::span<const CompClusterExt> clusters,
                            gsl::span<const unsigned char>::iterator& pattIt,
                            std::vector<o2::BaseCluster<float>>& output,
                            const TopologyDictionary* dict);

template <class iterator, typename T>
o2::math_utils::Point3D<T> extractClusterData(const CompClusterExt& c, iterator& iter, const TopologyDictionary* dict, T& sig2Row, T& sig2Col, unsigned int* clusterSize = nullptr)
{
  auto pattID = c.getPatternID();
  sig2Row = DefClusError2Row;
  sig2Col = DefClusError2Col; // Dummy COG errors (about half pixel size)
  if (pattID != CompCluster::InvalidPatternID) {
    sig2Row = dict->getErr2X(pattID);
    sig2Col = dict->getErr2Z(pattID);
    if (!dict->isGroup(pattID)) {
      if (clusterSize != nullptr) {
        *clusterSize = dict->getNpixels(pattID);
      }
      return dict->getClusterCoordinates<T>(c);
    }
    ClusterPattern patt(iter);
    if (clusterSize != nullptr) {
      *clusterSize = patt.getNPixels();
    }
    return dict->getClusterCoordinates<T>(c, patt);
  }
  ClusterPattern patt(iter);
  if (clusterSize != nullptr) {
    *clusterSize = patt.getNPixels();
  }
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

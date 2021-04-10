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
/// \file TrackerTraitsCPU.cxx
/// \brief
///

#include "ITStracking/TrackerTraitsCPU.h"

#include "CommonConstants/MathConstants.h"
#include "ITStracking/Cell.h"
#include "ITStracking/Constants.h"
#include "ITStracking/IndexTableUtils.h"
#include "ITStracking/Tracklet.h"
#include <fmt/format.h>
#include "ReconstructionDataFormats/Track.h"
#include <cassert>
#include <iostream>

#include "GPUCommonMath.h"

namespace o2
{
namespace its
{

void TrackerTraitsCPU::computeLayerTracklets()
{
  TimeFrame* tf = mTimeFrame;

  for (int rof0{0}; rof0 < tf->getNrof(); ++rof0) {
    gsl::span<const float3> primaryVertices = tf->getPrimaryVertices(rof0);
    for (int iLayer{0}; iLayer < mTrkParams.TrackletsPerRoad(); ++iLayer) {
      gsl::span<const Cluster> layer0 = tf->getClustersOnLayer(rof0, iLayer);
      if (layer0.empty()) {
        return;
      }

      const int currentLayerClustersNum{static_cast<int>(layer0.size())};

      for (int iCluster{0}; iCluster < currentLayerClustersNum; ++iCluster) {
        const Cluster& currentCluster{layer0.at(iCluster)};

        if (tf->isClusterUsed(iLayer, currentCluster.clusterId)) {
          continue;
        }

        for (auto& primaryVertex : primaryVertices) {
          const float tanLambda{(currentCluster.zCoordinate - primaryVertex.z) / currentCluster.radius};

          const float zAtRmin{tanLambda * (tf->getMinR(iLayer + 1) -
                                           currentCluster.radius) +
                              currentCluster.zCoordinate};
          const float zAtRmax{tanLambda * (tf->getMaxR(iLayer + 1) -
                                           currentCluster.radius) +
                              currentCluster.zCoordinate};

          const int4 selectedBinsRect{getBinsRect(currentCluster, iLayer, zAtRmin, zAtRmax,
                                                  mTrkParams.TrackletMaxDeltaZ[iLayer], mTrkParams.TrackletMaxDeltaPhi)};

          if (selectedBinsRect.x == 0 && selectedBinsRect.y == 0 && selectedBinsRect.z == 0 && selectedBinsRect.w == 0) {
            continue;
          }

          int phiBinsNum{selectedBinsRect.w - selectedBinsRect.y + 1};

          if (phiBinsNum < 0) {
            phiBinsNum += mTrkParams.PhiBins;
          }

          int minRof = (rof0 > 0) ? rof0 - 1 : 0;
          int maxRof = (rof0 == tf->getNrof() - 1) ? rof0 : rof0 + 1;
          for (int rof1{minRof}; rof1 < maxRof; ++rof1) {
            gsl::span<const Cluster> layer1 = tf->getClustersOnLayer(rof1, iLayer + 1);
            if (layer1.empty()) {
              continue;
            }

            for (int iPhiBin{selectedBinsRect.y}, iPhiCount{0}; iPhiCount < phiBinsNum;
                 iPhiBin = (++iPhiBin == tf->mIndexTableUtils.getNphiBins()) ? 0 : iPhiBin, iPhiCount++) {
              const int firstBinIndex{tf->mIndexTableUtils.getBinIndex(selectedBinsRect.x, iPhiBin)};
              const int maxBinIndex{firstBinIndex + selectedBinsRect.z - selectedBinsRect.x + 1};
              const int firstRowClusterIndex = tf->getIndexTables(rof1)[iLayer][firstBinIndex];
              const int maxRowClusterIndex = tf->getIndexTables(rof1)[iLayer][maxBinIndex];

              for (int iNextCluster{firstRowClusterIndex}; iNextCluster < maxRowClusterIndex; ++iNextCluster) {
                if (iNextCluster >= (int)tf->getClusters()[iLayer + 1].size()) {
                  break;
                }
                const Cluster& nextCluster{layer1.at(iNextCluster)};

                if (tf->isClusterUsed(iLayer + 1, nextCluster.clusterId)) {
                  continue;
                }

                const float deltaZ{gpu::GPUCommonMath::Abs(tanLambda * (nextCluster.radius - currentCluster.radius) +
                                                           currentCluster.zCoordinate - nextCluster.zCoordinate)};
                const float deltaPhi{gpu::GPUCommonMath::Abs(currentCluster.phi - nextCluster.phi)};

                if (deltaZ < mTrkParams.TrackletMaxDeltaZ[iLayer] &&
                    (deltaPhi < mTrkParams.TrackletMaxDeltaPhi ||
                     gpu::GPUCommonMath::Abs(deltaPhi - constants::math::TwoPi) < mTrkParams.TrackletMaxDeltaPhi)) {
                  const int currentSortedIndex{tf->getSortedIndex(rof0, iLayer, iCluster)};
                  if (iLayer > 0)
                    if (currentSortedIndex > tf->getTrackletsLookupTable()[iLayer - 1].size())
                      std::cout << "Issue with the tracklet LUT" << std::endl;
                  if (iLayer > 0 &&
                      tf->getTrackletsLookupTable()[iLayer - 1][currentSortedIndex] == constants::its::UnusedIndex) {
                    tf->getTrackletsLookupTable()[iLayer - 1][currentSortedIndex] =
                      tf->getTracklets()[iLayer].size();
                  }

                  tf->getTracklets()[iLayer].emplace_back(currentSortedIndex, tf->getSortedIndex(rof1, iLayer + 1, iNextCluster), currentCluster,
                                                          nextCluster, rof0, rof1);
                }
              }
            }
          }
        }
      }
    }
    // if (iLayer > 0 && iLayer < mTrkParams.TrackletsPerRoad() - 1 &&
    //     tf->getTracklets()[iLayer].size() > tf->getCellsLookupTable()[iLayer - 1].size()) {
    //   std::cout << "**** FATAL: not enough memory in the CellsLookupTable, increase the tracklet memory coefficients ****" << std::endl;
    //   exit(1);
    // }
  }
#ifdef CA_DEBUG
  std::cout << "+++ Number of tracklets per layer: ";
  for (int iLayer{0}; iLayer < mTrkParams.TrackletsPerRoad(); ++iLayer) {
    std::cout << primaryVertexContext->getTracklets()[iLayer].size() << "\t";
  }
#endif
  std::cout << "Number of tracklets " << std::endl;
  for (int iL{0}; iL < 6; ++iL) {
    std::cout << tf->getTracklets()[iL].size() << std::endl;
  }
}

void TrackerTraitsCPU::computeLayerCells()
{
  TimeFrame* tf = mTimeFrame;
  for (int iLayer{0}; iLayer < mTrkParams.CellsPerRoad(); ++iLayer) {

    if (tf->getTracklets()[iLayer + 1].empty() ||
        tf->getTracklets()[iLayer].empty()) {
      return;
    }

    const int currentLayerTrackletsNum{static_cast<int>(tf->getTracklets()[iLayer].size())};

    for (int iTracklet{0}; iTracklet < currentLayerTrackletsNum; ++iTracklet) {

      const Tracklet& currentTracklet{tf->getTracklets()[iLayer][iTracklet]};
      const int nextLayerClusterIndex{currentTracklet.secondClusterIndex};
      const int nextLayerFirstTrackletIndex{
        tf->getTrackletsLookupTable()[iLayer][nextLayerClusterIndex]};

      if (nextLayerFirstTrackletIndex == constants::its::UnusedIndex) {
        continue;
      }

      const Cluster& firstCellCluster{tf->getClusters()[iLayer][currentTracklet.firstClusterIndex]};
      const Cluster& secondCellCluster{
        tf->getClusters()[iLayer + 1][currentTracklet.secondClusterIndex]};
      const float firstCellClusterQuadraticRCoordinate{firstCellCluster.radius * firstCellCluster.radius};
      const float secondCellClusterQuadraticRCoordinate{secondCellCluster.radius *
                                                        secondCellCluster.radius};
      const float3 firstDeltaVector{secondCellCluster.xCoordinate - firstCellCluster.xCoordinate,
                                    secondCellCluster.yCoordinate - firstCellCluster.yCoordinate,
                                    secondCellClusterQuadraticRCoordinate - firstCellClusterQuadraticRCoordinate};
      const int nextLayerTrackletsNum{static_cast<int>(tf->getTracklets()[iLayer + 1].size())};

      for (int iNextLayerTracklet{nextLayerFirstTrackletIndex};
           iNextLayerTracklet < nextLayerTrackletsNum &&
           tf->getTracklets()[iLayer + 1][iNextLayerTracklet].firstClusterIndex ==
             nextLayerClusterIndex;
           ++iNextLayerTracklet) {

        const Tracklet& nextTracklet{tf->getTracklets()[iLayer + 1][iNextLayerTracklet]};
        const float deltaTanLambda{std::abs(currentTracklet.tanLambda - nextTracklet.tanLambda)};
        const float deltaPhi{std::abs(currentTracklet.phi - nextTracklet.phi)};

        if (deltaTanLambda < mTrkParams.CellMaxDeltaTanLambda &&
            (deltaPhi < mTrkParams.CellMaxDeltaPhi ||
             std::abs(deltaPhi - constants::math::TwoPi) < mTrkParams.CellMaxDeltaPhi)) {

          const float averageTanLambda{0.5f * (currentTracklet.tanLambda + nextTracklet.tanLambda)};
          const float directionZIntersection{-averageTanLambda * firstCellCluster.radius +
                                             firstCellCluster.zCoordinate};

          unsigned short romin = std::min(std::min(currentTracklet.rof[0], currentTracklet.rof[1]), nextTracklet.rof[1]);
          unsigned short romax = std::max(std::max(currentTracklet.rof[0], currentTracklet.rof[1]), nextTracklet.rof[1]);
          bool deltaZflag{false};
          gsl::span<const float3> primaryVertices{tf->getPrimaryVertices(romin, romax)};
          for (const auto& primaryVertex : primaryVertices)
            deltaZflag = std::abs(directionZIntersection - primaryVertex.z) < mTrkParams.CellMaxDeltaZ[iLayer];

          if (deltaZflag) {

            const Cluster& thirdCellCluster{
              tf->getClusters()[iLayer + 2][nextTracklet.secondClusterIndex]};

            const float thirdCellClusterQuadraticRCoordinate{thirdCellCluster.radius *
                                                             thirdCellCluster.radius};

            const float3 secondDeltaVector{thirdCellCluster.xCoordinate - firstCellCluster.xCoordinate,
                                           thirdCellCluster.yCoordinate - firstCellCluster.yCoordinate,
                                           thirdCellClusterQuadraticRCoordinate -
                                             firstCellClusterQuadraticRCoordinate};

            float3 cellPlaneNormalVector{math_utils::crossProduct(firstDeltaVector, secondDeltaVector)};

            const float vectorNorm{std::hypot(cellPlaneNormalVector.x, cellPlaneNormalVector.y, cellPlaneNormalVector.z)};

            if (vectorNorm < constants::math::FloatMinThreshold ||
                std::abs(cellPlaneNormalVector.z) < constants::math::FloatMinThreshold) {
              continue;
            }

            const float inverseVectorNorm{1.0f / vectorNorm};
            const float3 normVect{cellPlaneNormalVector.x * inverseVectorNorm,
                                  cellPlaneNormalVector.y * inverseVectorNorm,
                                  cellPlaneNormalVector.z * inverseVectorNorm};
            const float planeDistance{-normVect.x * (secondCellCluster.xCoordinate - tf->getBeamX()) -
                                      (normVect.y * secondCellCluster.yCoordinate - tf->getBeamY()) -
                                      normVect.z * secondCellClusterQuadraticRCoordinate};
            const float normVectZsquare{normVect.z * normVect.z};
            const float cellRadius{std::sqrt(
              (1.0f - normVectZsquare - 4.0f * planeDistance * normVect.z) /
              (4.0f * normVectZsquare))};
            const float2 circleCenter{-0.5f * normVect.x / normVect.z,
                                      -0.5f * normVect.y / normVect.z};
            const float dca{std::abs(cellRadius - std::hypot(circleCenter.x, circleCenter.y))};

            if (dca > mTrkParams.CellMaxDCA[iLayer]) {
              continue;
            }

            const float cellTrajectoryCurvature{1.0f / cellRadius};
            if (iLayer > 0 &&
                tf->getCellsLookupTable()[iLayer - 1][iTracklet] == constants::its::UnusedIndex) {

              tf->getCellsLookupTable()[iLayer - 1][iTracklet] =
                tf->getCells()[iLayer].size();
            }

            tf->getCells()[iLayer].emplace_back(
              currentTracklet.firstClusterIndex, nextTracklet.firstClusterIndex, nextTracklet.secondClusterIndex,
              iTracklet, iNextLayerTracklet, normVect, cellTrajectoryCurvature);
          }
        }
      }
    }
  }
  std::cout << "Number of cells " << std::endl;
  for (int iL{0}; iL < 5; ++iL) {
    std::cout << tf->getCells()[iL].size() << std::endl;
  }
}

void TrackerTraitsCPU::refitTracks(const std::vector<std::vector<TrackingFrameInfo>>& tf, std::vector<TrackITSExt>& tracks)
{
  std::vector<const Cell*> cells;
  for (int iLayer = 0; iLayer < mTrkParams.CellsPerRoad(); iLayer++) {
    cells.push_back(mTimeFrame->getCells()[iLayer].data());
  }
  std::vector<const Cluster*> clusters;
  for (int iLayer = 0; iLayer < mTrkParams.NLayers; iLayer++) {
    clusters.push_back(mTimeFrame->getClusters()[iLayer].data());
  }
  mChainRunITSTrackFit(*mChain, mTimeFrame->getRoads(), clusters, cells, tf, tracks);
}

} // namespace its
} // namespace o2

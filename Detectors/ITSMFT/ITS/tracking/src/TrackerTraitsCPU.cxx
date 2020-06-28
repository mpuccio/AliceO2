// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
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
    for (int iLayer{0}; iLayer < constants::its::TrackletsPerRoad; ++iLayer) {
      gsl::span<const Cluster> layer0 = tf->getClustersOnLayer(rof0, iLayer);
      if (layer0.empty()) {
        return;
      }

      gsl::span<const float3> primaryVertices = tf->getPrimaryVertices(rof0);
      const int currentLayerClustersNum{static_cast<int>(tf->getClustersOnLayer(rof0, iLayer).size())};

      for (int iCluster{0}; iCluster < currentLayerClustersNum; ++iCluster) {
        const Cluster& currentCluster{layer0[iCluster]};

        if (tf->isClusterUsed(iLayer, currentCluster.clusterId)) {
          continue;
        }

        for (auto& primaryVertex : primaryVertices) {
          const float tanLambda{(currentCluster.zCoordinate - primaryVertex.z) / currentCluster.radius};
          const float directionZIntersection{tanLambda * (constants::its::LayersRCoordinate()[iLayer + 1] -
                                                          currentCluster.radius) +
                                             currentCluster.zCoordinate};

          const int4 selectedBinsRect{getBinsRect(currentCluster, iLayer, directionZIntersection,
                                                  mTrkParams.TrackletMaxDeltaZ[iLayer], mTrkParams.TrackletMaxDeltaPhi)};

          if (selectedBinsRect.x == 0 && selectedBinsRect.y == 0 && selectedBinsRect.z == 0 && selectedBinsRect.w == 0) {
            continue;
          }

          int phiBinsNum{selectedBinsRect.w - selectedBinsRect.y + 1};

          if (phiBinsNum < 0) {
            phiBinsNum += constants::index_table::PhiBins;
          }

          int minRof = (rof0 > 0) ? rof0 - 1 : 0;
          for (int rof1{minRof}; rof1 < rof0 + 1; ++rof1) {
            gsl::span<const Cluster> layer1 = tf->getClustersOnLayer(rof1, iLayer + 1);
            if (layer1.empty()) {
              continue;
            }

            for (int iPhiBin{selectedBinsRect.y}, iPhiCount{0}; iPhiCount < phiBinsNum;
                 iPhiBin = ++iPhiBin == constants::index_table::PhiBins ? 0 : iPhiBin, iPhiCount++) {
              const int firstBinIndex{index_table_utils::getBinIndex(selectedBinsRect.x, iPhiBin)};
              const int maxBinIndex{firstBinIndex + selectedBinsRect.z - selectedBinsRect.x + 1};
              const int firstRowClusterIndex = tf->getIndexTables(rof1)[iLayer][firstBinIndex];
              const int maxRowClusterIndex = tf->getIndexTables(rof1)[iLayer][maxBinIndex];

              for (int iNextCluster{firstRowClusterIndex}; iNextCluster < maxRowClusterIndex; ++iNextCluster) {

                const Cluster& nextCluster{layer1[iNextCluster]};

                if (tf->isClusterUsed(iLayer + 1, nextCluster.clusterId)) {
                  continue;
                }

                const float deltaZ{gpu::GPUCommonMath::Abs(tanLambda * (nextCluster.radius - currentCluster.radius) +
                                                           currentCluster.zCoordinate - nextCluster.zCoordinate)};
                const float deltaPhi{gpu::GPUCommonMath::Abs(currentCluster.phi - nextCluster.phi)};

                if (deltaZ < mTrkParams.TrackletMaxDeltaZ[iLayer] &&
                    (deltaPhi < mTrkParams.TrackletMaxDeltaPhi ||
                     gpu::GPUCommonMath::Abs(deltaPhi - constants::math::TwoPi) < mTrkParams.TrackletMaxDeltaPhi)) {

                  if (iLayer > 0 &&
                      tf->getTrackletsLookupTable()[iLayer - 1][iCluster] == constants::its::UnusedIndex) {

                    tf->getTrackletsLookupTable()[iLayer - 1][iCluster] =
                      tf->getTracklets()[iLayer].size();
                  }

                  tf->getTracklets()[iLayer].emplace_back(tf->getSortedIndex(rof0, iLayer, iCluster), tf->getSortedIndex(rof1, iLayer + 1, iNextCluster), currentCluster,
                                                          nextCluster, rof0, rof1);
                }
              }
            }
          }
        }
      }
    }
  }
}

void TrackerTraitsCPU::computeLayerCells()
{
  TimeFrame* tf = mTimeFrame;
  for (int iLayer{0}; iLayer < constants::its::CellsPerRoad; ++iLayer) {

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
}

void TrackerTraitsCPU::refitTracks(const std::array<std::vector<TrackingFrameInfo>, 7>& tf, std::vector<TrackITSExt>& tracks)
{
  std::array<const Cell*, 5> cells;
  for (int iLayer = 0; iLayer < 5; iLayer++) {
    cells[iLayer] = mTimeFrame->getCells()[iLayer].data();
  }
  std::array<const Cluster*, 7> clusters;
  for (int iLayer = 0; iLayer < 7; iLayer++) {
    clusters[iLayer] = mTimeFrame->getClusters()[iLayer].data();
  }
  mChainRunITSTrackFit(*mChain, mTimeFrame->getRoads(), clusters, cells, tf, tracks);
}

} // namespace its
} // namespace o2

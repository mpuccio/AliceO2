// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Tracker.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "ITSReconstruction/CA/Tracker.h"

#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "DetectorsBase/Constants.h"
#include "ITSReconstruction/CA/Cell.h"
#include "ITSReconstruction/CA/Constants.h"
#include "ITSReconstruction/CA/Event.h"
#include "ITSReconstruction/CA/Layer.h"
#include "ITSReconstruction/CA/MathUtils.h"
#include "ITSReconstruction/CA/PrimaryVertexContext.h"
#include "ITSReconstruction/CA/Tracklet.h"

namespace o2
{
namespace ITS
{
namespace CA
{

/// Clusters are given from outside inward (cluster1 is the outermost). The innermost cluster is given in the tracking frame coordinates
/// whereas the others are referred to the global frame. This function is almost a clone of CookSeed, adapted to return a TrackParCov
Base::Track::TrackParCov Tracker::buildTrackSeed(const Cluster& cluster1, const Cluster& cluster2, const Cluster& cluster3, const TrackingFrameInfo& tf3) {
  const float ca = std::cos(tf3.alphaTrackingFrame), sa = std::sin(tf3.alphaTrackingFrame);
  const float x1 =  cluster1.xCoordinate * ca + cluster1.yCoordinate * sa;
  const float y1 = -cluster1.xCoordinate * sa + cluster1.yCoordinate * ca;
  const float z1 =  cluster1.zCoordinate;
  const float x2 =  cluster2.xCoordinate * ca + cluster2.yCoordinate * sa;
  const float y2 = -cluster2.xCoordinate * sa + cluster2.yCoordinate * ca;
  const float z2 =  cluster2.zCoordinate;
  const float x3 =  tf3.positionTrackingFrame[0];
  const float y3 =  tf3.positionTrackingFrame[1];
  const float z3 =  cluster3.zCoordinate;

  const float crv = MathUtils::computeCurvature(x1, y1, x2, y2, x3, y3);
  const float x0  = MathUtils::computeCurvatureCentreX(x1, y1, x2, y2, x3, y3);
  const float tgl12 = MathUtils::computeTanDipAngle(x1, y1, x2, y2, z1, z2);
  const float tgl23 = MathUtils::computeTanDipAngle(x2, y2, x3, y3, z2, z3);

  const float fy = 1. / (cluster2.rCoordinate - cluster3.rCoordinate);
  const float& tz = fy;
  const float cy = (MathUtils::computeCurvature(x1, y1, x2, y2 + Constants::ITS::Resolution, x3, y3) - crv) / \
    (Constants::ITS::Resolution * mBz * Base::Constants::kB2C) * 20.f; // FIXME: MS contribution to the cov[14] (*20 added)
  constexpr float s2 = Constants::ITS::Resolution * Constants::ITS::Resolution;

  return Base::Track::TrackParCov(
    tf3.xTrackingFrame,
    tf3.alphaTrackingFrame,
    {y3, z3, crv * (x3 - x0), 0.5f * (tgl12 + tgl23), std::abs(mBz) < Base::Constants::kAlmost0 ? Base::Constants::kAlmost0 : crv / (mBz * Base::Constants::kB2C)},
    {
      s2,
      0.f,     s2,
      s2 * fy, 0.f,     s2 * fy * fy,
      0.f,     s2 * tz, 0.f,          s2 * tz * tz,
      s2 * cy, 0.f,     s2 * fy * cy, 0.f,          s2 * cy * cy
    }
  );
}

Tracker::Tracker(const Event& event)
  : mBz{0.f}, mEvent{event}
{
  // Nothing to do
}

std::vector<std::vector<Road>> Tracker::clustersToTracks()
{
  const int verticesNum = mEvent.getPrimaryVerticesNum();
  std::vector<std::vector<Road>> roads;
  roads.reserve(verticesNum);

  for (int iVertex = 0; iVertex < verticesNum; ++iVertex) {
    PrimaryVertexContext primaryVertexContext{ mEvent, iVertex };

    computeTracklets(primaryVertexContext);
    computeCells(primaryVertexContext);
    findCellsNeighbours(primaryVertexContext);
    findTracks(primaryVertexContext);

    roads.emplace_back(primaryVertexContext.roads);
  }

  return roads;
}

void Tracker::computeTracklets(PrimaryVertexContext& primaryVertexContext)
{
  for (int iLayer = 0; iLayer < Constants::ITS::TrackletsPerRoad; ++iLayer) {
    const Layer& currentLayer = mEvent.getLayer(iLayer);
    const Layer& nextLayer = mEvent.getLayer(iLayer + 1);

    if (currentLayer.getClusters().empty() || nextLayer.getClusters().empty()) {
      continue;
    }

    const int currentLayerClustersNum = currentLayer.getClustersSize();

    if (iLayer < Constants::ITS::CellsPerRoad) {
      primaryVertexContext.trackletsLookupTable[iLayer].resize(nextLayer.getClustersSize(),
                                                               Constants::ITS::UnusedIndex);
    }

    for (int iCluster = 0; iCluster < currentLayerClustersNum; ++iCluster) {
      const Cluster& currentCluster = primaryVertexContext.clusters[iLayer][iCluster];

      if (primaryVertexContext.usedClusters[iLayer][currentCluster.clusterId]) {
        continue;
      }

      const float tanLambda =
        (currentCluster.zCoordinate - mEvent.getPrimaryVertex(primaryVertexContext.primaryVertexIndex)[2]) /
        currentCluster.rCoordinate;
      const float directionZIntersection =
        tanLambda * (Constants::ITS::LayersRCoordinate[iLayer + 1] - currentCluster.rCoordinate) +
        currentCluster.zCoordinate;

      const std::vector<int> nextLayerBinsSubset = primaryVertexContext.indexTables[iLayer].selectBins(
        directionZIntersection - 2 * Constants::Thresholds::ZCoordinateCut,
        directionZIntersection + 2 * Constants::Thresholds::ZCoordinateCut,
        currentCluster.phiCoordinate - Constants::Thresholds::PhiCoordinateCut,
        currentCluster.phiCoordinate + Constants::Thresholds::PhiCoordinateCut);

      if (nextLayerBinsSubset.empty()) {
        continue;
      }

      bool isFirstTrackletFromCurrentCluster = true;

      const int lastClusterBin = nextLayerBinsSubset.size() - 1;

      for (int iClusterBin = 0; iClusterBin <= lastClusterBin; ++iClusterBin) {
        const int currentBinIndex = nextLayerBinsSubset[iClusterBin];
        const int currentBinFirstCluster = primaryVertexContext.indexTables[iLayer].getBin(currentBinIndex);
        const int nextBinFirstClusterIndex = primaryVertexContext.indexTables[iLayer].getBin(currentBinIndex + 1);

        for (int iNextLayerCluster = currentBinFirstCluster; iNextLayerCluster < nextBinFirstClusterIndex;
             ++iNextLayerCluster) {
          const Cluster& nextCluster = primaryVertexContext.clusters[iLayer + 1][iNextLayerCluster];

          if (primaryVertexContext.usedClusters[iLayer + 1][currentCluster.clusterId]) {
            continue;
          }

          const float deltaZ = std::abs(tanLambda * (nextCluster.rCoordinate - currentCluster.rCoordinate) +
                                        currentCluster.zCoordinate - nextCluster.zCoordinate);
          const float deltaPhi = std::abs(currentCluster.phiCoordinate - nextCluster.phiCoordinate);

          if (deltaZ < Constants::Thresholds::TrackletMaxDeltaZ[iLayer] &&
              (deltaPhi < Constants::Thresholds::PhiCoordinateCut ||
               std::abs(deltaPhi - Base::Constants::k2PI) < Constants::Thresholds::PhiCoordinateCut)) {
            if (iLayer > 0 && isFirstTrackletFromCurrentCluster) {
              primaryVertexContext.trackletsLookupTable[iLayer - 1][iCluster] =
                primaryVertexContext.tracklets[iLayer].size();
              isFirstTrackletFromCurrentCluster = false;
            }

            const float trackletTanLambda = (currentCluster.zCoordinate - nextCluster.zCoordinate) /
                                            (currentCluster.rCoordinate - nextCluster.rCoordinate);
            const float trackletPhi = std::atan2(currentCluster.yCoordinate - nextCluster.yCoordinate,
                                                 currentCluster.xCoordinate - nextCluster.xCoordinate);

            primaryVertexContext.tracklets[iLayer].emplace_back(iCluster, iNextLayerCluster, trackletTanLambda,
                                                                trackletPhi);
          }
        }
      }
    }
  }
}

void Tracker::computeCells(PrimaryVertexContext& primaryVertexContext)
{
  const std::array<float, 3>& primaryVertex = mEvent.getPrimaryVertex(primaryVertexContext.primaryVertexIndex);

  for (int iLayer = 0; iLayer < Constants::ITS::CellsPerRoad; ++iLayer) {
    if (primaryVertexContext.tracklets[iLayer + 1].empty() || primaryVertexContext.tracklets[iLayer].empty()) {
      continue;
    }

    if (iLayer < Constants::ITS::CellsPerRoad - 1) {
      primaryVertexContext.cellsLookupTable[iLayer].resize(primaryVertexContext.tracklets[iLayer + 1].size(),
                                                           Constants::ITS::UnusedIndex);
    }

    const int currentLayerTrackletsNum = primaryVertexContext.tracklets[iLayer].size();

    for (int iTracklet = 0; iTracklet < currentLayerTrackletsNum; ++iTracklet) {
      const Tracklet& currentTracklet = primaryVertexContext.tracklets[iLayer][iTracklet];
      const int nextLayerClusterIndex = currentTracklet.secondClusterIndex;
      const int nextLayerFirstTrackletIndex = primaryVertexContext.trackletsLookupTable[iLayer][nextLayerClusterIndex];

      if (nextLayerFirstTrackletIndex == Constants::ITS::UnusedIndex) {
        continue;
      }

      const Cluster& firstCellCluster = primaryVertexContext.clusters[iLayer][currentTracklet.firstClusterIndex];
      const Cluster& secondCellCluster =
        primaryVertexContext.clusters[iLayer + 1][currentTracklet.secondClusterIndex];

      const float firstCellClusterQuadraticRCoordinate = firstCellCluster.rCoordinate * firstCellCluster.rCoordinate;
      const float secondCellClusterQuadraticRCoordinate = secondCellCluster.rCoordinate * secondCellCluster.rCoordinate;

      const std::array<float, 3> firstDeltaVector{ { secondCellCluster.xCoordinate - firstCellCluster.xCoordinate,
                                                     secondCellCluster.yCoordinate - firstCellCluster.yCoordinate,
                                                     secondCellClusterQuadraticRCoordinate -
                                                       firstCellClusterQuadraticRCoordinate } };

      bool isFirstCellForCurrentTracklet = true;
      const int nextLayerTrackletsNum = primaryVertexContext.tracklets[iLayer + 1].size();

      for (int iNextLayerTracklet = nextLayerFirstTrackletIndex;
           iNextLayerTracklet < nextLayerTrackletsNum &&
           primaryVertexContext.tracklets[iLayer + 1][iNextLayerTracklet].firstClusterIndex == nextLayerClusterIndex;
           ++iNextLayerTracklet) {
        const Tracklet& nextTracklet = primaryVertexContext.tracklets[iLayer + 1][iNextLayerTracklet];
        const float deltaTanLambda = std::abs(currentTracklet.tanLambda - nextTracklet.tanLambda);
        const float deltaPhi = std::abs(currentTracklet.phiCoordinate - nextTracklet.phiCoordinate);

        if (deltaTanLambda < Constants::Thresholds::CellMaxDeltaTanLambda &&
            (deltaPhi < Constants::Thresholds::CellMaxDeltaPhiThreshold ||
             std::abs(deltaPhi - Base::Constants::k2PI) < Constants::Thresholds::CellMaxDeltaPhiThreshold)) {
          const float averageTanLambda = 0.5f * (currentTracklet.tanLambda + nextTracklet.tanLambda);
          const float directionZIntersection =
            -averageTanLambda * firstCellCluster.rCoordinate + firstCellCluster.zCoordinate;
          const float deltaZ = std::abs(directionZIntersection - primaryVertex[2]);

          if (deltaZ < Constants::Thresholds::CellMaxDeltaZ[iLayer]) {
            const Cluster& thirdCellCluster =
              primaryVertexContext.clusters[iLayer + 2][nextTracklet.secondClusterIndex];

            const float thirdCellClusterQuadraticRCoordinate =
              thirdCellCluster.rCoordinate * thirdCellCluster.rCoordinate;

            const std::array<float, 3> secondDeltaVector{ { thirdCellCluster.xCoordinate - firstCellCluster.xCoordinate,
                                                            thirdCellCluster.yCoordinate - firstCellCluster.yCoordinate,
                                                            thirdCellClusterQuadraticRCoordinate -
                                                              firstCellClusterQuadraticRCoordinate } };

            std::array<float, 3> cellPlaneNormalVector{ MathUtils::crossProduct(firstDeltaVector,
                                                                                  secondDeltaVector) };

            const float vectorNorm = std::sqrt(cellPlaneNormalVector[0] * cellPlaneNormalVector[0] +
                                               cellPlaneNormalVector[1] * cellPlaneNormalVector[1] +
                                               cellPlaneNormalVector[2] * cellPlaneNormalVector[2]);

            if (vectorNorm < Base::Constants::kAlmost0 ||
                std::abs(cellPlaneNormalVector[2]) < Base::Constants::kAlmost0) {
              continue;
            }

            const float inverseVectorNorm = 1.0f / vectorNorm;

            const std::array<float, 3> normalizedPlaneVector{ { cellPlaneNormalVector[0] * inverseVectorNorm,
                                                                cellPlaneNormalVector[1] * inverseVectorNorm,
                                                                cellPlaneNormalVector[2] * inverseVectorNorm } };

            const float planeDistance = -normalizedPlaneVector[0] * (secondCellCluster.xCoordinate - primaryVertex[0]) -
                                        (normalizedPlaneVector[1] * secondCellCluster.yCoordinate - primaryVertex[1]) -
                                        normalizedPlaneVector[2] * secondCellClusterQuadraticRCoordinate;

            const float normalizedPlaneVectorQuadraticZCoordinate = normalizedPlaneVector[2] * normalizedPlaneVector[2];

            const float cellTrajectoryRadius = std::sqrt(
              (1.0f - normalizedPlaneVectorQuadraticZCoordinate - 4.0f * planeDistance * normalizedPlaneVector[2]) /
              (4.0f * normalizedPlaneVectorQuadraticZCoordinate));

            const std::array<float, 2> circleCenter{ { -0.5f * normalizedPlaneVector[0] / normalizedPlaneVector[2],
                                                       -0.5f * normalizedPlaneVector[1] / normalizedPlaneVector[2] } };

            const float distanceOfClosestApproach = std::abs(
              cellTrajectoryRadius - std::sqrt(circleCenter[0] * circleCenter[0] + circleCenter[1] * circleCenter[1]));

            if (distanceOfClosestApproach >
                Constants::Thresholds::CellMaxDCA[iLayer]) {
              continue;
            }

            const float cellTrajectoryCurvature = 1.0f / cellTrajectoryRadius;

            if (isFirstCellForCurrentTracklet && iLayer > 0) {
              primaryVertexContext.cellsLookupTable[iLayer - 1][iTracklet] = primaryVertexContext.cells[iLayer].size();
              isFirstCellForCurrentTracklet = false;
            }

            primaryVertexContext.cells[iLayer].emplace_back(
              currentTracklet.firstClusterIndex, nextTracklet.firstClusterIndex, nextTracklet.secondClusterIndex,
              iTracklet, iNextLayerTracklet, normalizedPlaneVector, cellTrajectoryCurvature);
          }
        }
      }
    }
  }
}

void Tracker::findCellsNeighbours(PrimaryVertexContext& primaryVertexContext)
{
  for (int iLayer = 0; iLayer < Constants::ITS::CellsPerRoad - 1; ++iLayer) {
    if (primaryVertexContext.cells[iLayer + 1].empty() || primaryVertexContext.cellsLookupTable[iLayer].empty()) {
      continue;
    }

    int layerCellsNum = primaryVertexContext.cells[iLayer].size();

    for (int iCell = 0; iCell < layerCellsNum; ++iCell) {
      const Cell& currentCell = primaryVertexContext.cells[iLayer][iCell];

      const int nextLayerTrackletIndex = currentCell.getTrackletIndex(1);
      const int nextLayerFirstCellIndex = primaryVertexContext.cellsLookupTable[iLayer][nextLayerTrackletIndex];

      if (nextLayerFirstCellIndex == Constants::ITS::UnusedIndex) {
        continue;
      }

      const int nextLayerCellsNum = primaryVertexContext.cells[iLayer + 1].size();

      for (int iNextLayerCell = nextLayerFirstCellIndex;
           iNextLayerCell < nextLayerCellsNum &&
           primaryVertexContext.cells[iLayer + 1][iNextLayerCell].getTrackletIndex(0) == nextLayerTrackletIndex;
           ++iNextLayerCell) {
        Cell& nextCell = primaryVertexContext.cells[iLayer + 1][iNextLayerCell];
        const std::array<float, 3> currentCellNormalVector = currentCell.getNormalVectorCoordinates();
        const std::array<float, 3> nextCellNormalVector = nextCell.getNormalVectorCoordinates();

        const std::array<float, 3> normalVectorsDeltaVector = {
          { currentCellNormalVector[0] - nextCellNormalVector[0], currentCellNormalVector[1] - nextCellNormalVector[1],
            currentCellNormalVector[2] - nextCellNormalVector[2] }
        };

        const float deltaNormalVectorsModulus = (normalVectorsDeltaVector[0] * normalVectorsDeltaVector[0]) +
                                                (normalVectorsDeltaVector[1] * normalVectorsDeltaVector[1]) +
                                                (normalVectorsDeltaVector[2] * normalVectorsDeltaVector[2]);

        const float deltaCurvature = std::abs(currentCell.getCurvature() - nextCell.getCurvature());

        if (deltaNormalVectorsModulus < Constants::Thresholds::NeighbourCellMaxNormalVectorsDelta[iLayer] &&
            deltaCurvature < Constants::Thresholds::NeighbourCellMaxCurvaturesDelta[iLayer]) {
          nextCell.combineCells(currentCell, iCell);
        }
      }
    }
  }
}

void Tracker::findRoads(PrimaryVertexContext& primaryVertexContext)
{
  for (int iLevel = Constants::ITS::CellsPerRoad; iLevel >= Constants::Thresholds::CellsMinLevel; --iLevel) {
    const int minimumLevel = iLevel - 1;

    for (int iLayer = Constants::ITS::CellsPerRoad - 1; iLayer >= minimumLevel; --iLayer) {
      const int levelCellsNum = primaryVertexContext.cells[iLayer].size();

      for (int iCell = 0; iCell < levelCellsNum; ++iCell) {
        Cell& currentCell = primaryVertexContext.cells[iLayer][iCell];

        if (currentCell.getLevel() != iLevel) {
          continue;
        }

        primaryVertexContext.roads.emplace_back(iLayer, iCell);

        const int cellNeighboursNum = currentCell.getNumberOfNeighbours();
        int isFirstValidNeighbour = true;

        for (int iNeighbourCell = 0; iNeighbourCell < cellNeighboursNum; ++iNeighbourCell) {
          const int neighbourCellId = currentCell.getNeighbourCellId(iNeighbourCell);
          const Cell& neighbourCell = primaryVertexContext.cells[iLayer - 1][neighbourCellId];

          if (iLevel - 1 != neighbourCell.getLevel()) {
            continue;
          }

          if (isFirstValidNeighbour) {
            isFirstValidNeighbour = false;
          } else {
            primaryVertexContext.roads.emplace_back(iLayer, iCell);
          }

          traverseCellsTree(primaryVertexContext, neighbourCellId, iLayer - 1);
        }

        // TODO: crosscheck for short track iterations
        // currentCell.setLevel(0);
      }
    }
  }
}

void Tracker::findTracks(PrimaryVertexContext& vertexContext) {
  vertexContext.tracks.reserve(vertexContext.roads.size());
  std::vector<Track> tracks;
  tracks.reserve(vertexContext.roads.size());
  for (auto& road : vertexContext.roads) {
    std::array<int, 7> clusters {Constants::ITS::UnusedIndex};
    int lastCellLevel = -1;
    for (int iCell{0}; iCell < Constants::ITS::CellsPerRoad; ++iCell) {
      const int cellIndex = road[iCell];
      if (cellIndex == Constants::ITS::UnusedIndex) {
        continue;
      } else {
        clusters[iCell] = vertexContext.cells[iCell][cellIndex].getClusterIndex(0);
        clusters[iCell + 1] = vertexContext.cells[iCell][cellIndex].getClusterIndex(1);
        clusters[iCell + 2] = vertexContext.cells[iCell][cellIndex].getClusterIndex(2);
        lastCellLevel = iCell;
      }
    }

    /// Track seed preparation. Clusters are numbered progressively from the outermost to the innermost.
    const auto& cluster1_glo = mEvent.getLayer(lastCellLevel + 2).getCluster(clusters[lastCellLevel + 2]);
    const auto& cluster2_glo = mEvent.getLayer(lastCellLevel + 1).getCluster(clusters[lastCellLevel + 1]);
    const auto& cluster3_glo = mEvent.getLayer(lastCellLevel).getCluster(clusters[lastCellLevel]);
    /// From global cluster id to tracking cluster id (== the one used as input of the tracking code)
    for (int iC{0}; iC < clusters.size(); iC++) {
      clusters[iC] = mEvent.getLayer(iC).getCluster(clusters[iC]).clusterId;
    }
    const auto& cluster3_tf = mEvent.getLayer(lastCellLevel).getTrackingFrameInfo(clusters[lastCellLevel]);
    
    Track temporaryTrack{
      buildTrackSeed(cluster1_glo, cluster2_glo, cluster3_glo, cluster3_tf),
      0.f,
      clusters
    };

    bool fitSuccess = true;
    for (int iCluster{Constants::ITS::LayersNumber-3}; iCluster--; ) {
      if (temporaryTrack.mClusters[iCluster] == Constants::ITS::UnusedIndex) {
        continue;
      }

      const TrackingFrameInfo& trackingHit = mEvent.getTrackingFrameInfo(iCluster, temporaryTrack.mClusters[iCluster]);

      fitSuccess = temporaryTrack.mParam.Rotate(trackingHit.alphaTrackingFrame);
      if (!fitSuccess) {
        break;
      }

      fitSuccess = temporaryTrack.mParam.PropagateTo(trackingHit.xTrackingFrame, mBz);
      if (!fitSuccess) {
        break;
      }

      temporaryTrack.mChi2 += temporaryTrack.mParam.GetPredictedChi2(trackingHit.positionTrackingFrame,trackingHit.covarianceTrackingFrame);
      fitSuccess = temporaryTrack.mParam.Update(trackingHit.positionTrackingFrame,trackingHit.covarianceTrackingFrame);
      if (!fitSuccess) {
        break;
      }

      const float xx0 = (iCluster > 2) ? 0.008f : 0.003f;            // Rough layer thickness
      constexpr float radiationLength = 9.36f; // Radiation length of Si [cm]
      constexpr float density = 2.33f;         // Density of Si [g/cm^3]
      fitSuccess = temporaryTrack.mParam.CorrectForMaterial(xx0, xx0 * radiationLength * density, true);
      if (!fitSuccess) {
        break;
      }
    }

    if (!fitSuccess) {
      continue;
    }
    tracks.emplace_back(temporaryTrack);
  }

  std::sort(tracks.begin(),tracks.end(),[](Track& track1, Track& track2) {
    return track1.mChi2 > track2.mChi2;
  });

  for (auto& track : tracks) {
    bool sharingCluster = false;
    for (int iCluster{0}; iCluster < Constants::ITS::LayersNumber; ++iCluster) {
      if (track.mClusters[iCluster] == Constants::ITS::UnusedIndex) {
        continue;
      }
      sharingCluster |= vertexContext.usedClusters[iCluster][track.mClusters[iCluster]];
    }

    if (sharingCluster) {
      continue;
    }
    for (int iCluster{0}; iCluster < Constants::ITS::LayersNumber; ++iCluster) {
      if (track.mClusters[iCluster] == Constants::ITS::UnusedIndex) {
        continue;
      }
      vertexContext.usedClusters[iCluster][track.mClusters[iCluster]] = true;
    }
    vertexContext.tracks.emplace_back(track);
  }
}

void Tracker::traverseCellsTree(PrimaryVertexContext& primaryVertexContext, const int currentCellId,
                                  const int currentLayerId)
{
  if (currentLayerId < 0) {
    return;
  }

  Cell& currentCell = primaryVertexContext.cells[currentLayerId][currentCellId];
  const int currentCellLevel = currentCell.getLevel();
  const int cellNeighboursNum = currentCell.getNumberOfNeighbours();
  int isFirstValidNeighbour = true;

  primaryVertexContext.roads.back().addCell(currentLayerId, currentCellId);

  for (int iNeighbourCell = 0; iNeighbourCell < cellNeighboursNum; ++iNeighbourCell) {
    const int neighbourCellId = currentCell.getNeighbourCellId(iNeighbourCell);
    const Cell& neighbourCell = primaryVertexContext.cells[currentLayerId - 1][neighbourCellId];

    if (currentCellLevel - 1 != neighbourCell.getLevel()) {
      continue;
    }

    if (isFirstValidNeighbour) {
      isFirstValidNeighbour = false;

    } else {
      primaryVertexContext.roads.push_back(primaryVertexContext.roads.back());
    }

    traverseCellsTree(primaryVertexContext, neighbourCellId, currentLayerId - 1);
  }

  currentCell.setLevel(0);
}


}
}
}

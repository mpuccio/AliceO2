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
/// \file Tracker.cxx
/// \brief 
///

#include "ITSReconstruction/CA/Tracker.h"

#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

#include "ITSReconstruction/CA/Cell.h"
#include "ITSReconstruction/CA/Constants.h"
#include "ITSReconstruction/CA/Definitions.h"
#include "ITSReconstruction/CA/Event.h"
#include "ITSReconstruction/CA/IndexTableUtils.h"
#include "ITSReconstruction/CA/Layer.h"
#include "ITSReconstruction/CA/MathUtils.h"
#include "ITSReconstruction/CA/PrimaryVertexContext.h"
#include "ITSReconstruction/CA/Tracklet.h"
#include "ITSReconstruction/CA/TrackingUtils.h"

namespace o2
{
namespace ITS
{
namespace CA
{

#if !TRACKINGITSU_GPU_MODE
template<>
void TrackerTraits<false>::computeLayerTracklets(PrimaryVertexContext& primaryVertexContext)
{
  for (int iLayer { 0 }; iLayer < Constants::ITS::TrackletsPerRoad; ++iLayer) {

    if (primaryVertexContext.getClusters()[iLayer].empty() || primaryVertexContext.getClusters()[iLayer + 1].empty()) {

      return;
    }

    const float3 &primaryVertex = primaryVertexContext.getPrimaryVertex();
    const int currentLayerClustersNum { static_cast<int>(primaryVertexContext.getClusters()[iLayer].size()) };

    for (int iCluster { 0 }; iCluster < currentLayerClustersNum; ++iCluster) {

      const Cluster& currentCluster { primaryVertexContext.getClusters()[iLayer][iCluster] };

      const float tanLambda { (currentCluster.zCoordinate - primaryVertex.z) / currentCluster.rCoordinate };
      const float directionZIntersection { tanLambda
          * (Constants::ITS::LayersRCoordinate()[iLayer + 1] - currentCluster.rCoordinate)
          + currentCluster.zCoordinate };

      const int4 selectedBinsRect { TrackingUtils::getBinsRect(currentCluster, iLayer, directionZIntersection) };

      if (selectedBinsRect.x == 0 && selectedBinsRect.y == 0 && selectedBinsRect.z == 0 && selectedBinsRect.w == 0) {

        continue;
      }

      int phiBinsNum { selectedBinsRect.w - selectedBinsRect.y + 1 };

      if (phiBinsNum < 0) {

        phiBinsNum += Constants::IndexTable::PhiBins;
      }

      for (int iPhiBin { selectedBinsRect.y }, iPhiCount { 0 }; iPhiCount < phiBinsNum;
          iPhiBin = ++iPhiBin == Constants::IndexTable::PhiBins ? 0 : iPhiBin, iPhiCount++) {

        const int firstBinIndex { IndexTableUtils::getBinIndex(selectedBinsRect.x, iPhiBin) };
        const int maxBinIndex { firstBinIndex + selectedBinsRect.z - selectedBinsRect.x + 1 };
        const int firstRowClusterIndex = primaryVertexContext.getIndexTables()[iLayer][firstBinIndex];
        const int maxRowClusterIndex = primaryVertexContext.getIndexTables()[iLayer][maxBinIndex];

        for (int iNextLayerCluster { firstRowClusterIndex }; iNextLayerCluster <= maxRowClusterIndex;
            ++iNextLayerCluster) {

          const Cluster& nextCluster { primaryVertexContext.getClusters()[iLayer + 1][iNextLayerCluster] };

          const float deltaZ { MATH_ABS(
              tanLambda * (nextCluster.rCoordinate - currentCluster.rCoordinate) + currentCluster.zCoordinate
                  - nextCluster.zCoordinate) };
          const float deltaPhi { MATH_ABS(currentCluster.phiCoordinate - nextCluster.phiCoordinate) };

          if (deltaZ < Constants::Thresholds::TrackletMaxDeltaZThreshold()[iLayer]
              && (deltaPhi < Constants::Thresholds::PhiCoordinateCut
                  || MATH_ABS(deltaPhi - Constants::Math::TwoPi) < Constants::Thresholds::PhiCoordinateCut)) {

            if (iLayer > 0
                && primaryVertexContext.getTrackletsLookupTable()[iLayer - 1][iCluster]
                    == Constants::ITS::UnusedIndex) {

              primaryVertexContext.getTrackletsLookupTable()[iLayer - 1][iCluster] =
                  primaryVertexContext.getTracklets()[iLayer].size();
            }

            primaryVertexContext.getTracklets()[iLayer].emplace_back(iCluster, iNextLayerCluster, currentCluster,
                nextCluster);
          }
        }
      }
    }
  }
}

template<>
void TrackerTraits<false>::computeLayerCells(PrimaryVertexContext& primaryVertexContext)
{
  for (int iLayer { 0 }; iLayer < Constants::ITS::CellsPerRoad; ++iLayer) {

    if (primaryVertexContext.getTracklets()[iLayer + 1].empty()
        || primaryVertexContext.getTracklets()[iLayer].empty()) {

      return;
    }

    const float3 &primaryVertex = primaryVertexContext.getPrimaryVertex();
    const int currentLayerTrackletsNum { static_cast<int>(primaryVertexContext.getTracklets()[iLayer].size()) };

    for (int iTracklet { 0 }; iTracklet < currentLayerTrackletsNum; ++iTracklet) {

      const Tracklet& currentTracklet { primaryVertexContext.getTracklets()[iLayer][iTracklet] };
      const int nextLayerClusterIndex { currentTracklet.secondClusterIndex };
      const int nextLayerFirstTrackletIndex {
          primaryVertexContext.getTrackletsLookupTable()[iLayer][nextLayerClusterIndex] };

      if (nextLayerFirstTrackletIndex == Constants::ITS::UnusedIndex) {

        continue;
      }

      const Cluster& firstCellCluster { primaryVertexContext.getClusters()[iLayer][currentTracklet.firstClusterIndex] };
      const Cluster& secondCellCluster {
          primaryVertexContext.getClusters()[iLayer + 1][currentTracklet.secondClusterIndex] };
      const float firstCellClusterQuadraticRCoordinate { firstCellCluster.rCoordinate * firstCellCluster.rCoordinate };
      const float secondCellClusterQuadraticRCoordinate { secondCellCluster.rCoordinate * secondCellCluster.rCoordinate };
      const float3 firstDeltaVector { secondCellCluster.xCoordinate - firstCellCluster.xCoordinate,
          secondCellCluster.yCoordinate - firstCellCluster.yCoordinate, secondCellClusterQuadraticRCoordinate
              - firstCellClusterQuadraticRCoordinate };
      const int nextLayerTrackletsNum { static_cast<int>(primaryVertexContext.getTracklets()[iLayer + 1].size()) };

      for (int iNextLayerTracklet { nextLayerFirstTrackletIndex };
          iNextLayerTracklet < nextLayerTrackletsNum
              && primaryVertexContext.getTracklets()[iLayer + 1][iNextLayerTracklet].firstClusterIndex
                  == nextLayerClusterIndex; ++iNextLayerTracklet) {

        const Tracklet& nextTracklet { primaryVertexContext.getTracklets()[iLayer + 1][iNextLayerTracklet] };
        const float deltaTanLambda { std::abs(currentTracklet.tanLambda - nextTracklet.tanLambda) };
        const float deltaPhi { std::abs(currentTracklet.phiCoordinate - nextTracklet.phiCoordinate) };

        if (deltaTanLambda < Constants::Thresholds::CellMaxDeltaTanLambdaThreshold
            && (deltaPhi < Constants::Thresholds::CellMaxDeltaPhiThreshold
                || std::abs(deltaPhi - Constants::Math::TwoPi) < Constants::Thresholds::CellMaxDeltaPhiThreshold)) {

          const float averageTanLambda { 0.5f * (currentTracklet.tanLambda + nextTracklet.tanLambda) };
          const float directionZIntersection { -averageTanLambda * firstCellCluster.rCoordinate
              + firstCellCluster.zCoordinate };
          const float deltaZ { std::abs(directionZIntersection - primaryVertex.z) };

          if (deltaZ < Constants::Thresholds::CellMaxDeltaZThreshold()[iLayer]) {

            const Cluster& thirdCellCluster {
                primaryVertexContext.getClusters()[iLayer + 2][nextTracklet.secondClusterIndex] };

            const float thirdCellClusterQuadraticRCoordinate { thirdCellCluster.rCoordinate
                * thirdCellCluster.rCoordinate };

            const float3 secondDeltaVector { thirdCellCluster.xCoordinate - firstCellCluster.xCoordinate,
                thirdCellCluster.yCoordinate - firstCellCluster.yCoordinate, thirdCellClusterQuadraticRCoordinate
                    - firstCellClusterQuadraticRCoordinate };

            float3 cellPlaneNormalVector { MathUtils::crossProduct(firstDeltaVector, secondDeltaVector) };

            const float vectorNorm { std::sqrt(
                cellPlaneNormalVector.x * cellPlaneNormalVector.x + cellPlaneNormalVector.y * cellPlaneNormalVector.y
                    + cellPlaneNormalVector.z * cellPlaneNormalVector.z) };

            if (vectorNorm < Constants::Math::FloatMinThreshold
                || std::abs(cellPlaneNormalVector.z) < Constants::Math::FloatMinThreshold) {

              continue;
            }

            const float inverseVectorNorm { 1.0f / vectorNorm };
            const float3 normalizedPlaneVector { cellPlaneNormalVector.x * inverseVectorNorm, cellPlaneNormalVector.y
                * inverseVectorNorm, cellPlaneNormalVector.z * inverseVectorNorm };
            const float planeDistance { -normalizedPlaneVector.x * (secondCellCluster.xCoordinate - primaryVertex.x)
                - (normalizedPlaneVector.y * secondCellCluster.yCoordinate - primaryVertex.y)
                - normalizedPlaneVector.z * secondCellClusterQuadraticRCoordinate };
            const float normalizedPlaneVectorQuadraticZCoordinate { normalizedPlaneVector.z * normalizedPlaneVector.z };
            const float cellTrajectoryRadius { std::sqrt(
                (1.0f - normalizedPlaneVectorQuadraticZCoordinate - 4.0f * planeDistance * normalizedPlaneVector.z)
                    / (4.0f * normalizedPlaneVectorQuadraticZCoordinate)) };
            const float2 circleCenter { -0.5f * normalizedPlaneVector.x / normalizedPlaneVector.z, -0.5f
                * normalizedPlaneVector.y / normalizedPlaneVector.z };
            const float distanceOfClosestApproach { std::abs(
                cellTrajectoryRadius - std::sqrt(circleCenter.x * circleCenter.x + circleCenter.y * circleCenter.y)) };

            if (distanceOfClosestApproach
                > Constants::Thresholds::CellMaxDistanceOfClosestApproachThreshold()[iLayer]) {

              continue;
            }

            const float cellTrajectoryCurvature { 1.0f / cellTrajectoryRadius };

            if (iLayer > 0
                && primaryVertexContext.getCellsLookupTable()[iLayer - 1][iTracklet] == Constants::ITS::UnusedIndex) {

              primaryVertexContext.getCellsLookupTable()[iLayer - 1][iTracklet] =
                  primaryVertexContext.getCells()[iLayer].size();
            }

            primaryVertexContext.getCells()[iLayer].emplace_back(currentTracklet.firstClusterIndex,
                nextTracklet.firstClusterIndex, nextTracklet.secondClusterIndex, iTracklet, iNextLayerTracklet,
                normalizedPlaneVector, cellTrajectoryCurvature);
          }
        }
      }
    }
  }
}
#endif

template<bool IsGPU>
Tracker<IsGPU>::Tracker()
{
  // Nothing to do
}

template<bool IsGPU>
std::vector<std::vector<Road>> Tracker<IsGPU>::clustersToTracks(const Event& event)
{
  const int verticesNum { event.getPrimaryVerticesNum() };
  std::vector<std::vector<Road>> roads { };
  roads.reserve(verticesNum);

  for (int iVertex { 0 }; iVertex < verticesNum; ++iVertex) {

    mPrimaryVertexContext.initialize(event, iVertex);

    computeTracklets();
    computeCells();
    findCellsNeighbours();
    findTracks();
    computeMontecarloLabels();

    roads.emplace_back(mPrimaryVertexContext.getRoads());
  }

  return roads;
}

template<bool IsGPU>
std::vector<std::vector<Road>> Tracker<IsGPU>::clustersToTracksVerbose(const Event& event)
{
  const int verticesNum { event.getPrimaryVerticesNum() };
  std::vector<std::vector<Road>> roads { };
  roads.reserve(verticesNum);

  for (int iVertex { 0 }; iVertex < verticesNum; ++iVertex) {

    clock_t t1 { }, t2 { };
    float diff { };

    t1 = clock();

    mPrimaryVertexContext.initialize(event, iVertex);

    t2 = clock();
    diff = ((float) t2 - (float) t1) / (CLOCKS_PER_SEC / 1000);
    std::cout << std::setw(2) << " - Context initialized in: " << diff << "ms" << std::endl;

    evaluateTask(&Tracker<IsGPU>::computeTracklets, "Tracklets Finding");
    evaluateTask(&Tracker<IsGPU>::computeCells, "Cells Finding");
    evaluateTask(&Tracker<IsGPU>::findCellsNeighbours, "Neighbours Finding");
    evaluateTask(&Tracker<IsGPU>::findTracks, "Tracks Finding");
    evaluateTask(&Tracker<IsGPU>::computeMontecarloLabels, "Computing Montecarlo Labels");

    t2 = clock();
    diff = ((float) t2 - (float) t1) / (CLOCKS_PER_SEC / 1000);
    std::cout << std::setw(2) << " - Vertex " << iVertex + 1 << " completed in: " << diff << "ms" << std::endl;

    roads.emplace_back(mPrimaryVertexContext.getRoads());
  }

  return roads;
}

template<bool IsGPU>
std::vector<std::vector<Road>> Tracker<IsGPU>::clustersToTracksMemoryBenchmark(
    const Event& event, std::ofstream & memoryBenchmarkOutputStream)
{
  const int verticesNum { event.getPrimaryVerticesNum() };
  std::vector<std::vector<Road>> roads { };
  roads.reserve(verticesNum);

  for (int iVertex { 0 }; iVertex < verticesNum; ++iVertex) {

    mPrimaryVertexContext.initialize(event, iVertex);

    for (int iLayer { 0 }; iLayer < Constants::ITS::LayersNumber; ++iLayer) {

      memoryBenchmarkOutputStream << mPrimaryVertexContext.getClusters()[iLayer].size() << "\t";
    }

    memoryBenchmarkOutputStream << std::endl;

#if !TRACKINGITSU_GPU_MODE
    for (int iLayer { 0 }; iLayer < Constants::ITS::TrackletsPerRoad; ++iLayer) {

      memoryBenchmarkOutputStream << mPrimaryVertexContext.getTracklets()[iLayer].capacity() << "\t";
    }

    memoryBenchmarkOutputStream << std::endl;

    computeTracklets();

    for (int iLayer { 0 }; iLayer < Constants::ITS::TrackletsPerRoad; ++iLayer) {

      memoryBenchmarkOutputStream << mPrimaryVertexContext.getTracklets()[iLayer].size() << "\t";
    }

    memoryBenchmarkOutputStream << std::endl;
#endif

    for (int iLayer { 0 }; iLayer < Constants::ITS::CellsPerRoad; ++iLayer) {

      memoryBenchmarkOutputStream << mPrimaryVertexContext.getCells()[iLayer].capacity() << "\t";
    }

    memoryBenchmarkOutputStream << std::endl;

    computeCells();

    for (int iLayer { 0 }; iLayer < Constants::ITS::CellsPerRoad; ++iLayer) {

      memoryBenchmarkOutputStream << mPrimaryVertexContext.getCells()[iLayer].size() << "\t";
    }

    memoryBenchmarkOutputStream << std::endl;

    findCellsNeighbours();
    findTracks();
    computeMontecarloLabels();

    roads.emplace_back(mPrimaryVertexContext.getRoads());

    memoryBenchmarkOutputStream << mPrimaryVertexContext.getRoads().size() << std::endl;
  }

  return roads;
}

template<bool IsGPU>
std::vector<std::vector<Road>> Tracker<IsGPU>::clustersToTracksTimeBenchmark(
    const Event& event, std::ofstream& timeBenchmarkOutputStream)
{
  const int verticesNum = event.getPrimaryVerticesNum();
  std::vector<std::vector<Road>> roads;
  roads.reserve(verticesNum);

  for (int iVertex = 0; iVertex < verticesNum; ++iVertex) {

    clock_t t1, t2;
    float diff, total = .0f;

    t1 = clock();

    mPrimaryVertexContext.initialize(event, iVertex);

    t2 = clock();
    diff = ((float) t2 - (float) t1) / (CLOCKS_PER_SEC / 1000);
    total += diff;
    timeBenchmarkOutputStream << diff << "\t";

    total += evaluateTask(&Tracker<IsGPU>::computeTracklets, nullptr, timeBenchmarkOutputStream);
    total += evaluateTask(&Tracker<IsGPU>::computeCells, nullptr, timeBenchmarkOutputStream);
    total += evaluateTask(&Tracker<IsGPU>::findCellsNeighbours, nullptr, timeBenchmarkOutputStream);
    total += evaluateTask(&Tracker<IsGPU>::findTracks, nullptr, timeBenchmarkOutputStream);
    total += evaluateTask(&Tracker<IsGPU>::computeMontecarloLabels, nullptr, timeBenchmarkOutputStream);

    timeBenchmarkOutputStream << total << std::endl;

    roads.emplace_back(mPrimaryVertexContext.getRoads());
  }

  return roads;
}

template<bool IsGPU>
void Tracker<IsGPU>::computeTracklets()
{
  Trait::computeLayerTracklets(mPrimaryVertexContext);
}

template<bool IsGPU>
void Tracker<IsGPU>::computeCells()
{
  Trait::computeLayerCells(mPrimaryVertexContext);
}

template<bool IsGPU>
void Tracker<IsGPU>::findCellsNeighbours()
{
  for (int iLayer { 0 }; iLayer < Constants::ITS::CellsPerRoad - 1; ++iLayer) {

    if (mPrimaryVertexContext.getCells()[iLayer + 1].empty()
        || mPrimaryVertexContext.getCellsLookupTable()[iLayer].empty()) {

      continue;
    }

    int layerCellsNum { static_cast<int>(mPrimaryVertexContext.getCells()[iLayer].size()) };

    for (int iCell { 0 }; iCell < layerCellsNum; ++iCell) {

      const Cell& currentCell { mPrimaryVertexContext.getCells()[iLayer][iCell] };
      const int nextLayerTrackletIndex { currentCell.getSecondTrackletIndex() };
      const int nextLayerFirstCellIndex { mPrimaryVertexContext.getCellsLookupTable()[iLayer][nextLayerTrackletIndex] };

      if (nextLayerFirstCellIndex != Constants::ITS::UnusedIndex
          && mPrimaryVertexContext.getCells()[iLayer + 1][nextLayerFirstCellIndex].getFirstTrackletIndex()
              == nextLayerTrackletIndex) {

        const int nextLayerCellsNum { static_cast<int>(mPrimaryVertexContext.getCells()[iLayer + 1].size()) };
        mPrimaryVertexContext.getCellsNeighbours()[iLayer].resize(nextLayerCellsNum);

        for (int iNextLayerCell { nextLayerFirstCellIndex };
            iNextLayerCell < nextLayerCellsNum
                && mPrimaryVertexContext.getCells()[iLayer + 1][iNextLayerCell].getFirstTrackletIndex()
                    == nextLayerTrackletIndex; ++iNextLayerCell) {

          Cell& nextCell { mPrimaryVertexContext.getCells()[iLayer + 1][iNextLayerCell] };
          const float3 currentCellNormalVector { currentCell.getNormalVectorCoordinates() };
          const float3 nextCellNormalVector { nextCell.getNormalVectorCoordinates() };
          const float3 normalVectorsDeltaVector { currentCellNormalVector.x - nextCellNormalVector.x,
              currentCellNormalVector.y - nextCellNormalVector.y, currentCellNormalVector.z - nextCellNormalVector.z };

          const float deltaNormalVectorsModulus { (normalVectorsDeltaVector.x * normalVectorsDeltaVector.x)
              + (normalVectorsDeltaVector.y * normalVectorsDeltaVector.y)
              + (normalVectorsDeltaVector.z * normalVectorsDeltaVector.z) };
          const float deltaCurvature { std::abs(currentCell.getCurvature() - nextCell.getCurvature()) };

          if (deltaNormalVectorsModulus < Constants::Thresholds::NeighbourCellMaxNormalVectorsDelta[iLayer]
              && deltaCurvature < Constants::Thresholds::NeighbourCellMaxCurvaturesDelta[iLayer]) {

            mPrimaryVertexContext.getCellsNeighbours()[iLayer][iNextLayerCell].push_back(iCell);

            const int currentCellLevel { currentCell.getLevel() };

            if (currentCellLevel >= nextCell.getLevel()) {

              nextCell.setLevel(currentCellLevel + 1);
            }
          }
        }
      }
    }
  }
}

template<bool IsGPU>
void Tracker<IsGPU>::findTracks()
{
  for (int iLevel { Constants::ITS::CellsPerRoad }; iLevel >= Constants::Thresholds::CellsMinLevel; --iLevel) {

    const int minimumLevel { iLevel - 1 };

    for (int iLayer { Constants::ITS::CellsPerRoad - 1 }; iLayer >= minimumLevel; --iLayer) {

      const int levelCellsNum { static_cast<int>(mPrimaryVertexContext.getCells()[iLayer].size()) };

      for (int iCell { 0 }; iCell < levelCellsNum; ++iCell) {

        Cell& currentCell { mPrimaryVertexContext.getCells()[iLayer][iCell] };

        if (currentCell.getLevel() != iLevel) {

          continue;
        }

        mPrimaryVertexContext.getRoads().emplace_back(iLayer, iCell);

        const int cellNeighboursNum {
            static_cast<int>(mPrimaryVertexContext.getCellsNeighbours()[iLayer - 1][iCell].size()) };
        bool isFirstValidNeighbour = true;

        for (int iNeighbourCell { 0 }; iNeighbourCell < cellNeighboursNum; ++iNeighbourCell) {

          const int neighbourCellId = mPrimaryVertexContext.getCellsNeighbours()[iLayer - 1][iCell][iNeighbourCell];
          const Cell& neighbourCell = mPrimaryVertexContext.getCells()[iLayer - 1][neighbourCellId];

          if (iLevel - 1 != neighbourCell.getLevel()) {

            continue;
          }

          if (isFirstValidNeighbour) {

            isFirstValidNeighbour = false;

          } else {

            mPrimaryVertexContext.getRoads().emplace_back(iLayer, iCell);
          }

          traverseCellsTree(neighbourCellId, iLayer - 1);
        }

        //TODO: crosscheck for short track iterations
        //currentCell.setLevel(0);
      }
    }
  }
}

template<bool IsGPU>
void Tracker<IsGPU>::traverseCellsTree(const int currentCellId, const int currentLayerId)
{
  Cell& currentCell { mPrimaryVertexContext.getCells()[currentLayerId][currentCellId] };
  const int currentCellLevel = currentCell.getLevel();

  mPrimaryVertexContext.getRoads().back().addCell(currentLayerId, currentCellId);

  if (currentLayerId > 0) {

    const int cellNeighboursNum {
        static_cast<int>(mPrimaryVertexContext.getCellsNeighbours()[currentLayerId - 1][currentCellId].size()) };
    bool isFirstValidNeighbour = true;

    for (int iNeighbourCell { 0 }; iNeighbourCell < cellNeighboursNum; ++iNeighbourCell) {

      const int neighbourCellId =
          mPrimaryVertexContext.getCellsNeighbours()[currentLayerId - 1][currentCellId][iNeighbourCell];
      const Cell& neighbourCell = mPrimaryVertexContext.getCells()[currentLayerId - 1][neighbourCellId];

      if (currentCellLevel - 1 != neighbourCell.getLevel()) {

        continue;
      }

      if (isFirstValidNeighbour) {

        isFirstValidNeighbour = false;

      } else {

        mPrimaryVertexContext.getRoads().push_back(mPrimaryVertexContext.getRoads().back());
      }

      traverseCellsTree(neighbourCellId, currentLayerId - 1);
    }
  }

  //TODO: crosscheck for short track iterations
  //currentCell.setLevel(0);
}

template<bool IsGPU>
void Tracker<IsGPU>::computeMontecarloLabels()
{
/// Moore’s Voting Algorithm

  int roadsNum { static_cast<int>(mPrimaryVertexContext.getRoads().size()) };

  for (int iRoad { 0 }; iRoad < roadsNum; ++iRoad) {

    Road& currentRoad { mPrimaryVertexContext.getRoads()[iRoad] };
    int maxOccurrencesValue { Constants::ITS::UnusedIndex };
    int count { 0 };
    bool isFakeRoad { false };
    bool isFirstRoadCell { true };

    for (int iCell { 0 }; iCell < Constants::ITS::CellsPerRoad; ++iCell) {

      const int currentCellIndex { currentRoad[iCell] };

      if (currentCellIndex == Constants::ITS::UnusedIndex) {

        if (isFirstRoadCell) {

          continue;

        } else {

          break;
        }
      }

      const Cell& currentCell { mPrimaryVertexContext.getCells()[iCell][currentCellIndex] };

      if (isFirstRoadCell) {

        maxOccurrencesValue =
            mPrimaryVertexContext.getClusters()[iCell][currentCell.getFirstClusterIndex()].monteCarloId;
        count = 1;

        const int secondMonteCarlo {
          mPrimaryVertexContext.getClusters()[iCell + 1][currentCell.getSecondClusterIndex()].monteCarloId };

        if (secondMonteCarlo == maxOccurrencesValue) {

          ++count;

        } else {

          maxOccurrencesValue = secondMonteCarlo;
          count = 1;
          isFakeRoad = true;
        }

        isFirstRoadCell = false;
      }

      const int currentMonteCarlo {
        mPrimaryVertexContext.getClusters()[iCell + 2][currentCell.getThirdClusterIndex()].monteCarloId };

      if (currentMonteCarlo == maxOccurrencesValue) {

        ++count;

      } else {

        --count;
        isFakeRoad = true;
      }

      if (count == 0) {

        maxOccurrencesValue = currentMonteCarlo;
        count = 1;
      }
    }

    currentRoad.setLabel(maxOccurrencesValue);
    currentRoad.setFakeRoad(isFakeRoad);
  }
}

template<bool IsGPU>
float Tracker<IsGPU>::evaluateTask(void (Tracker<IsGPU>::*task)(void), const char *taskName)
{
  return evaluateTask(task, taskName, std::cout);
}

template<bool IsGPU>
float Tracker<IsGPU>::evaluateTask(void (Tracker<IsGPU>::*task)(void), const char *taskName,
    std::ostream& ostream)
{
  clock_t t1, t2;
  float diff;

  t1 = clock();

  (this->*task)();

  t2 = clock();
  diff = ((float) t2 - (float) t1) / (CLOCKS_PER_SEC / 1000);

  if (taskName == nullptr) {

    ostream << diff << "\t";

  } else {

    ostream << std::setw(2) << " - " << taskName << " completed in: " << diff << "ms" << std::endl;
  }

  return diff;
}

template class Tracker<TRACKINGITSU_GPU_MODE> ;

}
}
}

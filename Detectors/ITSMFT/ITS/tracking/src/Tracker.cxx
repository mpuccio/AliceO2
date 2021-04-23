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

#include "ITStracking/Tracker.h"

#include "ITStracking/Cell.h"
#include "ITStracking/Constants.h"
#include "ITStracking/IndexTableUtils.h"
#include "ITStracking/Tracklet.h"
#include "ITStracking/TrackerTraits.h"
#include "ITStracking/TrackerTraitsCPU.h"
#include "ITStracking/TrackingConfigParam.h"

#include "ReconstructionDataFormats/Track.h"
#include <cassert>
#include <iostream>
#include <dlfcn.h>
#include <cstdlib>
#include <string>

namespace o2
{
namespace its
{

Tracker::Tracker(o2::its::TrackerTraits* traits)
{
  /// Initialise standard configuration with 1 iteration
  mTrkParams.resize(1);
  mMemParams.resize(1);
  mTraits = traits;
  mTimeFrame = mTraits->getTimeFrame();
#ifdef CA_DEBUG
  mDebugger = new StandaloneDebugger("dbg_ITSTrackerCPU.root");
#endif
}
#ifdef CA_DEBUG
Tracker::~Tracker()
{
  delete mDebugger;
}
#else
Tracker::~Tracker() = default;
#endif

void Tracker::clustersToTracks(std::ostream& timeBenchmarkOutputStream)
{
  mTracks.clear();
  mTrackLabels.clear();

  double total{0};
  for (int iteration = 0; iteration < mTrkParams.size(); ++iteration) {
    mTraits->UpdateTrackingParameters(mTrkParams[iteration]);

    total += evaluateTask(&Tracker::initialiseTimeFrame, "Context initialisation",
                          timeBenchmarkOutputStream, iteration, mMemParams[iteration], mTrkParams[iteration]);
    total += evaluateTask(&Tracker::computeTracklets, "Tracklet finding", timeBenchmarkOutputStream);
    mTimeFrame->printTrackletLUTs();
    total += evaluateTask(&Tracker::computeCells, "Cell finding", timeBenchmarkOutputStream);
    mTimeFrame->printCellLUTs();
    total += evaluateTask(&Tracker::findCellsNeighbours, "Neighbour finding", timeBenchmarkOutputStream, iteration);
    total += evaluateTask(&Tracker::findRoads, "Road finding", timeBenchmarkOutputStream, iteration);
    total += evaluateTask(&Tracker::findTracks, "Track finding", timeBenchmarkOutputStream);
  }

  if (constants::DoTimeBenchmarks)
    timeBenchmarkOutputStream << std::setw(2) << " - "
                              << "Vertex processing completed in: " << total << "ms" << std::endl;

  if (mTimeFrame->hasMCinformation()) {
    computeTracksMClabels();
  } else {
    rectifyClusterIndices();
  }
}

void Tracker::computeTracklets()
{
  mTraits->computeLayerTracklets();
}

void Tracker::computeCells()
{
  mTraits->computeLayerCells();
}

void Tracker::findCellsNeighbours(int& iteration)
{
  for (int iLayer{0}; iLayer < mTrkParams[iteration].CellsPerRoad() - 1; ++iLayer) {

    if (mTimeFrame->getCells()[iLayer + 1].empty() ||
        mTimeFrame->getCellsLookupTable()[iLayer].empty()) {
      continue;
    }

    int layerCellsNum{static_cast<int>(mTimeFrame->getCells()[iLayer].size())};
    const int nextLayerCellsNum{static_cast<int>(mTimeFrame->getCells()[iLayer + 1].size())};
    mTimeFrame->getCellsNeighbours()[iLayer].resize(nextLayerCellsNum);

    for (int iCell{0}; iCell < layerCellsNum; ++iCell) {

      const Cell& currentCell{mTimeFrame->getCells()[iLayer][iCell]};
      const int nextLayerTrackletIndex{currentCell.getSecondTrackletIndex()};
      std::cout << iLayer << "\t" << mTimeFrame->getCellsLookupTable()[iLayer].size() << "\t" << nextLayerTrackletIndex << std::endl;
      const int nextLayerFirstCellIndex{mTimeFrame->getCellsLookupTable()[iLayer][nextLayerTrackletIndex]};
      const int nextLayerLastCellIndex{mTimeFrame->getCellsLookupTable()[iLayer][nextLayerTrackletIndex + 1]};
      for (int iNextLayerCell{nextLayerFirstCellIndex}; iNextLayerCell < nextLayerLastCellIndex; ++iNextLayerCell) {

        Cell& nextCell{mTimeFrame->getCells()[iLayer + 1][iNextLayerCell]};
        if (nextCell.getFirstTrackletIndex() != nextLayerTrackletIndex) {
          std::cout << "Problem with the Cell LUT" << std::endl;
          break;
        }

        const float3 currentCellNormalVector{currentCell.getNormalVectorCoordinates()};
        const float3 nextCellNormalVector{nextCell.getNormalVectorCoordinates()};
        const float3 normalVectorsDeltaVector{currentCellNormalVector.x - nextCellNormalVector.x,
                                              currentCellNormalVector.y - nextCellNormalVector.y,
                                              currentCellNormalVector.z - nextCellNormalVector.z};

        const float deltaNormalVectorsModulus{(normalVectorsDeltaVector.x * normalVectorsDeltaVector.x) +
                                              (normalVectorsDeltaVector.y * normalVectorsDeltaVector.y) +
                                              (normalVectorsDeltaVector.z * normalVectorsDeltaVector.z)};
        const float deltaCurvature{std::abs(currentCell.getCurvature() - nextCell.getCurvature())};

        if (deltaNormalVectorsModulus < mTrkParams[iteration].NeighbourMaxDeltaN[iLayer] &&
            deltaCurvature < mTrkParams[iteration].NeighbourMaxDeltaCurvature[iLayer]) {

          mTimeFrame->getCellsNeighbours()[iLayer][iNextLayerCell].push_back(iCell);

          const int currentCellLevel{currentCell.getLevel()};

          if (currentCellLevel >= nextCell.getLevel()) {
            nextCell.setLevel(currentCellLevel + 1);
          }
        }
      }
    }
  }
}

void Tracker::findRoads(int& iteration)
{
  for (int iLevel{mTrkParams[iteration].CellsPerRoad()}; iLevel >= mTrkParams[iteration].CellMinimumLevel(); --iLevel) {
    CA_DEBUGGER(int nRoads = -mTimeFrame->getRoads().size());
    const int minimumLevel{iLevel - 1};

    for (int iLayer{mTrkParams[iteration].CellsPerRoad() - 1}; iLayer >= minimumLevel; --iLayer) {

      const int levelCellsNum{static_cast<int>(mTimeFrame->getCells()[iLayer].size())};

      for (int iCell{0}; iCell < levelCellsNum; ++iCell) {

        Cell& currentCell{mTimeFrame->getCells()[iLayer][iCell]};

        if (currentCell.getLevel() != iLevel) {
          continue;
        }

        mTimeFrame->getRoads().emplace_back(iLayer, iCell);

        /// For 3 clusters roads (useful for cascades and hypertriton) we just store the single cell
        /// and we do not do the candidate tree traversal
        if (iLevel == 1) {
          continue;
        }

        const int cellNeighboursNum{static_cast<int>(
          mTimeFrame->getCellsNeighbours()[iLayer - 1][iCell].size())};
        bool isFirstValidNeighbour = true;

        for (int iNeighbourCell{0}; iNeighbourCell < cellNeighboursNum; ++iNeighbourCell) {

          const int neighbourCellId = mTimeFrame->getCellsNeighbours()[iLayer - 1][iCell][iNeighbourCell];
          const Cell& neighbourCell = mTimeFrame->getCells()[iLayer - 1][neighbourCellId];

          if (iLevel - 1 != neighbourCell.getLevel()) {
            continue;
          }

          if (isFirstValidNeighbour) {

            isFirstValidNeighbour = false;

          } else {

            mTimeFrame->getRoads().emplace_back(iLayer, iCell);
          }

          traverseCellsTree(neighbourCellId, iLayer - 1);
        }

        // TODO: crosscheck for short track iterations
        // currentCell.setLevel(0);
      }
    }
#ifdef CA_DEBUG
    nRoads += mTimeFrame->getRoads().size();
    std::cout << "+++ Roads with " << iLevel + 2 << " clusters: " << nRoads << " / " << mTimeFrame->getRoads().size() << std::endl;
#endif
  }
  std::cout << "Number of roads: " << mTimeFrame->getRoads().size() << std::endl;
}

void Tracker::findTracks()
{
  mTracks.resize(mTimeFrame->getNrof());
  std::vector<TrackITSExt> tracks;
  tracks.reserve(mTimeFrame->getRoads().size());
#ifdef CA_DEBUG
  std::vector<int> roadCounters(mTrkParams[0].NLayers - 3, 0);
  std::vector<int> fitCounters(mTrkParams[0].NLayers - 3, 0);
  std::vector<int> backpropagatedCounters(mTrkParams[0].NLayers - 3, 0);
  std::vector<int> refitCounters(mTrkParams[0].NLayers - 3, 0);
  std::vector<int> nonsharingCounters(mTrkParams[0].NLayers - 3, 0);
#endif

  for (auto& road : mTimeFrame->getRoads()) {
    std::vector<int> clusters(mTrkParams[0].NLayers, constants::its::UnusedIndex);
    int lastCellLevel = constants::its::UnusedIndex;
    CA_DEBUGGER(int nClusters = 2);
    int firstTracklet{constants::its::UnusedIndex};
    std::vector<int> tracklets(mTrkParams[0].TrackletsPerRoad(), constants::its::UnusedIndex);

    for (int iCell{0}; iCell < mTrkParams[0].CellsPerRoad(); ++iCell) {
      const int cellIndex = road[iCell];
      if (cellIndex == constants::its::UnusedIndex) {
        continue;
      } else {
        if (firstTracklet == constants::its::UnusedIndex)
          firstTracklet = iCell;
        tracklets[iCell] = mTimeFrame->getCells()[iCell][cellIndex].getFirstTrackletIndex();
        tracklets[iCell + 1] = mTimeFrame->getCells()[iCell][cellIndex].getSecondTrackletIndex();
        clusters[iCell] = mTimeFrame->getCells()[iCell][cellIndex].getFirstClusterIndex();
        clusters[iCell + 1] = mTimeFrame->getCells()[iCell][cellIndex].getSecondClusterIndex();
        clusters[iCell + 2] = mTimeFrame->getCells()[iCell][cellIndex].getThirdClusterIndex();
        assert(clusters[iCell] != constants::its::UnusedIndex &&
               clusters[iCell + 1] != constants::its::UnusedIndex &&
               clusters[iCell + 2] != constants::its::UnusedIndex);
        lastCellLevel = iCell;
        CA_DEBUGGER(nClusters++);
      }
    }

    CA_DEBUGGER(assert(nClusters >= mTrkParams[0].MinTrackLength));
    int count{1};
    unsigned short rof{mTimeFrame->getTracklets()[firstTracklet][tracklets[firstTracklet]].rof[0]};
    for (int iT = firstTracklet; iT < 6; ++iT) {
      if (tracklets[iT] == constants::its::UnusedIndex)
        continue;
      if (rof == mTimeFrame->getTracklets()[iT][tracklets[iT]].rof[1])
        count++;
      else {
        if (count == 1)
          rof = mTimeFrame->getTracklets()[iT][tracklets[iT]].rof[1];
        else
          count--;
      }
    }

    assert(nClusters >= mTrkParams[0].MinTrackLength);
    CA_DEBUGGER(roadCounters[nClusters - 4]++);

    if (lastCellLevel == constants::its::UnusedIndex) {
      continue;
    }

    /// From primary vertex context index to event index (== the one used as input of the tracking code)
    for (int iC{0}; iC < clusters.size(); iC++) {
      if (clusters[iC] != constants::its::UnusedIndex) {
        clusters[iC] = mTimeFrame->getClusters()[iC][clusters[iC]].clusterId;
      }
    }

    /// Track seed preparation. Clusters are numbered progressively from the outermost to the innermost.
    const auto& cluster1_glo = mTimeFrame->getUnsortedClusters()[lastCellLevel + 2].at(clusters[lastCellLevel + 2]);
    const auto& cluster2_glo = mTimeFrame->getUnsortedClusters()[lastCellLevel + 1].at(clusters[lastCellLevel + 1]);
    const auto& cluster3_glo = mTimeFrame->getUnsortedClusters()[lastCellLevel].at(clusters[lastCellLevel]);

    const auto& cluster3_tf = mTimeFrame->getTrackingFrameInfoOnLayer(lastCellLevel).at(clusters[lastCellLevel]);

    /// FIXME!
    TrackITSExt temporaryTrack{buildTrackSeed(cluster1_glo, cluster2_glo, cluster3_glo, cluster3_tf)};
    for (size_t iC = 0; iC < clusters.size(); ++iC) {
      temporaryTrack.setExternalClusterIndex(iC, clusters[iC], clusters[iC] != constants::its::UnusedIndex);
    }
    bool fitSuccess = fitTrack(temporaryTrack, mTrkParams[0].NLayers - 4, -1, -1);
    if (!fitSuccess) {
      continue;
    }
    CA_DEBUGGER(fitCounters[nClusters - 4]++);
    temporaryTrack.resetCovariance();
    fitSuccess = fitTrack(temporaryTrack, 0, mTrkParams[0].NLayers, 1, mTrkParams[0].FitIterationMaxChi2[0]);
    if (!fitSuccess) {
      continue;
    }
    CA_DEBUGGER(backpropagatedCounters[nClusters - 4]++);
    temporaryTrack.getParamOut() = temporaryTrack;
    temporaryTrack.resetCovariance();
    fitSuccess = fitTrack(temporaryTrack, mTrkParams[0].NLayers - 1, -1, -1, mTrkParams[0].FitIterationMaxChi2[1]);
#ifdef CA_DEBUG
    mDebugger->dumpTrackToBranchWithInfo("testBranch", temporaryTrack, event, mPrimaryVertexContext, true);
#endif
    if (!fitSuccess) {
      continue;
    }
    CA_DEBUGGER(refitCounters[nClusters - 4]++);
    // temporaryTrack.setROFrame(rof);
    tracks.emplace_back(temporaryTrack);
    CA_DEBUGGER(assert(nClusters == temporaryTrack.getNumberOfClusters()));
  }
  std::cout << "Fitted tracks: " << tracks.size() << std::endl;
  //mTraits->refitTracks(event.getTrackingFrameInfo(), tracks);

  std::sort(tracks.begin(), tracks.end(),
            [](TrackITSExt& track1, TrackITSExt& track2) { return track1.isBetter(track2, 1.e6f); });

#ifdef CA_DEBUG
  // std::array<int, 26> sharingMatrix{0};
  // int prevNclusters = 7;
  // auto cumulativeIndex = [](int ncl) -> int {
  //   constexpr int idx[5] = {0, 5, 11, 18, 26};
  //   return idx[ncl - 4];
  // };
  // std::array<int, 4> xcheckCounters{0};
#endif

  for (auto& track : tracks) {
    CA_DEBUGGER(int nClusters = 0);
    int nShared = 0;
    for (int iLayer{0}; iLayer < mTrkParams[0].NLayers; ++iLayer) {
      if (track.getClusterIndex(iLayer) == constants::its::UnusedIndex) {
        continue;
      }
      nShared += int(mTimeFrame->isClusterUsed(iLayer, track.getClusterIndex(iLayer)));
      CA_DEBUGGER(nClusters++);
    }

    // #ifdef CA_DEBUG
    //     assert(nClusters == track.getNumberOfClusters());
    //     xcheckCounters[nClusters - 4]++;
    //     assert(nShared <= nClusters);
    //     sharingMatrix[cumulativeIndex(nClusters) + nShared]++;
    // #endif

    if (nShared > mTrkParams[0].ClusterSharing) {
      continue;
    }

    // #ifdef CA_DEBUG
    //     nonsharingCounters[nClusters - 4]++;
    //     assert(nClusters <= prevNclusters);
    //     prevNclusters = nClusters;
    // #endif

    for (int iLayer{0}; iLayer < mTrkParams[0].NLayers; ++iLayer) {
      if (track.getClusterIndex(iLayer) == constants::its::UnusedIndex) {
        continue;
      }
      mTimeFrame->markUsedCluster(iLayer, track.getClusterIndex(iLayer));
    }
    mTracks.emplace_back(track);
  }

#ifdef CA_DEBUG
  std::cout << "+++ Found candidates with 4, 5, 6 and 7 clusters:\t";
  for (int count : roadCounters)
    std::cout << count << "\t";
  std::cout << std::endl;

  std::cout << "+++ Fitted candidates with 4, 5, 6 and 7 clusters:\t";
  for (int count : fitCounters)
    std::cout << count << "\t";
  std::cout << std::endl;

  std::cout << "+++ Backprop candidates with 4, 5, 6 and 7 clusters:\t";
  for (int count : backpropagatedCounters)
    std::cout << count << "\t";
  std::cout << std::endl;

  std::cout << "+++ Refitted candidates with 4, 5, 6 and 7 clusters:\t";
  for (int count : refitCounters)
    std::cout << count << "\t";
  std::cout << std::endl;

  // std::cout << "+++ Cross check counters for 4, 5, 6 and 7 clusters:\t";
  // for (size_t iCount = 0; iCount < refitCounters.size(); ++iCount) {
  //   std::cout << xcheckCounters[iCount] << "\t";
  //   //assert(refitCounters[iCount] == xcheckCounters[iCount]);
  // }
  // std::cout << std::endl;

  // std::cout << "+++ Nonsharing candidates with 4, 5, 6 and 7 clusters:\t";
  // for (int count : nonsharingCounters)
  //   std::cout << count << "\t";
  // std::cout << std::endl;

  // std::cout << "+++ Sharing matrix:\n";
  // for (int iCl = 4; iCl <= 7; ++iCl) {
  //   std::cout << "+++ ";
  //   for (int iSh = cumulativeIndex(iCl); iSh < cumulativeIndex(iCl + 1); ++iSh) {
  //     std::cout << sharingMatrix[iSh] << "\t";
  //   }
  //   std::cout << std::endl;
  // }
#endif
}

bool Tracker::fitTrack(TrackITSExt& track, int start, int end, int step, const float chi2cut)
{
  auto propInstance = o2::base::Propagator::Instance();
  track.setChi2(0);
  for (int iLayer{start}; iLayer != end; iLayer += step) {
    if (track.getClusterIndex(iLayer) == constants::its::UnusedIndex) {
      continue;
    }
    const TrackingFrameInfo& trackingHit = mTimeFrame->getTrackingFrameInfoOnLayer(iLayer).at(track.getClusterIndex(iLayer));

    if (!track.rotate(trackingHit.alphaTrackingFrame)) {
      return false;
    }

    if (!propInstance->propagateToX(track, trackingHit.xTrackingFrame, getBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, mCorrType)) {
      return false;
    }

    auto predChi2{track.getPredictedChi2(trackingHit.positionTrackingFrame, trackingHit.covarianceTrackingFrame)};
    if (predChi2 > chi2cut) {
      return false;
    }
    track.setChi2(track.getChi2() + predChi2);
    if (!track.o2::track::TrackParCov::update(trackingHit.positionTrackingFrame, trackingHit.covarianceTrackingFrame)) {
      return false;
    }
  }
  return true;
}

void Tracker::traverseCellsTree(const int currentCellId, const int currentLayerId)
{
  Cell& currentCell{mTimeFrame->getCells()[currentLayerId][currentCellId]};
  const int currentCellLevel = currentCell.getLevel();

  mTimeFrame->getRoads().back().addCell(currentLayerId, currentCellId);

  if (currentLayerId > 0 && currentCellLevel > 1) {
    const int cellNeighboursNum{static_cast<int>(
      mTimeFrame->getCellsNeighbours()[currentLayerId - 1][currentCellId].size())};
    bool isFirstValidNeighbour = true;

    for (int iNeighbourCell{0}; iNeighbourCell < cellNeighboursNum; ++iNeighbourCell) {

      const int neighbourCellId =
        mTimeFrame->getCellsNeighbours()[currentLayerId - 1][currentCellId][iNeighbourCell];
      const Cell& neighbourCell = mTimeFrame->getCells()[currentLayerId - 1][neighbourCellId];

      if (currentCellLevel - 1 != neighbourCell.getLevel()) {
        continue;
      }

      if (isFirstValidNeighbour) {
        isFirstValidNeighbour = false;
      } else {
        mTimeFrame->getRoads().push_back(mTimeFrame->getRoads().back());
      }

      traverseCellsTree(neighbourCellId, currentLayerId - 1);
    }
  }

  // TODO: crosscheck for short track iterations
  // currentCell.setLevel(0);
}

void Tracker::computeRoadsMClabels()
{
  /// Moore's Voting Algorithm
  if (!mTimeFrame->hasMCinformation()) {
    return;
  }

  mTimeFrame->initialiseRoadLabels();

  int roadsNum{static_cast<int>(mTimeFrame->getRoads().size())};

  for (int iRoad{0}; iRoad < roadsNum; ++iRoad) {

    Road& currentRoad{mTimeFrame->getRoads()[iRoad]};
    unsigned long long maxOccurrencesValue{0xffffffffffffffff};
    int count{0};
    bool isFakeRoad{false};
    bool isFirstRoadCell{true};

    for (int iCell{0}; iCell < mTrkParams[0].CellsPerRoad(); ++iCell) {
      const int currentCellIndex{currentRoad[iCell]};

      if (currentCellIndex == constants::its::UnusedIndex) {
        if (isFirstRoadCell) {
          continue;
        } else {
          break;
        }
      }

      const Cell& currentCell{mTimeFrame->getCells()[iCell][currentCellIndex]};

      if (isFirstRoadCell) {

        const int cl0index{mTimeFrame->getClusters()[iCell][currentCell.getFirstClusterIndex()].clusterId};
        auto& cl0labs{mTimeFrame->getClusterLabels(iCell, cl0index)};
        maxOccurrencesValue = cl0labs;
        count = 1;

        const int cl1index{mTimeFrame->getClusters()[iCell + 1][currentCell.getSecondClusterIndex()].clusterId};
        auto& cl1labs{mTimeFrame->getClusterLabels(iCell + 1, cl1index)};

        if (cl1labs == maxOccurrencesValue) {
          ++count;
        } else {
          maxOccurrencesValue = cl1labs;
          count = 1;
          isFakeRoad = true;
        }

        isFirstRoadCell = false;
      }

      const int cl2index{mTimeFrame->getClusters()[iCell + 2][currentCell.getThirdClusterIndex()].clusterId};
      auto& cl2labs{mTimeFrame->getClusterLabels(iCell + 2, cl2index)};

      if (cl2labs == maxOccurrencesValue) {
        ++count;
      } else {
        --count;
        isFakeRoad = true;
      }

      if (count == 0) {
        maxOccurrencesValue = cl2labs;
        count = 1;
      }
    }

    mTimeFrame->setRoadLabel(iRoad, maxOccurrencesValue, isFakeRoad);
  }
}

void Tracker::computeTracksMClabels()
{
  /// Moore's Voting Algorithm
  if (!mTimeFrame->hasMCinformation()) {
    return;
  }

  int tracksNum{static_cast<int>(mTracks.size())};
  mTrackLabels.resize(mTracks.size());

  for (auto& track : mTracks) {

    unsigned long long maxOccurrencesValue{0xffffffffffffffff};
    int count{0};
    bool isFakeTrack{false};

    for (int iCluster = 0; iCluster < TrackITSExt::MaxClusters; ++iCluster) {
      const int index = track.getClusterIndex(iCluster);
      if (index == constants::its::UnusedIndex) {
        continue;
      }
      unsigned long long currentLabel = mTimeFrame->getClusterLabels(iCluster, index);
      if (currentLabel == maxOccurrencesValue) {
        ++count;
      } else {
        if (count != 0) { // only in the first iteration count can be 0 at this point
          isFakeTrack = true;
          --count;
        }

        unsigned long long currentLabel = mTimeFrame->getClusterLabels(iCluster, index);
        if (currentLabel == maxOccurrencesValue) {
          ++count;
        } else {
          if (count != 0) { // only in the first iteration count can be 0 at this point
            isFakeTrack = true;
            --count;
          }
          if (count == 0) {
            maxOccurrencesValue = currentLabel;
            count = 1;
          }
        }
        track.setExternalClusterIndex(iCluster, mTimeFrame->getClusterExternalIndex(iCluster, index));
      }

      if (isFakeTrack) {
        maxOccurrencesValue |= static_cast<unsigned long long>(0x1) << 63;
      }
    }
    mTrackLabels.emplace_back(maxOccurrencesValue);
  }
}

void Tracker::rectifyClusterIndices()
{
  int tracksNum{static_cast<int>(mTracks.size())};
  for (auto& track : mTracks) {
    for (int iCluster = 0; iCluster < TrackITSExt::MaxClusters; ++iCluster) {
      const int index = track.getClusterIndex(iCluster);
      if (index != constants::its::UnusedIndex) {
        track.setExternalClusterIndex(iCluster, mTimeFrame->getClusterExternalIndex(iCluster, index));
      }
    }
  }
}

/// Clusters are given from outside inward (cluster1 is the outermost). The innermost cluster is given in the tracking
/// frame coordinates
/// whereas the others are referred to the global frame. This function is almost a clone of CookSeed, adapted to return
/// a TrackParCov
track::TrackParCov Tracker::buildTrackSeed(const Cluster& cluster1, const Cluster& cluster2,
                                           const Cluster& cluster3, const TrackingFrameInfo& tf3)
{
  const float ca = std::cos(tf3.alphaTrackingFrame), sa = std::sin(tf3.alphaTrackingFrame);
  const float x1 = cluster1.xCoordinate * ca + cluster1.yCoordinate * sa;
  const float y1 = -cluster1.xCoordinate * sa + cluster1.yCoordinate * ca;
  const float z1 = cluster1.zCoordinate;
  const float x2 = cluster2.xCoordinate * ca + cluster2.yCoordinate * sa;
  const float y2 = -cluster2.xCoordinate * sa + cluster2.yCoordinate * ca;
  const float z2 = cluster2.zCoordinate;
  const float x3 = tf3.xTrackingFrame;
  const float y3 = tf3.positionTrackingFrame[0];
  const float z3 = tf3.positionTrackingFrame[1];

  const float crv = math_utils::computeCurvature(x1, y1, x2, y2, x3, y3);
  const float x0 = math_utils::computeCurvatureCentreX(x1, y1, x2, y2, x3, y3);
  const float tgl12 = math_utils::computeTanDipAngle(x1, y1, x2, y2, z1, z2);
  const float tgl23 = math_utils::computeTanDipAngle(x2, y2, x3, y3, z2, z3);

  const float fy = 1. / (cluster2.radius - cluster3.radius);
  const float& tz = fy;
  const float cy = (math_utils::computeCurvature(x1, y1, x2, y2 + constants::its::Resolution, x3, y3) - crv) /
                   (constants::its::Resolution * getBz() * o2::constants::math::B2C) *
                   20.f; // FIXME: MS contribution to the cov[14] (*20 added)
  constexpr float s2 = constants::its::Resolution * constants::its::Resolution;

  return track::TrackParCov(tf3.xTrackingFrame, tf3.alphaTrackingFrame,
                            {y3, z3, crv * (x3 - x0), 0.5f * (tgl12 + tgl23),
                             std::abs(getBz()) < o2::constants::math::Almost0 ? o2::constants::math::Almost0
                                                                              : crv / (getBz() * o2::constants::math::B2C)},
                            {s2, 0.f, s2, s2 * fy, 0.f, s2 * fy * fy, 0.f, s2 * tz, 0.f, s2 * tz * tz, s2 * cy, 0.f,
                             s2 * fy * cy, 0.f, s2 * cy * cy});
}

void Tracker::getGlobalConfiguration()
{
  auto& tc = o2::its::TrackerParamConfig::Instance();
  if (tc.useMatCorrTGeo) {
    setCorrType(o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrTGeo);
  }
}

} // namespace its
} // namespace o2

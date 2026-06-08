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
/// \file CATracker.cxx
/// \brief
///

#include "ITSMFTTracking/CATrackTypes.h"
#include "ITSMFTTracking/CATracker.h"

#include <algorithm>
#include <numeric>
#include <vector>

#include "Framework/Logger.h"
#include "ITStracking/Constants.h"
#include "ITStracking/BoundedAllocator.h"
#include "MFTTracking/MFTTrackingParam.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::itsmft::tracking
{

namespace
{
template <int NLayers>
void computeTracksMClabels(TimeFrame<NLayers>& tf)
{
  auto& trackLabels = tf.getTracksLabel();
  trackLabels.clear();
  trackLabels.reserve(tf.getNumberOfTracks());

  for (auto& track : tf.getTracks()) {
    std::vector<std::pair<MCCompLabel, size_t>> occurrences;
    for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
      const int index = track.getClusterIndex(iLayer);
      if (index == constants::UnusedIndex) {
        continue;
      }
      const auto labels = tf.getClusterLabels(iLayer, index);
      bool found{false};
      for (auto& occurrence : occurrences) {
        for (const auto& label : labels) {
          if (label == occurrence.first) {
            ++occurrence.second;
            found = true;
          }
        }
      }
      if (!found) {
        for (const auto& label : labels) {
          occurrences.emplace_back(label, 1);
        }
      }
    }

    MCCompLabel maxOccurrencesValue;
    if (!occurrences.empty()) {
      std::sort(occurrences.begin(), occurrences.end(), [](const auto& e1, const auto& e2) {
        return e1.second > e2.second;
      });
      maxOccurrencesValue = occurrences[0].first;
      if constexpr (NLayers == constants::MFTNLayers) {
        const float threshold = o2::mft::MFTTrackingParam::Instance().TrueTrackMCThreshold;
        if (static_cast<float>(occurrences[0].second) / track.getNumberOfClusters() < threshold) {
          maxOccurrencesValue.setFakeFlag();
        }
      } else if (occurrences[0].second < static_cast<size_t>(track.getNumberOfClusters())) {
        maxOccurrencesValue.setFakeFlag();
      }
    } else {
      maxOccurrencesValue.setFakeFlag();
    }
    trackLabels.emplace_back(maxOccurrencesValue);
  }
}
} // namespace

template <int NLayers>
Tracker<NLayers>::Tracker(TrackerTraits<NLayers>* traits) : mTraits(traits)
{
}

template <int NLayers>
void Tracker<NLayers>::adoptTimeFrame(TimeFrame<NLayers>& tf)
{
  mTimeFrame = &tf;
  mTraits->adoptTimeFrame(&tf);
}

template <int NLayers>
float Tracker<NLayers>::clustersToTracks()
{
  mTraits->updateTrackingParameters(mTrkParams);

  int maxNvertices{-1};
  if (mTrkParams[0].PerPrimaryVertexProcessing) {
    maxNvertices = mTimeFrame->getROFVertexLookupTableView().getMaxVerticesPerROF();
  }

  float total{0.f};
  try {
    for (int iteration = 0; iteration < static_cast<int>(mTrkParams.size()); ++iteration) {
      mMemoryPool->setMaxMemory(mTrkParams[iteration].MaxMemory);
      if (mTrkParams[iteration].PassFlags[IterationStep::UseUPCMask]) {
        mTimeFrame->useUPCMask();
      }

      int iVertex = std::min(maxNvertices, 0);
      initialiseTimeFrame(iteration);
      do {
        computeTracklets(iteration, iVertex);
        computeCells(iteration);
        findCellsNeighbours(iteration);
        findRoads(iteration);
      } while (++iVertex < maxNvertices);
    }
  } catch (const BoundedMemoryResource::MemoryLimitExceeded& err) {
    LOGP(error, "CA tracker exceeded memory limit: {}", err.what());
    if (mTrkParams[0].DropTFUponFailure) {
      mTimeFrame->wipe();
      return -1.f;
    }
    throw;
  } catch (const std::exception& err) {
    LOGP(error, "CA tracker failed: {}", err.what());
    mTimeFrame->getTracks().clear();
    return -1.f;
  }

  if (mTimeFrame->hasMCinformation()) {
    computeTracksMClabels(*mTimeFrame);
  }
  rectifyClusterIndices();
  sortTracks();
  return total;
}

template <int NLayers>
void Tracker<NLayers>::rectifyClusterIndices()
{
  for (auto& track : mTimeFrame->getTracks()) {
    for (int iCluster = 0; iCluster < CATrackType<NLayers>::MaxClusters; ++iCluster) {
      const int index = track.getClusterIndex(iCluster);
      if (index == constants::UnusedIndex) {
        continue;
      }
      track.setExternalClusterIndex(iCluster, mTimeFrame->getClusterExternalIndex(iCluster, index));
    }
  }
}

template <int NLayers>
void Tracker<NLayers>::sortTracks()
{
  auto& tracks = mTimeFrame->getTracks();
  bounded_vector<size_t> indices(tracks.size(), mMemoryPool.get());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&tracks](size_t i, size_t j) {
    const auto& a = tracks[i];
    const auto& b = tracks[j];
    const auto aLower = a.getTimeStamp().getTimeStamp() - a.getTimeStamp().getTimeStampError();
    const auto bLower = b.getTimeStamp().getTimeStamp() - b.getTimeStamp().getTimeStampError();
    if (aLower != bLower) {
      return aLower < bLower;
    }
    return a.getChi2() < b.getChi2();
  });

  bounded_vector<CATrackType<NLayers>> sortedTracks(mMemoryPool.get());
  sortedTracks.reserve(tracks.size());
  for (size_t idx : indices) {
    sortedTracks.push_back(tracks[idx]);
  }
  tracks.swap(sortedTracks);

  if (mTimeFrame->hasMCinformation()) {
    auto& trackLabels = mTimeFrame->getTracksLabel();
    bounded_vector<MCCompLabel> sortedLabels(mMemoryPool.get());
    sortedLabels.reserve(trackLabels.size());
    for (size_t idx : indices) {
      sortedLabels.push_back(trackLabels[idx]);
    }
    trackLabels.swap(sortedLabels);
  }
}

template class Tracker<7>;
template class Tracker<10>;

} // namespace o2::itsmft::tracking

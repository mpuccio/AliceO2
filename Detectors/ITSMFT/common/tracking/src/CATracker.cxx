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

#include "ITSMFTTracking/CATracker.h"

#include <algorithm>
#include <numeric>

#include "Framework/Logger.h"
#include "ITStracking/Constants.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITSMFTTracking/MCLabelAccumulator.h"
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
    MCLabelAccumulator labels;
    for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
      const int index = track.getClusterIndex(iLayer);
      if (index == constants::UnusedIndex) {
        continue;
      }
      labels.addCluster(tf.getClusterLabels(iLayer, index));
    }
    trackLabels.emplace_back(labels.finalize());
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
  } catch (const TraversalException& err) {
    // Structural/configuration failure (bad layout, stale layout, policy or
    // index mismatch): never a per-TF data problem, so DropTFUponFailure
    // never applies. Always wipe before propagating -- see class-level
    // comment: never rely on "the process is going down anyway".
    LOGP(error, "CA tracker hit a structural traversal failure: {}", err.what());
    mTimeFrame->wipe();
    throw;
  } catch (const BoundedMemoryResource::MemoryLimitExceeded& err) {
    // Recoverable, per-TF resource failure: the bounded pool's configured
    // budget was exceeded for this TimeFrame's data volume.
    LOGP(error, "CA tracker exceeded memory limit: {}", err.what());
    mTimeFrame->wipe();
    if (mTrkParams[0].DropTFUponFailure) {
      return kDroppedTimeFrameResult;
    }
    throw;
  } catch (const std::bad_alloc& err) {
    // Also recoverable/per-TF: several CA scratch containers on the hot path
    // (e.g. TrackerTraits::findCellsNeighbours' cellsNeighboursByTarget and
    // its per-thread TBB storage) allocate from the plain heap rather than
    // the bounded pool, so genuine memory pressure surfaces here as a plain
    // bad_alloc rather than MemoryLimitExceeded. Handled identically.
    LOGP(error, "CA tracker allocation failed: {}", err.what());
    mTimeFrame->wipe();
    if (mTrkParams[0].DropTFUponFailure) {
      return kDroppedTimeFrameResult;
    }
    throw;
  } catch (const std::exception& err) {
    // Unclassified: not a recognized recoverable-resource failure, so it is
    // treated as structural and always propagates, regardless of
    // DropTFUponFailure. A future explicitly typed
    // RecoverableTimeFrameException may extend the recoverable set; until
    // then, recoverability is never inferred from std::exception alone.
    LOGP(error, "CA tracker failed with an unclassified exception; treating as structural: {}", err.what());
    mTimeFrame->wipe();
    throw;
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

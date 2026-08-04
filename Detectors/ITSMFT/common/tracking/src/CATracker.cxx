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
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::itsmft::tracking
{

namespace
{
template <int NLayers, typename ScratchT>
void computeTracksMClabels(ScratchT& scratch)
{
  auto& trackLabels = scratch.getTracksLabel();
  trackLabels.clear();
  trackLabels.reserve(scratch.getNumberOfTracks());

  for (auto& track : scratch.getTracks()) {
    MCLabelAccumulator labels;
    for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
      const int index = track.getClusterIndex(iLayer);
      if (index == constants::UnusedIndex) {
        continue;
      }
      labels.addCluster(scratch.getClusterLabels(iLayer, index));
    }
    trackLabels.emplace_back(labels.finalize());
  }
}
} // namespace

template <int NLayers, typename ScratchT, typename BindingT>
Tracker<NLayers, ScratchT, BindingT>::Tracker(TrackerTraitsN* traits) : mTraits(traits)
{
}

template <int NLayers, typename ScratchT, typename BindingT>
void Tracker<NLayers, ScratchT, BindingT>::adoptScratch(ScratchT& scratch)
{
  mScratch = &scratch;
  mTraits->adoptScratch(&scratch);
}

template <int NLayers, typename ScratchT, typename BindingT>
void Tracker<NLayers, ScratchT, BindingT>::adoptFrame(TimeFrame& frame)
{
  mFrame = &frame;
  mTraits->adoptFrame(&frame);
}

template <int NLayers, typename ScratchT, typename BindingT>
TrackingResult Tracker<NLayers, ScratchT, BindingT>::clustersToTracks()
{
  mTraits->updateTrackingParameters(mTrkParams);

  int maxNvertices{-1};
  if (mTrkParams[0].PerPrimaryVertexProcessing) {
    maxNvertices = mScratch->getROFVertexLookupTableView().getMaxVerticesPerROF();
  }

  float total{0.f};
  try {
    for (int iteration = 0; iteration < static_cast<int>(mTrkParams.size()); ++iteration) {
      mMemoryPool->setMaxMemory(mTrkParams[iteration].MaxMemory);
      if (mTrkParams[iteration].PassFlags[IterationStep::UseUPCMask]) {
        mScratch->useUPCMask();
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
    // never applies. Always reset before propagating -- see class-level
    // comment: never rely on "the process is going down anyway".
    LOGP(error, "CA tracker hit a structural traversal failure: {}", err.what());
    if (mMFTPublicationCompatibility != nullptr) {
      mMFTPublicationCompatibility->clear();
    }
    if (mITSSharedClusterCompatibility != nullptr) {
      mITSSharedClusterCompatibility->clear();
    }
    resetTimeFrameEvent(*mFrame, *mScratch);
    throw;
  } catch (const BoundedMemoryResource::MemoryLimitExceeded& err) {
    // Recoverable, per-TF resource failure: the bounded pool's configured
    // budget was exceeded for this TimeFrame's data volume.
    LOGP(error, "CA tracker exceeded memory limit: {}", err.what());
    if (mMFTPublicationCompatibility != nullptr) {
      mMFTPublicationCompatibility->clear();
    }
    if (mITSSharedClusterCompatibility != nullptr) {
      mITSSharedClusterCompatibility->clear();
    }
    resetTimeFrameEvent(*mFrame, *mScratch);
    if (mTrkParams[0].DropTFUponFailure) {
      return TrackingResult{TrackingOutcome::RecoverableDropped, 0.f};
    }
    throw;
  } catch (const std::bad_alloc& err) {
    // Also recoverable/per-TF: several CA scratch containers on the hot path
    // (e.g. TrackerTraits::findCellsNeighbours' cellsNeighboursByTarget and
    // its per-thread TBB storage) allocate from the plain heap rather than
    // the bounded pool, so genuine memory pressure surfaces here as a plain
    // bad_alloc rather than MemoryLimitExceeded. Handled identically.
    LOGP(error, "CA tracker allocation failed: {}", err.what());
    if (mMFTPublicationCompatibility != nullptr) {
      mMFTPublicationCompatibility->clear();
    }
    if (mITSSharedClusterCompatibility != nullptr) {
      mITSSharedClusterCompatibility->clear();
    }
    resetTimeFrameEvent(*mFrame, *mScratch);
    if (mTrkParams[0].DropTFUponFailure) {
      return TrackingResult{TrackingOutcome::RecoverableDropped, 0.f};
    }
    throw;
  } catch (const std::exception& err) {
    // Unclassified: not a recognized recoverable-resource failure, so it is
    // treated as structural and always propagates, regardless of
    // DropTFUponFailure. A future explicitly typed
    // RecoverableTimeFrameException may extend the recoverable set; until
    // then, recoverability is never inferred from std::exception alone.
    LOGP(error, "CA tracker failed with an unclassified exception; treating as structural: {}", err.what());
    if (mMFTPublicationCompatibility != nullptr) {
      mMFTPublicationCompatibility->clear();
    }
    if (mITSSharedClusterCompatibility != nullptr) {
      mITSSharedClusterCompatibility->clear();
    }
    resetTimeFrameEvent(*mFrame, *mScratch);
    throw;
  }

  if (mScratch->hasMCinformation()) {
    computeTracksMClabels<NLayers>(*mScratch);
  }
  rectifyClusterIndices();
  sortTracks();
  return TrackingResult{TrackingOutcome::Success, total};
}

template <int NLayers, typename ScratchT, typename BindingT>
void Tracker<NLayers, ScratchT, BindingT>::rectifyClusterIndices()
{
  for (auto& track : mScratch->getTracks()) {
    for (int iCluster = 0; iCluster < CATrackType<NLayers>::MaxClusters; ++iCluster) {
      const int index = track.getClusterIndex(iCluster);
      if (index == constants::UnusedIndex) {
        continue;
      }
      // Capture the packed cluster size onto the track while `index` is
      // still this layer's own local identity (the domain mClusterSize is
      // stored in, see LegacyTrackerScratch::getClusterSize()): the very
      // next call overwrites track's cluster index in place with the
      // external/global identity, so this is the last point at which the
      // layer-local index needed to address mClusterSize[iCluster] is
      // recoverable from the track. Downstream publication
      // (TrackITSExt -> TrackITS, MFTCATrack) must read the size already
      // stored here rather than re-deriving it from the (by-then external)
      // cluster index.
      track.setClusterSize(iCluster, mScratch->getClusterSize(iCluster, index));
      track.setExternalClusterIndex(iCluster, mScratch->getClusterExternalIndex(iCluster, index));
    }
  }
}

template <int NLayers, typename ScratchT, typename BindingT>
void Tracker<NLayers, ScratchT, BindingT>::sortTracks()
{
  auto& tracks = mScratch->getTracks();
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

  if (mScratch->hasMCinformation()) {
    auto& trackLabels = mScratch->getTracksLabel();
    bounded_vector<MCCompLabel> sortedLabels(mMemoryPool.get());
    sortedLabels.reserve(trackLabels.size());
    for (size_t idx : indices) {
      sortedLabels.push_back(trackLabels[idx]);
    }
    trackLabels.swap(sortedLabels);
  }
}

// M6d: three explicit instantiations -- ITS, the standalone-MFT-workflow's
// own ITSMFTTrackingInterface<MFTNLayers> path (unchanged
// LegacyTrackerScratch<10>/DetectorTraversalBinding), and the combined
// workflow's new SurfaceTrackingScratch/SurfacePlanBinding instantiation.
// The two anonymous-namespace computeTracksMClabels<NLayers, ScratchT>
// instantiations these need are triggered implicitly by each Tracker<...>
// explicit instantiation below (via clustersToTracks()'s own call), not
// declared separately.
template class Tracker<7, LegacyTrackerScratch<7>, DetectorTraversalBinding>;
template class Tracker<10, LegacyTrackerScratch<10>, DetectorTraversalBinding>;
template class Tracker<10, SurfaceTrackingScratch, SurfacePlanBinding>;

} // namespace o2::itsmft::tracking

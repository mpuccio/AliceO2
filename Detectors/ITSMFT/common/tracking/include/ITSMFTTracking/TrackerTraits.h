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
/// \file TrackerTraits.h
/// \brief Shared CA tracker traits: same ITS-style tracklet/cell/road logic; MFT uses x-y LUT and forward refit
///

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKERTRAITS_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKERTRAITS_H_

#include <oneapi/tbb.h>

#include "ITSMFTTracking/CATrackTypes.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITStracking/BoundedAllocator.h"

namespace o2::itsmft::tracking
{

template <int NLayers>
class TrackerTraits
{
 public:
  using TimeFrameN = TimeFrame<NLayers>;
  using IndexTableUtilsN = o2::itsmft::IndexTableUtils<NLayers>;
  using CellSeedN = typename TimeFrameN::CellSeedN;
  using TrackSeedN = typename TimeFrameN::TrackSeedN;

  virtual ~TrackerTraits() = default;
  virtual void adoptTimeFrame(TimeFrameN* tf) { mTimeFrame = tf; }
  virtual void initialiseTimeFrame(const int iteration)
  {
    mTimeFrame->initialise(mTrkParams[iteration], mTrkParams[iteration].NLayers, iteration);
  }

  virtual void computeLayerTracklets(const int iteration, int iVertex);
  virtual void computeLayerCells(const int iteration);
  virtual void findCellsNeighbours(const int iteration);
  virtual void findRoads(const int iteration);

  template <typename InputSeed>
  void processNeighbours(int iteration, int defaultCellTopologyId, int iLevel, const bounded_vector<InputSeed>& currentCellSeed, const bounded_vector<int>& currentCellId, const bounded_vector<int>& currentCellTopologyId, bounded_vector<TrackSeedN>& updatedCellSeed, bounded_vector<int>& updatedCellId, bounded_vector<int>& updatedCellTopologyId);

  void acceptTracks(int iteration, bounded_vector<CATrackType<NLayers>>& tracks, bounded_vector<bounded_vector<int>>& firstClusters);
  void markTracks(int iteration);

  void updateTrackingParameters(const std::vector<TrackingParameters>& trkPars) { mTrkParams = trkPars; }
  TimeFrameN* getTimeFrame() { return mTimeFrame; }

  virtual void setBz(float bz);
  float getBz() const { return mBz; }
  virtual const char* getName() const noexcept { return "CPU"; }
  virtual bool isGPU() const noexcept { return false; }
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool) noexcept { mMemoryPool = pool; }
  auto getMemoryPool() const noexcept { return mMemoryPool; }

  void setNThreads(int n, std::shared_ptr<tbb::task_arena>& arena);
  int getNThreads() { return mTaskArena->max_concurrency(); }

  int getTFNumberOfClusters() const { return mTimeFrame->getNumberOfClusters(); }
  int getTFNumberOfTracklets() const { return mTimeFrame->getNumberOfTracklets(); }
  int getTFNumberOfCells() const { return mTimeFrame->getNumberOfCells(); }

 private:
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  std::shared_ptr<tbb::task_arena> mTaskArena;

 protected:
  TimeFrameN* mTimeFrame = nullptr;
  std::vector<TrackingParameters> mTrkParams;
  float mBz{-999.f};
};

using TrackerTraitsITS = TrackerTraits<constants::ITSNLayers>;
using TrackerTraitsMFT = TrackerTraits<constants::MFTNLayers>;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKERTRAITS_H_ */

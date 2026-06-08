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
/// \file CATracker.h
/// \brief Shared CA tracker orchestrator (same role as ITStracking/Tracker.h)
///

#ifndef ALICEO2_ITSMFT_TRACKING_CATRACKER_H_
#define ALICEO2_ITSMFT_TRACKING_CATRACKER_H_

#include <memory>
#include <vector>

#include <oneapi/tbb/task_arena.h>

#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"

namespace o2::itsmft::tracking
{

template <int NLayers>
class Tracker
{
 public:
  using TimeFrameN = TimeFrame<NLayers>;
  using TrackerTraitsN = TrackerTraits<NLayers>;

  explicit Tracker(TrackerTraitsN* traits);

  void adoptTimeFrame(TimeFrameN& tf);
  void setParameters(const std::vector<TrackingParameters>& p) { mTrkParams = p; }
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool) { mMemoryPool = pool; }
  void setBz(float bz) { mTraits->setBz(bz); }
  void setNThreads(int n, std::shared_ptr<tbb::task_arena>& arena) { mTraits->setNThreads(n, arena); }

  /// Run all configured iterations; returns elapsed ms or -1 on failure.
  float clustersToTracks();

  const TimeFrameN& getTimeFrame() const { return *mTimeFrame; }
  TimeFrameN& getTimeFrame() { return *mTimeFrame; }

 private:
  void initialiseTimeFrame(int iteration) { mTraits->initialiseTimeFrame(iteration); }
  void computeTracklets(int iteration, int iVertex) { mTraits->computeLayerTracklets(iteration, iVertex); }
  void computeCells(int iteration) { mTraits->computeLayerCells(iteration); }
  void findCellsNeighbours(int iteration) { mTraits->findCellsNeighbours(iteration); }
  void findRoads(int iteration) { mTraits->findRoads(iteration); }
  void rectifyClusterIndices();
  void sortTracks();

  TrackerTraitsN* mTraits = nullptr;
  TimeFrameN* mTimeFrame = nullptr;
  std::vector<TrackingParameters> mTrkParams;
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
};

template <int NLayers>
using CATracker = Tracker<NLayers>;

using TrackerITS = Tracker<constants::ITSNLayers>;
using TrackerMFT = Tracker<constants::MFTNLayers>;
using CATrackerITS = TrackerITS;
using CATrackerMFT = TrackerMFT;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_CATRACKER_H_ */

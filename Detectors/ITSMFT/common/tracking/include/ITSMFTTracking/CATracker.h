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

/// Gate 3 common-CA compatibility sentinel: the single float value
/// clustersToTracks() may return to mean "this TimeFrame was a recoverable
/// per-TF failure, DropTFUponFailure was set, and the TimeFrame has already
/// been fully wiped -- do not publish anything for it, and it is safe to
/// continue with the next TimeFrame." isDroppedTimeFrame() tests this exact
/// value, never a sign check, so no other negative result, NaN, or infinity
/// can be mistaken for a drop. Every other failure (structural or
/// unclassified) throws instead of returning a sentinel.
///
/// This is a bounded compatibility slice: a typed tracking outcome (success /
/// dropped / structural-failure, with a reason) should replace this float
/// sentinel before the common-CA failure contract is considered final.
inline constexpr float kDroppedTimeFrameResult = -1.f;

/// Exact-match test for the drop sentinel above. Deliberately not `result <
/// 0.f`: only the literal kDroppedTimeFrameResult value means "dropped",
/// so callers cannot silently widen the contract to other negative values.
inline bool isDroppedTimeFrame(float result)
{
  return result == kDroppedTimeFrameResult;
}

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

  /// Run all configured iterations. Returns elapsed ms on success, or the
  /// exact kDroppedTimeFrameResult sentinel when a recoverable per-TF
  /// failure was dropped (DropTFUponFailure=true); the TimeFrame is always
  /// fully wiped before that return. Any structural or unclassified failure,
  /// and any recoverable failure with DropTFUponFailure=false, throws
  /// instead of returning -- the TimeFrame is fully wiped before the
  /// exception propagates.
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

using TrackerITS = Tracker<ITSNLayers>;
using TrackerMFT = Tracker<o2::mft::constants::mft::LayersNumber>;
using CATrackerITS = TrackerITS;
using CATrackerMFT = TrackerMFT;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_CATRACKER_H_ */

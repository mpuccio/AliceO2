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

/// \file NeighbourStage.h
/// \brief The neighbour stage's per-triplet function.

#ifndef GTRK_TRACKING_NEIGHBOURSTAGE_H_
#define GTRK_TRACKING_NEIGHBOURSTAGE_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/Neighbour.h"
#include "GenericTracking/Tracking/NeighbourCut.h"
#include "GenericTracking/Tracking/OrderedRoadCounts.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/TrackletPaths.h"
#include "GenericTracking/Tracking/Triplet.h"
#include "GenericTracking/Tracking/TripletEndIndex.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include <cstddef>

namespace o2::gtrk
{

/// Emit every predecessor of one triplet `T2 = (b -> c, c -> d)`: the
/// triplets `T1 = (a -> b, b -> c)` ending with `T2.first` that pass the
/// neighbour cut, in the end index's order, which is the arena's.
///
/// The stage is a pure function of its inputs: a backend may run an element
/// twice, once to count and once to emit. The measurements `b`, `c` and `d`
/// are read once per element; each candidate adds `a`.
/// With \a readyLevels (always, in the trackers) generation follows the
/// layers' order: predecessor levels are final, only the highest accepted
/// stratum is emitted, and this element's level is written; without (the
/// tests' all-edges graph) every predecessor that passes is. \a roadCounts
/// needs \a readyLevels; its flags are initialized once before generation.
/// A repeated element gives the same result: counts are overwritten and marks
/// are idempotent.
template <DetectorLayout Det>
class NeighbourStage
{
 public:
  using Index = typename StageOutput<NeighbourT>::Index;

  // The delegating constructor retains this output for later level writes.
  // NOLINTBEGIN(readability-non-const-parameter)
  NeighbourStage(const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const StageOutput<TrackletT>& tracklets,
                 const StageOutput<TripletT>& triplets, const TripletEndIndex& ending, TrackletPathBounds paths, std::uint32_t* readyLevels = nullptr,
                 OrderedRoadCounts roadCounts = {})
    : NeighbourStage(&parameters, &frame, tracklets.view(), triplets.view(), ending.offsets(), ending.triplets(), paths, readyLevels, roadCounts)
  {
  }
  // NOLINTEND(readability-non-const-parameter)

  /// The same stage over inputs wherever they are -- a device's copies:
  /// every argument is only stored. The end
  /// index is `trackletCount + 1` offsets and the triplet indices they delimit.
  NeighbourStage(const IterationParameters<Det>* parameters, const TimeFrame<Det>* frame, TrackletConstView tracklets, TripletConstView triplets,
                 const Index* endOffsets, const Index* endTriplets, TrackletPathBounds paths, std::uint32_t* readyLevels = nullptr,
                 OrderedRoadCounts roadCounts = {})
    : mParameters(parameters),
      mFrame(frame),
      mTracklets(tracklets),
      mTriplets(triplets),
      mEndOffsets(endOffsets),
      mEndTriplets(endTriplets),
      mPaths(paths),
      mReadyLevels(readyLevels),
      mRoadCounts(roadCounts)
  {
  }

  template <typename Emitter>
  GPUdi() void operator()(Index element, Emitter& emitter) const
  {
    if (mReadyLevels != nullptr) {
      mReadyLevels[element] = 1U;
    }
    if (mRoadCounts.roads != nullptr) {
      mRoadCounts.base(element, mTracklets, mTriplets);
    }
    const Index shared = mTriplets.first[element];
    const Index begin = mEndOffsets[shared];
    const Index end = mEndOffsets[shared + 1U];
    countWork(WorkCounter::kNeighbourTriplets);
    countWork(WorkCounter::kNeighbourCandidates, end - begin);
    if (begin == end) {
      return;
    }
    const Index last = mTriplets.second[element];
    // Only triplets on a triplet path of minClusters, and only their
    // neighbours on one too.
    if (!mPaths.onPath(shared, last)) {
      return;
    }
    const std::size_t secondLink = mTracklets.link[shared];
    const LayerID secondLayer = mParameters->linkSources[secondLink];
    const LayerID thirdLayer = mParameters->linkTargets[secondLink];
    const Scalar curvature = mTriplets.curvature[element];
    const Scalar varCurvature = mTriplets.varCurvature[element];
    // A kinked triplet's middle point takes the kink allowance in the adjacent
    // fit; two kinked triplets never neighbour.
    const bool kinked = mTriplets.kinked[element] != 0U;
    const Scalar thirdScale = kinked ? mParameters->graphKinkScale2 : 1.F;
    const auto columns = mFrame->measurements();
    const auto b = timedPointAt(columns, mTracklets.source[shared]);
    const auto c = timedPointAt(columns, mTracklets.target[shared]);
    const auto d = timedPointAt(columns, mTracklets.target[last]);

    // This triplet's kink terms, made at its first candidate to need them.
    KinkTerms terms;
    bool termsMade = false;
    bool termsValid = false;

    std::uint32_t highest = 0U;
    if (mReadyLevels != nullptr) {
      for (Index candidate = begin; candidate < end; ++candidate) {
        const auto level = mReadyLevels[mEndTriplets[candidate]];
        highest = level > highest ? level : highest;
      }
    }
    for (auto stratum = highest;; --stratum) {
      bool accepted = false;
      std::uint64_t roads = 0U;
      for (Index candidate = begin; candidate < end; ++candidate) {
        const Index from = mEndTriplets[candidate];
        if (mReadyLevels != nullptr && mReadyLevels[from] != stratum) {
          continue;
        }
        if (std::uint32_t{mPaths.prefix[mTriplets.first[from]]} + 1U + mPaths.suffix[last] < mPaths.minClusters) {
          continue;
        }
        const bool fromKinked = mTriplets.kinked[from] != 0U;
        if (kinked && fromKinked) {
          countWork(WorkCounter::kNeighbourBothKinked);
          continue;
        }
        if (!neighbourCurvatureCompatible(*mParameters, mTriplets.curvature[from], mTriplets.varCurvature[from], curvature, varCurvature)) {
          countWork(WorkCounter::kNeighbourCurvature);
          continue;
        }
        const Index first = mTriplets.first[from];
        const auto a = timedPointAt(columns, mTracklets.source[first]);
        TripletFitResult fit;
        if (!neighbourCompatible(*mParameters, secondLayer, thirdLayer, a, b, c, d, terms, termsMade, termsValid, fit,
                                 fromKinked ? mParameters->graphKinkScale2 : 1.F, thirdScale)) {
          continue;
        }
        emitter.emit(Neighbour{.from = from});
        if (mRoadCounts.roads != nullptr) {
          mRoadCounts.add(from, roads);
        }
        accepted = true;
      }
      if (mReadyLevels == nullptr) {
        break;
      }
      if (accepted) {
        mReadyLevels[element] = stratum + 1U;
        if (mRoadCounts.roads != nullptr) {
          mRoadCounts.finish(element, roads);
        }
        break;
      }
      if (stratum <= 1U) {
        break;
      }
    }
  }

 private:
  const IterationParameters<Det>* mParameters = nullptr;
  const TimeFrame<Det>* mFrame = nullptr;
  TrackletConstView mTracklets;
  TripletConstView mTriplets;
  const Index* mEndOffsets = nullptr;
  const Index* mEndTriplets = nullptr;
  TrackletPathBounds mPaths;
  std::uint32_t* mReadyLevels = nullptr; ///< only with topologically ordered generation
  OrderedRoadCounts mRoadCounts;         ///< optional, only with readyLevels
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_NEIGHBOURSTAGE_H_

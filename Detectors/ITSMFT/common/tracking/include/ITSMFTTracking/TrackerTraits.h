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

#include <array>
#include <optional>
#include <stdexcept>
#include <string>

#include <oneapi/tbb.h>

#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TransitionPolicyBinding.h"
#include "ITSMFTTracking/TransitionPolicyDispatch.h"
#include "ITSMFTTracking/TransitionPolicyOperations.h"
#include "ITSMFTTracking/TransitionPolicyState.h"
#include "ITStracking/BoundedAllocator.h"

namespace o2::itsmft::tracking
{

enum class TraversalFailureReason : uint8_t {
  MissingLayout,
  StaleLayout,
  IterationOutOfRange,
  LegacyIndexMismatch,
  InvalidTraversalSchedule,
  MixedPolicyLayout,
  StateFamilyMismatch,
  InvalidPolicyParameters
};

class TraversalException final : public std::runtime_error
{
 public:
  TraversalException(int iteration, TraversalFailureReason reason)
    : std::runtime_error{"CA traversal initialization failed at iteration " + std::to_string(iteration) + " (reason=" + std::to_string(static_cast<int>(reason)) + ")"},
      mIteration{iteration}, mReason{reason}
  {
  }

  int getIteration() const noexcept { return mIteration; }
  TraversalFailureReason getReason() const noexcept { return mReason; }

 private:
  int mIteration{-1};
  TraversalFailureReason mReason{TraversalFailureReason::MissingLayout};
};

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
  virtual void initialiseTimeFrame(const int iteration);

  virtual void computeLayerTracklets(const int iteration, int iVertex);
  virtual void computeLayerCells(const int iteration);
  virtual void findCellsNeighbours(const int iteration);
  virtual void findRoads(const int iteration);

  template <TransitionPolicyTag Tag, typename InputSeed>
  void processNeighbours(int iteration, int defaultCellTopologyId, int iLevel, const bounded_vector<InputSeed>& currentCellSeed, const bounded_vector<int>& currentCellId, const bounded_vector<int>& currentCellTopologyId, bounded_vector<TrackSeedN>& updatedCellSeed, bounded_vector<int>& updatedCellId, bounded_vector<int>& updatedCellTopologyId, const typename TransitionPolicyTraits<Tag>::Params& params);

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

  int getTraversalGroupingCount() const noexcept { return mTraversalGroupingCount; }
  int getPolicyBindingCount(TransitionPolicyTag tag) const noexcept;
  bool hasTraversalCache() const noexcept { return mTraversalGrouping.has_value(); }

 private:
  void resetTraversalCache() noexcept;
  void validateLegacyParity(int iteration, const DetectorLayoutView& layout, TransitionPolicyTag& activeTag, bool& mixedPolicy) const;

  template <TransitionPolicyTag Tag>
  void computeLayerTrackletsForPolicy(int iteration,
                                      int iVertex,
                                      gsl::span<const TransitionId> transitionIds,
                                      const typename TransitionPolicyTraits<Tag>::Params& params);

  template <TransitionPolicyTag Tag>
  void computeLayerCellsForPolicy(int iteration,
                                  const typename TimeFrameN::TrackingTopologyN::View& topology,
                                  const typename TransitionPolicyTraits<Tag>::Params& params);

  template <TransitionPolicyTag Tag>
  void findCellsNeighboursForPolicy(int iteration,
                                    gsl::span<const CellTopologyId> scheduledCells,
                                    const typename TransitionPolicyTraits<Tag>::Params& params);

  template <TransitionPolicyTag Tag>
  void findRoadsForPolicy(int iteration, const typename TransitionPolicyTraits<Tag>::Params& params);

  // Gate 3 transition-preparation slice: relocated from TimeFrame::initialise()
  // (Architecture.md Sec 10/10.1). Called from initialiseTimeFrame() only
  // after all existing and new fallible validation for this iteration has
  // already succeeded (activeTag, cylinder/disk params, attachHitConfig,
  // geometryConfig, and -- for DiskDisk -- referenceCoordinateView); this
  // method itself never throws. Fills TimeFrame's already-sized
  // mTransitionMSAngles/mTransitionPhiCuts (container sizing stays in
  // TimeFrame::initialise(), unchanged) by iterating legacy transitionIds
  // 0..nTransitions-1 in increasing order directly off
  // TimeFrame::getTrackingTopologyView() -- not through mTraversalGrouping's
  // per-tag span -- so the loop-carried oneOverR ratchet threads in exactly
  // the same order the frozen legacy code used, with no dependency on
  // grouping-span ordering.
  template <TransitionPolicyTag Tag>
  void prepareTransitionScatteringAndBendingForPolicy(int iteration,
                                                      const LayerGeometryConfigView& geometryConfig,
                                                      const DiskDiskReferenceCoordinateView& referenceCoordinateView);

  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  std::shared_ptr<tbb::task_arena> mTaskArena;
  DetectorLayoutView mTraversalLayout{};
  std::optional<TransitionPolicyGrouping> mTraversalGrouping;
  std::optional<CylinderCylinderPolicyParams> mCylinderPolicyParams;
  std::optional<DiskDiskPolicyParams> mDiskPolicyParams;
  AttachHitPolicyConfigView mAttachHitConfig;
  int mTraversalGroupingCount{0};
  std::array<int, 2> mPolicyBindingCounts{};

 protected:
  TimeFrameN* mTimeFrame = nullptr;
  std::vector<TrackingParameters> mTrkParams;
  float mBz{-999.f};
};

using TrackerTraitsITS = TrackerTraits<ITSNLayers>;
using TrackerTraitsMFT = TrackerTraits<o2::mft::constants::mft::LayersNumber>;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKERTRAITS_H_ */

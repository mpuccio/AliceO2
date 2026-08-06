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
/// \file Tracker.h
/// \brief Shared runtime-plan tracker orchestrator.
///

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKER_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKER_H_

#include <cstdint>
#include <memory>
#include <vector>

#include <gsl/span>

#include <oneapi/tbb/task_arena.h>

#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/SurfaceGraph.h"
#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"

namespace o2::itsmft::tracking
{

/// Gate 4 C2 Slice 2: the typed outcome Tracker::clustersToTracks()
/// returns in place of the old float+sentinel compatibility contract
/// (kDroppedTimeFrameResult/isDroppedTimeFrame below, retained only because
/// ITSMFTTrackingInterface/CATrackerSpec.cxx still externally consume that
/// exact float contract -- see TrackingInterface.cxx::runTracking(), the
/// sole place that now translates between the two).
///
/// clustersToTracks() itself only ever *returns* Success or
/// RecoverableDropped: a recoverable, per-TF resource failure
/// (BoundedMemoryResource::MemoryLimitExceeded, std::bad_alloc) with
/// DropTFUponFailure=true. Every other failure -- structural/configuration
/// (TraversalException, any reason, including Gate 4 C2 Slice 1's
/// TraversalBindingMismatch), unclassified (any other std::exception), or a
/// recoverable failure with DropTFUponFailure=false -- retains its existing,
/// already-tested contract of propagating as a thrown C++ exception past
/// clustersToTracks()'s own boundary; this is deliberate ("retain exceptions
/// where that is the established contract"), not an oversight: reusing the
/// existing exception-based classification means a mismatched/invalid
/// binding can never be silently reclassified as a dropped, recoverable
/// result merely because DropTFUponFailure happens to be true. Structural is
/// part of this enum's vocabulary -- a future combined coordinator that
/// wraps clustersToTracks() in its own try/catch is expected to construct it
/// from the caught exception -- but no code in this slice constructs it via
/// a normal return.
enum class TrackingOutcome : uint8_t {
  Success,
  RecoverableDropped,
  Structural
};

/// clustersToTracks()'s complete return value on every path it does not
/// throw past its own boundary. `elapsedMs` is only meaningful when
/// `outcome == Success` (0.f otherwise, matching the old sentinel's implicit
/// contract of never being read on a drop).
struct TrackingResult {
  TrackingOutcome outcome{TrackingOutcome::Success};
  float elapsedMs{0.f};
};

struct TrackerIterationConfiguration {
  std::vector<SurfaceGraphSubgraph> graphSubgraphs;
  std::vector<SurfacePlanBinding::Declaration> bindings;
  std::vector<TrackingParameters> parameters;
};

struct TrackerInitialization {
  SurfaceCatalogView catalog;
  std::vector<TrackerIterationConfiguration> iterations;
  std::shared_ptr<BoundedMemoryResource> memoryPool;
};

enum class TrackerInitializationError : uint8_t {
  None,
  EmptyConfiguration,
  MissingCatalog,
  MissingMemoryPool,
  GraphBuildFailed,
  BindingCountMismatch,
  BindingBuildFailed,
  DuplicateSource,
  CapacityMismatch
};

struct TrackerInitializationResult {
  TrackerInitializationError error{TrackerInitializationError::None};
  std::size_t failedIteration{static_cast<std::size_t>(-1)};
  SurfaceGraphBuildError graphError{SurfaceGraphBuildError::None};
  SurfacePlanBindingError bindingError{SurfacePlanBindingError::None};
  bool ok() const noexcept { return error == TrackerInitializationError::None; }
};

/// Gate 3 common-CA compatibility sentinel, retained for
/// ITSMFTTrackingInterface's own external float contract (the ITS/MFT
/// workflow adapters still call isDroppedTimeFrame() on
/// processTimeFrame()'s return value) -- the single float value that means
/// "this TimeFrame was a recoverable per-TF failure, DropTFUponFailure was
/// set, and the TimeFrame has already been fully wiped -- do not publish
/// anything for it, and it is safe to continue with the next TimeFrame."
/// isDroppedTimeFrame() tests this exact value, never a sign check, so no
/// other negative result, NaN, or infinity can be mistaken for a drop.
inline constexpr float kDroppedTimeFrameResult = -1.f;

/// Exact-match test for the drop sentinel above. Deliberately not `result <
/// 0.f`: only the literal kDroppedTimeFrameResult value means "dropped",
/// so callers cannot silently widen the contract to other negative values.
inline bool isDroppedTimeFrame(float result) noexcept
{
  return result == kDroppedTimeFrameResult;
}

class Tracker
{
 public:
  explicit Tracker(TrackerTraits* traits);

  TrackerInitializationResult initialize(TimeFrame& frame, const TrackerInitialization& configuration);

  // Binds this tracker's two collaborators, each an independent bind-once
  // pointer -- neither owns nor stores a reference to the other.
  // `scratch`/`frame` must each outlive every subsequent clustersToTracks()
  // call.
  void adoptScratch(SurfaceTrackingScratch& scratch);
  void adoptFrame(TimeFrame& frame);
  void setSource(ClusterSourceId source) noexcept { mSource = source; }
  void setBz(float bz) { mTraits->setBz(bz); }
  void setNThreads(int n, std::shared_ptr<tbb::task_arena>& arena) { mTraits->setNThreads(n, arena); }

  /// Run all configured iterations. Returns {Success, elapsed ms} on
  /// success, or {RecoverableDropped, 0.f} when a recoverable per-TF failure
  /// was dropped (DropTFUponFailure=true); the event is always fully reset
  /// before that return.
  /// Any structural or unclassified failure, and any recoverable failure
  /// with DropTFUponFailure=false, throws instead of returning -- see
  /// TrackingOutcome's own doc for why this is deliberate -- the event is
  /// always fully reset before the exception propagates. This tracker never
  /// decides *which* reset a combined future owner with several
  /// participating scratches would want (see resetTimeFrameEvent()'s own
  /// doc) -- for this single-detector bridge, every recoverable failure here
  /// resets both its own scratch and the shared TimeFrame.
  TrackingResult clustersToTracks(TrackingOperationAdapter& operationAdapter);

  const SurfaceTrackingScratch& getScratch() const { return *mScratch; }
  SurfaceTrackingScratch& getScratch() { return *mScratch; }

 private:
  void initialiseTimeFrame(int iteration)
  {
    mTraits->adoptSurfacePlanBinding(mFrame->getBinding(iteration, mSource));
    mTraits->initialiseTimeFrame(iteration, mFrame->getGraphs());
  }
  void computeTracklets(int iteration, int iVertex) { mTraits->computeLayerTracklets(iteration, iVertex); }
  void computeCells(int iteration) { mTraits->computeLayerCells(iteration); }
  void findCellsNeighbours(int iteration) { mTraits->findCellsNeighbours(iteration); }
  void findRoads(int iteration, TrackingOperationAdapter& operationAdapter) { mTraits->findRoads(iteration, operationAdapter); }
  TrackerTraits* mTraits = nullptr;
  SurfaceTrackingScratch* mScratch = nullptr;
  TimeFrame* mFrame = nullptr;
  ClusterSourceId mSource{};
};
} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKER_H_ */

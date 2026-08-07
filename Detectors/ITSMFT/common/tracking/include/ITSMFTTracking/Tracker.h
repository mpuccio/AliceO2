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
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"

namespace o2::itsmft::tracking
{

/// Tracker::run() returns a typed outcome while the standalone interface keeps
/// its existing float publication contract at the adapter edge.
///
/// run() itself only ever *returns* Success or
/// RecoverableDropped: a recoverable, per-TF resource failure
/// (BoundedMemoryResource::MemoryLimitExceeded, std::bad_alloc) with
/// DropTFUponFailure=true. Every other failure -- structural/configuration
/// (TraversalException, any reason, including Gate 4 C2 Slice 1's
/// TraversalBindingMismatch), unclassified (any other std::exception), or a
/// recoverable failure with DropTFUponFailure=false -- retains its existing,
/// already-tested contract of propagating as a thrown C++ exception past
/// run()'s own boundary; this is deliberate ("retain exceptions
/// where that is the established contract"), not an oversight: reusing the
/// existing exception-based classification means a mismatched/invalid
/// binding can never be silently reclassified as a dropped, recoverable
/// result merely because DropTFUponFailure happens to be true. Structural is
/// part of this enum's vocabulary; workflow/adapter callers convert propagated
/// failures at their own boundary.
enum class TrackingOutcome : uint8_t {
  Success,
  RecoverableDropped,
  Structural
};

/// run()'s complete return value on every path it does not
/// throw past its own boundary. `elapsedMs` is only meaningful when
/// `outcome == Success` (0.f otherwise, matching the old sentinel's implicit
/// contract of never being read on a drop).
struct TrackingResult {
  TrackingOutcome outcome{TrackingOutcome::Success};
  float elapsedMs{0.f};
  // Cumulative accepted-result boundaries let application publication adapters
  // reproduce per-iteration staging without crossing the refit seam.
  std::vector<std::size_t> acceptedTrackCounts;
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

class Tracker
{
 public:
  explicit Tracker(TrackingOperationAdapter* operationAdapter = nullptr, ClusterSourceId source = {})
    : mOperationAdapter(operationAdapter), mSource(source)
  {
  }

  TrackerInitializationResult initialize(TimeFrame& frame, const TrackerInitialization& configuration);

  // The workflow/application owns this adapter and keeps it alive for every
  // run. Tracker stores no frame, graph, workspace, or event state.
  void setOperationAdapter(TrackingOperationAdapter* operationAdapter) noexcept { mOperationAdapter = operationAdapter; }
  void setSource(ClusterSourceId source) noexcept { mSource = source; }

  /// Run all configured iterations. Returns {Success, elapsed ms} on
  /// success, or {RecoverableDropped, 0.f} when a recoverable per-TF failure
  /// was dropped (DropTFUponFailure=true); the event is always fully reset
  /// before that return.
  /// Any structural or unclassified failure, and any recoverable failure
  /// with DropTFUponFailure=false, throws instead of returning -- see
  /// TrackingOutcome's own doc for why this is deliberate -- the event is
  /// always fully reset before the exception propagates.
  TrackingResult run(TimeFrame& frame, TrackerTraits& traits);

 private:
  TrackingOperationAdapter* mOperationAdapter = nullptr;
  ClusterSourceId mSource{};
};
} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKER_H_ */

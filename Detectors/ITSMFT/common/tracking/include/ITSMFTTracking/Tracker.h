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

namespace o2::itsmft::tracking
{

struct TrackerTestAccess;

/// `run()` returns `Success`, or `RecoverableDropped` for a recoverable
/// per-TimeFrame resource failure (`MemoryLimitExceeded` or `std::bad_alloc`)
/// when `DropTFUponFailure` is enabled.
/// Structural and unclassified failures, and recoverable failures with
/// dropping disabled, propagate as exceptions.
enum class TrackingOutcome : uint8_t {
  Success,
  RecoverableDropped,
  Structural
};

/// Complete return value for paths that do not throw. `elapsedMs` is meaningful
/// only when `outcome == Success`; it is 0.f otherwise.
struct TrackingResult {
  TrackingOutcome outcome{TrackingOutcome::Success};
  float elapsedMs{0.f};
  // Accepted-result counts are indexed by configured iteration.
  std::vector<std::size_t> acceptedTrackCounts;
};

struct TrackerIterationConfiguration {
  SurfaceGraphDefinition graph;
  TrackingParameters parameters;
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
  TraversalPlanBuildFailed,
  DuplicateSource,
  CapacityMismatch
};

struct TrackerInitializationResult {
  TrackerInitializationError error{TrackerInitializationError::None};
  std::size_t failedIteration{static_cast<std::size_t>(-1)};
  SurfaceGraphBuildError graphError{SurfaceGraphBuildError::None};
  bool ok() const noexcept { return error == TrackerInitializationError::None; }
};

class Tracker
{
 public:
  explicit Tracker(SeedRefitFunction refitFunction = nullptr)
    : mRefitFunction(refitFunction)
  {
  }

  TrackerInitializationResult initialize(TimeFrame& frame, const TrackerInitialization& configuration);

  // Tracker stores no frame, configuration, workspace, graph, or event state.
  void setSeedRefitFunction(SeedRefitFunction refitFunction) noexcept { mRefitFunction = refitFunction; }

  /// Run all configured iterations. Returns `Success` on success or
  /// `RecoverableDropped` when an allowed recoverable per-TF failure is
  /// dropped. The event is reset before a dropped return or propagated error.
  TrackingResult run(TimeFrame& frame, TrackerTraits& traits);

 private:
  friend struct TrackerTestAccess;
  void initializeTraversalWorkspace(TraversalWorkspaceView& view) const;
  SeedRefitFunction mRefitFunction = nullptr;
};
} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKER_H_ */

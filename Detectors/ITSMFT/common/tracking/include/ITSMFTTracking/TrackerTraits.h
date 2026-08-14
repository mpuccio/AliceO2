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
#include <vector>

#include <gsl/span>
#include <oneapi/tbb.h>

#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/GenericTrack.h"
#include "ITSMFTTracking/SurfaceGraph.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/detail/TrackingKernelParameters.h"
#include "ITStracking/BoundedAllocator.h"

namespace o2::itsmft::tracking
{

#ifndef GPUCA_GPUCODE

// Call-scoped native refit supplied by the detector/workflow edge.
// Publication and reset remain outside the tracking transaction.
using SeedRefitFunction = bool (*)(const TrackSeed& seed,
                                   const TrackingParameters& params,
                                   float bz,
                                   SurfaceTrackingScratch& scratch,
                                   gsl::span<const gsl::span<const GlobalMeasurement>> layerGlobals,
                                   gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                                   SurfaceCatalogView surfaceCatalog,
                                   gsl::span<const SurfaceId> orderedSurfaces,
                                   TrackingCandidate& candidate);

#endif

struct DiskReferenceCoordinateView;
struct TrackerTestAccess;

struct TraversalWorkspaceView {
  int iteration{-1};
  TimeFrame& frame;
  SurfaceTrackingScratch& scratch;
  SurfaceGraphView graph{};
  gsl::span<const TrackingParameters> parameters;
  float bz{0.f};
  TraversalWorkspace& workspace;
};

enum class TraversalFailureReason : uint8_t {
  MissingLayout,
  StaleLayout,
  IterationOutOfRange,
  SparseTopologyMismatch,
  InvalidTraversalSchedule,
  MixedSurfaceKindLayout,
  SurfaceKindMismatch,
  InvalidSurfaceParameters,
  // The iteration's index-table configuration is structurally invalid.
  InvalidIndexTableConfiguration,
  // A non-FirstPass configuration disagrees with the TimeFrame's configuration or LUT.
  IndexTableConfigurationMismatch,
  // LayerxX0 disagrees with mapped authoritative material, or the mapping is invalid.
  // Raised before tracking state is touched; the descriptor is never overwritten.
  LegacyMaterialMismatch,
  // The active SurfaceKind does not support the configured MatCorrType.
  // This structural error is reset and rethrown regardless of drop policy.
  // An unrecognized CorrType is reported separately by AttachHitConfigView::isValid().
  UnsupportedMaterialCorrectionMode,
  // Per-position normalized measurements disagree with the loaded frame or
  // compatibility data. Raised before tracking state is touched; spans commit only on success.
  NormalizedMeasurementMismatch,
  // orderedSurfaces does not map plan positions bijectively to global SurfaceId.
  SurfaceLayerMappingMismatch,
  // The adopted binding cannot translate a traversal ID to a compact scratch slot.
  // This is a binding/layout mismatch, detected before scratch access.
  TraversalBindingMismatch
};

class TraversalException final : public std::runtime_error
{
 public:
  TraversalException(int iteration, TraversalFailureReason reason)
    : std::runtime_error{"CA traversal initialization failed at iteration " + std::to_string(iteration) + " (reason=" + std::to_string(static_cast<int>(reason)) + ")"},
      mIteration{iteration},
      mReason{reason}
  {
  }

  int getIteration() const noexcept { return mIteration; }
  TraversalFailureReason getReason() const noexcept { return mReason; }

 private:
  int mIteration{-1};
  TraversalFailureReason mReason{TraversalFailureReason::MissingLayout};
};

// Backend implementation of a traversal supplied explicitly by Tracker.
class TrackerTraits
{
 public:
  virtual ~TrackerTraits() = default;
  // The production caller supplies all event and iteration state explicitly.
  void runTraversal(TraversalWorkspaceView view, SeedRefitFunction refitFunction);

  virtual const char* getName() const noexcept { return "CPU"; }
  virtual bool isGPU() const noexcept { return false; }
  void setNThreads(int n, std::shared_ptr<tbb::task_arena>& arena);
  int getNThreads() { return mTaskArena->max_concurrency(); }

 private:
  friend struct TrackerTestAccess;

  void acceptTracks(TraversalWorkspaceView& context, int iteration,
                    bounded_vector<TrackingCandidate>& tracks,
                    bounded_vector<bounded_vector<int>>& firstClusters);

  // Tracklet and cell enumeration are common; coordinate selection is owned
  // by their operation leaves.
  void computeLayerTracklets(TraversalWorkspaceView& context, int iteration, int iVertex);
  void computeLayerTrackletsImpl(TraversalWorkspaceView& context, int iteration, int iVertex,
                                 gsl::span<const EdgeId> edgeIds);
  void computeLayerCells(TraversalWorkspaceView& context, int iteration);
  void computeLayerCellsImpl(TraversalWorkspaceView& context, int iteration,
                             gsl::span<const CellTopologyId> cellIds);
  void findCellsNeighbours(TraversalWorkspaceView& context, int iteration);
  void findCellsNeighboursForSchedule(TraversalWorkspaceView& context, int iteration,
                                      gsl::span<const CellTopologyId> scheduledCells,
                                      const TrackingKernelParameters& params);

  void findRoads(TraversalWorkspaceView& context, int iteration, SeedRefitFunction refitFunction);
  void findRoadsImpl(TraversalWorkspaceView& context, int iteration, SeedRefitFunction refitFunction);

  // Neighbour processing helper; it does not encode a detector layer count.
  template <typename InputSeed>
  void processNeighbours(TraversalWorkspaceView& context, int iteration, int defaultCellTopologyId, int iLevel, const bounded_vector<InputSeed>& currentCellSeed, const bounded_vector<int>& currentCellId, const bounded_vector<int>& currentCellTopologyId, bounded_vector<TrackSeed>& updatedCellSeed, bounded_vector<int>& updatedCellId, bounded_vector<int>& updatedCellTopologyId, const TrackingKernelParameters& params);

  std::shared_ptr<tbb::task_arena> mTaskArena;
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKERTRAITS_H_ */

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
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
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

struct LayerGeometryConfigView;
struct DiskReferenceCoordinateView;

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

// The traits borrow the TimeFrame workspace, plan binding, parameters, and
// refit function. Plan-sized bounds come from the adopted plan and scratch.
class TrackerTraits
{
 public:
  virtual ~TrackerTraits() = default;
  // Borrowed call-scoped pointers; neither object owns the other.
  virtual void adoptScratch(SurfaceTrackingScratch* scratch) { mScratch = scratch; }
  virtual void adoptFrame(TimeFrame* frame) { mFrame = frame; }
  // `binding` is borrowed, never owned or copied, and must match the graph iteration.
  void adoptSurfacePlanBinding(const SurfacePlanBinding* binding) noexcept { mBinding = binding; }
  // `graphs` is the caller-owned immutable graph vector for this call.
  virtual void initialiseTimeFrame(const int iteration, const std::vector<SurfaceGraph>& graphs);

  virtual void computeLayerTracklets(const int iteration, int iVertex);
  virtual void computeLayerCells(const int iteration);
  virtual void findCellsNeighbours(const int iteration);
  virtual void findRoads(const int iteration, SeedRefitFunction refitFunction);

  void acceptTracks(int iteration,
                    bounded_vector<TrackingCandidate>& tracks,
                    bounded_vector<bounded_vector<int>>& firstClusters);

  // Keeps accepted GenericTrack/TrackSeed pairs until the workflow consumes the result.
  bounded_vector<TrackingCandidate>& acceptedTracksForSharedStatus();
  void clearAcceptedTracksForSharedStatus();

  void updateTrackingParameters(gsl::span<const TrackingParameters> trkPars) { mTrkParams = trkPars; }
  SurfaceTrackingScratch* getScratch() { return mScratch; }

  virtual void setBz(float bz);
  float getBz() const { return mBz; }
  virtual const char* getName() const noexcept { return "CPU"; }
  virtual bool isGPU() const noexcept { return false; }
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool) noexcept { mMemoryPool = pool; }
  auto getMemoryPool() const noexcept { return mMemoryPool; }

  void setNThreads(int n, std::shared_ptr<tbb::task_arena>& arena);
  int getNThreads() { return mTaskArena->max_concurrency(); }

  int getTFNumberOfClusters() const { return mScratch->getNumberOfClusters(); }
  int getTFNumberOfTracklets() const { return mScratch->getNumberOfTracklets(); }
  int getTFNumberOfCells() const { return mScratch->getNumberOfCells(); }

  bool hasTraversalCache() const noexcept { return mTraversalCacheValid; }
  // Nominal material resolved from SurfaceDescriptor::material and orderedSurfaces
  // by the last successful initialiseTimeFrame(). It is committed only on success
  // and reset otherwise. Read-only inspection; production uses mAttachHitConfig.
  gsl::span<const NominalSurfaceMaterial> getLayerMaterial() const noexcept { return {mLayerMaterial.data(), mLayerMaterial.size()}; }
  // Normalized measurements resolved once by the last successful
  // initialiseTimeFrame(), using orderedSurfaces. Spans are non-owning, reset
  // on failure, and valid only while the TimeFrame's normalized frame is alive.
  // Read-only inspection; production uses mLayerMeasurements.
  gsl::span<const gsl::span<const SurfaceMeasurement>> getLayerMeasurements() const noexcept { return {mLayerMeasurements.data(), mLayerMeasurements.size()}; }
  gsl::span<const gsl::span<const GlobalMeasurement>> getLayerGlobalMeasurements() const noexcept { return {mLayerGlobalMeasurements.data(), mLayerGlobalMeasurements.size()}; }

 private:
  void resetTraversalCache() noexcept;
  void validateSparsePlan(int iteration, const SurfaceGraphView& graph) const;
  int requireSurfacePosition(int iteration, SurfaceId id) const;

  // Global-ID to compact-scratch-slot translation used by this class.
  int requireScratchTransitionSlot(int iteration, TransitionId id) const;
  int requireScratchCellSlot(int iteration, CellTopologyId id) const;

  // Tracklet and cell enumeration are common; coordinate selection is owned
  // by their operation leaves.
  void computeLayerTrackletsImpl(int iteration, int iVertex,
                                 gsl::span<const TransitionId> transitionIds);
  void computeLayerCellsImpl(int iteration,
                             gsl::span<const CellTopologyId> cellIds);
  void findRoadsForGraph(int iteration, SeedRefitFunction refitFunction);

  void findCellsNeighboursForSchedule(int iteration,
                                      gsl::span<const CellTopologyId> scheduledCells,
                                      const TrackingKernelParameters& params);

  void findRoadsForSchedule(int iteration,
                            const TrackingKernelParameters& params,
                            SeedRefitFunction refitFunction,
                            gsl::span<const CellTopologyId> roadStartCells);

  // Neighbour processing helper; it does not encode a detector layer count.
  template <typename InputSeed>
  void processNeighbours(int iteration, int defaultCellTopologyId, int iLevel, const bounded_vector<InputSeed>& currentCellSeed, const bounded_vector<int>& currentCellId, const bounded_vector<int>& currentCellTopologyId, bounded_vector<TrackSeed>& updatedCellSeed, bounded_vector<int>& updatedCellId, bounded_vector<int>& updatedCellTopologyId, const TrackingKernelParameters& params);

  // Fills scratch transition arrays in binding order after validation. No throw.
  void prepareTransitionScatteringAndBending(int iteration,
                                             const LayerGeometryConfigView& geometryConfig,
                                             const DiskReferenceCoordinateView& referenceCoordinateView,
                                             gsl::span<const TransitionId> transitionIds);

  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  std::shared_ptr<tbb::task_arena> mTaskArena;
  SurfaceGraphView mTraversalGraph{};
  bool mTraversalCacheValid{false};
  // One scalar policy record shared by all surface operations.
  TrackingKernelParameters mKernelParameters{};
  // Borrowed disk reference coordinates, bound once per iteration. Empty for
  // Cylinder iterations; the cylinder path never reads them.
  gsl::span<const float> mDiskLayerReferenceZ{};
  std::vector<float> mDiskLayerReferenceZStorage;
  AttachHitConfigView mAttachHitConfig;
  // Per-position material resolved from SurfaceDescriptor::material through
  // orderedSurfaces. This compatibility cache is not written back to the descriptor.
  //
  // Staged locally and committed with the other traversal caches only when
  // initialiseTimeFrame() succeeds; resetTraversalCache() clears it otherwise.
  // Host-only plan-sized cache indexed by the binding's ordered surface position.
  std::vector<NominalSurfaceMaterial> mLayerMaterial;
  // Non-owning spans into the TimeFrame normalized frame, staged with the
  // traversal cache and cleared on failed initialization.
  std::vector<gsl::span<const SurfaceMeasurement>> mLayerMeasurements;
  std::vector<gsl::span<const GlobalMeasurement>> mLayerGlobalMeasurements;
  // Accepted candidates retained until shared-cluster marking and publication finish.
  bounded_vector<TrackingCandidate> mAcceptedTracksForSharedStatus;
  // Non-owning pointer to the adopted plan binding.
  const SurfacePlanBinding* mBinding = nullptr;

 protected:
  SurfaceTrackingScratch* mScratch = nullptr;
  TimeFrame* mFrame = nullptr;
  gsl::span<const TrackingParameters> mTrkParams;
  float mBz{-999.f};
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKERTRAITS_H_ */

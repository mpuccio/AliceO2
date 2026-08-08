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
#include "ITSMFTTracking/CommonTrack.h"
#include "ITSMFTTracking/SurfaceGraph.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "ITSMFTTracking/detail/CellFinding.h"
#include "ITSMFTTracking/detail/TrackingKernelParameters.h"
#include "ITStracking/BoundedAllocator.h"

namespace o2::itsmft::tracking
{

#ifndef GPUCA_GPUCODE

// Call-scoped native refit function supplied by the detector/workflow edge.
// The generic road stage invokes this function only; publication and reset
// remain outside the tracking transaction.
using SeedRefitFunction = bool (*)(const TrackSeed& seed,
                                   const TrackingParameters& params,
                                   float bz,
                                   SurfaceTrackingScratch& scratch,
                                   gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                                   SurfaceCatalogView surfaceCatalog,
                                   ClusterSourceId expectedSource,
                                   TrackingCandidate& candidate);

#endif

// Full definition lives in TrackletFinding.h (included by
// TrackerTraits.cxx, where the operation itself is called). Only used here
// by const reference in a private method declaration, so a forward
// declaration is sufficient and keeps this public header's dependency
// surface narrow.
struct DiskReferenceCoordinateView;

enum class TraversalFailureReason : uint8_t {
  MissingLayout,
  StaleLayout,
  IterationOutOfRange,
  SparseTopologyMismatch,
  InvalidTraversalSchedule,
  MixedSurfaceKindLayout,
  StateFamilyMismatch,
  InvalidSurfaceParameters,
  // The iteration's index-table configuration failed structural validation.
  InvalidIndexTableConfiguration,
  // A non-FirstPass configuration disagrees with the configuration and LUT
  // already owned by the TimeFrame.
  IndexTableConfigurationMismatch,
  // Legacy LayerxX0 disagrees with the authoritative material for the mapped
  // surface range, or that mapping is invalid/incomplete. Raised before
  // TimeFrame tracking state is touched; the descriptor is never overwritten.
  LegacyMaterialMismatch,
  // The active SurfaceKind does not support the configured MatCorrType.
  // This is structural, so it is reset and rethrown regardless of drop
  // policy. An unrecognized CorrType remains the distinct invalid-mode
  // failure reported by AttachHitConfigView::isValid().
  UnsupportedMaterialCorrectionMode,
  // The per-position normalized-measurement binding disagrees with the
  // loaded frame or legacy compatibility data: counts, surface identity,
  // ClusterRef/source, or external index. Raised before TimeFrame state is
  // touched; the non-owning spans are committed only on success.
  NormalizedMeasurementMismatch,
  // The iteration's orderedSurfaces does not define a bijection from plan
  // positions to global SurfaceId (for example, a duplicate entry).
  SurfaceLayerMappingMismatch,
  // An adopted SurfacePlanBinding cannot translate a traversal
  // TransitionId/CellTopologyId into a compact scratch slot. This indicates
  // a binding/layout mismatch, not event data, and is detected before the
  // offending index addresses scratch storage.
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

// The traits borrow the TimeFrame-owned workspace, plan binding, parameters,
// and refit function for each call. All plan-sized bounds come from the
// adopted runtime plan and scratch.
class TrackerTraits
{
 public:
  virtual ~TrackerTraits() = default;
  // Borrowed call-scoped pointers; neither object owns the other.
  virtual void adoptScratch(SurfaceTrackingScratch* scratch) { mScratch = scratch; }
  virtual void adoptFrame(TimeFrame* frame) { mFrame = frame; }
  // `binding` is borrowed for subsequent Tracker::run() calls and is never
  // owned or copied. It must describe the same graph iteration; a mismatch
  // is reported as TraversalBindingMismatch.
  void adoptSurfacePlanBinding(const SurfacePlanBinding* binding) noexcept { mBinding = binding; }
  // `graphs` is the caller-owned immutable graph vector for this call.
  virtual void initialiseTimeFrame(const int iteration, const std::vector<SurfaceGraph>& graphs);

  virtual void computeLayerTracklets(const int iteration, int iVertex);
  virtual void computeLayerCells(const int iteration);
  virtual void findCellsNeighbours(const int iteration);
  virtual void findRoads(const int iteration, SeedRefitFunction refitFunction);

  void acceptTracks(const TrackingParameters& parameters,
                    bounded_vector<TrackingCandidate>& tracks,
                    bounded_vector<bounded_vector<int>>& firstClusters);

  // The generic result path keeps accepted CommonTrack/TrackSeed pairs until
  // the workflow consumes the successful result. It contains no typed
  // accepted-track vector.
  bounded_vector<TrackingCandidate>& acceptedTracksForSharedStatus();
  void clearAcceptedTracksForSharedStatus();

  void updateTrackingParameters(gsl::span<const TrackingParameters> trkPars)
  {
    mTrkParams = trkPars;
    mTrkParamsByKind = {};
  }
  void updateTrackingParameters(gsl::span<const TrackingParameters> trkPars,
                                gsl::span<const std::array<TrackingParameters, 2>> trkParsByKind)
  {
    mTrkParams = trkPars;
    mTrkParamsByKind = trkParsByKind;
  }
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

  int getTraversalGroupingCount() const noexcept { return mTraversalGroupingCount; }
  bool hasTraversalCache() const noexcept { return mTraversalCacheValid; }
  // Authoritative per-surface-position nominal material resolved once by the
  // most recent successful initialiseTimeFrame() call, from
  // SurfaceDescriptor::material via this iteration's orderedSurfaces mapping
  // -- never inferred from a detector identity or geometry. Valid
  // (and committed) only after initialiseTimeFrame() returns without
  // throwing: a failed call -- for any reason, including one raised after
  // material validation itself already passed -- leaves this exactly as
  // resetTraversalCache() left it (reset/zero-filled), like every other
  // traversal cache; hasTraversalCache() is the existing, single source of
  // truth for whether the most recent call actually succeeded. Read-only
  // test/inspection accessor; production consumption is through
  // mAttachHitConfig, not this span directly.
  gsl::span<const NominalSurfaceMaterial> getLayerMaterial() const noexcept { return {mLayerMaterial.data(), mLayerMaterial.size()}; }
  // Authoritative per-surface-position normalized SurfaceMeasurement span,
  // resolved once by the most recent successful initialiseTimeFrame() call
  // from the already-loaded, already-validated
  // TimeFrame::getNormalizedFrame() via this iteration's orderedSurfaces
  // mapping -- never re-decoded, re-projected, or re-validated per candidate.
  // Same commit/reset contract as getLayerMaterial() above: valid only after
  // a successful initialiseTimeFrame() call, reset to empty spans by
  // resetTraversalCache() otherwise. Every span is non-owning and remains
  // valid only while the owning TimeFrame's normalized frame is alive and
  // has not been cleared/reloaded (see MultiSourceFrame's own lifetime
  // documentation). Read-only test/inspection accessor; production
  // consumption is through the private mLayerMeasurements member directly in
  // computeLayerCellsForKind/processNeighbours.
  gsl::span<const gsl::span<const SurfaceMeasurement>> getLayerMeasurements() const noexcept { return {mLayerMeasurements.data(), mLayerMeasurements.size()}; }

 private:
  const TrackingParameters& parametersForKind(int iteration, SurfaceKind kind) const noexcept
  {
    return mTrkParamsByKind.empty() ? mTrkParams[iteration]
                                    : mTrkParamsByKind[iteration][kind == SurfaceKind::Cylinder ? 0u : 1u];
  }
  void resetTraversalCache() noexcept;
  void validateSparsePlan(int iteration, const SurfaceGraphView& graph, std::optional<SurfaceKind>& activeKind, bool& mixedKind) const;
  int requireSurfacePosition(int iteration, SurfaceId id) const;

  // Sole global-ID to compact-scratch-slot translation used by this class.
  int requireScratchTransitionSlot(int iteration, TransitionId id) const;
  int requireScratchCellSlot(int iteration, CellTopologyId id) const;

  // Non-template operation targets forward to the kind-specific leaves using
  // pre-partitioned ids and kernel parameters.
  void computeLayerTrackletsCylinder(int iteration, int iVertex);
  void computeLayerTrackletsDisk(int iteration, int iVertex);
  void computeLayerCellsCylinder(int iteration);
  void computeLayerCellsDisk(int iteration);
  void findCellsNeighboursCylinder(int iteration);
  void findCellsNeighboursDisk(int iteration);
  void findRoadsCylinder(int iteration, SeedRefitFunction refitFunction);
  void findRoadsDisk(int iteration, SeedRefitFunction refitFunction);

  template <SurfaceKind Kind>
  void computeLayerTrackletsForKind(int iteration,
                                    int iVertex,
                                    gsl::span<const TransitionId> transitionIds,
                                    const TrackingKernelParameters& params);

  template <SurfaceKind Kind>
  void computeLayerCellsForKind(int iteration,
                                gsl::span<const CellTopologyId> cellIds,
                                const TrackingKernelParameters& params);

  template <SurfaceKind Kind>
  void findCellsNeighboursForKind(int iteration,
                                  gsl::span<const CellTopologyId> scheduledCells,
                                  const TrackingKernelParameters& params);

  template <SurfaceKind Kind>
  void findRoadsForKind(int iteration,
                        const TrackingKernelParameters& params,
                        SeedRefitFunction refitFunction,
                        gsl::span<const CellTopologyId> roadStartCells,
                        ClusterSourceId expectedSource);

  // Neighbour processing helper; it does not encode a detector layer count.
  template <SurfaceKind Kind, typename InputSeed>
  void processNeighbours(int iteration, int defaultCellTopologyId, int iLevel, const bounded_vector<InputSeed>& currentCellSeed, const bounded_vector<int>& currentCellId, const bounded_vector<int>& currentCellTopologyId, bounded_vector<TrackSeed>& updatedCellSeed, bounded_vector<int>& updatedCellId, bounded_vector<int>& updatedCellTopologyId, const TrackingKernelParameters& params);

  // Fills scratch-owned transition arrays in the binding's ordered sparse
  // transition order after fallible validation. This method does not throw.
  template <SurfaceKind Kind>
  void prepareTransitionScatteringAndBendingForKind(int iteration,
                                                    const LayerGeometryConfigView& geometryConfig,
                                                    const DiskReferenceCoordinateView& referenceCoordinateView,
                                                    gsl::span<const TransitionId> transitionIds);

  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  std::shared_ptr<tbb::task_arena> mTaskArena;
  SurfaceGraphView mTraversalGraph{};
  bool mTraversalCacheValid{false};
  // One committed record for the active endpoint SurfaceKind. It is reset
  // together with the other traversal caches and published only after all
  // fallible initialization checks succeed.
  std::array<TrackingKernelParameters, 2> mKernelParameters{};
  // Borrowed disk reference coordinates, bound once per iteration. Empty for
  // Cylinder iterations; the cylinder path never reads them.
  gsl::span<const float> mDiskLayerReferenceZ{};
  std::vector<float> mDiskLayerReferenceZStorage;
  AttachHitConfigView mAttachHitConfig;
  // Authoritative per-surface-position nominal material, resolved once per
  // initialiseTimeFrame() from SurfaceDescriptor::material via this
  // iteration's orderedSurfaces mapping (never inferred from detector identity,
  // radius, z, or numeric ordering). SurfaceDescriptor::material is
  // authoritative and is never overwritten here; this is a compatibility
  // cache for the current leaf operations.
  //
  // Commit contract: resolved into a local scratch array first and only
  // copied here at the very end of initialiseTimeFrame(), alongside every
  // other traversal cache -- never as soon as
  // material validation itself passes. A later fallible check in the same
  // call (index-table binding, legacy topology parity, state-family, or any
  // other surface-kind/geometry validation) must not leave this populated from an
  // iteration that ultimately failed; resetTraversalCache() zero-fills it at
  // the top of every call, and it stays that way unless the call returns
  // normally. See getLayerMaterial()'s doc for the read-side contract.
  // Host-only plan-sized cache. The position is the adopted binding's ordered
  // surface position; current leaf operations consume these per-position
  // parameter spans directly.
  std::vector<NominalSurfaceMaterial> mLayerMaterial;
  std::optional<SurfaceKind> mActiveKind;
  std::array<std::vector<TransitionId>, 2> mTransitionsByKind;
  std::array<std::vector<CellTopologyId>, 2> mCellsByKind;
  std::array<std::vector<CellTopologyId>, 2> mScheduledCellsByKind;
  std::array<std::vector<CellTopologyId>, 2> mRoadStartCellsByKind;
  // Non-owning per-surface-position spans into the TimeFrame-owned normalized
  // frame. They are staged and committed with the traversal cache, then
  // cleared on failed initialization.
  // Host-only non-owning view, sized to the adopted ordered surface span.
  std::vector<gsl::span<const SurfaceMeasurement>> mLayerMeasurements;
  int mTraversalGroupingCount{0};
  // Generic accepted candidates are retained only until shared-cluster
  // marking and the final adapter-owned publication seal complete.
  bounded_vector<TrackingCandidate> mAcceptedTracksForSharedStatus;
  // Non-owning pointer to the adopted plan binding.
  const SurfacePlanBinding* mBinding = nullptr;

 protected:
  SurfaceTrackingScratch* mScratch = nullptr;
  TimeFrame* mFrame = nullptr;
  gsl::span<const TrackingParameters> mTrkParams;
  gsl::span<const std::array<TrackingParameters, 2>> mTrkParamsByKind;
  float mBz{-999.f};
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKERTRAITS_H_ */

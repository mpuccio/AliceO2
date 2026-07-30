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

#include <gsl/span>
#include <oneapi/tbb.h>

#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TransitionPolicyBinding.h"
#include "ITSMFTTracking/TransitionPolicyDispatch.h"
#include "ITSMFTTracking/TransitionPolicyState.h"
#include "ITStracking/BoundedAllocator.h"

namespace o2::itsmft::tracking
{

// Full definition lives in TransitionPolicyOperations.h (included by
// TrackerTraits.cxx, where the operation itself is called). Only used here
// by const reference in a private method declaration, so a forward
// declaration is sufficient and keeps this public header's dependency
// surface narrow.
struct DiskDiskReferenceCoordinateView;

enum class TraversalFailureReason : uint8_t {
  MissingLayout,
  StaleLayout,
  IterationOutOfRange,
  LegacyIndexMismatch,
  InvalidTraversalSchedule,
  MixedPolicyLayout,
  StateFamilyMismatch,
  InvalidPolicyParameters,
  // The index-table configuration bound from this iteration's
  // TrackingParameters (IndexTableConfiguration.h) failed its own
  // structural validation.
  InvalidIndexTableConfiguration,
  // A non-FirstPass iteration's freshly bound index-table configuration
  // disagrees with the configuration (and populated LUT) the TimeFrame
  // already owns, which this iteration would otherwise silently reuse or
  // resort clusters into.
  IndexTableConfigurationMismatch,
  // Temporary compatibility check (Stage-B surface-material fast-fit slice):
  // this iteration's legacy TrackingParameters::LayerxX0 does not match the
  // authoritative SurfaceDescriptor::material resolved for the same layer
  // range, or the layer-to-surface mapping needed to resolve it is itself
  // invalid/incomplete. Raised before any TimeFrame tracking state is
  // touched; SurfaceDescriptor::material is never overwritten as a result.
  LegacyMaterialMismatch,
  // Stage-B activation constraint: policy preflight, wired into
  // initialiseTimeFrame() (see checkMaterialCorrectionModeSupport() in
  // TransitionPolicyBinding.h, called once per iteration right after
  // grouping/mixed-policy validation and before material staging/index-table
  // binding/TimeFrame::initialise()). The active transition policy's Stage-A
  // native SurfaceKinematicState path does not support this iteration's
  // configured MatCorrType: currently only CylinderCylinder with a non-NONE,
  // otherwise-recognized CorrType (USEMatCorrLUT/USEMatCorrTGeo remain
  // legacy-ITS-only). This is a structural/configuration failure like every
  // other TraversalException reason above: it is always wiped and rethrown
  // regardless of DropTFUponFailure -- it must never be reported as the
  // dropped-TF sentinel. A recognized-invalid CorrType value is a distinct,
  // pre-existing configuration-validation failure
  // (AttachHitPolicyConfigView::isValid()) and must never be reported as
  // this reason.
  UnsupportedMaterialCorrectionMode,
  // Stage-B normalized-CA-measurements slice: this iteration's one-time
  // normalized-measurement binding (mLayerMeasurements) failed its own
  // compatibility validation against the already-loaded normalized frame and
  // legacy compatibility structures -- covers a normalized measurement count
  // that disagrees with the corresponding legacy unsorted-cluster count, a
  // TrackingFrameInfo count that disagrees with that same count, a
  // measurement whose SurfaceId does not match the expected mapped surface,
  // an invalid measurement.cluster ClusterRef, a nonzero measurement source
  // under the current single-source TimeFrame contract, a negative legacy
  // external index, or a measurement.cluster.index that disagrees with the
  // corresponding legacy external index. Raised before any TimeFrame
  // tracking state is touched, alongside the LegacyMaterialMismatch check
  // this mirrors; mLayerMeasurements is never partially populated as a
  // result -- see mLayerMeasurements' own doc for the commit contract.
  NormalizedMeasurementMismatch,
  // Gate 4 Slice 0a (sparse-topology tracklet migration): this iteration's
  // orderedSurfaces does not define a bijection from legacy layer index onto
  // global SurfaceId -- i.e. two distinct legacy layers map to the same
  // SurfaceId. Detected in the same orderedSurfaces walk that resolves
  // mSurfaceToLegacyLayer (see that member's doc), immediately alongside the
  // existing per-entry LegacyMaterialMismatch validity/range check, before
  // any TimeFrame tracking state is touched; mSurfaceToLegacyLayer is never
  // partially populated as a result.
  SurfaceLayerMappingMismatch
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
  // Authoritative per-(legacy-)layer nominal material resolved once by the
  // most recent successful initialiseTimeFrame() call, from
  // SurfaceDescriptor::material via this iteration's orderedSurfaces mapping
  // -- never inferred from legacy index/detector identity/geometry. Valid
  // (and committed) only after initialiseTimeFrame() returns without
  // throwing: a failed call -- for any reason, including one raised after
  // material validation itself already passed -- leaves this exactly as
  // resetTraversalCache() left it (reset/zero-filled), like every other
  // traversal cache; hasTraversalCache() is the existing, single source of
  // truth for whether the most recent call actually succeeded. Read-only
  // test/inspection accessor; production consumption is through
  // mAttachHitConfig (TransitionPolicyBinding.h), not this span directly.
  gsl::span<const NominalSurfaceMaterial> getLayerMaterial() const noexcept { return {mLayerMaterial.data(), mLayerMaterial.size()}; }
  // Authoritative per-(legacy-)layer normalized SurfaceMeasurement span,
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
  // computeLayerCellsForPolicy/processNeighbours.
  gsl::span<const gsl::span<const SurfaceMeasurement>> getLayerMeasurements() const noexcept { return {mLayerMeasurements.data(), mLayerMeasurements.size()}; }

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
  // Gate 3 cell-road pre-cut slice: the legacy MFT reference-z span bound
  // once per iteration by bindLegacyMFTReferenceCoordinates() (static-storage
  // duration, never dangles), reused by passesCellRoadPrecut<DiskDisk> across
  // every candidate. Stored as the raw span -- not the full
  // DiskDiskReferenceCoordinateView -- so this header does not need that
  // type complete (see the forward declaration above). Left empty (default)
  // for CylinderCylinder iterations; passesCellRoadPrecut<CylinderCylinder>
  // never reads it.
  gsl::span<const float> mDiskLayerReferenceZ{};
  AttachHitPolicyConfigView mAttachHitConfig;
  // Authoritative per-(legacy-)layer nominal material, resolved once per
  // initialiseTimeFrame() from SurfaceDescriptor::material via this
  // iteration's orderedSurfaces mapping (never inferred from legacy index,
  // detector identity, radius, z, or numeric ordering). Compatibility-only:
  // SurfaceDescriptor::material is authoritative and never overwritten here;
  // this cache is temporary duplication that disappears once the final ITS
  // refit migrates off TrackingParameters::LayerxX0.
  //
  // Commit contract: resolved into a local scratch array first and only
  // copied here at the very end of initialiseTimeFrame(), alongside every
  // other traversal cache (mTraversalGrouping et al.) -- never as soon as
  // material validation itself passes. A later fallible check in the same
  // call (index-table binding, legacy topology parity, state-family, or any
  // other policy/geometry validation) must not leave this populated from an
  // iteration that ultimately failed; resetTraversalCache() zero-fills it at
  // the top of every call, and it stays that way unless the call returns
  // normally. See getLayerMaterial()'s doc for the read-side contract.
  std::array<NominalSurfaceMaterial, NLayers> mLayerMaterial{};
  // Gate 4 Slice 0a (sparse-topology tracklet migration): temporary bridge
  // from a global SurfaceId back to this TrackerTraits<NLayers>'s own
  // legacy layout-local layer index, so the migrated hot loops can resolve
  // a sparse SurfaceTransition's `from`/`to` endpoints to the legacy layer
  // indices TimeFrame's per-layer storage (clusters, index tables,
  // TrackingParameters::LayerRadii, mLayerMaterial, mLayerMeasurements) is
  // still keyed by. Sized to the full global SurfaceId domain
  // (MaxLayoutSurfaces, SurfaceId.h) -- never NLayers -- because SurfaceId
  // numbering is global, not per-detector; only the (at most) NLayers
  // entries this iteration's orderedSurfaces actually maps are ever valid
  // for this bridge, every other slot stays kInvalidLegacyLayer. Same
  // staged-then-committed contract as mLayerMaterial immediately above:
  // resolved into a local scratch array first, alongside mLayerMaterial, in
  // the same orderedSurfaces walk (see initialiseTimeFrame()'s step 2.5),
  // and committed here only in the final traversal-cache commit block.
  // resetTraversalCache() sentinel-fills every element at the top of every
  // call, and it stays that way unless the call returns normally.
  static constexpr uint8_t kInvalidLegacyLayer = 0xFFu;
  std::array<uint8_t, MaxLayoutSurfaces> mSurfaceToLegacyLayer{};
  // One-time normalized-measurement binding (Stage-B normalized-CA-
  // measurements slice): non-owning per-(legacy-)layer span into the
  // TimeFrame-owned normalized frame, resolved and validated once per
  // initialiseTimeFrame() call from this iteration's orderedSurfaces mapping
  // -- never an owning container (Architecture.md Sec 7: the common
  // TimeFrame/TrackerTraits own no compact clusters or duplicated
  // measurement storage). Same staged-then-committed contract as
  // mLayerMaterial above: resolved into a local scratch array first and
  // committed here only in the final traversal-cache commit block, alongside
  // mLayerMaterial and every other successfully staged state; a later
  // fallible check in the same call must not leave this populated from an
  // iteration that ultimately failed. resetTraversalCache() resets every
  // element to an empty span at the top of every call, and it stays that way
  // unless the call returns normally.
  std::array<gsl::span<const SurfaceMeasurement>, NLayers> mLayerMeasurements{};
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

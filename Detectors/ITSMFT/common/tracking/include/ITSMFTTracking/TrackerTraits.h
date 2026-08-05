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
#include "ITSMFTTracking/AcceptedTrackShadowPublisher.h"
#include "ITSMFTTracking/DetectorLayoutSet.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "ITSMFTTracking/detail/TransitionPolicyBinding.h"
#include "ITSMFTTracking/detail/TransitionPolicyDispatch.h"
#include "ITSMFTTracking/detail/TransitionPolicyState.h"
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
  SurfaceLayerMappingMismatch,
  // Gate 4 C2 Slice 1: an adopted SurfacePlanBinding could not translate a
  // global TransitionId/CellTopologyId encountered during traversal into a
  // compact scratch slot -- SurfacePlanBinding::getScratchTransitionSlot()/
  // getScratchCellSlot() returned std::nullopt. Since a correctly-built
  // binding structurally cannot hand a hot loop a foreign/unmapped id, this
  // can only mean the binding disagrees
  // with the layout/topology this iteration is actually running against
  // (e.g. a stale binding, or a binding built for the wrong DetectorLayoutSet
  // iteration) -- always a structural/configuration failure, never a per-TF
  // data problem, so DropTFUponFailure never applies. Raised before the
  // offending index is ever used to address scratch storage.
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

// The generic tracker owns one detector-neutral workspace and one immutable
// plan binding. NLayers remains an algorithm/storage parameter until M6g.
template <int NLayers>
class TrackerTraits
{
 public:
  using IndexTableUtilsN = o2::itsmft::IndexTableUtils<NLayers>;

  virtual ~TrackerTraits() = default;
  // Two independent bind-once pointers -- neither owns nor stores a
  // reference to the other.
  virtual void adoptScratch(SurfaceTrackingScratch* scratch) { mScratch = scratch; }
  virtual void adoptFrame(TimeFrame* frame) { mFrame = frame; }
  void adoptMFTPublicationCompatibility(MFTPublicationCompatibility* compatibility) noexcept { mAcceptedTrackShadowPublisher.adoptMFTPublicationCompatibility(compatibility); }
  void adoptITSSharedClusterCompatibility(ITSSharedClusterCompatibility* compatibility) noexcept { mAcceptedTrackShadowPublisher.adoptITSSharedClusterCompatibility(compatibility); }
  // M6f: bind the one SurfacePlanBinding used by the common-CA hot loops.
  // `binding` must outlive every subsequent clustersToTracks() call and is
  // never owned or copied. Direct algorithm tests may omit it when their
  // topology is already identity-indexed; production adapters always adopt
  // the binding built for the same DetectorLayoutSet iteration. A non-identity
  // binding mismatch surfaces as TraversalFailureReason::TraversalBindingMismatch,
  // never as a silent misread.
  void adoptSurfacePlanBinding(const SurfacePlanBinding* binding) noexcept { mBinding = binding; }
  // `layouts` is the owner's (ITSMFTTrackingInterface's) one immutable plan,
  // supplied explicitly by the caller (Gate 4 B2 Slice 2) -- this no longer
  // reads any layout/catalog state off TimeFrame.
  virtual void initialiseTimeFrame(const int iteration, const DetectorLayoutSet& layouts);

  virtual void computeLayerTracklets(const int iteration, int iVertex);
  virtual void computeLayerCells(const int iteration);
  virtual void findCellsNeighbours(const int iteration);
  virtual void findRoads(const int iteration);

  void acceptTracks(int iteration, bounded_vector<CATrackType<NLayers>>& tracks, bounded_vector<bounded_vector<int>>& firstClusters);
  void markTracks(int iteration);

  // The Surface path keeps accepted tracks only until ITS shared-cluster
  // compatibility is sealed. This is not publication staging: writers read
  // TimeFrame CommonTrack data exclusively. Legacy scratch users retain
  // their existing Group-C members until M6f.
  bounded_vector<CATrackType<NLayers>>& acceptedTracksForSharedStatus();
  void clearAcceptedTracksForSharedStatus();

  void updateTrackingParameters(const std::vector<TrackingParameters>& trkPars) { mTrkParams = trkPars; }
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

  // Gate 4 C2 Slice 1: the sole global-TransitionId/CellTopologyId-to-
  // compact-scratch-slot translation used anywhere in this class. Called
  // exactly once per traversal loop head or per dynamically-discovered
  // neighbour id (never inside a per-candidate inner loop -- see call sites
  // in TrackerTraits.cxx). With no binding adopted (mBinding == nullptr),
  // returns the id's raw .value() unchanged, i.e. an identity mapping --
  // this is exactly today's Gate 3 behavior, where the detector-only static
  // catalog already makes global ids dense/compact per detector. With a
  // binding adopted, returns its checked scratch slot or throws
  // TraversalFailureReason::TraversalBindingMismatch before the caller can
  // use an unmapped id to index scratch.
  int requireScratchTransitionSlot(int iteration, TransitionId id) const;
  int requireScratchCellSlot(int iteration, CellTopologyId id) const;

  // Gate 4 C2 Slice 1: outer-loop dispatch wrapper used everywhere this class
  // previously called dispatchTransitionPolicies(grouping, visitor) directly.
  // With no binding adopted, forwards unchanged to dispatchTransitionPolicies
  // (identical behavior to before this slice). With a binding adopted,
  // `grouping` (built from the possibly-multi-detector `layout`) is not used
  // for tag selection at all: the active tag is this NLayers instantiation's
  // own compile-time family (a binding only ever owns transitions of
  // the matching tag, enforced by SurfacePlanBinding::build()), so the
  // visitor is invoked exactly once, fed `binding`'s own filtered
  // transition/cell spans instead of `grouping`'s whole-layout ones -- this
  // is what prevents a combined (both-tags-present) `grouping` from firing
  // the visitor twice, once per tag, for a single-detector TrackerTraits.
  template <typename Visitor>
  void dispatchActivePolicy(const TransitionPolicyGrouping& grouping, Visitor&& visitor) const;

  // M5c: the compact, operation-local replacement for the four
  // detector-family/TransitionPolicyTag runtime branches that used to live
  // directly in computeLayerTracklets()/computeLayerCells()/
  // findCellsNeighbours()/findRoads() (one dispatchActivePolicy() call each,
  // every call). Holds only already-bound member-function pointers -- never
  // the TransitionPolicyTag (or the family it maps to) that selected them --
  // so those four shared hot-loop entry points below consume it with no Tag/
  // family branch of their own. Plain pointers-to-member, deliberately not a
  // type-erasing callable wrapper: each callable's target is one of the
  // eight non-template wrapper methods below (two per operation, one per
  // TransitionPolicyTag), whose signature is fixed and identical for both
  // tags -- so a plain pointer-to-member suffices, with no type erasure, no
  // capture storage, and no possibility of a heap allocation. The wrapper it
  // points to forwards to the existing
  // Tag-templated *ForPolicy leaf implementation, unchanged, using the ids
  // this same binding stores below and the corresponding
  // mCylinderPolicyParams/mDiskPolicyParams. bindTraversalOperation()
  // (TrackerTraits.cxx) is this struct's only producer: it fills every
  // member exactly once per successful initialiseTimeFrame() call, from
  // that call's already-validated activeTag/params/grouping (activeTag itself
  // derived earlier in the same call from actual endpoint SurfaceDescriptor
  // kinds via validateLegacyParity()'s tagOf(), never from NLayers or
  // detector identity). resetTraversalCache() clears it back to unbound,
  // alongside every other traversal cache, so a hot loop can never observe a
  // binding left over from a failed or unrelated iteration.
  struct TraversalOperationBinding {
    using ComputeTrackletsFn = void (TrackerTraits::*)(int iteration, int iVertex);
    using ComputeCellsFn = void (TrackerTraits::*)(int iteration);
    using FindNeighboursFn = void (TrackerTraits::*)(int iteration);
    using FindRoadsFn = void (TrackerTraits::*)(int iteration);

    ComputeTrackletsFn computeTracklets = nullptr;
    ComputeCellsFn computeCells = nullptr;
    FindNeighboursFn findNeighbours = nullptr;
    FindRoadsFn findRoads = nullptr;
    // The ids computeTracklets/computeCells/findNeighbours were bound
    // against -- resolved once by bindTraversalOperation(), from the same
    // grouping/binding lookup the removed per-call dispatch used to redo on
    // every call. Not a tag/family/detector-id/SurfaceKindPair: plain sparse
    // topology ids, non-owning, valid for exactly the traversal cache's own
    // lifetime (see resetTraversalCache()). findRoads needs none: its target
    // wrapper reads mTraversalGrouping/mBinding's road-start cells directly,
    // exactly as findRoadsForPolicy<Tag> already did before this slice.
    gsl::span<const TransitionId> boundTransitionIds{};
    gsl::span<const CellTopologyId> boundCellIds{};
    gsl::span<const CellTopologyId> boundScheduledCellIds{};
    bool bound = false;
  };

  // Builds mTraversalOperation for this iteration. Called once, at the very
  // end of initialiseTimeFrame(), after every other fallible validation and
  // traversal-cache commit in that call has already succeeded -- see that
  // method and TraversalOperationBinding's own doc above.
  void bindTraversalOperation(int iteration);

  // The eight non-template targets TraversalOperationBinding's member-
  // function pointers may point to -- one per (operation, TransitionPolicyTag)
  // pair, selected only inside bindTraversalOperation(). Each is a thin,
  // non-template forwarder to the corresponding *ForPolicy<Tag> leaf
  // implementation below, reading mTraversalOperation's bound ids and the
  // matching mCylinderPolicyParams/mDiskPolicyParams. Never called directly
  // from the four public hot-loop entry points -- only through
  // mTraversalOperation's bound pointer.
  void computeLayerTrackletsCylinderCylinder(int iteration, int iVertex);
  void computeLayerTrackletsDiskDisk(int iteration, int iVertex);
  void computeLayerCellsCylinderCylinder(int iteration);
  void computeLayerCellsDiskDisk(int iteration);
  void findCellsNeighboursCylinderCylinder(int iteration);
  void findCellsNeighboursDiskDisk(int iteration);
  void findRoadsCylinderCylinder(int iteration);
  void findRoadsDiskDisk(int iteration);

  template <TransitionPolicyTag Tag>
  void computeLayerTrackletsForPolicy(int iteration,
                                      int iVertex,
                                      gsl::span<const TransitionId> transitionIds,
                                      const typename TransitionPolicyTraits<Tag>::Params& params);

  template <TransitionPolicyTag Tag>
  void computeLayerCellsForPolicy(int iteration,
                                  gsl::span<const CellTopologyId> cellIds,
                                  const typename TransitionPolicyTraits<Tag>::Params& params);

  template <TransitionPolicyTag Tag>
  void findCellsNeighboursForPolicy(int iteration,
                                    gsl::span<const CellTopologyId> scheduledCells,
                                    const typename TransitionPolicyTraits<Tag>::Params& params);

  template <TransitionPolicyTag Tag>
  void findRoadsForPolicy(int iteration, const typename TransitionPolicyTraits<Tag>::Params& params);

  // M4 (GenericTrackingEngineMigration.md; ADR 0007 decision 7): moved from
  // public to private -- findRoadsForPolicy() is this method's only caller
  // (TrackerTraits.cxx), recursively, so there was never an external need
  // for it to be publicly callable. Explicit template instantiation of a
  // private member is unaffected by access control (TrackerTraits.cxx still
  // instantiates it for every (NLayers, Tag) pair it needs).
  template <TransitionPolicyTag Tag, typename InputSeed>
  void processNeighbours(int iteration, int defaultCellTopologyId, int iLevel, const bounded_vector<InputSeed>& currentCellSeed, const bounded_vector<int>& currentCellId, const bounded_vector<int>& currentCellTopologyId, bounded_vector<TrackSeed>& updatedCellSeed, bounded_vector<int>& updatedCellId, bounded_vector<int>& updatedCellTopologyId, const typename TransitionPolicyTraits<Tag>::Params& params);

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
  // M5c: see TraversalOperationBinding's own doc above.
  TraversalOperationBinding mTraversalOperation;
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
  // Host-only plan-sized cache.  The legacy layer index remains an adapter
  // coordinate until M7c, but its extent is supplied by the adopted plan,
  // never by the TrackerTraits template argument.
  std::vector<NominalSurfaceMaterial> mLayerMaterial;
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
  // Host-only non-owning view, sized to the adopted ordered surface span.
  std::vector<gsl::span<const SurfaceMeasurement>> mLayerMeasurements;
  int mTraversalGroupingCount{0};
  // A template-specialized serial accepted-track hook. Its ITS instantiation
  // is empty and has no MFT compatibility allocation or lookup.
  AcceptedTrackShadowPublisher<NLayers> mAcceptedTrackShadowPublisher;
  bounded_vector<CATrackType<NLayers>> mAcceptedTracksForSharedStatus;
  // M6f: non-owning pointer to the adopted common-CA plan binding. It is
  // populated by production adapters before traversal; nullptr is retained
  // only for identity-indexed direct algorithm tests (see above).
  const SurfacePlanBinding* mBinding = nullptr;

 protected:
  SurfaceTrackingScratch* mScratch = nullptr;
  TimeFrame* mFrame = nullptr;
  std::vector<TrackingParameters> mTrkParams;
  float mBz{-999.f};
};

using TrackerTraitsITS = TrackerTraits<ITSNLayers>;
using TrackerTraitsMFT = TrackerTraits<o2::mft::constants::mft::LayersNumber>;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKERTRAITS_H_ */

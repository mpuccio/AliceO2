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
/// \file TrackerTraits.cxx
/// \brief
///

#include <algorithm>
#include <array>
#include <iterator>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/enumerable_thread_specific.h>

#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "GPUCommonMath.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/CommonTrackShadow.h"
#include "ITStracking/Constants.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/DetectorTraits.h"
#include "ITSMFTTracking/IndexTableConfiguration.h"
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/SurfaceKinematicStateLegacyAdapters.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/LayerMask.h"
#include "ITStracking/ROFLookupTables.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TransitionPolicyBinding.h"
#include "ITSMFTTracking/TransitionPolicyOperations.h"
#include "ITStracking/TrackHelpers.h"
#include "ITStracking/Tracklet.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::itsmft::tracking
{

namespace math_utils = o2::its::math_utils;
using o2::its::deepVectorClear;
using o2::its::TimeEstBC;

struct PassMode {
  using OnePass = std::integral_constant<int, 0>;
  using TwoPassCount = std::integral_constant<int, 1>;
  using TwoPassInsert = std::integral_constant<int, 2>;
};

namespace
{
// A static "diamond" vertex (ITSCommonCATrackerParam.useDiamond /
// TrackerTraits::computeLayerTrackletsForPolicy) carries no genuine
// per-event timing: it stands in for every real primary vertex at once, so
// it must compare as time-compatible with whichever ROF it is being tested
// against. Vertex::getTimeStamp() is a fixed-width TimeEstBC whose error
// field is only 16 bits (DataFormatsITS/TimeEstBC.h) -- sized for one ROF's
// own uncertainty, not for an entire TimeFrame's real BC span (order
// 1e5-1e6 BC for a realistic TF) -- so one static timestamp value literally
// cannot overlap every ROF in a real TimeFrame at once; the correct
// full-TimeFrame envelope has to be re-derived at each (layer, rofId) this
// vertex is actually tested against, from that ROF's own real interval
// (ROFOverlapTableView::getLayer(layer).getROFTimeBounds(rofId, true) -- the
// exact same, already-configured per-TF bounds every real-vertex ROF
// compatibility check in this file already uses; not an independently
// invented arithmetic formula). A vertex timestamp built this way is
// therefore provably -- not just observed to be -- compatible with the ROF
// it was derived from: see this function's two call sites below for the two
// provable cases (an exact-match window, and a tracklet-time window that is
// always a subset of one of its two source ROF windows).
template <typename ROFOverlapView>
Vertex diamondVertexForROF(const Vertex& base, const ROFOverlapView& rofOverlapView, int layer, int rofId)
{
  Vertex v = base;
  v.setTimeStamp(rofOverlapView.getLayer(layer).getROFTimeBounds(rofId, true));
  return v;
}

template <int NLayers>
bool makeAcceptedTrackShadow(const CATrackType<NLayers>& track,
                             const LayerMeasurementSpans<NLayers>& layerMeasurements,
                             const TimeEstBC& selectedTimestamp,
                             CommonTrackShadowRecord& record)
{
  record = {};
  if constexpr (DetectorTraits<NLayers>::DetId == o2::detectors::DetID::ITS) {
    if (!legacy::importBarrelTrackParCov(track.getParamIn(), record.track.innerState) ||
        !legacy::importBarrelTrackParCov(track.getParamOut(), record.track.outerState)) {
      return false;
    }
  } else {
    if (!legacy::importLegacyForwardTrackParCov(track.getTrack(), record.track.innerState) ||
        !legacy::importLegacyForwardTrackParCov(track.getTrack().getOutParam(), record.track.outerState)) {
      return false;
    }
  }
  record.track.chi2 = track.getChi2();
  record.track.timestamp = {static_cast<TFBC>(selectedTimestamp.lower()), static_cast<TFBC>(selectedTimestamp.upper())};
  record.references.reserve(NLayers);
  for (int layer = 0; layer < NLayers; ++layer) {
    const int localIndex = track.getClusterIndex(layer);
    if (localIndex == constants::UnusedIndex) {
      continue;
    }
    if (localIndex < 0 || static_cast<size_t>(localIndex) >= layerMeasurements[layer].size()) {
      return false;
    }
    const auto& measurement = layerMeasurements[layer][localIndex];
    const TrackClusterReference reference{measurement.surface, SurfaceMeasurementIndex{static_cast<uint32_t>(localIndex)}};
    if (!reference.surface.isValid() || measurement.surface != reference.surface || !measurement.cluster.isValid()) {
      return false;
    }
    record.references.push_back(reference);
    record.track.hitSurfaces.set(reference.surface);
  }
  return !record.references.empty();
}
} // namespace

template <int NLayers>
void TrackerTraits<NLayers>::resetTraversalCache() noexcept
{
  mTraversalLayout = {};
  mTraversalGrouping.reset();
  mCylinderPolicyParams.reset();
  mDiskPolicyParams.reset();
  mDiskLayerReferenceZ = {};
  mAttachHitConfig = {};
  mLayerMaterial.fill(NominalSurfaceMaterial{});
  mSurfaceToLegacyLayer.fill(kInvalidLegacyLayer);
  mLayerMeasurements.fill(gsl::span<const SurfaceMeasurement>{});
  mTraversalGroupingCount = 0;
  mPolicyBindingCounts.fill(0);
}

template <int NLayers>
int TrackerTraits<NLayers>::getPolicyBindingCount(TransitionPolicyTag tag) const noexcept
{
  switch (tag) {
    case TransitionPolicyTag::CylinderCylinder:
      return mPolicyBindingCounts[0];
    case TransitionPolicyTag::DiskDisk:
      return mPolicyBindingCounts[1];
    case TransitionPolicyTag::Invalid:
      return 0;
  }
  return 0;
}

template <int NLayers>
void TrackerTraits<NLayers>::validateLegacyParity(int iteration,
                                                  const DetectorLayoutView& layout,
                                                  TransitionPolicyTag& activeTag,
                                                  bool& mixedPolicy) const
{
  const auto fail = [iteration]() { throw TraversalException{iteration, TraversalFailureReason::LegacyIndexMismatch}; };
  const auto sparse = layout.topology;
  const auto legacy = mScratch->getTrackingTopologyView();
  using LegacyId = typename ScratchN::TrackingTopologyN::Id;
  if (sparse.nTransitions != legacy.nTransitions || sparse.nCells != legacy.nCells ||
      (sparse.nTransitions != 0 && (sparse.transitions == nullptr || sparse.cellsByFirstTransitionOffsets == nullptr)) ||
      (sparse.nCells != 0 && (sparse.cells == nullptr || sparse.cellsByFirstTransition == nullptr))) {
    fail();
  }

  activeTag = TransitionPolicyTag::Invalid;
  mixedPolicy = false;
  for (uint32_t id = 0; id < sparse.nTransitions; ++id) {
    const auto sparseId = TransitionId{static_cast<uint16_t>(id)};
    const auto& current = sparse.getTransition(sparseId);
    const auto& reference = legacy.getTransition(static_cast<LegacyId>(id));
    if (current.from.value() != reference.fromLayer || current.to.value() != reference.toLayer ||
        current.skippedSurfaces.value() != LayerMask::skipped(reference.fromLayer, reference.toLayer).value()) {
      fail();
    }
    if (activeTag == TransitionPolicyTag::Invalid) {
      activeTag = current.policyTag;
    } else if (current.policyTag != activeTag) {
      mixedPolicy = true;
    }
  }

  for (uint32_t id = 0; id < sparse.nCells; ++id) {
    const auto sparseId = CellTopologyId{static_cast<uint16_t>(id)};
    const auto& current = sparse.getCell(sparseId);
    const auto& reference = legacy.getCell(static_cast<LegacyId>(id));
    if (current.firstTransition.value() != reference.firstTransition ||
        current.secondTransition.value() != reference.secondTransition ||
        current.hitSurfaces.value() != reference.hitLayerMask.value()) {
      fail();
    }
    const auto firstTag = sparse.getTransition(current.firstTransition).policyTag;
    const auto secondTag = sparse.getTransition(current.secondTransition).policyTag;
    if (firstTag != secondTag || (activeTag != TransitionPolicyTag::Invalid && firstTag != activeTag)) {
      mixedPolicy = true;
    }
  }

  if (sparse.seedingSurfaces.value() != mTrkParams[iteration].StartLayerMask.value() ||
      sparse.nCells != legacy.nCellsByFirstTransition) {
    fail();
  }
  for (uint32_t transition = 0; transition < sparse.nTransitions; ++transition) {
    const auto sparseRange = sparse.getCellsStartingWithTransition(TransitionId{static_cast<uint16_t>(transition)});
    const auto legacyRange = legacy.getCellsStartingWithTransition(static_cast<LegacyId>(transition));
    if (sparse.cellsByFirstTransitionOffsets[transition] != sparseRange.getFirstEntry() ||
        sparse.cellsByFirstTransitionOffsets[transition + 1] != sparseRange.getEntriesBound() ||
        sparseRange.getFirstEntry() != legacyRange.getFirstEntry() ||
        sparseRange.getEntries() != legacyRange.getEntries()) {
      fail();
    }
  }
  for (uint32_t entry = 0; entry < sparse.nCells; ++entry) {
    if (sparse.cellsByFirstTransition[entry].value() != legacy.cellsByFirstTransition[entry]) {
      fail();
    }
  }
}

template <int NLayers>
void TrackerTraits<NLayers>::initialiseTimeFrame(const int iteration, const DetectorLayoutSet& layouts)
{
  resetTraversalCache();

  // 1. Iteration bounds, checked before any index-table (or other)
  // configuration is touched. `layouts` is the caller's own immutable plan
  // (Gate 4 B2 Slice 2): its mere presence as a valid reference here is the
  // caller's guarantee that a plan exists, so there is no separate
  // missing/stale-layout state left to check.
  if (iteration < 0 || static_cast<size_t>(iteration) >= layouts.size()) {
    throw TraversalException{iteration, TraversalFailureReason::IterationOutOfRange};
  }
  const auto layout = layouts.getLayoutView(iteration);

  // 2. Grouping + single active tag, resolved from `layout` alone -- no
  // dependency on TimeFrame::initialise() having run, since neither
  // TransitionPolicyGrouping nor dispatchTransitionPolicies read the legacy
  // topology view that call populates.
  TransitionPolicyGrouping grouping{layout};
  ++mTraversalGroupingCount;
  if (!grouping.valid()) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidTraversalSchedule};
  }
  if (grouping.hasTag(TransitionPolicyTag::CylinderCylinder) && grouping.hasTag(TransitionPolicyTag::DiskDisk)) {
    throw TraversalException{iteration, TraversalFailureReason::MixedPolicyLayout};
  }

  // 2.1 (Stage-B activation, Architecture.md Sec 11): material-correction-mode
  // preflight for the active policy, once per iteration -- after layout/
  // grouping validation and mixed-policy resolution, before any material
  // staging, index-table binding, or TimeFrame tracking-state mutation.
  // `grouping` alone (not yet `activeTag`, which validateLegacyParity() only
  // resolves later) already tells us which single tag is active, exactly
  // like step 3's index-table dispatch below. An `Unsupported` result is a
  // structural/configuration failure (TraversalFailureReason doc); an
  // `InvalidMode` result is deliberately not raised here -- it defers to the
  // existing AttachHitPolicyConfigView::isValid() check further below, which
  // remains the single source of truth for "this CorrType value is not
  // recognized at all".
  MaterialCorrectionModeSupport materialModeSupport = MaterialCorrectionModeSupport::Supported;
  dispatchTransitionPolicies(grouping, [&](auto traits, auto /*transitionIds*/, auto /*cellIds*/) {
    using Traits = decltype(traits);
    materialModeSupport = checkMaterialCorrectionModeSupport<Traits::Tag>(mTrkParams[iteration].CorrType);
  });
  if (materialModeSupport == MaterialCorrectionModeSupport::Unsupported) {
    throw TraversalException{iteration, TraversalFailureReason::UnsupportedMaterialCorrectionMode};
  }

  // 2.5. Resolve and validate this iteration's authoritative per-layer
  // nominal material, entirely from `layouts`/`layout` (already resolved
  // above) and `mTrkParams[iteration]` -- before any TimeFrame tracking
  // state is touched (see TrackerTraits.h's TraversalFailureReason doc).
  // The layer-to-surface mapping is this DetectorLayoutSet's own validated
  // orderedSurfaces (never inferred from legacy index, detector identity,
  // radius, z, or numeric ordering); a size mismatch, an invalid/out-of-range
  // mapped SurfaceId, or a numeric disagreement against the temporary legacy
  // TrackingParameters::LayerxX0 all reject with the same
  // LegacyMaterialMismatch reason -- this whole block is one compatibility
  // precondition, not several. SurfaceDescriptor::material itself is never
  // written here. `stagedLayerMaterial` is kept local (not yet committed to
  // the mLayerMaterial member) until every remaining fallible check in this
  // function has also succeeded -- see the final commit block below -- so a
  // later failure leaves mLayerMaterial exactly as resetTraversalCache() left
  // it (reset/zero-filled), never partially populated from a failed
  // iteration.
  //
  // Gate 4 Slice 0a: the same walk also stages mSurfaceToLegacyLayer, the
  // reverse SurfaceId->legacyLayer bridge the migrated tracklet loop needs
  // (see that member's doc in TrackerTraits.h). `surfaceId.value() <
  // layout.nSurfaces <= MaxLayoutSurfaces` is already established by the
  // check immediately below, so indexing the MaxLayoutSurfaces-wide staged
  // array is always in-bounds. A duplicate SurfaceId across two legacyLayer
  // entries -- orderedSurfaces failing to be a bijection onto this
  // iteration's legacy layers -- rejects with SurfaceLayerMappingMismatch,
  // before any TimeFrame tracking state is touched, exactly like every
  // other check in this block.
  const auto& orderedSurfaces = layouts.getConfigurationKey().orderedSurfaces;
  if (orderedSurfaces.size() < static_cast<size_t>(NLayers)) {
    throw TraversalException{iteration, TraversalFailureReason::LegacyMaterialMismatch};
  }
  std::array<NominalSurfaceMaterial, NLayers> stagedLayerMaterial{};
  std::array<uint8_t, MaxLayoutSurfaces> stagedSurfaceToLegacyLayer;
  stagedSurfaceToLegacyLayer.fill(kInvalidLegacyLayer);
  for (int legacyLayer = 0; legacyLayer < NLayers; ++legacyLayer) {
    const auto surfaceId = orderedSurfaces[legacyLayer];
    if (!surfaceId.isValid() || surfaceId.value() >= layout.nSurfaces) {
      throw TraversalException{iteration, TraversalFailureReason::LegacyMaterialMismatch};
    }
    if (stagedSurfaceToLegacyLayer[surfaceId.value()] != kInvalidLegacyLayer) {
      throw TraversalException{iteration, TraversalFailureReason::SurfaceLayerMappingMismatch};
    }
    stagedSurfaceToLegacyLayer[surfaceId.value()] = static_cast<uint8_t>(legacyLayer);
    stagedLayerMaterial[legacyLayer] = layout.getSurface(surfaceId).material;
  }
  if (mTrkParams[iteration].LayerxX0.size() != static_cast<size_t>(NLayers)) {
    throw TraversalException{iteration, TraversalFailureReason::LegacyMaterialMismatch};
  }
  for (int legacyLayer = 0; legacyLayer < NLayers; ++legacyLayer) {
    if (mTrkParams[iteration].LayerxX0[legacyLayer] != stagedLayerMaterial[legacyLayer].xOverX0) {
      throw TraversalException{iteration, TraversalFailureReason::LegacyMaterialMismatch};
    }
  }

  // 2.6. One-time normalized-measurement binding (Stage-B normalized-CA-
  // measurements slice): resolve and validate this iteration's authoritative
  // per-(legacy-)layer normalized SurfaceMeasurement span, entirely from
  // `orderedSurfaces` (already resolved and validated above) and the
  // already-loaded TimeFrame normalized frame / legacy compatibility
  // structures -- before any TimeFrame tracking state is touched.
  // `stagedLayerMeasurements` is kept local (not yet committed to the
  // mLayerMeasurements member) until every remaining fallible check in this
  // function has also succeeded -- see the final commit block below -- so a
  // later failure leaves mLayerMeasurements exactly as resetTraversalCache()
  // left it (reset/empty spans), never partially populated from a failed
  // iteration.
  LayerMeasurementSpans<NLayers> stagedLayerMeasurements{};
  for (int legacyLayer = 0; legacyLayer < NLayers; ++legacyLayer) {
    const auto surfaceId = orderedSurfaces[legacyLayer];
    const auto measurements = mFrame->getNormalizedFrame().getSurfaceMeasurements(surfaceId);
    const auto& legacyClusters = mScratch->getUnsortedClusters()[legacyLayer];
    const auto& legacyHits = mScratch->getTrackingFrameInfoOnLayer(legacyLayer);
    if (measurements.size() != legacyClusters.size() || legacyHits.size() != legacyClusters.size() ||
        measurements.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
    }
    for (size_t i = 0; i < measurements.size(); ++i) {
      const auto& measurement = measurements[i];
      if (measurement.surface != surfaceId || !measurement.cluster.isValid() ||
          measurement.cluster.source != ClusterSourceId{0} ||
          measurement.cluster.index > static_cast<uint32_t>(std::numeric_limits<int>::max()) ||
          legacyClusters[i].clusterId != static_cast<int>(i)) {
        throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
      }
      const int legacyExternalIndex = mScratch->getClusterExternalIndex(legacyLayer, static_cast<int>(i));
      if (legacyExternalIndex < 0 || static_cast<uint32_t>(legacyExternalIndex) != measurement.cluster.index) {
        throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
      }
    }
    const auto rofBoundaries = mScratch->getROFrameClusters(legacyLayer);
    if (rofBoundaries.empty() || rofBoundaries.front() != 0 || rofBoundaries.back() != static_cast<int>(measurements.size())) {
      throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
    }
    for (size_t rof = 0; rof + 1 < rofBoundaries.size(); ++rof) {
      const int first = rofBoundaries[rof];
      const int last = rofBoundaries[rof + 1];
      if (first < 0 || last < first || last > static_cast<int>(measurements.size())) {
        throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
      }
      for (int index = first; index < last; ++index) {
        if (measurements[index].sourceROF != rof) {
          throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
        }
      }
    }
    stagedLayerMeasurements[legacyLayer] = measurements;
  }

  // 3. Bind + validate index-table configuration into a local scratch value,
  // dispatched on the single active tag via dispatchTransitionPolicies --
  // the same idiom computeLayerTracklets() uses below. The scratch is not
  // touched yet, so a failure here leaves it completely unchanged.
  typename ScratchN::IndexTableUtilsN stagedIndexTableConfig{};
  IndexTableConfigError indexTableConfigError = IndexTableConfigError::None;
  bool activePolicyTagResolved = false;
  bool stateFamilyMismatch = false;
  dispatchTransitionPolicies(grouping, [&](auto traits, auto /*transitionIds*/, auto /*cellIds*/) {
    using Traits = decltype(traits);
    activePolicyTagResolved = true;
    // Discards the call (and therefore the need for a
    // bindIndexTableConfiguration<Traits::Tag, NLayers> instantiation)
    // whenever this policy family's seed state cannot possibly match this
    // TrackerTraits<NLayers> instantiation's own CellSeedN -- the identical
    // compile-time compatibility guard computeLayerTracklets() already uses
    // below. A mismatch here is a genuine misconfiguration (this layout's
    // active transitions do not belong to this NLayers/state family), not a
    // NLayers-to-Tag policy selection: the active Tag itself still comes
    // exclusively from `grouping`/`layout` above.
    if constexpr (stateFamilyFromNLayers<NLayers>() != Traits::Family) {
      stateFamilyMismatch = true;
    } else {
      indexTableConfigError = bindIndexTableConfiguration<Traits::Tag, NLayers>(stagedIndexTableConfig, mTrkParams[iteration]);
    }
  });
  if (!activePolicyTagResolved || stateFamilyMismatch) {
    throw TraversalException{iteration, TraversalFailureReason::StateFamilyMismatch};
  }
  if (indexTableConfigError != IndexTableConfigError::None) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidIndexTableConfiguration};
  }

  // 4. LUT-reuse invariant: a non-FirstPass iteration will not have its
  // index-table configuration (re)committed or its LUT storage reallocated
  // by TimeFrame::initialise() below (TimeFrame.cxx, PassFlags[FirstPass]
  // gate) -- whether or not RebuildClusterLUT resorts clusters into the
  // existing LUT using whatever configuration TimeFrame already owns (e.g.
  // legacy ITS's async iteration 3, which sets RebuildClusterLUT without
  // FirstPass; ITS/tracking/src/Configuration.cxx). The freshly bound
  // configuration for such an iteration must therefore already match the
  // owned one exactly, checked before TimeFrame is touched at all.
  if (!mTrkParams[iteration].PassFlags[IterationStep::FirstPass] &&
      !indexTableConfigurationsMatch(stagedIndexTableConfig, mScratch->getIndexTableUtils())) {
    throw TraversalException{iteration, TraversalFailureReason::IndexTableConfigurationMismatch};
  }

  // 5. Only now is the scratch touched: it receives an already-validated
  // configuration by value and never inspects a tag or detector ID.
  mScratch->initialise(*mFrame, mTrkParams[iteration], mTrkParams[iteration].NLayers, iteration,
                       stagedIndexTableConfig, stagedLayerMeasurements);

  // A sorted Cluster is a locator/navigation cache only. Validate each
  // enabled ROF that can participate in a configured transition after every
  // TimeFrame initialisation, including LUT reuse and non-FirstPass paths.
  // The spans remain local until this check and all subsequent policy setup
  // have succeeded, so a structural failure cannot publish traversal caches.
  std::array<bool, NLayers> candidateReachableLayers{};
  for (int transitionId = 0; transitionId < mScratch->getTrackingTopologyView().nTransitions; ++transitionId) {
    const auto& transition = mScratch->getTrackingTopologyView().getTransition(transitionId);
    candidateReachableLayers[transition.fromLayer] = true;
    candidateReachableLayers[transition.toLayer] = true;
  }
  for (int layer = 0; layer < NLayers; ++layer) {
    if (!candidateReachableLayers[layer]) {
      continue;
    }
    const auto measurements = stagedLayerMeasurements[layer];
    const auto rofBoundaries = mScratch->getROFrameClusters(layer);
    const auto rofMask = mScratch->getROFMaskView();
    // Orchestration-only users can intentionally omit the mask altogether;
    // without it no ROF is candidate-reachable. Do not validate allocated
    // spans until a later configuration actually enables one.
    if (rofMask.mFlatMask == nullptr || rofMask.mLayerROFOffsets == nullptr) {
      continue;
    }
    for (int rof = 0; rof < mScratch->getNrof(layer); ++rof) {
      const auto sorted = mScratch->getClustersOnLayer(rof, layer);
      if (sorted.empty()) {
        continue;
      }
      if (!rofMask.isROFEnabled(layer, rof)) {
        continue;
      }
      const int first = rofBoundaries[rof];
      const int last = rofBoundaries[rof + 1];
      if (first < 0 || last < first || last > static_cast<int>(measurements.size()) ||
          sorted.size() != static_cast<size_t>(last - first)) {
        throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
      }
      std::vector<uint8_t> seen(static_cast<size_t>(last - first), uint8_t{0});
      for (const auto& locator : sorted) {
        const int clusterId = locator.clusterId;
        if (clusterId < first || clusterId >= last) {
          throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
        }
        const size_t localId = static_cast<size_t>(clusterId - first);
        if (seen[localId] != 0) {
          throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
        }
        seen[localId] = 1;
        const auto& measurement = measurements[clusterId];
        const int externalIndex = mScratch->getClusterExternalIndex(layer, clusterId);
        if (measurement.surface != orderedSurfaces[layer] || measurement.sourceROF != static_cast<uint32_t>(rof) || !measurement.cluster.isValid() ||
            measurement.cluster.source != ClusterSourceId{0} || externalIndex < 0 ||
            static_cast<uint32_t>(externalIndex) != measurement.cluster.index) {
          throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
        }
      }
      if (std::find(seen.begin(), seen.end(), uint8_t{0}) != seen.end()) {
        throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
      }
    }
  }

  // 6. validateLegacyParity needs mScratch->getTrackingTopologyView(),
  // which TimeFrame::initialise() just populated, so it necessarily still
  // runs after that call -- exactly as before, just with index-table
  // binding now happening ahead of it instead of buried inside it. Both
  // this resolution and step 2/3's are pure functions of the identical
  // `layout` value, so they agree by construction.
  TransitionPolicyTag activeTag = TransitionPolicyTag::Invalid;
  bool mixedPolicy = false;
  validateLegacyParity(iteration, layout, activeTag, mixedPolicy);
  if (mixedPolicy) {
    throw TraversalException{iteration, TraversalFailureReason::MixedPolicyLayout};
  }

  constexpr StateFamily cellStateFamily = stateFamilyFromNLayers<NLayers>();
  if (stateFamilyOf(activeTag) != cellStateFamily) {
    throw TraversalException{iteration, TraversalFailureReason::StateFamilyMismatch};
  }

  std::optional<CylinderCylinderPolicyParams> cylinderParams;
  std::optional<DiskDiskPolicyParams> diskParams;
  // Bound from the still-local stagedLayerMaterial, not the mLayerMaterial
  // member (not committed yet): attachHitConfig.layerMaterial is rebound to
  // point at mLayerMaterial once that member is actually populated, at the
  // final commit below, before it escapes into mAttachHitConfig.
  auto attachHitConfig = bindAttachHitPolicyConfig(
    gsl::span<const NominalSurfaceMaterial>(stagedLayerMaterial.data(), stagedLayerMaterial.size()), mTrkParams[iteration]);
  if (!attachHitConfig.isValid(NLayers)) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidPolicyParameters};
  }
  const auto geometryConfig = bindLayerGeometryConfig(mTrkParams[iteration], attachHitConfig);
  if (!geometryConfig.isValid(NLayers)) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidPolicyParameters};
  }
  DiskDiskReferenceCoordinateView referenceCoordinateView{};
  if (activeTag == TransitionPolicyTag::DiskDisk) {
    referenceCoordinateView = bindLegacyMFTReferenceCoordinates();
    if (!referenceCoordinateView.isValid(NLayers)) {
      throw TraversalException{iteration, TraversalFailureReason::InvalidPolicyParameters};
    }
  }
  if (activeTag == TransitionPolicyTag::CylinderCylinder) {
    cylinderParams = bindTransitionPolicyParams<TransitionPolicyTag::CylinderCylinder>(mTrkParams[iteration]);
    ++mPolicyBindingCounts[0];
    if (!cylinderParams->isValid()) {
      throw TraversalException{iteration, TraversalFailureReason::InvalidPolicyParameters};
    }
  } else if (activeTag == TransitionPolicyTag::DiskDisk) {
    diskParams = bindTransitionPolicyParams<TransitionPolicyTag::DiskDisk>(mTrkParams[iteration]);
    ++mPolicyBindingCounts[1];
    if (!diskParams->isValid()) {
      throw TraversalException{iteration, TraversalFailureReason::InvalidPolicyParameters};
    }
  } else {
    throw TraversalException{iteration, TraversalFailureReason::StateFamilyMismatch};
  }

  mTraversalLayout = layout;
  mTraversalGrouping.emplace(std::move(grouping));
  mCylinderPolicyParams = cylinderParams;
  mDiskPolicyParams = diskParams;
  mDiskLayerReferenceZ = referenceCoordinateView.perLayerReferenceZ;
  // Commit the material cache itself only now, alongside every other
  // traversal cache, then rebind attachHitConfig's span off the
  // about-to-be-destroyed local stagedLayerMaterial and onto the
  // now-populated, function-lifetime-independent mLayerMaterial member
  // before mAttachHitConfig (which outlives this call) retains it.
  mLayerMaterial = stagedLayerMaterial;
  // Gate 4 Slice 0a: committed alongside mLayerMaterial, from the same
  // orderedSurfaces walk above -- see mSurfaceToLegacyLayer's own doc.
  mSurfaceToLegacyLayer = stagedSurfaceToLegacyLayer;
  attachHitConfig.layerMaterial = gsl::span<const NominalSurfaceMaterial>(mLayerMaterial.data(), mLayerMaterial.size());
  mAttachHitConfig = attachHitConfig;
  // One-time normalized-measurement binding: committed here, alongside every
  // other successfully staged traversal cache -- see mLayerMeasurements' own
  // doc for the commit contract.
  mLayerMeasurements = stagedLayerMeasurements;

  // All fallible validation for this iteration (layout/grouping, legacy
  // parity, state-family, and every policy/geometry binding above) has now
  // succeeded. What follows is the relocated, total (non-throwing) per-layer
  // and per-transition scattering/bending preparation -- see the method doc
  // for the ordering/failure contract this relies on.
  if (activeTag == TransitionPolicyTag::CylinderCylinder) {
    prepareTransitionScatteringAndBendingForPolicy<TransitionPolicyTag::CylinderCylinder>(iteration, geometryConfig, referenceCoordinateView);
  } else {
    prepareTransitionScatteringAndBendingForPolicy<TransitionPolicyTag::DiskDisk>(iteration, geometryConfig, referenceCoordinateView);
  }
}

template <int NLayers>
template <TransitionPolicyTag Tag>
void TrackerTraits<NLayers>::prepareTransitionScatteringAndBendingForPolicy(
  int iteration,
  const LayerGeometryConfigView& geometryConfig,
  const DiskDiskReferenceCoordinateView& referenceCoordinateView)
{
  const auto& trkParam = mTrkParams[iteration];

  // Per-layer step: genuinely policy-specific (typed operation), preserving
  // the exact legacy loop bound (compile-time NLayers, not trkParam.NLayers --
  // matches TimeFrame.cxx's original "estimate MS per layer" loop verbatim).
  std::array<float, NLayers> msAngles{};
  for (unsigned int iLayer{0}; iLayer < NLayers; ++iLayer) {
    if constexpr (Tag == TransitionPolicyTag::CylinderCylinder) {
      msAngles[iLayer] = layerMultipleScatteringAngle<Tag>(
        LayerScatteringInputs<Tag>{geometryConfig.layerMaterial[iLayer].xOverX0}, trkParam.TrackletMinPt);
    } else {
      msAngles[iLayer] = layerMultipleScatteringAngle<Tag>(
        LayerScatteringInputs<Tag>{geometryConfig.layerMaterial[iLayer].xOverX0, geometryConfig.layerRadii[iLayer],
                                   referenceCoordinateView.perLayerReferenceZ[iLayer]},
        trkParam.TrackletMinPt);
    }
  }

  // Per-transition step: shared/Tag-independent post-clamp arithmetic behind
  // a Tag-specific curvature clamp. Iterates legacy transitionIds
  // 0..nTransitions-1 directly off the legacy-shaped TrackingTopology view (in
  // increasing order, matching the frozen code) rather than through
  // mTraversalGrouping's per-tag span, so the loop-carried oneOverR ratchet
  // is threaded in exactly the same order every time -- ordering does not
  // depend on, and is not proven against, grouping-span order.
  const auto& topology = mScratch->getTrackingTopologyView();
  auto& transitionMSAngles = mScratch->getTransitionMSAngles();
  auto& transitionPhiCuts = mScratch->getTransitionPhiCuts();
  float oneOverR{0.001f * 0.3f * std::abs(getBz()) / trkParam.TrackletMinPt};
  for (int transitionId{0}; transitionId < static_cast<int>(topology.nTransitions); ++transitionId) {
    const auto& transition = topology.getTransition(transitionId);
    const float r1 = trkParam.LayerRadii[transition.fromLayer];
    const float r2 = trkParam.LayerRadii[transition.toLayer];
    oneOverR = clampTransitionCurvature<Tag>(oneOverR, r2);
    const float res1 = o2::gpu::CAMath::Hypot(trkParam.PVres, mScratch->getPositionResolution(transition.fromLayer));
    const float res2 = o2::gpu::CAMath::Hypot(trkParam.PVres, mScratch->getPositionResolution(transition.toLayer));
    const auto prep = prepareTransitionScatteringAndBending(
      gsl::span<const float>(msAngles.data(), msAngles.size()), transition.fromLayer, transition.toLayer, r1, r2, oneOverR, res1, res2);
    transitionMSAngles[transitionId] = prep.msAngle;
    transitionPhiCuts[transitionId] = prep.phiCut;
  }
}

template <int NLayers>
void TrackerTraits<NLayers>::computeLayerTracklets(const int iteration, int iVertex)
{
  // Gate 4 Slice 0a: driven by the sparse topology cached on
  // initialiseTimeFrame() (mTraversalLayout), not a fresh legacy
  // getTrackingTopologyView() fetch. This clear/allocation pass is
  // tag-agnostic and runs once over the full sparse transition count,
  // regardless of which single policy tag is active for this layout (see
  // computeLayerTrackletsForPolicy() below for the per-tag-filtered body).
  const auto& topology = mTraversalLayout.topology;
  for (uint32_t transitionId = 0; transitionId < topology.nTransitions; ++transitionId) {
    mScratch->getTracklets()[transitionId].clear();
    mScratch->getTrackletsLabel(transitionId).clear();
    std::fill(mScratch->getTrackletsLookupTable()[transitionId].begin(), mScratch->getTrackletsLookupTable()[transitionId].end(), 0);
  }

  if (!mTraversalGrouping.has_value()) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidTraversalSchedule};
  }

  dispatchTransitionPolicies(*mTraversalGrouping, [&](auto traits, auto transitionIds, auto) {
    using Traits = decltype(traits);
    if constexpr (stateFamilyFromNLayers<NLayers>() != Traits::Family) {
      throw TraversalException{iteration, TraversalFailureReason::StateFamilyMismatch};
    } else if constexpr (Traits::Tag == TransitionPolicyTag::CylinderCylinder) {
      if (!mCylinderPolicyParams.has_value()) {
        throw TraversalException{iteration, TraversalFailureReason::InvalidPolicyParameters};
      }
      computeLayerTrackletsForPolicy<Traits::Tag>(iteration, iVertex, transitionIds, *mCylinderPolicyParams);
    } else if constexpr (Traits::Tag == TransitionPolicyTag::DiskDisk) {
      if (!mDiskPolicyParams.has_value()) {
        throw TraversalException{iteration, TraversalFailureReason::InvalidPolicyParameters};
      }
      computeLayerTrackletsForPolicy<Traits::Tag>(iteration, iVertex, transitionIds, *mDiskPolicyParams);
    }
  });
}

template <int NLayers>
template <TransitionPolicyTag Tag>
void TrackerTraits<NLayers>::computeLayerTrackletsForPolicy(
  const int iteration,
  const int iVertex,
  gsl::span<const TransitionId> transitionIds,
  const typename TransitionPolicyTraits<Tag>::Params& params)
{
  // Gate 4 Slice 0a: sparse topology (cached, not re-fetched from the
  // legacy view). `transitionIds` remains the caller-filtered, ascending
  // per-tag span from TransitionPolicyGrouping -- every loop below still
  // iterates exclusively over it (or over the tracklet storage it already
  // populated), never over the raw sparse transition count.
  const auto& topology = mTraversalLayout.topology;
  const Vertex diamondVert(mTrkParams[iteration].Diamond, mTrkParams[iteration].DiamondCov, 1, 1.f);

  mTaskArena->execute([&] {
    // Resolves one sparse SurfaceTransition's endpoints to this
    // TrackerTraits<NLayers>'s legacy layer indices through
    // mSurfaceToLegacyLayer (built/validated once per iteration in
    // initialiseTimeFrame(), see that member's doc). Called exactly once
    // per transitionId, outside every candidate (ROF/cluster/vertex) loop
    // below -- never re-derived per candidate.
    auto resolveTransitionLayers = [&](int transitionId) -> std::pair<int, int> {
      const auto& transition = topology.getTransition(TransitionId{static_cast<uint16_t>(transitionId)});
      return {static_cast<int>(mSurfaceToLegacyLayer[transition.from.value()]),
              static_cast<int>(mSurfaceToLegacyLayer[transition.to.value()])};
    };

    auto makeTransitionState = [&](int transitionId, int fromLayer, int toLayer) {
      if constexpr (Tag == TransitionPolicyTag::CylinderCylinder) {
        return TrackletProjectionState<Tag>{fromLayer,
                                            toLayer,
                                            mTrkParams[iteration].LayerRadii[toLayer] - mTrkParams[iteration].LayerRadii[fromLayer],
                                            mScratch->getMinR(toLayer),
                                            mScratch->getMaxR(toLayer),
                                            mScratch->getPositionResolution(fromLayer),
                                            mScratch->getTransitionMSAngle(transitionId),
                                            mScratch->getTransitionPhiCut(transitionId)};
      } else {
        const float fromZ = detail::mftLayerZ(fromLayer);
        const float toZ = detail::mftLayerZ(toLayer);
        return TrackletProjectionState<Tag>{fromLayer,
                                            toLayer,
                                            fromZ,
                                            toZ,
                                            toZ - fromZ,
                                            mTrkParams[iteration].LayerRadii[fromLayer],
                                            mScratch->getTransitionMSAngle(transitionId),
                                            mScratch->getTransitionPhiCut(transitionId)};
      }
    };

    auto forTracklets = [&](auto Mode, int transitionId, int fromLayer, int toLayer, const TrackletProjectionState<Tag>& transitionState, int pivotROF, int base, int& offset) -> int {
      if (!mScratch->getROFMaskView().isROFEnabled(fromLayer, pivotROF)) {
        return 0;
      }
      // A diamond vertex valid for this specific pivotROF is derived fresh
      // here (see diamondVertexForROF's doc comment above) rather than
      // reused from a single shared instance: the derivation is cheap
      // (two field reads + arithmetic already computed for every real ROF
      // timing check) and each forTracklets invocation owns its own stack
      // frame, so this is safe under the tbb::parallel_for dispatch below.
      Vertex diamondForROF{};
      gsl::span<const Vertex> primaryVertices;
      if (mTrkParams[iteration].UseDiamond) {
        diamondForROF = diamondVertexForROF(diamondVert, mScratch->getROFOverlapTableView(), fromLayer, pivotROF);
        primaryVertices = gsl::span<const Vertex>(&diamondForROF, 1);
      } else {
        primaryVertices = mScratch->getPrimaryVertices(*mFrame, fromLayer, pivotROF);
      }
      if (primaryVertices.empty()) {
        return 0;
      }
      const int startVtx = iVertex >= 0 ? iVertex : 0;
      const int endVtx = iVertex >= 0 ? o2::gpu::CAMath::Min(iVertex + 1, int(primaryVertices.size())) : int(primaryVertices.size());
      if (endVtx <= startVtx || (iVertex + 1) > primaryVertices.size()) {
        return 0;
      }

      const auto& rofOverlap = mScratch->getROFOverlapTableView().getOverlap(fromLayer, toLayer, pivotROF);
      if (!rofOverlap.getEntries()) {
        return 0;
      }

      int localCount = 0;
      auto& tracklets = mScratch->getTracklets()[transitionId];
      auto layer0 = mScratch->getClustersOnLayer(pivotROF, fromLayer);
      if (layer0.empty()) {
        return 0;
      }

      for (int iCluster = 0; iCluster < int(layer0.size()); ++iCluster) {
        const Cluster& currentCluster = layer0[iCluster];
        const int currentSortedIndex = mScratch->getSortedIndex(pivotROF, fromLayer, iCluster);
        if (mScratch->isClusterUsed(fromLayer, currentCluster.clusterId)) {
          continue;
        }
        const auto& sourceMeasurement = mLayerMeasurements[fromLayer][currentCluster.clusterId];

        for (int iV = startVtx; iV < endVtx; ++iV) {
          const auto& pv = primaryVertices[iV];
          if (!mScratch->getROFVertexLookupTableView().isVertexCompatible(fromLayer, pivotROF, pv)) {
            continue;
          }
          if (pv.isFlagSet(Vertex::Flags::UPCMode) != mTrkParams[iteration].PassFlags[IterationStep::SelectUPCVertices]) {
            continue;
          }
          TrackletSearchWindow<Tag> searchWindow{};
          if (!projectSearchWindow<Tag, NLayers>(sourceMeasurement, currentCluster, pv, transitionState, getBz(),
                                                 mScratch->getIndexTableUtils(), params, searchWindow)) {
            continue;
          }
          const auto bins = searchWindow.bins;
          const int rowBinsCount = mScratch->getIndexTableUtils().getNrowBins();
          int rowBinsNum = bins.w - bins.y + 1;
          if (rowBinsNum < 0) {
            rowBinsNum += rowBinsCount;
          }

          for (int targetROF = rofOverlap.getFirstEntry(); targetROF < rofOverlap.getEntriesBound(); ++targetROF) {
            if (!mScratch->getROFMaskView().isROFEnabled(toLayer, targetROF)) {
              continue;
            }
            auto layer1 = mScratch->getClustersOnLayer(targetROF, toLayer);
            if (layer1.empty()) {
              continue;
            }
            const auto ts = mScratch->getROFOverlapTableView().getTimeStamp(fromLayer, pivotROF, toLayer, targetROF);
            if (!ts.isCompatible(pv.getTimeStamp())) {
              continue;
            }
            const auto& targetIndexTable = mScratch->getIndexTable(targetROF, toLayer);
            const int colBinRange = (bins.z - bins.x) + 1;
            for (int iRow = 0; iRow < rowBinsNum; ++iRow) {
              int iRowBin = bins.y + iRow;
              if constexpr (Tag == TransitionPolicyTag::DiskDisk) {
                if (iRowBin >= rowBinsCount) {
                  break;
                }
              } else {
                iRowBin %= rowBinsCount;
              }
              const int firstBinIdx = mScratch->getIndexTableUtils().getBinIndex(bins.x, iRowBin);
              const int maxBinIdx = firstBinIdx + colBinRange;
              const int firstRow = targetIndexTable[firstBinIdx];
              const int lastRow = targetIndexTable[maxBinIdx];
              for (int iNext = firstRow; iNext < lastRow; ++iNext) {
                if (iNext >= int(layer1.size())) {
                  break;
                }
                const Cluster& nextCluster = layer1[iNext];
                if (mScratch->isClusterUsed(toLayer, nextCluster.clusterId)) {
                  continue;
                }
                const auto& targetMeasurement = mLayerMeasurements[toLayer][nextCluster.clusterId];

                float tanL = 0.f;
                const bool accepted = [&] {
                  if constexpr (Tag == TransitionPolicyTag::CylinderCylinder) {
                    return searchWindow.acceptCandidate(sourceMeasurement, currentCluster, targetMeasurement, nextCluster, tanL);
                  } else {
                    return searchWindow.acceptCandidate(sourceMeasurement, targetMeasurement, tanL);
                  }
                }();
                if (accepted) {
                  const float phi{o2::gpu::GPUCommonMath::ATan2(sourceMeasurement.global.y - targetMeasurement.global.y,
                                                                sourceMeasurement.global.x - targetMeasurement.global.x)};
                  if constexpr (decltype(Mode)::value == PassMode::OnePass::value) {
                    tracklets.emplace_back(currentSortedIndex, mScratch->getSortedIndex(targetROF, toLayer, iNext), tanL, phi, ts);
                  } else if constexpr (decltype(Mode)::value == PassMode::TwoPassCount::value) {
                    ++localCount;
                  } else if constexpr (decltype(Mode)::value == PassMode::TwoPassInsert::value) {
                    const int idx = base + offset++;
                    tracklets[idx] = Tracklet(currentSortedIndex, mScratch->getSortedIndex(targetROF, toLayer, iNext), tanL, phi, ts);
                  }
                }
              }
            }
          }
        }
      }
      return localCount;
    };

    int dummy{0};
    if (mTaskArena->max_concurrency() <= 1) {
      for (const auto typedTransitionId : transitionIds) {
        const int transitionId = typedTransitionId.value();
        const auto [fromLayer, toLayer] = resolveTransitionLayers(transitionId);
        const auto transitionState = makeTransitionState(transitionId, fromLayer, toLayer);
        const int startROF = 0, endROF = mScratch->getROFOverlapTableView().getLayer(fromLayer).mNROFsTF;
        for (int pivotROF{startROF}; pivotROF < endROF; ++pivotROF) {
          forTracklets(PassMode::OnePass{}, transitionId, fromLayer, toLayer, transitionState, pivotROF, 0, dummy);
        }
      }
    } else {
      tbb::parallel_for(0, static_cast<int>(transitionIds.size()), [&](const int transitionIndex) {
        const int transitionId = transitionIds[transitionIndex].value();
        const auto [fromLayer, toLayer] = resolveTransitionLayers(transitionId);
        const auto transitionState = makeTransitionState(transitionId, fromLayer, toLayer);
        const int startROF = 0, endROF = mScratch->getROFOverlapTableView().getLayer(fromLayer).mNROFsTF;
        bounded_vector<int> perROFCount((endROF - startROF) + 1, mMemoryPool.get());
        tbb::parallel_for(startROF, endROF, [&](const int pivotROF) {
          perROFCount[pivotROF - startROF] = forTracklets(PassMode::TwoPassCount{}, transitionId, fromLayer, toLayer, transitionState, pivotROF, 0, dummy);
        });
        std::exclusive_scan(perROFCount.begin(), perROFCount.end(), perROFCount.begin(), 0);
        const int nTracklets = perROFCount.back();
        mScratch->getTracklets()[transitionId].resize(nTracklets);
        if (nTracklets == 0) {
          return;
        }
        tbb::parallel_for(startROF, endROF, [&](const int pivotROF) {
          int baseIdx = perROFCount[pivotROF - startROF];
          if (baseIdx == perROFCount[pivotROF + 1 - startROF]) {
            return;
          }
          int localIdx = 0;
          forTracklets(PassMode::TwoPassInsert{}, transitionId, fromLayer, toLayer, transitionState, pivotROF, baseIdx, localIdx);
        });
      });
    }

    tbb::parallel_for(0, static_cast<int>(transitionIds.size()), [&](const int transitionIndex) {
      const int transitionId = transitionIds[transitionIndex].value();
      /// Sort tracklets & remove duplicates
      // duplicates can exist simply since we evaluate per vertex
      auto& trkl{mScratch->getTracklets()[transitionId]};
      std::sort(trkl.begin(), trkl.end());
      trkl.erase(std::unique(trkl.begin(), trkl.end()), trkl.end());
      trkl.shrink_to_fit();
      auto& lut{mScratch->getTrackletsLookupTable()[transitionId]};
      if (!trkl.empty()) {
        for (const auto& tkl : trkl) {
          lut[tkl.firstClusterIndex + 1]++;
        }
        std::inclusive_scan(lut.begin(), lut.end(), lut.begin());
      }
    });

    /// Create tracklets labels
    if (mScratch->hasMCinformation() && mTrkParams[iteration].CreateArtefactLabels) {
      tbb::parallel_for(0, static_cast<int>(transitionIds.size()), [&](const int transitionIndex) {
        const int transitionId = transitionIds[transitionIndex].value();
        const auto [fromLayer, toLayer] = resolveTransitionLayers(transitionId);
        for (auto& trk : mScratch->getTracklets()[transitionId]) {
          MCCompLabel label;
          int currentId{mScratch->getClusters()[fromLayer][trk.firstClusterIndex].clusterId};
          int nextId{mScratch->getClusters()[toLayer][trk.secondClusterIndex].clusterId};
          for (const auto& lab1 : mScratch->getClusterLabels(fromLayer, currentId)) {
            for (const auto& lab2 : mScratch->getClusterLabels(toLayer, nextId)) {
              if (lab1 == lab2 && lab1.isValid()) {
                label = lab1;
                break;
              }
            }
            if (label.isValid()) {
              break;
            }
          }
          mScratch->getTrackletsLabel(transitionId).emplace_back(label);
        }
      });
    }
  });
}

template <int NLayers>
void TrackerTraits<NLayers>::computeLayerCells(const int iteration)
{
  // Gate 4 Slice 0b: driven by the sparse topology cached on
  // initialiseTimeFrame() (mTraversalLayout), not a fresh legacy
  // getTrackingTopologyView() fetch.
  const auto& topology = mTraversalLayout.topology;
  // Defensive size-consistency check, mirroring findCellsNeighboursForPolicy's
  // own precedent: topology.nCells is no longer definitionally
  // mScratch->getCells().size() once the loop bound below is sparse-derived,
  // so check it explicitly, fail-closed, before the clear loop touches either
  // container -- this is the only point that can actually intercept the
  // hazard, since the clear loop itself indexes both unconditionally.
  if (mScratch->getCells().size() != topology.nCells ||
      mScratch->getCellsLookupTable().size() != topology.nCells) {
    throw TraversalException{iteration, TraversalFailureReason::LegacyIndexMismatch};
  }
  // This clear/allocation pass is tag-agnostic and runs once over the full
  // sparse cell/transition count, regardless of which single policy tag is
  // active for this layout (see computeLayerCellsForPolicy() below for the
  // per-tag-filtered body).
  for (uint32_t cellTopologyId = 0; cellTopologyId < topology.nCells; ++cellTopologyId) {
    deepVectorClear(mScratch->getCells()[cellTopologyId]);
    deepVectorClear(mScratch->getCellsLookupTable()[cellTopologyId]);
    if (mScratch->hasMCinformation() && mTrkParams[iteration].CreateArtefactLabels) {
      deepVectorClear(mScratch->getCellsLabel(cellTopologyId));
    }
  }

  if (!mTraversalGrouping.has_value()) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidTraversalSchedule};
  }
  if (!mAttachHitConfig.isValid(NLayers)) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidPolicyParameters};
  }

  dispatchTransitionPolicies(*mTraversalGrouping, [&](auto traits, auto, auto cellIds) {
    using Traits = decltype(traits);
    if constexpr (stateFamilyFromNLayers<NLayers>() != Traits::Family) {
      throw TraversalException{iteration, TraversalFailureReason::StateFamilyMismatch};
    } else if constexpr (Traits::Tag == TransitionPolicyTag::CylinderCylinder) {
      if (!mCylinderPolicyParams.has_value()) {
        throw TraversalException{iteration, TraversalFailureReason::InvalidPolicyParameters};
      }
      computeLayerCellsForPolicy<Traits::Tag>(iteration, cellIds, *mCylinderPolicyParams);
    } else if constexpr (Traits::Tag == TransitionPolicyTag::DiskDisk) {
      // The size check is a defensive invariant, not independently reachable
      // through the public API today: mDiskLayerReferenceZ and
      // mDiskPolicyParams are always committed together, at the same point,
      // gated by the same activeTag validation in initialiseTimeFrame() (see
      // that method). It guards against a future refactor accidentally
      // decoupling the two commits, not a state this call site can currently
      // observe on its own.
      if (!mDiskPolicyParams.has_value() || mDiskLayerReferenceZ.size() < static_cast<size_t>(NLayers)) {
        throw TraversalException{iteration, TraversalFailureReason::InvalidPolicyParameters};
      }
      computeLayerCellsForPolicy<Traits::Tag>(iteration, cellIds, *mDiskPolicyParams);
    }
  });

  for (uint32_t transitionId = 0; transitionId < topology.nTransitions; ++transitionId) {
    deepVectorClear(mScratch->getTracklets()[transitionId]);
    deepVectorClear(mScratch->getTrackletsLabel(transitionId));
  }
}

template <int NLayers>
template <TransitionPolicyTag Tag>
void TrackerTraits<NLayers>::computeLayerCellsForPolicy(
  const int iteration,
  gsl::span<const CellTopologyId> cellIds,
  const typename TransitionPolicyTraits<Tag>::Params& params)
{
  // Gate 4 Slice 0b: sparse topology (cached, not re-fetched from the
  // legacy view). `cellIds` remains the caller-filtered, ascending per-tag
  // span from TransitionPolicyGrouping -- the main loop below iterates
  // exclusively over it, never over the raw sparse cell count. The
  // defensive size-consistency check for mScratch's per-cellTopologyId
  // storage lives in computeLayerCells() (the caller), not here -- see that
  // function's own doc: its clear loop is the first thing to touch those
  // containers, so that is the only point that can actually intercept a
  // size mismatch before it becomes undefined behaviour.
  const auto& topology = mTraversalLayout.topology;

  mTaskArena->execute([&] {
    // Resolves one sparse SurfaceCellTopology's three hit surfaces to this
    // TrackerTraits<NLayers>'s legacy layer indices through
    // mSurfaceToLegacyLayer (built/validated once per iteration in
    // initialiseTimeFrame(), see that member's doc). Called exactly once
    // per CellTopologyId, outside the per-tracklet candidate loop below --
    // never re-derived per candidate. Every endpoint SurfaceId is checked
    // valid and in-range before it is ever used to index
    // mSurfaceToLegacyLayer, and the mapped entry itself is checked against
    // the invalid sentinel before use as a legacy-layer index -- fails
    // closed with the same LegacyIndexMismatch reason
    // findCellsNeighboursForPolicy already uses for its own analogous
    // sparse-derived-index defensive checks.
    auto resolveCellHitLayers = [&](const auto& cellTopology) -> std::array<int, 3> {
      const auto& firstTransition = topology.getTransition(cellTopology.firstTransition);
      const auto& secondTransition = topology.getTransition(cellTopology.secondTransition);
      const std::array<SurfaceId, 3> surfaces{firstTransition.from, firstTransition.to, secondTransition.to};
      std::array<int, 3> layers{};
      for (int i = 0; i < 3; ++i) {
        const auto surfaceId = surfaces[i];
        if (!surfaceId.isValid() || surfaceId.value() >= MaxLayoutSurfaces) {
          throw TraversalException{iteration, TraversalFailureReason::LegacyIndexMismatch};
        }
        const auto mapped = mSurfaceToLegacyLayer[surfaceId.value()];
        if (mapped == kInvalidLegacyLayer) {
          throw TraversalException{iteration, TraversalFailureReason::LegacyIndexMismatch};
        }
        layers[i] = static_cast<int>(mapped);
      }
      return layers;
    };

    auto forTrackletCells = [&](auto Mode, int firstTransitionId, int secondTransitionId, const std::array<int, 3>& hitLayers, bounded_vector<CellSeedN>& layerCells, int iTracklet, int offset = 0) -> int {
      const Tracklet& currentTracklet{mScratch->getTracklets()[firstTransitionId][iTracklet]};
      const int nextLayerClusterIndex{currentTracklet.secondClusterIndex};
      const int nextLayerFirstTrackletIndex{mScratch->getTrackletsLookupTable()[secondTransitionId][nextLayerClusterIndex]};
      const int nextLayerLastTrackletIndex{mScratch->getTrackletsLookupTable()[secondTransitionId][nextLayerClusterIndex + 1]};
      int foundCells{0};
      for (int iNextTracklet{nextLayerFirstTrackletIndex}; iNextTracklet < nextLayerLastTrackletIndex; ++iNextTracklet) {
        const Tracklet& nextTracklet{mScratch->getTracklets()[secondTransitionId][iNextTracklet]};
        if (nextTracklet.firstClusterIndex != nextLayerClusterIndex) {
          break;
        }
        if (!currentTracklet.getTimeStamp().isCompatible(nextTracklet.getTimeStamp())) {
          continue;
        }

        const float deltaTanLambdaSigma = std::abs(currentTracklet.tanLambda - nextTracklet.tanLambda) / mTrkParams[iteration].CellDeltaTanLambdaSigma;
        if (deltaTanLambdaSigma < mTrkParams[iteration].NSigmaCut) {

          /// Track seed preparation. Clusters are numbered progressively from the innermost going outward.
          const int clusId[3]{
            mScratch->getClusters()[hitLayers[0]][currentTracklet.firstClusterIndex].clusterId,
            mScratch->getClusters()[hitLayers[1]][nextTracklet.firstClusterIndex].clusterId,
            mScratch->getClusters()[hitLayers[2]][nextTracklet.secondClusterIndex].clusterId};

          const auto& measurementInner = mLayerMeasurements[hitLayers[0]][clusId[0]];
          const auto& measurementMiddle = mLayerMeasurements[hitLayers[1]][clusId[1]];
          const auto& measurementOuter = mLayerMeasurements[hitLayers[2]][clusId[2]];

          // MFT geometric road pre-cut: TrackerTraits-owned, outside buildCellSeed
          // (Architecture.md Sec 10 / TransitionPolicyOperations.h doc on buildCellSeed).
          // One unconditional call for both families -- no detector-ID/Tag
          // branch here; CylinderCylinder's specialization is an inline
          // no-op returning true.
          const GlobalPoint3F pointInner = measurementInner.global;
          const GlobalPoint3F pointMiddle = measurementMiddle.global;
          const GlobalPoint3F pointOuter = measurementOuter.global;
          if (!passesCellRoadPrecut<Tag>(pointInner, pointMiddle, pointOuter,
                                         hitLayers[0], hitLayers[1], hitLayers[2],
                                         mDiskLayerReferenceZ, params)) {
            continue;
          }

          // Strictly {inner, middle, outer}: CylinderCylinder reads [1] then
          // [0] (outer slot unused), DiskDisk reads [2], [1], [0].
          const std::array<NominalSurfaceMaterial, 3> material{
            mAttachHitConfig.layerMaterial[hitLayers[0]],
            mAttachHitConfig.layerMaterial[hitLayers[1]],
            mAttachHitConfig.layerMaterial[hitLayers[2]]};

          SurfaceKinematicState state{};
          float chi2{0.f};
          OperationFailureReason buildReason{};
          const bool good = o2::itsmft::tracking::buildCellSeed<Tag>(
            measurementInner, measurementMiddle, measurementOuter,
            material, getBz(), kCompatibilityAbsCharge, kCompatibilityPID, state, chi2, params, buildReason);

          if (good) {
            TimeEstBC ts = currentTracklet.getTimeStamp();
            ts += nextTracklet.getTimeStamp();
            // Gate 4 Slice 0b: reconstructed directly from the three
            // already-resolved legacy layer indices via LayerMask's
            // existing 3-int constructor -- mirrors legacy
            // TrackingTopology::init()'s own hitLayerMask construction
            // (Mask{first.fromLayer, first.toLayer, second.toLayer}), and
            // is provably identical to it under validateLegacyParity's
            // existing hitSurfaces==hitLayerMask cross-check. No generic
            // SurfaceMask->LayerMask converter is introduced.
            const LayerMask hitLayerMask{hitLayers[0], hitLayers[1], hitLayers[2]};
            if constexpr (decltype(Mode)::value == PassMode::OnePass::value) {
              layerCells.emplace_back(hitLayerMask, clusId[0], clusId[1], clusId[2], iTracklet, iNextTracklet, state, chi2, ts);
              ++foundCells;
            } else if constexpr (decltype(Mode)::value == PassMode::TwoPassCount::value) {
              ++foundCells;
            } else if constexpr (decltype(Mode)::value == PassMode::TwoPassInsert::value) {
              layerCells[offset++] = CellSeedN(hitLayerMask, clusId[0], clusId[1], clusId[2], iTracklet, iNextTracklet, state, chi2, ts);
              ++foundCells;
            } else {
              static_assert(false, "Unknown mode!");
            }
          }
        }
      }
      return foundCells;
    };

    for (const auto typedCellId : cellIds) {
      const int cellTopologyId = typedCellId.value();
      const auto& cellTopology = topology.getCell(typedCellId);
      const int firstTransitionId = cellTopology.firstTransition.value();
      const int secondTransitionId = cellTopology.secondTransition.value();
      if (mScratch->getTracklets()[firstTransitionId].empty() ||
          mScratch->getTracklets()[secondTransitionId].empty()) {
        continue;
      }
      const auto hitLayers = resolveCellHitLayers(cellTopology);

      auto& layerCells = mScratch->getCells()[cellTopologyId];
      const int currentLayerTrackletsNum{static_cast<int>(mScratch->getTracklets()[firstTransitionId].size())};
      bounded_vector<int> perTrackletCount(currentLayerTrackletsNum + 1, 0, mMemoryPool.get());
      if (mTaskArena->max_concurrency() <= 1) {
        for (int iTracklet{0}; iTracklet < currentLayerTrackletsNum; ++iTracklet) {
          perTrackletCount[iTracklet] = forTrackletCells(PassMode::OnePass{}, firstTransitionId, secondTransitionId, hitLayers, layerCells, iTracklet);
        }
        std::exclusive_scan(perTrackletCount.begin(), perTrackletCount.end(), perTrackletCount.begin(), 0);
      } else {
        tbb::parallel_for(0, currentLayerTrackletsNum, [&](const int iTracklet) {
          perTrackletCount[iTracklet] = forTrackletCells(PassMode::TwoPassCount{}, firstTransitionId, secondTransitionId, hitLayers, layerCells, iTracklet);
        });

        std::exclusive_scan(perTrackletCount.begin(), perTrackletCount.end(), perTrackletCount.begin(), 0);
        auto totalCells{perTrackletCount.back()};
        if (totalCells == 0) {
          auto& lut = mScratch->getCellsLookupTable()[cellTopologyId];
          lut.resize(currentLayerTrackletsNum + 1);
          std::fill(lut.begin(), lut.end(), 0);
          continue;
        }
        layerCells.resize(totalCells);

        tbb::parallel_for(0, currentLayerTrackletsNum, [&](const int iTracklet) {
          int offset = perTrackletCount[iTracklet];
          if (offset == perTrackletCount[iTracklet + 1]) {
            return;
          }
          forTrackletCells(PassMode::TwoPassInsert{}, firstTransitionId, secondTransitionId, hitLayers, layerCells, iTracklet, offset);
        });
      }

      auto& lut = mScratch->getCellsLookupTable()[cellTopologyId];
      lut.resize(currentLayerTrackletsNum + 1);
      std::copy_n(perTrackletCount.begin(), currentLayerTrackletsNum + 1, lut.begin());

      if (mScratch->hasMCinformation() && mTrkParams[iteration].CreateArtefactLabels) {
        auto& labels = mScratch->getCellsLabel(cellTopologyId);
        labels.reserve(layerCells.size());
        for (const auto& cell : layerCells) {
          MCCompLabel currentLab{mScratch->getTrackletsLabel(firstTransitionId)[cell.getFirstTrackletIndex()]};
          MCCompLabel nextLab{mScratch->getTrackletsLabel(secondTransitionId)[cell.getSecondTrackletIndex()]};
          labels.emplace_back(currentLab == nextLab ? currentLab : MCCompLabel());
        }
      }
    }
  });
}

template <int NLayers>
void TrackerTraits<NLayers>::findCellsNeighbours(const int iteration)
{
  if (!mTraversalGrouping.has_value()) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidTraversalSchedule};
  }
  dispatchTransitionPolicies(*mTraversalGrouping, [&](auto traits, auto, auto) {
    using Traits = decltype(traits);
    if constexpr (stateFamilyFromNLayers<NLayers>() != Traits::Family) {
      throw TraversalException{iteration, TraversalFailureReason::StateFamilyMismatch};
    } else if constexpr (Traits::Tag == TransitionPolicyTag::CylinderCylinder) {
      if (!mCylinderPolicyParams.has_value()) {
        throw TraversalException{iteration, TraversalFailureReason::InvalidPolicyParameters};
      }
      findCellsNeighboursForPolicy<Traits::Tag>(iteration, mTraversalGrouping->scheduledCellsForTag(Traits::Tag), *mCylinderPolicyParams);
    } else if constexpr (Traits::Tag == TransitionPolicyTag::DiskDisk) {
      if (!mDiskPolicyParams.has_value()) {
        throw TraversalException{iteration, TraversalFailureReason::InvalidPolicyParameters};
      }
      findCellsNeighboursForPolicy<Traits::Tag>(iteration, mTraversalGrouping->scheduledCellsForTag(Traits::Tag), *mDiskPolicyParams);
    }
  });
}

template <int NLayers>
template <TransitionPolicyTag Tag>
void TrackerTraits<NLayers>::findCellsNeighboursForPolicy(
  int iteration,
  gsl::span<const CellTopologyId> scheduledCells,
  const typename TransitionPolicyTraits<Tag>::Params& params)
{
  const auto topology = mTraversalLayout.topology;
  if (mScratch->getCells().size() != topology.nCells ||
      mScratch->getCellsLookupTable().size() != topology.nCells ||
      mScratch->getCellsNeighbours().size() != topology.nCells ||
      mScratch->getCellsNeighboursTopology().size() != topology.nCells ||
      mScratch->getCellsNeighboursLUT().size() != topology.nCells) {
    throw TraversalException{iteration, TraversalFailureReason::LegacyIndexMismatch};
  }
  mTaskArena->execute([&] {
    std::vector<bounded_vector<CellNeighbour>> cellsNeighboursByTarget;
    cellsNeighboursByTarget.reserve(topology.nCells);
    for (uint32_t cellTopologyId = 0; cellTopologyId < topology.nCells; ++cellTopologyId) {
      deepVectorClear(mScratch->getCellsNeighbours()[cellTopologyId]);
      deepVectorClear(mScratch->getCellsNeighboursTopology()[cellTopologyId]);
      deepVectorClear(mScratch->getCellsNeighboursLUT()[cellTopologyId]);
      cellsNeighboursByTarget.emplace_back(mMemoryPool.get());
    }

    for (const auto scheduledId : scheduledCells) {
      const auto cellTopologyId = scheduledId.value();
      if (cellTopologyId >= topology.nCells || cellTopologyId >= mScratch->getCells().size() ||
          cellTopologyId >= mScratch->getCellsLookupTable().size()) {
        throw TraversalException{iteration, TraversalFailureReason::LegacyIndexMismatch};
      }
      const auto& cellTopology = topology.getCell(scheduledId);
      if (mScratch->getCells()[cellTopologyId].empty()) {
        continue;
      }
      const auto successors = topology.getCellsStartingWithTransition(cellTopology.secondTransition);
      if (!successors.getEntries()) {
        continue;
      }

      tbb::enumerable_thread_specific<bounded_vector<CellNeighbour>> sourceNeighbours([&]() { return bounded_vector<CellNeighbour>{mMemoryPool.get()}; });
      tbb::parallel_for(0, static_cast<int>(mScratch->getCells()[cellTopologyId].size()), [&](const int iCell) {
        auto& localNeighbours = sourceNeighbours.local();
        const auto& currentCellSeed{mScratch->getCells()[cellTopologyId][iCell]};
        const int nextLayerTrackletIndex{currentCellSeed.getSecondTrackletIndex()};
        for (uint32_t iSuccessor = 0; iSuccessor < successors.getEntries(); ++iSuccessor) {
          const auto nextTopologyId = topology.cellsByFirstTransition[successors.getFirstEntry() + iSuccessor];
          const auto nextCellTopologyId = nextTopologyId.value();
          if (nextCellTopologyId >= topology.nCells || nextCellTopologyId >= mScratch->getCells().size() ||
              nextCellTopologyId >= mScratch->getCellsLookupTable().size()) {
            throw TraversalException{iteration, TraversalFailureReason::LegacyIndexMismatch};
          }
          if (mScratch->getCells()[nextCellTopologyId].empty() ||
              mScratch->getCellsLookupTable()[nextCellTopologyId].empty()) {
            continue;
          }
          const auto& nextCellLUT = mScratch->getCellsLookupTable()[nextCellTopologyId];
          if (nextLayerTrackletIndex < 0 || nextLayerTrackletIndex + 1 >= static_cast<int>(nextCellLUT.size())) {
            continue;
          }
          const int nextLayerFirstCellIndex{nextCellLUT[nextLayerTrackletIndex]};
          const int nextLayerLastCellIndex{nextCellLUT[nextLayerTrackletIndex + 1]};
          if (nextLayerFirstCellIndex < 0 || nextLayerLastCellIndex < nextLayerFirstCellIndex ||
              nextLayerLastCellIndex > static_cast<int>(mScratch->getCells()[nextCellTopologyId].size())) {
            throw TraversalException{iteration, TraversalFailureReason::LegacyIndexMismatch};
          }
          for (int iNextCell{nextLayerFirstCellIndex}; iNextCell < nextLayerLastCellIndex; ++iNextCell) {
            const auto& nextCellSeedRef{mScratch->getCells()[nextCellTopologyId][iNextCell]};
            if (nextCellSeedRef.getFirstTrackletIndex() != nextLayerTrackletIndex || !currentCellSeed.getTimeStamp().isCompatible(nextCellSeedRef.getTimeStamp())) {
              break;
            }

            if (!o2::itsmft::tracking::cellsAreCompatible<Tag>(currentCellSeed.state(), nextCellSeedRef.state(),
                                                               currentCellSeed.getSecondClusterIndex(), nextCellSeedRef.getFirstClusterIndex(),
                                                               getBz(), params)) {
              continue;
            }

            const int nextLevel = currentCellSeed.getLevel() + 1;
            localNeighbours.emplace_back(cellTopologyId, iCell, nextCellTopologyId, iNextCell, nextLevel);
          }
        }
      });

      bounded_vector<size_t> count(topology.nCells, 0, mMemoryPool.get());
      for (const auto& localNeighbours : sourceNeighbours) {
        for (const auto& neigh : localNeighbours) {
          ++count[neigh.nextCellTopology];
        }
      }
      for (size_t i{0}; i < topology.nCells; ++i) {
        cellsNeighboursByTarget[i].reserve(count[i]);
      }
      for (const auto& localNeighbours : sourceNeighbours) {
        for (const auto& neigh : localNeighbours) {
          cellsNeighboursByTarget[neigh.nextCellTopology].emplace_back(neigh);
          if (neigh.level > mScratch->getCells()[neigh.nextCellTopology][neigh.nextCell].getLevel()) {
            mScratch->getCells()[neigh.nextCellTopology][neigh.nextCell].setLevel(neigh.level);
          }
        }
      }
    }

    for (uint32_t cellTopologyId = 0; cellTopologyId < topology.nCells; ++cellTopologyId) {
      auto& cellsNeighbours = cellsNeighboursByTarget[cellTopologyId];
      if (cellsNeighbours.empty()) {
        continue;
      }

      std::sort(cellsNeighbours.begin(), cellsNeighbours.end(), [](const auto& a, const auto& b) {
        return a.nextCell < b.nextCell;
      });

      auto& cellsNeighbourLUT = mScratch->getCellsNeighboursLUT()[cellTopologyId];
      cellsNeighbourLUT.assign(mScratch->getCells()[cellTopologyId].size(), 0);
      for (const auto& neigh : cellsNeighbours) {
        ++cellsNeighbourLUT[neigh.nextCell];
      }
      std::inclusive_scan(cellsNeighbourLUT.begin(), cellsNeighbourLUT.end(), cellsNeighbourLUT.begin());

      mScratch->getCellsNeighbours()[cellTopologyId].reserve(cellsNeighbours.size());
      mScratch->getCellsNeighboursTopology()[cellTopologyId].reserve(cellsNeighbours.size());
      std::ranges::transform(cellsNeighbours, std::back_inserter(mScratch->getCellsNeighbours()[cellTopologyId]), [](const auto& neigh) { return neigh.cell; });
      std::ranges::transform(cellsNeighbours, std::back_inserter(mScratch->getCellsNeighboursTopology()[cellTopologyId]), [](const auto& neigh) { return neigh.cellTopology; });
    }

    // clean up LUTs
    for (auto& cellLUT : mScratch->getCellsLookupTable()) {
      deepVectorClear(cellLUT);
    }
  });
}

template <int NLayers>
template <TransitionPolicyTag Tag, typename InputSeed>
void TrackerTraits<NLayers>::processNeighbours(int iteration, int defaultCellTopologyId, int iLevel, const bounded_vector<InputSeed>& currentCellSeed, const bounded_vector<int>& currentCellId, const bounded_vector<int>& currentCellTopologyId, bounded_vector<TrackSeedN>& updatedCellSeeds, bounded_vector<int>& updatedCellsIds, bounded_vector<int>& updatedCellsTopologyIds, const typename TransitionPolicyTraits<Tag>::Params& params)
{
  const auto layerMaterial = mAttachHitConfig.layerMaterial;

  mTaskArena->execute([&] {
    auto forCellNeighbours = [&](auto Mode, int iCell, int offset = 0) -> int {
      const auto& currentCell{currentCellSeed[iCell]};
      const int cellTopologyId = currentCellTopologyId.empty() ? defaultCellTopologyId : currentCellTopologyId[iCell];

      if constexpr (decltype(Mode)::value != PassMode::TwoPassInsert::value) {
        if (currentCell.getLevel() != iLevel) {
          return 0;
        }
        if (currentCellId.empty()) {
          for (int layer = 0; layer < NLayers; ++layer) {
            const int clusterIndex = currentCell.getCluster(layer);
            if (clusterIndex != constants::UnusedIndex && mScratch->isClusterUsed(layer, clusterIndex)) {
              return 0; /// this we do only on the first iteration, hence the check on currentCellId
            }
          }
        }
      }

      const int cellId = currentCellId.empty() ? iCell : currentCellId[iCell];
      if (cellTopologyId < 0 || mScratch->getCellsNeighboursLUT()[cellTopologyId].empty()) {
        return 0;
      }
      const int startNeighbourId{cellId ? mScratch->getCellsNeighboursLUT()[cellTopologyId][cellId - 1] : 0};
      const int endNeighbourId{mScratch->getCellsNeighboursLUT()[cellTopologyId][cellId]};
      int foundSeeds{0};
      for (int iNeighbourCell{startNeighbourId}; iNeighbourCell < endNeighbourId; ++iNeighbourCell) {
        const int neighbourCellTopologyId = mScratch->getCellsNeighboursTopology()[cellTopologyId][iNeighbourCell];
        const int neighbourCellId = mScratch->getCellsNeighbours()[cellTopologyId][iNeighbourCell];
        const auto& neighbourCell = mScratch->getCells()[neighbourCellTopologyId][neighbourCellId];
        if (neighbourCell.getSecondTrackletIndex() != currentCell.getFirstTrackletIndex()) {
          continue;
        }
        if (!currentCell.getTimeStamp().isCompatible(neighbourCell.getTimeStamp())) {
          continue;
        }
        if (currentCell.getLevel() - 1 != neighbourCell.getLevel()) {
          continue;
        }
        const int neighbourLayer = neighbourCell.getInnerLayer();
        const int neighbourCluster = neighbourCell.getFirstClusterIndex();
        if (mScratch->isClusterUsed(neighbourLayer, neighbourCluster)) {
          continue;
        }

        /// Let's start the fitting procedure
        TrackSeedN seed{currentCell};
        seed.getTimeStamp() = currentCell.getTimeStamp();
        seed.getTimeStamp() += neighbourCell.getTimeStamp();

        const auto& measurement = mLayerMeasurements[neighbourLayer][neighbourCluster];
        float chi2 = seed.getChi2();
        OperationFailureReason attachReason{};
        if (!o2::itsmft::tracking::attachHit<Tag>(seed.state(), measurement, layerMaterial[neighbourLayer], getBz(), chi2, params, attachReason)) {
          continue;
        }
        seed.setChi2(chi2);

        if constexpr (decltype(Mode)::value != PassMode::TwoPassCount::value) {
          seed.getClusters()[neighbourLayer] = neighbourCluster;
          auto mask = seed.getHitLayerMask();
          mask.set(neighbourLayer);
          seed.setHitLayerMask(mask);
          seed.setLevel(neighbourCell.getLevel());
          seed.setFirstTrackletIndex(neighbourCell.getFirstTrackletIndex());
          seed.setSecondTrackletIndex(neighbourCell.getSecondTrackletIndex());
        }

        if constexpr (decltype(Mode)::value == PassMode::OnePass::value) {
          updatedCellSeeds.push_back(seed);
          updatedCellsIds.push_back(neighbourCellId);
          updatedCellsTopologyIds.push_back(neighbourCellTopologyId);
        } else if constexpr (decltype(Mode)::value == PassMode::TwoPassCount::value) {
          ++foundSeeds;
        } else if constexpr (decltype(Mode)::value == PassMode::TwoPassInsert::value) {
          updatedCellSeeds[offset] = seed;
          updatedCellsIds[offset] = neighbourCellId;
          updatedCellsTopologyIds[offset++] = neighbourCellTopologyId;
        } else {
          static_assert(false, "Unknown mode!");
        }
      }
      return foundSeeds;
    };

    const int nCells = static_cast<int>(currentCellSeed.size());
    if (mTaskArena->max_concurrency() <= 1) {
      for (int iCell{0}; iCell < nCells; ++iCell) {
        forCellNeighbours(PassMode::OnePass{}, iCell);
      }
    } else {
      bounded_vector<int> perCellCount(nCells + 1, 0, mMemoryPool.get());
      tbb::parallel_for(0, nCells, [&](const int iCell) {
        perCellCount[iCell] = forCellNeighbours(PassMode::TwoPassCount{}, iCell);
      });

      std::exclusive_scan(perCellCount.begin(), perCellCount.end(), perCellCount.begin(), 0);
      auto totalNeighbours{perCellCount.back()};
      if (totalNeighbours == 0) {
        return;
      }
      updatedCellSeeds.resize(totalNeighbours);
      updatedCellsIds.resize(totalNeighbours);
      updatedCellsTopologyIds.resize(totalNeighbours);

      tbb::parallel_for(0, nCells, [&](const int iCell) {
        int offset = perCellCount[iCell];
        if (offset == perCellCount[iCell + 1]) {
          return;
        }
        forCellNeighbours(PassMode::TwoPassInsert{}, iCell, offset);
      });
    }
  });
}

template <int NLayers>
void TrackerTraits<NLayers>::findRoads(const int iteration)
{
  if (!mTraversalGrouping.has_value()) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidTraversalSchedule};
  }
  // Defensive parity check (Architecture.md Sec 9 / validateLegacyParity()):
  // findRoadsForPolicy() below indexes mScratch->getCells() -- a legacy,
  // per-legacy-CellTopologyId TimeFrame array -- with the .value() of sparse
  // CellTopologyIds returned by roadStartCellsForTag(). That is only safe
  // because validateLegacyParity(), run once per initialiseTimeFrame() before
  // mTraversalLayout/mTraversalGrouping are committed, already proves the
  // sparse and legacy topologies have identical cell count and per-index
  // correspondence for this iteration. Re-checking the count here, at the
  // road-finding boundary and before any indexing, turns a future desync
  // between that cached proof and TimeFrame's own containers into an
  // explicit failure instead of an out-of-bounds/misaligned read.
  if (mScratch->getCells().size() != mTraversalLayout.topology.nCells) {
    throw TraversalException{iteration, TraversalFailureReason::LegacyIndexMismatch};
  }
  dispatchTransitionPolicies(*mTraversalGrouping, [&](auto traits, auto, auto) {
    using Traits = decltype(traits);
    if constexpr (stateFamilyFromNLayers<NLayers>() != Traits::Family) {
      throw TraversalException{iteration, TraversalFailureReason::StateFamilyMismatch};
    } else if constexpr (Traits::Tag == TransitionPolicyTag::CylinderCylinder) {
      if (!mCylinderPolicyParams.has_value()) {
        throw TraversalException{iteration, TraversalFailureReason::InvalidPolicyParameters};
      }
      findRoadsForPolicy<Traits::Tag>(iteration, *mCylinderPolicyParams);
    } else if constexpr (Traits::Tag == TransitionPolicyTag::DiskDisk) {
      if (!mDiskPolicyParams.has_value()) {
        throw TraversalException{iteration, TraversalFailureReason::InvalidPolicyParameters};
      }
      findRoadsForPolicy<Traits::Tag>(iteration, *mDiskPolicyParams);
    }
  });
}

template <int NLayers>
template <TransitionPolicyTag Tag>
void TrackerTraits<NLayers>::findRoadsForPolicy(const int iteration, const typename TransitionPolicyTraits<Tag>::Params& params)
{
  bounded_vector<bounded_vector<int>> firstClusters(mTrkParams[iteration].NLayers, bounded_vector<int>(mMemoryPool.get()), mMemoryPool.get());
  firstClusters.resize(mTrkParams[iteration].NLayers);
  const auto propagator = o2::base::Propagator::Instance();
  const TrackingFrameInfo* tfInfos[NLayers]{};
  const Cluster* unsortedClusters[NLayers]{};
  for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
    tfInfos[iLayer] = mScratch->getTrackingFrameInfoOnLayer(iLayer).data();
    unsortedClusters[iLayer] = mScratch->getUnsortedClusters()[iLayer].data();
  }
  // Road-start selection is topology-derived, not a StartLayerMask/LayerMask
  // runtime decision (Architecture.md Sec 10, D003): mTraversalGrouping's
  // roadStartCellsForTag(Tag) (TransitionPolicyDispatch.h) is the deterministic
  // ascending-CellTopologyId subsequence of the sparse layout's cells whose
  // traversal endpoint is a seeding SurfaceId, cached once per
  // initialiseTimeFrame() call and reused unchanged across every startLevel
  // pass below and across every repeated findRoads() call in the
  // PerPrimaryVertexProcessing loop (CATracker.cxx). StartLayerMask itself
  // remains a legacy configuration/layout-construction input (see
  // positionalSurfaceMask() in TimeFrame.cxx and validateLegacyParity()) --
  // it is simply no longer read here. Each returned CellTopologyId is a
  // sparse identifier; only its numeric .value() is used, and only to index
  // mScratch->getCells(), the parity-validated legacy per-cell-topology
  // TimeFrame array (see the defensive count check in findRoads()). No
  // SurfaceId is used as a legacy layer/vector index anywhere in this
  // function.
  // Road-length filter bound: maximum absolute q/pT, in the same (GeV/c)^-1
  // units as SurfaceKinematicState::parameters[4]. Applied identically to
  // both families via getQOverPt() (raw signed value, never squared); no
  // NLayers/DetID/state-family branch. std::abs() of a NaN/+-Inf q/pT is
  // never <= this finite bound (standard IEEE-754 comparison semantics), so
  // non-finite seeds are rejected deterministically without extra checks --
  // see testCellRepresentation.cxx's dedicated road-filter tests for focused
  // coverage of that behavior, per the correction's requirement not to rely
  // on it without a test.
  constexpr float maxAbsQOverPt = 1.e3f;
  for (int startLevel{mTrkParams[iteration].CellsPerRoad()}; startLevel >= mTrkParams[iteration].CellMinimumLevel(); --startLevel) {

    auto seedFilter = [&](const auto& seed) {
      return seed.getHitLayerMask().isAllowed(mTrkParams[iteration].MaxHoles, mTrkParams[iteration].HoleLayerMask) &&
             seed.getHitLayerMask().length() >= mTrkParams[iteration].MinTrackLength &&
             std::abs(seed.getQOverPt()) <= maxAbsQOverPt && seed.getChi2() <= mTrkParams[iteration].MaxChi2NDF * ((startLevel + 2) * 2 - 5);
    };

    bounded_vector<TrackSeedN> trackSeeds(mMemoryPool.get());
    for (const auto startId : mTraversalGrouping->roadStartCellsForTag(Tag)) {
      const int startCellTopologyId = startId.value();
      // Cell population is per-event/per-vertex data, never cached in
      // TransitionPolicyGrouping: this check must stay here, evaluated at
      // runtime against the current iVertex's TimeFrame content.
      if (mScratch->getCells()[startCellTopologyId].empty()) {
        continue;
      }

      bounded_vector<int> lastCellId(mMemoryPool.get()), updatedCellId(mMemoryPool.get());
      bounded_vector<int> lastCellTopologyId(mMemoryPool.get()), updatedCellTopologyId(mMemoryPool.get());
      bounded_vector<TrackSeedN> lastCellSeed(mMemoryPool.get()), updatedCellSeed(mMemoryPool.get());

      processNeighbours<Tag>(iteration, startCellTopologyId, startLevel, mScratch->getCells()[startCellTopologyId], lastCellId, lastCellTopologyId, updatedCellSeed, updatedCellId, updatedCellTopologyId, params);

      int level = startLevel;
      while (level > 2 && !updatedCellSeed.empty()) {
        lastCellSeed.swap(updatedCellSeed);
        lastCellId.swap(updatedCellId);
        lastCellTopologyId.swap(updatedCellTopologyId);
        deepVectorClear(updatedCellSeed); /// tame the memory peaks
        deepVectorClear(updatedCellId);   /// tame the memory peaks
        deepVectorClear(updatedCellTopologyId);
        processNeighbours<Tag>(iteration, constants::UnusedIndex, --level, lastCellSeed, lastCellId, lastCellTopologyId, updatedCellSeed, updatedCellId, updatedCellTopologyId, params);
      }
      deepVectorClear(lastCellId);         /// tame the memory peaks
      deepVectorClear(lastCellTopologyId); /// tame the memory peaks
      deepVectorClear(lastCellSeed);       /// tame the memory peaks

      if (!updatedCellSeed.empty()) {
        trackSeeds.reserve(trackSeeds.size() + std::count_if(updatedCellSeed.begin(), updatedCellSeed.end(), seedFilter));
        std::copy_if(updatedCellSeed.begin(), updatedCellSeed.end(), std::back_inserter(trackSeeds), seedFilter);
      }
    }

    if (trackSeeds.empty()) {
      continue;
    }

    using TrackT = typename DetectorTraits<NLayers>::TrackType;
    bounded_vector<TrackT> tracks(mMemoryPool.get());
    mTaskArena->execute([&] {
      auto forSeed = [&](auto Mode, int iSeed, int offset = 0) {
        TrackT temporaryTrack;
        const bool refitSuccess = DetectorTraits<NLayers>::refitSeed(trackSeeds[iSeed],
                                                                     temporaryTrack,
                                                                     mTrkParams[iteration],
                                                                     mBz,
                                                                     *mScratch,
                                                                     tfInfos,
                                                                     unsortedClusters,
                                                                     propagator,
                                                                     mLayerMeasurements);
        if (refitSuccess) {
          DetectorTraits<NLayers>::copySeedPatternToTrack(temporaryTrack, trackSeeds[iSeed]);
          if constexpr (decltype(Mode)::value == PassMode::OnePass::value) {
            tracks.push_back(temporaryTrack);
          } else if constexpr (decltype(Mode)::value == PassMode::TwoPassCount::value) {
            // nothing to do
          } else if constexpr (decltype(Mode)::value == PassMode::TwoPassInsert::value) {
            tracks[offset] = temporaryTrack;
          } else {
            static_assert(false, "Unknown mode!");
          }
          return 1;
        }
        return 0;
      };

      const int nSeeds = static_cast<int>(trackSeeds.size());
      if (mTaskArena->max_concurrency() <= 1) {
        for (int iSeed{0}; iSeed < nSeeds; ++iSeed) {
          forSeed(PassMode::OnePass{}, iSeed);
        }
      } else {
        bounded_vector<int> perSeedCount(nSeeds + 1, 0, mMemoryPool.get());
        tbb::parallel_for(0, nSeeds, [&](const int iSeed) {
          perSeedCount[iSeed] = forSeed(PassMode::TwoPassCount{}, iSeed);
        });

        std::exclusive_scan(perSeedCount.begin(), perSeedCount.end(), perSeedCount.begin(), 0);
        auto totalTracks{perSeedCount.back()};
        if (totalTracks == 0) {
          return;
        }
        tracks.resize(totalTracks);

        tbb::parallel_for(0, nSeeds, [&](const int iSeed) {
          if (perSeedCount[iSeed] == perSeedCount[iSeed + 1]) {
            return;
          }
          forSeed(PassMode::TwoPassInsert{}, iSeed, perSeedCount[iSeed]);
        });
      }

      deepVectorClear(trackSeeds);
    });

    // Same ordering as o2::its::track::isBetter (longer track, then lower chi2).
    std::sort(tracks.begin(), tracks.end(), [](const TrackT& a, const TrackT& b) {
      const auto ncla = a.getNumberOfClusters();
      const auto nclb = b.getNumberOfClusters();
      return (ncla == nclb) ? (a.getChi2() < b.getChi2()) : ncla > nclb;
    });
    acceptTracks(iteration, tracks, firstClusters);
  }
  markTracks(iteration);
}

template <int NLayers>
void TrackerTraits<NLayers>::acceptTracks(int iteration, bounded_vector<CATrackType<NLayers>>& tracks, bounded_vector<bounded_vector<int>>& firstClusters)
{
  auto& trks = mScratch->getTracks();
  trks.reserve(trks.size() + tracks.size());
  const float smallestROFHalf = mScratch->getROFOverlapTableView().getClockLayer().mROFLength * 0.5f;
  for (auto& track : tracks) {
    int nShared = 0;
    bool isFirstShared{false};
    int firstLayer{-1}, firstCluster{-1};
    for (int iLayer{0}; iLayer < mTrkParams[iteration].NLayers; ++iLayer) {
      if (track.getClusterIndex(iLayer) == constants::UnusedIndex) {
        continue;
      }
      bool isShared = mScratch->isClusterUsed(iLayer, track.getClusterIndex(iLayer));
      nShared += int(isShared);
      if (firstLayer < 0) {
        firstCluster = track.getClusterIndex(iLayer);
        isFirstShared = isShared && mTrkParams[iteration].AllowSharingFirstCluster && std::find(firstClusters[iLayer].begin(), firstClusters[iLayer].end(), firstCluster) != firstClusters[iLayer].end();
        firstLayer = iLayer;
      }
    }

    /// do not account for the first cluster in the shared clusters number if it is allowed
    if (nShared - int(isFirstShared && mTrkParams[iteration].AllowSharingFirstCluster) > mTrkParams[iteration].SharedMaxClusters) {
      continue;
    }

    bool firstCls{true}, nominalCompatible{true};
    TimeEstBC nominalTS, expandedTS;
    for (int iLayer{0}; iLayer < mTrkParams[iteration].NLayers; ++iLayer) {
      if (track.getClusterIndex(iLayer) == constants::UnusedIndex) {
        continue;
      }
      mScratch->markUsedCluster(iLayer, track.getClusterIndex(iLayer));
      int currentROF = mScratch->getClusterROF(iLayer, track.getClusterIndex(iLayer));
      const auto nominalROFTS = mScratch->getROFOverlapTableView().getLayer(iLayer).getROFTimeBounds(currentROF);
      const auto expandedROFTS = mScratch->getROFOverlapTableView().getLayer(iLayer).getROFTimeBounds(currentROF, true);
      if (firstCls) {
        firstCls = false;
        nominalTS = nominalROFTS;
        expandedTS = expandedROFTS;
      } else {
        if (nominalCompatible) {
          if (nominalTS.isCompatible(nominalROFTS)) {
            nominalTS += nominalROFTS;
          } else {
            nominalCompatible = false;
          }
        }
        if (!expandedTS.isCompatible(expandedROFTS)) {
          LOGP(fatal, "TS {}+/-{} are incompatible with {}+/-{}, this should not happen!", expandedROFTS.getTimeStamp(), expandedROFTS.getTimeStampError(), expandedTS.getTimeStamp(), expandedTS.getTimeStampError());
        }
        expandedTS += expandedROFTS;
      }
    }
    const auto selectedTimestamp = nominalCompatible ? nominalTS : expandedTS;
    track.getTimeStamp() = selectedTimestamp.makeSymmetrical();
    // this is a sanity clamp
    // we cannot be worse than the clock so we clamp to this
    if (track.getTimeStamp().getTimeStampError() > smallestROFHalf) {
      track.getTimeStamp().setTimeStampError(smallestROFHalf);
    }
    DetectorTraits<NLayers>::clearTransientLayerPattern(track);
    trks.emplace_back(track);

    // Shadow-only owner-thread publication. The candidate/refit loops above
    // may be parallel, but acceptTracks() runs after they have joined and is
    // the one deterministic legacy accepted-track order. Do not move this
    // into a task-arena worker or GPU path.
    CommonTrackShadowRecord shadow;
    if (!makeAcceptedTrackShadow<NLayers>(track, mLayerMeasurements, selectedTimestamp, shadow) ||
        !mAcceptedTrackShadowPublisher.publish(*mFrame, shadow, track)) {
      LOGP(fatal, "CommonTrack shadow construction failed for an accepted {} CA track", DetectorTraits<NLayers>::DetId == o2::detectors::DetID::ITS ? "ITS" : "MFT");
    }

    if (mTrkParams[iteration].AllowSharingFirstCluster) {
      firstClusters[firstLayer].push_back(firstCluster);
    }
  }
}

template <int NLayers>
void TrackerTraits<NLayers>::markTracks(int iteration)
{
  if (mTrkParams[iteration].AllowSharingFirstCluster) {
    /// Now we have to set the shared cluster flag
    auto& tracks = mScratch->getTracks();

    bounded_vector<int> fclusSort(tracks.size(), mMemoryPool.get());
    std::iota(fclusSort.begin(), fclusSort.end(), 0);
    std::sort(fclusSort.begin(), fclusSort.end(), [&tracks](int a, int b) {
      return tracks[a].getClusterIndex(tracks[a].getFirstClusterLayer()) < tracks[b].getClusterIndex(tracks[b].getFirstClusterLayer());
    });

    auto areTracksSelected = [this, iteration](const auto& t1, const auto& t2) {
      const auto t1FirstLayer{t1.getFirstClusterLayer()}, t2FirstLayer{t2.getFirstClusterLayer()};
      if (t1FirstLayer != t2FirstLayer) {
        return false;
      }
      if (mScratch->getClusterROF(t1FirstLayer, t1.getClusterIndex(t1FirstLayer)) != mScratch->getClusterROF(t2FirstLayer, t2.getClusterIndex(t2FirstLayer))) {
        return false;
      }
      if (!math_utils::isPhiDifferenceBelow(t1.getPhi(), t2.getPhi(), mTrkParams[iteration].SharedClusterMaxDeltaPhi)) {
        return false;
      }
      if (std::abs(t1.getEta() - t2.getEta()) > mTrkParams[iteration].SharedClusterMaxDeltaEta) {
        return false;
      }
      if (mTrkParams[iteration].SharedClusterOppositeSign) {
        if (DetectorTraits<NLayers>::haveSamePolarity(t1, t2)) {
          return false;
        }
      }
      return true;
    };

    for (int i{0}; i < static_cast<int>(fclusSort.size()); ++i) {
      auto& track = tracks[fclusSort[i]];
      for (int j{i + 1}; j < static_cast<int>(fclusSort.size()) && tracks[fclusSort[j]].getClusterIndex(tracks[fclusSort[j]].getFirstClusterLayer()) == track.getClusterIndex(track.getFirstClusterLayer()); ++j) {
        auto& track2 = tracks[fclusSort[j]];
        if (areTracksSelected(track, track2)) {
          track.setSharedClusters();
          track2.setSharedClusters();
        }
      }
    }
  }
}

template <int NLayers>
void TrackerTraits<NLayers>::setBz(float bz)
{
  mBz = bz;
  mFrame->setBz(bz);
}

template <int NLayers>
void TrackerTraits<NLayers>::setNThreads(int n, std::shared_ptr<tbb::task_arena>& arena)
{
#if defined(OPTIMISATION_OUTPUT)
  mTaskArena = std::make_shared<tbb::task_arena>(1);
#else
  if (arena == nullptr) {
    mTaskArena = std::make_shared<tbb::task_arena>(std::abs(n));
    LOGP(info, "Setting tracker with {} threads.", n);
  } else {
    mTaskArena = arena;
  }
#endif
}

template class TrackerTraits<7>;
template class TrackerTraits<10>;
template void TrackerTraits<7>::processNeighbours<TransitionPolicyTag::CylinderCylinder, typename TrackerTraits<7>::CellSeedN>(int, int, int, const bounded_vector<typename TrackerTraits<7>::CellSeedN>&, const bounded_vector<int>&, const bounded_vector<int>&, bounded_vector<typename TrackerTraits<7>::TrackSeedN>&, bounded_vector<int>&, bounded_vector<int>&, const CylinderCylinderPolicyParams&);
template void TrackerTraits<7>::processNeighbours<TransitionPolicyTag::CylinderCylinder, typename TrackerTraits<7>::TrackSeedN>(int, int, int, const bounded_vector<typename TrackerTraits<7>::TrackSeedN>&, const bounded_vector<int>&, const bounded_vector<int>&, bounded_vector<typename TrackerTraits<7>::TrackSeedN>&, bounded_vector<int>&, bounded_vector<int>&, const CylinderCylinderPolicyParams&);
template void TrackerTraits<10>::processNeighbours<TransitionPolicyTag::DiskDisk, typename TrackerTraits<10>::CellSeedN>(int, int, int, const bounded_vector<typename TrackerTraits<10>::CellSeedN>&, const bounded_vector<int>&, const bounded_vector<int>&, bounded_vector<typename TrackerTraits<10>::TrackSeedN>&, bounded_vector<int>&, bounded_vector<int>&, const DiskDiskPolicyParams&);
template void TrackerTraits<10>::processNeighbours<TransitionPolicyTag::DiskDisk, typename TrackerTraits<10>::TrackSeedN>(int, int, int, const bounded_vector<typename TrackerTraits<10>::TrackSeedN>&, const bounded_vector<int>&, const bounded_vector<int>&, bounded_vector<typename TrackerTraits<10>::TrackSeedN>&, bounded_vector<int>&, bounded_vector<int>&, const DiskDiskPolicyParams&);

} // namespace o2::itsmft::tracking

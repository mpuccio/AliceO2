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
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/enumerable_thread_specific.h>

#include "Framework/Logger.h"
#include "GPUCommonMath.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/detail/CommonTrackShadow.h"
#include "ITStracking/Constants.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IndexTableConfiguration.h"
#include "ITSMFTTracking/detail/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/SurfaceMask.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "ITSMFTTracking/detail/CandidateFinding.h"
#include "ITStracking/Tracklet.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::itsmft::tracking
{

namespace math_utils = o2::its::math_utils;
using o2::its::deepVectorClear;
using o2::its::TimeEstBC;

namespace detail
{
// SurfacePlanBinding's compact owned-surface index is the one translation
// used by the tracker while resolving a plan into layer-local storage.
inline std::optional<uint16_t> ownedSurfacePosition(const SurfacePlanBinding& binding, SurfaceId id) noexcept
{
  return binding.getOwnedSurfaceIndex(id);
}
} // namespace detail

struct PassMode {
  using OnePass = std::integral_constant<int, 0>;
  using TwoPassCount = std::integral_constant<int, 1>;
  using TwoPassInsert = std::integral_constant<int, 2>;
};

namespace
{
constexpr uint8_t kCompatibilityAbsCharge = 1;
const o2::track::PID kCompatibilityPID = o2::track::PID::Pion;

constexpr std::size_t kindIndex(SurfaceKind kind) noexcept
{
  return kind == SurfaceKind::Cylinder ? 0u : 1u;
}

bool makeCandidateShadow(const TrackingCandidate& candidate,
                         gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                         CommonTrackShadowRecord& record) noexcept
{
  record = {};
  record.track = candidate.track;
  record.track.hitSurfaces = {};
  record.references.reserve(layerMeasurements.size());
  for (std::size_t position = 0; position < layerMeasurements.size(); ++position) {
    const int localIndex = candidate.getClusterIndex(static_cast<int>(position));
    if (localIndex == o2::its::constants::UnusedIndex) {
      continue;
    }
    if (localIndex < 0 || static_cast<std::size_t>(localIndex) >= layerMeasurements[position].size()) {
      return false;
    }
    const auto& measurement = layerMeasurements[position][localIndex];
    const TrackClusterReference reference{measurement.surface, SurfaceMeasurementIndex{static_cast<uint32_t>(localIndex)}};
    if (!reference.surface.isValid() || measurement.surface != reference.surface || !measurement.cluster.isValid()) {
      return false;
    }
    record.references.push_back(reference);
    record.track.hitSurfaces.set(reference.surface);
  }
  return !record.references.empty();
}

// A static "diamond" vertex (ITSCommonCATrackerParam.useDiamond /
// TrackerTraits::computeLayerTracklets) carries no genuine
// per-event timing: it stands in for every real primary vertex at once, so
// it must compare as time-compatible with whichever ROF it is being tested
// against. Vertex::getTimeStamp() is a fixed-width TimeEstBC whose error
// field is only 16 bits (DataFormatsITS/TimeEstBC.h) -- sized for one ROF's
// own uncertainty, not for an entire TimeFrame's real BC span (order
// 1e5-1e6 BC for a realistic TF) -- so one static timestamp value literally
// cannot overlap every ROF in a real TimeFrame at once; the correct
// full-TimeFrame envelope has to be re-derived at each (layer, rofId) this
// vertex is actually tested against, from that ROF's own real interval
// (RuntimeROFOverlapView::getLayer(layer).getROFTimeBounds(rofId, true) -- the
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

// Host-only conversion from ROOT-visible TrackingParameters into the one
// device-portable record. This is intentionally private to the initialization
// seam and is called once per iteration.
TrackingKernelParameters bindTrackingKernelParameters(const TrackingParameters& params) noexcept
{
  TrackingKernelParameters out;
  out.trackletMinPt = params.TrackletMinPt;
  out.nSigmaCut = params.NSigmaCut;
  out.maxChi2ClusterAttachment = params.MaxChi2ClusterAttachment;
  out.maxChi2NDF = params.MaxChi2NDF;
  out.pvResolution = params.PVres;
  return out;
}

} // namespace

void TrackerTraits::resetTraversalCache() noexcept
{
  mTraversalGraph = {};
  mTraversalCacheValid = false;
  mKernelParameters = {};
  mDiskLayerReferenceZ = {};
  mDiskLayerReferenceZStorage.clear();
  mAttachHitConfig = {};
  const auto resetSurfaceCount = mScratch == nullptr ? std::size_t{0} : mScratch->getNOwnedSurfaces();
  mLayerMaterial.assign(resetSurfaceCount, NominalSurfaceMaterial{});
  mActiveKind.reset();
  mLayerMeasurements.assign(resetSurfaceCount, gsl::span<const SurfaceMeasurement>{});
  for (auto& group : mTransitionsByKind) {
    group.clear();
  }
  for (auto& group : mScheduledCellsByKind) {
    group.clear();
  }
  for (auto& group : mRoadStartCellsByKind) {
    group.clear();
  }
  mTraversalGroupingCount = 0;
}

int TrackerTraits::requireScratchTransitionSlot(int iteration, TransitionId id) const
{
  if (mBinding == nullptr) {
    throw TraversalException{iteration, TraversalFailureReason::MissingLayout};
  }
  const auto slot = mBinding->getScratchTransitionSlot(id);
  if (!slot) {
    throw TraversalException{iteration, TraversalFailureReason::TraversalBindingMismatch};
  }
  return static_cast<int>(*slot);
}

int TrackerTraits::requireScratchCellSlot(int iteration, CellTopologyId id) const
{
  if (mBinding == nullptr) {
    throw TraversalException{iteration, TraversalFailureReason::MissingLayout};
  }
  const auto slot = mBinding->getScratchCellSlot(id);
  if (!slot) {
    throw TraversalException{iteration, TraversalFailureReason::TraversalBindingMismatch};
  }
  return static_cast<int>(*slot);
}

int TrackerTraits::requireSurfacePosition(int iteration, SurfaceId id) const
{
  if (!id.isValid()) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  if (mBinding == nullptr) {
    throw TraversalException{iteration, TraversalFailureReason::MissingLayout};
  }
  const auto position = mBinding->getOwnedSurfaceIndex(id);
  if (!position || *position >= mScratch->getNOwnedSurfaces()) {
    throw TraversalException{iteration, TraversalFailureReason::TraversalBindingMismatch};
  }
  return static_cast<int>(*position);
}

void TrackerTraits::findCellsNeighboursCylinder(int iteration)
{
  findCellsNeighboursForKind<SurfaceKind::Cylinder>(iteration, mScheduledCellsByKind[0], mKernelParameters);
}

void TrackerTraits::findCellsNeighboursDisk(int iteration)
{
  findCellsNeighboursForKind<SurfaceKind::Disk>(iteration, mScheduledCellsByKind[1], mKernelParameters);
}

void TrackerTraits::findRoadsCylinder(int iteration, SeedRefitFunction refitFunction)
{
  ClusterSourceId source{};
  for (const auto transitionId : mTransitionsByKind[0]) {
    const auto position = requireSurfacePosition(iteration, mTraversalGraph.getTransition(transitionId).from);
    if (!mLayerMeasurements[position].empty()) {
      source = mLayerMeasurements[position].front().cluster.source;
      break;
    }
  }
  findRoadsForKind<SurfaceKind::Cylinder>(iteration, mKernelParameters, refitFunction, mRoadStartCellsByKind[0], source);
}

void TrackerTraits::findRoadsDisk(int iteration, SeedRefitFunction refitFunction)
{
  ClusterSourceId source{};
  for (const auto transitionId : mTransitionsByKind[1]) {
    const auto position = requireSurfacePosition(iteration, mTraversalGraph.getTransition(transitionId).from);
    if (!mLayerMeasurements[position].empty()) {
      source = mLayerMeasurements[position].front().cluster.source;
      break;
    }
  }
  findRoadsForKind<SurfaceKind::Disk>(iteration, mKernelParameters, refitFunction, mRoadStartCellsByKind[1], source);
}

void TrackerTraits::validateSparsePlan(int iteration,
                                       const SurfaceGraphView& layout,
                                       std::optional<SurfaceKind>& activeKind,
                                       bool& mixedKind) const
{
  const auto fail = [iteration]() { throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch}; };
  const auto& topology = layout;
  if (layout.surfaces == nullptr || layout.nSurfaces == 0 ||
      (topology.nTransitions != 0 && (topology.transitions == nullptr || topology.cellsByFirstTransitionOffsets == nullptr)) ||
      (topology.nCells != 0 && (topology.cells == nullptr || topology.cellsByFirstTransition == nullptr))) {
    fail();
  }

  const auto kindOf = [&layout](SurfaceId surface) {
    if (!surface.isValid() || surface.value() >= layout.nSurfaces) {
      return std::optional<SurfaceKind>{};
    }
    return std::optional<SurfaceKind>{layout.getSurface(surface).kind};
  };
  const auto observeKind = [&](std::optional<SurfaceKind> kind) {
    if (!kind) {
      fail();
    }
    if (!activeKind) {
      activeKind = kind;
    } else if (*activeKind != *kind) {
      mixedKind = true;
    }
  };

  activeKind.reset();
  mixedKind = false;
  if (mBinding == nullptr) {
    fail();
  }
  const auto transitions = mBinding->getGlobalTransitions();
  const auto cells = mBinding->getGlobalCells();
  if (transitions.empty() || transitions.size() > topology.nTransitions || cells.size() > topology.nCells) {
    fail();
  }
  for (const auto id : transitions) {
    if (!id.isValid() || id.value() >= topology.nTransitions || !mBinding->getScratchTransitionSlot(id)) {
      fail();
    }
    const auto& transition = topology.getTransition(id);
    if (!mBinding->getOwnedSurfaceIndex(transition.from) || !mBinding->getOwnedSurfaceIndex(transition.to) ||
        !transition.skippedSurfaces.isSubsetOf(mBinding->getOwnedSurfaces())) {
      fail();
    }
    observeKind(kindOf(transition.from));
    if (kindOf(transition.to) != activeKind) {
      mixedKind = true;
    }
  }
  for (const auto id : cells) {
    if (!id.isValid() || id.value() >= topology.nCells || !mBinding->getScratchCellSlot(id)) {
      fail();
    }
    const auto& cell = topology.getCell(id);
    if (!mBinding->getScratchTransitionSlot(cell.firstTransition) ||
        !mBinding->getScratchTransitionSlot(cell.secondTransition) ||
        !cell.hitSurfaces.isSubsetOf(mBinding->getOwnedSurfaces())) {
      fail();
    }
    const auto firstKind = kindOf(topology.getTransition(cell.firstTransition).from);
    const auto secondKind = kindOf(topology.getTransition(cell.secondTransition).from);
    observeKind(firstKind);
    if (secondKind != firstKind) {
      mixedKind = true;
    }
  }
  for (const auto id : mBinding->getGlobalScheduledCells()) {
    if (!mBinding->getScratchCellSlot(id)) {
      fail();
    }
  }
  for (const auto id : mBinding->getGlobalRoadStartCells()) {
    if (!mBinding->getScratchCellSlot(id)) {
      fail();
    }
  }
  (void)mTrkParams[iteration];
}

void TrackerTraits::initialiseTimeFrame(const int iteration, const std::vector<SurfaceGraph>& graphs)
{
  resetTraversalCache();

  // 1. Check iteration bounds before accessing configuration.
  if (iteration < 0 || static_cast<size_t>(iteration) >= graphs.size()) {
    throw TraversalException{iteration, TraversalFailureReason::IterationOutOfRange};
  }
  const auto layout = graphs[iteration].getView();
  if (mTrkParams[iteration].PassFlags[IterationStep::FirstPass]) {
    clearAcceptedTracksForSharedStatus();
  }

  // 2. Resolve the one participant binding and its active family. Schedule
  // validation and all filtered ordering are performed by SurfacePlanBinding.
  if (mBinding == nullptr) {
    throw TraversalException{iteration, TraversalFailureReason::MissingLayout};
  }
  ++mTraversalGroupingCount;

  // Resolve the active SurfaceKind from the sparse-plan endpoint before
  // any kind-specific binding runs. A combined layout may contain both
  // families, so a participant binding supplies the one kind this tracker
  // owns; no layer-count selection is used here.
  const auto boundTransitions = mBinding->getGlobalTransitions();
  if (boundTransitions.empty() || boundTransitions.front().value() >= layout.nTransitions) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  mActiveKind = layout.getSurface(layout.getTransition(boundTransitions.front()).from).kind;

  // 2.1 Validate material-correction support before staging or mutating
  // TimeFrame tracking state. Unsupported modes are structural; unrecognized
  // values remain classified by AttachHitConfigView::isValid().
  MaterialCorrectionModeSupport materialModeSupport = MaterialCorrectionModeSupport::Supported;
  if (mActiveKind) {
    materialModeSupport = materialCorrectionModeSupport(*mActiveKind, mTrkParams[iteration].CorrType);
  }
  if (materialModeSupport == MaterialCorrectionModeSupport::Unsupported) {
    throw TraversalException{iteration, TraversalFailureReason::UnsupportedMaterialCorrectionMode};
  }

  // 2.5. Resolve and validate per-surface-position material from the layout
  // and iteration parameters before mutating TimeFrame tracking state.
  // The surface-position mapping is this graph vector's own validated
  // orderedSurfaces (never inferred from a detector identity,
  // radius, z, or numeric ordering); a size mismatch, an invalid/out-of-range
  // mapped SurfaceId, or a numeric disagreement against the adapter-owned
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
  const auto orderedSurfaces = mBinding->getOrderedSurfaces();
  const int activeSurfaceCount = static_cast<int>(mScratch->getNOwnedSurfaces());
  // This is the remaining adapter-edge compatibility check: the application
  // configuration must agree with the adopted plan count. No shared
  // traversal loop uses the compatibility field as its extent.
  if (activeSurfaceCount <= 0 || activeSurfaceCount > MaxLayoutSurfaces ||
      orderedSurfaces.size() != static_cast<std::size_t>(activeSurfaceCount) ||
      mTrkParams[iteration].NLayers != activeSurfaceCount) {
    throw TraversalException{iteration, TraversalFailureReason::LegacyMaterialMismatch};
  }
  std::vector<NominalSurfaceMaterial> stagedLayerMaterial(static_cast<std::size_t>(activeSurfaceCount));
  for (int surfacePosition = 0; surfacePosition < activeSurfaceCount; ++surfacePosition) {
    const auto surfaceId = orderedSurfaces[surfacePosition];
    if (!surfaceId.isValid() || surfaceId.value() >= layout.nSurfaces) {
      throw TraversalException{iteration, TraversalFailureReason::LegacyMaterialMismatch};
    }
    if (std::find(orderedSurfaces.begin(), orderedSurfaces.begin() + surfacePosition, surfaceId) != orderedSurfaces.begin() + surfacePosition) {
      throw TraversalException{iteration, TraversalFailureReason::SurfaceLayerMappingMismatch};
    }
    stagedLayerMaterial[surfacePosition] = layout.getSurface(surfaceId).material;
  }
  if (mTrkParams[iteration].LayerxX0.size() != static_cast<size_t>(activeSurfaceCount)) {
    throw TraversalException{iteration, TraversalFailureReason::LegacyMaterialMismatch};
  }
  for (int surfacePosition = 0; surfacePosition < activeSurfaceCount; ++surfacePosition) {
    if (mTrkParams[iteration].LayerxX0[surfacePosition] != stagedLayerMaterial[surfacePosition].xOverX0) {
      throw TraversalException{iteration, TraversalFailureReason::LegacyMaterialMismatch};
    }
  }

  // 2.6. Resolve and validate this iteration's normalized-measurement binding:
  // per-surface-position normalized SurfaceMeasurement span, entirely from
  // `orderedSurfaces` (already resolved and validated above) and the
  // already-loaded TimeFrame normalized frame / legacy compatibility
  // structures -- before any TimeFrame tracking state is touched.
  // `stagedLayerMeasurements` is kept local (not yet committed to the
  // mLayerMeasurements member) until every remaining fallible check in this
  // function has also succeeded -- see the final commit block below -- so a
  // later failure leaves mLayerMeasurements exactly as resetTraversalCache()
  // left it (reset/empty spans), never partially populated from a failed
  // iteration.
  std::vector<gsl::span<const SurfaceMeasurement>> stagedLayerMeasurements(static_cast<std::size_t>(activeSurfaceCount));
  for (int surfacePosition = 0; surfacePosition < activeSurfaceCount; ++surfacePosition) {
    const auto surfaceId = orderedSurfaces[surfacePosition];
    const auto measurements = mFrame->getNormalizedFrame().getSurfaceMeasurements(surfaceId);
    const auto& clusters = mScratch->getUnsortedClusters()[surfacePosition];
    const auto& hits = mScratch->getTrackingFrameInfoOnLayer(surfacePosition);
    if (measurements.size() != clusters.size() || hits.size() != clusters.size() ||
        measurements.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
    }
    for (size_t i = 0; i < measurements.size(); ++i) {
      const auto& measurement = measurements[i];
      if (measurement.surface != surfaceId || !measurement.cluster.isValid() ||
          measurement.cluster.index > static_cast<uint32_t>(std::numeric_limits<int>::max()) ||
          clusters[i].clusterId != static_cast<int>(i)) {
        throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
      }
      const int externalIndex = mScratch->getClusterExternalIndex(surfacePosition, static_cast<int>(i));
      if (externalIndex < 0 || static_cast<uint32_t>(externalIndex) != measurement.cluster.index) {
        throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
      }
    }
    const auto rofBoundaries = mScratch->getROFrameClusters(surfacePosition);
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
    stagedLayerMeasurements[surfacePosition] = measurements;
  }

  // 3. Bind and validate index-table configuration in a local value. The
  // scratch is not touched until this validation succeeds.
  std::array<IndexTableUtilsCore, 2> indexTableConfigByKind{};
  std::vector<IndexTableUtilsCore> stagedIndexTableConfigs(activeSurfaceCount);
  for (const auto kind : {SurfaceKind::Cylinder, SurfaceKind::Disk}) {
    bool present = false;
    for (const auto surface : orderedSurfaces) {
      present = present || layout.getSurface(surface).kind == kind;
    }
    if (!present) {
      continue;
    }
    if (bindIndexTableConfiguration(indexTableConfigByKind[kindIndex(kind)], mTrkParams[iteration],
                                    activeSurfaceCount, kind) != IndexTableConfigError::None) {
      throw TraversalException{iteration, TraversalFailureReason::InvalidIndexTableConfiguration};
    }
  }
  for (int position = 0; position < activeSurfaceCount; ++position) {
    stagedIndexTableConfigs[position] = indexTableConfigByKind[kindIndex(layout.getSurface(orderedSurfaces[position]).kind)];
  }

  // 4. LUT-reuse invariant: a non-FirstPass iteration will not have its
  // index-table configuration (re)committed or its LUT storage reallocated
  // by the scratch initialization below (PassFlags[FirstPass]
  // gate) -- whether or not RebuildClusterLUT resorts clusters into the
  // existing LUT using whatever configuration TimeFrame already owns (e.g.
  // configuration for such an iteration must already match the
  // owned one exactly, checked before TimeFrame is touched at all.
  if (!mTrkParams[iteration].PassFlags[IterationStep::FirstPass]) {
    for (int position = 0; position < activeSurfaceCount; ++position) {
      if (!indexTableConfigurationsMatch(stagedIndexTableConfigs[position], mScratch->getIndexTableUtils(position), activeSurfaceCount)) {
        throw TraversalException{iteration, TraversalFailureReason::IndexTableConfigurationMismatch};
      }
    }
  }

  // 5. Only now is the scratch touched: it receives an already-validated
  // configuration and the sparse plan's actual ordered ids. It never inspects
  // a layer count or detector ID.
  const auto transitionIds = mBinding->getGlobalTransitions();
  const auto cellIds = mBinding->getGlobalCells();
  mScratch->initialise(*mFrame, mTrkParams[iteration], activeSurfaceCount, iteration,
                       stagedIndexTableConfigs, layout, transitionIds, cellIds,
                       orderedSurfaces, stagedLayerMeasurements);

  // A sorted Cluster is a locator/navigation cache only. Validate each
  // enabled ROF that can participate in a configured transition after every
  // TimeFrame initialisation, including LUT reuse and non-FirstPass paths.
  // The spans remain local until this check and all subsequent kind setup
  // have succeeded, so a structural failure cannot publish traversal caches.
  std::vector<bool> candidateReachableLayers(static_cast<std::size_t>(activeSurfaceCount), false);
  for (const auto transitionId : transitionIds) {
    const auto& transition = layout.getTransition(transitionId);
    const auto fromSlot = mBinding->getOwnedSurfaceIndex(transition.from);
    const auto toSlot = mBinding->getOwnedSurfaceIndex(transition.to);
    if (!fromSlot || !toSlot || *fromSlot >= candidateReachableLayers.size() || *toSlot >= candidateReachableLayers.size()) {
      throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
    }
    candidateReachableLayers[*fromSlot] = true;
    candidateReachableLayers[*toSlot] = true;
  }
  for (int layer = 0; layer < activeSurfaceCount; ++layer) {
    if (!candidateReachableLayers[layer]) {
      continue;
    }
    const auto measurements = stagedLayerMeasurements[layer];
    const auto rofBoundaries = mScratch->getROFrameClusters(layer);
    const auto rofMask = mScratch->getROFViews(layer).mask;
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
      if (!mScratch->isROFEnabled(layer, rof)) {
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
            externalIndex < 0 ||
            static_cast<uint32_t>(externalIndex) != measurement.cluster.index) {
          throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
        }
      }
      if (std::find(seen.begin(), seen.end(), uint8_t{0}) != seen.end()) {
        throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
      }
    }
  }

  // 6. Validate the sparse plan/binding after scratch initialization and before
  // committing traversal state. No layer-indexed topology oracle is consulted.
  std::optional<SurfaceKind> activeKind;
  bool mixedKind = false;
  validateSparsePlan(iteration, layout, activeKind, mixedKind);
  (void)mixedKind;
  for (const auto transitionId : mBinding->getGlobalTransitions()) {
    const auto kind = layout.getSurface(layout.getTransition(transitionId).from).kind;
    mTransitionsByKind[kindIndex(kind)].push_back(transitionId);
  }
  for (const auto cellId : mBinding->getGlobalScheduledCells()) {
    const auto& cell = layout.getCell(cellId);
    const auto kind = layout.getSurface(layout.getTransition(cell.firstTransition).from).kind;
    mScheduledCellsByKind[kindIndex(kind)].push_back(cellId);
  }
  for (const auto cellId : mBinding->getGlobalRoadStartCells()) {
    const auto& cell = layout.getCell(cellId);
    const auto kind = layout.getSurface(layout.getTransition(cell.firstTransition).from).kind;
    mRoadStartCellsByKind[kindIndex(kind)].push_back(cellId);
  }

  // Bound from the still-local stagedLayerMaterial, not the mLayerMaterial
  // member (not committed yet): attachHitConfig.layerMaterial is rebound to
  // point at mLayerMaterial once that member is actually populated, at the
  // final commit below, before it escapes into mAttachHitConfig.
  auto attachHitConfig = bindAttachHitConfig(
    gsl::span<const NominalSurfaceMaterial>(stagedLayerMaterial.data(), stagedLayerMaterial.size()), mTrkParams[iteration]);
  if (!attachHitConfig.isValid(activeSurfaceCount)) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
  }
  const auto geometryConfig = bindLayerGeometryConfig(mTrkParams[iteration], attachHitConfig);
  if (!geometryConfig.isValid(activeSurfaceCount)) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
  }
  DiskReferenceCoordinateView referenceCoordinateView{};
  if (!mTransitionsByKind[1].empty()) {
    const auto legacyDiskReference = bindLegacyMFTReferenceCoordinates();
    mDiskLayerReferenceZStorage.assign(activeSurfaceCount, 0.f);
    std::size_t diskLayer = 0;
    for (int position = 0; position < activeSurfaceCount; ++position) {
      if (layout.getSurface(orderedSurfaces[position]).kind != SurfaceKind::Disk) {
        continue;
      }
      if (diskLayer >= legacyDiskReference.perLayerReferenceZ.size()) {
        throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
      }
      mDiskLayerReferenceZStorage[position] = legacyDiskReference.perLayerReferenceZ[diskLayer++];
    }
    referenceCoordinateView = {mDiskLayerReferenceZStorage};
    if (!referenceCoordinateView.isValid(activeSurfaceCount)) {
      throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
    }
  }
  if (!activeKind || (*activeKind != SurfaceKind::Cylinder && *activeKind != SurfaceKind::Disk)) {
    throw TraversalException{iteration, TraversalFailureReason::StateFamilyMismatch};
  }
  mKernelParameters = bindTrackingKernelParameters(mTrkParams[iteration]);
  if (!mKernelParameters.isValid()) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
  }

  mTraversalGraph = layout;
  mTraversalCacheValid = true;
  mDiskLayerReferenceZ = referenceCoordinateView.perLayerReferenceZ;
  // Commit the material cache itself only now, alongside every other
  // traversal cache, then rebind attachHitConfig's span off the
  // about-to-be-destroyed local stagedLayerMaterial and onto the
  // now-populated, function-lifetime-independent mLayerMaterial member
  // before mAttachHitConfig (which outlives this call) retains it.
  mLayerMaterial = std::move(stagedLayerMaterial);
  attachHitConfig.layerMaterial = gsl::span<const NominalSurfaceMaterial>(mLayerMaterial.data(), mLayerMaterial.size());
  mAttachHitConfig = attachHitConfig;
  // One-time normalized-measurement binding: committed here, alongside every
  // other successfully staged traversal cache -- see mLayerMeasurements' own
  // doc for the commit contract.
  mLayerMeasurements = std::move(stagedLayerMeasurements);

  // All fallible validation for this iteration (layout/grouping, legacy
  // parity, state-family, and every kind/geometry binding above) has now
  // succeeded. What follows is the relocated, total (non-throwing) per-layer
  // and per-transition scattering/bending preparation -- see the method doc
  // for the ordering/failure contract this relies on.
  if (!mTransitionsByKind[0].empty()) {
    prepareTransitionScatteringAndBendingForKind<SurfaceKind::Cylinder>(iteration, geometryConfig, referenceCoordinateView,
                                                                        mTransitionsByKind[0]);
  }
  if (!mTransitionsByKind[1].empty()) {
    prepareTransitionScatteringAndBendingForKind<SurfaceKind::Disk>(iteration, geometryConfig, referenceCoordinateView,
                                                                    mTransitionsByKind[1]);
  }
}

template <SurfaceKind Kind>
void TrackerTraits::prepareTransitionScatteringAndBendingForKind(
  int iteration,
  const LayerGeometryConfigView& geometryConfig,
  const DiskReferenceCoordinateView& referenceCoordinateView,
  gsl::span<const TransitionId> transitionIds)
{
  const auto& trkParam = mTrkParams[iteration];
  const auto& topology = mTraversalGraph;

  // Per-layer step: genuinely kind-specific (typed operation), but the
  // extent is the adopted plan, not the TrackerTraits compatibility width.
  const int activeSurfaceCount = static_cast<int>(mScratch->getNOwnedSurfaces());
  std::vector<float> msAngles(static_cast<std::size_t>(activeSurfaceCount));
  for (int iLayer{0}; iLayer < activeSurfaceCount; ++iLayer) {
    if constexpr (Kind == SurfaceKind::Cylinder) {
      msAngles[iLayer] = cylinderLayerMultipleScatteringAngle(
        CylinderLayerScatteringInputs{geometryConfig.layerMaterial[iLayer].xOverX0}, trkParam.TrackletMinPt);
    } else {
      msAngles[iLayer] = diskLayerMultipleScatteringAngle(
        DiskLayerScatteringInputs{geometryConfig.layerMaterial[iLayer].xOverX0, geometryConfig.layerRadii[iLayer],
                                  referenceCoordinateView.perLayerReferenceZ[iLayer]},
        trkParam.TrackletMinPt);
    }
  }

  // Per-transition step: shared/Kind-independent post-clamp arithmetic behind
  // a Kind-specific curvature clamp. Iterate the binding's validated sparse
  // transition order, which is the compact scratch order and preserves the
  // existing global topology ordering without consulting a layer topology.
  // The operation binding is installed after this preparation step. Resolve
  // the already-validated plan/binding span directly here; using the
  // not-yet-bound operation cache would silently leave every transition
  // preparation entry at its reset value.
  auto& transitionMSAngles = mScratch->getTransitionMSAngles();
  auto& transitionPhiCuts = mScratch->getTransitionPhiCuts();
  float oneOverR{0.001f * 0.3f * std::abs(getBz()) / trkParam.TrackletMinPt};
  for (const auto transitionId : transitionIds) {
    const auto transitionSlot = requireScratchTransitionSlot(iteration, transitionId);
    const auto& transition = topology.getTransition(transitionId);
    const int fromLayer = requireSurfacePosition(iteration, transition.from);
    const int toLayer = requireSurfacePosition(iteration, transition.to);
    const float r1 = trkParam.LayerRadii[fromLayer];
    const float r2 = trkParam.LayerRadii[toLayer];
    if constexpr (Kind == SurfaceKind::Cylinder) {
      oneOverR = clampCylinderTransitionCurvature(oneOverR, r2);
    } else {
      oneOverR = clampDiskTransitionCurvature(oneOverR, r2);
    }
    const float res1 = o2::gpu::CAMath::Hypot(trkParam.PVres, mScratch->getPositionResolution(fromLayer));
    const float res2 = o2::gpu::CAMath::Hypot(trkParam.PVres, mScratch->getPositionResolution(toLayer));
    const auto prep = prepareTransitionScatteringAndBending(
      gsl::span<const float>(msAngles.data(), msAngles.size()), fromLayer, toLayer, r1, r2, oneOverR, res1, res2);
    transitionMSAngles[transitionSlot] = prep.msAngle;
    transitionPhiCuts[transitionSlot] = prep.phiCut;
  }
}

void TrackerTraits::computeLayerTracklets(const int iteration, int iVertex)
{
  // Clear the scratch-owned compact transition storage before dispatch.
  const auto scratchTransitionCount = mScratch->getTracklets().size();
  for (size_t transitionId = 0; transitionId < scratchTransitionCount; ++transitionId) {
    mScratch->getTracklets()[transitionId].clear();
    mScratch->getTrackletsLabel(transitionId).clear();
    std::fill(mScratch->getTrackletsLookupTable()[transitionId].begin(), mScratch->getTrackletsLookupTable()[transitionId].end(), 0);
  }

  if (!mTraversalCacheValid) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidTraversalSchedule};
  }

  computeLayerTrackletsImpl(iteration, iVertex, mBinding->getGlobalTransitions());
}

void TrackerTraits::computeLayerTrackletsImpl(
  const int iteration,
  const int iVertex,
  gsl::span<const TransitionId> transitionIds)
{
  const auto& topology = mTraversalGraph;
  const Vertex diamondVert(mTrkParams[iteration].Diamond, mTrkParams[iteration].DiamondCov, 1, 1.f);

  mTaskArena->execute([&] {
    // Resolves one sparse SurfaceTransition's endpoints to runtime-plan
    // positions through the immutable binding. Called exactly once per
    // transitionId, outside every candidate (ROF/cluster/vertex) loop below.
    auto resolveTransitionLayers = [&](int transitionId) -> std::pair<int, int> {
      const auto& transition = topology.getTransition(TransitionId{static_cast<uint16_t>(transitionId)});
      return {requireSurfacePosition(iteration, transition.from),
              requireSurfacePosition(iteration, transition.to)};
    };

    auto makeTransitionState = [&](SurfaceKind kind, int transitionId, int fromLayer, int toLayer) {
      TrackletProjectionState state{};
      const auto layerRadii = gsl::span<const float>{mTrkParams[iteration].LayerRadii.data(),
                                                     mTrkParams[iteration].LayerRadii.size()};
      if (!bindTrackletProjectionState(kind, fromLayer, toLayer, layerRadii, mDiskLayerReferenceZ,
                                       mScratch->getMinR(toLayer), mScratch->getMaxR(toLayer),
                                       mScratch->getPositionResolution(fromLayer),
                                       mScratch->getTransitionMSAngle(transitionId),
                                       mScratch->getTransitionPhiCut(transitionId), state)) {
        throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
      }
      return state;
    };

    auto forTracklets = [&](auto Mode, int transitionId, int fromLayer, int toLayer, const auto& transitionState, int pivotROF, int base, int& offset) -> int {
      if (!mScratch->isROFEnabled(fromLayer, pivotROF)) {
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
        diamondForROF = diamondVertexForROF(diamondVert, mScratch->getROFViews(fromLayer).overlap,
                                            mScratch->getROFLocalLayer(fromLayer), pivotROF);
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

      const auto& rofOverlap = mScratch->getROFOverlap(fromLayer, toLayer, pivotROF);
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
        const o2::its::Cluster& currentCluster = layer0[iCluster];
        const int currentSortedIndex = mScratch->getSortedIndex(pivotROF, fromLayer, iCluster);
        if (mScratch->isClusterUsed(fromLayer, currentCluster.clusterId)) {
          continue;
        }
        const auto& sourceMeasurement = mLayerMeasurements[fromLayer][currentCluster.clusterId];

        for (int iV = startVtx; iV < endVtx; ++iV) {
          const auto& pv = primaryVertices[iV];
          if (!mScratch->isVertexCompatible(fromLayer, pivotROF, pv)) {
            continue;
          }
          if (pv.isFlagSet(Vertex::Flags::UPCMode) != mTrkParams[iteration].PassFlags[IterationStep::SelectUPCVertices]) {
            continue;
          }
          TrackletSearchWindow searchWindow{};
          const bool projected = projectTrackletSearchWindow(sourceMeasurement, currentCluster, pv,
                                                             transitionState, getBz(),
                                                             mScratch->getIndexTableUtils(toLayer),
                                                             mKernelParameters, searchWindow);
          if (!projected) {
            continue;
          }
          const auto bins = trackletSearchBins(searchWindow);
          const auto& indexTableUtils = mScratch->getIndexTableUtils(toLayer);
          const int rowBinsNum = trackletSearchRowCount(searchWindow, indexTableUtils);

          for (int targetROF = rofOverlap.getFirstEntry(); targetROF < rofOverlap.getEntriesBound(); ++targetROF) {
            if (!mScratch->isROFEnabled(toLayer, targetROF)) {
              continue;
            }
            auto layer1 = mScratch->getClustersOnLayer(targetROF, toLayer);
            if (layer1.empty()) {
              continue;
            }
            const auto ts = mScratch->getROFTimeStamp(fromLayer, pivotROF, toLayer, targetROF);
            if (!ts.isCompatible(pv.getTimeStamp())) {
              continue;
            }
            const auto& targetIndexTable = mScratch->getIndexTable(targetROF, toLayer);
            const int colBinRange = (bins.z - bins.x) + 1;
            for (int iRow = 0; iRow < rowBinsNum; ++iRow) {
              const int iRowBin = trackletSearchRowBin(searchWindow, iRow, indexTableUtils);
              if (iRowBin < 0) {
                break;
              }
              const int firstBinIdx = indexTableUtils.getBinIndex(bins.x, iRowBin);
              const int maxBinIdx = firstBinIdx + colBinRange;
              const int firstRow = targetIndexTable[firstBinIdx];
              const int lastRow = targetIndexTable[maxBinIdx];
              for (int iNext = firstRow; iNext < lastRow; ++iNext) {
                if (iNext >= int(layer1.size())) {
                  break;
                }
                const o2::its::Cluster& nextCluster = layer1[iNext];
                if (mScratch->isClusterUsed(toLayer, nextCluster.clusterId)) {
                  continue;
                }
                const auto& targetMeasurement = mLayerMeasurements[toLayer][nextCluster.clusterId];

                float tanL = 0.f;
                const bool accepted = acceptTrackletCandidate(searchWindow, sourceMeasurement, currentCluster,
                                                              targetMeasurement, nextCluster, tanL);
                if (accepted) {
                  const float phi{o2::gpu::GPUCommonMath::ATan2(sourceMeasurement.global.y - targetMeasurement.global.y,
                                                                sourceMeasurement.global.x - targetMeasurement.global.x)};
                  if constexpr (decltype(Mode)::value == PassMode::OnePass::value) {
                    tracklets.emplace_back(currentSortedIndex, mScratch->getSortedIndex(targetROF, toLayer, iNext), tanL, phi, ts);
                  } else if constexpr (decltype(Mode)::value == PassMode::TwoPassCount::value) {
                    ++localCount;
                  } else if constexpr (decltype(Mode)::value == PassMode::TwoPassInsert::value) {
                    const int idx = base + offset++;
                    tracklets[idx] = o2::its::Tracklet(currentSortedIndex, mScratch->getSortedIndex(targetROF, toLayer, iNext), tanL, phi, ts);
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
        const int scratchTransitionId = requireScratchTransitionSlot(iteration, typedTransitionId);
        const auto [fromLayer, toLayer] = resolveTransitionLayers(transitionId);
        const auto kind = topology.getSurface(topology.getTransition(typedTransitionId).from).kind;
        const auto transitionState = makeTransitionState(kind, scratchTransitionId, fromLayer, toLayer);
        const int startROF = 0, endROF = mScratch->getROFTiming(fromLayer).mNROFsTF;
        for (int pivotROF{startROF}; pivotROF < endROF; ++pivotROF) {
          forTracklets(PassMode::OnePass{}, scratchTransitionId, fromLayer, toLayer, transitionState, pivotROF, 0, dummy);
        }
      }
    } else {
      tbb::parallel_for(0, static_cast<int>(transitionIds.size()), [&](const int transitionIndex) {
        const auto typedTransitionId = transitionIds[transitionIndex];
        const int transitionId = typedTransitionId.value();
        const int scratchTransitionId = requireScratchTransitionSlot(iteration, typedTransitionId);
        const auto [fromLayer, toLayer] = resolveTransitionLayers(transitionId);
        const auto kind = topology.getSurface(topology.getTransition(typedTransitionId).from).kind;
        const auto transitionState = makeTransitionState(kind, scratchTransitionId, fromLayer, toLayer);
        const int startROF = 0, endROF = mScratch->getROFTiming(fromLayer).mNROFsTF;
        bounded_vector<int> perROFCount((endROF - startROF) + 1, mMemoryPool.get());
        tbb::parallel_for(startROF, endROF, [&](const int pivotROF) {
          perROFCount[pivotROF - startROF] = forTracklets(PassMode::TwoPassCount{}, scratchTransitionId, fromLayer, toLayer, transitionState, pivotROF, 0, dummy);
        });
        std::exclusive_scan(perROFCount.begin(), perROFCount.end(), perROFCount.begin(), 0);
        const int nTracklets = perROFCount.back();
        mScratch->getTracklets()[scratchTransitionId].resize(nTracklets);
        if (nTracklets == 0) {
          return;
        }
        tbb::parallel_for(startROF, endROF, [&](const int pivotROF) {
          int baseIdx = perROFCount[pivotROF - startROF];
          if (baseIdx == perROFCount[pivotROF + 1 - startROF]) {
            return;
          }
          int localIdx = 0;
          forTracklets(PassMode::TwoPassInsert{}, scratchTransitionId, fromLayer, toLayer, transitionState, pivotROF, baseIdx, localIdx);
        });
      });
    }

    tbb::parallel_for(0, static_cast<int>(transitionIds.size()), [&](const int transitionIndex) {
      const int scratchTransitionId = requireScratchTransitionSlot(iteration, transitionIds[transitionIndex]);
      /// Sort tracklets & remove duplicates
      // duplicates can exist simply since we evaluate per vertex
      auto& trkl{mScratch->getTracklets()[scratchTransitionId]};
      std::sort(trkl.begin(), trkl.end());
      trkl.erase(std::unique(trkl.begin(), trkl.end()), trkl.end());
      trkl.shrink_to_fit();
      auto& lut{mScratch->getTrackletsLookupTable()[scratchTransitionId]};
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
        const auto typedTransitionId = transitionIds[transitionIndex];
        const int transitionId = typedTransitionId.value();
        const int scratchTransitionId = requireScratchTransitionSlot(iteration, typedTransitionId);
        const auto [fromLayer, toLayer] = resolveTransitionLayers(transitionId);
        for (auto& trk : mScratch->getTracklets()[scratchTransitionId]) {
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
          mScratch->getTrackletsLabel(scratchTransitionId).emplace_back(label);
        }
      });
    }
  });
}

void TrackerTraits::computeLayerCells(const int iteration)
{
  // Defensive size-consistency check against this scratch's own
  // already-allocated compact cell count, never topology.nCells directly --
  // see computeLayerTracklets()'s identical reasoning for
  // scratchTransitionCount above. The two parallel per-cell-topology
  // containers must still agree with each other before the clear loop
  // touches either.
  const auto scratchCellCount = mScratch->getCells().size();
  if (mScratch->getCellsLookupTable().size() != scratchCellCount) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  // This clear/allocation pass runs once over the full sparse cell count.
  for (size_t cellTopologyId = 0; cellTopologyId < scratchCellCount; ++cellTopologyId) {
    deepVectorClear(mScratch->getCells()[cellTopologyId]);
    deepVectorClear(mScratch->getCellsLookupTable()[cellTopologyId]);
    if (mScratch->hasMCinformation() && mTrkParams[iteration].CreateArtefactLabels) {
      deepVectorClear(mScratch->getCellsLabel(cellTopologyId));
    }
  }

  if (!mTraversalCacheValid) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidTraversalSchedule};
  }
  if (!mAttachHitConfig.isValid(static_cast<int>(mScratch->getNOwnedSurfaces()))) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
  }

  computeLayerCellsImpl(iteration, mBinding->getGlobalCells());

  const auto scratchTransitionCount = mScratch->getTracklets().size();
  for (size_t transitionId = 0; transitionId < scratchTransitionCount; ++transitionId) {
    deepVectorClear(mScratch->getTracklets()[transitionId]);
    deepVectorClear(mScratch->getTrackletsLabel(transitionId));
  }
}

void TrackerTraits::computeLayerCellsImpl(
  const int iteration,
  gsl::span<const CellTopologyId> cellIds)
{
  // The defensive size-consistency check for mScratch's per-cellTopologyId
  // storage lives in computeLayerCells(), before this traversal begins.
  const auto& topology = mTraversalGraph;

  mTaskArena->execute([&] {
    struct CellHitBinding {
      std::array<int, 3> layers{};
      std::array<SurfaceDescriptor, 3> surfaces{};
    };

    // Resolves one sparse SurfaceCellTopology's three hit surfaces to
    // runtime-plan positions through the immutable binding. Called exactly
    // once per CellTopologyId, outside the per-tracklet candidate loop.
    auto resolveCellHitBinding = [&](const auto& cellTopology) -> CellHitBinding {
      const auto& firstTransition = topology.getTransition(cellTopology.firstTransition);
      const auto& secondTransition = topology.getTransition(cellTopology.secondTransition);
      const std::array<SurfaceId, 3> surfaces{firstTransition.from, firstTransition.to, secondTransition.to};
      CellHitBinding binding;
      for (int i = 0; i < 3; ++i) {
        const auto surfaceId = surfaces[i];
        binding.layers[i] = requireSurfacePosition(iteration, surfaceId);
        binding.surfaces[i] = topology.getSurface(surfaceId);
      }
      return binding;
    };

    auto forTrackletCells = [&](auto Mode, SurfaceKind kind, int firstTransitionId, int secondTransitionId, const CellHitBinding& hitBinding, bounded_vector<CellSeed>& layerCells, int iTracklet, int offset = 0) -> int {
      const auto& hitLayers = hitBinding.layers;
      const o2::its::Tracklet& currentTracklet{mScratch->getTracklets()[firstTransitionId][iTracklet]};
      const int nextLayerClusterIndex{currentTracklet.secondClusterIndex};
      const int nextLayerFirstTrackletIndex{mScratch->getTrackletsLookupTable()[secondTransitionId][nextLayerClusterIndex]};
      const int nextLayerLastTrackletIndex{mScratch->getTrackletsLookupTable()[secondTransitionId][nextLayerClusterIndex + 1]};
      int foundCells{0};
      for (int iNextTracklet{nextLayerFirstTrackletIndex}; iNextTracklet < nextLayerLastTrackletIndex; ++iNextTracklet) {
        const o2::its::Tracklet& nextTracklet{mScratch->getTracklets()[secondTransitionId][iNextTracklet]};
        if (nextTracklet.firstClusterIndex != nextLayerClusterIndex) {
          break;
        }
        if (!currentTracklet.getTimeStamp().isCompatible(nextTracklet.getTimeStamp())) {
          continue;
        }

        /// Track seed preparation. Clusters are numbered progressively from the innermost going outward.
        const int clusId[3]{
          mScratch->getClusters()[hitLayers[0]][currentTracklet.firstClusterIndex].clusterId,
          mScratch->getClusters()[hitLayers[1]][nextTracklet.firstClusterIndex].clusterId,
          mScratch->getClusters()[hitLayers[2]][nextTracklet.secondClusterIndex].clusterId};

        const auto& measurementInner = mLayerMeasurements[hitLayers[0]][clusId[0]];
        const auto& measurementMiddle = mLayerMeasurements[hitLayers[1]][clusId[1]];
        const auto& measurementOuter = mLayerMeasurements[hitLayers[2]][clusId[2]];
        std::array<DirectionObservation, 3> directionObservations{};
        if (!makeDirectionObservation(hitBinding.surfaces[0], measurementInner, directionObservations[0]) ||
            !makeDirectionObservation(hitBinding.surfaces[1], measurementMiddle, directionObservations[1]) ||
            !makeDirectionObservation(hitBinding.surfaces[2], measurementOuter, directionObservations[2])) {
          continue;
        }
        CellDirectionCompatibility directionCompatibility{};
        // The outgoing transition begins at the middle hit; its precomputed
        // TrackletMinPt scattering is the equivalent thin-kick uncertainty.
        const double transitionMSAngle = static_cast<double>(mScratch->getTransitionMSAngle(secondTransitionId));
        const DirectionProcessNoise directionProcessNoise{transitionMSAngle * transitionMSAngle};
        if (cellDirectionsAreCompatible(directionObservations, directionProcessNoise, mKernelParameters.nSigmaCut,
                                        directionCompatibility)) {

          // Strictly {inner, middle, outer}: Cylinder reads [1] then
          // [0] (outer slot unused), Disk reads [2], [1], [0].
          const std::array<NominalSurfaceMaterial, 3> material{
            mAttachHitConfig.layerMaterial[hitLayers[0]],
            mAttachHitConfig.layerMaterial[hitLayers[1]],
            mAttachHitConfig.layerMaterial[hitLayers[2]]};

          SurfaceKinematicState state{};
          float chi2{0.f};
          OperationFailureReason buildReason{};
          const bool good = buildCellSeed(kind, measurementInner, measurementMiddle, measurementOuter,
                                          material, getBz(), kCompatibilityAbsCharge, kCompatibilityPID,
                                          state, chi2, mKernelParameters, buildReason);

          if (good) {
            TimeEstBC ts = currentTracklet.getTimeStamp();
            ts += nextTracklet.getTimeStamp();
            // Reconstructed directly from the three already-resolved plan
            // positions via LayerMask's existing
            // 3-int constructor. The positions are validated against the
            // sparse cell's hit-surface mask during plan validation, so no
            // generic SurfaceMask->LayerMask converter is needed.
            const LayerMask hitLayerMask{hitLayers[0], hitLayers[1], hitLayers[2]};
            if constexpr (decltype(Mode)::value == PassMode::OnePass::value) {
              layerCells.emplace_back(hitLayerMask, clusId[0], clusId[1], clusId[2], iTracklet, iNextTracklet, state, chi2, ts);
              ++foundCells;
            } else if constexpr (decltype(Mode)::value == PassMode::TwoPassCount::value) {
              ++foundCells;
            } else if constexpr (decltype(Mode)::value == PassMode::TwoPassInsert::value) {
              layerCells[offset++] = CellSeed(hitLayerMask, clusId[0], clusId[1], clusId[2], iTracklet, iNextTracklet, state, chi2, ts);
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
      // Sole translation point for this loop iteration: every scratch access
      // below (getCells(), getCellsLookupTable(),
      // getCellsLabel(), getTracklets(), getTrackletsLookupTable(),
      // getTrackletsLabel()) is addressed exclusively through these three
      // already-translated compact slots, never through typedCellId/
      // cellTopology.firstTransition/secondTransition's own .value() again.
      const int cellTopologyId = requireScratchCellSlot(iteration, typedCellId);
      const auto& cellTopology = topology.getCell(typedCellId);
      const auto kind = topology.getSurface(topology.getTransition(cellTopology.firstTransition).from).kind;
      const int firstTransitionId = requireScratchTransitionSlot(iteration, cellTopology.firstTransition);
      const int secondTransitionId = requireScratchTransitionSlot(iteration, cellTopology.secondTransition);
      if (mScratch->getTracklets()[firstTransitionId].empty() ||
          mScratch->getTracklets()[secondTransitionId].empty()) {
        continue;
      }
      const auto hitBinding = resolveCellHitBinding(cellTopology);

      auto& layerCells = mScratch->getCells()[cellTopologyId];
      const int currentLayerTrackletsNum{static_cast<int>(mScratch->getTracklets()[firstTransitionId].size())};
      bounded_vector<int> perTrackletCount(currentLayerTrackletsNum + 1, 0, mMemoryPool.get());
      if (mTaskArena->max_concurrency() <= 1) {
        for (int iTracklet{0}; iTracklet < currentLayerTrackletsNum; ++iTracklet) {
          perTrackletCount[iTracklet] = forTrackletCells(PassMode::OnePass{}, kind, firstTransitionId, secondTransitionId, hitBinding, layerCells, iTracklet);
        }
        std::exclusive_scan(perTrackletCount.begin(), perTrackletCount.end(), perTrackletCount.begin(), 0);
      } else {
        tbb::parallel_for(0, currentLayerTrackletsNum, [&](const int iTracklet) {
          perTrackletCount[iTracklet] = forTrackletCells(PassMode::TwoPassCount{}, kind, firstTransitionId, secondTransitionId, hitBinding, layerCells, iTracklet);
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
          forTrackletCells(PassMode::TwoPassInsert{}, kind, firstTransitionId, secondTransitionId, hitBinding, layerCells, iTracklet, offset);
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

void TrackerTraits::findCellsNeighbours(const int iteration)
{
  if (!mTraversalCacheValid) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidTraversalSchedule};
  }
  for (std::size_t slot = 0; slot < mScratch->getCellsNeighbours().size(); ++slot) {
    deepVectorClear(mScratch->getCellsNeighbours()[slot]);
    deepVectorClear(mScratch->getCellsNeighboursTopology()[slot]);
    deepVectorClear(mScratch->getCellsNeighboursLUT()[slot]);
  }
  if (!mScheduledCellsByKind[0].empty()) {
    findCellsNeighboursCylinder(iteration);
  }
  if (!mScheduledCellsByKind[1].empty()) {
    findCellsNeighboursDisk(iteration);
  }
  for (auto& cellLUT : mScratch->getCellsLookupTable()) {
    deepVectorClear(cellLUT);
  }
}

template <SurfaceKind Kind>
void TrackerTraits::findCellsNeighboursForKind(
  int iteration,
  gsl::span<const CellTopologyId> scheduledCells,
  const TrackingKernelParameters& params)
{
  const auto topology = mTraversalGraph;
  // Check against this scratch's allocated compact cell count.
  const auto scratchCellCount = mScratch->getCells().size();
  if (mScratch->getCellsLookupTable().size() != scratchCellCount ||
      mScratch->getCellsNeighbours().size() != scratchCellCount ||
      mScratch->getCellsNeighboursTopology().size() != scratchCellCount ||
      mScratch->getCellsNeighboursLUT().size() != scratchCellCount) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  mTaskArena->execute([&] {
    std::vector<bounded_vector<CellNeighbour>> cellsNeighboursByTarget;
    cellsNeighboursByTarget.reserve(scratchCellCount);
    for (size_t cellTopologyId = 0; cellTopologyId < scratchCellCount; ++cellTopologyId) {
      cellsNeighboursByTarget.emplace_back(mMemoryPool.get());
    }

    for (const auto scheduledId : scheduledCells) {
      // Sole translation point for this loop iteration: every
      // mScratch->getCells...() access below, for this
      // source cell, is addressed through this one already-translated slot.
      const int cellTopologyId = requireScratchCellSlot(iteration, scheduledId);
      if (static_cast<size_t>(cellTopologyId) >= scratchCellCount ||
          static_cast<size_t>(cellTopologyId) >= mScratch->getCellsLookupTable().size()) {
        throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
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
          // Dynamically discovered neighbours are translated here;
          // `nextTopologyId` is read from the global
          // topology's own CSR array (topology.cellsByFirstTransition), not
          // from any precomputed per-kind span, so it is translated to its
          // compact scratch slot here, exactly once, before any scratch
          // access below uses it.
          const auto nextTopologyId = topology.cellsByFirstTransition[successors.getFirstEntry() + iSuccessor];
          const int nextCellTopologyId = requireScratchCellSlot(iteration, nextTopologyId);
          if (static_cast<size_t>(nextCellTopologyId) >= mScratch->getCells().size() ||
              static_cast<size_t>(nextCellTopologyId) >= mScratch->getCellsLookupTable().size()) {
            throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
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
            throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
          }
          for (int iNextCell{nextLayerFirstCellIndex}; iNextCell < nextLayerLastCellIndex; ++iNextCell) {
            const auto& nextCellSeedRef{mScratch->getCells()[nextCellTopologyId][iNextCell]};
            if (nextCellSeedRef.getFirstTrackletIndex() != nextLayerTrackletIndex || !currentCellSeed.getTimeStamp().isCompatible(nextCellSeedRef.getTimeStamp())) {
              break;
            }

            const bool compatible = [&] {
              if constexpr (Kind == SurfaceKind::Cylinder) {
                return cellsCylinderAreCompatible(currentCellSeed.state(), nextCellSeedRef.state(),
                                                  currentCellSeed.getSecondClusterIndex(), nextCellSeedRef.getFirstClusterIndex(),
                                                  getBz(), params);
              } else {
                return cellsDiskAreCompatible(currentCellSeed.state(), nextCellSeedRef.state(),
                                              currentCellSeed.getSecondClusterIndex(), nextCellSeedRef.getFirstClusterIndex(),
                                              getBz(), params);
              }
            }();
            if (!compatible) {
              continue;
            }

            const int nextLevel = currentCellSeed.getLevel() + 1;
            localNeighbours.emplace_back(cellTopologyId, iCell, nextCellTopologyId, iNextCell, nextLevel);
          }
        }
      });

      bounded_vector<size_t> count(scratchCellCount, 0, mMemoryPool.get());
      for (const auto& localNeighbours : sourceNeighbours) {
        for (const auto& neigh : localNeighbours) {
          ++count[neigh.nextCellTopology];
        }
      }
      for (size_t i{0}; i < scratchCellCount; ++i) {
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

    for (size_t cellTopologyId = 0; cellTopologyId < scratchCellCount; ++cellTopologyId) {
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

  });
}

template <SurfaceKind Kind, typename InputSeed>
void TrackerTraits::processNeighbours(int iteration, int defaultCellTopologyId, int iLevel, const bounded_vector<InputSeed>& currentCellSeed, const bounded_vector<int>& currentCellId, const bounded_vector<int>& currentCellTopologyId, bounded_vector<TrackSeed>& updatedCellSeeds, bounded_vector<int>& updatedCellsIds, bounded_vector<int>& updatedCellsTopologyIds, const TrackingKernelParameters& params)
{
  const auto layerMaterial = mAttachHitConfig.layerMaterial;
  const int activeSurfaceCount = static_cast<int>(mScratch->getNOwnedSurfaces());

  mTaskArena->execute([&] {
    auto forCellNeighbours = [&](auto Mode, int iCell, int offset = 0) -> int {
      const auto& currentCell{currentCellSeed[iCell]};
      const int cellTopologyId = currentCellTopologyId.empty() ? defaultCellTopologyId : currentCellTopologyId[iCell];

      if constexpr (decltype(Mode)::value != PassMode::TwoPassInsert::value) {
        if (currentCell.getLevel() != iLevel) {
          return 0;
        }
        if (currentCellId.empty()) {
          for (int layer = 0; layer < activeSurfaceCount; ++layer) {
            const int clusterIndex = currentCell.getCluster(layer);
            if (clusterIndex != o2::its::constants::UnusedIndex && mScratch->isClusterUsed(layer, clusterIndex)) {
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
        TrackSeed seed{currentCell};
        seed.getTimeStamp() = currentCell.getTimeStamp();
        seed.getTimeStamp() += neighbourCell.getTimeStamp();

        const auto& measurement = mLayerMeasurements[neighbourLayer][neighbourCluster];
        float chi2 = seed.getChi2();
        OperationFailureReason attachReason{};
        const bool attached = [&] {
          if constexpr (Kind == SurfaceKind::Cylinder) {
            return attachCylinderHit(seed.state(), measurement, layerMaterial[neighbourLayer], getBz(), chi2, params, attachReason);
          } else {
            return attachDiskHit(seed.state(), measurement, layerMaterial[neighbourLayer], getBz(), chi2, params, attachReason);
          }
        }();
        if (!attached) {
          continue;
        }
        seed.setChi2(chi2);

        if constexpr (decltype(Mode)::value != PassMode::TwoPassCount::value) {
          seed.setCluster(neighbourLayer, neighbourCluster);
          if (neighbourLayer < 0 || neighbourLayer >= activeSurfaceCount) {
            throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
          }
          // TrackSeed::SurfaceMask is the fixed-capacity compact plan
          // position space used by the CA acceptance/refit loops. Global
          // SurfaceIds stay on normalized measurements and CommonTrack
          // references; mixing them into this local mask would make a
          // combined sparse binding look like extra holes. The binding's
          // ordered position is already `neighbourLayer` here.
          auto surfaceMask = seed.getSurfaceMask();
          surfaceMask.set(SurfaceId{static_cast<uint16_t>(neighbourLayer)});
          seed.setSurfaceMask(surfaceMask);
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

void TrackerTraits::findRoads(const int iteration, SeedRefitFunction refitFunction)
{
  if (!mTraversalCacheValid) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidTraversalSchedule};
  }
  // Defensive sparse-plan check:
  // findRoadsForKind() below indexes the scratch cell vectors with the
  // compact slots returned by the bound sparse plan. The count is checked
  // before indexing so a future plan/storage desync becomes an explicit
  // failure instead of an out-of-bounds read; no legacy topology parity is
  // required by this path.
  // The expected count is the binding's own owned-cell count.
  if (mScratch->getCells().size() != mBinding->getGlobalCells().size()) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  if (!mRoadStartCellsByKind[0].empty()) {
    findRoadsCylinder(iteration, refitFunction);
  }
  if (!mRoadStartCellsByKind[1].empty()) {
    findRoadsDisk(iteration, refitFunction);
  }
}

template <SurfaceKind Kind>
void TrackerTraits::findRoadsForKind(const int iteration,
                                     const TrackingKernelParameters& params,
                                     SeedRefitFunction refitFunction,
                                     gsl::span<const CellTopologyId> roadStartCells,
                                     ClusterSourceId expectedSource)
{
  const int activeSurfaceCount = static_cast<int>(mScratch->getNOwnedSurfaces());
  bounded_vector<bounded_vector<int>> firstClusters(activeSurfaceCount, bounded_vector<int>(mMemoryPool.get()), mMemoryPool.get());
  firstClusters.resize(activeSurfaceCount);
  // Road starts are the binding's deterministic sparse-plan subsequence whose
  // endpoint is seeding-eligible. CellTopologyId values are resolved through
  // compact slots; SurfaceId is never used as a vector index here.
  // Road-length filter bound: maximum absolute q/pT, in the same (GeV/c)^-1
  // units as SurfaceKinematicState::parameters[4]. Applied identically to
  // both families via getQOverPt() (raw signed value, never squared); no
  // layer-count/family branch. std::abs() of a NaN/+-Inf q/pT is
  // never <= this finite bound (standard IEEE-754 comparison semantics), so
  // non-finite seeds are rejected deterministically without extra checks --
  // see testCellRepresentation.cxx's dedicated road-filter tests for focused
  // coverage of that behavior, per the correction's requirement not to rely
  // on it without a test.
  constexpr float maxAbsQOverPt = 1.e3f;
  const int cellsPerRoad = static_cast<int>(mScratch->getNOwnedSurfaces()) - 2;
  for (int startLevel{cellsPerRoad}; startLevel >= mTrkParams[iteration].CellMinimumLevel(); --startLevel) {

    auto seedFilter = [&](const auto& seed) {
      return seed.getHitLayerMask().isAllowed(mTrkParams[iteration].MaxHoles, mTrkParams[iteration].HoleLayerMask) &&
             seed.getHitLayerMask().length() >= mTrkParams[iteration].MinTrackLength &&
             std::abs(seed.getQOverPt()) <= maxAbsQOverPt && seed.getChi2() <= mTrkParams[iteration].MaxChi2NDF * ((startLevel + 2) * 2 - 5);
    };

    bounded_vector<TrackSeed> trackSeeds(mMemoryPool.get());
    // The binding supplies the ownership-filtered road-start span.
    for (const auto startId : roadStartCells) {
      // Translate the road-start cell to its compact slot once.
      const int startCellTopologyId = requireScratchCellSlot(iteration, startId);
      // Cell population is per-event/per-vertex data, never cached in
      // SurfacePlanBinding: this check must stay here, evaluated at
      // runtime against the current iVertex's TimeFrame content.
      if (mScratch->getCells()[startCellTopologyId].empty()) {
        continue;
      }

      bounded_vector<int> lastCellId(mMemoryPool.get()), updatedCellId(mMemoryPool.get());
      bounded_vector<int> lastCellTopologyId(mMemoryPool.get()), updatedCellTopologyId(mMemoryPool.get());
      bounded_vector<TrackSeed> lastCellSeed(mMemoryPool.get()), updatedCellSeed(mMemoryPool.get());

      processNeighbours<Kind>(iteration, startCellTopologyId, startLevel, mScratch->getCells()[startCellTopologyId], lastCellId, lastCellTopologyId, updatedCellSeed, updatedCellId, updatedCellTopologyId, params);

      int level = startLevel;
      while (level > 2 && !updatedCellSeed.empty()) {
        lastCellSeed.swap(updatedCellSeed);
        lastCellId.swap(updatedCellId);
        lastCellTopologyId.swap(updatedCellTopologyId);
        deepVectorClear(updatedCellSeed); /// tame the memory peaks
        deepVectorClear(updatedCellId);   /// tame the memory peaks
        deepVectorClear(updatedCellTopologyId);
        processNeighbours<Kind>(iteration, o2::its::constants::UnusedIndex, --level, lastCellSeed, lastCellId, lastCellTopologyId, updatedCellSeed, updatedCellId, updatedCellTopologyId, params);
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

    bounded_vector<TrackingCandidate> tracks(mMemoryPool.get());
    mTaskArena->execute([&] {
      auto forSeed = [&](auto Mode, int iSeed, int offset = 0) {
        TrackingCandidate temporaryTrack;
        temporaryTrack.seed = trackSeeds[iSeed];
        const bool refitSuccess = refitFunction(trackSeeds[iSeed],
                                                mTrkParams[iteration],
                                                mBz,
                                                *mScratch,
                                                mLayerMeasurements,
                                                mTraversalGraph.getSurfaceCatalogView(),
                                                expectedSource,
                                                temporaryTrack);
        if (refitSuccess) {
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
    std::sort(tracks.begin(), tracks.end(), [](const TrackingCandidate& a, const TrackingCandidate& b) {
      const auto ncla = a.getNumberOfClusters();
      const auto nclb = b.getNumberOfClusters();
      return (ncla == nclb) ? (a.track.chi2 < b.track.chi2) : ncla > nclb;
    });
    acceptTracks(iteration, tracks, firstClusters);
  }
}

void TrackerTraits::acceptTracks(int iteration,
                                 bounded_vector<TrackingCandidate>& tracks,
                                 bounded_vector<bounded_vector<int>>& firstClusters)
{
  auto& trks = acceptedTracksForSharedStatus();
  trks.reserve(trks.size() + tracks.size());
  const int activeSurfaceCount = static_cast<int>(mScratch->getNOwnedSurfaces());
  for (auto& track : tracks) {
    int nShared = 0;
    bool isFirstShared{false};
    int firstLayer{-1}, firstCluster{-1};
    for (int iLayer{0}; iLayer < activeSurfaceCount; ++iLayer) {
      if (track.getClusterIndex(iLayer) == o2::its::constants::UnusedIndex) {
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
    float smallestROFHalf = std::numeric_limits<float>::max();
    for (int iLayer{0}; iLayer < activeSurfaceCount; ++iLayer) {
      if (track.getClusterIndex(iLayer) == o2::its::constants::UnusedIndex) {
        continue;
      }
      smallestROFHalf = std::min(smallestROFHalf, mScratch->getROFTiming(iLayer).mROFLength * 0.5f);
      mScratch->markUsedCluster(iLayer, track.getClusterIndex(iLayer));
      int currentROF = mScratch->getClusterROF(iLayer, track.getClusterIndex(iLayer));
      const auto nominalROFTS = mScratch->getROFTiming(iLayer).getROFTimeBounds(currentROF);
      const auto expandedROFTS = mScratch->getROFTiming(iLayer).getROFTimeBounds(currentROF, true);
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
    const auto selectedTimestampSymmetric = selectedTimestamp.makeSymmetrical();
    // This is the same sanity clamp as the legacy symmetric timestamp, but
    // committed directly to the detector-neutral CommonTrack interval.
    const float selectedTimestampError = std::min(selectedTimestampSymmetric.getTimeStampError(), smallestROFHalf);
    track.track.timestamp = {static_cast<TFBC>(selectedTimestampSymmetric.getTimeStamp() - selectedTimestampError),
                             static_cast<TFBC>(selectedTimestampSymmetric.getTimeStamp() + selectedTimestampError)};
    trks.emplace_back(track);

    // Generic owner-thread publication. The candidate/refit loops above may
    // be parallel, but acceptTracks() runs after they have joined and is the
    // one deterministic accepted-result order. Typed sidecars are completed
    // later by the application adapter from the same generic results.
    CommonTrackShadowRecord shadow;
    if (!makeCandidateShadow(track, mLayerMeasurements, shadow)) {
      LOGP(fatal, "CommonTrack publication failed for an accepted CA track");
    }
    const auto commonTrackIndex = publishCommonTrackShadow(*mFrame, shadow);
    if (!commonTrackIndex) {
      LOGP(fatal, "CommonTrack publication failed for an accepted CA track");
    }
    trks.back().commonTrackIndex = *commonTrackIndex;

    if (mTrkParams[iteration].AllowSharingFirstCluster) {
      firstClusters[firstLayer].push_back(firstCluster);
    }
  }
}

bounded_vector<TrackingCandidate>& TrackerTraits::acceptedTracksForSharedStatus()
{
  return mAcceptedTracksForSharedStatus;
}

void TrackerTraits::clearAcceptedTracksForSharedStatus()
{
  deepVectorClear(mAcceptedTracksForSharedStatus, mMemoryPool.get());
}

void TrackerTraits::setBz(float bz)
{
  mBz = bz;
  mFrame->setBz(bz);
}

void TrackerTraits::setNThreads(int n, std::shared_ptr<tbb::task_arena>& arena)
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

} // namespace o2::itsmft::tracking

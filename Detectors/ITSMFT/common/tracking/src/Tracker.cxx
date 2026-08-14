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
/// \file Tracker.cxx
/// \brief
///

#include "ITSMFTTracking/Tracker.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <ranges>
#include <utility>

#include "Framework/Logger.h"
#include "GPUCommonMath.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITSMFTTracking/IndexTableConfiguration.h"
#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/detail/TrackerTraversalPreparation.h"

namespace o2::itsmft::tracking
{

namespace
{
constexpr std::size_t kindIndex(SurfaceKind kind) noexcept
{
  return kind == SurfaceKind::Cylinder ? 0u : 1u;
}

#define BIND_TRAVERSAL_CONTEXT(context)                                         \
  auto* mScratch = &(context).scratch;                                          \
  auto* mFrame = &(context).frame;                                              \
  const auto& mMemoryPool = mScratch->getMemoryPool();                          \
  const auto* mBinding = &(context).binding;                                    \
  const auto mTrkParams = (context).parameters;                                 \
  const auto mBz = (context).bz;                                                \
  const auto& mTraversalGraph = (context).graph;                                \
  auto& mTraversalCacheValid = (context).workspace.valid;                       \
  auto& mKernelParameters = (context).workspace.kernelParameters;               \
  auto& mDiskLayerReferenceZ = (context).workspace.diskLayerReferenceZView;     \
  auto& mDiskLayerReferenceZStorage = (context).workspace.diskLayerReferenceZ;  \
  auto& mAttachHitConfig = (context).workspace.attachHitConfig;                 \
  auto& mLayerMaterial = (context).workspace.layerMaterial;                     \
  auto& mLayerMeasurements = (context).workspace.layerMeasurements;             \
  auto& mLayerGlobalMeasurements = (context).workspace.layerGlobalMeasurements; \
  auto& mAcceptedTracksForSharedStatus = (context).workspace.acceptedTracks

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

void validateSparsePlan(const TraversalWorkspaceView& context, int iteration, const SurfaceGraphView& layout)
{
  BIND_TRAVERSAL_CONTEXT(context);
  const auto fail = [iteration]() { throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch}; };
  const auto& topology = layout;
  if (layout.surfaces == nullptr || layout.nSurfaces == 0 ||
      (topology.nTransitions != 0 && (topology.transitions == nullptr || topology.cellsByFirstTransitionOffsets == nullptr)) ||
      (topology.nCells != 0 && (topology.cells == nullptr || topology.cellsByFirstTransition == nullptr))) {
    fail();
  }

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

void prepareTransitionScatteringAndBending(
  TraversalWorkspaceView& context,
  int iteration,
  const DiskReferenceCoordinateView& referenceCoordinateView,
  gsl::span<const TransitionId> transitionIds)
{
  BIND_TRAVERSAL_CONTEXT(context);
  const auto& trkParam = mTrkParams[iteration];
  const auto& topology = mTraversalGraph;

  const int activeSurfaceCount = static_cast<int>(mScratch->getNOwnedSurfaces());
  std::array<float, MaxLayoutSurfaces> msAngles{};
  for (int iLayer{0}; iLayer < activeSurfaceCount; ++iLayer) {
    const auto surface = mBinding->getOrderedSurfaces()[iLayer];
    if (topology.getSurface(surface).kind == SurfaceKind::Cylinder) {
      msAngles[iLayer] = cylinderLayerMultipleScatteringAngle(
        CylinderLayerScatteringInputs{mAttachHitConfig.layerMaterial[iLayer].xOverX0}, trkParam.TrackletMinPt);
    } else {
      msAngles[iLayer] = diskLayerMultipleScatteringAngle(
        DiskLayerScatteringInputs{mAttachHitConfig.layerMaterial[iLayer].xOverX0, trkParam.LayerRadii[iLayer],
                                  referenceCoordinateView.perLayerReferenceZ[iLayer]},
        trkParam.TrackletMinPt);
    }
  }

  auto& transitionMSAngles = mScratch->getTransitionMSAngles();
  auto& transitionPhiCuts = mScratch->getTransitionPhiCuts();
  const float oneOverR{0.001f * 0.3f * std::abs(mBz) / trkParam.TrackletMinPt};
  for (const auto transitionId : transitionIds) {
    const auto transitionSlot = mBinding->getScratchTransitionSlot(transitionId);
    if (!transitionSlot) {
      throw TraversalException{iteration, TraversalFailureReason::TraversalBindingMismatch};
    }
    const auto& transition = topology.getTransition(transitionId);
    const auto fromPosition = mBinding->getOwnedSurfaceIndex(transition.from);
    const auto toPosition = mBinding->getOwnedSurfaceIndex(transition.to);
    if (!fromPosition || !toPosition || *fromPosition >= static_cast<std::size_t>(activeSurfaceCount) ||
        *toPosition >= static_cast<std::size_t>(activeSurfaceCount)) {
      throw TraversalException{iteration, TraversalFailureReason::TraversalBindingMismatch};
    }
    const int fromLayer = static_cast<int>(*fromPosition);
    const int toLayer = static_cast<int>(*toPosition);
    const float r1 = std::min(trkParam.LayerRadii[fromLayer], trkParam.LayerRadii[toLayer]);
    const float r2 = std::max(trkParam.LayerRadii[fromLayer], trkParam.LayerRadii[toLayer]);
    const float transitionOneOverR = clampTransitionCurvature(oneOverR, r2);
    const float res1 = o2::gpu::CAMath::Hypot(trkParam.PVres, mScratch->getPositionResolution(fromLayer));
    const float res2 = o2::gpu::CAMath::Hypot(trkParam.PVres, mScratch->getPositionResolution(toLayer));
    const auto prep = ::o2::itsmft::tracking::prepareTransitionScatteringAndBending(
      gsl::span<const float>(msAngles.data(), static_cast<std::size_t>(activeSurfaceCount)), fromLayer, toLayer, r1, r2, transitionOneOverR, res1, res2);
    transitionMSAngles[*transitionSlot] = prep.msAngle;
    transitionPhiCuts[*transitionSlot] = prep.phiCut;
  }
}
} // namespace

void Tracker::initializeTraversalWorkspace(TraversalWorkspaceView& context) const
{
  BIND_TRAVERSAL_CONTEXT(context);
  const int iteration = context.iteration;
  // Invalidate before preflight so no failed preparation can leave a reusable
  // traversal view. The workspace itself is reset only after preflight.
  mTraversalCacheValid = false;

  if (iteration < 0 || static_cast<size_t>(iteration) >= mTrkParams.size()) {
    throw TraversalException{iteration, TraversalFailureReason::IterationOutOfRange};
  }
  const auto layout = mTraversalGraph;
  if (mBinding == nullptr) {
    throw TraversalException{iteration, TraversalFailureReason::MissingLayout};
  }

  const auto boundTransitions = mBinding->getGlobalTransitions();
  if (boundTransitions.empty() || boundTransitions.front().value() >= layout.nTransitions) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }

  for (const auto surface : mBinding->getOrderedSurfaces()) {
    if (materialCorrectionModeSupport(layout.getSurface(surface).kind, mTrkParams[iteration].CorrType) ==
        MaterialCorrectionModeSupport::Unsupported) {
      throw TraversalException{iteration, TraversalFailureReason::UnsupportedMaterialCorrectionMode};
    }
  }

  // 2.5. Resolve and validate material through this graph's orderedSurfaces
  // before touching TimeFrame state. LayerxX0 mismatches and invalid mappings
  // use LegacyMaterialMismatch. Keep stagedLayerMaterial local until all
  // fallible checks succeed, so failures leave the cache reset.
  //
  const auto orderedSurfaces = mBinding->getOrderedSurfaces();
  const int activeSurfaceCount = static_cast<int>(mScratch->getNOwnedSurfaces());
  // The application configuration must match the adopted plan count.
  if (activeSurfaceCount <= 0 || activeSurfaceCount > MaxLayoutSurfaces ||
      orderedSurfaces.size() != static_cast<std::size_t>(activeSurfaceCount) ||
      mTrkParams[iteration].NLayers != activeSurfaceCount) {
    throw TraversalException{iteration, TraversalFailureReason::LegacyMaterialMismatch};
  }
  const auto activeCount = static_cast<std::size_t>(activeSurfaceCount);
  std::array<NominalSurfaceMaterial, MaxLayoutSurfaces> stagedLayerMaterial{};
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

  // 2.6. Resolve and validate normalized measurements from orderedSurfaces
  // and the loaded TimeFrame before touching tracking state. Keep staged spans
  // local until all fallible checks succeed; failures leave the cache empty.
  std::array<gsl::span<const SurfaceMeasurement>, MaxLayoutSurfaces> stagedLayerMeasurements{};
  std::array<gsl::span<const GlobalMeasurement>, MaxLayoutSurfaces> stagedLayerGlobalMeasurements{};
  for (int surfacePosition = 0; surfacePosition < activeSurfaceCount; ++surfacePosition) {
    const auto surfaceId = orderedSurfaces[surfacePosition];
    const auto measurements = mFrame->getSurfaceMeasurements(surfaceId);
    const auto globals = mFrame->getGlobalMeasurements(surfaceId);
    const auto& clusters = mScratch->getUnsortedClusters()[surfacePosition];
    const auto& hits = mScratch->getTrackingFrameInfoOnLayer(surfacePosition);
    if (measurements.size() != globals.size() || globals.size() != clusters.size() || hits.size() != clusters.size() ||
        measurements.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
    }
    for (size_t i = 0; i < measurements.size(); ++i) {
      const auto& global = globals[i];
      if (global.surface != surfaceId || !global.cluster.isValid() ||
          global.cluster.index > static_cast<uint32_t>(std::numeric_limits<int>::max()) ||
          clusters[i].clusterId != static_cast<int>(i)) {
        throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
      }
      const int externalIndex = mScratch->getClusterExternalIndex(surfacePosition, static_cast<int>(i));
      if (externalIndex < 0 || static_cast<uint32_t>(externalIndex) != global.cluster.index) {
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
        if (globals[index].sourceROF != rof) {
          throw TraversalException{iteration, TraversalFailureReason::NormalizedMeasurementMismatch};
        }
      }
    }
    stagedLayerMeasurements[surfacePosition] = measurements;
    stagedLayerGlobalMeasurements[surfacePosition] = globals;
  }

  std::array<IndexTableUtilsCore, 2> indexTableConfigByKind{};
  std::array<IndexTableUtilsCore, MaxLayoutSurfaces> stagedIndexTableConfigs{};
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

  if (!mTrkParams[iteration].PassFlags[IterationStep::FirstPass]) {
    for (int position = 0; position < activeSurfaceCount; ++position) {
      if (!indexTableConfigurationsMatch(stagedIndexTableConfigs[position], mScratch->getIndexTableUtils(position), activeSurfaceCount)) {
        throw TraversalException{iteration, TraversalFailureReason::IndexTableConfigurationMismatch};
      }
    }
  }

  // All inputs needed for scratch sizing and traversal-derived state now
  // passed validation. Reset the selected iteration workspace as one unit;
  // it also clears any prior accepted candidates for this iteration.
  context.workspace.reset(mMemoryPool.get());

  const auto transitionIds = mBinding->getGlobalTransitions();
  const auto cellIds = mBinding->getGlobalCells();
  const gsl::span<const IndexTableUtilsCore> stagedIndexTableConfigView{stagedIndexTableConfigs.data(), activeCount};
  const gsl::span<const gsl::span<const GlobalMeasurement>> stagedLayerGlobalMeasurementView{
    stagedLayerGlobalMeasurements.data(), activeCount};
  // initialise() builds allocator-backed locator containers needed by the
  // remaining sorted-ROF validation. If that validation fails, Tracker::run()
  // resets the whole TimeFrame before returning or propagating the error.
  mScratch->initialise(*mFrame, mTrkParams[iteration], activeSurfaceCount, iteration,
                       stagedIndexTableConfigView, layout, transitionIds, cellIds,
                       orderedSurfaces, stagedLayerGlobalMeasurementView);

  // Sorted clusters are a locator cache. Validate every enabled ROF that can
  // participate in a configured transition, including LUT-reuse paths.
  // Keep spans local until validation and kind setup complete.
  std::array<bool, MaxLayoutSurfaces> candidateReachableLayers{};
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
    // Orchestration-only users may omit the mask; without it no ROF is reachable.
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
        const auto& measurement = stagedLayerGlobalMeasurements[layer][clusterId];
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

  validateSparsePlan(context, iteration, layout);

  // Bind the bounded preflight material; rebind to frame-owned workspace storage at commit.
  auto attachHitConfig = bindAttachHitConfig(
    gsl::span<const NominalSurfaceMaterial>(stagedLayerMaterial.data(), activeCount), mTrkParams[iteration]);
  if (!attachHitConfig.isValid(activeSurfaceCount)) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
  }
  if (mTrkParams[iteration].LayerRadii.size() < static_cast<size_t>(activeSurfaceCount)) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
  }
  DiskReferenceCoordinateView referenceCoordinateView{};
  std::array<float, MaxLayoutSurfaces> stagedDiskLayerReferenceZ{};
  const auto hasDiskSurface = std::ranges::any_of(orderedSurfaces, [&layout](SurfaceId id) {
    return layout.getSurface(id).kind == SurfaceKind::Disk;
  });
  if (hasDiskSurface) {
    const auto legacyDiskReference = bindLegacyMFTReferenceCoordinates();
    std::size_t diskLayer = 0;
    for (int position = 0; position < activeSurfaceCount; ++position) {
      if (layout.getSurface(orderedSurfaces[position]).kind != SurfaceKind::Disk) {
        continue;
      }
      if (diskLayer >= legacyDiskReference.perLayerReferenceZ.size()) {
        throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
      }
      stagedDiskLayerReferenceZ[position] = legacyDiskReference.perLayerReferenceZ[diskLayer++];
    }
    referenceCoordinateView = {gsl::span<const float>{stagedDiskLayerReferenceZ.data(), activeCount}};
    if (!referenceCoordinateView.isValid(activeSurfaceCount)) {
      throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
    }
  }
  mKernelParameters = bindTrackingKernelParameters(mTrkParams[iteration]);
  if (!mKernelParameters.isValid()) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
  }

  mDiskLayerReferenceZStorage.assign(referenceCoordinateView.perLayerReferenceZ.begin(),
                                     referenceCoordinateView.perLayerReferenceZ.end());
  mDiskLayerReferenceZ = mDiskLayerReferenceZStorage;
  // Commit material with the other traversal caches, then rebind attachHitConfig
  // from the local staging span to the persistent member.
  mLayerMaterial.assign(stagedLayerMaterial.begin(), stagedLayerMaterial.begin() + activeSurfaceCount);
  attachHitConfig.layerMaterial = gsl::span<const NominalSurfaceMaterial>(mLayerMaterial.data(), mLayerMaterial.size());
  mAttachHitConfig = attachHitConfig;
  // Commit normalized measurements with the other traversal caches.
  mLayerMeasurements.assign(stagedLayerMeasurements.begin(), stagedLayerMeasurements.begin() + activeSurfaceCount);
  mLayerGlobalMeasurements.assign(stagedLayerGlobalMeasurements.begin(), stagedLayerGlobalMeasurements.begin() + activeSurfaceCount);
  mTraversalCacheValid = true;

  // All fallible validation is complete. The remaining per-layer and
  // per-transition scattering/bending preparation is non-throwing.
  prepareTransitionScatteringAndBending(context, iteration, referenceCoordinateView,
                                        mBinding->getGlobalTransitions());
}

#undef BIND_TRAVERSAL_CONTEXT

TrackerInitializationResult Tracker::initialize(TimeFrame& frame, const TrackerInitialization& configuration)
{
  TrackerInitializationResult result;
  if (configuration.iterations.empty()) {
    result.error = TrackerInitializationError::EmptyConfiguration;
    return result;
  }
  if (configuration.catalog.surfaces == nullptr || configuration.catalog.nSurfaces == 0) {
    result.error = TrackerInitializationError::MissingCatalog;
    return result;
  }
  if (!configuration.memoryPool) {
    result.error = TrackerInitializationError::MissingMemoryPool;
    return result;
  }

  std::vector<SurfaceGraph> graphs;
  std::vector<TrackingParameters> parameters;
  std::vector<std::unique_ptr<SurfacePlanBinding>> bindings;
  std::vector<TrackingWorkspaceCapacity> capacities;
  graphs.reserve(configuration.iterations.size());
  parameters.reserve(configuration.iterations.size());
  bindings.reserve(configuration.iterations.size());
  capacities.reserve(configuration.iterations.size());

  for (std::size_t iteration = 0; iteration < configuration.iterations.size(); ++iteration) {
    const auto& input = configuration.iterations[iteration];
    const auto graphResult = SurfaceGraphBuilder{configuration.catalog, input.graph}.build();
    if (!graphResult.ok()) {
      result.error = TrackerInitializationError::GraphBuildFailed;
      result.failedIteration = iteration;
      result.graphError = graphResult.error;
      return result;
    }

    if (input.graph.orderedSurfaces.empty()) {
      result.error = TrackerInitializationError::BindingBuildFailed;
      result.failedIteration = iteration;
      return result;
    }
    SurfaceMask ownedSurfaces;
    for (const auto surface : input.graph.orderedSurfaces) {
      ownedSurfaces.set(surface);
    }
    auto bindingResult = SurfacePlanBinding::build(graphResult.graph->getView(), ownedSurfaces,
                                                   input.graph.orderedSurfaces);
    if (!bindingResult.ok()) {
      result.error = TrackerInitializationError::BindingBuildFailed;
      result.failedIteration = iteration;
      result.bindingError = bindingResult.error;
      return result;
    }
    if (input.parameters.NLayers != 0 && input.parameters.NLayers != input.graph.orderedSurfaces.size()) {
      result.error = TrackerInitializationError::CapacityMismatch;
      result.failedIteration = iteration;
      return result;
    }
    capacities.push_back(TrackingWorkspaceCapacity{
      input.graph.orderedSurfaces.size(), bindingResult.binding->getGlobalTransitions().size(),
      bindingResult.binding->getGlobalCells().size()});
    graphs.push_back(*graphResult.graph);
    parameters.push_back(input.parameters);
    bindings.push_back(std::move(bindingResult.binding));
  }

  if (!frame.commitConfiguration(std::move(graphs), std::move(parameters), std::move(bindings),
                                 std::move(capacities), configuration.memoryPool)) {
    result.error = TrackerInitializationError::CapacityMismatch;
    return result;
  }
  return result;
}

TrackingResult Tracker::run(TimeFrame& frame, TrackerTraits& traits)
{
  if (mRefitFunction == nullptr || !frame.isConfigured() ||
      frame.getNIterations() == 0 || frame.getBinding(0) == nullptr) {
    throw TraversalException{-1, TraversalFailureReason::MissingLayout};
  }
  auto& scratch = frame.getWorkspace();
  const auto& trkParams = frame.getTrackingParameters();
  const auto& memoryPool = frame.getMemoryPool();

  float total{0.f};
  std::vector<std::size_t> acceptedTrackCounts;
  acceptedTrackCounts.reserve(trkParams.size());
  try {
    for (int iteration = 0; iteration < static_cast<int>(trkParams.size()); ++iteration) {
      // Apply a tighter event-local limit when configured; this also lets
      // workflows and tests inject a resource failure after loading.
      if (trkParams[iteration].MaxMemory != std::numeric_limits<size_t>::max() &&
          memoryPool->getMaxMemory() > trkParams[iteration].MaxMemory) {
        memoryPool->setMaxMemory(trkParams[iteration].MaxMemory);
      }
      if (trkParams[iteration].PassFlags[IterationStep::UseUPCMask]) {
        scratch.useUPCMask();
      }

      const auto* binding = frame.getBinding(iteration);
      if (binding == nullptr) {
        throw TraversalException{iteration, TraversalFailureReason::MissingLayout};
      }
      auto& workspace = scratch.getTraversalWorkspace(static_cast<std::size_t>(iteration));
      TraversalWorkspaceView view{iteration, frame, scratch, frame.getGraph(iteration).getView(), *binding,
                                  trkParams, frame.getBz(), workspace};
      initializeTraversalWorkspace(view);
      traits.runTraversal(view, mRefitFunction);
      acceptedTrackCounts.push_back(workspace.acceptedTracks.size());
    }
  } catch (const TraversalException& err) {
    // Structural/configuration failures are not per-TF data failures, so
    // DropTFUponFailure does not apply. Reset before propagating.
    LOGP(error, "CA tracker hit a structural traversal failure: {}", err.what());
    frame.resetEvent();
    throw;
  } catch (const BoundedMemoryResource::MemoryLimitExceeded& err) {
    // Recoverable per-TF resource failure: the bounded pool budget was
    // exceeded for this TimeFrame.
    LOGP(error, "CA tracker exceeded memory limit: {}", err.what());
    frame.resetEvent();
    if (trkParams[0].DropTFUponFailure) {
      return TrackingResult{TrackingOutcome::RecoverableDropped, 0.f};
    }
    throw;
  } catch (const std::bad_alloc& err) {
    // Some CA scratch containers use the plain heap instead of the bounded
    // pool, so memory pressure can surface as bad_alloc. Handle it likewise.
    LOGP(error, "CA tracker allocation failed: {}", err.what());
    frame.resetEvent();
    if (trkParams[0].DropTFUponFailure) {
      return TrackingResult{TrackingOutcome::RecoverableDropped, 0.f};
    }
    throw;
  } catch (const std::exception& err) {
    // Unclassified exceptions are treated as structural and always propagate;
    // recoverability is not inferred from std::exception alone.
    LOGP(error, "CA tracker failed with an unclassified exception; treating as structural: {}", err.what());
    frame.resetEvent();
    throw;
  }

  return TrackingResult{TrackingOutcome::Success, total, std::move(acceptedTrackCounts)};
}

} // namespace o2::itsmft::tracking

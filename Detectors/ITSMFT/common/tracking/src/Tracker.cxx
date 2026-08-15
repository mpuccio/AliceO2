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
#include <numeric>
#include <queue>
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

SurfaceLayoutDefinition makeSurfaceLayoutDefinition(const SurfaceGraphDefinition& graph)
{
  SurfaceLayoutDefinition layout;
  layout.orderedSurfaces = graph.orderedSurfaces;
  layout.maxHoles = graph.maxHoles;
  layout.holeSurfaces = graph.holeSurfaces;
  layout.seedingSurfaces = graph.seedingSurfaces;
  layout.componentOffsets = {0};
  for (uint16_t cut = 0; cut + 1 < graph.orderedSurfaces.size(); ++cut) {
    const auto crossesCut = std::any_of(graph.basePairs.begin(), graph.basePairs.end(), [cut](const auto pair) {
      return pair.fromIndex <= cut && pair.toIndex > cut;
    });
    if (!crossesCut) {
      layout.componentOffsets.push_back(static_cast<uint16_t>(cut + 1));
    }
  }
  return layout;
}

SurfaceMask disabledSurfaceMask(const SurfaceLayout& layout, LayerMask inactiveLayers)
{
  SurfaceMask disabled;
  const auto ordered = layout.getOrderedSurfaces();
  for (uint16_t position = 0; position < ordered.size(); ++position) {
    if (inactiveLayers.has(static_cast<int>(position))) {
      disabled.set(ordered[position]);
    }
  }
  return disabled;
}

} // namespace

void buildTraversalPlanImpl(TraversalWorkspace& workspace, const SurfaceGraphView& graph, LayerMask inactiveLayers, int iteration)
{
  const auto fail = [iteration]() {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  };
  if (graph.surfaces == nullptr || graph.nSurfaces == 0 || graph.orderedSurfaces == nullptr || graph.nOrderedSurfaces == 0 ||
      (graph.nEdges != 0 && graph.edges == nullptr) || (graph.nCells != 0 && graph.cells == nullptr)) {
    fail();
  }

  workspace.orderedSurfaces.clear();
  workspace.orderedSurfaces.reserve(graph.nOrderedSurfaces);
  for (uint32_t position = 0; position < graph.nOrderedSurfaces; ++position) {
    workspace.orderedSurfaces.push_back(graph.orderedSurfaces[position]);
  }
  if (workspace.orderedSurfaces.empty()) {
    fail();
  }
  workspace.topologyCatalog = graph.getSurfaceCatalogView();
  workspace.topology.orderedSurfaces = workspace.orderedSurfaces;
  workspace.topology.surfacePositionById.assign(MaxLayoutSurfaces, -1);
  for (uint16_t position = 0; position < workspace.orderedSurfaces.size(); ++position) {
    workspace.topology.surfacePositionById[workspace.orderedSurfaces[position].value()] = static_cast<int16_t>(position);
  }
  workspace.surfaceSlotById.assign(graph.nSurfaces, -1);
  for (uint16_t position = 0; position < workspace.orderedSurfaces.size(); ++position) {
    const auto surface = workspace.orderedSurfaces[position];
    if (!surface.isValid() || surface.value() >= graph.nSurfaces || workspace.activeSurfaces.has(surface) ||
        workspace.surfaceSlotById[surface.value()] >= 0) {
      fail();
    }
    if (!inactiveLayers.has(static_cast<int>(position))) {
      workspace.activeSurfaces.set(surface);
    }
    workspace.surfaceSlotById[surface.value()] = static_cast<int16_t>(position);
  }
  if (workspace.activeSurfaces.empty()) {
    fail();
  }
  workspace.topology.activeSurfaces = workspace.activeSurfaces;
  workspace.topology.activeSurfaceList.clear();
  for (const auto surface : workspace.orderedSurfaces) {
    if (workspace.activeSurfaces.has(surface)) {
      workspace.topology.activeSurfaceList.push_back(surface);
    }
  }
  workspace.topology.seedingSurfaces = graph.seedingSurfaces;

  std::vector<uint32_t> indegree(graph.nSurfaces, 0);
  std::vector<uint32_t> rank(graph.nSurfaces, 0);
  for (uint32_t rawId = 0; rawId < graph.nEdges; ++rawId) {
    const auto& edge = graph.getEdge(EdgeId{static_cast<uint16_t>(rawId)});
    if (!edge.from.isValid() || !edge.to.isValid() || edge.from.value() >= graph.nSurfaces || edge.to.value() >= graph.nSurfaces) {
      fail();
    }
    ++indegree[edge.to.value()];
  }
  const auto laterSurface = [](SurfaceId lhs, SurfaceId rhs) { return rhs < lhs; };
  std::priority_queue<SurfaceId, std::vector<SurfaceId>, decltype(laterSurface)> ready{laterSurface};
  for (uint16_t surface = 0; surface < graph.nSurfaces; ++surface) {
    if (indegree[surface] == 0) {
      ready.push(SurfaceId{surface});
    }
  }
  uint32_t visited = 0;
  while (!ready.empty()) {
    const auto surface = ready.top();
    ready.pop();
    ++visited;
    for (uint32_t rawId = 0; rawId < graph.nEdges; ++rawId) {
      const auto& edge = graph.getEdge(EdgeId{static_cast<uint16_t>(rawId)});
      if (edge.from != surface) {
        continue;
      }
      rank[edge.to.value()] = std::max(rank[edge.to.value()], rank[surface.value()] + 1);
      if (--indegree[edge.to.value()] == 0) {
        ready.push(edge.to);
      }
    }
  }
  if (visited != graph.nSurfaces) {
    fail();
  }
  std::vector<uint16_t> componentParent(graph.nSurfaces);
  std::iota(componentParent.begin(), componentParent.end(), uint16_t{0});
  const auto componentRoot = [&componentParent](uint16_t surface) {
    while (componentParent[surface] != surface) {
      componentParent[surface] = componentParent[componentParent[surface]];
      surface = componentParent[surface];
    }
    return surface;
  };
  for (uint32_t rawId = 0; rawId < graph.nEdges; ++rawId) {
    const auto& edge = graph.getEdge(EdgeId{static_cast<uint16_t>(rawId)});
    const auto fromRoot = componentRoot(edge.from.value());
    const auto toRoot = componentRoot(edge.to.value());
    if (fromRoot != toRoot) {
      componentParent[toRoot] = fromRoot;
    }
  }
  std::vector<uint32_t> componentOrder(graph.nSurfaces, std::numeric_limits<uint32_t>::max());
  for (uint32_t position = 0; position < workspace.orderedSurfaces.size(); ++position) {
    const auto root = componentRoot(workspace.orderedSurfaces[position].value());
    componentOrder[root] = std::min(componentOrder[root], position);
  }

  workspace.edgeSlotById.clear();
  workspace.cellSlotById.clear();
  std::vector<int16_t> selectedEdgeByGraphId(graph.nEdges, -1);
  for (uint32_t rawId = 0; rawId < graph.nEdges; ++rawId) {
    const auto id = EdgeId{static_cast<uint16_t>(rawId)};
    const auto& edge = graph.getEdge(id);
    const bool fromActive = workspace.activeSurfaces.has(edge.from);
    const bool toActive = workspace.activeSurfaces.has(edge.to);
    if (!fromActive || !toActive) {
      continue;
    }
    const auto selectedId = EdgeId{static_cast<uint16_t>(workspace.topology.edges.size())};
    selectedEdgeByGraphId[id.value()] = static_cast<int16_t>(selectedId.value());
    workspace.topology.edges.push_back(edge);
    workspace.edgeSlotById.push_back(static_cast<int16_t>(workspace.edges.size()));
    workspace.edges.push_back(selectedId);
  }
  if (workspace.edges.empty()) {
    fail();
  }
  for (uint32_t rawId = 0; rawId < graph.nCells; ++rawId) {
    const auto id = CellTopologyId{static_cast<uint16_t>(rawId)};
    const auto& cell = graph.getCell(id);
    if (!cell.firstEdge.isValid() || !cell.secondEdge.isValid() || cell.firstEdge.value() >= graph.nEdges || cell.secondEdge.value() >= graph.nEdges) {
      fail();
    }
    const auto firstSelected = selectedEdgeByGraphId[cell.firstEdge.value()];
    const auto secondSelected = selectedEdgeByGraphId[cell.secondEdge.value()];
    const bool firstActive = firstSelected >= 0;
    const bool secondActive = secondSelected >= 0;
    if (!firstActive || !secondActive || !cell.hitSurfaces.isSubsetOf(workspace.activeSurfaces)) {
      continue;
    }
    const auto selectedId = CellTopologyId{static_cast<uint16_t>(workspace.topology.paths.size())};
    workspace.topology.paths.push_back(CellPath{EdgeId{static_cast<uint16_t>(firstSelected)}, EdgeId{static_cast<uint16_t>(secondSelected)}});
    workspace.cellSlotById.push_back(static_cast<int16_t>(workspace.cells.size()));
    workspace.cells.push_back(selectedId);
  }
  const auto scheduleOrder = [&](CellTopologyId lhs, CellTopologyId rhs) {
    const auto lhsTarget = workspace.topology.edges[workspace.topology.paths[lhs.value()].second.value()].to;
    const auto rhsTarget = workspace.topology.edges[workspace.topology.paths[rhs.value()].second.value()].to;
    const auto lhsComponent = componentOrder[componentRoot(lhsTarget.value())];
    const auto rhsComponent = componentOrder[componentRoot(rhsTarget.value())];
    if (lhsComponent != rhsComponent) {
      return lhsComponent < rhsComponent;
    }
    return rank[lhsTarget.value()] != rank[rhsTarget.value()] ? rank[lhsTarget.value()] < rank[rhsTarget.value()] : lhs < rhs;
  };
  workspace.scheduledCells = workspace.cells;
  std::sort(workspace.scheduledCells.begin(), workspace.scheduledCells.end(), scheduleOrder);
  for (const auto id : workspace.cells) {
    if (workspace.topology.seedingSurfaces.has(workspace.topology.edges[workspace.topology.paths[id.value()].second.value()].to)) {
      workspace.roadStartCells.push_back(id);
    }
  }
  std::sort(workspace.roadStartCells.begin(), workspace.roadStartCells.end(), scheduleOrder);
  workspace.roadStartComponentOffsets.push_back(0);
  uint32_t previousComponent = std::numeric_limits<uint32_t>::max();
  for (uint32_t offset = 0; offset < workspace.roadStartCells.size(); ++offset) {
    const auto target = workspace.topology.edges[workspace.topology.paths[workspace.roadStartCells[offset].value()].second.value()].to;
    const auto component = componentOrder[componentRoot(target.value())];
    if (component != previousComponent && offset != 0) {
      workspace.roadStartComponentOffsets.push_back(offset);
    }
    previousComponent = component;
  }
  workspace.roadStartComponentOffsets.push_back(static_cast<uint32_t>(workspace.roadStartCells.size()));

  workspace.topology.pathsByFirstEdgeOffsets.assign(workspace.topology.edges.size() + 1, 0);
  for (const auto& path : workspace.topology.paths) {
    ++workspace.topology.pathsByFirstEdgeOffsets[path.first.value() + 1];
  }
  for (size_t offset = 1; offset < workspace.topology.pathsByFirstEdgeOffsets.size(); ++offset) {
    workspace.topology.pathsByFirstEdgeOffsets[offset] += workspace.topology.pathsByFirstEdgeOffsets[offset - 1];
  }
  workspace.topology.pathsByFirstEdge.resize(workspace.topology.paths.size());
  auto cursor = workspace.topology.pathsByFirstEdgeOffsets;
  for (uint32_t path = 0; path < workspace.topology.paths.size(); ++path) {
    workspace.topology.pathsByFirstEdge[cursor[workspace.topology.paths[path].first.value()]++] = CellTopologyId{static_cast<uint16_t>(path)};
  }
  workspace.topology.scheduledPaths = workspace.scheduledCells;
  workspace.topology.roadStartPaths = workspace.roadStartCells;
  workspace.topology.roadStartComponentOffsets = workspace.roadStartComponentOffsets;
}

namespace
{
void validateSparsePlan(const TraversalWorkspaceView& context, int iteration, const TraversalTopologyView& layout)
{
  const auto& workspace = context.workspace;
  const auto parameters = context.parameters;
  const auto fail = [iteration]() { throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch}; };
  const auto& topology = layout;
  if (layout.catalog.surfaces == nullptr || layout.catalog.nSurfaces == 0 ||
      (topology.nEdges != 0 && (topology.edges == nullptr || topology.pathsByFirstEdgeOffsets == nullptr)) ||
      (topology.nPaths != 0 && (topology.paths == nullptr || topology.pathsByFirstEdge == nullptr))) {
    fail();
  }

  const auto& edges = workspace.edges;
  const auto& cells = workspace.cells;
  if (edges.empty() || edges.size() > topology.nEdges || cells.size() > topology.nPaths) {
    fail();
  }
  for (const auto id : edges) {
    if (!id.isValid() || id.value() >= topology.nEdges || !workspace.getEdgeSlot(id)) {
      fail();
    }
    const auto& edge = topology.getEdge(id);
    if (!workspace.getSurfaceSlot(edge.from) || !workspace.getSurfaceSlot(edge.to)) {
      fail();
    }
  }
  for (const auto id : cells) {
    if (!id.isValid() || id.value() >= topology.nPaths || !workspace.getCellSlot(id)) {
      fail();
    }
    const auto& path = topology.getPath(id);
    const auto& firstEdge = topology.getEdge(path.first);
    const auto& secondEdge = topology.getEdge(path.second);
    if (!workspace.getEdgeSlot(path.first) || !workspace.getEdgeSlot(path.second) ||
        !workspace.activeSurfaces.has(firstEdge.from) || !workspace.activeSurfaces.has(firstEdge.to) ||
        !workspace.activeSurfaces.has(secondEdge.to)) {
      fail();
    }
  }
  for (const auto id : workspace.scheduledCells) {
    if (!workspace.getCellSlot(id)) {
      fail();
    }
  }
  for (const auto id : workspace.roadStartCells) {
    if (!workspace.getCellSlot(id)) {
      fail();
    }
  }
  (void)parameters[iteration];
}

void prepareTraversalEdgeTolerances(
  TraversalWorkspaceView& context,
  int iteration,
  const DiskReferenceCoordinateView& referenceCoordinateView,
  gsl::span<const EdgeId> edgeIds)
{
  auto& scratch = context.scratch;
  const auto parameters = context.parameters;
  const auto& graph = context.topology;
  const auto& attachHitConfig = context.workspace.attachHitConfig;
  const auto& trkParam = parameters[iteration];
  const auto& topology = graph;

  const int activeSurfaceCount = static_cast<int>(context.workspace.orderedSurfaces.size());
  std::array<float, MaxLayoutSurfaces> msAngles{};
  for (int iLayer{0}; iLayer < activeSurfaceCount; ++iLayer) {
    const auto surface = context.workspace.orderedSurfaces[iLayer];
    if (topology.getSurface(surface).kind == SurfaceKind::Cylinder) {
      msAngles[iLayer] = cylinderLayerMultipleScatteringAngle(
        CylinderLayerScatteringInputs{attachHitConfig.layerMaterial[iLayer].xOverX0}, trkParam.TrackletMinPt);
    } else {
      msAngles[iLayer] = diskLayerMultipleScatteringAngle(
        DiskLayerScatteringInputs{attachHitConfig.layerMaterial[iLayer].xOverX0, trkParam.LayerRadii[iLayer],
                                  referenceCoordinateView.perLayerReferenceZ[iLayer]},
        trkParam.TrackletMinPt);
    }
  }

  auto& edgeMSAngles = scratch.getEdgeMSAngles();
  auto& edgePhiCuts = scratch.getEdgePhiCuts();
  const float oneOverR{0.001f * 0.3f * std::abs(context.bz) / trkParam.TrackletMinPt};
  for (const auto edgeId : edgeIds) {
    const auto edgeSlot = context.workspace.getEdgeSlot(edgeId);
    if (!edgeSlot) {
      throw TraversalException{iteration, TraversalFailureReason::TraversalBindingMismatch};
    }
    const auto& edge = topology.getEdge(edgeId);
    const auto fromPosition = context.workspace.getSurfaceSlot(edge.from);
    const auto toPosition = context.workspace.getSurfaceSlot(edge.to);
    if (!fromPosition || !toPosition || *fromPosition >= static_cast<std::size_t>(activeSurfaceCount) ||
        *toPosition >= static_cast<std::size_t>(activeSurfaceCount)) {
      throw TraversalException{iteration, TraversalFailureReason::TraversalBindingMismatch};
    }
    const int fromLayer = static_cast<int>(*fromPosition);
    const int toLayer = static_cast<int>(*toPosition);
    const float r1 = std::min(trkParam.LayerRadii[fromLayer], trkParam.LayerRadii[toLayer]);
    const float r2 = std::max(trkParam.LayerRadii[fromLayer], trkParam.LayerRadii[toLayer]);
    const float edgeOneOverR = clampEdgeCurvature(oneOverR, r2);
    const float res1 = o2::gpu::CAMath::Hypot(trkParam.PVres, scratch.getPositionResolution(fromLayer));
    const float res2 = o2::gpu::CAMath::Hypot(trkParam.PVres, scratch.getPositionResolution(toLayer));
    const auto prep = ::o2::itsmft::tracking::prepareEdgeScatteringAndBending(
      gsl::span<const float>(msAngles.data(), static_cast<std::size_t>(activeSurfaceCount)), fromLayer, toLayer, r1, r2, edgeOneOverR, res1, res2);
    edgeMSAngles[*edgeSlot] = prep.msAngle;
    edgePhiCuts[*edgeSlot] = prep.phiCut;
  }
}
} // namespace

void Tracker::buildTraversalPlan(TraversalWorkspace& workspace, const SurfaceLayout& layout, SurfaceMask disabledSurfaces, int iteration) const
{
  TraversalWorkspace staged;
  const auto derived = deriveTraversalTopology(layout, disabledSurfaces);
  if (!derived.ok()) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  staged.topology = *derived.topology;
  staged.topologyCatalog = layout.getSurfaceCatalog();
  staged.activeSurfaces = staged.topology.activeSurfaces;
  staged.orderedSurfaces = staged.topology.orderedSurfaces;
  staged.surfaceSlotById = staged.topology.surfacePositionById;
  staged.edgeSlotById.reserve(staged.topology.edges.size());
  staged.edges.reserve(staged.topology.edges.size());
  for (uint16_t edge = 0; edge < staged.topology.edges.size(); ++edge) {
    staged.edgeSlotById.push_back(static_cast<int16_t>(edge));
    staged.edges.push_back(EdgeId{edge});
  }
  staged.cellSlotById.reserve(staged.topology.paths.size());
  staged.cells.reserve(staged.topology.paths.size());
  for (uint16_t path = 0; path < staged.topology.paths.size(); ++path) {
    staged.cellSlotById.push_back(static_cast<int16_t>(path));
    staged.cells.push_back(CellTopologyId{path});
  }
  staged.scheduledCells = staged.topology.scheduledPaths;
  staged.roadStartCells = staged.topology.roadStartPaths;
  staged.roadStartComponentOffsets = staged.topology.roadStartComponentOffsets;

  workspace.activeSurfaces = staged.activeSurfaces;
  workspace.topologyCatalog = staged.topologyCatalog;
  workspace.topology = std::move(staged.topology);
  workspace.orderedSurfaces = std::move(staged.orderedSurfaces);
  workspace.surfaceSlotById = std::move(staged.surfaceSlotById);
  workspace.edgeSlotById = std::move(staged.edgeSlotById);
  workspace.cellSlotById = std::move(staged.cellSlotById);
  workspace.edges = std::move(staged.edges);
  workspace.cells = std::move(staged.cells);
  workspace.roadStartCells = std::move(staged.roadStartCells);
  workspace.roadStartComponentOffsets = std::move(staged.roadStartComponentOffsets);
  workspace.scheduledCells = std::move(staged.scheduledCells);
}

void Tracker::initializeTraversalWorkspace(TraversalWorkspaceView& context) const
{
  auto* mScratch = &context.scratch;
  auto* mFrame = &context.frame;
  const auto& mMemoryPool = mScratch->getMemoryPool();
  const auto mTrkParams = context.parameters;
  const auto& surfaceLayout = mFrame->getLayout(static_cast<std::size_t>(context.iteration));
  auto& mTraversalCacheValid = context.workspace.valid;
  auto& mKernelParameters = context.workspace.kernelParameters;
  auto& mDiskLayerReferenceZ = context.workspace.diskLayerReferenceZView;
  auto& mDiskLayerReferenceZStorage = context.workspace.diskLayerReferenceZ;
  auto& mAttachHitConfig = context.workspace.attachHitConfig;
  auto& mLayerMaterial = context.workspace.layerMaterial;
  auto& mLayerMeasurements = context.workspace.layerMeasurements;
  auto& mLayerGlobalMeasurements = context.workspace.layerGlobalMeasurements;
  const int iteration = context.iteration;
  // Invalidate before preflight so no failed preparation can leave a reusable
  // traversal view. The workspace itself is reset only after preflight.
  mTraversalCacheValid = false;

  if (iteration < 0 || static_cast<size_t>(iteration) >= mTrkParams.size()) {
    throw TraversalException{iteration, TraversalFailureReason::IterationOutOfRange};
  }
  context.workspace.reset(mMemoryPool.get());
  buildTraversalPlan(context.workspace, surfaceLayout, disabledSurfaceMask(surfaceLayout, mTrkParams[iteration].InactiveLayerMask), iteration);
  context.topology = context.workspace.getTopologyView();
  const auto& layout = context.topology;
  const auto& boundEdges = context.workspace.edges;

  for (const auto surface : context.workspace.orderedSurfaces) {
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
  const auto orderedSurfaces = context.workspace.orderedSurfaces;
  const int activeSurfaceCount = static_cast<int>(context.workspace.orderedSurfaces.size());
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
    if (!surfaceId.isValid() || surfaceId.value() >= layout.catalog.nSurfaces) {
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

  std::array<IndexTableUtilsCore, MaxLayoutSurfaces> stagedIndexTableConfigs{};
  std::array<SurfaceChartRange, MaxLayoutSurfaces> chartRanges{};
  for (int position = 0; position < activeSurfaceCount; ++position) {
    chartRanges[position] = layout.getSurface(orderedSurfaces[position]).chartRange;
  }
  const gsl::span<const SurfaceChartRange> chartRangeView{chartRanges.data(), activeCount};
  for (int position = 0; position < activeSurfaceCount; ++position) {
    if (bindIndexTableConfiguration(stagedIndexTableConfigs[position], mTrkParams[iteration], activeSurfaceCount,
                                    layout.getSurface(orderedSurfaces[position]).kind, chartRangeView) != IndexTableConfigError::None) {
      throw TraversalException{iteration, TraversalFailureReason::InvalidIndexTableConfiguration};
    }
  }

  if (!mTrkParams[iteration].PassFlags[IterationStep::FirstPass]) {
    for (int position = 0; position < activeSurfaceCount; ++position) {
      if (!indexTableConfigurationsMatch(stagedIndexTableConfigs[position], mScratch->getIndexTableUtils(position), activeSurfaceCount)) {
        throw TraversalException{iteration, TraversalFailureReason::IndexTableConfigurationMismatch};
      }
    }
  }

  const auto& edgeIds = context.workspace.edges;
  const auto& cellIds = context.workspace.cells;
  const gsl::span<const IndexTableUtilsCore> stagedIndexTableConfigView{stagedIndexTableConfigs.data(), activeCount};
  const gsl::span<const gsl::span<const GlobalMeasurement>> stagedLayerGlobalMeasurementView{
    stagedLayerGlobalMeasurements.data(), activeCount};
  // initialise() builds allocator-backed locator containers needed by the
  // remaining sorted-ROF validation. If that validation fails, Tracker::run()
  // resets the whole TimeFrame before returning or propagating the error.
  mScratch->initialise(*mFrame, mTrkParams[iteration], activeSurfaceCount, iteration,
                       stagedIndexTableConfigView, context.topology, edgeIds, cellIds,
                       orderedSurfaces, stagedLayerGlobalMeasurementView);

  // Sorted clusters are a locator cache. Validate every enabled ROF that can
  // participate in a configured edge, including LUT-reuse paths.
  // Keep spans local until validation and kind setup complete.
  std::array<bool, MaxLayoutSurfaces> candidateReachableLayers{};
  for (const auto edgeId : edgeIds) {
    const auto& edge = layout.getEdge(edgeId);
    const auto fromSlot = context.workspace.getSurfaceSlot(edge.from);
    const auto toSlot = context.workspace.getSurfaceSlot(edge.to);
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

  validateSparsePlan(context, iteration, context.topology);

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
  // per-edge scattering/bending preparation is non-throwing.
  prepareTraversalEdgeTolerances(context, iteration, referenceCoordinateView,
                                 context.workspace.edges);
}

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
  std::vector<SurfaceLayout> layouts;
  std::vector<TrackingParameters> parameters;
  std::vector<TrackingWorkspaceCapacity> capacities;
  graphs.reserve(configuration.iterations.size());
  layouts.reserve(configuration.iterations.size());
  parameters.reserve(configuration.iterations.size());
  capacities.reserve(configuration.iterations.size());

  for (std::size_t iteration = 0; iteration < configuration.iterations.size(); ++iteration) {
    const auto& input = configuration.iterations[iteration];
    SurfaceLayout layout{gsl::span<const SurfaceDescriptor>{configuration.catalog.surfaces, configuration.catalog.nSurfaces},
                         makeSurfaceLayoutDefinition(input.graph)};
    if (!layout.valid()) {
      result.error = TrackerInitializationError::TraversalPlanBuildFailed;
      result.failedIteration = iteration;
      return result;
    }
    const auto graphResult = SurfaceGraphBuilder{configuration.catalog, input.graph}.build();
    if (!graphResult.ok()) {
      result.error = TrackerInitializationError::GraphBuildFailed;
      result.failedIteration = iteration;
      result.graphError = graphResult.error;
      return result;
    }

    if (input.graph.orderedSurfaces.empty()) {
      result.error = TrackerInitializationError::TraversalPlanBuildFailed;
      result.failedIteration = iteration;
      return result;
    }
    if (input.parameters.NLayers != 0 && input.parameters.NLayers != input.graph.orderedSurfaces.size()) {
      result.error = TrackerInitializationError::CapacityMismatch;
      result.failedIteration = iteration;
      return result;
    }
    capacities.push_back(TrackingWorkspaceCapacity{
      input.graph.orderedSurfaces.size(), graphResult.graph->getEdges().size(),
      graphResult.graph->getCells().size()});
    layouts.push_back(std::move(layout));
    graphs.push_back(*graphResult.graph);
    parameters.push_back(input.parameters);
  }

  if (!frame.commitConfiguration(std::move(layouts), std::move(graphs), std::move(parameters),
                                 std::move(capacities), configuration.memoryPool)) {
    result.error = TrackerInitializationError::CapacityMismatch;
    return result;
  }
  return result;
}

TrackingResult Tracker::run(TimeFrame& frame, TrackerTraits& traits)
{
  if (mRefitFunction == nullptr || !frame.isConfigured() || frame.getNIterations() == 0) {
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

      auto& workspace = scratch.getTraversalWorkspace(static_cast<std::size_t>(iteration));
      TraversalWorkspaceView view{iteration, frame, scratch, frame.getLayout(iteration).getSurfaceCatalog(),
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

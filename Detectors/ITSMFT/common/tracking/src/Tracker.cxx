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
#include <utility>

#include "Framework/Logger.h"
#include "ITStracking/BoundedAllocator.h"

namespace o2::itsmft::tracking
{

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
  std::vector<std::array<TrackingParameters, 2>> parametersByKind;
  std::vector<std::unique_ptr<SurfacePlanBinding>> bindings;
  std::vector<TrackingWorkspaceCapacity> capacities;
  graphs.reserve(configuration.iterations.size());
  parameters.reserve(configuration.iterations.size());
  parametersByKind.reserve(configuration.iterations.size());
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
    std::array<TrackingParameters, 2> resolvedByKind{input.parameters, input.parameters};
    for (std::size_t kind = 0; kind < resolvedByKind.size(); ++kind) {
      if (input.parametersByKind[kind]) {
        resolvedByKind[kind] = *input.parametersByKind[kind];
      }
      if (resolvedByKind[kind].NLayers != 0 && resolvedByKind[kind].NLayers != input.graph.orderedSurfaces.size()) {
        result.error = TrackerInitializationError::CapacityMismatch;
        result.failedIteration = iteration;
        return result;
      }
    }
    capacities.push_back(TrackingWorkspaceCapacity{
      input.graph.orderedSurfaces.size(), bindingResult.binding->getGlobalTransitions().size(),
      bindingResult.binding->getGlobalCells().size()});
    graphs.push_back(*graphResult.graph);
    parameters.push_back(input.parameters);
    parametersByKind.push_back(std::move(resolvedByKind));
    bindings.push_back(std::move(bindingResult.binding));
  }

  if (!frame.commitConfiguration(std::move(graphs), std::move(parameters), std::move(parametersByKind), std::move(bindings),
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
  const auto& trkParamsByKind = frame.getTrackingParametersByKind();
  const auto& memoryPool = frame.getMemoryPool();
  traits.adoptFrame(&frame);
  traits.adoptScratch(&scratch);
  traits.setBz(frame.getBz());
  traits.updateTrackingParameters(trkParams, trkParamsByKind);

  int maxNvertices{-1};
  if (trkParams[0].PerPrimaryVertexProcessing) {
    maxNvertices = scratch.getMaxVerticesPerROF();
  }

  float total{0.f};
  std::vector<std::size_t> acceptedTrackCounts;
  acceptedTrackCounts.reserve(trkParams.size());
  try {
    for (int iteration = 0; iteration < static_cast<int>(trkParams.size()); ++iteration) {
      // Keep a deliberately tightened event-local bound. This is also the
      // only way a workflow/test can inject a resource failure after loading;
      // configuration still supplies the normal upper bound.
      if (trkParams[iteration].MaxMemory != std::numeric_limits<size_t>::max() &&
          memoryPool->getMaxMemory() > trkParams[iteration].MaxMemory) {
        memoryPool->setMaxMemory(trkParams[iteration].MaxMemory);
      }
      if (trkParams[iteration].PassFlags[IterationStep::UseUPCMask]) {
        scratch.useUPCMask();
      }

      int iVertex = std::min(maxNvertices, 0);
      traits.adoptSurfacePlanBinding(frame.getBinding(iteration));
      traits.initialiseTimeFrame(iteration, frame.getGraphs());
      do {
        traits.computeLayerTracklets(iteration, iVertex);
        traits.computeLayerCells(iteration);
        traits.findCellsNeighbours(iteration);
        traits.findRoads(iteration, mRefitFunction);
      } while (++iVertex < maxNvertices);
      acceptedTrackCounts.push_back(traits.acceptedTracksForSharedStatus().size());
    }
  } catch (const TraversalException& err) {
    // Structural/configuration failure (bad layout, stale layout, surface-kind or
    // index mismatch): never a per-TF data problem, so DropTFUponFailure
    // never applies. Always reset before propagating -- see class-level
    // comment: never rely on "the process is going down anyway".
    LOGP(error, "CA tracker hit a structural traversal failure: {}", err.what());
    frame.resetEvent();
    throw;
  } catch (const BoundedMemoryResource::MemoryLimitExceeded& err) {
    // Recoverable, per-TF resource failure: the bounded pool's configured
    // budget was exceeded for this TimeFrame's data volume.
    LOGP(error, "CA tracker exceeded memory limit: {}", err.what());
    frame.resetEvent();
    if (trkParams[0].DropTFUponFailure) {
      return TrackingResult{TrackingOutcome::RecoverableDropped, 0.f};
    }
    throw;
  } catch (const std::bad_alloc& err) {
    // Also recoverable/per-TF: several CA scratch containers on the hot path
    // (e.g. TrackerTraits::findCellsNeighbours' cellsNeighboursByTarget and
    // its per-thread TBB storage) allocate from the plain heap rather than
    // the bounded pool, so genuine memory pressure surfaces here as a plain
    // bad_alloc rather than MemoryLimitExceeded. Handled identically.
    LOGP(error, "CA tracker allocation failed: {}", err.what());
    frame.resetEvent();
    if (trkParams[0].DropTFUponFailure) {
      return TrackingResult{TrackingOutcome::RecoverableDropped, 0.f};
    }
    throw;
  } catch (const std::exception& err) {
    // Unclassified: not a recognized recoverable-resource failure, so it is
    // treated as structural and always propagates, regardless of
    // DropTFUponFailure. A future explicitly typed
    // RecoverableTimeFrameException may extend the recoverable set; until
    // then, recoverability is never inferred from std::exception alone.
    LOGP(error, "CA tracker failed with an unclassified exception; treating as structural: {}", err.what());
    frame.resetEvent();
    throw;
  }

  return TrackingResult{TrackingOutcome::Success, total, std::move(acceptedTrackCounts)};
}

} // namespace o2::itsmft::tracking

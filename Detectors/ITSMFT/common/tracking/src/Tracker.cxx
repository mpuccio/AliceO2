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

#include "Framework/Logger.h"
#include "ITStracking/BoundedAllocator.h"

namespace o2::itsmft::tracking
{

Tracker::Tracker(TrackerTraits* traits) : mTraits(traits)
{
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
  std::vector<std::vector<TrackingParameters>> parameters;
  std::vector<TimeFrame::BindingSet> bindings;
  std::vector<std::vector<TrackingWorkspaceCapacity>> capacities;
  graphs.reserve(configuration.iterations.size());
  parameters.reserve(configuration.iterations.size());
  bindings.reserve(configuration.iterations.size());
  capacities.reserve(configuration.iterations.size());

  for (std::size_t iteration = 0; iteration < configuration.iterations.size(); ++iteration) {
    const auto& input = configuration.iterations[iteration];
    SurfaceGraphBuilder builder{configuration.catalog};
    for (const auto& subgraph : input.graphSubgraphs) {
      builder.addSubgraph(subgraph);
    }
    const auto graphResult = builder.build();
    if (!graphResult.ok()) {
      result.error = TrackerInitializationError::GraphBuildFailed;
      result.failedIteration = iteration;
      result.graphError = graphResult.error;
      return result;
    }

    TimeFrame::BindingSet iterationBindings;
    std::vector<TrackingWorkspaceCapacity> iterationCapacities;
    if (input.bindings.empty() || input.bindings.size() != input.graphSubgraphs.size() ||
        input.parameters.size() != input.bindings.size()) {
      result.error = TrackerInitializationError::BindingCountMismatch;
      result.failedIteration = iteration;
      return result;
    }
    for (const auto& declaration : input.bindings) {
      if (std::any_of(iterationBindings.begin(), iterationBindings.end(), [&](const auto& binding) {
            return binding && binding->getSource() == declaration.source;
          })) {
        result.error = TrackerInitializationError::DuplicateSource;
        result.failedIteration = iteration;
        return result;
      }
      auto bindingResult = SurfacePlanBinding::build(graphResult.graph->getView(), declaration);
      if (!bindingResult.ok()) {
        result.error = TrackerInitializationError::BindingBuildFailed;
        result.failedIteration = iteration;
        result.bindingError = bindingResult.error;
        return result;
      }
      if (input.parameters[iterationBindings.size()].NLayers != 0 &&
          input.parameters[iterationBindings.size()].NLayers != declaration.orderedSurfaces.size()) {
        result.error = TrackerInitializationError::CapacityMismatch;
        result.failedIteration = iteration;
        return result;
      }
      iterationCapacities.push_back(TrackingWorkspaceCapacity{
        declaration.orderedSurfaces.size(), bindingResult.binding->getGlobalTransitions().size(),
        bindingResult.binding->getGlobalCells().size()});
      iterationBindings.push_back(std::move(bindingResult.binding));
    }
    graphs.push_back(*graphResult.graph);
    parameters.push_back(input.parameters);
    bindings.push_back(std::move(iterationBindings));
    capacities.push_back(std::move(iterationCapacities));
  }

  for (const auto& iterationBindings : bindings) {
    if (std::none_of(iterationBindings.begin(), iterationBindings.end(), [&](const auto& binding) {
          return binding && binding->getSource() == mSource;
        })) {
      result.error = TrackerInitializationError::BindingCountMismatch;
      return result;
    }
  }

  if (!frame.commitConfiguration(std::move(graphs), std::move(parameters), std::move(bindings),
                                 std::move(capacities), configuration.memoryPool)) {
    result.error = TrackerInitializationError::CapacityMismatch;
    return result;
  }
  adoptFrame(frame);
  const auto& sourceParameters = frame.getTrackingParameters(mSource);
  mTraits->updateTrackingParameters(sourceParameters);
  return result;
}

void Tracker::adoptFrame(TimeFrame& frame)
{
  mFrame = &frame;
  mTraits->adoptFrame(&frame);
  if (frame.isConfigured()) {
    mTraits->adoptScratch(&frame.getWorkspace(mSource));
  }
}

TrackingResult Tracker::clustersToTracks(TrackingOperationAdapter& operationAdapter)
{
  if (mFrame == nullptr || !mFrame->isConfigured() ||
      mFrame->getNIterations() == 0 || mFrame->getBinding(0, mSource) == nullptr) {
    throw TraversalException{-1, TraversalFailureReason::MissingLayout};
  }
  auto& scratch = mFrame->getWorkspace(mSource);
  const auto& trkParams = mFrame->getTrackingParameters(mSource);
  const auto& memoryPool = mFrame->getMemoryPool();
  mTraits->updateTrackingParameters(trkParams);

  int maxNvertices{-1};
  if (trkParams[0].PerPrimaryVertexProcessing) {
    maxNvertices = scratch.getROFVertexLookupView().getMaxVerticesPerROF();
  }

  float total{0.f};
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
      initialiseTimeFrame(iteration);
      do {
        computeTracklets(iteration, iVertex);
        computeCells(iteration);
        findCellsNeighbours(iteration);
        findRoads(iteration, operationAdapter);
      } while (++iVertex < maxNvertices);
    }
  } catch (const TraversalException& err) {
    // Structural/configuration failure (bad layout, stale layout, policy or
    // index mismatch): never a per-TF data problem, so DropTFUponFailure
    // never applies. Always reset before propagating -- see class-level
    // comment: never rely on "the process is going down anyway".
    LOGP(error, "CA tracker hit a structural traversal failure: {}", err.what());
    operationAdapter.resetAdapterState();
    mFrame->resetEvent();
    throw;
  } catch (const BoundedMemoryResource::MemoryLimitExceeded& err) {
    // Recoverable, per-TF resource failure: the bounded pool's configured
    // budget was exceeded for this TimeFrame's data volume.
    LOGP(error, "CA tracker exceeded memory limit: {}", err.what());
    operationAdapter.resetAdapterState();
    mFrame->resetEvent();
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
    operationAdapter.resetAdapterState();
    mFrame->resetEvent();
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
    operationAdapter.resetAdapterState();
    mFrame->resetEvent();
    throw;
  }

  return TrackingResult{TrackingOutcome::Success, total};
}

} // namespace o2::itsmft::tracking

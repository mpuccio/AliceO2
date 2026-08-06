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

void Tracker::adoptScratch(SurfaceTrackingScratch& scratch)
{
  mScratch = &scratch;
  mTraits->adoptScratch(&scratch);
}

void Tracker::adoptFrame(TimeFrame& frame)
{
  mFrame = &frame;
  mTraits->adoptFrame(&frame);
}

TrackingResult Tracker::clustersToTracks(TrackingOperationAdapter& operationAdapter)
{
  mTraits->updateTrackingParameters(mTrkParams);

  int maxNvertices{-1};
  if (mTrkParams[0].PerPrimaryVertexProcessing) {
    maxNvertices = mScratch->getROFVertexLookupView().getMaxVerticesPerROF();
  }

  float total{0.f};
  try {
    for (int iteration = 0; iteration < static_cast<int>(mTrkParams.size()); ++iteration) {
      mMemoryPool->setMaxMemory(mTrkParams[iteration].MaxMemory);
      if (mTrkParams[iteration].PassFlags[IterationStep::UseUPCMask]) {
        mScratch->useUPCMask();
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
    resetTimeFrameEvent(*mFrame, *mScratch);
    throw;
  } catch (const BoundedMemoryResource::MemoryLimitExceeded& err) {
    // Recoverable, per-TF resource failure: the bounded pool's configured
    // budget was exceeded for this TimeFrame's data volume.
    LOGP(error, "CA tracker exceeded memory limit: {}", err.what());
    operationAdapter.resetAdapterState();
    resetTimeFrameEvent(*mFrame, *mScratch);
    if (mTrkParams[0].DropTFUponFailure) {
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
    resetTimeFrameEvent(*mFrame, *mScratch);
    if (mTrkParams[0].DropTFUponFailure) {
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
    resetTimeFrameEvent(*mFrame, *mScratch);
    throw;
  }

  return TrackingResult{TrackingOutcome::Success, total};
}

} // namespace o2::itsmft::tracking

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_TEST_TRAVERSALTESTSUPPORT_H_
#define ALICEO2_ITSMFT_TRACKING_TEST_TRAVERSALTESTSUPPORT_H_

#include <stdexcept>

#include "ITSMFTTracking/Tracker.h"

namespace o2::itsmft::tracking
{

// Test-only access to the Tracker-owned initialization transaction and
// explicit backend stages. It never persists the resulting view.
struct TrackerTestAccess {
  static void useDiskXYReference(Tracker& tracker, bool enabled) noexcept
  {
    tracker.mUseDiskXYReferenceForTesting = enabled;
  }

  static TraversalWorkspaceView prepare(Tracker& tracker, TimeFrame& frame, int iteration)
  {
    if (iteration < 0 || static_cast<std::size_t>(iteration) >= frame.getTrackingParameters().size()) {
      throw std::out_of_range{"test traversal iteration"};
    }
    const auto* binding = frame.getBinding(static_cast<std::size_t>(iteration));
    if (binding == nullptr) {
      throw std::runtime_error{"test traversal requires configured binding"};
    }
    auto& scratch = frame.getWorkspace();
    TraversalWorkspaceView view{iteration,
                                frame,
                                scratch,
                                frame.getGraph(static_cast<std::size_t>(iteration)).getView(),
                                *binding,
                                frame.getTrackingParameters(),
                                frame.getBz(),
                                scratch.getTraversalWorkspace(static_cast<std::size_t>(iteration))};
    tracker.initializeTraversalWorkspace(view);
    return view;
  }

  static void computeTracklets(TrackerTraits& traits, TraversalWorkspaceView& view, int vertex)
  {
    traits.computeLayerTracklets(view, view.iteration, vertex);
  }

  static void computeCells(TrackerTraits& traits, TraversalWorkspaceView& view)
  {
    traits.computeLayerCells(view, view.iteration);
  }

  static void findNeighbours(TrackerTraits& traits, TraversalWorkspaceView& view)
  {
    traits.findCellsNeighbours(view, view.iteration);
  }

  static void findRoads(TrackerTraits& traits, TraversalWorkspaceView& view, SeedRefitFunction refit)
  {
    traits.findRoads(view, view.iteration, refit);
  }
};

} // namespace o2::itsmft::tracking

#endif

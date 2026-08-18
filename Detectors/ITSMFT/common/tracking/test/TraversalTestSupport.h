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
  static IterationContext prepare(Tracker& tracker, TimeFrame& frame, int iteration)
  {
    const auto* configuration = iteration < 0 ? nullptr : tracker.getIterationConfiguration(static_cast<std::size_t>(iteration));
    if (configuration == nullptr || !tracker.isConfiguredFor(frame)) {
      throw std::out_of_range{"test traversal iteration"};
    }
    auto& scratch = frame.getScratch();
    auto& iterationScratch = scratch.getIteration(static_cast<std::size_t>(iteration));
    auto layerGlobalMeasurements = tracker.prepareTimeFrame(frame);
    IterationContext view{iteration,
                          frame,
                          scratch,
                          configuration->getTopologyView(frame.getLayout().getSurfaceCatalog()),
                          *configuration,
                          tracker.mDetectorConfiguration,
                          std::move(layerGlobalMeasurements),
                          frame.getBz(),
                          iterationScratch};
    tracker.initializeIterationScratch(view);
    return view;
  }

  static void computeTracklets(TrackerTraits& traits, IterationContext& view, int vertex)
  {
    traits.computeLayerTracklets(view, view.iteration, vertex);
  }

  static void computeCells(TrackerTraits& traits, IterationContext& view)
  {
    traits.computeLayerCells(view, view.iteration);
  }

  static void findNeighbours(TrackerTraits& traits, IterationContext& view)
  {
    traits.findCellsNeighbours(view, view.iteration);
  }

  static void findRoads(TrackerTraits& traits, IterationContext& view, SeedRefitFunction refit)
  {
    traits.findRoads(view, view.iteration, refit);
  }
};

} // namespace o2::itsmft::tracking

#endif

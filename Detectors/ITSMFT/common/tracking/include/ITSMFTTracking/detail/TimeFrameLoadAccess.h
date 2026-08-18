// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.

#ifndef ALICEO2_ITSMFT_TRACKING_DETAIL_TIMEFRAMELOADACCESS_H_
#define ALICEO2_ITSMFT_TRACKING_DETAIL_TIMEFRAMELOADACCESS_H_

#include "ITSMFTTracking/TimeFrame.h"

namespace o2::itsmft::tracking::detail
{

// Narrow internal access for input adapters. TimeFrame remains independent
// of detector input representations while adapters fill its owned storage.
struct TimeFrameLoadAccess {
  static void prepareMeasurements(TimeFrame& frame, std::size_t nSurfaces)
  {
    frame.mLayerGlobalMeasurements.resize(nSurfaces);
    frame.mLayerSurfaceMeasurements.resize(nSurfaces);
    frame.mLayerClusterLabels.resize(nSurfaces);
    frame.mLayerUsedClusters.resize(nSurfaces);
  }

  static auto& globalMeasurements(TimeFrame& frame) { return frame.mLayerGlobalMeasurements; }
  static auto& surfaceMeasurements(TimeFrame& frame) { return frame.mLayerSurfaceMeasurements; }
  static auto& clusterLabels(TimeFrame& frame) { return frame.mLayerClusterLabels; }

  static void finishMeasurements(TimeFrame& frame, bool hasMCInformation)
  {
    for (std::size_t surface = 0; surface < frame.mLayerSurfaceMeasurements.size(); ++surface) {
      frame.mLayerUsedClusters[surface].assign(frame.mLayerSurfaceMeasurements[surface].size(), uint8_t{0});
    }
    frame.mHasMCInformation = hasMCInformation;
  }

  static void setMeasurements(
    TimeFrame& frame, std::vector<std::vector<GlobalMeasurement>>&& globals,
    std::vector<std::vector<SurfaceMeasurement>>&& measurements,
    std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>&& labels,
    bool hasMCInformation)
  {
    frame.mLayerGlobalMeasurements = std::move(globals);
    frame.mLayerSurfaceMeasurements = std::move(measurements);
    frame.mLayerClusterLabels = std::move(labels);
    frame.mLayerUsedClusters.resize(frame.mLayerSurfaceMeasurements.size());
    finishMeasurements(frame, hasMCInformation);
  }

  static void setNavigation(TimeFrame& frame, std::vector<std::vector<int>>&& rofBoundaries,
                            RuntimeROFViews defaultViews,
                            std::vector<RuntimeROFViews>&& viewsBySurface,
                            std::vector<uint16_t>&& localLayerBySurface)
  {
    frame.mROFramesClusters = std::move(rofBoundaries);
    frame.mROFViews = defaultViews;
    frame.mROFViewsBySurface = std::move(viewsBySurface);
    frame.mROFLocalLayerBySurface = std::move(localLayerBySurface);
    frame.mUseUPC = false;
  }
};

} // namespace o2::itsmft::tracking::detail

#endif

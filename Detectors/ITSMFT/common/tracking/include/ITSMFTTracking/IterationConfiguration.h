// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.

#ifndef ALICEO2_ITSMFT_TRACKING_ITERATIONCONFIGURATION_H_
#define ALICEO2_ITSMFT_TRACKING_ITERATIONCONFIGURATION_H_

#include <optional>
#include <vector>

#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/TraversalTopology.h"
#include "ITSMFTTracking/detail/TrackingKernelParameters.h"

namespace o2::itsmft::tracking
{

// Tracker-owned data derived once from the invariant detector layout.
struct DetectorConfiguration {
  std::vector<NominalSurfaceMaterial> layerMaterial;
  std::vector<float> layerRadii;
  std::vector<IndexTableUtilsCore> indexTableConfigs;
  std::vector<float> positionResolutions;
  std::vector<uint32_t> addTimeError;
  std::vector<float> layerResolution;
  std::vector<float> systError2Row;
  std::vector<float> systError2Col;
};

// Tracker-owned, immutable instructions for one tracking iteration.
struct IterationConfiguration {
  IterationParameters parameters;
  TraversalTopology topology;
  std::vector<EdgeId> edges;
  std::vector<CellPathId> cells;
  TrackingKernelParameters kernelParameters{};

  std::optional<uint16_t> getSurfaceSlot(LayerId id) const noexcept
  {
    if (!id.isValid() || id.value() >= topology.surfacePositionById.size() || topology.surfacePositionById[id.value()] < 0) {
      return std::nullopt;
    }
    return static_cast<uint16_t>(topology.surfacePositionById[id.value()]);
  }

  std::optional<uint16_t> getEdgeSlot(EdgeId id) const noexcept
  {
    return id.isValid() && id.value() < edges.size() ? std::optional<uint16_t>{id.value()} : std::nullopt;
  }

  std::optional<uint16_t> getCellSlot(CellPathId id) const noexcept
  {
    return id.isValid() && id.value() < cells.size() ? std::optional<uint16_t>{id.value()} : std::nullopt;
  }

  TraversalTopologyView getTopologyView(SurfaceCatalogView catalog) const noexcept
  {
    return topology.getView(catalog);
  }
};

} // namespace o2::itsmft::tracking

#endif

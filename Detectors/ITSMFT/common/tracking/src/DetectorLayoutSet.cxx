// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file DetectorLayoutSet.cxx
/// \brief One-shot DetectorLayoutSet construction from a borrowed static catalog
///

#include "ITSMFTTracking/DetectorLayoutSet.h"

namespace o2::itsmft::tracking
{

DetectorLayoutSetBuildResult buildDetectorLayoutSet(SurfaceCatalogView catalog,
                                                    gsl::span<const SurfaceId> orderedSurfaces,
                                                    TransitionPolicyTag policyTag,
                                                    gsl::span<const o2::itsmft::TrackingParameters> trackingParameters)
{
  DetectorLayoutConfigurationKey key;
  key.orderedSurfaces.assign(orderedSurfaces.begin(), orderedSurfaces.end());
  key.policyTag = policyTag;
  key.iterations.reserve(trackingParameters.size());
  for (const auto& parameters : trackingParameters) {
    if (parameters.NLayers < 0 || static_cast<size_t>(parameters.NLayers) > orderedSurfaces.size()) {
      return {.error = DetectorLayoutSetBuildError::InvalidActiveCount,
              .failedIteration = key.iterations.size()};
    }
    key.iterations.push_back(DetectorLayoutIterationConfiguration{
      static_cast<uint32_t>(parameters.NLayers), parameters.MaxHoles,
      parameters.HoleLayerMask, parameters.StartLayerMask});
  }

  std::vector<DetectorLayout> staging;
  staging.reserve(key.iterations.size());
  for (size_t iteration = 0; iteration < key.iterations.size(); ++iteration) {
    const auto& configuration = key.iterations[iteration];
    std::vector<SurfaceId> activeSurfaces(orderedSurfaces.begin(), orderedSurfaces.begin() + configuration.activeCount);
    DetectorLayoutSubgraph subgraph;
    subgraph.orderedSurfaces = std::move(activeSurfaces);
    subgraph.maxHoles = configuration.maxHoles;
    subgraph.holeSurfaces = positionalSurfaceMask(configuration.holeLayerMask, orderedSurfaces, configuration.activeCount);
    subgraph.seedingSurfaces = positionalSurfaceMask(configuration.startLayerMask, orderedSurfaces, configuration.activeCount);
    subgraph.policyTag = policyTag;

    DetectorLayoutBuilder builder{catalog};
    auto buildResult = builder.addSubgraph(std::move(subgraph)).build();
    if (!buildResult.ok()) {
      return {.error = DetectorLayoutSetBuildError::LayoutBuilderFailure,
              .failedIteration = iteration,
              .layoutBuildError = buildResult.error,
              .topologyError = buildResult.topologyError,
              .layoutError = buildResult.layoutError};
    }
    staging.push_back(std::move(*buildResult.layout));
  }

  DetectorLayoutSetBuildResult result;
  result.rebuilt = true;
  result.layout.emplace(std::move(key), catalog, std::move(staging));
  return result;
}

} // namespace o2::itsmft::tracking

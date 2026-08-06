// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"

#include <algorithm>
#include <memory>
#include <vector>

#include "ITSMFTTracking/TimeFrame.h"

namespace o2::itsmft::tracking
{

LoadSourcesResult MultiSourceTimeFrameLoader::load(TimeFrame& frame, gsl::span<const ClusterSourceInput> sources,
                                                   SurfaceCatalogView catalog, const o2::InteractionRecord& origin)
{
  if (!frame.isConfigured()) {
    return {MultiSourceLoadError::FrameNotConfigured};
  }
  // Stage normalized ownership first, generically, over every source at
  // once. The existing source-level loader remains the authoritative decoder
  // and timing validator.
  MultiSourceFrame normalized;
  const auto normalizedResult = loadSources(normalized, catalog, sources, origin);
  if (!normalizedResult.ok()) {
    return normalizedResult;
  }

  if (sources.size() != frame.getNConfiguredSources()) {
    return {MultiSourceLoadError::NonDenseSourceIds};
  }

  std::vector<std::unique_ptr<SurfaceTrackingScratch>> stagedWorkspaces;
  stagedWorkspaces.reserve(sources.size());
  for (const auto& source : sources) {
    const auto* binding = frame.getBinding(0, source.id);
    if (binding == nullptr || binding->getSource() != source.id ||
        source.layerToSurface.size() != binding->getOrderedSurfaces().size() ||
        !std::equal(source.layerToSurface.begin(), source.layerToSurface.end(), binding->getOrderedSurfaces().begin())) {
      return {MultiSourceLoadError::InvalidLayerMapping, source.id};
    }

    auto staged = std::make_unique<SurfaceTrackingScratch>();
    auto& live = frame.getWorkspace(source.id);
    if (live.hasFrameworkAllocator()) {
      staged->setFrameworkAllocator(live.getFrameworkAllocator());
    }
    staged->setMemoryPool(frame.getMemoryPool());
    staged->adoptPlan(live.getNOwnedSurfaces(), 0, 0);
    staged->setROFViews(source.rofViews);

    TimeFrame disposable;
    ClusterSourceInput one = source;
    one.id = ClusterSourceId{0};
    auto result = staged->loadNormalizedSource(disposable, *one.decoder, origin, one.timing,
                                               one.clusters, one.patterns, one.rofs,
                                               one.dictionary, one.labels, one.detector,
                                               one.layerToSurface, catalog, one.applySysErrors);
    if (!result.ok()) {
      result.source = source.id;
      return result;
    }
    stagedWorkspaces.push_back(std::move(staged));
  }

  std::vector<ClusterSourceId> sourceOrder;
  sourceOrder.reserve(sources.size());
  for (const auto& source : sources) {
    sourceOrder.push_back(source.id);
  }
  if (!frame.commitLoadedEvent(std::move(normalized), sourceOrder, std::move(stagedWorkspaces))) {
    return {MultiSourceLoadError::OtherMalformedInput};
  }
  return {};
}

} // namespace o2::itsmft::tracking

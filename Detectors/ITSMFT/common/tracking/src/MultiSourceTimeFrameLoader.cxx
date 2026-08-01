// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"

#include <array>

namespace o2::itsmft::tracking
{
template <int NLayers>
void MultiSourceTimeFrameLoader::prepareStagingScratch(LegacyTrackerScratch<NLayers>& staged, const LegacyTrackerScratch<NLayers>& live)
{
  // Keep allocator identity, not merely allocator configuration: pmr swaps
  // in the final commit are defined only for equal resources.
  staged.mExternalAllocator = live.mExternalAllocator;
  staged.mExtMemoryPool = live.mExtMemoryPool;
  staged.setMemoryPool(live.mMemoryPool);
}

template <int NLayers>
void MultiSourceTimeFrameLoader::commitNormalizedBackfill(LegacyTrackerScratch<NLayers>& live, LegacyTrackerScratch<NLayers>& staged) noexcept
{
  // This is the same load-owned subset committed by loadNormalizedSource().
  // All swaps are same-allocator and therefore allocation-free/no-throw.
  for (int layer = 0; layer < NLayers; ++layer) {
    live.mUnsortedClusters[layer].swap(staged.mUnsortedClusters[layer]);
    live.mTrackingFrameInfo[layer].swap(staged.mTrackingFrameInfo[layer]);
    live.mClusterExternalIndices[layer].swap(staged.mClusterExternalIndices[layer]);
    live.mClusterSize[layer].swap(staged.mClusterSize[layer]);
    live.mROFramesClusters[layer].swap(staged.mROFramesClusters[layer]);
    live.mClusterLabels[layer] = staged.mClusterLabels[layer];
  }
  for (int i = 0; i < 2; ++i) {
    live.mNTrackletsPerCluster[i].swap(staged.mNTrackletsPerCluster[i]);
    live.mNTrackletsPerClusterSum[i].swap(staged.mNTrackletsPerClusterSum[i]);
  }
}

template <int NLayers>
LoadSourcesResult stageScratch(LegacyTrackerScratch<NLayers>& staged,
                               const ClusterSourceInput& source,
                               SurfaceCatalogView catalog,
                               const o2::InteractionRecord& origin)
{
  // The old bridge's one-source API correctly owns the legacy conversion.
  // It is invoked only on a disposable TimeFrame; the shared owner's
  // normalized frame was already staged once, above, from both sources.
  TimeFrame disposable;
  ClusterSourceInput one = source;
  one.id = ClusterSourceId{0};
  return staged.loadNormalizedSource(disposable, *one.decoder, origin, one.timing,
                                     one.clusters, one.patterns, one.rofs,
                                     one.dictionary, one.labels, one.detector,
                                     one.layerToSurface, catalog, one.applySysErrors);
}
LoadSourcesResult MultiSourceTimeFrameLoader::loadITSAndMFT(TimeFrame& frame,
                                                            LegacyTrackerScratchITS& itsScratch,
                                                            LegacyTrackerScratchMFT& mftScratch,
                                                            const ClusterSourceInput& itsSource,
                                                            const ClusterSourceInput& mftSource,
                                                            SurfaceCatalogView catalog,
                                                            const o2::InteractionRecord& origin)
{
  // This owner is deliberately ITS+MFT-shaped, rather than a generic
  // heterogeneous tracker entry point. Keeping the source positions fixed
  // makes source metadata and every ClusterRef stable and makes it impossible
  // to accidentally backfill a detector into the other scratch.
  if (itsSource.id != ClusterSourceId{0} || itsSource.detector != o2::detectors::DetID::ITS) {
    return {MultiSourceLoadError::UnsupportedDetector, itsSource.id};
  }
  if (mftSource.id != ClusterSourceId{1} || mftSource.detector != o2::detectors::DetID::MFT) {
    return {MultiSourceLoadError::UnsupportedDetector, mftSource.id};
  }

  // Stage normalized ownership first with the existing multi-source loader.
  // It preserves per-surface ownership and global SurfaceId mappings.
  MultiSourceFrame normalized;
  const std::array<ClusterSourceInput, 2> sources{itsSource, mftSource};
  const auto normalizedResult = loadSources(normalized, catalog, sources, origin);
  if (!normalizedResult.ok()) {
    return normalizedResult;
  }

  LegacyTrackerScratchITS stagedITS;
  LegacyTrackerScratchMFT stagedMFT;
  prepareStagingScratch(stagedITS, itsScratch);
  prepareStagingScratch(stagedMFT, mftScratch);
  if (const auto result = stageScratch(stagedITS, itsSource, catalog, origin); !result.ok()) {
    return result;
  }
  if (const auto result = stageScratch(stagedMFT, mftSource, catalog, origin); !result.ok()) {
    return result;
  }

  // Nothing above mutates a caller-owned object. This final sequence cannot
  // throw; once it starts, all three owners become the one staged event.
  frame.commitNormalizedFrame(std::move(normalized));
  commitNormalizedBackfill(itsScratch, stagedITS);
  commitNormalizedBackfill(mftScratch, stagedMFT);
  return normalizedResult;
}

void MultiSourceTimeFrameLoader::resetITSAndMFTEvent(TimeFrame& frame,
                                                     LegacyTrackerScratchITS& itsScratch,
                                                     LegacyTrackerScratchMFT& mftScratch) noexcept
{
  itsScratch.resetScratch();
  mftScratch.resetScratch();
  frame.wipe();
}

} // namespace o2::itsmft::tracking

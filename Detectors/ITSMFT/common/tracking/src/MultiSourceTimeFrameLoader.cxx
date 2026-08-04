// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"

#include <array>
#include <vector>

#include "ITSMFTTracking/TimeFrame.h"

namespace o2::itsmft::tracking
{

template <int NLayers>
LoadSourcesResult MultiSourceTimeFrameLoader::LoadTargetImpl<NLayers>::stage(
  const ClusterSourceInput& source, SurfaceCatalogView catalog, const o2::InteractionRecord& origin)
{
  // Keep allocator identity, not merely allocator configuration: pmr swaps
  // in commit() are defined only for equal resources.
  mStaged.mExternalAllocator = mLive.mExternalAllocator;
  mStaged.mExtMemoryPool = mLive.mExtMemoryPool;
  mStaged.setMemoryPool(mLive.mMemoryPool);

  // The single-source legacy conversion, invoked only on a disposable
  // TimeFrame -- the shared owner's normalized frame is staged once,
  // generically, by loadEvent() itself, from every binding's source.
  TimeFrame disposable;
  ClusterSourceInput one = source;
  one.id = ClusterSourceId{0};
  return mStaged.loadNormalizedSource(disposable, *one.decoder, origin, one.timing,
                                      one.clusters, one.patterns, one.rofs,
                                      one.dictionary, one.labels, one.detector,
                                      one.layerToSurface, catalog, one.applySysErrors);
}

template <int NLayers>
void MultiSourceTimeFrameLoader::LoadTargetImpl<NLayers>::commit() noexcept
{
  // Same-allocator swaps only -- never throws (see stage()'s allocator-
  // identity binding above).
  for (int layer = 0; layer < NLayers; ++layer) {
    mLive.mUnsortedClusters[layer].swap(mStaged.mUnsortedClusters[layer]);
    mLive.mTrackingFrameInfo[layer].swap(mStaged.mTrackingFrameInfo[layer]);
    mLive.mClusterExternalIndices[layer].swap(mStaged.mClusterExternalIndices[layer]);
    mLive.mClusterSize[layer].swap(mStaged.mClusterSize[layer]);
    mLive.mROFramesClusters[layer].swap(mStaged.mROFramesClusters[layer]);
    mLive.mClusterLabels[layer] = mStaged.mClusterLabels[layer];
  }
  for (int i = 0; i < 2; ++i) {
    mLive.mNTrackletsPerCluster[i].swap(mStaged.mNTrackletsPerCluster[i]);
    mLive.mNTrackletsPerClusterSum[i].swap(mStaged.mNTrackletsPerClusterSum[i]);
  }
}

template class MultiSourceTimeFrameLoader::LoadTargetImpl<ITSNLayers>;
template class MultiSourceTimeFrameLoader::LoadTargetImpl<MFTNLayers>;

// M6d: LoadTargetImplSurface mirrors LoadTargetImpl<NLayers>::stage()/
// commit() exactly (same allocator-identity-preservation pattern, same
// mStaged.loadNormalizedSource() call, same swap-commit loop), just over a
// runtime orderedSurfaces.size() bound instead of a fixed NLayers one.
LoadSourcesResult MultiSourceTimeFrameLoader::LoadTargetImplSurface::stage(
  const ClusterSourceInput& source, SurfaceCatalogView catalog, const o2::InteractionRecord& origin)
{
  mStaged.mExternalAllocator = mLive.mExternalAllocator;
  mStaged.mExtMemoryPool = mLive.mExtMemoryPool;
  mStaged.setMemoryPool(mLive.getMemoryPool());
  // mStaged is a fresh SurfaceTrackingScratch, never itself through
  // adoptPlan() (unlike LegacyTrackerScratch<NLayers>'s own mStaged, whose
  // Group A member is a fixed-size std::array<T, NLayers> valid at every
  // index from default construction) -- its Group A outer
  // std::vector<bounded_vector<T>> members must be sized to mLive's own
  // owned-surface count before loadNormalizedSource() below indexes them.
  // Group B/D counts are irrelevant here (loadNormalizedSource() never
  // touches them), so 0/0 is deliberate, not a placeholder oversight.
  mStaged.adoptPlan(mLive.getNOwnedSurfaces(), 0, 0);

  TimeFrame disposable;
  ClusterSourceInput one = source;
  one.id = ClusterSourceId{0};
  return mStaged.loadNormalizedSource(disposable, *one.decoder, origin, one.timing,
                                      one.clusters, one.patterns, one.rofs,
                                      one.dictionary, one.labels, one.detector,
                                      one.layerToSurface, catalog, one.applySysErrors);
}

void MultiSourceTimeFrameLoader::LoadTargetImplSurface::commit() noexcept
{
  // Deliberately NOT SurfaceTrackingScratch::swap() (a whole-object swap
  // that would also exchange Group B/D and the adopted plan size):
  // `mStaged` was never through adoptPlan(), so a whole-object swap would
  // clobber `mLive`'s already-populated Group B/D/plan-size state with
  // `mStaged`'s empty defaults. Same-allocator swaps only, same narrow
  // Group-A scope as LoadTargetImpl<NLayers>::commit() above -- never
  // throws (see stage()'s allocator-identity binding above).
  for (std::size_t layer = 0; layer < mLive.getNOwnedSurfaces(); ++layer) {
    mLive.mUnsortedClusters[layer].swap(mStaged.mUnsortedClusters[layer]);
    mLive.mTrackingFrameInfo[layer].swap(mStaged.mTrackingFrameInfo[layer]);
    mLive.mClusterExternalIndices[layer].swap(mStaged.mClusterExternalIndices[layer]);
    mLive.mClusterSize[layer].swap(mStaged.mClusterSize[layer]);
    mLive.mROFramesClusters[layer].swap(mStaged.mROFramesClusters[layer]);
    mLive.mClusterLabels[layer] = mStaged.mClusterLabels[layer];
  }
  for (int i = 0; i < 2; ++i) {
    mLive.mNTrackletsPerCluster[i].swap(mStaged.mNTrackletsPerCluster[i]);
    mLive.mNTrackletsPerClusterSum[i].swap(mStaged.mNTrackletsPerClusterSum[i]);
  }
}

LoadSourcesResult MultiSourceTimeFrameLoader::loadEvent(TimeFrame& frame, gsl::span<const AtomicLoadBinding> bindings,
                                                        SurfaceCatalogView catalog, const o2::InteractionRecord& origin)
{
  std::vector<ClusterSourceInput> sources;
  sources.reserve(bindings.size());
  for (const auto& binding : bindings) {
    sources.push_back(binding.source);
  }

  // Stage normalized ownership first, generically, over every source at
  // once: loadSources() already preserves per-surface ownership and
  // global SurfaceId mappings regardless of source count -- untouched by
  // this milestone.
  MultiSourceFrame normalized;
  const auto normalizedResult = loadSources(normalized, catalog, sources, origin);
  if (!normalizedResult.ok()) {
    return normalizedResult;
  }

  for (const auto& binding : bindings) {
    if (const auto result = binding.target.stage(binding.source, catalog, origin); !result.ok()) {
      // Nothing has been committed yet: `frame` and every binding's live
      // scratch (including any binding staged successfully before this
      // one) are all still exactly as they were before this call.
      return result;
    }
  }

  // Nothing above mutated a caller-owned object. This final sequence
  // cannot throw; once it starts, every owner becomes the one staged
  // event.
  frame.commitNormalizedFrame(std::move(normalized));
  for (const auto& binding : bindings) {
    binding.target.commit();
  }
  return normalizedResult;
}

LoadSourcesResult MultiSourceTimeFrameLoader::loadITSAndMFT(TimeFrame& frame,
                                                            LegacyTrackerScratchITS& itsScratch,
                                                            LegacyTrackerScratchMFT& mftScratch,
                                                            const ClusterSourceInput& itsSource,
                                                            const ClusterSourceInput& mftSource,
                                                            SurfaceCatalogView catalog,
                                                            const o2::InteractionRecord& origin)
{
  // Thin compatibility wrapper (M2b): the fixed ITS=0/MFT=1 position
  // contract lives only here, not in the generic loadEvent() transaction
  // this forwards to. Keeping the source positions fixed makes source
  // metadata and every ClusterRef stable and makes it impossible to
  // accidentally backfill a detector into the other scratch.
  if (itsSource.id != ClusterSourceId{0} || itsSource.detector != o2::detectors::DetID::ITS) {
    return {MultiSourceLoadError::UnsupportedDetector, itsSource.id};
  }
  if (mftSource.id != ClusterSourceId{1} || mftSource.detector != o2::detectors::DetID::MFT) {
    return {MultiSourceLoadError::UnsupportedDetector, mftSource.id};
  }

  LoadTargetImpl<ITSNLayers> itsTarget{itsScratch};
  LoadTargetImpl<MFTNLayers> mftTarget{mftScratch};
  const std::array<AtomicLoadBinding, 2> bindings{
    AtomicLoadBinding{itsSource, itsTarget},
    AtomicLoadBinding{mftSource, mftTarget}};
  return loadEvent(frame, gsl::span<const AtomicLoadBinding>{bindings}, catalog, origin);
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

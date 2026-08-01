// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
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
/// \file LegacyTrackerScratch.cxx
/// \brief
///

#include <cmath>
#include <limits>
#include <new>
#include <numeric>
#include <stdexcept>
#include <string>

#include "Framework/Logger.h"
#include "ITSMFTTracking/IndexTableConfiguration.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/LegacyTrackerScratch.h"
#include "ITStracking/MathUtils.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "MFTTracking/Constants.h"

namespace
{
struct ClusterHelper {
  float rowCoord;
  float r;
  int bin;
  int ind;
};

// Implementation detail of LegacyTrackerScratch::initialise()/
// prepareClusters()'s checked index-table allocation-size arithmetic
// (checkedIndexTableSizeProduct, ITSMFTTracking/IndexTableConfiguration.h).
// Deriving from std::bad_alloc -- rather than TrackerTraits.h's
// TraversalException -- is deliberate: an unrepresentable allocation size
// here is a function of this scratch's actual loaded ROF volume combined
// with an otherwise-valid bin configuration, the same per-TF resource-
// failure category as a genuine bad_alloc/BoundedMemoryResource::
// MemoryLimitExceeded, not a structural TrackingParameters problem. This
// keeps it caught (and DropTFUponFailure-gated) by the existing
// `catch (const std::bad_alloc&)` handling in CATracker.cxx/
// TrackingInterface.cxx without this file depending on TrackerTraits.h. Only
// its std::bad_alloc base-class behavior is part of the failure contract;
// the type itself is not part of any public API.
class IndexTableAllocationSizeError final : public std::bad_alloc
{
 public:
  const char* what() const noexcept override { return "index-table allocation size is unrepresentable"; }
};
} // namespace

namespace o2::itsmft::tracking
{

#ifndef GPUCA_GPUCODE
namespace
{
// Internal invariant violation: loadNormalizedSource() staged a legacy
// backfill vector (mUnsortedClusters/mTrackingFrameInfo/
// mClusterExternalIndices/mClusterSize/mROFramesClusters, or the non-per-
// layer mNTrackletsPerCluster/mNTrackletsPerClusterSum, reported with the
// sentinel layer -1) with a memory-resource pointer different from its
// corresponding live scratch vector. Correctly configured operation always
// derives both the staged and live pointers from the same mMemoryPool.get()
// calls, so this can never be triggered by caller input; it is a defensive
// internal-invariant gate checked unconditionally right before commit, not a
// validation outcome reported through LoadSourcesResult. Internal to this
// translation unit: not part of LegacyTrackerScratch's public API.
class NormalizedBackfillAllocatorMismatch final : public std::logic_error
{
 public:
  explicit NormalizedBackfillAllocatorMismatch(int layer)
    : std::logic_error("LegacyTrackerScratch::loadNormalizedSource(): staged/live memory-resource mismatch on layer " + std::to_string(layer))
  {
  }
};

} // namespace
#endif

using o2::its::clearResizeBoundedVector;
using o2::its::deepVectorClear;
namespace math_utils = o2::its::math_utils;
using o2::itsmft::IndexTableCoordType;
using o2::itsmft::IterationStep;
using o2::itsmft::TrackingParameters;

template <int NLayers>
void LegacyTrackerScratch<NLayers>::loadROFrameData(gsl::span<const o2::itsmft::ROFRecord> rofs,
                                                    gsl::span<const itsmft::CompClusterExt> clusters,
                                                    gsl::span<const unsigned char>::iterator& pattIt,
                                                    const itsmft::TopologyDictionary* dict,
                                                    int layer,
                                                    const dataformats::MCTruthContainer<MCCompLabel>* mcLabels,
                                                    o2::detectors::DetID::ID detId)
{
  if (NLayers != constants::nLayersForDet(detId)) {
    LOGP(fatal, "LegacyTrackerScratch<{}> is incompatible with detector {} (expected {} layers)",
         NLayers, static_cast<int>(detId), constants::nLayersForDet(detId));
  }
  ioutils::fillMatrixCache(detId);
  resetROFrameData(layer);
  prepareROFrameData(clusters, layer, detId);

  // check for missing/empty/unset rofs
  // the code requires consistent monotonically increasing input without gaps
  const auto& timing = mROFOverlapTableView.getLayer(layer >= 0 ? layer : 0);
  if (timing.mNROFsTF != rofs.size()) {
    LOGP(fatal, "Received inconsistent number of rofs on layer:{} expected:{} received:{}", layer, timing.mNROFsTF, rofs.size());
  }

  for (int32_t iRof{0}; iRof < rofs.size(); ++iRof) {
    const auto& rof = rofs[iRof];
    for (int clusterId{rof.getFirstEntry()}; clusterId < rof.getFirstEntry() + rof.getNEntries(); ++clusterId) {
      const auto& c = clusters[clusterId];
      int lay{0};
      unsigned int clusterSize{0};
      TrackingFrameInfo tfInfo;
      if (detId == o2::detectors::DetID::MFT) {
        o2::itsmft::ioutils::loadClusterTrackingFrameInfo<o2::detectors::DetID::MFT>(c, pattIt, dict, lay, clusterSize, tfInfo);
      } else {
        o2::itsmft::ioutils::loadClusterTrackingFrameInfo<o2::detectors::DetID::ITS>(c, pattIt, dict, lay, clusterSize, tfInfo);
      }
      mClusterSize[layer >= 0 ? layer : 0][clusterId] = std::clamp(clusterSize, 0u, 255u);
      addTrackingFrameInfoToLayer(lay, tfInfo);
      addClusterToLayer(lay, tfInfo.xCoordinate, tfInfo.yCoordinate, tfInfo.zCoordinate, mUnsortedClusters[lay].size());
      addClusterExternalIndexToLayer(lay, clusterId);
    }
    // effectively calculating an exclusive sum
    if (layer >= 0) {
      mROFramesClusters[layer][iRof + 1] = mUnsortedClusters[layer].size();
    } else {
      for (unsigned int iL{0}; iL < mUnsortedClusters.size(); ++iL) {
        mROFramesClusters[iL][iRof + 1] = mUnsortedClusters[iL].size();
      }
    }
  }

  if (layer == 1 || layer == -1) {
    for (auto i = 0; i < mNTrackletsPerCluster.size(); ++i) {
      mNTrackletsPerCluster[i].resize(mUnsortedClusters[1].size());
      mNTrackletsPerClusterSum[i].resize(mUnsortedClusters[1].size() + 1);
    }
  }

  if (mcLabels != nullptr) {
    mClusterLabels[layer >= 0 ? layer : 0] = mcLabels;
  } else {
    mClusterLabels[layer >= 0 ? layer : 0] = nullptr;
  }
}

#ifndef GPUCA_GPUCODE
template <int NLayers>
LoadSourcesResult LegacyTrackerScratch<NLayers>::loadNormalizedSource(
  TimeFrame& frame,
  const ClusterDecoder& decoder,
  const o2::InteractionRecord& origin,
  const ROFTimingConfig& timing,
  gsl::span<const itsmft::CompClusterExt> clusters,
  gsl::span<const unsigned char> patterns,
  gsl::span<const o2::itsmft::ROFRecord> rofs,
  const itsmft::TopologyDictionary* dictionary,
  const dataformats::MCTruthContainer<MCCompLabel>* labels,
  o2::detectors::DetID::ID detId,
  gsl::span<const SurfaceId> orderedSurfaces,
  SurfaceCatalogView catalogView,
  bool applySysErrors)
{
  // Exactly one source is ever submitted here, and loadSources() requires
  // dense, zero-based source IDs, so ClusterSourceId{0} is the only value
  // that could ever succeed; it is fixed internally rather than accepted as
  // a parameter that would misleadingly suggest other IDs are supported.
  constexpr ClusterSourceId kSourceId{0};

  // nLayersForDet() maps every non-MFT detector to ITSNLayers, so it alone
  // cannot distinguish ITS from an unsupported detector that happens to
  // share NLayers==7; the detector identity must be checked explicitly and
  // first, before nLayersForDet() or the caller-supplied catalog/mapping.
  if (detId != o2::detectors::DetID::ITS && detId != o2::detectors::DetID::MFT) {
    return {MultiSourceLoadError::UnsupportedDetector, kSourceId};
  }
  if (NLayers != constants::nLayersForDet(detId)) {
    return {MultiSourceLoadError::UnsupportedDetector, kSourceId};
  }
  // Neither this scratch nor `frame` owns a catalog/plan of its own (Gate 4
  // B2 Slice 2): the caller (ITSMFTTrackingInterface, the plan's one owner)
  // supplies both explicitly. An empty view means the caller has no plan at
  // all.
  if (catalogView.surfaces == nullptr || catalogView.nSurfaces == 0) {
    return {MultiSourceLoadError::SurfaceCatalogNotConfigured, kSourceId};
  }
  if (orderedSurfaces.size() != static_cast<size_t>(NLayers)) {
    return {MultiSourceLoadError::InvalidLayerMapping, kSourceId};
  }
  // Fixed-size, non-allocating duplicate check: the caller-supplied catalog
  // is already guaranteed (by construction -- see buildDetectorLayoutSet())
  // to fit within the fixed 32-bit SurfaceMask capacity, so this preflight
  // step still performs no allocation before decoding/mutation.
  SurfaceMask mappedSurfaceSeen{};
  for (const auto& surfaceId : orderedSurfaces) {
    if (!surfaceId.isValid() || surfaceId.value() >= catalogView.nSurfaces) {
      return {MultiSourceLoadError::InvalidLayerMapping, kSourceId};
    }
    if (mappedSurfaceSeen.has(surfaceId)) {
      return {MultiSourceLoadError::InvalidLayerMapping, kSourceId};
    }
    mappedSurfaceSeen.set(surfaceId);
  }
  for (const auto& surfaceId : orderedSurfaces) {
    if (catalogView.getSurface(surfaceId).detectorId != static_cast<uint8_t>(detId)) {
      return {MultiSourceLoadError::DetectorSurfaceMismatch, kSourceId};
    }
  }

  const gsl::span<const SurfaceId> layerToSurface = orderedSurfaces;

  // Stage into a scratch owner: loadSources() itself never mutates its
  // `frame` argument on failure, but staging separately also protects
  // `frame`'s *existing* normalized owner (and, by construction below, every
  // legacy compatibility structure here) from a failed reload.
  MultiSourceFrame staged;
  ClusterSourceInput src;
  src.id = kSourceId;
  src.detector = detId;
  src.clusters = clusters;
  src.patterns = patterns;
  src.rofs = rofs;
  src.dictionary = dictionary;
  src.labels = labels;
  src.layerToSurface = layerToSurface;
  src.timing = timing;
  src.decoder = &decoder;
  src.applySysErrors = applySysErrors; // Matches loadROFrameData()'s own default (loadClusterTrackingFrameInfo<DetId>(..., applySysErrors=true)).

  const auto result = loadSources(staged, catalogView, gsl::span<const ClusterSourceInput>(&src, 1), origin);
  if (!result.ok()) {
    // Nothing below has run: `frame`'s normalized owner and every legacy
    // compatibility structure here are exactly as they were before this
    // call.
    return result;
  }

  // From here on, loadSources() has already succeeded. Strong exception
  // safety for the legacy backfill: every field below is first built on a
  // complete set of *local* staged owners -- constructed with the same
  // memory-resource pointers as the live vectors they will replace -- and
  // measurements are read from the local `staged` normalized frame, never
  // from `frame`'s live one. `frame` and every legacy member here are left
  // untouched until the allocator invariant gate below has passed, so an
  // allocation failure anywhere in this staging (the only way this can
  // throw, e.g. BoundedMemoryResource::MemoryLimitExceeded) unwinds the stack
  // through plain local-variable destruction and leaves every live
  // representation exactly as it was before this call.
  const bool isMFT = (detId == o2::detectors::DetID::MFT);
  auto* pool = mMemoryPool.get();
  const auto nROFs = static_cast<size_t>(rofs.size());

  // Staged owners, one entry per layer, built below with std::vector<T>::
  // emplace_back so each element is constructed in place with its final
  // memory-resource pointer directly -- never by default-constructing (with
  // whatever default allocator that implies) and then swapping in a
  // differently-allocated replacement. bounded_vector is std::pmr::vector,
  // whose swap() is only well-defined when both vectors' allocators compare
  // equal (std::pmr::polymorphic_allocator neither propagates on swap nor is
  // ever always_equal); swapping unequal-allocator vectors is undefined
  // behavior, and in practice this libc++ resolves it by keeping each side's
  // own original allocator, silently discarding the intended one. Each
  // vector is reserved to NLayers up front so growing it one emplace_back
  // per layer can never reallocate (and thus never move/invalidate) an
  // already-staged element.
  std::vector<bounded_vector<Cluster>> stagedUnsortedClusters;
  std::vector<bounded_vector<TrackingFrameInfo>> stagedTrackingFrameInfo;
  std::vector<bounded_vector<int>> stagedClusterExternalIndices;
  std::vector<bounded_vector<uint8_t>> stagedClusterSize;
  std::vector<bounded_vector<int>> stagedROFramesClusters;
  std::array<const dataformats::MCTruthContainer<MCCompLabel>*, NLayers> stagedClusterLabels{nullptr};
  stagedUnsortedClusters.reserve(NLayers);
  stagedTrackingFrameInfo.reserve(NLayers);
  stagedClusterExternalIndices.reserve(NLayers);
  stagedClusterSize.reserve(NLayers);
  stagedROFramesClusters.reserve(NLayers);

  for (int layer = 0; layer < NLayers; ++layer) {
    const auto measurements = staged.getSurfaceMeasurements(layerToSurface[layer]);

    // mr/pool are nullptr whenever no framework allocator and no memory pool
    // have been configured (setMemoryPool() was never called); a nullptr
    // std::pmr::memory_resource* is not a valid polymorphic_allocator target
    // (constructing a non-empty vector through it would dereference the
    // nullptr). Falling back to the *live* vector's own current allocator
    // resource -- exactly what deepVectorClear()/clearResizeBoundedVector()
    // do for the same nullptr case -- avoids that and, since the
    // corresponding live vector has not been touched yet, is exactly the
    // resource the invariant gate below requires the staged vector to
    // match.
    auto* mr = getMaybeFrameworkHostResource();
    auto* unsortedClustersMr = mr != nullptr ? mr : mUnsortedClusters[layer].get_allocator().resource();
    auto* trackingFrameInfoMr = mr != nullptr ? mr : mTrackingFrameInfo[layer].get_allocator().resource();
    auto* clusterExternalIndicesMr = pool != nullptr ? pool : mClusterExternalIndices[layer].get_allocator().resource();
    auto* clusterSizeMr = pool != nullptr ? pool : mClusterSize[layer].get_allocator().resource();
    auto* rofFramesClustersMr = mr != nullptr ? mr : mROFramesClusters[layer].get_allocator().resource();

    stagedUnsortedClusters.emplace_back(std::pmr::polymorphic_allocator<Cluster>{unsortedClustersMr});
    stagedTrackingFrameInfo.emplace_back(std::pmr::polymorphic_allocator<TrackingFrameInfo>{trackingFrameInfoMr});
    stagedClusterExternalIndices.emplace_back(std::pmr::polymorphic_allocator<int>{clusterExternalIndicesMr});
    stagedClusterSize.emplace_back(measurements.size(), uint8_t{0}, std::pmr::polymorphic_allocator<uint8_t>{clusterSizeMr});
    stagedROFramesClusters.emplace_back(nROFs + 1, 0, std::pmr::polymorphic_allocator<int>{rofFramesClustersMr});

    size_t mi{0};
    for (const auto& m : measurements) {
      TrackingFrameInfo tfInfo;
      if (isMFT) {
        // Recreate the established synthetic legacy MFT representation
        // (TrackingFrameInfoAdapters.h::makeTrackingFrameInfo<MFT>) from the
        // normalized global position and row/column covariance. This is
        // deliberately not m.frame, which for a Disk-kind SurfaceMeasurement
        // holds the (z, x, y) disk-frame projection, not the legacy
        // synthetic layout existing production code consumes. The normalized
        // covariance is authoritative: decoding/loading has already applied
        // detector systematics when requested, and this compatibility
        // backfill copies that covariance without adding to it.
        tfInfo = TrackingFrameInfo{
          m.global.x, m.global.y, m.global.z,
          m.global.x, 0.f,
          std::array<float, 2>{m.global.y, m.global.z},
          std::array<float, 3>{m.covariance.uu, m.covariance.uv, m.covariance.vv}};
      } else {
        // As for MFT above, compatibility storage copies the authoritative
        // normalized covariance produced by decoding/loading.
        tfInfo = TrackingFrameInfo{
          m.global.x, m.global.y, m.global.z,
          m.frame.q, m.frame.frameAngle,
          std::array<float, 2>{m.frame.u, m.frame.v},
          std::array<float, 3>{m.covariance.uu, m.covariance.uv, m.covariance.vv}};
      }
      stagedTrackingFrameInfo[layer].push_back(tfInfo);
      stagedUnsortedClusters[layer].emplace_back(m.global.x, m.global.y, m.global.z, static_cast<int>(stagedUnsortedClusters[layer].size()));
      stagedClusterExternalIndices[layer].push_back(static_cast<int>(m.cluster.index));
      stagedClusterSize[layer][mi] = static_cast<uint8_t>(std::clamp(m.shape.nPixels, 0u, 255u));
      ++mi;
    }

    size_t mj{0};
    for (size_t r = 0; r < nROFs; ++r) {
      while (mj < measurements.size() && measurements[mj].sourceROF == static_cast<uint32_t>(r)) {
        ++mj;
      }
      stagedROFramesClusters[layer][r + 1] = static_cast<int>(mj);
    }

    stagedClusterLabels[layer] = labels;
  }

  // mNTrackletsPerCluster/mNTrackletsPerClusterSum are not per-layer (see
  // their declaration: std::array<bounded_vector<int>, 2>, sized from layer
  // 1's cluster count, consumed later by computeTrackletsPerROFScans() on
  // the CA hot path). loadROFrameData() has always sized them as part of its
  // own commit; loadNormalizedSource() must too, or a scratch loaded through
  // this path silently carries stale/empty tracklet-per-cluster scratch into
  // tracking. Sized from the *staged* (not yet committed) layer-1 count, so
  // this stays consistent with everything else staged above.
  const size_t nClustersLayer1 = stagedUnsortedClusters[1].size();
  auto* nTrackletsPerClusterMr = pool != nullptr ? pool : mNTrackletsPerCluster[0].get_allocator().resource();
  auto* nTrackletsPerClusterSumMr = pool != nullptr ? pool : mNTrackletsPerClusterSum[0].get_allocator().resource();
  // List-initialized at declaration, not assigned: each element is
  // *constructed* directly with its final allocator. std::pmr::vector
  // move-assignment falls back to an element-wise move when allocators
  // differ (polymorphic_allocator neither propagates on move-assignment nor
  // is ever always_equal), which would silently discard the intended pool --
  // the same hazard the per-layer vectors above avoid via emplace_back.
  std::array<bounded_vector<int>, 2> stagedNTrackletsPerCluster{
    bounded_vector<int>(nClustersLayer1, 0, std::pmr::polymorphic_allocator<int>{nTrackletsPerClusterMr}),
    bounded_vector<int>(nClustersLayer1, 0, std::pmr::polymorphic_allocator<int>{nTrackletsPerClusterMr})};
  std::array<bounded_vector<int>, 2> stagedNTrackletsPerClusterSum{
    bounded_vector<int>(nClustersLayer1 + 1, 0, std::pmr::polymorphic_allocator<int>{nTrackletsPerClusterSumMr}),
    bounded_vector<int>(nClustersLayer1 + 1, 0, std::pmr::polymorphic_allocator<int>{nTrackletsPerClusterSumMr})};

  // Invariant gate: every staged vector must share its live counterpart's
  // memory resource before anything is committed. Construction above always
  // derives both pointers from the same mMemoryPool.get() calls, so this is
  // unreachable through caller input; a mismatch here means the staging
  // above is internally inconsistent, which is reported as an internal logic
  // error rather than a LoadSourcesResult error code. This check -- and
  // everything above it -- must run before `frame` or any legacy member here
  // is touched.
  for (int layer = 0; layer < NLayers; ++layer) {
    if (mUnsortedClusters[layer].get_allocator().resource() != stagedUnsortedClusters[layer].get_allocator().resource() ||
        mTrackingFrameInfo[layer].get_allocator().resource() != stagedTrackingFrameInfo[layer].get_allocator().resource() ||
        mClusterExternalIndices[layer].get_allocator().resource() != stagedClusterExternalIndices[layer].get_allocator().resource() ||
        mClusterSize[layer].get_allocator().resource() != stagedClusterSize[layer].get_allocator().resource() ||
        mROFramesClusters[layer].get_allocator().resource() != stagedROFramesClusters[layer].get_allocator().resource()) {
      throw NormalizedBackfillAllocatorMismatch(layer);
    }
  }
  // mNTrackletsPerCluster/mNTrackletsPerClusterSum are not indexed by layer:
  // every one of their two elements is checked individually, not just
  // element zero, since each was constructed as an independent bounded_vector.
  for (int i = 0; i < 2; ++i) {
    if (mNTrackletsPerCluster[i].get_allocator().resource() != stagedNTrackletsPerCluster[i].get_allocator().resource() ||
        mNTrackletsPerClusterSum[i].get_allocator().resource() != stagedNTrackletsPerClusterSum[i].get_allocator().resource()) {
      throw NormalizedBackfillAllocatorMismatch(-1); // -1: not a per-layer container; see declaration comment above.
    }
  }

  // Commit: nothing past this point can throw. frame.commitNormalizedFrame()
  // is documented noexcept (its own move-assignment is statically nothrow,
  // asserted there); every legacy container exchange below is a same-
  // allocator bounded_vector::swap (the gate above just proved the
  // allocators equal, which is the vector swap's own precondition for
  // well-defined, allocation-free behavior); the per-layer label pointers
  // are raw-pointer assignments.
  static_assert(noexcept(std::declval<bounded_vector<Cluster>&>().swap(std::declval<bounded_vector<Cluster>&>())));
  // mNTrackletsPerCluster/mNTrackletsPerClusterSum are bounded_vector<int>,
  // a distinct instantiation from bounded_vector<Cluster> above; asserted
  // separately rather than assumed to swap noexcept by analogy.
  static_assert(noexcept(std::declval<bounded_vector<int>&>().swap(std::declval<bounded_vector<int>&>())));

  frame.commitNormalizedFrame(std::move(staged));
  for (int layer = 0; layer < NLayers; ++layer) {
    mUnsortedClusters[layer].swap(stagedUnsortedClusters[layer]);
    mTrackingFrameInfo[layer].swap(stagedTrackingFrameInfo[layer]);
    mClusterExternalIndices[layer].swap(stagedClusterExternalIndices[layer]);
    mClusterSize[layer].swap(stagedClusterSize[layer]);
    mROFramesClusters[layer].swap(stagedROFramesClusters[layer]);
    mClusterLabels[layer] = stagedClusterLabels[layer];
  }
  for (int i = 0; i < 2; ++i) {
    mNTrackletsPerCluster[i].swap(stagedNTrackletsPerCluster[i]);
    mNTrackletsPerClusterSum[i].swap(stagedNTrackletsPerClusterSum[i]);
  }

  return result;
}
#endif

template <int NLayers>
void LegacyTrackerScratch<NLayers>::resetROFrameData(int layer)
{
  if (layer >= 0) {
    deepVectorClear(mUnsortedClusters[layer], getMaybeFrameworkHostResource());
    deepVectorClear(mTrackingFrameInfo[layer], getMaybeFrameworkHostResource());
    deepVectorClear(mClusterExternalIndices[layer], mMemoryPool.get());
    clearResizeBoundedVector(mROFramesClusters[layer], mROFOverlapTableView.getLayer(layer).mNROFsTF + 1, getMaybeFrameworkHostResource());
  } else {
    for (int iLayer{0}; iLayer < NLayers; ++iLayer) {
      deepVectorClear(mUnsortedClusters[iLayer], getMaybeFrameworkHostResource());
      deepVectorClear(mTrackingFrameInfo[iLayer], getMaybeFrameworkHostResource());
      deepVectorClear(mClusterExternalIndices[iLayer], mMemoryPool.get());
      clearResizeBoundedVector(mROFramesClusters[iLayer], mROFOverlapTableView.getLayer(iLayer).mNROFsTF + 1, getMaybeFrameworkHostResource());
    }
  }
}

template <int NLayers>
void LegacyTrackerScratch<NLayers>::prepareROFrameData(gsl::span<const itsmft::CompClusterExt> clusters, int layer, o2::detectors::DetID::ID detId)
{
  if (layer >= 0) {
    mUnsortedClusters[layer].reserve(clusters.size());
    mTrackingFrameInfo[layer].reserve(clusters.size());
    mClusterExternalIndices[layer].reserve(clusters.size());
    clearResizeBoundedVector(mClusterSize[layer], clusters.size(), mMemoryPool.get());
  } else {
    clearResizeBoundedVector(mClusterSize[0], clusters.size(), mMemoryPool.get());
    std::array<size_t, NLayers> clusterCountPerLayer{0};
    for (const auto& cls : clusters) {
      ++clusterCountPerLayer[ioutils::getClusterLayer(detId, cls)];
    }
    for (int iLayer{0}; iLayer < NLayers; ++iLayer) {
      mUnsortedClusters[iLayer].reserve(clusterCountPerLayer[iLayer]);
      mTrackingFrameInfo[iLayer].reserve(clusterCountPerLayer[iLayer]);
      mClusterExternalIndices[iLayer].reserve(clusterCountPerLayer[iLayer]);
    }
  }
}

template <int NLayers>
void LegacyTrackerScratch<NLayers>::prepareClusters(const TimeFrame& frame, const TrackingParameters& trkParam, const int maxLayers,
                                                    const LayerMeasurementSpans<NLayers>& layerMeasurements)
{
  // numBins/stride come from the already-committed, already-validated
  // mIndexTableUtils -- never again from trkParam.RowBins/ColBins -- via
  // checked size_t arithmetic (see the identical reasoning in initialise()).
  const int colBinsCount = mIndexTableUtils.getNcolBins();
  std::size_t numBins = 0;
  if (!checkedIndexTableSizeProduct(static_cast<std::size_t>(mIndexTableUtils.getNrowBins()),
                                    static_cast<std::size_t>(colBinsCount), numBins) ||
      numBins == std::numeric_limits<std::size_t>::max()) {
    throw IndexTableAllocationSizeError{};
  }
  const std::size_t stride{numBins + 1};
  bounded_vector<ClusterHelper> cHelper(mMemoryPool.get());
  bounded_vector<int> clsPerBin(numBins, 0, mMemoryPool.get());
  bounded_vector<int> lutPerBin(numBins, 0, mMemoryPool.get());
  const std::array<float, 2> beamXY{frame.getBeamX(), frame.getBeamY()};
  for (int iLayer{0}, stopLayer = std::min(trkParam.NLayers, maxLayers); iLayer < stopLayer; ++iLayer) {
    for (int rof{0}; rof < getNrof(iLayer); ++rof) {
      if (!mROFMaskView.isROFEnabled(iLayer, rof)) {
        continue;
      }
      const auto& unsortedClusters{getUnsortedClustersOnLayer(rof, iLayer)};
      const int clustersNum{static_cast<int>(unsortedClusters.size())};
      auto* tableBase = mIndexTables[iLayer].data() + rof * stride;

      cHelper.resize(clustersNum);

      const bool useXYBinning = mIndexTableUtils.getCoordType() == IndexTableCoordType::XY;
      for (int iCluster{0}; iCluster < clustersNum; ++iCluster) {
        const Cluster& c = unsortedClusters[iCluster];
        const auto& measurement = layerMeasurements[iLayer][c.clusterId];
        ClusterHelper& h = cHelper[iCluster];

        const float x = measurement.global.x - (useXYBinning ? 0.f : beamXY[0]);
        const float y = measurement.global.y - (useXYBinning ? 0.f : beamXY[1]);
        const float z = measurement.global.z;

        const float rowCoord = useXYBinning ? measurement.global.y : math_utils::computePhi(x, y);
        const float colCoord = useXYBinning ? measurement.global.x : z;
        int colBin{mIndexTableUtils.getColBinIndex(iLayer, colCoord)};
        if (colBin < 0 || colBin >= colBinsCount) {
          colBin = std::clamp(colBin, 0, colBinsCount - 1);
          mBogusClusters[iLayer]++;
        }
        int bin = mIndexTableUtils.getBinIndex(colBin, mIndexTableUtils.getRowBinIndex(rowCoord));
        h.rowCoord = rowCoord;
        h.r = math_utils::hypot(x, y);
        mMinR[iLayer] = o2::gpu::GPUCommonMath::Min(h.r, mMinR[iLayer]);
        mMaxR[iLayer] = o2::gpu::GPUCommonMath::Max(h.r, mMaxR[iLayer]);
        h.bin = bin;
        h.ind = clsPerBin[bin]++;
      }
      std::exclusive_scan(clsPerBin.begin(), clsPerBin.end(), lutPerBin.begin(), 0);

      auto clusters2beSorted{getClustersOnLayer(rof, iLayer)};
      for (int iCluster{0}; iCluster < clustersNum; ++iCluster) {
        const ClusterHelper& h = cHelper[iCluster];
        Cluster& c = clusters2beSorted[lutPerBin[h.bin] + h.ind];

        c = unsortedClusters[iCluster];
        const auto& measurement = layerMeasurements[iLayer][c.clusterId];
        c.phi = useXYBinning ? math_utils::computePhi(measurement.global.x, measurement.global.y) : h.rowCoord;
        c.radius = h.r;
        c.indexTableBinIndex = h.bin;
      }
      std::copy_n(lutPerBin.data(), clsPerBin.size(), tableBase);
      std::fill_n(tableBase + clsPerBin.size(), stride - clsPerBin.size(), clustersNum);

      std::fill(clsPerBin.begin(), clsPerBin.end(), 0);
      cHelper.clear();
    }
  }
}

template <int NLayers>
void LegacyTrackerScratch<NLayers>::initVertexingTopology(const TrackingParameters& trkParam)
{
  mVertexingTopology.init(3, trkParam.MaxHoles, LayerMask{trkParam.HoleLayerMask});
}

template <int NLayers>
void LegacyTrackerScratch<NLayers>::initDefaultTrackingTopology(const TrackingParameters& trkParam, const int maxLayers)
{
  mDefaultTrackingTopology.init(maxLayers, trkParam.MaxHoles, LayerMask{trkParam.HoleLayerMask});
}

template <int NLayers>
void LegacyTrackerScratch<NLayers>::initTrackerTopologies(gsl::span<const TrackingParameters> trkParams, const int maxLayers)
{
  mTrackerTopologies.resize(trkParams.size());
  for (size_t iteration = 0; iteration < trkParams.size(); ++iteration) {
    const int iterationMaxLayers = std::min(maxLayers, trkParams[iteration].NLayers);
    mTrackerTopologies[iteration].init(iterationMaxLayers, trkParams[iteration].MaxHoles, LayerMask{trkParams[iteration].HoleLayerMask});
  }
}

template <int NLayers>
void LegacyTrackerScratch<NLayers>::initialise(const TimeFrame& frame, const TrackingParameters& trkParam, const int maxLayers, const int iteration,
                                               const IndexTableUtilsN& indexTableConfig,
                                               const LayerMeasurementSpans<NLayers>& layerMeasurements)
{
  mTrackingTopologyView = iteration != constants::UnusedIndex ? mTrackerTopologies[iteration].getView() : (maxLayers == 3 ? mVertexingTopology.getView() : mDefaultTrackingTopology.getView());

  if (trkParam.PassFlags[IterationStep::FirstPass]) {
    deepVectorClear(mTracks);
    deepVectorClear(mTracksLabel);
    deepVectorClear(mLines);
    deepVectorClear(mLinesLabels);
    clearResizeBoundedVector(mLinesLabels, getNrof(1), mMemoryPool.get());
    // The caller (TrackerTraits::initialiseTimeFrame) has already bound and
    // validated this configuration -- from the transition-policy tag
    // resolved through the validated DetectorLayout/TransitionPolicyGrouping,
    // never from NLayers or DetID -- before this call. LegacyTrackerScratch
    // commits it by value and never inspects a tag or detector ID itself.
    mIndexTableUtils = indexTableConfig;
    clearResizeBoundedVector(mPositionResolution, trkParam.NLayers, mMemoryPool.get());
    clearResizeBoundedVector(mBogusClusters, trkParam.NLayers, mMemoryPool.get());
    deepVectorClear(mTrackletClusters);
    for (unsigned int iLayer{0}; iLayer < std::min((int)mClusters.size(), maxLayers); ++iLayer) {
      clearResizeBoundedVector(mClusters[iLayer], mUnsortedClusters[iLayer].size(), getMaybeFrameworkHostResource(maxLayers != NLayers));
      clearResizeBoundedVector(mUsedClusters[iLayer], mUnsortedClusters[iLayer].size(), getMaybeFrameworkHostResource(maxLayers != NLayers));
      mPositionResolution[iLayer] = o2::gpu::CAMath::Sqrt((0.5f * (trkParam.SystError2Col[iLayer] + trkParam.SystError2Row[iLayer])) + (trkParam.LayerResolution[iLayer] * trkParam.LayerResolution[iLayer]));
    }
    clearResizeBoundedVector(mLines, getNrof(1), mMemoryPool.get());
    clearResizeBoundedVector(mTrackletClusters, getNrof(1), mMemoryPool.get());

    // Sized from the just-committed, already-validated mIndexTableUtils --
    // never again from trkParam.RowBins/ColBins directly -- with checked
    // size_t arithmetic at every multiplication (RowBins*ColBins is already
    // bind-time bounded to fit int, but nROFs is a runtime quantity
    // independent of that bound, so the final product is checked here too
    // rather than assumed safe).
    std::size_t indexTableStride = 0;
    if (!checkedIndexTableSizeProduct(static_cast<std::size_t>(mIndexTableUtils.getNrowBins()),
                                      static_cast<std::size_t>(mIndexTableUtils.getNcolBins()), indexTableStride) ||
        indexTableStride == std::numeric_limits<std::size_t>::max()) {
      throw IndexTableAllocationSizeError{};
    }
    ++indexTableStride;
    for (int iLayer{0}; iLayer < NLayers; ++iLayer) {
      std::size_t indexTableSize = 0;
      if (!checkedIndexTableSizeProduct(static_cast<std::size_t>(getNrof(iLayer)), indexTableStride, indexTableSize)) {
        throw IndexTableAllocationSizeError{};
      }
      clearResizeBoundedVector(mIndexTables[iLayer], indexTableSize, getMaybeFrameworkHostResource());
    }

    // Cluster covariance is immutable after decoding/loading:
    // SurfaceMeasurement owns the authoritative value (including detector
    // systematics exactly once when requested), and TrackingFrameInfo is only
    // its compatibility copy. initialise() never mutates either
    // representation, so FirstPass, later iterations, LUT rebuild/reuse,
    // repeated initialisation, wipe and reload cannot accumulate systematic
    // covariance.

    mMinR.fill(std::numeric_limits<float>::max());
    mMaxR.fill(std::numeric_limits<float>::min());
  }
  clearResizeBoundedVector(mCells, mTrackingTopologyView.nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsLookupTable, mTrackingTopologyView.nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighbours, mTrackingTopologyView.nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighboursTopology, mTrackingTopologyView.nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighboursLUT, mTrackingTopologyView.nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellLabels, mTrackingTopologyView.nCells, mMemoryPool.get());
  clearResizeBoundedVector(mTracklets, mTrackingTopologyView.nTransitions, mMemoryPool.get());
  clearResizeBoundedVector(mTrackletLabels, mTrackingTopologyView.nTransitions, mMemoryPool.get());
  clearResizeBoundedVector(mTrackletsLookupTable, mTrackingTopologyView.nTransitions, mMemoryPool.get());
  clearResizeBoundedVector(mTransitionPhiCuts, mTrackingTopologyView.nTransitions, mMemoryPool.get());
  clearResizeBoundedVector(mTransitionMSAngles, mTrackingTopologyView.nTransitions, mMemoryPool.get());
  mNTrackletsPerROF.resize(2);
  for (auto& v : mNTrackletsPerROF) {
    v = bounded_vector<int>(getNrof(1) + 1, 0, mMemoryPool.get());
  }
  if (trkParam.PassFlags[IterationStep::RebuildClusterLUT]) {
    prepareClusters(frame, trkParam, maxLayers, layerMeasurements);
  }
  mTotalTracklets = {0, 0};
  if (maxLayers < trkParam.NLayers) { // Vertexer only, but in both iterations
    for (size_t iLayer{0}; iLayer < maxLayers; ++iLayer) {
      deepVectorClear(mUsedClusters[iLayer]);
      clearResizeBoundedVector(mUsedClusters[iLayer], mUnsortedClusters[iLayer].size(), mMemoryPool.get());
    }
  }

  // Per-layer position resolution (unchanged; multiple-scattering angle and
  // per-transition phi-cut/bending preparation are relocated to
  // TrackerTraits::initialiseTimeFrame() -- see TransitionPolicyOperations.h
  // layerMultipleScatteringAngle<Tag>/prepareTransitionScatteringAndBending).
  for (unsigned int iLayer{0}; iLayer < NLayers; ++iLayer) {
    mPositionResolution[iLayer] = o2::gpu::CAMath::Sqrt((0.5f * (trkParam.SystError2Col[iLayer] + trkParam.SystError2Row[iLayer])) + (trkParam.LayerResolution[iLayer] * trkParam.LayerResolution[iLayer]));
  }

  // Transition tracklet/label/LUT container clearing and resizing (unchanged;
  // both legacy branches were already identical here, so no NLayers/detector
  // branch is needed once the (relocated) value computation is removed).
  for (int transitionId{0}; transitionId < (int)mTracklets.size(); ++transitionId) {
    const auto& transition = mTrackingTopologyView.getTransition(transitionId);
    deepVectorClear(mTracklets[transitionId]);
    deepVectorClear(mTrackletLabels[transitionId]);
    deepVectorClear(mTrackletsLookupTable[transitionId]);
    mTrackletsLookupTable[transitionId].resize(mClusters[transition.fromLayer].size() + 1, 0);
  }

  for (int cellId{0}; cellId < (int)mCells.size(); ++cellId) {
    deepVectorClear(mCells[cellId]);
    deepVectorClear(mCellsLookupTable[cellId]);
    deepVectorClear(mCellsNeighbours[cellId]);
    deepVectorClear(mCellsNeighboursTopology[cellId]);
    deepVectorClear(mCellsNeighboursLUT[cellId]);
    deepVectorClear(mCellLabels[cellId]);
  }
}

template <int NLayers>
unsigned long LegacyTrackerScratch<NLayers>::getArtefactsMemory() const
{
  unsigned long size{0};
  for (const auto& trkl : mTracklets) {
    size += sizeof(Tracklet) * trkl.size();
  }
  for (const auto& cells : mCells) {
    size += sizeof(CellSeedN) * cells.size();
  }
  for (const auto& cellsN : mCellsNeighbours) {
    size += sizeof(int) * cellsN.size();
  }
  for (const auto& cellsN : mCellsNeighboursTopology) {
    size += sizeof(int) * cellsN.size();
  }
  return size;
}

template <int NLayers>
void LegacyTrackerScratch<NLayers>::printArtefactsMemory() const
{
  LOGP(info, "LegacyTrackerScratch: Artefacts occupy {:.2f} MB", getArtefactsMemory() / constants::MB);
}

template <int NLayers>
void LegacyTrackerScratch<NLayers>::computeTrackletsPerROFScans()
{
  for (ushort iLayer = 0; iLayer < 2; ++iLayer) {
    for (unsigned int iRof{0}; iRof < getNrof(1); ++iRof) {
      if (mROFMaskView.isROFEnabled(1, iRof)) {
        mTotalTracklets[iLayer] += mNTrackletsPerROF[iLayer][iRof];
      }
    }
    std::exclusive_scan(mNTrackletsPerROF[iLayer].begin(), mNTrackletsPerROF[iLayer].end(), mNTrackletsPerROF[iLayer].begin(), 0);
    std::exclusive_scan(mNTrackletsPerCluster[iLayer].begin(), mNTrackletsPerCluster[iLayer].end(), mNTrackletsPerClusterSum[iLayer].begin(), 0);
  }
}

template <int NLayers>
void LegacyTrackerScratch<NLayers>::setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
{
  mMemoryPool = pool;

  auto initVector = [&]<typename T>(bounded_vector<T>& vec, bool useExternal = false) {
    std::pmr::memory_resource* mr = (useExternal) ? mExtMemoryPool.get() : mMemoryPool.get();
    deepVectorClear(vec, mr);
  };

  auto initContainers = [&]<typename Container>(Container& container, bool useExternal = false) {
    for (auto& v : container) {
      initVector(v, useExternal);
    }
  };

  // these will only reside on the host for the cpu part
  initContainers(mClusterExternalIndices);
  initContainers(mNTrackletsPerCluster);
  initContainers(mNTrackletsPerClusterSum);
  initContainers(mNClustersPerROF);
  initVector(mTransitionPhiCuts);
  initVector(mTransitionMSAngles);
  initVector(mPositionResolution);
  initContainers(mClusterSize);
  initVector(mPValphaX);
  initVector(mBogusClusters);
  initContainers(mTrackletsIndexROF);
  initVector(mTracks);
  initContainers(mTracklets);
  initContainers(mCells);
  initContainers(mCellsNeighbours);
  initContainers(mCellsLookupTable);
  // MC info (we don't know if we have MC)
  initContainers(mLinesLabels);
  initContainers(mTrackletLabels);
  initContainers(mCellLabels);
  initVector(mTracksLabel);
  // these will use possibly an externally provided allocator
  initContainers(mClusters, hasFrameworkAllocator());
  initContainers(mUsedClusters, hasFrameworkAllocator());
  initContainers(mUnsortedClusters, hasFrameworkAllocator());
  initContainers(mIndexTables, hasFrameworkAllocator());
  initContainers(mTrackingFrameInfo, hasFrameworkAllocator());
  initContainers(mROFramesClusters, hasFrameworkAllocator());
}

template <int NLayers>
void LegacyTrackerScratch<NLayers>::setFrameworkAllocator(ExternalAllocator* ext)
{
  mExternalAllocator = ext;
  mExtMemoryPool = std::make_shared<BoundedMemoryResource>(mExternalAllocator);
}

template <int NLayers>
void LegacyTrackerScratch<NLayers>::resetScratch()
{
  deepVectorClear(mTracks);
  deepVectorClear(mTracklets);
  deepVectorClear(mCells);
  deepVectorClear(mCellsNeighbours);
  deepVectorClear(mCellsNeighboursTopology);
  deepVectorClear(mCellsLookupTable);
  deepVectorClear(mTrackletsLookupTable);
  deepVectorClear(mClusterExternalIndices);
  deepVectorClear(mNTrackletsPerCluster);
  deepVectorClear(mNTrackletsPerClusterSum);
  deepVectorClear(mNClustersPerROF);
  deepVectorClear(mTransitionPhiCuts);
  deepVectorClear(mTransitionMSAngles);
  deepVectorClear(mPositionResolution);
  deepVectorClear(mClusterSize);
  deepVectorClear(mPValphaX);
  deepVectorClear(mBogusClusters);
  deepVectorClear(mTrackletsIndexROF);
  deepVectorClear(mTrackletClusters);
  deepVectorClear(mLines);
  // if we use the external host allocator then the assumption is that we
  // don't clear the memory ourself
  if (!hasFrameworkAllocator()) {
    deepVectorClear(mClusters);
    deepVectorClear(mUsedClusters);
    deepVectorClear(mUnsortedClusters);
    deepVectorClear(mIndexTables);
    deepVectorClear(mTrackingFrameInfo);
    deepVectorClear(mROFramesClusters);
  }
  // only needed to clear if we have MC info
  if (hasMCinformation()) {
    deepVectorClear(mLinesLabels);
    deepVectorClear(mTrackletLabels);
    deepVectorClear(mCellLabels);
    deepVectorClear(mTracksLabel);
  }
  // mClusterLabels holds non-owning pointers into caller-supplied MC label
  // containers, not owned storage -- reset to nullptr, not deepVectorClear'd.
  mClusterLabels.fill(nullptr);
}

template class LegacyTrackerScratch<ITSNLayers>;
template class LegacyTrackerScratch<o2::mft::constants::mft::LayersNumber>;

} // namespace o2::itsmft::tracking

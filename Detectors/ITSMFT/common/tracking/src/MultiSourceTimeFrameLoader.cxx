// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <memory>
#include <vector>

#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/SurfaceMask.h"
#include "ITSMFTTracking/TimeFrame.h"

namespace o2::itsmft::tracking
{

LoadSourcesResult MultiSourceTimeFrameLoader::load(TimeFrame& frame, gsl::span<const ClusterSourceInput> sources,
                                                   SurfaceCatalogView catalog, const o2::InteractionRecord& origin)
{
  if (!frame.isConfigured()) {
    return {MultiSourceLoadError::FrameNotConfigured};
  }
  // Stage measurement data for every source before committing the event. The
  // source-level loader remains the decoder and timing-validation boundary.
  TimeFrame stagedMeasurements;
  const auto loadResult = loadSources(stagedMeasurements, catalog, sources, origin);
  if (!loadResult.ok()) {
    return loadResult;
  }

  const auto* binding = frame.getBinding(0);
  if (binding == nullptr) {
    return {MultiSourceLoadError::FrameNotConfigured};
  }
  SurfaceMask mappedSurfaces;
  for (const auto& source : sources) {
    for (const auto surface : source.layerToSurface) {
      if (!surface.isValid() || mappedSurfaces.has(surface) || !binding->getOwnedSurfaceIndex(surface)) {
        return {MultiSourceLoadError::InvalidLayerMapping, source.id};
      }
      mappedSurfaces.set(surface);
    }
  }
  if (mappedSurfaces != binding->getOwnedSurfaces()) {
    // Attribute an omitted surface when exactly one input stream owns its
    // detector family. Ambiguous same-detector splits remain unattributed.
    for (const auto surface : binding->getOrderedSurfaces()) {
      if (mappedSurfaces.has(surface)) {
        continue;
      }
      ClusterSourceId owner;
      for (const auto& source : sources) {
        if (static_cast<uint8_t>(source.detector) != catalog.getSurface(surface).detectorId) {
          continue;
        }
        if (owner.isValid()) {
          return {MultiSourceLoadError::InvalidLayerMapping};
        }
        owner = source.id;
      }
      return {MultiSourceLoadError::InvalidLayerMapping, owner};
    }
    return {MultiSourceLoadError::InvalidLayerMapping};
  }

  auto staged = std::make_unique<SurfaceTrackingScratch>();
  auto& live = frame.getWorkspace();
  if (live.hasFrameworkAllocator()) {
    staged->setFrameworkAllocator(live.getFrameworkAllocator());
  }
  staged->setMemoryPool(frame.getMemoryPool());
  staged->adoptPlan(live.getNOwnedSurfaces(), 0, 0);
  const auto backfillResult = staged->backfillNormalizedSources(stagedMeasurements, sources, binding->getOrderedSurfaces());
  if (!backfillResult.ok()) {
    return backfillResult;
  }

  if (!frame.commitLoadedEvent(std::move(stagedMeasurements), std::move(staged))) {
    return {MultiSourceLoadError::OtherMalformedInput};
  }
  return {};
}

namespace
{
bool isSupportedDetector(o2::detectors::DetID::ID det) noexcept
{
  return det == o2::detectors::DetID::ITS || det == o2::detectors::DetID::MFT;
}

bool covariance2DIsPositiveSemidefinite(float varianceFirst, float covariance,
                                        float varianceSecond) noexcept
{
  if (!std::isfinite(varianceFirst) || !std::isfinite(covariance) ||
      !std::isfinite(varianceSecond) || varianceFirst < 0.f || varianceSecond < 0.f) {
    return false;
  }
  const double diagonalProduct = static_cast<double>(varianceFirst) * varianceSecond;
  const double covarianceSquared = static_cast<double>(covariance) * covariance;
  const double tolerance = 16. * std::numeric_limits<float>::epsilon() *
                           std::max(diagonalProduct, covarianceSquared);
  return diagonalProduct - covarianceSquared >= -tolerance;
}

bool globalCovarianceIsPositiveSemidefinite(const GlobalCovariance3F& covariance) noexcept
{
  if (!covariance2DIsPositiveSemidefinite(covariance.xx, covariance.xy, covariance.yy) ||
      !covariance2DIsPositiveSemidefinite(covariance.xx, covariance.xz, covariance.zz) ||
      !covariance2DIsPositiveSemidefinite(covariance.yy, covariance.yz, covariance.zz)) {
    return false;
  }
  const double determinant =
    static_cast<double>(covariance.xx) * covariance.yy * covariance.zz +
    2. * static_cast<double>(covariance.xy) * covariance.xz * covariance.yz -
    static_cast<double>(covariance.xx) * covariance.yz * covariance.yz -
    static_cast<double>(covariance.yy) * covariance.xz * covariance.xz -
    static_cast<double>(covariance.zz) * covariance.xy * covariance.xy;
  const double scale = std::max({std::abs(static_cast<double>(covariance.xx) * covariance.yy * covariance.zz),
                                 std::abs(2. * static_cast<double>(covariance.xy) * covariance.xz * covariance.yz),
                                 std::abs(static_cast<double>(covariance.xx) * covariance.yz * covariance.yz),
                                 std::abs(static_cast<double>(covariance.yy) * covariance.xz * covariance.xz),
                                 std::abs(static_cast<double>(covariance.zz) * covariance.xy * covariance.xy)});
  return std::isfinite(determinant) &&
         determinant >= -32. * std::numeric_limits<float>::epsilon() * scale;
}

bool decodedMeasurementIsValid(const GlobalMeasurement& global,
                               const SurfaceMeasurement& local) noexcept
{
  return std::isfinite(global.position.x) && std::isfinite(global.position.y) &&
         std::isfinite(global.position.z) &&
         globalCovarianceIsPositiveSemidefinite(global.covariance) &&
         std::isfinite(local.frame.q) && std::isfinite(local.frame.u) &&
         std::isfinite(local.frame.v) && std::isfinite(local.frame.frameAngle) &&
         covariance2DIsPositiveSemidefinite(local.covariance.uu,
                                            local.covariance.uv,
                                            local.covariance.vv);
}

MultiSourceLoadError mapDecodeError(ClusterDecodeError error) noexcept
{
  switch (error) {
    case ClusterDecodeError::None:
      return MultiSourceLoadError::None;
    case ClusterDecodeError::MissingDictionary:
      return MultiSourceLoadError::MissingDictionary;
    case ClusterDecodeError::TruncatedExplicitPattern:
      return MultiSourceLoadError::TruncatedExplicitPattern;
    case ClusterDecodeError::MalformedExplicitPattern:
      return MultiSourceLoadError::MalformedExplicitPattern;
    case ClusterDecodeError::InvalidPatternId:
      return MultiSourceLoadError::InvalidPatternId;
    case ClusterDecodeError::InvalidSensor:
      return MultiSourceLoadError::InvalidSensor;
    case ClusterDecodeError::InvalidLayer:
      return MultiSourceLoadError::InvalidDecodedLayer;
    case ClusterDecodeError::InvalidLayerMapping:
      return MultiSourceLoadError::InvalidLayerMapping;
    case ClusterDecodeError::GeometryUnavailable:
      return MultiSourceLoadError::GeometryUnavailable;
    case ClusterDecodeError::OtherMalformedInput:
      return MultiSourceLoadError::OtherMalformedInput;
  }
  return MultiSourceLoadError::OtherMalformedInput;
}
} // namespace

LoadSourcesResult loadSources(TimeFrame& frame,
                              const SurfaceCatalogView& catalog,
                              gsl::span<const ClusterSourceInput> sources,
                              const o2::InteractionRecord& origin)
{
  const auto nSources = static_cast<uint32_t>(sources.size());

  std::vector<bool> seen(nSources, false);
  for (const auto& src : sources) {
    if (!src.id.isValid() || src.id.value() >= nSources) {
      return {MultiSourceLoadError::NonDenseSourceIds, src.id};
    }
    if (seen[src.id.value()]) {
      return {MultiSourceLoadError::DuplicateSourceId, src.id};
    }
    seen[src.id.value()] = true;
    if (!isSupportedDetector(src.detector)) {
      return {MultiSourceLoadError::UnsupportedDetector, src.id};
    }
    if (src.decoder == nullptr) {
      return {MultiSourceLoadError::MissingDecoder, src.id};
    }
    if (!src.clusters.empty() && src.dictionary == nullptr) {
      return {MultiSourceLoadError::MissingDictionary, src.id, 0, 0};
    }
  }

  std::vector<std::vector<GlobalMeasurement>> perSurfaceGlobal(catalog.nSurfaces);
  std::vector<std::vector<SurfaceMeasurement>> perSurface(catalog.nSurfaces);
  std::vector<SourceROFInfo> sourceInfo(nSources);
  std::vector<std::vector<ROFIntervalBC>> perSourceIntervals(nSources);
  std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*> labelSources(nSources, nullptr);

  for (const auto& src : sources) {
    const auto srcIdx = src.id.value();
    sourceInfo[srcIdx] = SourceROFInfo{src.id, static_cast<uint32_t>(src.rofs.size())};
    labelSources[srcIdx] = src.labels;

    int64_t expectedNext = 0;
    for (uint32_t r = 0; r < src.rofs.size(); ++r) {
      const auto& rof = src.rofs[r];
      const int64_t first = rof.getFirstEntry();
      const int64_t n = rof.getNEntries();
      if (n < 0 || first != expectedNext) {
        return {MultiSourceLoadError::InvalidROFRange, src.id, r};
      }
      expectedNext = first + n;
      if (expectedNext > static_cast<int64_t>(src.clusters.size())) {
        return {MultiSourceLoadError::InvalidROFRange, src.id, r};
      }
    }
    if (expectedNext != static_cast<int64_t>(src.clusters.size())) {
      return {MultiSourceLoadError::InvalidROFRange, src.id, static_cast<uint32_t>(src.rofs.size())};
    }

    auto& intervals = perSourceIntervals[srcIdx];
    intervals.reserve(src.rofs.size());
    for (uint32_t r = 0; r < src.rofs.size(); ++r) {
      const auto built = computeROFIntervalBC(src.rofs[r].getBCData(), origin, src.timing, r);
      if (!built.ok()) {
        return LoadSourcesResult{.error = MultiSourceLoadError::TimingError, .source = src.id, .rof = r, .timingDetail = built.error};
      }
      intervals.push_back(built.interval);
    }

    src.decoder->prepare();
    BoundedPatternCursor patterns{src.patterns};
    for (uint32_t r = 0; r < src.rofs.size(); ++r) {
      const auto& rof = src.rofs[r];
      const auto firstEntry = rof.getFirstEntry();
      const auto nEntries = rof.getNEntries();
      for (int32_t clusterId = firstEntry; clusterId < firstEntry + nEntries; ++clusterId) {
        const auto& cluster = src.clusters[clusterId];
        const auto externalIndex = static_cast<uint32_t>(clusterId);
        const auto decoded = src.decoder->decode(cluster, patterns, src.dictionary, src.layerToSurface,
                                                 src.id, externalIndex, r, src.applySysErrors);
        if (!decoded.ok()) {
          return {mapDecodeError(decoded.error), src.id, r, externalIndex};
        }
        if (!decoded.layerMapped || decoded.layer < 0 ||
            static_cast<size_t>(decoded.layer) >= src.layerToSurface.size()) {
          return {MultiSourceLoadError::InvalidLayerMapping, src.id, r, externalIndex};
        }
        const auto expectedSurface = src.layerToSurface[decoded.layer];
        if (!expectedSurface.isValid() || expectedSurface.value() >= catalog.nSurfaces) {
          return {MultiSourceLoadError::InvalidLayerMapping, src.id, r, externalIndex};
        }
        auto global = decoded.global;
        if (global.surface != expectedSurface ||
            global.cluster.source != src.id ||
            global.cluster.index != externalIndex ||
            global.sourceROF != r ||
            global.sensor.detector != static_cast<uint32_t>(src.detector)) {
          return {MultiSourceLoadError::InconsistentDecoderMetadata, src.id, r, externalIndex};
        }
        const auto& surfaceDescriptor = catalog.getSurface(expectedSurface);
        if (surfaceDescriptor.detectorId != static_cast<uint8_t>(src.detector)) {
          return {MultiSourceLoadError::DetectorSurfaceMismatch, src.id, r, externalIndex};
        }
        if (surfaceDescriptor.kind != decoded.kind) {
          return {MultiSourceLoadError::SurfaceKindMismatch, src.id, r, externalIndex};
        }
        if (!decodedMeasurementIsValid(global, decoded.measurement)) {
          return {MultiSourceLoadError::OtherMalformedInput, src.id, r, externalIndex};
        }
        global.radius = std::hypot(global.position.x, global.position.y);
        perSurfaceGlobal[expectedSurface.value()].push_back(global);
        perSurface[expectedSurface.value()].push_back(decoded.measurement);
      }
    }
    if (!patterns.empty()) {
      return {MultiSourceLoadError::TrailingPatternData, src.id,
              static_cast<uint32_t>(src.rofs.size()), static_cast<uint32_t>(src.clusters.size())};
    }
  }

  std::vector<uint32_t> sourceOffsets(nSources + 1, 0);
  for (uint32_t s = 0; s < nSources; ++s) {
    sourceOffsets[s + 1] = sourceOffsets[s] + static_cast<uint32_t>(perSourceIntervals[s].size());
  }
  std::vector<ROFIntervalBC> flatIntervals;
  flatIntervals.reserve(sourceOffsets.back());
  for (uint32_t s = 0; s < nSources; ++s) {
    flatIntervals.insert(flatIntervals.end(), perSourceIntervals[s].begin(), perSourceIntervals[s].end());
  }

  frame.assignLoadedMeasurements(std::move(perSurfaceGlobal), std::move(perSurface), std::move(sourceInfo),
                                 std::move(flatIntervals), std::move(sourceOffsets), std::move(labelSources));
  return {};
}

namespace
{
class NormalizedBackfillAllocatorMismatch final : public std::logic_error
{
 public:
  explicit NormalizedBackfillAllocatorMismatch(int layer)
    : std::logic_error("SurfaceTrackingScratch::loadNormalizedSource(): staged/live memory-resource mismatch on layer " + std::to_string(layer))
  {
  }
};
} // namespace

LoadSourcesResult SurfaceTrackingScratch::loadNormalizedSource(
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
  // The workspace is shared by both common-CA detector sources.
  constexpr ClusterSourceId kSourceId{0};
  if (detId != o2::detectors::DetID::MFT && detId != o2::detectors::DetID::ITS) {
    return {MultiSourceLoadError::UnsupportedDetector, kSourceId};
  }
  if (catalogView.surfaces == nullptr || catalogView.nSurfaces == 0) {
    return {MultiSourceLoadError::SurfaceCatalogNotConfigured, kSourceId};
  }
  if (orderedSurfaces.size() != mNOwnedSurfaces) {
    return {MultiSourceLoadError::InvalidLayerMapping, kSourceId};
  }
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
  const std::size_t nOwnedSurfaces = orderedSurfaces.size();

  TimeFrame staged;
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
  src.applySysErrors = applySysErrors;

  const auto result = loadSources(staged, catalogView, gsl::span<const ClusterSourceInput>(&src, 1), origin);
  if (!result.ok()) {
    return result;
  }

  const bool isMFT = (detId == o2::detectors::DetID::MFT);
  auto* pool = mMemoryPool.get();
  const auto nROFs = static_cast<size_t>(rofs.size());

  std::vector<bounded_vector<o2::its::Cluster>> stagedUnsortedClusters;
  std::vector<bounded_vector<o2::its::TrackingFrameInfo>> stagedTrackingFrameInfo;
  std::vector<bounded_vector<int>> stagedClusterExternalIndices;
  std::vector<bounded_vector<uint8_t>> stagedClusterSize;
  std::vector<bounded_vector<int>> stagedROFramesClusters;
  std::vector<const dataformats::MCTruthContainer<MCCompLabel>*> stagedClusterLabels(nOwnedSurfaces, nullptr);
  stagedUnsortedClusters.reserve(nOwnedSurfaces);
  stagedTrackingFrameInfo.reserve(nOwnedSurfaces);
  stagedClusterExternalIndices.reserve(nOwnedSurfaces);
  stagedClusterSize.reserve(nOwnedSurfaces);
  stagedROFramesClusters.reserve(nOwnedSurfaces);

  for (std::size_t layer = 0; layer < nOwnedSurfaces; ++layer) {
    const auto measurements = staged.getSurfaceMeasurements(layerToSurface[layer]);
    const auto globals = staged.getGlobalMeasurements(layerToSurface[layer]);

    auto* mr = getMaybeFrameworkHostResource();
    auto* unsortedClustersMr = mr != nullptr ? mr : mUnsortedClusters[layer].get_allocator().resource();
    auto* trackingFrameInfoMr = mr != nullptr ? mr : mTrackingFrameInfo[layer].get_allocator().resource();
    auto* clusterExternalIndicesMr = pool != nullptr ? pool : mClusterExternalIndices[layer].get_allocator().resource();
    auto* clusterSizeMr = pool != nullptr ? pool : mClusterSize[layer].get_allocator().resource();
    auto* rofFramesClustersMr = mr != nullptr ? mr : mROFramesClusters[layer].get_allocator().resource();

    stagedUnsortedClusters.emplace_back(std::pmr::polymorphic_allocator<o2::its::Cluster>{unsortedClustersMr});
    stagedTrackingFrameInfo.emplace_back(std::pmr::polymorphic_allocator<o2::its::TrackingFrameInfo>{trackingFrameInfoMr});
    stagedClusterExternalIndices.emplace_back(std::pmr::polymorphic_allocator<int>{clusterExternalIndicesMr});
    stagedClusterSize.emplace_back(measurements.size(), uint8_t{0}, std::pmr::polymorphic_allocator<uint8_t>{clusterSizeMr});
    stagedROFramesClusters.emplace_back(nROFs + 1, 0, std::pmr::polymorphic_allocator<int>{rofFramesClustersMr});

    size_t mi{0};
    for (std::size_t measurementIndex = 0; measurementIndex < measurements.size(); ++measurementIndex) {
      const auto& m = measurements[measurementIndex];
      const auto& global = globals[measurementIndex];
      o2::its::TrackingFrameInfo tfInfo;
      if (isMFT) {
        // Recreate the synthetic legacy MFT representation from
        // normalized global position and row/column covariance.
        tfInfo = o2::its::TrackingFrameInfo{
          global.position.x, global.position.y, global.position.z,
          global.position.x, 0.f,
          std::array<float, 2>{global.position.y, global.position.z},
          std::array<float, 3>{m.covariance.uu, m.covariance.uv, m.covariance.vv}};
      } else {
        // ITS compatibility representation, using the normalized measurement.
        tfInfo = o2::its::TrackingFrameInfo{
          global.position.x, global.position.y, global.position.z,
          m.frame.q, m.frame.frameAngle,
          std::array<float, 2>{m.frame.u, m.frame.v},
          std::array<float, 3>{m.covariance.uu, m.covariance.uv, m.covariance.vv}};
      }
      stagedTrackingFrameInfo[layer].push_back(tfInfo);
      stagedUnsortedClusters[layer].emplace_back(global.position.x, global.position.y, global.position.z, static_cast<int>(stagedUnsortedClusters[layer].size()));
      stagedClusterExternalIndices[layer].push_back(static_cast<int>(global.cluster.index));
      stagedClusterSize[layer][mi] = static_cast<uint8_t>(std::clamp(global.shape.nPixels, 0u, 255u));
      ++mi;
    }

    size_t mj{0};
    for (size_t r = 0; r < nROFs; ++r) {
      while (mj < globals.size() && globals[mj].sourceROF == static_cast<uint32_t>(r)) {
        ++mj;
      }
      stagedROFramesClusters[layer][r + 1] = static_cast<int>(mj);
    }

    stagedClusterLabels[layer] = labels;
  }

  const size_t nClustersLayer1 = stagedUnsortedClusters[1].size();
  auto* nTrackletsPerClusterMr = pool != nullptr ? pool : mNTrackletsPerCluster[0].get_allocator().resource();
  auto* nTrackletsPerClusterSumMr = pool != nullptr ? pool : mNTrackletsPerClusterSum[0].get_allocator().resource();
  std::array<bounded_vector<int>, 2> stagedNTrackletsPerCluster{
    bounded_vector<int>(nClustersLayer1, 0, std::pmr::polymorphic_allocator<int>{nTrackletsPerClusterMr}),
    bounded_vector<int>(nClustersLayer1, 0, std::pmr::polymorphic_allocator<int>{nTrackletsPerClusterMr})};
  std::array<bounded_vector<int>, 2> stagedNTrackletsPerClusterSum{
    bounded_vector<int>(nClustersLayer1 + 1, 0, std::pmr::polymorphic_allocator<int>{nTrackletsPerClusterSumMr}),
    bounded_vector<int>(nClustersLayer1 + 1, 0, std::pmr::polymorphic_allocator<int>{nTrackletsPerClusterSumMr})};
  std::vector<ClusterSourceId> stagedSourceBySurface(nOwnedSurfaces, kSourceId);

  for (std::size_t layer = 0; layer < nOwnedSurfaces; ++layer) {
    if (mUnsortedClusters[layer].get_allocator().resource() != stagedUnsortedClusters[layer].get_allocator().resource() ||
        mTrackingFrameInfo[layer].get_allocator().resource() != stagedTrackingFrameInfo[layer].get_allocator().resource() ||
        mClusterExternalIndices[layer].get_allocator().resource() != stagedClusterExternalIndices[layer].get_allocator().resource() ||
        mClusterSize[layer].get_allocator().resource() != stagedClusterSize[layer].get_allocator().resource() ||
        mROFramesClusters[layer].get_allocator().resource() != stagedROFramesClusters[layer].get_allocator().resource()) {
      throw NormalizedBackfillAllocatorMismatch(static_cast<int>(layer));
    }
  }
  for (int i = 0; i < 2; ++i) {
    if (mNTrackletsPerCluster[i].get_allocator().resource() != stagedNTrackletsPerCluster[i].get_allocator().resource() ||
        mNTrackletsPerClusterSum[i].get_allocator().resource() != stagedNTrackletsPerClusterSum[i].get_allocator().resource()) {
      throw NormalizedBackfillAllocatorMismatch(-1);
    }
  }

  static_assert(noexcept(std::declval<bounded_vector<o2::its::Cluster>&>().swap(std::declval<bounded_vector<o2::its::Cluster>&>())));
  static_assert(noexcept(std::declval<bounded_vector<int>&>().swap(std::declval<bounded_vector<int>&>())));

  // The view was constructed by the adapter for this staged event. The
  // frame commit resets the previous event, including its non-owning ROF
  // context, so reinstall this already-validated view only after the new
  // normalized frame is atomically live.
  const auto stagedROFViews = mROFViews;
  frame.commitMeasurements(std::move(staged));
  setROFViews(stagedROFViews);
  mSourceBySurface.swap(stagedSourceBySurface);
  for (std::size_t layer = 0; layer < nOwnedSurfaces; ++layer) {
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

LoadSourcesResult SurfaceTrackingScratch::backfillNormalizedSources(
  const TimeFrame& measurements,
  gsl::span<const ClusterSourceInput> sources,
  gsl::span<const SurfaceId> orderedSurfaces)
{
  if (orderedSurfaces.size() != mNOwnedSurfaces || sources.empty()) {
    return {MultiSourceLoadError::InvalidLayerMapping};
  }

  mROFViews = sources.front().rofViews;
  mROFViewsBySurface.assign(mNOwnedSurfaces, {});
  mROFLocalLayerBySurface.assign(mNOwnedSurfaces, 0);
  mSourceBySurface.assign(mNOwnedSurfaces, ClusterSourceId::invalid());
  SurfaceMask seenSurfaces;
  for (std::size_t layer = 0; layer < orderedSurfaces.size(); ++layer) {
    const auto surface = orderedSurfaces[layer];
    if (!surface.isValid() || seenSurfaces.has(surface)) {
      return {MultiSourceLoadError::InvalidLayerMapping};
    }
    seenSurfaces.set(surface);

    const ClusterSourceInput* owner = nullptr;
    uint16_t localLayer = 0;
    for (const auto& source : sources) {
      const auto it = std::find(source.layerToSurface.begin(), source.layerToSurface.end(), surface);
      if (it == source.layerToSurface.end()) {
        continue;
      }
      if (owner != nullptr) {
        return {MultiSourceLoadError::InvalidLayerMapping, source.id};
      }
      owner = &source;
      localLayer = static_cast<uint16_t>(std::distance(source.layerToSurface.begin(), it));
    }
    if (owner == nullptr) {
      return {MultiSourceLoadError::InvalidLayerMapping};
    }
    mSourceBySurface[layer] = owner->id;

    const auto localMeasurements = measurements.getSurfaceMeasurements(surface);
    const auto globals = measurements.getGlobalMeasurements(surface);
    auto& clusters = mUnsortedClusters[layer];
    auto& trackingInfo = mTrackingFrameInfo[layer];
    auto& externalIndices = mClusterExternalIndices[layer];
    auto& clusterSizes = mClusterSize[layer];
    auto& rofBoundaries = mROFramesClusters[layer];
    clusters.clear();
    trackingInfo.clear();
    externalIndices.clear();
    clusterSizes.assign(localMeasurements.size(), uint8_t{0});
    rofBoundaries.assign(owner->rofs.size() + 1, 0);
    clusters.reserve(localMeasurements.size());
    trackingInfo.reserve(localMeasurements.size());
    externalIndices.reserve(localMeasurements.size());

    const bool isMFT = owner->detector == o2::detectors::DetID::MFT;
    std::size_t measurementIndex = 0;
    for (std::size_t localIndex = 0; localIndex < localMeasurements.size(); ++localIndex) {
      const auto& measurement = localMeasurements[localIndex];
      const auto& global = globals[localIndex];
      if (global.surface != surface || !global.cluster.isValid() ||
          global.cluster.source != owner->id || global.sourceROF >= owner->rofs.size()) {
        return {MultiSourceLoadError::InconsistentDecoderMetadata, owner->id,
                global.sourceROF, global.cluster.index};
      }
      if (isMFT) {
        trackingInfo.emplace_back(
          global.position.x, global.position.y, global.position.z,
          global.position.x, 0.f,
          std::array<float, 2>{global.position.y, global.position.z},
          std::array<float, 3>{measurement.covariance.uu, measurement.covariance.uv, measurement.covariance.vv});
      } else {
        trackingInfo.emplace_back(
          global.position.x, global.position.y, global.position.z,
          measurement.frame.q, measurement.frame.frameAngle,
          std::array<float, 2>{measurement.frame.u, measurement.frame.v},
          std::array<float, 3>{measurement.covariance.uu, measurement.covariance.uv, measurement.covariance.vv});
      }
      clusters.emplace_back(global.position.x, global.position.y, global.position.z,
                            static_cast<int>(clusters.size()));
      externalIndices.push_back(static_cast<int>(global.cluster.index));
      clusterSizes[measurementIndex++] = static_cast<uint8_t>(std::clamp(global.shape.nPixels, 0u, 255u));
    }

    std::size_t cursor = 0;
    for (std::size_t rof = 0; rof < owner->rofs.size(); ++rof) {
      while (cursor < globals.size() && globals[cursor].sourceROF == rof) {
        ++cursor;
      }
      rofBoundaries[rof + 1] = static_cast<int>(cursor);
    }
    if (cursor != localMeasurements.size()) {
      return {MultiSourceLoadError::InconsistentDecoderMetadata, owner->id};
    }
    mClusterLabels[layer] = owner->labels;
    mROFViewsBySurface[layer] = owner->rofViews;
    mROFLocalLayerBySurface[layer] = localLayer;
  }

  const std::size_t pivotLayer = orderedSurfaces.size() > 1 ? 1 : 0;
  const std::size_t pivotClusters = orderedSurfaces.empty() ? 0 : mUnsortedClusters[pivotLayer].size();
  for (int i = 0; i < 2; ++i) {
    mNTrackletsPerCluster[i].assign(pivotClusters, 0);
    mNTrackletsPerClusterSum[i].assign(pivotClusters + 1, 0);
  }
  mUseUPC = false;
  return {};
}

namespace
{
std::string formatLoadSourcesResult(const char* label, const LoadSourcesResult& result)
{
  return std::format("{}: error={} source={} rof={} clusterIndex={} timingDetail={}",
                     label, static_cast<int>(result.error), result.source.value(),
                     result.rof, result.clusterIndex, static_cast<int>(result.timingDetail));
}
} // namespace

RecoverableLoadFailure::RecoverableLoadFailure(const LoadSourcesResult& result)
  : std::runtime_error(formatLoadSourcesResult("TimeFrame loading boundary: recoverable data failure", result)),
    mResult(result)
{
}

TimeFrameLoadException::TimeFrameLoadException(TimeFrameLoadFailureReason reason, std::string message)
  : std::runtime_error(std::move(message)), mReason(reason)
{
}

TimeFrameLoadException::TimeFrameLoadException(const LoadSourcesResult& result)
  : std::runtime_error(formatLoadSourcesResult("TimeFrame loading boundary: structural failure", result)),
    mReason(TimeFrameLoadFailureReason::LoadSourcesFailure),
    mLoadResult(result)
{
}

bool isRecoverableLoadError(MultiSourceLoadError error, TimingBuildError timingDetail) noexcept
{
  switch (error) {
    case MultiSourceLoadError::InvalidROFRange:
    case MultiSourceLoadError::TruncatedExplicitPattern:
    case MultiSourceLoadError::MalformedExplicitPattern:
    case MultiSourceLoadError::InvalidPatternId:
    case MultiSourceLoadError::InvalidSensor:
    case MultiSourceLoadError::InvalidDecodedLayer:
    case MultiSourceLoadError::OtherMalformedInput:
    case MultiSourceLoadError::TrailingPatternData:
      return true;
    case MultiSourceLoadError::TimingError:
      return timingDetail == TimingBuildError::Overflow;
    case MultiSourceLoadError::None:
    case MultiSourceLoadError::NonDenseSourceIds:
    case MultiSourceLoadError::DuplicateSourceId:
    case MultiSourceLoadError::UnsupportedDetector:
    case MultiSourceLoadError::MissingDecoder:
    case MultiSourceLoadError::InvalidLayerMapping:
    case MultiSourceLoadError::DetectorSurfaceMismatch:
    case MultiSourceLoadError::InconsistentDecoderMetadata:
    case MultiSourceLoadError::SurfaceKindMismatch:
    case MultiSourceLoadError::SurfaceCatalogNotConfigured:
    case MultiSourceLoadError::SurfaceCatalogStale:
    case MultiSourceLoadError::MissingDictionary:
    case MultiSourceLoadError::GeometryUnavailable:
    case MultiSourceLoadError::FrameNotConfigured:
      return false;
  }
  return false;
}

} // namespace o2::itsmft::tracking

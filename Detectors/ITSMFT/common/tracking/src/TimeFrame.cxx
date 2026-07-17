// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
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
/// \file TimeFrame.cxx
/// \brief
///

#include <numeric>

#include "Framework/Logger.h"
#include "ITSMFTTracking/DetectorTraits.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITStracking/MathUtils.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"

namespace
{
struct ClusterHelper {
  float rowCoord;
  float r;
  int bin;
  int ind;
};
} // namespace

namespace o2::itsmft::tracking
{

#ifndef GPUCA_GPUCODE
namespace
{
SurfaceMask positionalSurfaceMask(LayerMask layerMask, gsl::span<const SurfaceId> orderedSurfaces, uint32_t activeCount)
{
  SurfaceMask result;
  for (uint32_t position = 0; position < activeCount; ++position) {
    if (layerMask.has(position)) {
      result.set(orderedSurfaces[position]);
    }
  }
  return result;
}

DetectorSurfaceCatalogValidationError validateSurfaceCatalog(const DetectorSurfaceCatalogRequest& request,
                                                             gsl::span<const SurfaceDescriptor> catalog)
{
  if (request.detector < o2::detectors::DetID::First || request.detector > o2::detectors::DetID::Last) {
    return DetectorSurfaceCatalogValidationError::InvalidDetector;
  }
  if (!request.firstSurface.isValid()) {
    return DetectorSurfaceCatalogValidationError::InvalidFirstSurface;
  }
  if (request.detectorSurfaceCount == 0) {
    return DetectorSurfaceCatalogValidationError::EmptyDetector;
  }
  const uint64_t expectedSize = static_cast<uint64_t>(request.firstSurface.value()) + request.detectorSurfaceCount;
  if (expectedSize > MaxLayoutSurfaces) {
    return DetectorSurfaceCatalogValidationError::TooManySurfaces;
  }
  if (catalog.size() != expectedSize) {
    return DetectorSurfaceCatalogValidationError::SizeMismatch;
  }
  for (size_t globalIndex = 0; globalIndex < catalog.size(); ++globalIndex) {
    if (catalog[globalIndex].id != SurfaceId{static_cast<uint16_t>(globalIndex)}) {
      return DetectorSurfaceCatalogValidationError::NonDenseGlobalSurfaceIds;
    }
  }

  std::vector<bool> detectorSurfaceIndices(request.detectorSurfaceCount, false);
  for (size_t globalIndex = request.firstSurface.value(); globalIndex < expectedSize; ++globalIndex) {
    const auto& surface = catalog[globalIndex];
    if (surface.detectorId != static_cast<uint8_t>(request.detector)) {
      return DetectorSurfaceCatalogValidationError::DetectorMismatch;
    }
    if (surface.detectorSurfaceIndex == std::numeric_limits<uint16_t>::max()) {
      return DetectorSurfaceCatalogValidationError::MissingDetectorSurfaceIndex;
    }
    if (surface.detectorSurfaceIndex >= request.detectorSurfaceCount) {
      return DetectorSurfaceCatalogValidationError::DetectorSurfaceIndexOutOfRange;
    }
    if (detectorSurfaceIndices[surface.detectorSurfaceIndex]) {
      return DetectorSurfaceCatalogValidationError::DuplicateDetectorSurfaceIndex;
    }
    detectorSurfaceIndices[surface.detectorSurfaceIndex] = true;
  }
  if (std::find(detectorSurfaceIndices.begin(), detectorSurfaceIndices.end(), false) != detectorSurfaceIndices.end()) {
    return DetectorSurfaceCatalogValidationError::MissingDetectorSurfaceIndex;
  }
  return DetectorSurfaceCatalogValidationError::None;
}
} // namespace

template <int NLayers>
bool TimeFrame<NLayers>::detectorLayoutsCurrent() const noexcept
{
  return mDetectorLayouts.has_value() && mRequiredDetectorLayoutConfiguration.has_value() &&
         mDetectorLayouts->getConfigurationKey() == *mRequiredDetectorLayoutConfiguration &&
         mDetectorLayouts->getConfigurationKey().geometryEpoch == mRequiredDetectorGeometryEpoch;
}

template <int NLayers>
void TimeFrame<NLayers>::invalidateDetectorLayouts() noexcept
{
  const auto previousEpoch = mRequiredDetectorGeometryEpoch;
  mRequiredDetectorGeometryEpoch = nextDetectorGeometryEpoch(previousEpoch);
  if (previousEpoch == std::numeric_limits<DetectorGeometryEpoch>::max()) {
    // Epoch one may have existed before the counter wrapped. Dropping the old
    // owner is the only way to prevent that ancient catalog becoming current.
    mDetectorLayouts.reset();
  }
  if (mRequiredDetectorLayoutConfiguration) {
    mRequiredDetectorLayoutConfiguration->geometryEpoch = mRequiredDetectorGeometryEpoch;
  }
}

template <int NLayers>
DetectorLayoutSetBuildResult TimeFrame<NLayers>::ensureDetectorLayouts(const DetectorSurfaceCatalogProvider* provider,
                                                                       const DetectorSurfaceCatalogRequest& catalogRequest,
                                                                       gsl::span<const SurfaceId> orderedSurfaces,
                                                                       TransitionPolicyTag policyTag,
                                                                       gsl::span<const TrackingParameters> trackingParameters)
{
  DetectorLayoutConfigurationKey key;
  key.geometryEpoch = mRequiredDetectorGeometryEpoch;
  key.catalogRequest = catalogRequest;
  key.orderedSurfaces.assign(orderedSurfaces.begin(), orderedSurfaces.end());
  key.policyTag = policyTag;
  key.iterations.reserve(trackingParameters.size());
  for (const auto& parameters : trackingParameters) {
    if (parameters.NLayers < 0 || static_cast<size_t>(parameters.NLayers) > orderedSurfaces.size()) {
      const auto failedIteration = key.iterations.size();
      mRequiredDetectorLayoutConfiguration.reset();
      return {.error = DetectorLayoutSetBuildError::InvalidActiveCount,
              .failedIteration = failedIteration};
    }
    key.iterations.push_back(DetectorLayoutIterationConfiguration{
      static_cast<uint32_t>(parameters.NLayers), parameters.MaxHoles,
      parameters.HoleLayerMask, parameters.StartLayerMask});
  }

  mRequiredDetectorLayoutConfiguration = key;
  if (detectorLayoutsCurrent()) {
    return {};
  }
  if (provider == nullptr) {
    return {.error = DetectorLayoutSetBuildError::MissingProvider};
  }

  auto catalogResult = provider->buildCatalog(catalogRequest);
  if (!catalogResult.ok()) {
    return {.error = DetectorLayoutSetBuildError::CatalogProviderFailure,
            .catalogError = catalogResult.error};
  }
  const auto catalogValidationError = validateSurfaceCatalog(catalogRequest, catalogResult.catalog);
  if (catalogValidationError != DetectorSurfaceCatalogValidationError::None) {
    return {.error = DetectorLayoutSetBuildError::InvalidCatalog,
            .catalogValidationError = catalogValidationError};
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

    DetectorLayoutBuilder builder{catalogResult.catalog};
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

  DetectorLayoutSet stagedSet{std::move(key), std::move(catalogResult.catalog), std::move(staging)};
  static_assert(std::is_nothrow_move_constructible_v<DetectorLayoutSet>);
  mDetectorLayouts.emplace(std::move(stagedSet));
  return {.rebuilt = true};
}
#endif

using o2::its::clearResizeBoundedVector;
using o2::its::deepVectorClear;
namespace math_utils = o2::its::math_utils;
using o2::itsmft::IndexTableCoordType;
using o2::itsmft::IterationStep;
using o2::itsmft::TrackingParameters;

template <int NLayers>
void TimeFrame<NLayers>::addPrimaryVertex(const Vertex& vert)
{
  mPrimaryVertices.emplace_back(vert);
  if (!isBeamPositionOverridden) {
    const float w = vert.getNContributors();
    mBeamPos[0] = (mBeamPos[0] * mBeamPosWeight + vert.getX() * w) / (mBeamPosWeight + w);
    mBeamPos[1] = (mBeamPos[1] * mBeamPosWeight + vert.getY() * w) / (mBeamPosWeight + w);
    mBeamPosWeight += w;
  }
}

template <int NLayers>
void TimeFrame<NLayers>::loadROFrameData(gsl::span<const o2::itsmft::ROFRecord> rofs,
                                         gsl::span<const itsmft::CompClusterExt> clusters,
                                         gsl::span<const unsigned char>::iterator& pattIt,
                                         const itsmft::TopologyDictionary* dict,
                                         int layer,
                                         const dataformats::MCTruthContainer<MCCompLabel>* mcLabels,
                                         o2::detectors::DetID::ID detId)
{
  mDetId = detId;
  if (NLayers != constants::nLayersForDet(detId)) {
    LOGP(fatal, "TimeFrame<{}> is incompatible with detector {} (expected {} layers)",
         NLayers, static_cast<int>(detId), constants::nLayersForDet(detId));
  }
  ioutils::fillMatrixCache(detId);
  resetROFrameData(layer);
  prepareROFrameData(clusters, layer);

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

template <int NLayers>
LoadSourcesResult TimeFrame<NLayers>::loadNormalizedSource(
  const DetectorLayoutView& layout,
  gsl::span<const SurfaceId> layerToSurface,
  const ClusterDecoder& decoder,
  const o2::InteractionRecord& origin,
  const ROFTimingConfig& timing,
  gsl::span<const itsmft::CompClusterExt> clusters,
  gsl::span<const unsigned char> patterns,
  gsl::span<const o2::itsmft::ROFRecord> rofs,
  const itsmft::TopologyDictionary* dictionary,
  const dataformats::MCTruthContainer<MCCompLabel>* labels,
  o2::detectors::DetID::ID detId,
  bool applySysErrors)
{
  // Exactly one source is ever submitted here, and loadSources() requires
  // dense, zero-based source IDs, so ClusterSourceId{0} is the only value
  // that could ever succeed; it is fixed internally rather than accepted as
  // a parameter that would misleadingly suggest other IDs are supported.
  constexpr ClusterSourceId kSourceId{0};

  if (NLayers != constants::nLayersForDet(detId)) {
    return {MultiSourceLoadError::UnsupportedDetector, kSourceId};
  }
  if (layerToSurface.size() != static_cast<size_t>(NLayers)) {
    return {MultiSourceLoadError::InvalidLayerMapping, kSourceId};
  }

  // Stage into a scratch owner: loadSources() itself never mutates its
  // `frame` argument on failure, but staging separately also protects the
  // *existing* mNormalizedFrame (and, by construction below, every legacy
  // compatibility structure) from a failed reload.
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

  const auto result = loadSources(staged, layout, gsl::span<const ClusterSourceInput>(&src, 1), origin);
  if (!result.ok()) {
    // Nothing below has run: mNormalizedFrame and every legacy compatibility
    // structure are exactly as they were before this call.
    return result;
  }

  // From here on, loadSources() has already succeeded and every remaining
  // step is plain data transformation over its output; no further
  // *validation* failure is possible, so mNormalizedFrame and the legacy
  // compatibility structures below are always updated together whenever this
  // function returns an ok() result. This is not a strong exception-safety
  // guarantee: an allocation failure inside the backfill below could still
  // throw after mNormalizedFrame has already been committed, leaving the two
  // representations inconsistent; that guarantee is not required here.
  mNormalizedFrame = std::move(staged);
  mDetId = detId;

  const bool isMFT = (detId == o2::detectors::DetID::MFT);
  auto* mr = getMaybeFrameworkHostResource();
  const auto nROFs = static_cast<size_t>(rofs.size());

  for (int layer = 0; layer < NLayers; ++layer) {
    const auto measurements = mNormalizedFrame.getSurfaceMeasurements(layerToSurface[layer]);

    deepVectorClear(mUnsortedClusters[layer], mr);
    deepVectorClear(mTrackingFrameInfo[layer], mr);
    deepVectorClear(mClusterExternalIndices[layer], mMemoryPool.get());
    clearResizeBoundedVector(mClusterSize[layer], measurements.size(), mMemoryPool.get());
    clearResizeBoundedVector(mROFramesClusters[layer], nROFs + 1, mr, 0);

    size_t mi{0};
    for (const auto& m : measurements) {
      TrackingFrameInfo tfInfo;
      if (isMFT) {
        // Recreate the established synthetic legacy MFT representation
        // (TrackingFrameInfoAdapters.h::makeTrackingFrameInfo<MFT>) from the
        // normalized global position and row/column covariance. This is
        // deliberately not m.frame, which for a Disk-kind SurfaceMeasurement
        // holds the (z, x, y) disk-frame projection, not the legacy
        // synthetic layout existing production code consumes.
        tfInfo = TrackingFrameInfo{
          m.global.x, m.global.y, m.global.z,
          m.global.x, 0.f,
          std::array<float, 2>{m.global.y, m.global.z},
          std::array<float, 3>{m.covariance.uu, m.covariance.uv, m.covariance.vv}};
      } else {
        tfInfo = TrackingFrameInfo{
          m.global.x, m.global.y, m.global.z,
          m.frame.q, m.frame.frameAngle,
          std::array<float, 2>{m.frame.u, m.frame.v},
          std::array<float, 3>{m.covariance.uu, m.covariance.uv, m.covariance.vv}};
      }
      addTrackingFrameInfoToLayer(layer, tfInfo);
      addClusterToLayer(layer, m.global.x, m.global.y, m.global.z, mUnsortedClusters[layer].size());
      addClusterExternalIndexToLayer(layer, static_cast<int>(m.cluster.index));
      mClusterSize[layer][mi] = static_cast<uint8_t>(std::clamp(m.shape.nPixels, 0u, 255u));
      ++mi;
    }

    size_t mj{0};
    for (size_t r = 0; r < nROFs; ++r) {
      while (mj < measurements.size() && measurements[mj].sourceROF == static_cast<uint32_t>(r)) {
        ++mj;
      }
      mROFramesClusters[layer][r + 1] = static_cast<int>(mj);
    }

    mClusterLabels[layer] = labels;
  }

  return result;
}

template <int NLayers>
void TimeFrame<NLayers>::resetROFrameData(int layer)
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
void TimeFrame<NLayers>::prepareROFrameData(gsl::span<const itsmft::CompClusterExt> clusters, int layer)
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
      ++clusterCountPerLayer[ioutils::getClusterLayer(mDetId, cls)];
    }
    for (int iLayer{0}; iLayer < NLayers; ++iLayer) {
      mUnsortedClusters[iLayer].reserve(clusterCountPerLayer[iLayer]);
      mTrackingFrameInfo[iLayer].reserve(clusterCountPerLayer[iLayer]);
      mClusterExternalIndices[iLayer].reserve(clusterCountPerLayer[iLayer]);
    }
  }
}

template <int NLayers>
void TimeFrame<NLayers>::prepareClusters(const TrackingParameters& trkParam, const int maxLayers)
{
  const int numBins{trkParam.RowBins * trkParam.ColBins};
  const int stride{numBins + 1};
  bounded_vector<ClusterHelper> cHelper(mMemoryPool.get());
  bounded_vector<int> clsPerBin(numBins, 0, mMemoryPool.get());
  bounded_vector<int> lutPerBin(numBins, 0, mMemoryPool.get());
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
        ClusterHelper& h = cHelper[iCluster];

        const float x = c.xCoordinate - (useXYBinning ? 0.f : mBeamPos[0]);
        const float y = c.yCoordinate - (useXYBinning ? 0.f : mBeamPos[1]);
        const float z = c.zCoordinate;

        const float rowCoord = useXYBinning ? c.yCoordinate : math_utils::computePhi(x, y);
        const float colCoord = useXYBinning ? c.xCoordinate : z;
        int colBin{mIndexTableUtils.getColBinIndex(iLayer, colCoord)};
        if (colBin < 0 || colBin >= trkParam.ColBins) {
          colBin = std::clamp(colBin, 0, trkParam.ColBins - 1);
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
        c.phi = useXYBinning ? math_utils::computePhi(c.xCoordinate, c.yCoordinate) : h.rowCoord;
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
void TimeFrame<NLayers>::initVertexingTopology(const TrackingParameters& trkParam)
{
  mVertexingTopology.init(3, trkParam.MaxHoles, LayerMask{trkParam.HoleLayerMask});
}

template <int NLayers>
void TimeFrame<NLayers>::initDefaultTrackingTopology(const TrackingParameters& trkParam, const int maxLayers)
{
  mDefaultTrackingTopology.init(maxLayers, trkParam.MaxHoles, LayerMask{trkParam.HoleLayerMask});
}

template <int NLayers>
void TimeFrame<NLayers>::initTrackerTopologies(gsl::span<const TrackingParameters> trkParams, const int maxLayers)
{
  mTrackerTopologies.resize(trkParams.size());
  for (size_t iteration = 0; iteration < trkParams.size(); ++iteration) {
    const int iterationMaxLayers = std::min(maxLayers, trkParams[iteration].NLayers);
    mTrackerTopologies[iteration].init(iterationMaxLayers, trkParams[iteration].MaxHoles, LayerMask{trkParams[iteration].HoleLayerMask});
  }
}

template <int NLayers>
void TimeFrame<NLayers>::initialise(const TrackingParameters& trkParam, const int maxLayers, const int iteration)
{
  mTrackingTopologyView = iteration != constants::UnusedIndex ? mTrackerTopologies[iteration].getView() : (maxLayers == 3 ? mVertexingTopology.getView() : mDefaultTrackingTopology.getView());

  if (trkParam.PassFlags[IterationStep::FirstPass]) {
    deepVectorClear(mTracks);
    deepVectorClear(mTracksLabel);
    deepVectorClear(mLines);
    deepVectorClear(mLinesLabels);
    if (trkParam.PassFlags[IterationStep::ResetVertices]) {
      deepVectorClear(mPrimaryVertices);
      deepVectorClear(mPrimaryVerticesLabels);
    }
    clearResizeBoundedVector(mLinesLabels, getNrof(1), mMemoryPool.get());
    DetectorTraits<NLayers>::configureIndexTableUtils(mIndexTableUtils, trkParam);
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

    for (int iLayer{0}; iLayer < NLayers; ++iLayer) {
      clearResizeBoundedVector(mIndexTables[iLayer], getNrof(iLayer) * ((trkParam.ColBins * trkParam.RowBins) + 1), getMaybeFrameworkHostResource());
    }
    for (int iLayer{0}; iLayer < trkParam.NLayers; ++iLayer) {
      if (trkParam.SystError2Row[iLayer] > 0.f || trkParam.SystError2Col[iLayer] > 0.f) {
        for (auto& tfInfo : mTrackingFrameInfo[iLayer]) {
          /// Account for alignment systematics in the cluster covariance matrix
          tfInfo.covarianceTrackingFrame[0] += trkParam.SystError2Row[iLayer];
          tfInfo.covarianceTrackingFrame[2] += trkParam.SystError2Col[iLayer];
        }
      }
    }

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
    prepareClusters(trkParam, maxLayers);
  }
  mTotalTracklets = {0, 0};
  if (maxLayers < trkParam.NLayers) { // Vertexer only, but in both iterations
    for (size_t iLayer{0}; iLayer < maxLayers; ++iLayer) {
      deepVectorClear(mUsedClusters[iLayer]);
      clearResizeBoundedVector(mUsedClusters[iLayer], mUnsortedClusters[iLayer].size(), mMemoryPool.get());
    }
  }

  // estimate MS per layer
  std::array<float, NLayers> msAngles{};
  for (unsigned int iLayer{0}; iLayer < NLayers; ++iLayer) {
    if constexpr (NLayers == o2::mft::constants::mft::LayersNumber) {
      msAngles[iLayer] = detail::mftLayerMSAngle(iLayer, trkParam);
    } else {
      msAngles[iLayer] = math_utils::MSangle(0.14f, trkParam.TrackletMinPt, trkParam.LayerxX0[iLayer]);
    }
    mPositionResolution[iLayer] = o2::gpu::CAMath::Sqrt((0.5f * (trkParam.SystError2Col[iLayer] + trkParam.SystError2Row[iLayer])) + (trkParam.LayerResolution[iLayer] * trkParam.LayerResolution[iLayer]));
  }

  // for each transition calculate the phi-cuts + integrated MS
  if constexpr (NLayers == o2::mft::constants::mft::LayersNumber) {
    float oneOverR{0.001f * 0.3f * std::abs(mBz) / trkParam.TrackletMinPt};
    for (int transitionId{0}; transitionId < (int)mTracklets.size(); ++transitionId) {
      const auto& transition = mTrackingTopologyView.getTransition(transitionId);
      float ms2 = 0.f;
      for (int layer = transition.fromLayer; layer < transition.toLayer; ++layer) {
        ms2 += math_utils::Sq(msAngles[layer]);
      }
      mTransitionMSAngles[transitionId] = o2::gpu::CAMath::Sqrt(ms2);
      const float& r1 = trkParam.LayerRadii[transition.fromLayer];
      const float& r2 = trkParam.LayerRadii[transition.toLayer];
      oneOverR = (0.5f * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
      const float res1 = o2::gpu::CAMath::Hypot(trkParam.PVres, mPositionResolution[transition.fromLayer]);
      const float res2 = o2::gpu::CAMath::Hypot(trkParam.PVres, mPositionResolution[transition.toLayer]);
      const float cosTheta1half = o2::gpu::CAMath::Sqrt(1.f - math_utils::Sq(0.5f * r1 * oneOverR));
      const float cosTheta2half = o2::gpu::CAMath::Sqrt(1.f - math_utils::Sq(0.5f * r2 * oneOverR));
      const float x = (r2 * cosTheta1half) - (r1 * cosTheta2half);
      const float delta = o2::gpu::CAMath::Sqrt(1.f / (1.f - 0.25f * math_utils::Sq(x * oneOverR)) *
                                                (math_utils::Sq((0.25f * r1 * r2 * math_utils::Sq(oneOverR) / cosTheta2half) + cosTheta1half) * math_utils::Sq(res1) +
                                                 math_utils::Sq((0.25f * r1 * r2 * math_utils::Sq(oneOverR) / cosTheta1half) + cosTheta2half) * math_utils::Sq(res2)));
      // For MFT: bending angle (rad) at TrackletMinPt, used to widen x-y tracklet windows.
      mTransitionPhiCuts[transitionId] = o2::gpu::CAMath::Min(o2::gpu::CAMath::ASin(0.5f * x * oneOverR) + 2.f * mTransitionMSAngles[transitionId] + delta, o2::constants::math::PI * 0.5f);

      deepVectorClear(mTracklets[transitionId]);
      deepVectorClear(mTrackletLabels[transitionId]);
      deepVectorClear(mTrackletsLookupTable[transitionId]);
      mTrackletsLookupTable[transitionId].resize(mClusters[transition.fromLayer].size() + 1, 0);
    }
  } else {
    float oneOverR{0.001f * 0.3f * std::abs(mBz) / trkParam.TrackletMinPt};
    for (int transitionId{0}; transitionId < (int)mTracklets.size(); ++transitionId) {
      const auto& transition = mTrackingTopologyView.getTransition(transitionId);
      float ms2 = 0.f;
      for (int layer = transition.fromLayer; layer < transition.toLayer; ++layer) {
        ms2 += math_utils::Sq(msAngles[layer]);
      }
      mTransitionMSAngles[transitionId] = o2::gpu::CAMath::Sqrt(ms2);
      const float& r1 = trkParam.LayerRadii[transition.fromLayer];
      const float& r2 = trkParam.LayerRadii[transition.toLayer];
      oneOverR = (0.5 * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
      const float res1 = o2::gpu::CAMath::Hypot(trkParam.PVres, mPositionResolution[transition.fromLayer]);
      const float res2 = o2::gpu::CAMath::Hypot(trkParam.PVres, mPositionResolution[transition.toLayer]);
      const float cosTheta1half = o2::gpu::CAMath::Sqrt(1.f - math_utils::Sq(0.5f * r1 * oneOverR));
      const float cosTheta2half = o2::gpu::CAMath::Sqrt(1.f - math_utils::Sq(0.5f * r2 * oneOverR));
      float x = (r2 * cosTheta1half) - (r1 * cosTheta2half);
      float delta = o2::gpu::CAMath::Sqrt(1.f / (1.f - 0.25f * math_utils::Sq(x * oneOverR)) * (math_utils::Sq((0.25f * r1 * r2 * math_utils::Sq(oneOverR) / cosTheta2half) + cosTheta1half) * math_utils::Sq(res1) + math_utils::Sq((0.25f * r1 * r2 * math_utils::Sq(oneOverR) / cosTheta1half) + cosTheta2half) * math_utils::Sq(res2)));
      /// the expression std::asin(0.5f * x * oneOverR) is equivalent to std::aCos(0.5f * r1 * oneOverR) - std::acos(0.5 * r2 * oneOverR)
      mTransitionPhiCuts[transitionId] = o2::gpu::CAMath::Min(o2::gpu::CAMath::ASin(0.5f * x * oneOverR) + 2.f * mTransitionMSAngles[transitionId] + delta, o2::constants::math::PI * 0.5f);

      deepVectorClear(mTracklets[transitionId]);
      deepVectorClear(mTrackletLabels[transitionId]);
      deepVectorClear(mTrackletsLookupTable[transitionId]);
      mTrackletsLookupTable[transitionId].resize(mClusters[transition.fromLayer].size() + 1, 0);
    }
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
unsigned long TimeFrame<NLayers>::getArtefactsMemory() const
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
void TimeFrame<NLayers>::printArtefactsMemory() const
{
  LOGP(info, "TimeFrame: Artefacts occupy {:.2f} MB", getArtefactsMemory() / constants::MB);
}

template <int NLayers>
void TimeFrame<NLayers>::computeTrackletsPerROFScans()
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
void TimeFrame<NLayers>::setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
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
  initVector(mPrimaryVertices);
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
  initVector(mPrimaryVerticesLabels);
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
void TimeFrame<NLayers>::setFrameworkAllocator(ExternalAllocator* ext)
{
  mExternalAllocator = ext;
  mExtMemoryPool = std::make_shared<BoundedMemoryResource>(mExternalAllocator);
}

template <int NLayers>
void TimeFrame<NLayers>::wipe()
{
  deepVectorClear(mTracks);
  deepVectorClear(mTracklets);
  deepVectorClear(mCells);
  deepVectorClear(mCellsNeighbours);
  deepVectorClear(mCellsNeighboursTopology);
  deepVectorClear(mCellsLookupTable);
  deepVectorClear(mPrimaryVertices);
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
    deepVectorClear(mPrimaryVerticesLabels);
    deepVectorClear(mTrackletLabels);
    deepVectorClear(mCellLabels);
    deepVectorClear(mTracksLabel);
  }
}

template class TimeFrame<ITSNLayers>;
template class TimeFrame<o2::mft::constants::mft::LayersNumber>;

} // namespace o2::itsmft::tracking

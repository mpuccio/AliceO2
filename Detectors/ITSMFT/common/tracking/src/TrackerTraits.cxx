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
/// \file TrackerTraits.cxx
/// \brief
///

#include <algorithm>
#include <array>
#include <iterator>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/enumerable_thread_specific.h>

#include "Framework/Logger.h"
#include "GPUCommonMath.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITSMFTTracking/Cell.h"
#include "ITStracking/Constants.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IndexTableConfiguration.h"
#include "ITSMFTTracking/Propagator.h"
#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/detail/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/SurfaceMask.h"
#include "ITSMFTTracking/TripletFitting.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "ITSMFTTracking/detail/CandidateFinding.h"
#include "ITSMFTTracking/detail/DirectionCompatibility.h"
#include "ITStracking/Tracklet.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::itsmft::tracking
{

namespace math_utils = o2::its::math_utils;
using o2::its::deepVectorClear;
using o2::its::TimeEstBC;

struct PassMode {
  using OnePass = std::integral_constant<int, 0>;
  using TwoPassCount = std::integral_constant<int, 1>;
  using TwoPassInsert = std::integral_constant<int, 2>;
};

namespace
{
constexpr uint8_t kCompatibilityAbsCharge = 1;
const o2::track::PID kCompatibilityPID = o2::track::PID::Pion;

constexpr std::size_t kindIndex(SurfaceKind kind) noexcept
{
  return kind == SurfaceKind::Cylinder ? 0u : 1u;
}

std::optional<uint32_t> appendGenericTrack(TimeFrame& frame,
                                           const TrackingCandidate& candidate,
                                           gsl::span<const gsl::span<const GlobalMeasurement>> layerMeasurements)
{
  GenericTrack track = candidate.track;
  track.hitSurfaces = {};
  std::vector<TrackClusterReference> resolvedReferences;
  resolvedReferences.reserve(layerMeasurements.size());
  for (std::size_t position = 0; position < layerMeasurements.size(); ++position) {
    const int localIndex = candidate.getClusterIndex(static_cast<int>(position));
    if (localIndex == o2::its::constants::UnusedIndex) {
      continue;
    }
    if (localIndex < 0 || static_cast<std::size_t>(localIndex) >= layerMeasurements[position].size()) {
      return std::nullopt;
    }
    const auto& measurement = layerMeasurements[position][localIndex];
    const TrackClusterReference reference{measurement.surface, SurfaceMeasurementIndex{static_cast<uint32_t>(localIndex)}};
    const auto* resolved = frame.getGlobalMeasurement(reference.surface, reference.index);
    if (!reference.surface.isValid() || measurement.surface != reference.surface || !measurement.cluster.isValid() ||
        resolved == nullptr || resolved != &measurement || resolved->surface != reference.surface || !resolved->cluster.isValid()) {
      return std::nullopt;
    }
    resolvedReferences.push_back(reference);
    track.hitSurfaces.set(reference.surface);
  }
  if (!track.innerState.hasRecognizedKind() || !track.outerState.hasRecognizedKind() ||
      !track.timestamp.isValid() || resolvedReferences.empty()) {
    return std::nullopt;
  }

  auto& tracks = frame.getGenericTracks();
  auto& references = frame.getTrackClusterIndices();
  const auto oldTrackSize = tracks.size();
  const auto oldReferenceSize = references.size();
  if (oldTrackSize > std::numeric_limits<uint32_t>::max() || oldReferenceSize > std::numeric_limits<uint32_t>::max() ||
      resolvedReferences.size() > std::numeric_limits<uint32_t>::max() - oldReferenceSize) {
    return std::nullopt;
  }

  try {
    references.reserve(oldReferenceSize + resolvedReferences.size());
    tracks.reserve(oldTrackSize + 1);
    for (const auto& reference : resolvedReferences) {
      references.push_back(reference);
    }
    track.firstClusterRef = static_cast<uint32_t>(oldReferenceSize);
    track.clusterRefEnd = static_cast<uint32_t>(references.size());
    tracks.push_back(track);
  } catch (...) {
    references.resize(oldReferenceSize);
    tracks.resize(oldTrackSize);
    throw;
  }
  return static_cast<uint32_t>(oldTrackSize);
}

// A static diamond vertex represents all primary vertices and has no event
// timestamp. Derive its envelope from the tested ROF's configured bounds;
// TimeEstBC cannot represent a full TimeFrame. The resulting timestamp is
// compatible by construction with that ROF.
template <typename ROFOverlapView>
Vertex diamondVertexForROF(const Vertex& base, const ROFOverlapView& rofOverlapView, int layer, int rofId)
{
  Vertex v = base;
  v.setTimeStamp(rofOverlapView.getLayer(layer).getROFTimeBounds(rofId, true));
  return v;
}

// Convert ROOT-visible parameters to the device-portable record once per iteration.
} // namespace

void TrackerTraits::runTraversal(TraversalWorkspaceView view, SeedRefitFunction refitFunction)
{
  if (view.iteration < 0 || static_cast<std::size_t>(view.iteration) >= view.parameters.size()) {
    throw TraversalException{view.iteration, TraversalFailureReason::IterationOutOfRange};
  }
  if (!view.workspace.valid) {
    throw TraversalException{view.iteration, TraversalFailureReason::InvalidTraversalSchedule};
  }
  int maxNvertices{-1};
  if (view.parameters[0].PerPrimaryVertexProcessing) {
    maxNvertices = view.scratch.getMaxVerticesPerROF();
  }
  int iVertex = std::min(maxNvertices, 0);
  do {
    computeLayerTracklets(view, view.iteration, iVertex);
    computeLayerCells(view, view.iteration);
    findCellsNeighbours(view, view.iteration);
    findRoads(view, view.iteration, refitFunction);
  } while (++iVertex < maxNvertices);
}

int requireScratchLinkSlot(const TraversalWorkspaceView& context, int iteration, LinkId id)
{
  const auto slot = context.binding.getScratchLinkSlot(id);
  if (!slot) {
    throw TraversalException{iteration, TraversalFailureReason::TraversalBindingMismatch};
  }
  return static_cast<int>(*slot);
}

int requireScratchCellSlot(const TraversalWorkspaceView& context, int iteration, CellTopologyId id)
{
  const auto slot = context.binding.getScratchCellSlot(id);
  if (!slot) {
    throw TraversalException{iteration, TraversalFailureReason::TraversalBindingMismatch};
  }
  return static_cast<int>(*slot);
}

int requireSurfacePosition(const TraversalWorkspaceView& context, int iteration, SurfaceId id)
{
  if (!id.isValid()) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  const auto position = context.binding.getOwnedSurfaceIndex(id);
  if (!position || *position >= context.scratch.getNOwnedSurfaces()) {
    throw TraversalException{iteration, TraversalFailureReason::TraversalBindingMismatch};
  }
  return static_cast<int>(*position);
}

void TrackerTraits::computeLayerTracklets(TraversalWorkspaceView& context, const int iteration, int iVertex)
{
  auto& scratch = context.scratch;
  auto& workspace = context.workspace;
  const auto scratchLinkCount = scratch.getTracklets().size();
  for (size_t linkId = 0; linkId < scratchLinkCount; ++linkId) {
    scratch.getTracklets()[linkId].clear();
    scratch.getTrackletsLabel(linkId).clear();
    std::fill(scratch.getTrackletsLookupTable()[linkId].begin(), scratch.getTrackletsLookupTable()[linkId].end(), 0);
  }

  if (!workspace.valid) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidTraversalSchedule};
  }

  computeLayerTrackletsImpl(context, iteration, iVertex, context.binding.getGlobalLinks());
}

void TrackerTraits::computeLayerTrackletsImpl(
  TraversalWorkspaceView& context,
  const int iteration,
  const int iVertex,
  gsl::span<const LinkId> linkIds)
{
  auto* mScratch = &context.scratch;
  const auto& mMemoryPool = mScratch->getMemoryPool();
  auto* mFrame = &context.frame;
  const auto* mBinding = &context.binding;
  const auto mTrkParams = context.parameters;
  const auto mBz = context.bz;
  const auto& mTraversalGraph = context.graph;
  auto& mKernelParameters = context.workspace.kernelParameters;
  const auto& mDiskLayerReferenceZ = context.workspace.diskLayerReferenceZView;
  const auto& mLayerGlobalMeasurements = context.workspace.layerGlobalMeasurements;
  const auto& topology = mTraversalGraph;
  const Vertex diamondVert(mTrkParams[iteration].Diamond, mTrkParams[iteration].DiamondCov, 1, 1.f);

  mTaskArena->execute([&] {
    auto resolveLinkLayers = [&](int linkId) -> std::pair<int, int> {
      const auto& link = topology.getLink(LinkId{static_cast<uint16_t>(linkId)});
      return {requireSurfacePosition(context, iteration, link.from),
              requireSurfacePosition(context, iteration, link.to)};
    };

    auto forTracklets = [&](auto Mode, int linkId, int fromLayer, int toLayer, SurfaceKind kind,
                            const TrackletProjectionCache& linkCache, int pivotROF, int base, int& offset) -> int {
      if (!mScratch->isROFEnabled(fromLayer, pivotROF)) {
        return 0;
      }
      // Derive a diamond vertex for this pivot ROF; each invocation owns its
      // stack frame, so this is safe inside the parallel dispatch.
      Vertex diamondForROF{};
      gsl::span<const Vertex> primaryVertices;
      if (mTrkParams[iteration].UseDiamond) {
        diamondForROF = diamondVertexForROF(diamondVert, mScratch->getROFViews(fromLayer).overlap,
                                            mScratch->getROFLocalLayer(fromLayer), pivotROF);
        primaryVertices = gsl::span<const Vertex>(&diamondForROF, 1);
      } else {
        primaryVertices = mScratch->getPrimaryVertices(*mFrame, fromLayer, pivotROF);
      }
      if (primaryVertices.empty()) {
        return 0;
      }
      const int startVtx = iVertex >= 0 ? iVertex : 0;
      const int endVtx = iVertex >= 0 ? o2::gpu::CAMath::Min(iVertex + 1, int(primaryVertices.size())) : int(primaryVertices.size());
      if (endVtx <= startVtx || (iVertex + 1) > primaryVertices.size()) {
        return 0;
      }

      const auto& rofOverlap = mScratch->getROFOverlap(fromLayer, toLayer, pivotROF);
      if (!rofOverlap.getEntries()) {
        return 0;
      }

      int localCount = 0;
      auto& tracklets = mScratch->getTracklets()[linkId];
      auto layer0 = mScratch->getClustersOnLayer(pivotROF, fromLayer);
      if (layer0.empty()) {
        return 0;
      }

      for (int iCluster = 0; iCluster < int(layer0.size()); ++iCluster) {
        const o2::its::Cluster& currentCluster = layer0[iCluster];
        const int currentSortedIndex = mScratch->getSortedIndex(pivotROF, fromLayer, iCluster);
        if (mScratch->isClusterUsed(fromLayer, currentCluster.clusterId)) {
          continue;
        }
        const auto& sourceMeasurement = mLayerGlobalMeasurements[fromLayer][currentCluster.clusterId];

        for (int iV = startVtx; iV < endVtx; ++iV) {
          const auto& pv = primaryVertices[iV];
          if (!mScratch->isVertexCompatible(fromLayer, pivotROF, pv)) {
            continue;
          }
          if (pv.isFlagSet(Vertex::Flags::UPCMode) != mTrkParams[iteration].PassFlags[IterationStep::SelectUPCVertices]) {
            continue;
          }
          TrackletSearchWindow searchWindow{};
          const bool projected = projectTrackletSearchWindow(sourceMeasurement, currentCluster, pv,
                                                             kind, linkCache, mBz,
                                                             mScratch->getIndexTableUtils(toLayer),
                                                             mKernelParameters, searchWindow);
          if (!projected) {
            continue;
          }
          const auto bins = searchWindow.bins;
          const auto& indexTableUtils = mScratch->getIndexTableUtils(toLayer);
          int rowBinsNum = bins.w - bins.y + 1;
          if (indexTableUtils.getCoordType() == IndexTableCoordType::PhiZ && rowBinsNum < 0) {
            rowBinsNum += indexTableUtils.getNrowBins();
          }
          rowBinsNum = std::max(0, rowBinsNum);

          for (int targetROF = rofOverlap.getFirstEntry(); targetROF < rofOverlap.getEntriesBound(); ++targetROF) {
            if (!mScratch->isROFEnabled(toLayer, targetROF)) {
              continue;
            }
            auto layer1 = mScratch->getClustersOnLayer(targetROF, toLayer);
            if (layer1.empty()) {
              continue;
            }
            const auto ts = mScratch->getROFTimeStamp(fromLayer, pivotROF, toLayer, targetROF);
            if (!ts.isCompatible(pv.getTimeStamp())) {
              continue;
            }
            const auto& targetIndexTable = mScratch->getIndexTable(targetROF, toLayer);
            const int colBinRange = (bins.z - bins.x) + 1;
            for (int iRow = 0; iRow < rowBinsNum; ++iRow) {
              int iRowBin = bins.y + iRow;
              if (indexTableUtils.getCoordType() == IndexTableCoordType::PhiZ) {
                iRowBin %= indexTableUtils.getNrowBins();
              }
              if (iRowBin < 0 || iRowBin >= indexTableUtils.getNrowBins()) {
                break;
              }
              const int firstBinIdx = indexTableUtils.getBinIndex(bins.x, iRowBin);
              const int maxBinIdx = firstBinIdx + colBinRange;
              const int firstRow = targetIndexTable[firstBinIdx];
              const int lastRow = targetIndexTable[maxBinIdx];
              for (int iNext = firstRow; iNext < lastRow; ++iNext) {
                if (iNext >= int(layer1.size())) {
                  break;
                }
                const o2::its::Cluster& nextCluster = layer1[iNext];
                if (mScratch->isClusterUsed(toLayer, nextCluster.clusterId)) {
                  continue;
                }
                const auto& targetMeasurement = mLayerGlobalMeasurements[toLayer][nextCluster.clusterId];

                float tanL = 0.f;
                const bool accepted = acceptTrackletCandidate(searchWindow, sourceMeasurement, currentCluster,
                                                              targetMeasurement, nextCluster, kind,
                                                              mKernelParameters.nSigmaCut, tanL);
                if (accepted) {
                  const float phi{o2::gpu::GPUCommonMath::ATan2(sourceMeasurement.position.y - targetMeasurement.position.y,
                                                                sourceMeasurement.position.x - targetMeasurement.position.x)};
                  if constexpr (decltype(Mode)::value == PassMode::OnePass::value) {
                    tracklets.emplace_back(currentSortedIndex, mScratch->getSortedIndex(targetROF, toLayer, iNext), tanL, phi, ts);
                  } else if constexpr (decltype(Mode)::value == PassMode::TwoPassCount::value) {
                    ++localCount;
                  } else if constexpr (decltype(Mode)::value == PassMode::TwoPassInsert::value) {
                    const int idx = base + offset++;
                    tracklets[idx] = o2::its::Tracklet(currentSortedIndex, mScratch->getSortedIndex(targetROF, toLayer, iNext), tanL, phi, ts);
                  }
                }
              }
            }
          }
        }
      }
      return localCount;
    };

    int dummy{0};
    if (mTaskArena->max_concurrency() <= 1) {
      for (const auto typedLinkId : linkIds) {
        const int linkId = typedLinkId.value();
        const int scratchLinkId = requireScratchLinkSlot(context, iteration, typedLinkId);
        const auto [fromLayer, toLayer] = resolveLinkLayers(linkId);
        const auto kind = topology.getSurface(topology.getLink(typedLinkId).from).kind;
        TrackletProjectionCache linkCache{};
        const auto layerRadii = gsl::span<const float>{mTrkParams[iteration].LayerRadii.data(),
                                                       mTrkParams[iteration].LayerRadii.size()};
        if (!bindTrackletProjectionCache(fromLayer, toLayer, layerRadii, mDiskLayerReferenceZ,
                                         mScratch->getMinR(toLayer), mScratch->getMaxR(toLayer),
                                         mScratch->getMinZ(toLayer), mScratch->getMaxZ(toLayer),
                                         mScratch->getPositionResolution(fromLayer),
                                         mScratch->getLinkMSAngle(scratchLinkId),
                                         mScratch->getLinkPhiCut(scratchLinkId), linkCache)) {
          throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
        }
        const int startROF = 0, endROF = mScratch->getROFTiming(fromLayer).mNROFsTF;
        for (int pivotROF{startROF}; pivotROF < endROF; ++pivotROF) {
          forTracklets(PassMode::OnePass{}, scratchLinkId, fromLayer, toLayer, kind, linkCache, pivotROF, 0, dummy);
        }
      }
    } else {
      tbb::parallel_for(0, static_cast<int>(linkIds.size()), [&](const int linkIndex) {
        const auto typedLinkId = linkIds[linkIndex];
        const int linkId = typedLinkId.value();
        const int scratchLinkId = requireScratchLinkSlot(context, iteration, typedLinkId);
        const auto [fromLayer, toLayer] = resolveLinkLayers(linkId);
        const auto kind = topology.getSurface(topology.getLink(typedLinkId).from).kind;
        TrackletProjectionCache linkCache{};
        const auto layerRadii = gsl::span<const float>{mTrkParams[iteration].LayerRadii.data(),
                                                       mTrkParams[iteration].LayerRadii.size()};
        if (!bindTrackletProjectionCache(fromLayer, toLayer, layerRadii, mDiskLayerReferenceZ,
                                         mScratch->getMinR(toLayer), mScratch->getMaxR(toLayer),
                                         mScratch->getMinZ(toLayer), mScratch->getMaxZ(toLayer),
                                         mScratch->getPositionResolution(fromLayer),
                                         mScratch->getLinkMSAngle(scratchLinkId),
                                         mScratch->getLinkPhiCut(scratchLinkId), linkCache)) {
          throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
        }
        const int startROF = 0, endROF = mScratch->getROFTiming(fromLayer).mNROFsTF;
        bounded_vector<int> perROFCount((endROF - startROF) + 1, mMemoryPool.get());
        tbb::parallel_for(startROF, endROF, [&](const int pivotROF) {
          perROFCount[pivotROF - startROF] = forTracklets(PassMode::TwoPassCount{}, scratchLinkId, fromLayer, toLayer, kind, linkCache, pivotROF, 0, dummy);
        });
        std::exclusive_scan(perROFCount.begin(), perROFCount.end(), perROFCount.begin(), 0);
        const int nTracklets = perROFCount.back();
        mScratch->getTracklets()[scratchLinkId].resize(nTracklets);
        if (nTracklets == 0) {
          return;
        }
        tbb::parallel_for(startROF, endROF, [&](const int pivotROF) {
          int baseIdx = perROFCount[pivotROF - startROF];
          if (baseIdx == perROFCount[pivotROF + 1 - startROF]) {
            return;
          }
          int localIdx = 0;
          forTracklets(PassMode::TwoPassInsert{}, scratchLinkId, fromLayer, toLayer, kind, linkCache, pivotROF, baseIdx, localIdx);
        });
      });
    }

    tbb::parallel_for(0, static_cast<int>(linkIds.size()), [&](const int linkIndex) {
      const int scratchLinkId = requireScratchLinkSlot(context, iteration, linkIds[linkIndex]);
      /// Sort tracklets & remove duplicates
      // duplicates can exist simply since we evaluate per vertex
      auto& trkl{mScratch->getTracklets()[scratchLinkId]};
      std::sort(trkl.begin(), trkl.end());
      trkl.erase(std::unique(trkl.begin(), trkl.end()), trkl.end());
      trkl.shrink_to_fit();
      auto& lut{mScratch->getTrackletsLookupTable()[scratchLinkId]};
      if (!trkl.empty()) {
        for (const auto& tkl : trkl) {
          lut[tkl.firstClusterIndex + 1]++;
        }
        std::inclusive_scan(lut.begin(), lut.end(), lut.begin());
      }
    });

    /// Create tracklets labels
    if (mScratch->hasMCinformation() && mTrkParams[iteration].CreateArtefactLabels) {
      tbb::parallel_for(0, static_cast<int>(linkIds.size()), [&](const int linkIndex) {
        const auto typedLinkId = linkIds[linkIndex];
        const int linkId = typedLinkId.value();
        const int scratchLinkId = requireScratchLinkSlot(context, iteration, typedLinkId);
        const auto [fromLayer, toLayer] = resolveLinkLayers(linkId);
        for (auto& trk : mScratch->getTracklets()[scratchLinkId]) {
          MCCompLabel label;
          int currentId{mScratch->getClusters()[fromLayer][trk.firstClusterIndex].clusterId};
          int nextId{mScratch->getClusters()[toLayer][trk.secondClusterIndex].clusterId};
          for (const auto& lab1 : mScratch->getClusterLabels(fromLayer, currentId)) {
            for (const auto& lab2 : mScratch->getClusterLabels(toLayer, nextId)) {
              if (lab1 == lab2 && lab1.isValid()) {
                label = lab1;
                break;
              }
            }
            if (label.isValid()) {
              break;
            }
          }
          mScratch->getTrackletsLabel(scratchLinkId).emplace_back(label);
        }
      });
    }
  });
}

void TrackerTraits::computeLayerCells(TraversalWorkspaceView& context, const int iteration)
{
  auto& scratch = context.scratch;
  const auto& workspace = context.workspace;
  const auto scratchCellCount = scratch.getCells().size();
  if (scratch.getCellsLookupTable().size() != scratchCellCount) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  for (size_t cellTopologyId = 0; cellTopologyId < scratchCellCount; ++cellTopologyId) {
    deepVectorClear(scratch.getCells()[cellTopologyId]);
    deepVectorClear(scratch.getCellsLookupTable()[cellTopologyId]);
    if (scratch.hasMCinformation() && context.parameters[iteration].CreateArtefactLabels) {
      deepVectorClear(scratch.getCellsLabel(cellTopologyId));
    }
  }

  if (!workspace.valid) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidTraversalSchedule};
  }
  if (!workspace.attachHitConfig.isValid(static_cast<int>(scratch.getNOwnedSurfaces()))) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
  }

  computeLayerCellsImpl(context, iteration, context.binding.getGlobalCells());

  const auto scratchLinkCount = scratch.getTracklets().size();
  for (size_t linkId = 0; linkId < scratchLinkCount; ++linkId) {
    deepVectorClear(scratch.getTracklets()[linkId]);
    deepVectorClear(scratch.getTrackletsLabel(linkId));
  }
}

void TrackerTraits::computeLayerCellsImpl(
  TraversalWorkspaceView& context,
  const int iteration,
  gsl::span<const CellTopologyId> cellIds)
{
  auto* mScratch = &context.scratch;
  const auto& mMemoryPool = mScratch->getMemoryPool();
  const auto* mBinding = &context.binding;
  const auto mTrkParams = context.parameters;
  const auto mBz = context.bz;
  const auto& mTraversalGraph = context.graph;
  const auto& mKernelParameters = context.workspace.kernelParameters;
  const auto& mAttachHitConfig = context.workspace.attachHitConfig;
  const auto& mLayerMaterial = context.workspace.layerMaterial;
  const auto& mLayerMeasurements = context.workspace.layerMeasurements;
  const auto& mLayerGlobalMeasurements = context.workspace.layerGlobalMeasurements;
  const auto& topology = mTraversalGraph;

  mTaskArena->execute([&] {
    struct CellHitBinding {
      std::array<int, 3> layers{};
      std::array<SurfaceDescriptor, 3> surfaces{};
    };

    auto resolveCellHitBinding = [&](const auto& cellTopology) -> CellHitBinding {
      const auto& firstLink = topology.getLink(cellTopology.firstLink);
      const auto& secondLink = topology.getLink(cellTopology.secondLink);
      const std::array<SurfaceId, 3> surfaces{firstLink.from, firstLink.to, secondLink.to};
      CellHitBinding binding;
      for (int i = 0; i < 3; ++i) {
        const auto surfaceId = surfaces[i];
        binding.layers[i] = requireSurfacePosition(context, iteration, surfaceId);
        binding.surfaces[i] = topology.getSurface(surfaceId);
      }
      return binding;
    };

    auto forTrackletCells = [&](auto Mode, SurfaceKind kind, int firstLinkId, int secondLinkId, const CellHitBinding& hitBinding, bounded_vector<CellSeed>& layerCells, int iTracklet, int offset = 0) -> int {
      const auto& hitLayers = hitBinding.layers;
      const o2::its::Tracklet& currentTracklet{mScratch->getTracklets()[firstLinkId][iTracklet]};
      const int nextLayerClusterIndex{currentTracklet.secondClusterIndex};
      const int nextLayerFirstTrackletIndex{mScratch->getTrackletsLookupTable()[secondLinkId][nextLayerClusterIndex]};
      const int nextLayerLastTrackletIndex{mScratch->getTrackletsLookupTable()[secondLinkId][nextLayerClusterIndex + 1]};
      int foundCells{0};
      for (int iNextTracklet{nextLayerFirstTrackletIndex}; iNextTracklet < nextLayerLastTrackletIndex; ++iNextTracklet) {
        const o2::its::Tracklet& nextTracklet{mScratch->getTracklets()[secondLinkId][iNextTracklet]};
        if (nextTracklet.firstClusterIndex != nextLayerClusterIndex) {
          break;
        }
        if (!currentTracklet.getTimeStamp().isCompatible(nextTracklet.getTimeStamp())) {
          continue;
        }

        /// Prepare the track seed; clusters are numbered from inner to outer.
        const int clusId[3]{
          mScratch->getClusters()[hitLayers[0]][currentTracklet.firstClusterIndex].clusterId,
          mScratch->getClusters()[hitLayers[1]][nextTracklet.firstClusterIndex].clusterId,
          mScratch->getClusters()[hitLayers[2]][nextTracklet.secondClusterIndex].clusterId};

        const auto& measurementInner = mLayerMeasurements[hitLayers[0]][clusId[0]];
        const auto& measurementMiddle = mLayerMeasurements[hitLayers[1]][clusId[1]];
        const auto& measurementOuter = mLayerMeasurements[hitLayers[2]][clusId[2]];
        const auto& globalInner = mLayerGlobalMeasurements[hitLayers[0]][clusId[0]];
        const auto& globalMiddle = mLayerGlobalMeasurements[hitLayers[1]][clusId[1]];
        const auto& globalOuter = mLayerGlobalMeasurements[hitLayers[2]][clusId[2]];
        const double linkMSAngle = static_cast<double>(mScratch->getLinkMSAngle(secondLinkId));
        const DirectionProcessNoise directionProcessNoise{linkMSAngle * linkMSAngle};
        std::array<TransverseDirectionObservation, 3> transverseObservations{};
        if (!makeTransverseDirectionObservation(globalInner, transverseObservations[0]) ||
            !makeTransverseDirectionObservation(globalMiddle, transverseObservations[1]) ||
            !makeTransverseDirectionObservation(globalOuter, transverseObservations[2])) {
          continue;
        }
        TransverseDirectionCompatibility transverseCompatibility{};
        if (!trackletDirectionsAreTransverselyCompatible(
              transverseObservations, currentTracklet.phi, nextTracklet.phi,
              directionProcessNoise, mBz,
              mKernelParameters.trackletMinPt, mKernelParameters.nSigmaCut,
              transverseCompatibility)) {
          continue;
        }
        std::array<DirectionObservation, 3> directionObservations{};
        if (!makeDirectionObservation(globalInner, directionObservations[0]) ||
            !makeDirectionObservation(globalMiddle, directionObservations[1]) ||
            !makeDirectionObservation(globalOuter, directionObservations[2])) {
          continue;
        }
        CellDirectionCompatibility directionCompatibility{};
        if (cellDirectionsAreCompatible(directionObservations, directionProcessNoise, mKernelParameters.nSigmaCut,
                                        directionCompatibility)) {

          // Strictly {inner, middle, outer}: Cylinder reads [1] then
          // [0] (outer slot unused), Disk reads [2], [1], [0].
          const std::array<NominalSurfaceMaterial, 3> material{
            mAttachHitConfig.layerMaterial[hitLayers[0]],
            mAttachHitConfig.layerMaterial[hitLayers[1]],
            mAttachHitConfig.layerMaterial[hitLayers[2]]};

          SurfaceKinematicState state{};
          float chi2{0.f};
          OperationFailureReason buildReason{};
          const bool good = buildCellSeed(kind, globalInner, globalMiddle, globalOuter,
                                          measurementInner, measurementMiddle, measurementOuter,
                                          material, mBz, kCompatibilityAbsCharge, kCompatibilityPID,
                                          state, chi2, mKernelParameters, buildReason);

          if (good) {
            TimeEstBC ts = currentTracklet.getTimeStamp();
            ts += nextTracklet.getTimeStamp();
            // Build directly from the resolved plan positions; plan validation
            // already checked them against the cell's hit-surface mask.
            const LayerMask hitLayerMask{hitLayers[0], hitLayers[1], hitLayers[2]};
            TripletFitFactor tripletFactor{};
            if constexpr (decltype(Mode)::value != PassMode::TwoPassCount::value) {
              std::array<TripletFitObservation, 3> observations{};
              if (!makeTripletFitObservation(globalInner, observations[0]) ||
                  !makeTripletFitObservation(globalMiddle, observations[1]) ||
                  !makeTripletFitObservation(globalOuter, observations[2]) ||
                  !makeTripletFitFactor(observations, tripletFactor)) {
                tripletFactor = {};
              }
            }
            if constexpr (decltype(Mode)::value == PassMode::OnePass::value) {
              layerCells.emplace_back(hitLayerMask, clusId[0], clusId[1], clusId[2], iTracklet, iNextTracklet, state, chi2, ts);
              layerCells.back().tripletFactor() = tripletFactor;
              ++foundCells;
            } else if constexpr (decltype(Mode)::value == PassMode::TwoPassCount::value) {
              ++foundCells;
            } else if constexpr (decltype(Mode)::value == PassMode::TwoPassInsert::value) {
              layerCells[offset++] = CellSeed(hitLayerMask, clusId[0], clusId[1], clusId[2], iTracklet, iNextTracklet, state, chi2, ts);
              layerCells[offset - 1].tripletFactor() = tripletFactor;
              ++foundCells;
            } else {
              static_assert(false, "Unknown mode!");
            }
          }
        }
      }
      return foundCells;
    };

    for (const auto typedCellId : cellIds) {
      const int cellTopologyId = requireScratchCellSlot(context, iteration, typedCellId);
      const auto& cellTopology = topology.getCell(typedCellId);
      const auto kind = topology.getSurface(topology.getLink(cellTopology.firstLink).from).kind;
      const int firstLinkId = requireScratchLinkSlot(context, iteration, cellTopology.firstLink);
      const int secondLinkId = requireScratchLinkSlot(context, iteration, cellTopology.secondLink);
      if (mScratch->getTracklets()[firstLinkId].empty() ||
          mScratch->getTracklets()[secondLinkId].empty()) {
        continue;
      }
      const auto hitBinding = resolveCellHitBinding(cellTopology);

      auto& layerCells = mScratch->getCells()[cellTopologyId];
      const int currentLayerTrackletsNum{static_cast<int>(mScratch->getTracklets()[firstLinkId].size())};
      bounded_vector<int> perTrackletCount(currentLayerTrackletsNum + 1, 0, mMemoryPool.get());
      if (mTaskArena->max_concurrency() <= 1) {
        for (int iTracklet{0}; iTracklet < currentLayerTrackletsNum; ++iTracklet) {
          perTrackletCount[iTracklet] = forTrackletCells(PassMode::OnePass{}, kind, firstLinkId, secondLinkId, hitBinding, layerCells, iTracklet);
        }
        std::exclusive_scan(perTrackletCount.begin(), perTrackletCount.end(), perTrackletCount.begin(), 0);
      } else {
        tbb::parallel_for(0, currentLayerTrackletsNum, [&](const int iTracklet) {
          perTrackletCount[iTracklet] = forTrackletCells(PassMode::TwoPassCount{}, kind, firstLinkId, secondLinkId, hitBinding, layerCells, iTracklet);
        });

        std::exclusive_scan(perTrackletCount.begin(), perTrackletCount.end(), perTrackletCount.begin(), 0);
        auto totalCells{perTrackletCount.back()};
        if (totalCells == 0) {
          auto& lut = mScratch->getCellsLookupTable()[cellTopologyId];
          lut.resize(currentLayerTrackletsNum + 1);
          std::fill(lut.begin(), lut.end(), 0);
          continue;
        }
        layerCells.resize(totalCells);

        tbb::parallel_for(0, currentLayerTrackletsNum, [&](const int iTracklet) {
          int offset = perTrackletCount[iTracklet];
          if (offset == perTrackletCount[iTracklet + 1]) {
            return;
          }
          forTrackletCells(PassMode::TwoPassInsert{}, kind, firstLinkId, secondLinkId, hitBinding, layerCells, iTracklet, offset);
        });
      }

      auto& lut = mScratch->getCellsLookupTable()[cellTopologyId];
      lut.resize(currentLayerTrackletsNum + 1);
      std::copy_n(perTrackletCount.begin(), currentLayerTrackletsNum + 1, lut.begin());

      if (mScratch->hasMCinformation() && mTrkParams[iteration].CreateArtefactLabels) {
        auto& labels = mScratch->getCellsLabel(cellTopologyId);
        labels.reserve(layerCells.size());
        for (const auto& cell : layerCells) {
          MCCompLabel currentLab{mScratch->getTrackletsLabel(firstLinkId)[cell.getFirstTrackletIndex()]};
          MCCompLabel nextLab{mScratch->getTrackletsLabel(secondLinkId)[cell.getSecondTrackletIndex()]};
          labels.emplace_back(currentLab == nextLab ? currentLab : MCCompLabel());
        }
      }
    }
  });
}

void TrackerTraits::findCellsNeighbours(TraversalWorkspaceView& context, const int iteration)
{
  auto& scratch = context.scratch;
  const auto& workspace = context.workspace;
  if (!workspace.valid) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidTraversalSchedule};
  }
  for (std::size_t slot = 0; slot < scratch.getCellsNeighbours().size(); ++slot) {
    deepVectorClear(scratch.getCellsNeighbours()[slot]);
    deepVectorClear(scratch.getCellsNeighboursTopology()[slot]);
    deepVectorClear(scratch.getCellsNeighboursLUT()[slot]);
  }
  const auto scheduledCells = context.binding.getGlobalScheduledCells();
  if (!scheduledCells.empty()) {
    findCellsNeighboursForSchedule(context, iteration, scheduledCells, workspace.kernelParameters);
  }
  for (auto& cellLUT : scratch.getCellsLookupTable()) {
    deepVectorClear(cellLUT);
  }
}

void TrackerTraits::findCellsNeighboursForSchedule(
  TraversalWorkspaceView& context,
  int iteration,
  gsl::span<const CellTopologyId> scheduledCells,
  const TrackingKernelParameters& params)
{
  auto* mScratch = &context.scratch;
  const auto& mMemoryPool = mScratch->getMemoryPool();
  const auto& mTraversalGraph = context.graph;
  const auto& mLayerMaterial = context.workspace.layerMaterial;
  const auto& mLayerGlobalMeasurements = context.workspace.layerGlobalMeasurements;
  const auto topology = mTraversalGraph;
  const auto scratchCellCount = mScratch->getCells().size();
  if (mScratch->getCellsLookupTable().size() != scratchCellCount ||
      mScratch->getCellsNeighbours().size() != scratchCellCount ||
      mScratch->getCellsNeighboursTopology().size() != scratchCellCount ||
      mScratch->getCellsNeighboursLUT().size() != scratchCellCount) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  mTaskArena->execute([&] {
    std::vector<bounded_vector<CellNeighbour>> cellsNeighboursByTarget;
    cellsNeighboursByTarget.reserve(scratchCellCount);
    for (size_t cellTopologyId = 0; cellTopologyId < scratchCellCount; ++cellTopologyId) {
      cellsNeighboursByTarget.emplace_back(mMemoryPool.get());
    }

    for (const auto scheduledId : scheduledCells) {
      const int cellTopologyId = requireScratchCellSlot(context, iteration, scheduledId);
      if (static_cast<size_t>(cellTopologyId) >= scratchCellCount ||
          static_cast<size_t>(cellTopologyId) >= mScratch->getCellsLookupTable().size()) {
        throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
      }
      const auto& cellTopology = topology.getCell(scheduledId);
      const int currentProcessLink = requireScratchLinkSlot(context, iteration, cellTopology.secondLink);
      const double currentMSAngle = mScratch->getLinkMSAngle(currentProcessLink);
      const double currentAngularVariance = currentMSAngle * currentMSAngle;
      if (mScratch->getCells()[cellTopologyId].empty()) {
        continue;
      }
      const auto successors = topology.getCellsStartingWithLink(cellTopology.secondLink);
      if (!successors.getEntries()) {
        continue;
      }

      tbb::enumerable_thread_specific<bounded_vector<CellNeighbour>> sourceNeighbours([&]() { return bounded_vector<CellNeighbour>{mMemoryPool.get()}; });
      tbb::parallel_for(0, static_cast<int>(mScratch->getCells()[cellTopologyId].size()), [&](const int iCell) {
        auto& localNeighbours = sourceNeighbours.local();
        const auto& currentCellSeed{mScratch->getCells()[cellTopologyId][iCell]};
        const int nextLayerTrackletIndex{currentCellSeed.getSecondTrackletIndex()};
        for (uint32_t iSuccessor = 0; iSuccessor < successors.getEntries(); ++iSuccessor) {
          // Translate dynamically discovered neighbours from the global CSR
          // topology before accessing scratch.
          const auto nextTopologyId = topology.cellsByFirstLink[successors.getFirstEntry() + iSuccessor];
          const int nextCellTopologyId = requireScratchCellSlot(context, iteration, nextTopologyId);
          const auto& nextCellTopology = topology.getCell(nextTopologyId);
          const int nextProcessLink = requireScratchLinkSlot(context, iteration, nextCellTopology.secondLink);
          const double nextMSAngle = mScratch->getLinkMSAngle(nextProcessLink);
          const double nextAngularVariance = nextMSAngle * nextMSAngle;
          if (static_cast<size_t>(nextCellTopologyId) >= mScratch->getCells().size() ||
              static_cast<size_t>(nextCellTopologyId) >= mScratch->getCellsLookupTable().size()) {
            throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
          }
          if (mScratch->getCells()[nextCellTopologyId].empty() ||
              mScratch->getCellsLookupTable()[nextCellTopologyId].empty()) {
            continue;
          }
          const auto& nextCellLUT = mScratch->getCellsLookupTable()[nextCellTopologyId];
          if (nextLayerTrackletIndex < 0 || nextLayerTrackletIndex + 1 >= static_cast<int>(nextCellLUT.size())) {
            continue;
          }
          const int nextLayerFirstCellIndex{nextCellLUT[nextLayerTrackletIndex]};
          const int nextLayerLastCellIndex{nextCellLUT[nextLayerTrackletIndex + 1]};
          if (nextLayerFirstCellIndex < 0 || nextLayerLastCellIndex < nextLayerFirstCellIndex ||
              nextLayerLastCellIndex > static_cast<int>(mScratch->getCells()[nextCellTopologyId].size())) {
            throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
          }
          for (int iNextCell{nextLayerFirstCellIndex}; iNextCell < nextLayerLastCellIndex; ++iNextCell) {
            const auto& nextCellSeedRef{mScratch->getCells()[nextCellTopologyId][iNextCell]};
            if (nextCellSeedRef.getFirstTrackletIndex() != nextLayerTrackletIndex || !currentCellSeed.getTimeStamp().isCompatible(nextCellSeedRef.getTimeStamp())) {
              break;
            }

            const auto currentMiddle = currentCellSeed.getClusterReference(1);
            const auto currentOuter = currentCellSeed.getClusterReference(2);
            const auto nextInner = nextCellSeedRef.getClusterReference(0);
            const auto nextMiddle = nextCellSeedRef.getClusterReference(1);
            if (currentMiddle.surfacePosition != nextInner.surfacePosition ||
                currentMiddle.clusterIndex != nextInner.clusterIndex ||
                currentOuter.surfacePosition != nextMiddle.surfacePosition ||
                currentOuter.clusterIndex != nextMiddle.clusterIndex) {
              continue;
            }

            const std::array<CellClusterReference, 4> references{
              currentCellSeed.getClusterReference(0), currentMiddle,
              currentOuter, nextCellSeedRef.getClusterReference(2)};
            std::array<TripletFitObservation, 4> observations{};
            bool observationsValid = true;
            for (std::size_t hit = 0; hit < references.size(); ++hit) {
              const auto reference = references[hit];
              if (reference.surfacePosition < 0 ||
                  static_cast<std::size_t>(reference.surfacePosition) >= mLayerGlobalMeasurements.size() ||
                  reference.clusterIndex < 0 ||
                  static_cast<std::size_t>(reference.clusterIndex) >= mLayerGlobalMeasurements[reference.surfacePosition].size() ||
                  !makeTripletFitObservation(
                    mLayerGlobalMeasurements[reference.surfacePosition][reference.clusterIndex], observations[hit])) {
                observationsValid = false;
                break;
              }
            }
            AdjacentTripletFitResult adjacentFit{};
            const bool fitValid = observationsValid &&
                                  fitAdjacentTripletFactors(
                                    currentCellSeed.tripletFactor(), nextCellSeedRef.tripletFactor(), observations,
                                    {currentAngularVariance, nextAngularVariance}, adjacentFit);
            if (!fitValid || adjacentFit.chi2 > params.maxChi2ClusterAttachment) {
              continue;
            }

            const int nextLevel = currentCellSeed.getLevel() + 1;
            localNeighbours.emplace_back(cellTopologyId, iCell, nextCellTopologyId, iNextCell, nextLevel);
          }
        }
      });

      bounded_vector<size_t> count(scratchCellCount, 0, mMemoryPool.get());
      for (const auto& localNeighbours : sourceNeighbours) {
        for (const auto& neigh : localNeighbours) {
          ++count[neigh.nextCellTopology];
        }
      }
      for (size_t i{0}; i < scratchCellCount; ++i) {
        cellsNeighboursByTarget[i].reserve(count[i]);
      }
      for (const auto& localNeighbours : sourceNeighbours) {
        for (const auto& neigh : localNeighbours) {
          cellsNeighboursByTarget[neigh.nextCellTopology].emplace_back(neigh);
          if (neigh.level > mScratch->getCells()[neigh.nextCellTopology][neigh.nextCell].getLevel()) {
            mScratch->getCells()[neigh.nextCellTopology][neigh.nextCell].setLevel(neigh.level);
          }
        }
      }
    }

    for (size_t cellTopologyId = 0; cellTopologyId < scratchCellCount; ++cellTopologyId) {
      auto& cellsNeighbours = cellsNeighboursByTarget[cellTopologyId];
      if (cellsNeighbours.empty()) {
        continue;
      }

      std::sort(cellsNeighbours.begin(), cellsNeighbours.end(), [](const auto& a, const auto& b) {
        return a.nextCell < b.nextCell;
      });

      auto& cellsNeighbourLUT = mScratch->getCellsNeighboursLUT()[cellTopologyId];
      cellsNeighbourLUT.assign(mScratch->getCells()[cellTopologyId].size(), 0);
      for (const auto& neigh : cellsNeighbours) {
        ++cellsNeighbourLUT[neigh.nextCell];
      }
      std::inclusive_scan(cellsNeighbourLUT.begin(), cellsNeighbourLUT.end(), cellsNeighbourLUT.begin());

      mScratch->getCellsNeighbours()[cellTopologyId].reserve(cellsNeighbours.size());
      mScratch->getCellsNeighboursTopology()[cellTopologyId].reserve(cellsNeighbours.size());
      std::ranges::transform(cellsNeighbours, std::back_inserter(mScratch->getCellsNeighbours()[cellTopologyId]), [](const auto& neigh) { return neigh.cell; });
      std::ranges::transform(cellsNeighbours, std::back_inserter(mScratch->getCellsNeighboursTopology()[cellTopologyId]), [](const auto& neigh) { return neigh.cellTopology; });
    }
  });
}

template <typename InputSeed>
void TrackerTraits::processNeighbours(TraversalWorkspaceView& context, int iteration, int defaultCellTopologyId, int iLevel, const bounded_vector<InputSeed>& currentCellSeed, const bounded_vector<int>& currentCellId, const bounded_vector<int>& currentCellTopologyId, bounded_vector<TrackSeed>& updatedCellSeeds, bounded_vector<int>& updatedCellsIds, bounded_vector<int>& updatedCellsTopologyIds, const TrackingKernelParameters& params)
{
  auto* mScratch = &context.scratch;
  const auto& mMemoryPool = mScratch->getMemoryPool();
  const auto mBz = context.bz;
  const auto& mAttachHitConfig = context.workspace.attachHitConfig;
  const auto& mLayerMeasurements = context.workspace.layerMeasurements;
  const auto layerMaterial = mAttachHitConfig.layerMaterial;
  const int activeSurfaceCount = static_cast<int>(mScratch->getNOwnedSurfaces());

  mTaskArena->execute([&] {
    auto forCellNeighbours = [&](auto Mode, int iCell, int offset = 0) -> int {
      const auto& currentCell{currentCellSeed[iCell]};
      const int cellTopologyId = currentCellTopologyId.empty() ? defaultCellTopologyId : currentCellTopologyId[iCell];

      if constexpr (decltype(Mode)::value != PassMode::TwoPassInsert::value) {
        if (currentCell.getLevel() != iLevel) {
          return 0;
        }
        if (currentCellId.empty()) {
          for (int layer = 0; layer < activeSurfaceCount; ++layer) {
            const int clusterIndex = currentCell.getCluster(layer);
            if (clusterIndex != o2::its::constants::UnusedIndex && mScratch->isClusterUsed(layer, clusterIndex)) {
              return 0; /// this we do only on the first iteration, hence the check on currentCellId
            }
          }
        }
      }

      const int cellId = currentCellId.empty() ? iCell : currentCellId[iCell];
      if (cellTopologyId < 0 || mScratch->getCellsNeighboursLUT()[cellTopologyId].empty()) {
        return 0;
      }
      const int startNeighbourId{cellId ? mScratch->getCellsNeighboursLUT()[cellTopologyId][cellId - 1] : 0};
      const int endNeighbourId{mScratch->getCellsNeighboursLUT()[cellTopologyId][cellId]};
      int foundSeeds{0};
      for (int iNeighbourCell{startNeighbourId}; iNeighbourCell < endNeighbourId; ++iNeighbourCell) {
        const int neighbourCellTopologyId = mScratch->getCellsNeighboursTopology()[cellTopologyId][iNeighbourCell];
        const int neighbourCellId = mScratch->getCellsNeighbours()[cellTopologyId][iNeighbourCell];
        const auto& neighbourCell = mScratch->getCells()[neighbourCellTopologyId][neighbourCellId];
        if (neighbourCell.getSecondTrackletIndex() != currentCell.getFirstTrackletIndex()) {
          continue;
        }
        if (!currentCell.getTimeStamp().isCompatible(neighbourCell.getTimeStamp())) {
          continue;
        }
        if (currentCell.getLevel() - 1 != neighbourCell.getLevel()) {
          continue;
        }
        const int neighbourLayer = neighbourCell.getInnerLayer();
        const int neighbourCluster = neighbourCell.getFirstClusterIndex();
        if (mScratch->isClusterUsed(neighbourLayer, neighbourCluster)) {
          continue;
        }

        /// Let's start the fitting procedure
        TrackSeed seed{currentCell};
        seed.getTimeStamp() = currentCell.getTimeStamp();
        seed.getTimeStamp() += neighbourCell.getTimeStamp();

        const auto& measurement = mLayerMeasurements[neighbourLayer][neighbourCluster];
        float chi2 = seed.getChi2();
        OperationFailureReason attachReason{};
        const bool attached = Propagator::attachMeasurement(seed.state(), measurement, layerMaterial[neighbourLayer], mBz,
                                                            material::MaterialTraversalDirection::OppositeMomentum, true,
                                                            params.maxChi2ClusterAttachment, chi2, attachReason);
        if (!attached) {
          continue;
        }
        seed.setChi2(chi2);

        if constexpr (decltype(Mode)::value != PassMode::TwoPassCount::value) {
          seed.setCluster(neighbourLayer, neighbourCluster);
          if (neighbourLayer < 0 || neighbourLayer >= activeSurfaceCount) {
            throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
          }
          // TrackSeed::SurfaceMask is the fixed-capacity compact plan
          // position space used by the CA acceptance/refit loops. Global
          // SurfaceIds stay on normalized measurements and GenericTrack
          // references; mixing them into this local mask would make a
          // combined sparse binding look like extra holes. The binding's
          // ordered position is already `neighbourLayer` here.
          auto surfaceMask = seed.getSurfaceMask();
          surfaceMask.set(SurfaceId{static_cast<uint16_t>(neighbourLayer)});
          seed.setSurfaceMask(surfaceMask);
          seed.setLevel(neighbourCell.getLevel());
          seed.setFirstTrackletIndex(neighbourCell.getFirstTrackletIndex());
          seed.setSecondTrackletIndex(neighbourCell.getSecondTrackletIndex());
        }

        if constexpr (decltype(Mode)::value == PassMode::OnePass::value) {
          updatedCellSeeds.push_back(seed);
          updatedCellsIds.push_back(neighbourCellId);
          updatedCellsTopologyIds.push_back(neighbourCellTopologyId);
        } else if constexpr (decltype(Mode)::value == PassMode::TwoPassCount::value) {
          ++foundSeeds;
        } else if constexpr (decltype(Mode)::value == PassMode::TwoPassInsert::value) {
          updatedCellSeeds[offset] = seed;
          updatedCellsIds[offset] = neighbourCellId;
          updatedCellsTopologyIds[offset++] = neighbourCellTopologyId;
        } else {
          static_assert(false, "Unknown mode!");
        }
      }
      return foundSeeds;
    };

    const int nCells = static_cast<int>(currentCellSeed.size());
    if (mTaskArena->max_concurrency() <= 1) {
      for (int iCell{0}; iCell < nCells; ++iCell) {
        forCellNeighbours(PassMode::OnePass{}, iCell);
      }
    } else {
      bounded_vector<int> perCellCount(nCells + 1, 0, mMemoryPool.get());
      tbb::parallel_for(0, nCells, [&](const int iCell) {
        perCellCount[iCell] = forCellNeighbours(PassMode::TwoPassCount{}, iCell);
      });

      std::exclusive_scan(perCellCount.begin(), perCellCount.end(), perCellCount.begin(), 0);
      auto totalNeighbours{perCellCount.back()};
      if (totalNeighbours == 0) {
        return;
      }
      updatedCellSeeds.resize(totalNeighbours);
      updatedCellsIds.resize(totalNeighbours);
      updatedCellsTopologyIds.resize(totalNeighbours);

      tbb::parallel_for(0, nCells, [&](const int iCell) {
        int offset = perCellCount[iCell];
        if (offset == perCellCount[iCell + 1]) {
          return;
        }
        forCellNeighbours(PassMode::TwoPassInsert{}, iCell, offset);
      });
    }
  });
}

void TrackerTraits::findRoads(TraversalWorkspaceView& context, const int iteration, SeedRefitFunction refitFunction)
{
  if (!context.workspace.valid) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidTraversalSchedule};
  }
  if (context.scratch.getCells().size() != context.binding.getGlobalCells().size()) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  findRoadsImpl(context, iteration, refitFunction);
}

void TrackerTraits::findRoadsImpl(TraversalWorkspaceView& context, const int iteration, SeedRefitFunction refitFunction)
{
  auto* mScratch = &context.scratch;
  const auto& mMemoryPool = mScratch->getMemoryPool();
  const auto* mBinding = &context.binding;
  const auto mTrkParams = context.parameters;
  const auto mBz = context.bz;
  const auto& mTraversalGraph = context.graph;
  const auto& mKernelParameters = context.workspace.kernelParameters;
  const auto& mLayerMeasurements = context.workspace.layerMeasurements;
  const auto& mLayerGlobalMeasurements = context.workspace.layerGlobalMeasurements;
  const auto roadStartCells = mBinding->getGlobalRoadStartCells();
  const int activeSurfaceCount = static_cast<int>(mScratch->getNOwnedSurfaces());
  bounded_vector<bounded_vector<int>> firstClusters(activeSurfaceCount, bounded_vector<int>(mMemoryPool.get()), mMemoryPool.get());
  firstClusters.resize(activeSurfaceCount);
  // Road starts are the binding's seeding-eligible sparse-plan subsequence.
  // CellTopologyId values use compact slots; SurfaceId is never a vector index.
  // Filter roads by absolute q/pT in parameters[4]'s units, identically for
  // both families. Non-finite values fail the finite-bound comparison.
  constexpr float maxAbsQOverPt = 1.e3f;
  const int cellsPerRoad = static_cast<int>(mScratch->getNOwnedSurfaces()) - 2;
  const auto componentOffsets = mBinding->getRoadStartComponentOffsets();
  if (componentOffsets.empty() || componentOffsets.front() != 0 || componentOffsets.back() != roadStartCells.size()) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  for (size_t component = 0; component + 1 < componentOffsets.size(); ++component) {
    const auto componentRoadStarts = roadStartCells.subspan(componentOffsets[component],
                                                            componentOffsets[component + 1] - componentOffsets[component]);
    for (int startLevel{cellsPerRoad}; startLevel >= mTrkParams[iteration].CellMinimumLevel(); --startLevel) {

      auto seedFilter = [&](const auto& seed) {
        return seed.getHitLayerMask().isAllowed(mTrkParams[iteration].MaxHoles, mTrkParams[iteration].HoleLayerMask) &&
               seed.getHitLayerMask().length() >= mTrkParams[iteration].MinTrackLength &&
               std::abs(seed.getQOverPt()) <= maxAbsQOverPt && seed.getChi2() <= mTrkParams[iteration].MaxChi2NDF * ((startLevel + 2) * 2 - 5);
      };

      bounded_vector<TrackSeed> trackSeeds(mMemoryPool.get());
      // The binding supplies the ownership-filtered road-start span.
      for (const auto startId : componentRoadStarts) {
        // Translate the road-start cell to its compact slot once.
        const int startCellTopologyId = requireScratchCellSlot(context, iteration, startId);
        // Cell population is per-event/per-vertex data, so check it against
        // the current vertex rather than caching it in SurfacePlanBinding.
        if (mScratch->getCells()[startCellTopologyId].empty()) {
          continue;
        }

        bounded_vector<int> lastCellId(mMemoryPool.get()), updatedCellId(mMemoryPool.get());
        bounded_vector<int> lastCellTopologyId(mMemoryPool.get()), updatedCellTopologyId(mMemoryPool.get());
        bounded_vector<TrackSeed> lastCellSeed(mMemoryPool.get()), updatedCellSeed(mMemoryPool.get());

        processNeighbours(context, iteration, startCellTopologyId, startLevel, mScratch->getCells()[startCellTopologyId], lastCellId, lastCellTopologyId, updatedCellSeed, updatedCellId, updatedCellTopologyId, mKernelParameters);

        int level = startLevel;
        while (level > 2 && !updatedCellSeed.empty()) {
          lastCellSeed.swap(updatedCellSeed);
          lastCellId.swap(updatedCellId);
          lastCellTopologyId.swap(updatedCellTopologyId);
          deepVectorClear(updatedCellSeed); /// tame the memory peaks
          deepVectorClear(updatedCellId);   /// tame the memory peaks
          deepVectorClear(updatedCellTopologyId);
          processNeighbours(context, iteration, o2::its::constants::UnusedIndex, --level, lastCellSeed, lastCellId, lastCellTopologyId, updatedCellSeed, updatedCellId, updatedCellTopologyId, mKernelParameters);
        }
        deepVectorClear(lastCellId);         /// tame the memory peaks
        deepVectorClear(lastCellTopologyId); /// tame the memory peaks
        deepVectorClear(lastCellSeed);       /// tame the memory peaks

        if (!updatedCellSeed.empty()) {
          trackSeeds.reserve(trackSeeds.size() + std::count_if(updatedCellSeed.begin(), updatedCellSeed.end(), seedFilter));
          std::copy_if(updatedCellSeed.begin(), updatedCellSeed.end(), std::back_inserter(trackSeeds), seedFilter);
        }
      }

      if (trackSeeds.empty()) {
        continue;
      }

      bounded_vector<TrackingCandidate> tracks(mMemoryPool.get());
      mTaskArena->execute([&] {
        auto forSeed = [&](auto Mode, int iSeed, int offset = 0) {
          TrackingCandidate temporaryTrack;
          temporaryTrack.seed = trackSeeds[iSeed];
          const bool refitSuccess = refitFunction(trackSeeds[iSeed],
                                                  mTrkParams[iteration],
                                                  mBz,
                                                  *mScratch,
                                                  mLayerGlobalMeasurements,
                                                  mLayerMeasurements,
                                                  mTraversalGraph.getSurfaceCatalogView(),
                                                  mBinding->getOrderedSurfaces(),
                                                  temporaryTrack);
          if (refitSuccess) {
            if constexpr (decltype(Mode)::value == PassMode::OnePass::value) {
              tracks.push_back(temporaryTrack);
            } else if constexpr (decltype(Mode)::value == PassMode::TwoPassCount::value) {
              // nothing to do
            } else if constexpr (decltype(Mode)::value == PassMode::TwoPassInsert::value) {
              tracks[offset] = temporaryTrack;
            } else {
              static_assert(false, "Unknown mode!");
            }
            return 1;
          }
          return 0;
        };

        const int nSeeds = static_cast<int>(trackSeeds.size());
        if (mTaskArena->max_concurrency() <= 1) {
          for (int iSeed{0}; iSeed < nSeeds; ++iSeed) {
            forSeed(PassMode::OnePass{}, iSeed);
          }
        } else {
          bounded_vector<int> perSeedCount(nSeeds + 1, 0, mMemoryPool.get());
          tbb::parallel_for(0, nSeeds, [&](const int iSeed) {
            perSeedCount[iSeed] = forSeed(PassMode::TwoPassCount{}, iSeed);
          });

          std::exclusive_scan(perSeedCount.begin(), perSeedCount.end(), perSeedCount.begin(), 0);
          auto totalTracks{perSeedCount.back()};
          if (totalTracks == 0) {
            return;
          }
          tracks.resize(totalTracks);

          tbb::parallel_for(0, nSeeds, [&](const int iSeed) {
            if (perSeedCount[iSeed] == perSeedCount[iSeed + 1]) {
              return;
            }
            forSeed(PassMode::TwoPassInsert{}, iSeed, perSeedCount[iSeed]);
          });
        }

        deepVectorClear(trackSeeds);
      });

      // Same ordering as o2::its::track::isBetter (longer track, then lower chi2).
      std::sort(tracks.begin(), tracks.end(), [](const TrackingCandidate& a, const TrackingCandidate& b) {
        const auto ncla = a.getNumberOfClusters();
        const auto nclb = b.getNumberOfClusters();
        return (ncla == nclb) ? (a.track.chi2 < b.track.chi2) : ncla > nclb;
      });
      acceptTracks(context, iteration, tracks, firstClusters);
    }
  }
}

void TrackerTraits::acceptTracks(TraversalWorkspaceView& context, int iteration,
                                 bounded_vector<TrackingCandidate>& tracks,
                                 bounded_vector<bounded_vector<int>>& firstClusters)
{
  auto* mScratch = &context.scratch;
  auto* mFrame = &context.frame;
  const auto mTrkParams = context.parameters;
  auto& mAcceptedTracksForSharedStatus = context.workspace.acceptedTracks;
  const auto& mLayerGlobalMeasurements = context.workspace.layerGlobalMeasurements;
  auto& trks = mAcceptedTracksForSharedStatus;
  trks.reserve(trks.size() + tracks.size());
  const int activeSurfaceCount = static_cast<int>(mScratch->getNOwnedSurfaces());
  for (auto& track : tracks) {
    int nShared = 0;
    bool isFirstShared{false};
    int firstLayer{-1}, firstCluster{-1};
    for (int iLayer{0}; iLayer < activeSurfaceCount; ++iLayer) {
      if (track.getClusterIndex(iLayer) == o2::its::constants::UnusedIndex) {
        continue;
      }
      bool isShared = mScratch->isClusterUsed(iLayer, track.getClusterIndex(iLayer));
      nShared += int(isShared);
      if (firstLayer < 0) {
        firstCluster = track.getClusterIndex(iLayer);
        isFirstShared = isShared && mTrkParams[iteration].AllowSharingFirstCluster && std::find(firstClusters[iLayer].begin(), firstClusters[iLayer].end(), firstCluster) != firstClusters[iLayer].end();
        firstLayer = iLayer;
      }
    }

    /// do not account for the first cluster in the shared clusters number if it is allowed
    if (nShared - int(isFirstShared && mTrkParams[iteration].AllowSharingFirstCluster) > mTrkParams[iteration].SharedMaxClusters) {
      continue;
    }

    bool firstCls{true}, nominalCompatible{true};
    TimeEstBC nominalTS, expandedTS;
    float smallestROFHalf = std::numeric_limits<float>::max();
    for (int iLayer{0}; iLayer < activeSurfaceCount; ++iLayer) {
      if (track.getClusterIndex(iLayer) == o2::its::constants::UnusedIndex) {
        continue;
      }
      smallestROFHalf = std::min(smallestROFHalf, mScratch->getROFTiming(iLayer).mROFLength * 0.5f);
      mScratch->markUsedCluster(iLayer, track.getClusterIndex(iLayer));
      int currentROF = mScratch->getClusterROF(iLayer, track.getClusterIndex(iLayer));
      const auto nominalROFTS = mScratch->getROFTiming(iLayer).getROFTimeBounds(currentROF);
      const auto expandedROFTS = mScratch->getROFTiming(iLayer).getROFTimeBounds(currentROF, true);
      if (firstCls) {
        firstCls = false;
        nominalTS = nominalROFTS;
        expandedTS = expandedROFTS;
      } else {
        if (nominalCompatible) {
          if (nominalTS.isCompatible(nominalROFTS)) {
            nominalTS += nominalROFTS;
          } else {
            nominalCompatible = false;
          }
        }
        if (!expandedTS.isCompatible(expandedROFTS)) {
          LOGP(fatal, "TS {}+/-{} are incompatible with {}+/-{}, this should not happen!", expandedROFTS.getTimeStamp(), expandedROFTS.getTimeStampError(), expandedTS.getTimeStamp(), expandedTS.getTimeStampError());
        }
        expandedTS += expandedROFTS;
      }
    }
    const auto selectedTimestamp = nominalCompatible ? nominalTS : expandedTS;
    const auto selectedTimestampSymmetric = selectedTimestamp.makeSymmetrical();
    // This is the same sanity clamp as the legacy symmetric timestamp, but
    // committed directly to the detector-neutral GenericTrack interval.
    const float selectedTimestampError = std::min(selectedTimestampSymmetric.getTimeStampError(), smallestROFHalf);
    track.track.timestamp = {static_cast<TFBC>(selectedTimestampSymmetric.getTimeStamp() - selectedTimestampError),
                             static_cast<TFBC>(selectedTimestampSymmetric.getTimeStamp() + selectedTimestampError)};
    trks.emplace_back(track);

    // acceptTracks() is the serial owner-thread publication boundary. Typed
    // sidecars are completed later by the application adapter.
    const auto genericTrackIndex = appendGenericTrack(*mFrame, track, mLayerGlobalMeasurements);
    if (!genericTrackIndex) {
      LOGP(fatal, "GenericTrack publication failed for an accepted CA track");
    }
    trks.back().genericTrackIndex = *genericTrackIndex;

    if (mTrkParams[iteration].AllowSharingFirstCluster) {
      firstClusters[firstLayer].push_back(firstCluster);
    }
  }
}

void TrackerTraits::setNThreads(int n, std::shared_ptr<tbb::task_arena>& arena)
{
#if defined(OPTIMISATION_OUTPUT)
  mTaskArena = std::make_shared<tbb::task_arena>(1);
#else
  if (arena == nullptr) {
    mTaskArena = std::make_shared<tbb::task_arena>(std::abs(n));
    LOGP(info, "Setting tracker with {} threads.", n);
  } else {
    mTaskArena = arena;
  }
#endif
}

} // namespace o2::itsmft::tracking

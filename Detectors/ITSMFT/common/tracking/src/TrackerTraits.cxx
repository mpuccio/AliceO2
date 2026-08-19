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
#include "ITStracking/MathUtils.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IndexTableConfiguration.h"
#include "ITSMFTTracking/Propagator.h"
#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/detail/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/LayerMask.h"
#include "ITSMFTTracking/TripletFitting.h"
#include "ITSMFTTracking/detail/TimeFrameScratch.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/CandidateFinding.h"
#include "ITSMFTTracking/detail/DirectionCompatibility.h"
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

bool appendGenericTrack(TimeFrame& frame,
                        const TrackingCandidate& candidate,
                        gsl::span<const gsl::span<const GlobalMeasurement>> layerMeasurements,
                        gsl::span<const LayerId> orderedSurfaces)
{
  GenericTrack track = candidate.track;
  track.hitLayers = {};
  std::vector<TrackClusterReference> resolvedReferences;
  resolvedReferences.reserve(layerMeasurements.size());
  for (std::size_t position = 0; position < layerMeasurements.size(); ++position) {
    const int localIndex = candidate.getClusterIndex(static_cast<int>(position));
    if (localIndex == o2::its::constants::UnusedIndex) {
      continue;
    }
    if (localIndex < 0 || static_cast<std::size_t>(localIndex) >= layerMeasurements[position].size()) {
      return false;
    }
    if (position >= orderedSurfaces.size()) {
      return false;
    }
    const auto& measurement = layerMeasurements[position][localIndex];
    const TrackClusterReference reference{orderedSurfaces[position], 0, measurement.clusterId};
    if (!reference.isValid()) {
      return false;
    }
    resolvedReferences.push_back(reference);
    track.hitLayers.set(static_cast<int>(position));
  }
  if (!track.innerState.hasRecognizedKind() || !track.outerState.hasRecognizedKind() ||
      !track.timestamp.isValid() || resolvedReferences.empty()) {
    return false;
  }

  auto& tracks = frame.getGenericTracks();
  auto& references = frame.getTrackClusterIndices();
  const auto oldTrackSize = tracks.size();
  const auto oldReferenceSize = references.size();
  if (oldTrackSize > std::numeric_limits<uint32_t>::max() || oldReferenceSize > std::numeric_limits<uint32_t>::max() ||
      resolvedReferences.size() > std::numeric_limits<uint32_t>::max() - oldReferenceSize) {
    return false;
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
  return true;
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

void TrackerTraits::runTraversal(IterationContext view, SeedRefitFunction refitFunction)
{
  if (view.iteration < 0) {
    throw TraversalException{view.iteration, TraversalFailureReason::IterationOutOfRange};
  }
  int maxNvertices{-1};
  if (view.configuration.parameters.PerPrimaryVertexProcessing) {
    maxNvertices = view.frame.getMaxVerticesPerROF();
  }
  int iVertex = std::min(maxNvertices, 0);
  do {
    computeLayerTracklets(view, view.iteration, iVertex);
    computeLayerCells(view, view.iteration);
    findCellsNeighbours(view, view.iteration);
    findRoads(view, view.iteration, refitFunction);
  } while (++iVertex < maxNvertices);
}

int requireScratchEdgeSlot(const IterationContext& context, int iteration, EdgeId id)
{
  const auto slot = context.configuration.getEdgeSlot(id);
  if (!slot) {
    throw TraversalException{iteration, TraversalFailureReason::TraversalBindingMismatch};
  }
  return static_cast<int>(*slot);
}

int requireScratchCellSlot(const IterationContext& context, int iteration, CellPathId id)
{
  const auto slot = context.configuration.getCellSlot(id);
  if (!slot) {
    throw TraversalException{iteration, TraversalFailureReason::TraversalBindingMismatch};
  }
  return static_cast<int>(*slot);
}

int requireSurfacePosition(const IterationContext& context, int iteration, LayerId id)
{
  if (!id.isValid()) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  const auto position = context.configuration.getSurfaceSlot(id);
  if (!position || *position >= context.configuration.topology.orderedSurfaces.size()) {
    throw TraversalException{iteration, TraversalFailureReason::TraversalBindingMismatch};
  }
  return static_cast<int>(*position);
}

void TrackerTraits::computeLayerTracklets(IterationContext& context, const int iteration, int iVertex)
{
  auto& scratch = context.scratch;
  const auto scratchEdgeCount = scratch.getTracklets().size();
  for (size_t edgeId = 0; edgeId < scratchEdgeCount; ++edgeId) {
    scratch.getTracklets()[edgeId].clear();
    scratch.getTrackletsLabel(edgeId).clear();
    std::fill(scratch.getTrackletsLookupTable()[edgeId].begin(), scratch.getTrackletsLookupTable()[edgeId].end(), 0);
  }

  computeLayerTrackletsImpl(context, iteration, iVertex, context.configuration.edges);
}

void TrackerTraits::computeLayerTrackletsImpl(
  IterationContext& context,
  const int iteration,
  const int iVertex,
  gsl::span<const EdgeId> edgeIds)
{
  auto* mScratch = &context.scratch;
  const auto& mMemoryPool = mScratch->getMemoryPool();
  auto* mFrame = &context.frame;
  const auto& trkParam = context.configuration.parameters;
  const auto mBz = context.bz;
  const auto& mTraversalGraph = context.topology;
  const auto& mKernelParameters = context.configuration.kernelParameters;
  const auto& mLayerGlobalMeasurements = context.layerGlobalMeasurements;
  const auto& orderedSurfaces = context.configuration.topology.orderedSurfaces;
  const auto& topology = mTraversalGraph;
  const Vertex diamondVert(trkParam.Diamond, trkParam.DiamondCov, 1, 1.f);

  mTaskArena->execute([&] {
    auto resolveEdgeLayers = [&](int edgeId) -> std::pair<int, int> {
      const auto& edge = topology.getEdge(EdgeId{static_cast<uint16_t>(edgeId)});
      return {requireSurfacePosition(context, iteration, edge.from),
              requireSurfacePosition(context, iteration, edge.to)};
    };

    auto forTracklets = [&](auto Mode, int edgeId, int fromLayer, int toLayer, SurfaceKind kind,
                            const TrackletProjectionCache& edgeCache, int pivotROF, int base, int& offset) -> int {
      if (!mFrame->isROFEnabled(fromLayer, pivotROF)) {
        return 0;
      }
      // Derive a diamond vertex for this pivot ROF; each invocation owns its
      // stack frame, so this is safe inside the parallel dispatch.
      Vertex diamondForROF{};
      gsl::span<const Vertex> primaryVertices;
      if (trkParam.UseDiamond) {
        diamondForROF = diamondVertexForROF(diamondVert, mFrame->getROFViews(fromLayer).overlap,
                                            mFrame->getROFLocalLayer(fromLayer), pivotROF);
        primaryVertices = gsl::span<const Vertex>(&diamondForROF, 1);
      } else {
        primaryVertices = mFrame->getPrimaryVertices(fromLayer, pivotROF);
      }
      if (primaryVertices.empty()) {
        return 0;
      }
      const int startVtx = iVertex >= 0 ? iVertex : 0;
      const int endVtx = iVertex >= 0 ? o2::gpu::CAMath::Min(iVertex + 1, int(primaryVertices.size())) : int(primaryVertices.size());
      if (endVtx <= startVtx || (iVertex + 1) > primaryVertices.size()) {
        return 0;
      }

      const auto& rofOverlap = mFrame->getROFOverlap(fromLayer, toLayer, pivotROF);
      if (!rofOverlap.getEntries()) {
        return 0;
      }

      int localCount = 0;
      auto& tracklets = mScratch->getTracklets()[edgeId];
      auto layer0 = mFrame->getClustersOnLayer(pivotROF, fromLayer);
      if (layer0.empty()) {
        return 0;
      }

      for (int iCluster = 0; iCluster < int(layer0.size()); ++iCluster) {
        const GlobalMeasurement& sourceMeasurement = layer0[iCluster];
        const int currentSortedIndex = mFrame->getSortedIndex(pivotROF, fromLayer, iCluster);
        if (mFrame->isClusterUsed(fromLayer, sourceMeasurement.clusterId)) {
          continue;
        }

        for (int iV = startVtx; iV < endVtx; ++iV) {
          const auto& pv = primaryVertices[iV];
          if (!mFrame->isVertexCompatible(fromLayer, pivotROF, pv)) {
            continue;
          }
          if (pv.isFlagSet(Vertex::Flags::UPCMode) != trkParam.PassFlags[IterationStep::SelectUPCVertices]) {
            continue;
          }
          TrackletSearchWindow searchWindow{};
          const bool projected = projectTrackletSearchWindow(sourceMeasurement, pv, mFrame->getBeamX(), mFrame->getBeamY(),
                                                             mFrame->getBeamPositionVariance(),
                                                             kind, edgeCache, mBz,
                                                             mFrame->getIndexTableUtils(toLayer),
                                                             mKernelParameters, searchWindow);
          if (!projected) {
            continue;
          }
          const auto bins = searchWindow.bins;
          const auto& indexTableUtils = mFrame->getIndexTableUtils(toLayer);
          int rowBinsNum = bins.w - bins.y + 1;
          const bool periodicPhi = indexTableUtils.getCoordType() == IndexTableCoordType::PhiZ ||
                                   indexTableUtils.getCoordType() == IndexTableCoordType::PhiR;
          if (periodicPhi && rowBinsNum < 0) {
            rowBinsNum += indexTableUtils.getNrowBins();
          }
          rowBinsNum = std::max(0, rowBinsNum);

          for (int targetROF = rofOverlap.getFirstEntry(); targetROF < rofOverlap.getEntriesBound(); ++targetROF) {
            if (!mFrame->isROFEnabled(toLayer, targetROF)) {
              continue;
            }
            auto layer1 = mFrame->getClustersOnLayer(targetROF, toLayer);
            if (layer1.empty()) {
              continue;
            }
            const auto ts = mFrame->getROFTimeStamp(fromLayer, pivotROF, toLayer, targetROF);
            if (!ts.isCompatible(pv.getTimeStamp())) {
              continue;
            }
            const auto& targetIndexTable = mFrame->getIndexTable(targetROF, toLayer);
            const int colBinRange = (bins.z - bins.x) + 1;
            for (int iRow = 0; iRow < rowBinsNum; ++iRow) {
              int iRowBin = bins.y + iRow;
              if (periodicPhi) {
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
                const GlobalMeasurement& targetMeasurement = layer1[iNext];
                if (mFrame->isClusterUsed(toLayer, targetMeasurement.clusterId)) {
                  continue;
                }

                float tanL = 0.f;
                const bool accepted = acceptTrackletCandidate(searchWindow, sourceMeasurement,
                                                              targetMeasurement, kind,
                                                              mKernelParameters.nSigmaCut, tanL);
                if (accepted) {
                  const float phi{o2::gpu::GPUCommonMath::ATan2(sourceMeasurement.y - targetMeasurement.y,
                                                                sourceMeasurement.x - targetMeasurement.x)};
                  if constexpr (decltype(Mode)::value == PassMode::OnePass::value) {
                    tracklets.emplace_back(currentSortedIndex, mFrame->getSortedIndex(targetROF, toLayer, iNext), tanL, phi, ts);
                  } else if constexpr (decltype(Mode)::value == PassMode::TwoPassCount::value) {
                    ++localCount;
                  } else if constexpr (decltype(Mode)::value == PassMode::TwoPassInsert::value) {
                    const int idx = base + offset++;
                    tracklets[idx] = Tracklet(currentSortedIndex, mFrame->getSortedIndex(targetROF, toLayer, iNext), tanL, phi, ts);
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
      for (const auto typedEdgeId : edgeIds) {
        const int edgeId = typedEdgeId.value();
        const int scratchEdgeId = requireScratchEdgeSlot(context, iteration, typedEdgeId);
        const auto [fromLayer, toLayer] = resolveEdgeLayers(edgeId);
        const auto kind = topology.getSurface(topology.getEdge(typedEdgeId).from).kind;
        TrackletProjectionCache edgeCache{};
        const auto layerRadii = gsl::span<const float>{context.detectorConfiguration.layerRadii};
        if (!bindTrackletProjectionCache(fromLayer, toLayer, layerRadii,
                                         mFrame->getMinR(toLayer), mFrame->getMaxR(toLayer),
                                         mFrame->getMinZ(toLayer), mFrame->getMaxZ(toLayer),
                                         context.detectorConfiguration.positionResolutions[fromLayer],
                                         mScratch->getEdgeMSAngle(scratchEdgeId),
                                         mScratch->getEdgePhiCut(scratchEdgeId), edgeCache)) {
          throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
        }
        const int startROF = 0, endROF = mFrame->getROFTiming(fromLayer).mNROFsTF;
        for (int pivotROF{startROF}; pivotROF < endROF; ++pivotROF) {
          forTracklets(PassMode::OnePass{}, scratchEdgeId, fromLayer, toLayer, kind, edgeCache, pivotROF, 0, dummy);
        }
      }
    } else {
      tbb::parallel_for(0, static_cast<int>(edgeIds.size()), [&](const int edgeIndex) {
        const auto typedEdgeId = edgeIds[edgeIndex];
        const int edgeId = typedEdgeId.value();
        const int scratchEdgeId = requireScratchEdgeSlot(context, iteration, typedEdgeId);
        const auto [fromLayer, toLayer] = resolveEdgeLayers(edgeId);
        const auto kind = topology.getSurface(topology.getEdge(typedEdgeId).from).kind;
        TrackletProjectionCache edgeCache{};
        const auto layerRadii = gsl::span<const float>{context.detectorConfiguration.layerRadii};
        if (!bindTrackletProjectionCache(fromLayer, toLayer, layerRadii,
                                         mFrame->getMinR(toLayer), mFrame->getMaxR(toLayer),
                                         mFrame->getMinZ(toLayer), mFrame->getMaxZ(toLayer),
                                         context.detectorConfiguration.positionResolutions[fromLayer],
                                         mScratch->getEdgeMSAngle(scratchEdgeId),
                                         mScratch->getEdgePhiCut(scratchEdgeId), edgeCache)) {
          throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
        }
        const int startROF = 0, endROF = mFrame->getROFTiming(fromLayer).mNROFsTF;
        bounded_vector<int> perROFCount((endROF - startROF) + 1, mMemoryPool.get());
        tbb::parallel_for(startROF, endROF, [&](const int pivotROF) {
          perROFCount[pivotROF - startROF] = forTracklets(PassMode::TwoPassCount{}, scratchEdgeId, fromLayer, toLayer, kind, edgeCache, pivotROF, 0, dummy);
        });
        std::exclusive_scan(perROFCount.begin(), perROFCount.end(), perROFCount.begin(), 0);
        const int nTracklets = perROFCount.back();
        mScratch->getTracklets()[scratchEdgeId].resize(nTracklets);
        if (nTracklets == 0) {
          return;
        }
        tbb::parallel_for(startROF, endROF, [&](const int pivotROF) {
          int baseIdx = perROFCount[pivotROF - startROF];
          if (baseIdx == perROFCount[pivotROF + 1 - startROF]) {
            return;
          }
          int localIdx = 0;
          forTracklets(PassMode::TwoPassInsert{}, scratchEdgeId, fromLayer, toLayer, kind, edgeCache, pivotROF, baseIdx, localIdx);
        });
      });
    }

    tbb::parallel_for(0, static_cast<int>(edgeIds.size()), [&](const int edgeIndex) {
      const int scratchEdgeId = requireScratchEdgeSlot(context, iteration, edgeIds[edgeIndex]);
      /// Sort tracklets & remove duplicates
      // duplicates can exist simply since we evaluate per vertex
      auto& trkl{mScratch->getTracklets()[scratchEdgeId]};
      std::sort(trkl.begin(), trkl.end());
      trkl.erase(std::unique(trkl.begin(), trkl.end()), trkl.end());
      trkl.shrink_to_fit();
      auto& lut{mScratch->getTrackletsLookupTable()[scratchEdgeId]};
      if (!trkl.empty()) {
        for (const auto& tkl : trkl) {
          lut[tkl.firstClusterIndex + 1]++;
        }
        std::inclusive_scan(lut.begin(), lut.end(), lut.begin());
      }
    });

    /// Create tracklets labels
    if (mFrame->hasMCinformation() && trkParam.CreateArtefactLabels) {
      tbb::parallel_for(0, static_cast<int>(edgeIds.size()), [&](const int edgeIndex) {
        const auto typedEdgeId = edgeIds[edgeIndex];
        const int edgeId = typedEdgeId.value();
        const int scratchEdgeId = requireScratchEdgeSlot(context, iteration, typedEdgeId);
        const auto [fromLayer, toLayer] = resolveEdgeLayers(edgeId);
        for (auto& trk : mScratch->getTracklets()[scratchEdgeId]) {
          MCCompLabel label;
          const auto currentId = mFrame->getClusters()[fromLayer][trk.firstClusterIndex].clusterId;
          const auto nextId = mFrame->getClusters()[toLayer][trk.secondClusterIndex].clusterId;
          for (const auto& lab1 : mFrame->getLabels(orderedSurfaces[fromLayer], currentId)) {
            for (const auto& lab2 : mFrame->getLabels(orderedSurfaces[toLayer], nextId)) {
              if (lab1 == lab2 && lab1.isValid()) {
                label = lab1;
                break;
              }
            }
            if (label.isValid()) {
              break;
            }
          }
          mScratch->getTrackletsLabel(scratchEdgeId).emplace_back(label);
        }
      });
    }
  });
}

void TrackerTraits::computeLayerCells(IterationContext& context, const int iteration)
{
  auto& scratch = context.scratch;
  const auto scratchCellCount = scratch.getCells().size();
  if (scratch.getCellsLookupTable().size() != scratchCellCount) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  for (size_t cellPathId = 0; cellPathId < scratchCellCount; ++cellPathId) {
    deepVectorClear(scratch.getCells()[cellPathId]);
    deepVectorClear(scratch.getCellsLookupTable()[cellPathId]);
    if (context.frame.hasMCinformation() && context.configuration.parameters.CreateArtefactLabels) {
      deepVectorClear(scratch.getCellsLabel(cellPathId));
    }
  }

  if (!bindAttachHitConfig(context.detectorConfiguration.layerMaterial, context.configuration.parameters)
         .isValid(static_cast<int>(context.configuration.topology.orderedSurfaces.size()))) {
    throw TraversalException{iteration, TraversalFailureReason::InvalidSurfaceParameters};
  }

  computeLayerCellsImpl(context, iteration, context.configuration.cells);

  const auto scratchEdgeCount = scratch.getTracklets().size();
  for (size_t edgeId = 0; edgeId < scratchEdgeCount; ++edgeId) {
    deepVectorClear(scratch.getTracklets()[edgeId]);
    deepVectorClear(scratch.getTrackletsLabel(edgeId));
  }
}

void TrackerTraits::computeLayerCellsImpl(
  IterationContext& context,
  const int iteration,
  gsl::span<const CellPathId> cellIds)
{
  auto* mScratch = &context.scratch;
  const auto& mMemoryPool = mScratch->getMemoryPool();
  const auto& trkParam = context.configuration.parameters;
  const auto mBz = context.bz;
  const auto& mTraversalGraph = context.topology;
  const auto& mKernelParameters = context.configuration.kernelParameters;
  const auto mAttachHitConfig = bindAttachHitConfig(context.detectorConfiguration.layerMaterial, trkParam);
  const auto& mLayerMaterial = context.detectorConfiguration.layerMaterial;
  const auto& mLayerGlobalMeasurements = context.layerGlobalMeasurements;
  const auto& topology = mTraversalGraph;

  mTaskArena->execute([&] {
    struct CellHitBinding {
      std::array<int, 3> layers{};
      std::array<SurfaceDescriptor, 3> surfaces{};
    };

    auto resolveCellHitBinding = [&](const auto& cellTopology) -> CellHitBinding {
      const auto& firstEdge = topology.getEdge(cellTopology.first);
      const auto& secondEdge = topology.getEdge(cellTopology.second);
      const std::array<LayerId, 3> surfaces{firstEdge.from, firstEdge.to, secondEdge.to};
      CellHitBinding binding;
      for (int i = 0; i < 3; ++i) {
        const auto LayerId = surfaces[i];
        binding.layers[i] = requireSurfacePosition(context, iteration, LayerId);
        binding.surfaces[i] = topology.getSurface(LayerId);
      }
      return binding;
    };

    auto forTrackletCells = [&](auto Mode, SurfaceKind kind, int firstEdgeId, int secondEdgeId, const CellHitBinding& hitBinding, bounded_vector<CellSeed>& layerCells, int iTracklet, int offset = 0) -> int {
      const auto& hitLayers = hitBinding.layers;
      const Tracklet& currentTracklet{mScratch->getTracklets()[firstEdgeId][iTracklet]};
      const int nextLayerClusterIndex{currentTracklet.secondClusterIndex};
      const int nextLayerFirstTrackletIndex{mScratch->getTrackletsLookupTable()[secondEdgeId][nextLayerClusterIndex]};
      const int nextLayerLastTrackletIndex{mScratch->getTrackletsLookupTable()[secondEdgeId][nextLayerClusterIndex + 1]};
      int foundCells{0};
      for (int iNextTracklet{nextLayerFirstTrackletIndex}; iNextTracklet < nextLayerLastTrackletIndex; ++iNextTracklet) {
        const Tracklet& nextTracklet{mScratch->getTracklets()[secondEdgeId][iNextTracklet]};
        if (nextTracklet.firstClusterIndex != nextLayerClusterIndex) {
          break;
        }
        if (!currentTracklet.getTimeStamp().isCompatible(nextTracklet.getTimeStamp())) {
          continue;
        }

        /// Prepare the track seed; clusters are numbered from inner to outer.
        const int sortedId[3]{currentTracklet.firstClusterIndex, nextTracklet.firstClusterIndex, nextTracklet.secondClusterIndex};
        const auto& globalInner = mLayerGlobalMeasurements[hitLayers[0]][sortedId[0]];
        const auto& globalMiddle = mLayerGlobalMeasurements[hitLayers[1]][sortedId[1]];
        const auto& globalOuter = mLayerGlobalMeasurements[hitLayers[2]][sortedId[2]];
        const auto* measurementInner = context.frame.getSurfaceMeasurement(hitBinding.surfaces[0].id, globalInner.clusterId);
        const auto* measurementMiddle = context.frame.getSurfaceMeasurement(hitBinding.surfaces[1].id, globalMiddle.clusterId);
        const auto* measurementOuter = context.frame.getSurfaceMeasurement(hitBinding.surfaces[2].id, globalOuter.clusterId);
        if (measurementInner == nullptr || measurementMiddle == nullptr || measurementOuter == nullptr) {
          continue;
        }
        const double edgeMSAngle = static_cast<double>(mScratch->getEdgeMSAngle(secondEdgeId));
        const DirectionProcessNoise directionProcessNoise{edgeMSAngle * edgeMSAngle};
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
        if (cellDirectionsAreCompatible(directionObservations, directionProcessNoise, context.frame.getBeamPositionVariance(), mKernelParameters.nSigmaCut,
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
          const bool good = buildCellSeed(kind, globalInner, globalMiddle,
                                          *measurementInner, *measurementMiddle, *measurementOuter,
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
              layerCells.emplace_back(hitLayerMask, sortedId[0], sortedId[1], sortedId[2], iTracklet, iNextTracklet, state, chi2, ts);
              layerCells.back().tripletFactor() = tripletFactor;
              ++foundCells;
            } else if constexpr (decltype(Mode)::value == PassMode::TwoPassCount::value) {
              ++foundCells;
            } else if constexpr (decltype(Mode)::value == PassMode::TwoPassInsert::value) {
              layerCells[offset++] = CellSeed(hitLayerMask, sortedId[0], sortedId[1], sortedId[2], iTracklet, iNextTracklet, state, chi2, ts);
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
      const int cellPathId = requireScratchCellSlot(context, iteration, typedCellId);
      const auto& cellTopology = topology.getPath(typedCellId);
      const auto kind = topology.getSurface(topology.getEdge(cellTopology.first).from).kind;
      const int firstEdgeId = requireScratchEdgeSlot(context, iteration, cellTopology.first);
      const int secondEdgeId = requireScratchEdgeSlot(context, iteration, cellTopology.second);
      if (mScratch->getTracklets()[firstEdgeId].empty() ||
          mScratch->getTracklets()[secondEdgeId].empty()) {
        continue;
      }
      const auto hitBinding = resolveCellHitBinding(cellTopology);

      auto& layerCells = mScratch->getCells()[cellPathId];
      const int currentLayerTrackletsNum{static_cast<int>(mScratch->getTracklets()[firstEdgeId].size())};
      bounded_vector<int> perTrackletCount(currentLayerTrackletsNum + 1, 0, mMemoryPool.get());
      if (mTaskArena->max_concurrency() <= 1) {
        for (int iTracklet{0}; iTracklet < currentLayerTrackletsNum; ++iTracklet) {
          perTrackletCount[iTracklet] = forTrackletCells(PassMode::OnePass{}, kind, firstEdgeId, secondEdgeId, hitBinding, layerCells, iTracklet);
        }
        std::exclusive_scan(perTrackletCount.begin(), perTrackletCount.end(), perTrackletCount.begin(), 0);
      } else {
        tbb::parallel_for(0, currentLayerTrackletsNum, [&](const int iTracklet) {
          perTrackletCount[iTracklet] = forTrackletCells(PassMode::TwoPassCount{}, kind, firstEdgeId, secondEdgeId, hitBinding, layerCells, iTracklet);
        });

        std::exclusive_scan(perTrackletCount.begin(), perTrackletCount.end(), perTrackletCount.begin(), 0);
        auto totalCells{perTrackletCount.back()};
        if (totalCells == 0) {
          auto& lut = mScratch->getCellsLookupTable()[cellPathId];
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
          forTrackletCells(PassMode::TwoPassInsert{}, kind, firstEdgeId, secondEdgeId, hitBinding, layerCells, iTracklet, offset);
        });
      }

      auto& lut = mScratch->getCellsLookupTable()[cellPathId];
      lut.resize(currentLayerTrackletsNum + 1);
      std::copy_n(perTrackletCount.begin(), currentLayerTrackletsNum + 1, lut.begin());

      if (context.frame.hasMCinformation() && trkParam.CreateArtefactLabels) {
        auto& labels = mScratch->getCellsLabel(cellPathId);
        labels.reserve(layerCells.size());
        for (const auto& cell : layerCells) {
          MCCompLabel currentLab{mScratch->getTrackletsLabel(firstEdgeId)[cell.getFirstTrackletIndex()]};
          MCCompLabel nextLab{mScratch->getTrackletsLabel(secondEdgeId)[cell.getSecondTrackletIndex()]};
          labels.emplace_back(currentLab == nextLab ? currentLab : MCCompLabel());
        }
      }
    }
  });
}

void TrackerTraits::findCellsNeighbours(IterationContext& context, const int iteration)
{
  auto& scratch = context.scratch;
  for (std::size_t slot = 0; slot < scratch.getCellsNeighbours().size(); ++slot) {
    deepVectorClear(scratch.getCellsNeighbours()[slot]);
    deepVectorClear(scratch.getCellsNeighboursTopology()[slot]);
    deepVectorClear(scratch.getCellsNeighboursLUT()[slot]);
  }
  const auto& scheduledCells = context.configuration.topology.scheduledPaths;
  if (!scheduledCells.empty()) {
    findCellsNeighboursForSchedule(context, iteration, scheduledCells, context.configuration.kernelParameters);
  }
  for (auto& cellLUT : scratch.getCellsLookupTable()) {
    deepVectorClear(cellLUT);
  }
}

void TrackerTraits::findCellsNeighboursForSchedule(
  IterationContext& context,
  int iteration,
  gsl::span<const CellPathId> scheduledCells,
  const TrackingKernelParameters& params)
{
  auto* mScratch = &context.scratch;
  const auto& mMemoryPool = mScratch->getMemoryPool();
  const auto& mTraversalGraph = context.topology;
  const auto& mLayerMaterial = context.detectorConfiguration.layerMaterial;
  const auto& mLayerGlobalMeasurements = context.layerGlobalMeasurements;
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
    for (size_t cellPathId = 0; cellPathId < scratchCellCount; ++cellPathId) {
      cellsNeighboursByTarget.emplace_back(mMemoryPool.get());
    }

    for (const auto scheduledId : scheduledCells) {
      const int cellPathId = requireScratchCellSlot(context, iteration, scheduledId);
      if (static_cast<size_t>(cellPathId) >= scratchCellCount ||
          static_cast<size_t>(cellPathId) >= mScratch->getCellsLookupTable().size()) {
        throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
      }
      const auto& cellTopology = topology.getPath(scheduledId);
      const int currentProcessEdge = requireScratchEdgeSlot(context, iteration, cellTopology.second);
      const double currentMSAngle = mScratch->getEdgeMSAngle(currentProcessEdge);
      const double currentAngularVariance = currentMSAngle * currentMSAngle;
      if (mScratch->getCells()[cellPathId].empty()) {
        continue;
      }
      const auto successors = topology.getPathsStartingWithEdge(cellTopology.second);
      if (!successors.getEntries()) {
        continue;
      }

      tbb::enumerable_thread_specific<bounded_vector<CellNeighbour>> sourceNeighbours([&]() { return bounded_vector<CellNeighbour>{mMemoryPool.get()}; });
      tbb::parallel_for(0, static_cast<int>(mScratch->getCells()[cellPathId].size()), [&](const int iCell) {
        auto& localNeighbours = sourceNeighbours.local();
        const auto& currentCellSeed{mScratch->getCells()[cellPathId][iCell]};
        const int nextLayerTrackletIndex{currentCellSeed.getSecondTrackletIndex()};
        for (uint32_t iSuccessor = 0; iSuccessor < successors.getEntries(); ++iSuccessor) {
          // Translate dynamically discovered neighbours from the global CSR
          // topology before accessing scratch.
          const auto nextTopologyId = topology.pathsByFirstEdge[successors.getFirstEntry() + iSuccessor];
          const int nextCellPathId = requireScratchCellSlot(context, iteration, nextTopologyId);
          const auto& nextCellTopology = topology.getPath(nextTopologyId);
          const int nextProcessEdge = requireScratchEdgeSlot(context, iteration, nextCellTopology.second);
          const double nextMSAngle = mScratch->getEdgeMSAngle(nextProcessEdge);
          const double nextAngularVariance = nextMSAngle * nextMSAngle;
          if (static_cast<size_t>(nextCellPathId) >= mScratch->getCells().size() ||
              static_cast<size_t>(nextCellPathId) >= mScratch->getCellsLookupTable().size()) {
            throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
          }
          if (mScratch->getCells()[nextCellPathId].empty() ||
              mScratch->getCellsLookupTable()[nextCellPathId].empty()) {
            continue;
          }
          const auto& nextCellLUT = mScratch->getCellsLookupTable()[nextCellPathId];
          if (nextLayerTrackletIndex < 0 || nextLayerTrackletIndex + 1 >= static_cast<int>(nextCellLUT.size())) {
            continue;
          }
          const int nextLayerFirstCellIndex{nextCellLUT[nextLayerTrackletIndex]};
          const int nextLayerLastCellIndex{nextCellLUT[nextLayerTrackletIndex + 1]};
          if (nextLayerFirstCellIndex < 0 || nextLayerLastCellIndex < nextLayerFirstCellIndex ||
              nextLayerLastCellIndex > static_cast<int>(mScratch->getCells()[nextCellPathId].size())) {
            throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
          }
          for (int iNextCell{nextLayerFirstCellIndex}; iNextCell < nextLayerLastCellIndex; ++iNextCell) {
            const auto& nextCellSeedRef{mScratch->getCells()[nextCellPathId][iNextCell]};
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
            localNeighbours.emplace_back(cellPathId, iCell, nextCellPathId, iNextCell, nextLevel);
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

    for (size_t cellPathId = 0; cellPathId < scratchCellCount; ++cellPathId) {
      auto& cellsNeighbours = cellsNeighboursByTarget[cellPathId];
      if (cellsNeighbours.empty()) {
        continue;
      }

      std::sort(cellsNeighbours.begin(), cellsNeighbours.end(), [](const auto& a, const auto& b) {
        return a.nextCell < b.nextCell;
      });

      auto& cellsNeighbourLUT = mScratch->getCellsNeighboursLUT()[cellPathId];
      cellsNeighbourLUT.assign(mScratch->getCells()[cellPathId].size(), 0);
      for (const auto& neigh : cellsNeighbours) {
        ++cellsNeighbourLUT[neigh.nextCell];
      }
      std::inclusive_scan(cellsNeighbourLUT.begin(), cellsNeighbourLUT.end(), cellsNeighbourLUT.begin());

      mScratch->getCellsNeighbours()[cellPathId].reserve(cellsNeighbours.size());
      mScratch->getCellsNeighboursTopology()[cellPathId].reserve(cellsNeighbours.size());
      std::ranges::transform(cellsNeighbours, std::back_inserter(mScratch->getCellsNeighbours()[cellPathId]), [](const auto& neigh) { return neigh.cell; });
      std::ranges::transform(cellsNeighbours, std::back_inserter(mScratch->getCellsNeighboursTopology()[cellPathId]), [](const auto& neigh) { return neigh.cellTopology; });
    }
  });
}

template <typename InputSeed>
void TrackerTraits::processNeighbours(IterationContext& context, int iteration, int defaultCellPathId, int iLevel, const bounded_vector<InputSeed>& currentCellSeed, const bounded_vector<int>& currentCellId, const bounded_vector<int>& currentCellPathId, bounded_vector<TrackSeed>& updatedCellSeeds, bounded_vector<int>& updatedCellsIds, bounded_vector<int>& updatedCellsPathIds, const TrackingKernelParameters& params)
{
  auto* mScratch = &context.scratch;
  const auto& mMemoryPool = mScratch->getMemoryPool();
  const auto mBz = context.bz;
  const auto mAttachHitConfig = bindAttachHitConfig(context.detectorConfiguration.layerMaterial,
                                                    context.configuration.parameters);
  const auto layerMaterial = mAttachHitConfig.layerMaterial;
  const auto& mLayerGlobalMeasurements = context.layerGlobalMeasurements;
  const auto& orderedSurfaces = context.configuration.topology.orderedSurfaces;
  const int activeSurfaceCount = static_cast<int>(orderedSurfaces.size());

  mTaskArena->execute([&] {
    auto forCellNeighbours = [&](auto Mode, int iCell, int offset = 0) -> int {
      const auto& currentCell{currentCellSeed[iCell]};
      const int cellPathId = currentCellPathId.empty() ? defaultCellPathId : currentCellPathId[iCell];

      if constexpr (decltype(Mode)::value != PassMode::TwoPassInsert::value) {
        if (currentCell.getLevel() != iLevel) {
          return 0;
        }
        if (currentCellId.empty()) {
          for (int layer = 0; layer < activeSurfaceCount; ++layer) {
            const int clusterIndex = currentCell.getCluster(layer);
            if (clusterIndex != o2::its::constants::UnusedIndex &&
                context.frame.isClusterUsed(layer, mLayerGlobalMeasurements[layer][clusterIndex].clusterId)) {
              return 0; /// this we do only on the first iteration, hence the check on currentCellId
            }
          }
        }
      }

      const int cellId = currentCellId.empty() ? iCell : currentCellId[iCell];
      if (cellPathId < 0 || mScratch->getCellsNeighboursLUT()[cellPathId].empty()) {
        return 0;
      }
      const int startNeighbourId{cellId ? mScratch->getCellsNeighboursLUT()[cellPathId][cellId - 1] : 0};
      const int endNeighbourId{mScratch->getCellsNeighboursLUT()[cellPathId][cellId]};
      int foundSeeds{0};
      for (int iNeighbourCell{startNeighbourId}; iNeighbourCell < endNeighbourId; ++iNeighbourCell) {
        const int neighbourCellPathId = mScratch->getCellsNeighboursTopology()[cellPathId][iNeighbourCell];
        const int neighbourCellId = mScratch->getCellsNeighbours()[cellPathId][iNeighbourCell];
        const auto& neighbourCell = mScratch->getCells()[neighbourCellPathId][neighbourCellId];
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
        const auto& neighbourGlobal = mLayerGlobalMeasurements[neighbourLayer][neighbourCluster];
        if (context.frame.isClusterUsed(neighbourLayer, neighbourGlobal.clusterId)) {
          continue;
        }

        /// Let's start the fitting procedure
        TrackSeed seed{currentCell};
        seed.getTimeStamp() = currentCell.getTimeStamp();
        seed.getTimeStamp() += neighbourCell.getTimeStamp();

        const auto* measurement = context.frame.getSurfaceMeasurement(orderedSurfaces[neighbourLayer], neighbourGlobal.clusterId);
        if (measurement == nullptr) {
          continue;
        }
        float chi2 = seed.getChi2();
        OperationFailureReason attachReason{};
        const bool attached = Propagator::attachMeasurement(seed.state(), *measurement, layerMaterial[neighbourLayer], mBz,
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
          auto hitLayerMask = seed.getHitLayerMask();
          hitLayerMask.set(neighbourLayer);
          seed.setHitLayerMask(hitLayerMask);
          seed.setLevel(neighbourCell.getLevel());
          seed.setFirstTrackletIndex(neighbourCell.getFirstTrackletIndex());
          seed.setSecondTrackletIndex(neighbourCell.getSecondTrackletIndex());
        }

        if constexpr (decltype(Mode)::value == PassMode::OnePass::value) {
          updatedCellSeeds.push_back(seed);
          updatedCellsIds.push_back(neighbourCellId);
          updatedCellsPathIds.push_back(neighbourCellPathId);
        } else if constexpr (decltype(Mode)::value == PassMode::TwoPassCount::value) {
          ++foundSeeds;
        } else if constexpr (decltype(Mode)::value == PassMode::TwoPassInsert::value) {
          updatedCellSeeds[offset] = seed;
          updatedCellsIds[offset] = neighbourCellId;
          updatedCellsPathIds[offset++] = neighbourCellPathId;
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
      updatedCellsPathIds.resize(totalNeighbours);

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

void TrackerTraits::findRoads(IterationContext& context, const int iteration, SeedRefitFunction refitFunction)
{
  if (context.scratch.getCells().size() != context.configuration.cells.size()) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  findRoadsImpl(context, iteration, refitFunction);
}

void TrackerTraits::findRoadsImpl(IterationContext& context, const int iteration, SeedRefitFunction refitFunction)
{
  auto* mScratch = &context.scratch;
  const auto& mMemoryPool = mScratch->getMemoryPool();
  const auto& trkParam = context.configuration.parameters;
  const auto mBz = context.bz;
  const auto& mTraversalGraph = context.topology;
  const auto& mKernelParameters = context.configuration.kernelParameters;
  const auto& mLayerGlobalMeasurements = context.layerGlobalMeasurements;
  const gsl::span<const CellPathId> roadStartCells = context.configuration.topology.roadStartPaths;
  const auto& orderedSurfaces = context.configuration.topology.orderedSurfaces;
  const int activeSurfaceCount = static_cast<int>(orderedSurfaces.size());
  bounded_vector<bounded_vector<int>> firstClusters(activeSurfaceCount, bounded_vector<int>(mMemoryPool.get()), mMemoryPool.get());
  firstClusters.resize(activeSurfaceCount);
  // Road starts are the binding's seeding-eligible sparse-plan subsequence.
  // CellPathId values use compact slots; LayerId is never a vector index.
  // Filter roads by absolute q/pT in parameters[4]'s units, identically for
  // both families. Non-finite values fail the finite-bound comparison.
  constexpr float maxAbsQOverPt = 1.e3f;
  const auto seedingLayerMask = context.topology.seedingLayers;
  const auto nonSeedingLayerMask = ~seedingLayerMask;
  const int cellsPerRoad = seedingLayerMask.count() - 2;
  const auto& componentOffsets = context.configuration.topology.roadStartComponentOffsets;
  const auto holeLayerMask = context.frame.getLayout().getHoleLayers();
  if (componentOffsets.empty() || componentOffsets.front() != 0 || componentOffsets.back() != roadStartCells.size()) {
    throw TraversalException{iteration, TraversalFailureReason::SparseTopologyMismatch};
  }
  for (size_t component = 0; component + 1 < componentOffsets.size(); ++component) {
    const auto componentRoadStarts = roadStartCells.subspan(componentOffsets[component],
                                                            componentOffsets[component + 1] - componentOffsets[component]);
    for (int startLevel{cellsPerRoad}; startLevel >= trkParam.CellMinimumLevel(); --startLevel) {

      auto seedFilter = [&](const auto& seed) {
        const auto hitLayerMask = seed.getHitLayerMask();
        const int effectiveTrackLength = hitLayerMask.empty()
                                           ? 0
                                           : hitLayerMask.length() - (LayerMask::span(hitLayerMask.first(), hitLayerMask.last()) & nonSeedingLayerMask).count();
        const auto effectiveHoleMask = hitLayerMask.holeMask() & ~nonSeedingLayerMask;
        return effectiveHoleMask.isAllowedHoleMask(trkParam.MaxHoles, holeLayerMask) &&
               effectiveTrackLength >= trkParam.getMinSeedingClusters() &&
               std::abs(seed.getQOverPt()) <= maxAbsQOverPt && seed.getChi2() <= trkParam.MaxChi2NDF * ((startLevel + 2) * 2 - 5);
      };

      bounded_vector<TrackSeed> trackSeeds(mMemoryPool.get());
      // The binding supplies the ownership-filtered road-start span.
      for (const auto startId : componentRoadStarts) {
        // Translate the road-start cell to its compact slot once.
        const int startCellPathId = requireScratchCellSlot(context, iteration, startId);
        // Cell population is per-event/per-vertex data, so check it against
        // the current vertex rather than caching it in the pass plan.
        if (mScratch->getCells()[startCellPathId].empty()) {
          continue;
        }

        bounded_vector<int> lastCellId(mMemoryPool.get()), updatedCellId(mMemoryPool.get());
        bounded_vector<int> lastCellPathId(mMemoryPool.get()), updatedCellPathId(mMemoryPool.get());
        bounded_vector<TrackSeed> lastCellSeed(mMemoryPool.get()), updatedCellSeed(mMemoryPool.get());

        processNeighbours(context, iteration, startCellPathId, startLevel, mScratch->getCells()[startCellPathId], lastCellId, lastCellPathId, updatedCellSeed, updatedCellId, updatedCellPathId, mKernelParameters);

        int level = startLevel;
        while (level > 2 && !updatedCellSeed.empty()) {
          lastCellSeed.swap(updatedCellSeed);
          lastCellId.swap(updatedCellId);
          lastCellPathId.swap(updatedCellPathId);
          deepVectorClear(updatedCellSeed); /// tame the memory peaks
          deepVectorClear(updatedCellId);   /// tame the memory peaks
          deepVectorClear(updatedCellPathId);
          processNeighbours(context, iteration, o2::its::constants::UnusedIndex, --level, lastCellSeed, lastCellId, lastCellPathId, updatedCellSeed, updatedCellId, updatedCellPathId, mKernelParameters);
        }
        deepVectorClear(lastCellId);     /// tame the memory peaks
        deepVectorClear(lastCellPathId); /// tame the memory peaks
        deepVectorClear(lastCellSeed);   /// tame the memory peaks

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
                                                  context.frame,
                                                  trkParam,
                                                  mBz,
                                                  mLayerGlobalMeasurements,
                                                  mTraversalGraph.getSurfaceCatalogView(),
                                                  orderedSurfaces,
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

void TrackerTraits::acceptTracks(IterationContext& context, int iteration,
                                 bounded_vector<TrackingCandidate>& tracks,
                                 bounded_vector<bounded_vector<int>>& firstClusters)
{
  auto* mScratch = &context.scratch;
  auto* mFrame = &context.frame;
  const auto& trkParam = context.configuration.parameters;
  const auto& mLayerGlobalMeasurements = context.layerGlobalMeasurements;
  const auto& orderedSurfaces = context.configuration.topology.orderedSurfaces;
  const int activeSurfaceCount = static_cast<int>(orderedSurfaces.size());
  for (auto& track : tracks) {
    int nShared = 0;
    bool isFirstShared{false};
    int firstLayer{-1}, firstCluster{-1};
    for (int iLayer{0}; iLayer < activeSurfaceCount; ++iLayer) {
      if (track.getClusterIndex(iLayer) == o2::its::constants::UnusedIndex) {
        continue;
      }
      const auto clusterId = mLayerGlobalMeasurements[iLayer][track.getClusterIndex(iLayer)].clusterId;
      bool isShared = mFrame->isClusterUsed(iLayer, clusterId);
      nShared += int(isShared);
      if (firstLayer < 0) {
        firstCluster = track.getClusterIndex(iLayer);
        isFirstShared = isShared && trkParam.AllowSharingFirstCluster && std::find(firstClusters[iLayer].begin(), firstClusters[iLayer].end(), firstCluster) != firstClusters[iLayer].end();
        firstLayer = iLayer;
      }
    }

    /// do not account for the first cluster in the shared clusters number if it is allowed
    if (nShared - int(isFirstShared && trkParam.AllowSharingFirstCluster) > trkParam.SharedMaxClusters) {
      continue;
    }

    bool firstCls{true}, nominalCompatible{true};
    TimeEstBC nominalTS, expandedTS;
    float smallestROFHalf = std::numeric_limits<float>::max();
    for (int iLayer{0}; iLayer < activeSurfaceCount; ++iLayer) {
      if (track.getClusterIndex(iLayer) == o2::its::constants::UnusedIndex) {
        continue;
      }
      smallestROFHalf = std::min(smallestROFHalf, mFrame->getROFTiming(iLayer).mROFLength * 0.5f);
      const auto clusterId = mLayerGlobalMeasurements[iLayer][track.getClusterIndex(iLayer)].clusterId;
      mFrame->markUsedCluster(iLayer, clusterId);
      int currentROF = mFrame->getClusterROF(iLayer, track.getClusterIndex(iLayer));
      const auto nominalROFTS = mFrame->getROFTiming(iLayer).getROFTimeBounds(currentROF);
      const auto expandedROFTS = mFrame->getROFTiming(iLayer).getROFTimeBounds(currentROF, true);
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
    if (!appendGenericTrack(*mFrame, track, mLayerGlobalMeasurements, orderedSurfaces)) {
      LOGP(fatal, "GenericTrack publication failed for an accepted CA track");
    }

    if (trkParam.AllowSharingFirstCluster) {
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

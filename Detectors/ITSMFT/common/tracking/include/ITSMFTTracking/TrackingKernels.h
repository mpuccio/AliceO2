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

#ifndef ALICEO2_ITSMFT_TRACKINGKERNELS_H_
#define ALICEO2_ITSMFT_TRACKINGKERNELS_H_

#include <cmath>
#include <cstddef>
#include <span>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include "GPUCommonDef.h"
#include "GPUCommonMath.h"
#include "CommonConstants/MathConstants.h"
#include "ITSMFTTracking/BoundedAllocator.h"
#include "ITSMFTTracking/GenericTrack.h"
#include "ITSMFTTracking/TrackingPrimitives.h"
#include "ITSMFTTracking/Triplet.h"
#include "ITSMFTTracking/TrackSeed.h"
#include "ITSMFTTracking/Propagator.h"
#include "ITSMFTTracking/TripletFitting.h"
#include "ITSMFTTracking/RefitDriver.h"
#include "ITSMFTTracking/MathUtils.h"
#include "ITSMFTTracking/ROFViews.h"
#include "ITSMFTTracking/detail/CandidateFinding.h"
#include "ITSMFTTracking/detail/TrackingKernelParameters.h"

namespace o2::itsmft::tracking
{
// The CA steps shared by the CPU and GPU stages: each is written once, for
// one source element, and the stages only differ in how they run it over all
// sources and store what it emits.
struct TrackletTarget {
  float x, y, z, radius, phi;
  bool used;
};
struct TrackletJob {
  float x, y, z;
  float sourceReference, sourceProjected, slope;
  float varianceConstant, varianceLinear, varianceQuadratic;
  float phiPrediction, phiVariance, nSigmaCutSquared;
  bool cylinder;
};

GPUhdi() bool evaluateTracklet(const TrackletJob& job, const TrackletTarget& target, float& tanLambda, float& phi)
{
  if (target.used) {
    return false;
  }
  const float referenceDelta = (job.cylinder ? target.radius : target.z) - job.sourceReference;
  const float prediction = job.sourceProjected + job.slope * referenceDelta;
  const float variance = job.varianceConstant + referenceDelta * (job.varianceLinear + referenceDelta * job.varianceQuadratic);
  const float residual = prediction - (job.cylinder ? target.z : target.radius);
  const float phiResidual = ::remainderf(job.phiPrediction - target.phi, o2::constants::math::TwoPI);
  if (!(variance > 0.f && job.phiVariance > 0.f)) {
    return false;
  }
  // Two separate windows, as the legacy tracker cuts. The azimuthal one is an
  // acceptance bound (the bending at the minimum pT), not a Gaussian
  // resolution, and both residuals grow as 1/p: summing them as independent
  // pulls loses true low-momentum tracklets near the edge of either window.
  if (residual * residual >= job.nSigmaCutSquared * variance ||
      phiResidual * phiResidual >= job.nSigmaCutSquared * job.phiVariance) {
    return false;
  }
  const float chord = ::hypotf(target.x - job.x, target.y - job.y);
  if (!(chord > 1.e-6f)) {
    return false;
  }
  tanLambda = (target.z - job.z) / chord;
  phi = o2::gpu::GPUCommonMath::ATan2(job.y - target.y, job.x - target.x);
  return true;
}

// One edge of the tracklet search. The ranges are indexed by source ROF;
// vertexRanges and the frame's vertices are used without a diamond.
struct TrackletSearch {
  const RuntimeROFTableEntry* overlapRanges{};
  size_t nOverlapRanges{};
  const RuntimeROFTableEntry* vertexRanges{};
  size_t nVertexRanges{};
  ROFTimingLayer fromLayerTiming{};
  ROFTimingLayer toLayerTiming{};
  TrackletProjectionCache edgeCache{};
  o2::itsmft::IndexTableUtilsCore indexTableUtils{};
  float nSigmaCut{};
  float beamPositionVariance{};
  bool useDiamond{false};
  Vertex diamondBase{};
  bool selectUPC{false};
  int iVertex{-1};
  SurfaceKind kind{SurfaceKind::Cylinder};
};

// Every tracklet of one source cluster: over its compatible vertices, the
// overlapping target ROFs and the target bins of its search window.
template <typename Emit>
GPUhdi() void forEachTracklet(const FrameView& frame, const TrackletSearch& edge, int source, Emit emit)
{
  const auto& from = frame.layers[edge.edgeCache.fromLayer];
  const auto& to = frame.layers[edge.edgeCache.toLayer];
  const auto& cluster = from.clusters[source];
  const int rof = cluster.rof;
  const auto overlap = edge.overlapRanges[rof];
  if (!from.rofEnabled[rof] || !overlap.getEntries() || from.isUsed(cluster.clusterId)) {
    return;
  }
  const auto vertices = edge.useDiamond ? RuntimeROFTableEntry{0, 1} : edge.vertexRanges[rof];
  const int nVertices = vertices.getEntries();
  const int firstVertex = edge.iVertex >= 0 ? edge.iVertex : 0;
  const int endVertex = edge.iVertex >= 0 ? o2::gpu::CAMath::Min(edge.iVertex + 1, nVertices) : nVertices;
  const auto& utils = edge.indexTableUtils;
  for (int iVertex = firstVertex; iVertex < endVertex; ++iVertex) {
    Vertex vertex = edge.useDiamond ? edge.diamondBase : frame.vertices[vertices.getFirstEntry() + iVertex];
    if (edge.useDiamond) {
      vertex.setTimeStamp(edge.fromLayerTiming.getROFTimeBounds(rof, true));
    }
    TrackletSearchWindow window{};
    if (!isVertexCompatibleWithROFTiming(edge.fromLayerTiming, rof, vertex) || vertex.isFlagSet(Vertex::Flags::UPCMode) != edge.selectUPC ||
        !projectTrackletSearchWindow(cluster, vertex, edge.beamPositionVariance, edge.kind, edge.edgeCache, utils, edge.nSigmaCut, window)) {
      continue;
    }
    int rowBins = window.bins.w - window.bins.y + 1;
    if (rowBins < 0) {
      rowBins += utils.getNrowBins();
    }
    const TrackletJob job{cluster.x, cluster.y, cluster.z, window.sourceReferenceCoordinate, window.sourceProjectedCoordinate, window.slope,
                          window.varianceConstant, window.varianceLinear, window.varianceQuadratic, window.phiPrediction, window.phiVariance,
                          o2::its::math_utils::Sq(edge.nSigmaCut), edge.kind == SurfaceKind::Cylinder};
    for (int targetROF = overlap.getFirstEntry(); targetROF < overlap.getEntriesBound(); ++targetROF) {
      const auto time = edge.fromLayerTiming.getROFTimeBounds(rof, true) + edge.toLayerTiming.getROFTimeBounds(targetROF, true);
      if (!to.rofEnabled[targetROF] || !time.isCompatible(vertex.getTimeStamp())) {
        continue;
      }
      const int* table = to.indexTables + size_t(targetROF) * to.tableSize;
      const int firstTarget = to.rofClusters[targetROF];
      const int nTargets = to.rofClusters[targetROF + 1] - firstTarget;
      for (int row = 0; row < rowBins; ++row) {
        const int rowBin = (window.bins.y + row) % utils.getNrowBins();
        if (rowBin < 0) {
          break;
        }
        const int firstBin = utils.getBinIndex(window.bins.x, rowBin);
        const int endRow = o2::gpu::CAMath::Min(table[firstBin + window.bins.z - window.bins.x + 1], nTargets);
        for (int target = firstTarget + table[firstBin]; target < firstTarget + endRow; ++target) {
          const auto& next = to.clusters[target];
          float tanLambda, phi;
          if (evaluateTracklet(job, {next.x, next.y, next.z, next.radius, next.phi, to.isUsed(next.clusterId)}, tanLambda, phi)) {
            emit(Tracklet{source, target, tanLambda, phi, time});
          }
        }
      }
    }
  }
}

struct CellParameters {
  float angularTolerance;
  float maximumCurvature;
  bool disk;
};

GPUhdi() bool evaluateCell(const CellParameters& params, const Tracklet& first, const Tracklet& second,
                          const GlobalMeasurement& inner, const GlobalMeasurement& middle,
                          const GlobalMeasurement& outer, TripletFitFactor& factor)
{
  if (first.secondClusterIndex != second.firstClusterIndex || !first.isCompatible(second)) {
    return false;
  }
  const float lambda01 = std::atan(first.tanLambda);
  const float lambda12 = std::atan(second.tanLambda);
  const float sinTheta = std::max(std::abs(o2::its::math_utils::cosFloat(0.5f * (lambda01 + lambda12))), float(o2::constants::math::Almost0));
  // Disk scattering is expressed using pT_min as p: project with sin(theta).
  const float dipTolerance = params.disk ? params.angularTolerance * sinTheta : params.angularTolerance;
  if (std::abs(lambda01 - lambda12) > dipTolerance) {
    return false;
  }
  const float length01 = std::hypot(inner.x - middle.x, inner.y - middle.y);
  const float length12 = std::hypot(middle.x - outer.x, middle.y - outer.y);
  const float curvature = std::min(std::min(params.maximumCurvature, 2.f / length01), 2.f / length12);
  const float maximumBending = std::asin(std::clamp(0.5f * curvature * length01, 0.f, 1.f)) +
                               std::asin(std::clamp(0.5f * curvature * length12, 0.f, 1.f));
  const float deltaPhi = std::abs(std::remainder(first.phi - second.phi, o2::constants::math::TwoPI));
  // The same momentum correction cancels in the disk azimuthal projection.
  const float azimuthalTolerance = params.disk ? params.angularTolerance : params.angularTolerance / sinTheta;
  if (deltaPhi > maximumBending + azimuthalTolerance) {
    return false;
  }
  return makeTripletFitFactor({inner, middle, outer}, factor);
}

// One cell path: its two edges' tracklets, the second edge's lookup table
// (first tracklet of every middle cluster, plus the total) and its layers.
struct CellSearch {
  const Tracklet* first{};
  const Tracklet* second{};
  const int* secondLookup{};
  std::array<int, 3> layers{};
  CellParameters parameters{};
};

// Every cell of one first-edge tracklet, in second-tracklet order.
template <typename Emit>
GPUhdi() void forEachCell(const FrameView& frame, const CellSearch& path, int tracklet, Emit emit)
{
  const auto& first = path.first[tracklet];
  const int middle = first.secondClusterIndex;
  for (int next = path.secondLookup[middle]; next < path.secondLookup[middle + 1]; ++next) {
    const auto& second = path.second[next];
    if (second.firstClusterIndex != middle) {
      break;
    }
    TripletFitFactor factor{};
    if (evaluateCell(path.parameters, first, second, frame.layers[path.layers[0]].clusters[first.firstClusterIndex],
                     frame.layers[path.layers[1]].clusters[middle], frame.layers[path.layers[2]].clusters[second.secondClusterIndex], factor)) {
      Triplet cell{LayerMask{path.layers[0], path.layers[1], path.layers[2]}, first.firstClusterIndex, middle, second.secondClusterIndex,
                   tracklet, next, first.getTimeStamp() + second.getTimeStamp()};
      cell.tripletFactor() = factor;
      emit(cell);
    }
  }
}

inline bool cellFactorsEquivalent(const TripletFitFactor& a, const TripletFitFactor& b)
{
  const auto close = [](float x, float y) {
    return std::isfinite(x) && std::isfinite(y) && std::abs(x - y) <= 2.e-5f * std::max({1.f, std::abs(x), std::abs(y)});
  };
  if (!close(a.psi.theta, b.psi.theta) || !close(a.psi.phi, b.psi.phi) ||
      !close(a.rho.theta, b.rho.theta) || !close(a.rho.phi, b.rho.phi)) {
    return false;
  }
  for (int hit = 0; hit < 3; ++hit) {
    for (int coordinate = 0; coordinate < 3; ++coordinate) {
      if (!close(a.h[hit].theta[coordinate], b.h[hit].theta[coordinate]) ||
          !close(a.h[hit].phi[coordinate], b.h[hit].phi[coordinate])) {
        return false;
      }
    }
  }
  return true;
}

// A cell and its measurements; an invalid cluster reference rejects it.
struct NeighbourCell {
  Triplet triplet;
  std::array<GlobalMeasurement, 3> measurements{};
  bool measurementsValid{false};
};
GPUhdi() NeighbourCell resolveCell(const FrameView& frame, const Triplet& triplet)
{
  NeighbourCell cell{triplet, {}, true};
  for (int hit = 0; hit < 3; ++hit) {
    const auto reference = triplet.getClusterReference(hit);
    const auto* cluster = frame.cluster(reference.surfacePosition, reference.clusterIndex);
    if (!cluster) {
      cell.measurementsValid = false;
      break;
    }
    cell.measurements[hit] = *cluster;
  }
  return cell;
}

GPUhdi() bool evaluateNeighbour(const NeighbourCell& first, const NeighbourCell& second,
                                const std::array<float, 2>& angularVariance, float maxChi2)
{
  for (int hit = 0; hit < 2; ++hit) {
    const auto a = first.triplet.getClusterReference(hit + 1);
    const auto b = second.triplet.getClusterReference(hit);
    if (a.surfacePosition != b.surfacePosition || a.clusterIndex != b.clusterIndex) {
      return false;
    }
  }
  AdjacentTripletFitResult fit{};
  return first.measurementsValid && second.measurementsValid &&
         fitAdjacentTripletFactors(first.triplet.tripletFactor(), second.triplet.tripletFactor(),
                                   {first.measurements[0], first.measurements[1], first.measurements[2], second.measurements[2]},
                                   angularVariance, fit) && !(fit.chi2 > maxChi2);
}

// A source and a target path, with the target's lookup table (first cell of
// every first-edge tracklet, plus the total).
struct NeighbourSearch {
  const Triplet* sources{};
  int sourcePath{-1};
  const Triplet* targets{};
  const int* targetLookup{};
  size_t lookupSize{};
  std::array<float, 2> angularVariance{};
  float maxChi2{};
};

// Every neighbour of one source cell, stopping at the first target whose
// time is incompatible.
template <typename Emit>
GPUhdi() void forEachNeighbour(const FrameView& frame, const NeighbourSearch& pair, int source, Emit emit)
{
  const auto first = resolveCell(frame, pair.sources[source]);
  const int tracklet = first.triplet.getSecondTrackletIndex();
  if (tracklet < 0 || size_t(tracklet) + 1 >= pair.lookupSize) {
    return;
  }
  for (int target = pair.targetLookup[tracklet]; target < pair.targetLookup[tracklet + 1]; ++target) {
    const auto& second = pair.targets[target];
    if (second.getFirstTrackletIndex() != tracklet || !first.triplet.getTimeStamp().isCompatible(second.getTimeStamp())) {
      break;
    }
    if (evaluateNeighbour(first, resolveCell(frame, second), pair.angularVariance, pair.maxChi2)) {
      emit(CellNeighbour{source, pair.sourcePath, target});
    }
  }
}

struct SeedInput {
  Triplet cell;
  std::array<GlobalMeasurement, 3> globals;
  std::array<SurfaceMeasurement, 3> measurements;
  std::array<SurfaceDescriptor, 3> surfaces;
};

GPUdi() bool initializeTrackSeed(const SeedInput& input, float bz, float maxChi2, TrackSeed& output)
{
  SurfaceTrackState state{};
  float chi2{0.f};
  const auto& outer = input.measurements[2];
  const auto kind = input.surfaces[2].kind;

  float sinPhi = 0.f, cosPhi = 0.f, tanLambda = 0.f, qOverPt = 1.f / o2::track::kMostProbablePt;
  float curvatureSquared = 1.f;

  state.referenceCoordinate = outer.frame.q;
  state.alpha = (kind == SurfaceKind::Cylinder) ? outer.frame.frameAngle : 0.f;
  state.parameters[0] = outer.frame.u;
  state.parameters[1] = outer.frame.v;

  float cosAlpha, sinAlpha, x[3], y[3];
  o2::its::math_utils::sinCosFloat(state.alpha, sinAlpha, cosAlpha);
  for (int i{0}; i < 3; ++i) {
    const auto& pos = input.globals[i].position;
    x[i] = pos.x * cosAlpha + pos.y * sinAlpha;
    y[i] = -pos.x * sinAlpha + pos.y * cosAlpha;
  }
  const float dx = x[2] - x[1];
  const float dy = y[2] - y[1];
  const float chordLength = std::sqrt(dx * dx + dy * dy);
  const float inverseLength = 1.f / chordLength;

  const float chordCos = dx * inverseLength;
  const float chordSin = dy * inverseLength;
  tanLambda = -0.5f *
              (o2::its::math_utils::computeTanDipAngle(x[0], y[0], x[1], y[1], input.globals[0].position.z, input.globals[1].position.z) +
               o2::its::math_utils::computeTanDipAngle(x[1], y[1], x[2], y[2], input.globals[1].position.z, input.globals[2].position.z));

  if (std::abs(bz) < 0.01f) {
    cosPhi = chordCos;
    sinPhi = chordSin;
  } else {
    const float curvature =
      o2::its::math_utils::computeCurvature(
        x[2], y[2], x[1], y[1], x[0], y[0]);

    const float halfSin = 0.5f * curvature * chordLength;
    const float halfCos =
      std::sqrt((1.f - halfSin) * (1.f + halfSin));

    cosPhi = chordCos * halfCos - chordSin * halfSin;
    sinPhi = chordSin * halfCos + chordCos * halfSin;
    qOverPt = curvature /
              (bz * o2::constants::math::B2C);
    curvatureSquared = curvature * curvature;
  }

  float phi = o2::gpu::GPUCommonMath::ASin(sinPhi);
  if (cosPhi < 0.f) {
    phi = o2::constants::math::PI - phi;
  } else if (phi < 0.f) {
    phi += o2::constants::math::TwoPI;
  }

  state.parameters[2] = (kind == SurfaceKind::Cylinder) ? sinPhi : phi;
  state.parameters[3] = tanLambda;
  state.parameters[4] = qOverPt;
  state.covariance[packedCovarianceIndex(0, 0)] = outer.covariance.uu;
  state.covariance[packedCovarianceIndex(1, 0)] = outer.covariance.uv;
  state.covariance[packedCovarianceIndex(1, 1)] = outer.covariance.vv;
  state.covariance[packedCovarianceIndex(2, 2)] = (kind == SurfaceKind::Cylinder) ? o2::track::kCSnp2max : o2::track::kCSnp2max / (cosPhi * cosPhi);
  state.covariance[packedCovarianceIndex(3, 3)] = o2::track::kCTgl2max;
  state.covariance[packedCovarianceIndex(4, 4)] = o2::track::kC1Pt2max * std::clamp(curvatureSquared, 0.0005f, 1.f);

  state.kind = kind;
  state.flags = 0;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;

  const std::array<const SurfaceMeasurement*, 2> attachmentMeasurements{&input.measurements[1], &input.measurements[0]};
  const std::array<const SurfaceDescriptor*, 2> attachmentSurfaces{&input.surfaces[1], &input.surfaces[0]};
  for (int step = 0; step < 2; ++step) {
    const auto& targetSurface = *attachmentSurfaces[step];
    if (!Propagator::attachMeasurement(
          state, targetSurface, *attachmentMeasurements[step], bz,
          material::MaterialTraversalDirection::OppositeMomentum,
          step == 1,
          maxChi2,
          chi2)) {
      return false;
    }
  }

  output = TrackSeed{input.cell, state, chi2};
  return true;
}

struct RoadSeedEmission {
  TrackSeed seed;
  int cellId{-1};
  int cellPathId{-1};
};
struct RoadJob {
  TrackSeed seed;
  size_t firstTarget{}, endTarget{};
  SeedInput initial{};
  bool initialize{false};
  // Index of the input seed this job extends; results report it as their job.
  size_t source{};
};
struct RoadTarget {
  Triplet cell;
  SurfaceDescriptor surface;
  SurfaceMeasurement measurement;
  int cellId, cellPathId;
  bool available;
};
struct RoadResult {
  size_t job;
  RoadSeedEmission emission;
};
#ifndef GPUCA_GPUCODE
inline bool roadSeedsEquivalent(const RoadSeedEmission& expected, const RoadSeedEmission& actual)
{
  const auto& a = expected.seed;
  const auto& b = actual.seed;
  const auto close = [](float x, float y) {
    return std::isfinite(x) && std::isfinite(y) && std::abs(x - y) <= 1.e-4f + 2.e-3f * std::max(std::abs(x), std::abs(y));
  };
  if (expected.cellId != actual.cellId || expected.cellPathId != actual.cellPathId ||
      a.getClusters() != b.getClusters() || a.getHitLayerMask().value() != b.getHitLayerMask().value() ||
      a.getLevel() != b.getLevel() || a.getFirstTrackletIndex() != b.getFirstTrackletIndex() || a.getSecondTrackletIndex() != b.getSecondTrackletIndex() ||
      a.getTimeStamp().getTimeStamp() != b.getTimeStamp().getTimeStamp() ||
      a.getTimeStamp().getTimeStampError() != b.getTimeStamp().getTimeStampError() || !(std::isfinite(a.getChi2()) && std::isfinite(b.getChi2()) &&
        std::abs(a.getChi2() - b.getChi2()) <= 1.e-2f + 5.e-3f * std::max(std::abs(a.getChi2()), std::abs(b.getChi2())))) {
    return false;
  }
  const auto& x = a.state();
  const auto& y = b.state();
  if (x.kind != y.kind || x.pid.getID() != y.pid.getID() || x.absCharge != y.absCharge || x.flags != y.flags ||
      !close(x.referenceCoordinate, y.referenceCoordinate) || !close(x.alpha, y.alpha)) {
    return false;
  }
  for (int i = 0; i < 5; ++i) {
    if (!close(x.parameters[i], y.parameters[i])) {
      return false;
    }
  }
  for (int i = 0; i < 15; ++i) {
    if (!close(x.covariance[i], y.covariance[i])) {
      return false;
    }
  }
  return true;
}
#endif
// The seeds a road start level keeps for the refit. Missing layers may be
// allowed, but do not count toward the minimum track length; q/pT is cut in
// parameters[4]'s units, and non-finite values fail the bound.
struct RoadSeedSelector {
  LayerMask nonSeedingLayers;
  LayerMask holeLayers;
  int maxHoles{};
  int minTrackLength{};
  float maxAbsQOverPt{};
  float maxChi2{};
  GPUhdi() bool operator()(const TrackSeed& seed) const
  {
    const auto hitLayerMask = seed.getHitLayerMask();
    const auto effectiveHoleMask = hitLayerMask.holeMask() & ~nonSeedingLayers;
    return effectiveHoleMask.isAllowedHoleMask(maxHoles, holeLayers) && hitLayerMask.count() >= minTrackLength &&
           o2::gpu::CAMath::Abs(seed.getQOverPt()) <= maxAbsQOverPt && seed.getChi2() <= maxChi2;
  }
};

// A cell path's cells and neighbour graph as the road stage reads them. Only
// paths with neighbours have a lookup table: the first neighbour of every
// cell, plus the total.
struct RoadGraphPath {
  const Triplet* cells{};
  size_t nCells{};
  const int* lookup{};
  size_t nLookup{};
  const CellNeighbour* neighbours{};
  size_t nNeighbours{};
};
// One road extension: the graph, the traversal layers, the surfaces and the
// level of the seeds it extends.
struct RoadStep {
  const RoadGraphPath* paths{};
  size_t nPaths{};
  int nLayers{};
  int currentLevel{};
  SurfaceCatalogView catalog{};
};
enum RoadError : unsigned {
  RoadMissingMeasurement = 1,
  RoadInvalidLayer = 2,
  RoadInvalidGraph = 4,
};

// The seed input of a start cell, or with no input only whether it has one.
GPUhdi() unsigned makeSeedInput(const FrameView& frame, const SurfaceCatalogView& catalog, const NeighbourCell& cell, SeedInput* input)
{
  for (int hit = 0; hit < 3; ++hit) {
    const auto reference = cell.triplet.getClusterReference(hit);
    const LayerId surface{static_cast<uint16_t>(reference.surfacePosition)};
    const auto* measurement = cell.measurementsValid ? frame.layers[reference.surfacePosition].measurement(cell.measurements[hit].clusterId) : nullptr;
    if (!measurement || !catalog.hasSurface(surface)) {
      return RoadMissingMeasurement;
    }
    if (input) {
      input->globals[hit] = cell.measurements[hit];
      input->measurements[hit] = *measurement;
      input->surfaces[hit] = catalog.getSurface(surface);
    }
  }
  if (input) {
    input->cell = cell.triplet;
  }
  return 0;
}

// The neighbours [begin, end) a seed extends to: start cell cellId of path
// (startCell) or a previous extension ending on it. Seeds of another level
// and start cells with a used cluster extend to none.
GPUhdi() unsigned roadRange(const FrameView& frame, const RoadStep& step, int path, int cellId, int level, bool startCell, int& begin, int& end)
{
  begin = end = 0;
  if (level != step.currentLevel) {
    return 0;
  }
  if (path < 0 || size_t(path) >= step.nPaths) {
    return startCell ? RoadInvalidGraph : 0;
  }
  const auto& graph = step.paths[path];
  NeighbourCell cell;
  if (startCell) {
    cell = resolveCell(frame, graph.cells[cellId]);
    for (int hit = 0; hit < 3; ++hit) {
      const auto reference = cell.triplet.getClusterReference(hit);
      if (reference.surfacePosition >= step.nLayers || reference.clusterIndex == o2::its::constants::UnusedIndex) {
        continue;
      }
      if (!cell.measurementsValid) {
        return RoadInvalidGraph;
      }
      if (frame.layers[reference.surfacePosition].isUsed(cell.measurements[hit].clusterId)) {
        return 0;
      }
    }
  }
  if (!graph.nLookup) {
    return 0;
  }
  if (cellId < 0 || size_t(cellId) + 1 >= graph.nLookup || graph.lookup[cellId] < 0 || size_t(graph.lookup[cellId + 1]) > graph.nNeighbours) {
    return RoadInvalidGraph;
  }
  begin = graph.lookup[cellId];
  end = o2::gpu::CAMath::Max(begin, graph.lookup[cellId + 1]);
  return startCell && begin == end ? makeSeedInput(frame, step.catalog, cell, nullptr) : 0; // as for seeds that are extended
}

GPUhdi() unsigned makeRoadTarget(const FrameView& frame, const RoadStep& step, const CellNeighbour& neighbour, RoadTarget& target)
{
  target = {};
  if (neighbour.cellPath < 0 || size_t(neighbour.cellPath) >= step.nPaths || neighbour.cell < 0 ||
      size_t(neighbour.cell) >= step.paths[neighbour.cellPath].nCells) {
    return RoadInvalidGraph;
  }
  const auto cell = resolveCell(frame, step.paths[neighbour.cellPath].cells[neighbour.cell]);
  const int layer = cell.triplet.getInnerLayer();
  const LayerId surface{static_cast<uint16_t>(layer)};
  if (layer < 0 || layer >= step.nLayers || !step.catalog.hasSurface(surface)) {
    return RoadInvalidLayer;
  }
  if (!cell.measurementsValid) {
    return RoadInvalidGraph;
  }
  const auto clusterId = cell.measurements[0].clusterId;
  const auto* measurement = frame.layers[layer].measurement(clusterId);
  target = {cell.triplet, step.catalog.getSurface(surface), measurement ? *measurement : SurfaceMeasurement{}, neighbour.cell,
            neighbour.cellPath, measurement && !frame.layers[layer].isUsed(clusterId)};
  return 0;
}

inline void throwRoadError(unsigned error)
{
  if (error & RoadMissingMeasurement) {
    throw std::invalid_argument{"CA seed: missing surface measurement"};
  }
  if (error & RoadInvalidLayer) {
    throw std::invalid_argument{"CA road traversal: invalid neighbour layer"};
  }
  if (error) {
    throw std::runtime_error{"CA road traversal: invalid road graph"};
  }
}

static_assert(std::is_trivially_copyable_v<RoadJob>);
static_assert(std::is_trivially_copyable_v<RoadTarget>);
static_assert(std::is_trivially_copyable_v<RoadResult>);

GPUdi() bool extendRoad(const TrackSeed& current, const RoadTarget& target, float bz, float maxChi2, RoadSeedEmission& output)
{
  const auto& neighbour = target.cell;
  if (!target.available || neighbour.getSecondTrackletIndex() != current.getFirstTrackletIndex() ||
      !current.getTimeStamp().isCompatible(neighbour.getTimeStamp()) || current.getLevel() - 1 != neighbour.getLevel()) {
    return false;
  }
  TrackSeed seed{current};
  seed.getTimeStamp() += neighbour.getTimeStamp();
  float chi2 = seed.getChi2();
  if (!Propagator::attachMeasurement(seed.state(), target.surface, target.measurement, bz,
                                      material::MaterialTraversalDirection::OppositeMomentum, true, maxChi2, chi2)) {
    return false;
  }
  seed.setChi2(chi2);
  const int layer = neighbour.getInnerLayer();
  seed.setCluster(layer, neighbour.getFirstClusterIndex());
  auto mask = seed.getHitLayerMask();
  mask.set(layer);
  seed.setHitLayerMask(mask);
  seed.setLevel(neighbour.getLevel());
  seed.setFirstTrackletIndex(neighbour.getFirstTrackletIndex());
  seed.setSecondTrackletIndex(neighbour.getSecondTrackletIndex());
  output = {seed, target.cellId, target.cellPathId};
  return true;
}

// Every extension of a seed (initialised first if it is a start cell's) by
// the targets [first, end).
template <typename Emit>
GPUdi() void forEachRoad(TrackSeed seed, const SeedInput* initial, const RoadTarget* targets, size_t first, size_t end, float bz,
                          float maxChi2, Emit emit)
{
  if (first == end || (initial && !initializeTrackSeed(*initial, bz, maxChi2, seed))) {
    return;
  }
  for (size_t target = first; target < end; ++target) {
    RoadSeedEmission result;
    if (extendRoad(seed, targets[target], bz, maxChi2, result)) {
      emit(result);
    }
  }
}

struct RefitJob {
  TrackSeed seed;
  std::array<detail::RefitMeasurementSlot, MaxLayoutSurfaces> slots{};
  std::array<detail::CircleFitPoint, MaxLayoutSurfaces> points{};
  size_t nSlots{}, nPoints{};
  float minPt{};
};
struct RefitParameters {
  float bz{}, maxChi2ClusterAttachment{}, maxChi2NDF{};
  bool shiftReferenceToMeasurement{}, repeatRefitOut{};
};
struct RefitResult {
  SurfaceTrackState inner{}, outer{};
  float chi2{};
  bool accepted{false};
};
static_assert(std::is_trivially_copyable_v<RefitJob>);
static_assert(std::is_trivially_copyable_v<RefitResult>);
// A seed's refit job; minPt is the minimum pT for each number of missing
// layers.
GPUhdi() void prepareRefitJob(const FrameView& frame, const TrackSeed& seed, const float* minPt, size_t nMinPt, RefitJob& job)
{
  job.seed = seed;
  job.nSlots = prepareRefitSlots(frame, seed, job.slots.data(), job.points.data(), job.nPoints);
  const int missing = static_cast<int>(job.nSlots) - seed.getHitLayerMask().count();
  job.minPt = missing >= 0 && size_t(missing) < nMinPt ? minPt[missing] : 0.f;
}
GPUdi() RefitResult evaluateRefit(const RefitJob& job, SurfaceCatalogView catalog, const RefitParameters& p)
{
  RefitResult result{};
  result.accepted = fitPreparedTrackStateLegs(job.seed.state(), {job.slots.data(), job.nSlots}, {job.points.data(), job.nPoints},
                                            catalog, p.bz, p.shiftReferenceToMeasurement, p.maxChi2ClusterAttachment,
                                            p.maxChi2NDF, p.repeatRefitOut, job.minPt, result.inner, result.outer, result.chi2);
  return result;
}
#ifndef GPUCA_GPUCODE
inline bool refitsEquivalent(const RefitResult& a, const RefitResult& b)
{
  if (a.accepted != b.accepted) {
    return false;
  }
  if (!a.accepted) {
    return true;
  }
  const auto close = [](float x, float y, float absolute, float relative) {
    return std::isfinite(x) && std::isfinite(y) && std::abs(x - y) <= absolute + relative * std::max(std::abs(x), std::abs(y));
  };
  const auto stateEquivalent = [&](const auto& x, const auto& y) {
    if (x.kind != y.kind || x.flags != y.flags || x.absCharge != y.absCharge || x.pid.getID() != y.pid.getID() ||
        !close(x.referenceCoordinate, y.referenceCoordinate, 1.e-5f, 1.e-4f) || !close(x.alpha, y.alpha, 1.e-6f, 1.e-4f)) {
      return false;
    }
    for (int i = 0; i < 5; ++i) {
      if (!close(x.parameters[i], y.parameters[i], 1.e-5f, 2.e-3f)) {
        return false;
      }
    }
    for (int i = 0; i < 15; ++i) {
      if (!close(x.covariance[i], y.covariance[i], 1.e-8f, 2.e-2f)) {
        return false;
      }
    }
    return true;
  };
  return close(a.chi2, b.chi2, 1.e-2f, 5.e-3f) && stateEquivalent(a.inner, b.inner) && stateEquivalent(a.outer, b.outer);
}
#endif

} // namespace o2::itsmft::tracking
#endif

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

#include <boost/test/unit_test.hpp>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>
#include "ITSMFTTrackingGPU/TimeFrameGPU.h"

using namespace o2::itsmft::tracking;

namespace
{
// One path's search from its two edges' tracklets.
struct CellSearchInput {
  int firstEdge, secondEdge;
  std::array<int, 3> layers;
  CellParameters parameters;
};

// A timeframe with hand-made clusters and tracklets.
struct CellRig {
  int device{gpu::currentDevice()};
  gpu::TimeFrameGPU frame{device};

  // Clusters of every layer, then every edge's tracklets with their lookup
  // table (first tracklet of every source cluster, plus the total).
  void loadClusters(std::initializer_list<std::span<const GlobalMeasurement>> layers)
  {
    frame.initialise(layers.size());
    int index = 0;
    for (const auto& layer : layers) {
      gpu::LayerHostData host;
      host.clusters = layer;
      frame.loadLayer(index++, host);
    }
    frame.finishLoading();
  }
  void loadTracklets(std::initializer_list<std::pair<std::span<const Tracklet>, std::span<const int>>> edges, size_t nPaths)
  {
    frame.initialiseTracklets(edges.size());
    int index = 0;
    for (const auto& [tracklets, lookup] : edges) {
      frame.upload(frame.tracklets(index).items, tracklets);
      frame.upload(frame.tracklets(index++).lookup, lookup);
    }
    frame.initialiseCells(nPaths);
  }
  // Runs one path with the given buffer capacity; returns the number found.
  int search(const CellSearchInput& input, int path, int capacity)
  {
    frame.cells(path).reserve(capacity);
    CellSearch search;
    search.layers = input.layers;
    search.parameters = input.parameters;
    return gpu::computeCells(frame, search, input.firstEdge, input.secondEdge, path);
  }
  // As the traits do: grows the buffer until everything fits, then downloads.
  std::vector<Triplet> cells(const CellSearchInput& input, int path)
  {
    int capacity = 0;
    for (int found = search(input, path, capacity); found > capacity; found = search(input, path, capacity)) {
      capacity = found;
    }
    std::vector<Triplet> result;
    frame.download(frame.cells(path).items, result);
    frame.download(frame.cells(path).lookup, lookup);
    return result;
  }
  std::vector<int> lookup;
};

Tracklet segment(int first, int second, const GlobalMeasurement& a, const GlobalMeasurement& b, o2::its::TimeEstBC time)
{
  return Tracklet{first, second, (b.z - a.z) / std::hypot(b.x - a.x, b.y - a.y), std::atan2(a.y - b.y, a.x - b.x), time};
}
} // namespace

BOOST_AUTO_TEST_CASE(DeviceCellOrderingTimingDegeneracyAndCapacity)
{
  CellRig rig;
  std::array<GlobalMeasurement, 2> inner{}, middle{}, outer{};
  inner[0].position = {3.f, 0.1f, 0.9f};
  middle[0].position = {4.f, 0.15f, 1.05f};
  outer[0].position = {5.f, 0.201f, 1.2f};
  inner[1] = inner[0];
  middle[1] = middle[0];
  outer[1] = middle[0]; // coincident hits must not produce a fit factor
  std::vector<Tracklet> first;
  for (int i = 0; i < 257; ++i) {
    first.push_back(segment(i % 2, i % 2, inner[i % 2], middle[i % 2], {100, 20}));
  }
  const std::vector<Tracklet> second{
    segment(0, 0, middle[0], outer[0], {110, 20}),
    segment(0, 0, middle[0], outer[0], {120, 20}), // touching endpoints do not overlap
    Tracklet{1, 1, first[1].tanLambda, first[1].phi, {100, 20}}};
  const std::vector<int> firstLookup{0, 129, 257}, secondLookup{0, 2, 3};
  const CellSearchInput input{0, 1, {0, 1, 2}, {0.1f, 0.1f, false}};
  const auto load = [&] {
    rig.loadClusters({inner, middle, outer});
    rig.loadTracklets({{first, firstLookup}, {second, secondLookup}}, 1);
  };
  load();
  auto cells = rig.cells(input, 0);
  BOOST_REQUIRE_EQUAL(cells.size(), 129u);
  TripletFitFactor expected;
  BOOST_REQUIRE(makeTripletFitFactor({inner[0], middle[0], outer[0]}, expected));
  for (size_t i = 0; i < cells.size(); ++i) {
    // Sorted by (first, second) tracklet, as the host emits them.
    BOOST_CHECK_EQUAL(cells[i].getFirstTrackletIndex(), 2 * static_cast<int>(i));
    BOOST_CHECK_EQUAL(cells[i].getSecondTrackletIndex(), 0);
    BOOST_CHECK_EQUAL(cells[i].getLevel(), 1);
    BOOST_CHECK(cellFactorsEquivalent(cells[i].tripletFactor(), expected));
  }
  // First cell of every first-edge tracklet, plus the total.
  BOOST_REQUIRE_EQUAL(rig.lookup.size(), first.size() + 1);
  for (size_t i = 0; i < rig.lookup.size(); ++i) {
    BOOST_CHECK_EQUAL(rig.lookup[i], static_cast<int>((i + 1) / 2));
  }

  // Too small a buffer: everything is counted, nothing else is done.
  BOOST_CHECK_EQUAL(rig.search(input, 0, 1), 129);
  BOOST_CHECK_EQUAL(rig.frame.cells(0).items.size, 0u);
  BOOST_CHECK_EQUAL(rig.cells(input, 0).size(), 129u);

  // A new timeframe with the same layout but new contents.
  const auto saved = outer[0];
  outer[0] = middle[0];
  load();
  BOOST_CHECK(rig.cells(input, 0).empty());
  outer[0] = saved;
  load();
  BOOST_CHECK_EQUAL(rig.cells(input, 0).size(), 129u);
  auto tight = input;
  tight.parameters.angularTolerance = 0.f;
  BOOST_CHECK(rig.cells(tight, 0).empty());
  rig.loadTracklets({{{}, firstLookup}, {second, secondLookup}}, 1);
  BOOST_CHECK(rig.cells(input, 0).empty());
}

BOOST_AUTO_TEST_CASE(DeviceCellSmallCurvatureFactorDoesNotAmplifyFloatRounding)
{
  // Real-event regression: float intermediate asin rounding changed rho.theta
  // by about 1.9e-4 between CPU and CUDA for this triplet.
  CellRig rig;
  std::array<GlobalMeasurement, 1> inner{}, middle{}, outer{};
  inner[0].position = {-0x1.04af2ep+0f, -0x1.80e2a2p+1f, -0x1.4720bcp+2f};
  middle[0].position = {-0x1.5dbafcp+0f, -0x1.02dfb4p+2f, -0x1.b8332ep+2f};
  outer[0].position = {-0x1.a3422p+2f, -0x1.2b6c9ep+4f, -0x1.00047ap+5f};
  const std::array<Tracklet, 1> first{Tracklet{0, 0, 0.f, 0.f, {100, 20}}};
  const std::array<Tracklet, 1> second{Tracklet{0, 0, 0.f, 0.f, {110, 20}}};
  const std::array<int, 2> lookup{0, 1};
  rig.loadClusters({inner, middle, outer});
  rig.loadTracklets({{first, lookup}, {second, lookup}}, 1);
  const auto cells = rig.cells({0, 1, {0, 1, 2}, {1.f, 1.f, false}}, 0);
  BOOST_REQUIRE_EQUAL(cells.size(), 1u);
  TripletFitFactor reference;
  BOOST_REQUIRE(makeTripletFitFactor({inner[0], middle[0], outer[0]}, reference));
  BOOST_CHECK(cellFactorsEquivalent(cells[0].tripletFactor(), reference));
  BOOST_CHECK_SMALL(cells[0].tripletFactor().rho.theta - reference.rho.theta, 1.e-7f);
}

// The road stage builds jobs and targets from the resident cells and road
// graph exactly as the host does, and stale or inconsistent inputs fail.
BOOST_AUTO_TEST_CASE(DeviceBuiltRoadsMatchHostBuiltRoads)
{
  CellRig rig;
  const auto device = rig.device;
  auto& frame = rig.frame;
  // One cluster per cylinder layer at r = 3, 4, 5, 6, all with cluster id 0.
  constexpr int nLayers = 4;
  std::array<std::array<GlobalMeasurement, 1>, nLayers> globals{};
  std::array<std::array<SurfaceMeasurement, 1>, nLayers> surfaceMeasurements{};
  std::array<SurfaceDescriptor, nLayers> surfaces{};
  const std::array<float, nLayers> y{0.1f, 0.15f, 0.201f, 0.253f};
  for (int layer = 0; layer < nLayers; ++layer) {
    const float r = 3.f + layer;
    globals[layer][0].position = {r, y[layer], 0.15f * r + 0.45f};
    globals[layer][0].clusterId = 0;
    surfaceMeasurements[layer][0] = {{r, y[layer], 0.15f * r + 0.45f, 0.f}, {1.e-4f, 0.f, 1.e-4f}};
    surfaces[layer].kind = SurfaceKind::Cylinder;
    surfaces[layer].referenceCoordinate = r;
  }
  const auto layerSegment = [&](int layer) { return segment(0, 0, globals[layer][0], globals[layer + 1][0], {100, 20}); };
  constexpr size_t nSeeds = 64;
  const std::vector<Tracklet> t01(1, layerSegment(0)), t12(nSeeds, layerSegment(1)), t23(1, layerSegment(2));
  const std::vector<int> lookup01{0, 1}, lookup12{0, nSeeds}, lookup23{0, 1};

  // The timeframe as TrackerTraitsGPU uploads it.
  std::array<std::array<uint8_t, 1>, nLayers> used{};
  std::array<std::span<const uint8_t>, nLayers> usedSpans;
  const auto loadFrame = [&](int withoutMeasurements) {
    frame.initialise(nLayers);
    for (int layer = 0; layer < nLayers; ++layer) {
      usedSpans[layer] = used[layer];
      gpu::LayerHostData host;
      host.clusters = globals[layer];
      if (layer != withoutMeasurements) {
        host.surfaceMeasurements = surfaceMeasurements[layer];
      }
      frame.loadLayer(layer, host);
    }
    frame.loadUsedClusters(usedSpans);
    frame.finishLoading();
  };
  loadFrame(-1);
  // Path 0: cells on layers 0-2 (cell k uses t12[k]); path 1: layers 1-3.
  std::vector<Triplet> innerCells, outerCells;
  const auto buildCells = [&] {
    rig.loadTracklets({{t01, lookup01}, {t12, lookup12}, {t23, lookup23}}, 2);
    innerCells = rig.cells({0, 1, {0, 1, 2}, {0.1f, 0.1f, false}}, 0);
    outerCells = rig.cells({1, 2, {1, 2, 3}, {0.1f, 0.1f, false}}, 1);
  };
  buildCells();
  BOOST_REQUIRE_EQUAL(innerCells.size(), nSeeds);
  BOOST_REQUIRE_EQUAL(outerCells.size(), nSeeds);
  const auto withLevel = [](Triplet cell, int level) {
    cell.setLevel(level);
    return cell;
  };

  // The neighbour stage on the device: outer cell k's only neighbour is
  // inner cell k, which raises it to level 2.
  const auto findNeighbours = [&] {
    frame.initialiseNeighbours(2);
    frame.neighbours(1).reserve(nSeeds);
    const std::array<NeighbourSearch, 1> sources{NeighbourSearch{nullptr, 0, nullptr, nullptr, 0, {1.e-4f, 1.e-4f}, 1.e6f}};
    BOOST_REQUIRE_EQUAL(gpu::computeNeighbours(rig.frame, 1, sources), static_cast<int>(nSeeds));
  };
  findNeighbours();
  std::vector<CellNeighbour> graph;
  frame.download(frame.neighbours(1).items, graph);
  BOOST_REQUIRE_EQUAL(graph.size(), nSeeds);
  for (size_t k = 0; k < nSeeds; ++k) {
    BOOST_REQUIRE_EQUAL(innerCells[k].getSecondTrackletIndex(), static_cast<int>(k));
    BOOST_CHECK_EQUAL(graph[k].cell, static_cast<int>(k));
    BOOST_CHECK_EQUAL(graph[k].cellPath, 0);
    BOOST_CHECK_EQUAL(graph[k].nextCell, static_cast<int>(k));
  }
  frame.download(frame.cells(1).items, outerCells);
  for (const auto& cell : outerCells) {
    BOOST_CHECK_EQUAL(cell.getLevel(), 2);
  }

  // The host-built reference, as TrackerTraits::buildRoadJobs makes it.
  const auto hostJobs = [&](const std::array<std::array<uint8_t, 1>, nLayers>& used) {
    std::pair<bounded_vector<RoadJob>, bounded_vector<RoadTarget>> built;
    for (size_t k = 0; k < nSeeds; ++k) {
      if (used[1][0] || used[2][0] || used[3][0]) {
        continue;
      }
      RoadJob job{};
      job.source = k;
      job.firstTarget = built.second.size();
      job.endTarget = job.firstTarget + 1;
      job.initialize = true;
      job.initial.cell = withLevel(outerCells[k], 2);
      for (int hit = 0; hit < 3; ++hit) {
        job.initial.globals[hit] = globals[hit + 1][0];
        job.initial.measurements[hit] = surfaceMeasurements[hit + 1][0];
        job.initial.surfaces[hit] = surfaces[hit + 1];
      }
      built.first.push_back(job);
      built.second.push_back({withLevel(innerCells[k], 1), surfaces[0], surfaceMeasurements[0][0], static_cast<int>(k), 0, !used[0][0]});
    }
    return built;
  };
  constexpr float bz = 5.f, maxChi2 = 1.e6f;
  gpu::RoadExtension onDevice{1, false, 2, nLayers, bz, maxChi2, surfaces};
  const auto roads = [&](const gpu::RoadExtension& extension) {
    bounded_vector<RoadResult> output;
    gpu::extendRoads(frame, extension);
    frame.download(frame.roads[frame.lastRoads], output);
    return output;
  };
  // Not bound yet.
  BOOST_CHECK_THROW(roads(onDevice), std::runtime_error);
  gpu::bindRoadGraph(frame);
  const auto compare = [&] {
    frame.loadUsedClusters(usedSpans);
    const auto fromDevice = roads(onDevice);
    const auto [jobs, targets] = hostJobs(used);
    bounded_vector<RoadResult> fromHost;
    gpu::computeRoads(frame, {jobs, targets, bz, maxChi2, false, nSeeds}, fromHost);
    BOOST_REQUIRE_EQUAL(fromDevice.size(), fromHost.size());
    for (size_t i = 0; i < fromHost.size(); ++i) {
      BOOST_CHECK_EQUAL(fromDevice[i].job, fromHost[i].job);
      BOOST_CHECK(roadSeedsEquivalent(fromHost[i].emission, fromDevice[i].emission));
    }
    return fromDevice.size();
  };
  BOOST_CHECK_EQUAL(compare(), nSeeds);
  // Kept on the device, the same extension feeds the seed selection, which
  // appends the passing seeds in order.
  const auto expected = roads(onDevice);
  BOOST_REQUIRE_EQUAL(gpu::extendRoads(frame, onDevice), nSeeds);
  bounded_vector<RoadResult> kept;
  frame.download(frame.roads[frame.lastRoads], kept);
  BOOST_REQUIRE_EQUAL(kept.size(), expected.size());
  for (size_t i = 0; i < kept.size(); ++i) {
    BOOST_CHECK(roadSeedsEquivalent(expected[i].emission, kept[i].emission));
  }
  RoadSeedSelector selector{LayerMask{}, LayerMask{}, 0, 4, 1.e3f, 1.e30f};
  frame.trackSeeds.size = 0;
  BOOST_CHECK_EQUAL(gpu::selectRoadSeeds(frame, selector), nSeeds);
  BOOST_CHECK_EQUAL(gpu::selectRoadSeeds(frame, selector), 2 * nSeeds); // appended
  std::vector<TrackSeed> seeds;
  frame.download(frame.trackSeeds, seeds);
  for (size_t i = 0; i < seeds.size(); ++i) {
    BOOST_CHECK_EQUAL(seeds[i].getChi2(), kept[i % nSeeds].emission.seed.getChi2());
    BOOST_CHECK(selector(seeds[i]));
  }
  selector.maxChi2 = -1.f;
  BOOST_CHECK_EQUAL(gpu::selectRoadSeeds(frame, selector), 2 * nSeeds);

  // The seed refit prepares the same slots as the host and sorts the
  // accepted tracks longest first, then by chi2.
  const RefitParameters refitParameters{bz, 1.e6f, 1.e6f, false, false};
  const size_t nTracks = gpu::refitTrackSeeds(frame, surfaces, refitParameters, {});
  std::vector<RefitResult> results;
  frame.download(frame.refitResults, results);
  bounded_vector<RefitJob> refitJobs(seeds.size());
  for (size_t i = 0; i < seeds.size(); ++i) {
    auto& job = refitJobs[i];
    job.seed = seeds[i];
    for (int layer = 0; layer < nLayers; ++layer) {
      const auto& global = globals[layer][0];
      job.slots[layer] = {surfaceMeasurements[layer][0], LayerId{static_cast<uint16_t>(layer)}, true};
      job.points[layer] = {global.x, global.y, global.covariance.xx, global.covariance.xy, global.covariance.yy};
    }
    job.nSlots = job.nPoints = nLayers;
  }
  bounded_vector<RefitResult> reference;
  gpu::refitTracks(frame, {refitJobs, surfaces, refitParameters}, reference);
  size_t accepted = 0;
  for (size_t i = 0; i < seeds.size(); ++i) {
    BOOST_CHECK(refitsEquivalent(reference[i], results[i]));
    accepted += results[i].accepted;
  }
  BOOST_CHECK_EQUAL(nTracks, accepted);
  std::vector<TrackingCandidate> tracks;
  frame.download(frame.tracks, tracks);
  BOOST_CHECK_EQUAL(tracks.size(), nTracks);
  for (size_t i = 1; i < tracks.size(); ++i) {
    const int previous = tracks[i - 1].seed.getActiveLayerCount(), current = tracks[i].seed.getActiveLayerCount();
    BOOST_CHECK(previous > current || (previous == current && tracks[i - 1].track.chi2 <= tracks[i].track.chi2));
  }
  // The continuing expansion reads the device results; inner cells have no
  // neighbours, so it ends there.
  auto next = onDevice;
  next.reusePreviousSeeds = true;
  next.currentLevel = 1;
  BOOST_CHECK(roads(next).empty());

  used[0][0] = 1; // targets become unavailable
  BOOST_CHECK_EQUAL(compare(), 0u);
  used[0][0] = 0;
  used[2][0] = 1; // seeds are skipped
  BOOST_CHECK_EQUAL(compare(), 0u);
  used[2][0] = 0;
  BOOST_CHECK_EQUAL(compare(), nSeeds);

  onDevice.nLayers = 0; // every neighbour is outside the traversal
  BOOST_CHECK_THROW(roads(onDevice), std::invalid_argument);
  onDevice.nLayers = nLayers;
  loadFrame(3); // the seeds' outer measurement is missing
  BOOST_CHECK_THROW(roads(onDevice), std::invalid_argument);
  loadFrame(-1);
  BOOST_CHECK_EQUAL(roads(onDevice).size(), nSeeds);
  buildCells(); // a new cell stage makes the graph stale
  BOOST_CHECK_THROW(roads(onDevice), std::runtime_error);
  findNeighbours();
  gpu::bindRoadGraph(frame);
  BOOST_CHECK_EQUAL(roads(onDevice).size(), nSeeds);
}

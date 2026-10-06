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


#define BOOST_TEST_MODULE CommonTrackingNeighbourKernels
#include <boost/test/unit_test.hpp>
#include <vector>
#include "ITSMFTTrackingGPU/TimeFrameGPU.h"

using namespace o2::itsmft::tracking;

namespace
{
// Source cells on path 0 and target cells on path 1, with the targets'
// lookup table (first target cell of every first-edge tracklet); the cells'
// measurements are the clusters of the given layers.
struct NeighbourRig {
  int device{gpu::currentDevice()};
  gpu::TimeFrameGPU frame{device};
  void loadClusters(const std::vector<std::vector<GlobalMeasurement>>& layers)
  {
    frame.initialise(layers.size());
    for (size_t layer = 0; layer < layers.size(); ++layer) {
      gpu::LayerHostData host;
      host.clusters = layers[layer];
      frame.loadLayer(static_cast<int>(layer), host);
    }
    frame.finishLoading();
  }
  // Loads the cells, runs the search as the traits do (growing the buffer
  // until everything fits) and returns the target's neighbours.
  std::vector<CellNeighbour> run(const std::vector<Triplet>& sources, const std::vector<Triplet>& targets,
                                      std::span<const int> lookup, std::array<float, 2> variance, float maxChi2)
  {
    frame.initialiseCells(2);
    frame.upload(frame.cells(0).items, std::span<const Triplet>{sources});
    frame.upload(frame.cells(1).items, std::span<const Triplet>{targets});
    frame.upload(frame.cells(1).lookup, lookup);
    frame.initialiseNeighbours(2);
    int capacity = 0;
    for (;;) {
      frame.neighbours(1).reserve(capacity);
      const std::array<NeighbourSearch, 1> sources{NeighbourSearch{nullptr, 0, nullptr, nullptr, 0, variance, maxChi2}};
      const int found = gpu::computeNeighbours(frame, 1, sources);
      if (found <= capacity) {
        break;
      }
      capacity = found;
    }
    std::vector<CellNeighbour> neighbours;
    frame.download(frame.neighbours(1).items, neighbours);
    frame.download(frame.neighbours(1).lookup, lookupTable);
    std::vector<Triplet> cells;
    frame.download(frame.cells(1).items, cells);
    levels.clear();
    for (const auto& cell : cells) {
      levels.push_back(cell.getLevel());
    }
    return neighbours;
  }
  std::vector<int> lookupTable, levels;
};
} // namespace

BOOST_AUTO_TEST_CASE(NeighbourOrderingCutsAndLevels)
{
  NeighbourRig rig;
  std::array<GlobalMeasurement, 4> hits{};
  for (int hit = 0; hit < 4; ++hit) {
    hits[hit].position = {3.f + hit, 0.01f * hit * hit, 0.15f * hit};
    hits[hit].covariance = {1.e-4f, 0.f, 0.f, 1.e-4f, 0.f, 1.e-4f};
  }
  // One cluster per layer: hit l on layer l.
  rig.loadClusters({{hits[0]}, {hits[1]}, {hits[2]}, {hits[3]}});
  Triplet first{0, 0, 0, 0, 0, 0, {100, 20}};
  Triplet second{1, 0, 0, 0, 0, 0, {110, 20}};
  BOOST_REQUIRE(makeTripletFitFactor({hits[0], hits[1], hits[2]}, first.tripletFactor()));
  BOOST_REQUIRE(makeTripletFitFactor({hits[1], hits[2], hits[3]}, second.tripletFactor()));
  std::vector<Triplet> sources(257, first), targets(3, second);
  // One rejected candidate between accepted ones must not stop enumeration:
  // its outer cluster does not exist.
  targets[1].getClusters()[2] = 1;
  std::array<int, 2> lookup{0, 3};
  const std::array<float, 2> variance{1.e-4f, 1.e-4f};
  auto run = [&](float maxChi2 = 1.e6f) { return rig.run(sources, targets, lookup, variance, maxChi2); };
  auto output = run();
  // Sorted by target cell, then source cell.
  BOOST_REQUIRE_EQUAL(output.size(), 514u);
  for (size_t i = 0; i < output.size(); ++i) {
    BOOST_CHECK_EQUAL(output[i].nextCell, i < 257 ? 0 : 2);
    BOOST_CHECK_EQUAL(output[i].cell, static_cast<int>(i % 257));
    BOOST_CHECK_EQUAL(output[i].cellPath, 0);
  }
  const std::vector<int> expectedLookup{0, 257, 257, 514};
  BOOST_CHECK_EQUAL_COLLECTIONS(rig.lookupTable.begin(), rig.lookupTable.end(), expectedLookup.begin(), expectedLookup.end());
  // Targets reached from a level-1 source reach level 2; the rejected one keeps its own.
  const std::vector<int> expectedLevels{2, 1, 2};
  BOOST_CHECK_EQUAL_COLLECTIONS(rig.levels.begin(), rig.levels.end(), expectedLevels.begin(), expectedLevels.end());
  sources[5].setLevel(3);
  run();
  const std::vector<int> raisedLevels{4, 1, 4};
  BOOST_CHECK_EQUAL_COLLECTIONS(rig.levels.begin(), rig.levels.end(), raisedLevels.begin(), raisedLevels.end());
  sources[5] = first;

  // Preserve the existing break on an incompatible timestamp, even when a
  // later entry would be compatible. Endpoint contact does not overlap.
  targets[1].getTimeStamp() = {120, 20};
  output = run();
  BOOST_REQUIRE_EQUAL(output.size(), 257u);
  for (const auto& result : output) {
    BOOST_CHECK_EQUAL(result.nextCell, 0);
  }
  targets[1] = second;
  targets[1].setFirstTrackletIndex(1);
  BOOST_CHECK_EQUAL(run().size(), 257u);
  targets[1] = second;
  targets[0].setHitLayerMask(LayerMask{0, 2, 3});
  BOOST_CHECK_EQUAL(run().size(), 514u);
  targets[0] = second;
  sources[0].setSecondTrackletIndex(-1);
  sources[1].setSecondTrackletIndex(1);
  BOOST_CHECK_EQUAL(run().size(), 765u);
  sources[0] = sources[1] = first;
  BOOST_CHECK_EQUAL(run().size(), 771u);
  // Sparse paths may cross holes as long as the two shared references match.
  rig.loadClusters({{hits[0]}, {}, {hits[1]}, {}, {hits[2]}, {}, {hits[3]}});
  for (auto& source : sources) {
    source.setHitLayerMask(LayerMask{0, 2, 4});
  }
  for (auto& target : targets) {
    target.setHitLayerMask(LayerMask{2, 4, 6});
  }
  BOOST_CHECK_EQUAL(run().size(), 771u);
  // A positive chi2 threshold on either side of the computed fit exercises
  // the attachment cut, independently of candidate identity and timing.
  sources.resize(1);
  targets.resize(1);
  lookup[1] = 1;
  targets[0].tripletFactor().psi.theta += 0.1f;
  AdjacentTripletFitResult reference{};
  BOOST_REQUIRE(fitAdjacentTripletFactors(sources[0].tripletFactor(), targets[0].tripletFactor(),
                                          hits, variance, reference));
  BOOST_REQUIRE_GT(reference.chi2, 0.f);
  BOOST_CHECK(run(reference.chi2 * 0.5f).empty());
  BOOST_CHECK_EQUAL(run(reference.chi2 * 2.f).size(), 1u);
  BOOST_CHECK(run(-1.f).empty());
  targets[0].tripletFactor() = {};
  sources[0].tripletFactor() = {};
  BOOST_CHECK(run().empty());
  sources.clear();
  BOOST_CHECK(run().empty());

  // Paths must exist and differ.
  const std::array<NeighbourSearch, 1> source{NeighbourSearch{nullptr, 0, nullptr, nullptr, 0, variance, 1.f}};
  BOOST_CHECK_THROW(gpu::computeNeighbours(rig.frame, 0, source), std::out_of_range);
  BOOST_CHECK_THROW(gpu::computeNeighbours(rig.frame, 5, source), std::out_of_range);
}

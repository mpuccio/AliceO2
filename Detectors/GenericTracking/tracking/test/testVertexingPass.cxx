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

/// \file testVertexingPass.cxx
/// \brief The own vertexer before the tracking: one vertex per collision from the pass over the innermost
///        layers, the same list at any thread count, and the vertex-list
///        iterations finding with it the tracks they find with the true
///        vertices.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK VertexingPass
#include <boost/test/unit_test.hpp>
#include <numbers>

#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/Tracker.h"
#include "GenericTracking/Tracking/VertexingPass.h"

#include "ToyBarrel8Fixture.h"
#include "UniformMaterial.h"

#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

namespace o2::gtrk
{

namespace
{

using namespace toybarrel8;
using Pass = VertexingPass<Det, UniformMaterial>;

/// Two collisions of eight tracks each, 500 ns apart, every track 0.5 to 2.6
/// GeV from the origin.
constexpr std::size_t kCollisions = 2U;
constexpr std::size_t kTracksPerCollision = 8U;
constexpr std::size_t kTracks = kCollisions * kTracksPerCollision;
constexpr std::array<double, kCollisions> kStartNS{37., 537.};
constexpr double kDurationNS = 600.;

std::vector<Hit> hits()
{
  std::vector<Hit> result;
  for (std::size_t track = 0; track < kTracks; ++track) {
    const auto k = static_cast<double>(track);
    const double pt = 0.5 + (0.3 * static_cast<double>(track % kTracksPerCollision));
    const auto one =
      sensorHelixHits(std::fmod(0.3 + (2.1 * k), 2 * std::numbers::pi_v<double>), 0.03 * std::cos(1.3 * k), pt, kStartNS[track / kTracksPerCollision]);
    BOOST_REQUIRE(one.size() == kLayers);
    result.insert(result.end(), one.begin(), one.end());
  }
  return result;
}

/// The pass: the beam line, the three innermost layers, no holes, 0.3 GeV.
IterationSettings passSettings()
{
  auto settings = iterationSettings();
  settings.windows.minPtGeV = 0.3;
  settings.minClusters = 3;
  settings.maxKinks = 0;
  settings.skipLinks = false;
  settings.layerMask = 0b111U;
  settings.follow = false;
  return settings;
}

ScatteringTable<Det> noScattering()
{
  const std::array<double, kLayers> material{};
  return scatteringTable<Det>(material);
}

/// The pass's slices: three blocks, as the tracking's in these tests.
constexpr std::int32_t kSliceBlocks = 3;

Pass::Result findVertices(const Fixture& fixture, int threads, std::int32_t sliceBlocks = kSliceBlocks)
{
  Pass pass{passSettings(), VertexerSettings{}, Fixture::kBeamConstraint, noScattering(), std::make_shared<BoundedMemoryResource>(), 2U, threads};
  auto result = pass.find(*fixture.frame, sliceBlocks, 0);
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  return result;
}

/// A track as compared: its iteration, its rows and its chi2's bits.
struct Found {
  std::uint8_t iteration = 0U;
  std::vector<std::uint32_t> rows;
  std::uint32_t chi2 = 0U;
  bool operator==(const Found&) const = default;
};

/// The tracks of two vertex-list iterations, at 1 and 0.3 GeV, seeded from
/// \a vertices.
std::vector<Found> trackFrom(const Fixture& fixture, std::span<const VertexSeed> vertices)
{
  auto first = iterationSettings();
  first.vertexSource = VertexSource::kVertexList;
  auto second = first;
  second.windows.minPtGeV = 0.3;
  const std::vector<IterationSettings> settings{first, second};
  Tracker<Det, UniformMaterial> tracker{settings, Fixture::kBeamConstraint, noScattering(), std::make_shared<BoundedMemoryResource>(), 2U, 1};
  static const UniformMaterial none{};
  tracker.setTrackModel(trackingGeometry(), none, Propagator{});
  const auto result = tracker.run(*fixture.frame, tracker.planTimeframe(*fixture.frame, vertices, 3, 0.));
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  std::vector<Found> tracks;
  tracks.reserve(result.tracks->size());
  const auto view = result.tracks->view();
  for (std::size_t i = 0; i < result.tracks->size(); ++i) {
    tracks.push_back({.iteration = view.iteration[i],
                      .rows = {view.clusters[i].index, view.clusters[i].index + view.count[i]},
                      .chi2 = std::bit_cast<std::uint32_t>(view.chi2[i])});
  }
  return tracks;
}

} // namespace

/// After `endTimeframe` the pass's pool holds what it held before `find`, the
/// vertexer's work arrays included; the result is the caller's, off the pool.
BOOST_AUTO_TEST_CASE(EndTimeframeGivesTheMemoryBack)
{
  const Fixture fixture{hits(), kDurationNS};
  for (const int threads : {1, 4}) {
    BOOST_TEST_CONTEXT("threads " << threads)
    {
      const auto memory = std::make_shared<BoundedMemoryResource>();
      Pass pass{passSettings(), VertexerSettings{}, Fixture::kBeamConstraint, noScattering(), memory, 2U, threads};
      const auto before = memory->getUsedMemory();
      for (int timeframe = 0; timeframe < 2; ++timeframe) {
        const auto result = pass.find(*fixture.frame, kSliceBlocks, 0);
        BOOST_REQUIRE(result.status == StageStatus::kComplete);
        BOOST_TEST(result.found.vertices.size() == kCollisions);
        BOOST_TEST(memory->getUsedMemory() > before);
        BOOST_TEST(memory->getPeakMemory() > before);
        pass.endTimeframe();
        BOOST_TEST(memory->getUsedMemory() == before);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(OneVertexPerCollision)
{
  const Fixture fixture{hits(), kDurationNS};
  const auto result = findVertices(fixture, 1);
  const auto& found = result.found;
  // Every track's one triplet on layers 0-1-2 is a candidate; a combinatorial
  // one besides them fails the candidates' cuts.
  BOOST_TEST(found.triplets > kTracks);
  BOOST_TEST(found.candidates == kTracks);
  BOOST_REQUIRE(found.vertices.size() == kCollisions);
  for (std::size_t collision = 0; collision < kCollisions; ++collision) {
    BOOST_TEST_CONTEXT("collision " << collision)
    {
      const auto& vertex = found.vertices[collision];
      BOOST_TEST(found.contributors[collision] == kTracksPerCollision);
      BOOST_TEST(std::abs(vertex.z) < 5.e-4F);
      // Its box holds the collision's tick: the vertices come in time order.
      const auto startPS = static_cast<std::uint32_t>(kStartNS[collision] * 1000.);
      BOOST_TEST(vertex.timePS <= startPS);
      BOOST_TEST(startPS < vertex.timePS + vertex.widthPS);
    }
  }
}

BOOST_AUTO_TEST_CASE(TheListDoesNotDependOnTheThreads)
{
  const Fixture fixture{hits(), kDurationNS};
  const auto one = findVertices(fixture, 1).found;
  const auto four = findVertices(fixture, 4).found;
  BOOST_REQUIRE(one.vertices.size() == four.vertices.size());
  for (std::size_t v = 0; v < one.vertices.size(); ++v) {
    BOOST_TEST(one.vertices[v].timePS == four.vertices[v].timePS);
    BOOST_TEST(one.vertices[v].widthPS == four.vertices[v].widthPS);
    BOOST_TEST(std::bit_cast<std::uint32_t>(one.vertices[v].z) == std::bit_cast<std::uint32_t>(four.vertices[v].z));
    BOOST_TEST(std::bit_cast<std::uint32_t>(one.vertices[v].varZ) == std::bit_cast<std::uint32_t>(four.vertices[v].varZ));
    BOOST_TEST(one.contributors[v] == four.contributors[v]);
  }
}

/// Each slice finds its core's frames' vertices: the list and the claims
/// are the timeframe's in one slice, at any width.
BOOST_AUTO_TEST_CASE(TheListDoesNotDependOnTheSlices)
{
  const Fixture fixture{hits(), kDurationNS};
  const auto find = [&](std::int32_t width) {
    Pass pass{passSettings(), VertexerSettings{}, Fixture::kBeamConstraint, noScattering(), std::make_shared<BoundedMemoryResource>(), 2U, 2};
    auto result = pass.find(*fixture.frame, width, 0, true);
    BOOST_REQUIRE(result.status == StageStatus::kComplete);
    std::ranges::sort(result.claims, {}, [](const TripletClaim& claim) { return std::pair{claim.rows, claim.vertex}; });
    return result;
  };
  const auto whole = find(fixture.geometry.blockCount());
  BOOST_REQUIRE(whole.found.vertices.size() == kCollisions);
  for (const std::int32_t width : {1, 2, kSliceBlocks}) {
    BOOST_TEST_CONTEXT("slices of " << width << " blocks")
    {
      const auto sliced = find(width);
      BOOST_TEST(sliced.found.candidates == whole.found.candidates);
      BOOST_REQUIRE(sliced.found.vertices.size() == whole.found.vertices.size());
      for (std::size_t v = 0; v < whole.found.vertices.size(); ++v) {
        BOOST_TEST(sliced.found.vertices[v].timePS == whole.found.vertices[v].timePS);
        BOOST_TEST(sliced.found.vertices[v].widthPS == whole.found.vertices[v].widthPS);
        BOOST_TEST(std::bit_cast<std::uint32_t>(sliced.found.vertices[v].z) == std::bit_cast<std::uint32_t>(whole.found.vertices[v].z));
        BOOST_TEST(std::bit_cast<std::uint32_t>(sliced.found.vertices[v].varZ) == std::bit_cast<std::uint32_t>(whole.found.vertices[v].varZ));
        BOOST_TEST(sliced.found.contributors[v] == whole.found.contributors[v]);
      }
      BOOST_TEST(
        std::ranges::equal(sliced.claims, whole.claims, [](const TripletClaim& a, const TripletClaim& b) { return a.rows == b.rows && a.vertex == b.vertex; }));
    }
  }
}

/// The list seeds the vertex-list iterations: they find every track, the
/// same tracks as from the collisions' true vertices.
BOOST_AUTO_TEST_CASE(ItsListSeedsTheVertexListIterations)
{
  const Fixture fixture{hits(), kDurationNS};
  const auto found = findVertices(fixture, 1).found;
  std::vector<VertexSeed> truth;
  truth.reserve(kCollisions);
  for (const double startNS : kStartNS) {
    truth.push_back({.timePS = static_cast<std::uint32_t>(startNS * 1000.), .widthPS = 0U, .z = 0.F, .varZ = 1.e-4F, .varTransverse = 1.e-6F});
  }
  const auto fromPass = trackFrom(fixture, found.vertices);
  BOOST_TEST(fromPass.size() == kTracks);
  BOOST_TEST((fromPass == trackFrom(fixture, truth)));
}

} // namespace o2::gtrk

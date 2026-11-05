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

/// \file testTracker.cxx
/// \brief One CPU tracking iteration over a core and its halo.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Tracker
#include <boost/test/unit_test.hpp>

#include "ToyBarrel3.h"
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/MeasurementGridBuilder.h"
#include "GenericTracking/Tracking/Tracker.h"

#include "GridTestSupport.h"
#include "StageOutputComparison.h"
#include "SerialReferences.h"
#include "UniformMaterial.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

namespace o2::gtrk
{

namespace
{

using Det = detectors::ToyBarrel3;
using Frame = TimeFrame<Det>;
using Scratch = TimeFrameScratch<Det>;
using Sources = LayerRowMap<Det>;
using Output = StageOutput<TrackletT>;
using Controller = Tracker<Det, UniformMaterial>;
constexpr std::size_t kMeasurements = 4U;

struct InputPoint {
  LayerID layer = kInvalidLayerID;
  std::uint32_t timePS = 0U; ///< the time box
  std::uint32_t widthPS = 0U;
  Scalar r = 0.F;
};

SearchWindowParams windowParameters()
{
  SearchWindowParams result;
  result.minPtGeV = 1.;
  result.bzKiloGauss = 20.;
  result.positionSigma = 0.01;
  result.nSigma = 3.;
  return result;
}

GridParams gridParameters()
{
  GridParams result;
  result.blockWidth = 20'000;
  return test::withToyFrames(result);
}

class Fixture
{
 public:
  static constexpr VertexConstraint kBeamConstraint{.z = 0., .sigmaZ = 0.5, .sigmaTransverse = 0.001};

  Fixture()
    : geometry(test::cappedGrid<Det>(60., gridParameters(), windowParameters(), kBeamConstraint, 1, 1)),
      offsets(static_cast<std::size_t>(geometry.bucketCount()) + 1U)
  {
    // Rows 0 and 1 form a valid halo-only pair before the core and so inside
    // the slice gate. Row 2 is the only core source; its box ends where row 1's
    // ends, so they do not pair. Row 3 is a halo target whose wide box reaches
    // back into row 2's.
    constexpr std::array<InputPoint, kMeasurements> input{{
      {.layer = 0, .timePS = 0U, .widthPS = 20'000U, .r = 2.F},
      {.layer = 1, .timePS = 5'000U, .widthPS = 15'000U, .r = 4.F},
      {.layer = 0, .timePS = 20'100U, .widthPS = 24'900U, .r = 2.F},
      {.layer = 1, .timePS = 41'000U, .widthPS = 29'000U, .r = 4.F},
    }};

    frame = test::frameOf(geometry, input.size(), [&](std::size_t i) { return test::FramePoint{.layer = input[i].layer, .x = input[i].r, .r = input[i].r}; }, [&](std::size_t i) { return TimeEstPS{.start = input[i].timePS, .width = input[i].widthPS}; }, [&](std::size_t i) { return test::PointVariances{1e-4F, 1e-4F, 1e-4F}; }, columns, offsets);
    used.assign(frame->size(), kFree);
  }

  [[nodiscard]] std::uint32_t row(std::uint32_t cluster) const
  {
    return static_cast<std::uint32_t>(std::ranges::find(columns.clusterIndex, cluster) - columns.clusterIndex.begin());
  }

  GridGeometry<Det> geometry;
  std::vector<std::uint32_t> offsets;
  MeasurementStore columns;
  std::unique_ptr<Frame> frame;
  std::vector<std::uint8_t> used; ///< the used mask: every row free
  VertexSeed beamLine = vertexSeed(kBeamConstraint);
};

std::array<std::uint32_t, Sources::kLayerCount> candidateWidthBounds(const Fixture& fixture, const Frame::TimeSlice& slice)
{
  std::array<std::uint32_t, Sources::kLayerCount> result{};
  for (std::size_t layer = 0; layer < result.size(); ++layer) {
    for (auto row = slice.halo[layer].begin; row < slice.halo[layer].end; ++row) {
      result[layer] = std::max(result[layer], fixture.frame->measurement(row).widthPS);
    }
  }
  return result;
}

/// The toy's three layers make triplets and no road: its chains are 3
/// clusters, or no link would be taken (`takenLinks`).
IterationSettings iterationSettings() { return {.windows = windowParameters(), .vertexSource = VertexSource::kBeamLine, .minClusters = 3}; }

} // namespace

BOOST_AUTO_TEST_CASE(OneSliceMatchesSerialWithHaloSourcesInsideTheGate)
{
  const Fixture fixture;
  const Frame::TimeSlice slice{*fixture.frame, {.begin = 1, .end = 2}, 1'000};
  const auto coreSource = fixture.row(2U);
  const auto haloSource = fixture.row(0U);
  const auto haloTarget = fixture.row(3U);
  BOOST_REQUIRE(slice.core[0].begin <= coreSource && coreSource < slice.core[0].end);
  BOOST_REQUIRE(haloSource < slice.core[0].begin || haloSource >= slice.core[0].end);
  BOOST_REQUIRE(slice.halo[1].begin <= haloTarget && haloTarget < slice.halo[1].end);

  const std::array<double, Sources::kLayerCount> material{};
  const auto scattering = scatteringTable<Det>(material);
  const auto memory = std::make_shared<BoundedMemoryResource>();
  Controller tracker{iterationSettings(), Fixture::kBeamConstraint, scattering, memory, 1U};
  BOOST_REQUIRE(tracker.beginTimeframe(*fixture.frame) == StageStatus::kComplete);
  const auto result = test::runSlice(tracker, 0U, *fixture.frame, slice, {});
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  BOOST_REQUIRE(result.tracklets != nullptr);
  BOOST_REQUIRE(result.tracklets->size() > 0U);
  BOOST_TEST(tracker.trackletHighWater() == result.tracklets->size());
  BOOST_TEST(tracker.stagingHighWater() > 0U);

  const auto maxCandidateWidthPS = candidateWidthBounds(fixture, slice);
  BOOST_TEST(maxCandidateWidthPS[1] == 29'000U);
  const auto parameters = iterationParameters<Det>(iterationSettings(), Fixture::kBeamConstraint, scattering);
  const std::span<const VertexSeed> noVertices;
  const auto vertexWindows = vertexWindowParameters<Det>(parameters, iterationSettings().windows, noVertices);
  const Sources sources = Scratch::makeSources(slice.halo);
  // Both source layers contribute their halo rows. Layer 1 cannot emit here
  // because the fixture has no layer-2 rows; every source lies within the
  // gate, so the ungated serial stage is the reference.
  BOOST_REQUIRE(sources.size() == 4U);
  BOOST_TEST(sources.contains(0, haloSource));
  BOOST_TEST(sources.contains(1, fixture.row(1U)));
  BOOST_TEST(sources.contains(1, fixture.row(3U)));
  const TrackletStage<Det> stage{parameters, sources, *fixture.frame, slice.haloBlocks, fixture.beamLine, noVertices,
                                 maxCandidateWidthPS, vertexWindows, fixture.used.data()};
  Output serial;
  serial.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  serial.prepare(sources.size());
  BOOST_REQUIRE(runSerial(serial, stage) == StageStatus::kComplete);
  sortTrackletSegments(serial); // the tracker hands out each source's tracklets sorted
  BOOST_TEST(sameStageOutput(*result.tracklets, serial));

  BOOST_REQUIRE(result.tracklets->size() == 2U);
  const auto tracklets = result.tracklets->view();
  const bool haloFirst = tracklets[0].source == haloSource;
  const auto& fromHalo = tracklets[haloFirst ? 0 : 1];
  const auto& fromCore = tracklets[haloFirst ? 1 : 0];
  BOOST_TEST(fromHalo.source == haloSource);
  BOOST_TEST(fromHalo.target == fixture.row(1U));
  BOOST_TEST(fromCore.source == coreSource);
  BOOST_TEST(fromCore.target == haloTarget);
}

/// After `endTimeframe` the pool holds what it held before the timeframe
/// began, at one thread and at several, timeframe after timeframe.
BOOST_AUTO_TEST_CASE(EndTimeframeGivesTheMemoryBack)
{
  const Fixture fixture;
  const Frame::TimeSlice slice{*fixture.frame, {.begin = 1, .end = 2}, 1'000};
  const std::array<double, Sources::kLayerCount> material{};
  const auto scattering = scatteringTable<Det>(material);
  for (const int threads : {1, 4}) {
    BOOST_TEST_CONTEXT("threads " << threads)
    {
      const auto memory = std::make_shared<BoundedMemoryResource>();
      Controller tracker{iterationSettings(), Fixture::kBeamConstraint, scattering, memory, 1U, threads};
      const auto before = memory->getUsedMemory();
      for (int timeframe = 0; timeframe < 3; ++timeframe) {
        BOOST_REQUIRE(tracker.beginTimeframe(*fixture.frame) == StageStatus::kComplete);
        const auto result = test::runSlice(tracker, 0U, *fixture.frame, slice, {});
        BOOST_REQUIRE(result.status == StageStatus::kComplete);
        BOOST_TEST(result.tracklets->size() > 0U);
        BOOST_TEST(memory->getUsedMemory() > before);
        tracker.endTimeframe();
        BOOST_TEST(memory->getUsedMemory() == before);
      }
    }
  }
}

} // namespace o2::gtrk

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

/// \file testStageFunnel.cxx
/// \brief The stage funnel's chains, first losses and seeded iterations
///        (StageFunnel.h).

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK stage funnel
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Study/StageFunnel.h"

#include "ToyBarrel8Fixture.h"
#include "UniformMaterial.h"

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

namespace o2::gtrk::study
{

namespace
{
using Controller = toybarrel8::Controller;
using Funnel = StageFunnel<toybarrel8::Det, Controller>;

std::unique_ptr<Funnel> runFunnel(const toybarrel8::Fixture& fixture, std::span<const TruthRow> rows, const ParticleTable& truth, const std::vector<IterationSettings>& settings)
{
  using namespace toybarrel8;
  Controller tracker{settings, Fixture::kBeamConstraint, scatteringTable<Det>(std::array<double, kLayers>{}), std::make_shared<BoundedMemoryResource>(), 2U, 1};
  const UniformMaterial none{};
  tracker.setTrackModel(trackingGeometry(), none, Propagator{});
  std::vector<double> floors;
  floors.reserve(settings.size());
  for (const auto& setting : settings) {
    floors.push_back(setting.windows.minPtGeV);
  }
  auto funnel = std::make_unique<Funnel>(*fixture.frame, rows, truth, comparisonSettings<Det>(), vertexSeed(Fixture::kBeamConstraint), floors);
  tracker.setSliceObserver([&](const Controller& observed, const Controller::SliceObservation& slice) { funnel->observeSlice(observed, slice); });
  tracker.setIterationObserver([&](const Controller& observed, const Controller::IterationObservation& iteration) { funnel->observeIteration(observed, iteration); });
  const auto result = tracker.run(*fixture.frame, tracker.planTimeframe(*fixture.frame, {}, fixture.geometry.blockCount(), 0.));
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  return funnel;
}

std::vector<TruthRow> truthRows(const toybarrel8::Fixture& fixture)
{
  std::vector<TruthRow> rows(fixture.frame->size());
  for (std::uint32_t row = 0; row < rows.size(); ++row) {
    const auto hit = fixture.frame->measurement(row).clusterIndex;
    rows[row].layer = static_cast<LayerID>(hit % toybarrel8::kLayers);
    rows[row].particles[rows[row].count++] = (hit / toybarrel8::kLayers) + 1U;
  }
  return rows;
}
} // namespace

BOOST_AUTO_TEST_CASE(ChainsAreCountedInNodesWhateverTheEdgeOrder)
{
  // Rows 0-3 on layers 0-3, a branch 1 -> 4 on layer 2; listed outermost first.
  const std::array<int, 5> layer{0, 1, 2, 3, 2};
  const std::vector<std::pair<std::uint32_t, std::uint32_t>> edges{{2, 3}, {1, 4}, {1, 2}, {0, 1}};
  const auto layerOf = [&](std::uint32_t node) { return layer[node]; };
  BOOST_TEST(longestChain(edges, layerOf) == 4U);
  const auto lengths = chainLengths(edges, layerOf);
  BOOST_TEST(lengths.at(4) == 3U);
  BOOST_TEST(lengths.at(3) == 4U);
  BOOST_TEST(!lengths.contains(0));
  BOOST_TEST(longestChain({}, layerOf) == 0U);
}

BOOST_AUTO_TEST_CASE(AParticleIsLostAtTheFirstStageShortOfTheMinimum)
{
  BOOST_TEST((firstLoss({3, 3, 3, 3, 3}, 4, false) == FunnelStage::kTracklets));
  BOOST_TEST((firstLoss({7, 7, 7, 3, 0}, 4, false) == FunnelStage::kRoads));
  BOOST_TEST((firstLoss({7, 7, 7, 7, 4}, 4, false) == FunnelStage::kSelection));
  BOOST_TEST((firstLoss({0, 0, 0, 0, 0}, 4, true) == FunnelStage::kFound));
}

BOOST_AUTO_TEST_CASE(EachParticleIsSeededByTheIterationsWithARoadOfItsOwn)
{
  using namespace toybarrel8;
  // Particle 1 at 2 GeV, above both iterations' pT_min; particle 2 at
  // 0.45 GeV, below the first's 1 GeV, above the second's 0.3.
  std::vector<Hit> hits = sensorHelixHits(0.3, 0.03, 2., 37.);
  const auto soft = sensorHelixHits(1.2, 0.01, 0.45, 120.);
  BOOST_REQUIRE(hits.size() == kLayers);
  BOOST_REQUIRE(soft.size() == kLayers);
  hits.insert(hits.end(), soft.begin(), soft.end());
  const Fixture fixture{hits, 300.};
  // Each row's particle: its hit's track.
  const auto rows = truthRows(fixture);
  const ParticleTable truth{{1, {.pt = 2.F, .primary = true, .origin = ParticleOrigin::kPhysicalPrimary}}, {2, {.pt = 0.45F, .primary = true, .origin = ParticleOrigin::kPhysicalPrimary}}};
  auto first = iterationSettings();
  auto second = iterationSettings();
  second.windows.minPtGeV = 0.3;
  const std::vector<IterationSettings> settings{first, second};

  const auto funnel = runFunnel(fixture, rows, truth, settings);

  const auto seeded = funnel->seededBy();
  BOOST_REQUIRE(seeded.size() == 2U);
  BOOST_TEST(seeded.at(1) == 0b01U); // its measurements used after iteration 0
  BOOST_TEST(seeded.at(2) == 0b10U);
}

BOOST_AUTO_TEST_CASE(APairCountsAsLostOnlyIfNoIterationKeptIt)
{
  using namespace toybarrel8;
  // Particle 1 at 2 GeV, found by iteration 0, its rows used in iteration 1;
  // particle 2 at 0.45 GeV, below iteration 0's 1 GeV, found by iteration 1;
  // particle 3 at 0.45 GeV with its layer-3 measurement 100 ns late, so that
  // its pairs 2 -> 3 and 3 -> 4 fail the time test wherever they are tried.
  std::vector<Hit> hits = sensorHelixHits(0.3, 0.03, 2., 37.);
  const auto soft = sensorHelixHits(1.2, 0.01, 0.45, 120.);
  auto late = sensorHelixHits(2.4, -0.02, 0.45, 60.);
  BOOST_REQUIRE(hits.size() == kLayers);
  BOOST_REQUIRE(soft.size() == kLayers);
  BOOST_REQUIRE(late.size() == kLayers);
  late[3].timeNS += 100.;
  hits.insert(hits.end(), soft.begin(), soft.end());
  hits.insert(hits.end(), late.begin(), late.end());
  const Fixture fixture{hits, 300.};
  const auto rows = truthRows(fixture);
  const ParticleTable truth{{1, {.pt = 2.F, .primary = true, .origin = ParticleOrigin::kPhysicalPrimary}}, {2, {.pt = 0.45F, .primary = true, .origin = ParticleOrigin::kPhysicalPrimary}}, {3, {.pt = 0.45F, .primary = true, .origin = ParticleOrigin::kPhysicalPrimary}}};
  auto first = iterationSettings();
  auto second = iterationSettings();
  second.windows.minPtGeV = 0.3;
  const std::vector<IterationSettings> settings{first, second};

  const auto funnel = runFunnel(fixture, rows, truth, settings);

  const auto one = funnel->particleVerdict(1);
  const auto two = funnel->particleVerdict(2);
  const auto three = funnel->particleVerdict(3);
  BOOST_REQUIRE((one.has_value() && two.has_value() && three.has_value()));
  // NOLINTBEGIN(bugprone-unchecked-optional-access): BOOST_REQUIRE checks them, which the check does not model
  // Kept by iteration 0: iteration 1's "used" does not count.
  BOOST_TEST(!one->pairLoss.has_value());
  BOOST_TEST((one->lostAt == FunnelStage::kSelection)); // as if not found: every stage reached
  // Iteration 0 is above its pT, iteration 1 kept its pairs.
  BOOST_TEST(!two->pairLoss.has_value());
  BOOST_REQUIRE(three->pairLoss.has_value());
  BOOST_TEST((*three->pairLoss == PairLoss::kTime));
  // NOLINTEND(bugprone-unchecked-optional-access)
  BOOST_TEST(!funnel->particleVerdict(4).has_value());
}

BOOST_AUTO_TEST_CASE(AChainRaisedFromBelowIsLostAtTheRoads)
{
  using namespace toybarrel8;
  // One helix at 2 GeV whose layer-0 measurement belongs to another particle
  // (9, not findable): its triplet on layers 0-2 is a predecessor of the
  // particle's chain on layers 1-7, which it raises by a level, so the only
  // road through the chain carries the other's measurement.
  const auto hits = sensorHelixHits(0.3, 0.03, 2., 37.);
  BOOST_REQUIRE(hits.size() == kLayers);
  const Fixture fixture{hits, 300.};
  std::vector<TruthRow> rows(hits.size());
  for (std::uint32_t row = 0; row < rows.size(); ++row) {
    const auto hit = fixture.frame->measurement(row).clusterIndex;
    rows[row].layer = static_cast<LayerID>(hit % kLayers);
    rows[row].particles[rows[row].count++] = rows[row].layer == 0 ? 9U : 1U;
  }
  const ParticleTable truth{{1, {.pt = 2.F, .primary = true, .origin = ParticleOrigin::kPhysicalPrimary}}};
  const auto first = iterationSettings();
  const std::vector<IterationSettings> settings{first};

  const auto funnel = runFunnel(fixture, rows, truth, settings);

  const auto verdict = funnel->particleVerdict(1);
  BOOST_REQUIRE(verdict.has_value());
  // NOLINTBEGIN(bugprone-unchecked-optional-access): BOOST_REQUIRE checks them, which the check does not model
  BOOST_TEST((verdict->lostAt == FunnelStage::kRoads));
  BOOST_REQUIRE(verdict->roadLoss.has_value());
  BOOST_TEST((*verdict->roadLoss == RoadLoss::kRaisedBelow));
  BOOST_TEST(!verdict->raisedFit.has_value()); // no road fit set: nothing fitted
  // NOLINTEND(bugprone-unchecked-optional-access)
}

} // namespace o2::gtrk::study

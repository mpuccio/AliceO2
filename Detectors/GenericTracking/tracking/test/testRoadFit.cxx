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

/// \file testRoadFit.cxx
/// \brief The road fit: the outward pass with its early stop against the full
///        fit, its chi2 against its degrees of freedom, and the stage's
///        selection records, on ToyBarrel8's sensors.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK RoadFit
#include <boost/test/unit_test.hpp>
#include <numbers>

#include "GenericTracking/Study/Validation.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/RoadCount.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/Selection.h"
#include "GenericTracking/Tracking/TrackRefit.h"

#include "ToyBarrel8Fixture.h"
#include "SerialReferences.h"
#include "FixedPointReferences.h"
#include "UniformMaterial.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <random>
#include <span>
#include <vector>

namespace o2::gtrk
{

namespace
{

using namespace toybarrel8;
using Stage = RoadFitStage<Det, UniformMaterial, false, false, false, false>;
using RecordOutput = StageOutput<SelectionRecordOf<maxChainLayers<Det>>::template RecordT>;
constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
/// A road cut that passes everything, per degree of freedom.
const std::array<Scalar, (3U * kMaxClusters) - 5U> kNoCut = [] {
  std::array<Scalar, (3U * kMaxClusters) - 5U> cut{};
  cut.fill(std::numeric_limits<Scalar>::infinity());
  return cut;
}();

/// A track's direction, dip slope and pT: well separated azimuths.
struct Helix {
  double phi0 = 0.;
  double slope = 0.;
  double ptGeV = 0.;
};

Helix helix(std::size_t track)
{
  const auto k = static_cast<double>(track);
  return {
    .phi0 = std::fmod(0.4 + (2.3 * k), 2 * std::numbers::pi_v<double>),
    .slope = 0.04 * std::sin(1.7 * k),
    .ptGeV = 1.5 + (0.25 * static_cast<double>(track % 5U)),
  };
}

/// The tracks' hits, one per layer each.
std::vector<Hit> trackHits(std::size_t tracks)
{
  std::vector<Hit> hits;
  for (std::size_t track = 0; track < tracks; ++track) {
    const auto h = helix(track);
    const auto one = sensorHelixHits(h.phi0, h.slope, h.ptGeV);
    BOOST_REQUIRE(one.size() == kLayers);
    hits.insert(hits.end(), one.begin(), one.end());
  }
  return hits;
}

/// \a track's rows, innermost first.
std::array<std::uint32_t, kLayers> trackRows(const Fixture& fixture, std::size_t track)
{
  std::array<std::uint32_t, kLayers> rows{};
  for (std::size_t layer = 0; layer < kLayers; ++layer) {
    rows[layer] = fixture.row((kLayers * track) + layer);
  }
  return rows;
}

const Propagator kPropagator{};
const UniformMaterial kNoMaterial{};

IterationParameters<Det> parameters()
{
  const std::array<double, kLayers> material{};
  return iterationParameters<Det>(iterationSettings(), Fixture::kBeamConstraint, scatteringTable<Det>(std::array<double, kLayers>{}));
}

} // namespace

/// The beam-line update: a primary's road on
/// the outer four layers, its inward state updated with the beam line, has a
/// chi2 near zero -- the hits are exact -- and a narrower q/p. The beam line
/// moved across the track's direction raises the chi2 with the square of the
/// distance, and 5 mm, about ten of the road's d0 sigma here, refuses it.
BOOST_AUTO_TEST_CASE(TheBeamLineConstrainsAPrimarysOuterRoad)
{
  const Fixture fixture{trackHits(6U)};
  const auto iteration = parameters();
  const auto& tracking = trackingGeometry();
  for (std::size_t track = 0; track < 6U; ++track) {
    const auto all = trackRows(fixture, track);
    const std::array<std::uint32_t, 4> rows{all[4], all[5], all[6], all[7]};
    const auto fit = fitMeasurements(kPropagator, tracking, kNoMaterial, RowChain<Det>{.frame = fixture.frame.get(), .rows = rows.data()}, 4);
    BOOST_REQUIRE(fit.status == PropagationStatus::kOk);
    const auto innermost = pairPointAt(fixture.frame->measurements(), rows[0]);
    const auto layer = tracking.surface(fit.state.surface).layer;
    const auto& plane = tracking.planeOf(fit.state);
    auto constrained = fit.state;
    const auto update = updateWithBeamLine(kPropagator, iteration, plane, layer, innermost, constrained);
    BOOST_REQUIRE(update.ok);
    BOOST_TEST(update.chi2 < 0.01F);
    BOOST_TEST(constrained.variance(kQOverP) < 0.5F * fit.state.variance(kQOverP));
    BOOST_TEST(std::abs(constrained.qOverP() - fit.state.qOverP()) < 3.F * std::sqrt(constrained.variance(kQOverP)));
    const auto phi0 = static_cast<Scalar>(helix(track).phi0);
    const auto movedBy = [&](Scalar distance) {
      auto moved = iteration;
      moved.beamX = -distance * std::sin(phi0);
      moved.beamY = distance * std::cos(phi0);
      auto state = fit.state;
      const auto result = updateWithBeamLine(kPropagator, moved, plane, layer, innermost, state);
      BOOST_REQUIRE(result.ok);
      return result.chi2;
    };
    const Scalar one = movedBy(0.1F);
    BOOST_TEST(movedBy(0.2F) == 4.F * one, boost::test_tools::tolerance(0.02F));
    BOOST_TEST(movedBy(0.5F) > iteration.beamChi2Max);
  }
}

BOOST_AUTO_TEST_CASE(TheOutwardPassStopsWhereTheFullOneWouldReject)
{
  // Track 3's layer-4 hit moved by 10 sigma of loc0 along its sensor.
  auto hits = trackHits(6U);
  auto& moved = hits[(kLayers * 3U) + 4U];
  moved.loc0 += 10.F * std::sqrt(3.e-6F * static_cast<Scalar>(1U + ((kLayers * 3U + 4U) % 3U)));
  Scalar global[3]{};
  surfaceToGlobal(trackingGeometry().surface(moved.surface), moved.loc0, moved.loc1, global);
  moved.point.x = global[0];
  moved.point.y = global[1];
  moved.point.z = global[2];
  const Fixture fixture{hits};
  const auto iteration = parameters();
  const auto cut = study::roadChi2Cut(iteration, kLayers);
  BOOST_TEST(cut > 0.F);

  for (std::size_t track = 0; track < 6U; ++track) {
    const auto rows = trackRows(fixture, track);
    const RowChain<Det> chain{.frame = fixture.frame.get(), .rows = rows.data()};
    std::array<Measurement, kLayers> copied{};
    for (std::size_t i = 0; i < kLayers; ++i) {
      copied[i] = chain[static_cast<int>(i)];
    }
    const auto full = fitMeasurements(kPropagator, trackingGeometry(), kNoMaterial, copied.data(), static_cast<int>(kLayers));
    const auto outward = fitOutward(kPropagator, trackingGeometry(), kNoMaterial, chain, static_cast<int>(kLayers), kNoCut.data());
    const auto cutAt = fitOutward(kPropagator, trackingGeometry(), kNoMaterial, chain, static_cast<int>(kLayers), iteration.roadChi2MaxByNdf);
    // A kink cannot absorb a displaced cluster: the one after it misses.
    const auto kinked = fitOutward(kPropagator, trackingGeometry(), kNoMaterial, chain, static_cast<int>(kLayers), iteration.roadChi2MaxByNdf, kinkBudget(iteration));
    BOOST_TEST_CONTEXT("track " << track)
    {
      // Read by row or copied, the same pass, to the bit.
      BOOST_REQUIRE((full.status == PropagationStatus::kOk));
      BOOST_TEST((outward.status == RoadFitStatus::kAccepted));
      BOOST_TEST(outward.chi2 == full.outwardChi2);
      BOOST_TEST((cutAt.status == RoadFitStatus::kAccepted) == (full.outwardChi2 <= cut));
      BOOST_TEST((kinked.status == cutAt.status));
      if (track == 3U) {
        BOOST_TEST((cutAt.status == RoadFitStatus::kRejected));
        BOOST_TEST(full.outwardChi2 > 25.F);
      } else {
        BOOST_TEST(kinked.kinks == 0U);
        BOOST_TEST((cutAt.status == RoadFitStatus::kAccepted));
        BOOST_TEST(outward.chi2 < 1.F);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(SmearedRoadsHaveTheirDegreesOfFreedom)
{
  // Every hit smeared by its own errors in loc0 and loc1: the outward chi2 of
  // 8 clusters has 2 * 8 - 5 = 11 degrees of freedom.
  constexpr std::size_t kTracks = 400U;
  auto hits = trackHits(kTracks);
  std::mt19937_64 generator(20260926U); // NOLINT(cert-msc32-c,cert-msc51-cpp): the pulls are compared against fixed bounds
  std::normal_distribution<double> normal;
  for (std::size_t i = 0; i < hits.size(); ++i) {
    const auto scale = static_cast<double>(1U + (i % 3U));
    auto& hit = hits[i];
    hit.loc0 += static_cast<Scalar>(normal(generator) * std::sqrt(3.e-6 * scale));
    hit.loc1 += static_cast<Scalar>(normal(generator) * std::sqrt(1.e-6 * scale));
    Scalar global[3]{};
    surfaceToGlobal(trackingGeometry().surface(hit.surface), hit.loc0, hit.loc1, global);
    hit.point.x = global[0];
    hit.point.y = global[1];
    hit.point.z = global[2];
  }
  const Fixture fixture{hits, 100.};
  const auto cut = study::roadChi2Cut(parameters(), kLayers);
  double sum = 0.;
  std::size_t accepted = 0U;
  for (std::size_t track = 0; track < kTracks; ++track) {
    const auto rows = trackRows(fixture, track);
    const RowChain<Det> chain{.frame = fixture.frame.get(), .rows = rows.data()};
    const auto fit = fitOutward(kPropagator, trackingGeometry(), kNoMaterial, chain, static_cast<int>(kLayers), kNoCut.data());
    BOOST_REQUIRE((fit.status == RoadFitStatus::kAccepted));
    sum += static_cast<double>(fit.chi2);
    accepted += fit.chi2 <= cut ? 1U : 0U;
  }
  const double mean = sum / static_cast<double>(kTracks);
  BOOST_TEST_MESSAGE("mean chi2 " << mean << " (11), within the cut " << accepted << "/" << kTracks);
  // The mean's error over 400 fits is sqrt(2 * 11 / 400) = 0.23.
  BOOST_TEST(mean > 10.F);
  BOOST_TEST(mean < 12.F);
  BOOST_TEST(accepted >= 0.99 * kTracks);
}

/// A helix whose direction turns by \a turn on the sensor of layer \a kink --
/// a hard scatter there -- as noiseless measurements: layers up to the kink
/// from the helix, the rest propagated from the kink with the new direction.
static std::array<Measurement, kLayers> kinkedChain(std::size_t kink, double turn)
{
  const auto hits = sensorHelixHits(0.7, 0.03, 1.5);
  BOOST_REQUIRE(hits.size() == kLayers);
  std::array<Measurement, kLayers> chain{};
  const auto measurementOf = [](const BoundTrackState& state, std::size_t index) {
    return Measurement{
      .loc0 = state.loc0(),
      .loc1 = state.loc1(),
      .timePS = test::toyTimeBox(test::kToyTrackTimeNS).start,
      .widthPS = test::toyTimeBox(test::kToyTrackTimeNS).width,
      .varLoc0 = 3.e-6F,
      .covLoc01 = 0.F,
      .varLoc1 = 1.e-6F,
      .clusterIndex = static_cast<std::uint32_t>(index),
      .surface = state.surface,
    };
  };
  for (std::size_t layer = 0; layer <= kink; ++layer) {
    chain[layer] = Measurement{
      .loc0 = hits[layer].loc0,
      .loc1 = hits[layer].loc1,
      .timePS = test::toyTimeBox(hits[layer].timeNS).start,
      .widthPS = test::toyTimeBox(hits[layer].timeNS).width,
      .varLoc0 = 3.e-6F,
      .covLoc01 = 0.F,
      .varLoc1 = 1.e-6F,
      .clusterIndex = static_cast<std::uint32_t>(layer),
      .surface = hits[layer].surface,
    };
  }
  // The prefix's own fit gives the helix's state on the kink's sensor.
  auto state = fitMeasurements(kPropagator, trackingGeometry(), kNoMaterial, chain.data(), static_cast<int>(kink + 1U)).outwardState;
  state.par[kPhi] += static_cast<Scalar>(turn);
  for (std::size_t layer = kink + 1U; layer < kLayers; ++layer) {
    const auto& from = trackingGeometry().surface(state.surface);
    bool crossed = false;
    for (SensorID sensor = firstSensorOfLayer<Det>(static_cast<LayerID>(layer)); sensor < firstSensorOfLayer<Det>(static_cast<LayerID>(layer + 1U)) && !crossed; ++sensor) {
      BoundTrackState probe = state;
      if (propagate(kPropagator, probe, from, trackingGeometry().surface(sensor), PropagationDirection::kForward).status == PropagationStatus::kOk) {
        state = probe;
        crossed = true;
      }
    }
    BOOST_REQUIRE(crossed);
    chain[layer] = measurementOf(state, layer);
  }
  return chain;
}

BOOST_AUTO_TEST_CASE(AKinkRecoversAScatteredRoadWithTwoDegreesOfFreedom)
{
  constexpr std::size_t kKink = 3U;
  const auto iteration = parameters();
  const auto budget = kinkBudget(iteration);
  BOOST_REQUIRE(budget.max == 1);
  const auto chain = kinkedChain(kKink, 0.05);

  // 50 mrad: the toy's layers are centimetres apart. Without a kink the road
  // is beyond its cut; with one it passes, at two degrees of freedom fewer.
  const auto strict = fitOutward(kPropagator, trackingGeometry(), kNoMaterial, chain.data(), static_cast<int>(kLayers), iteration.roadChi2MaxByNdf);
  BOOST_TEST((strict.status == RoadFitStatus::kRejected));
  const auto kinked = fitOutward(kPropagator, trackingGeometry(), kNoMaterial, chain.data(), static_cast<int>(kLayers), iteration.roadChi2MaxByNdf, budget);
  BOOST_TEST((kinked.status == RoadFitStatus::kAccepted));
  BOOST_TEST(kinked.kinks == 1U);
  BOOST_TEST(kinked.chi2 <= study::roadChi2Cut(iteration, kLayers, 1U));

  // Each pass takes its own: the outward one at the cluster after the kink,
  // the inward one at the cluster before it, where the direction it brings
  // from outside stops fitting.
  const auto fit = fitMeasurements(kPropagator, trackingGeometry(), kNoMaterial, chain.data(), static_cast<int>(kLayers), budget);
  BOOST_REQUIRE((fit.status == PropagationStatus::kOk));
  BOOST_TEST(fit.outwardKinks == 1U);
  BOOST_TEST(fit.kinks == 1U);
  BOOST_TEST(fit.kinkMeasurement == static_cast<std::int16_t>(kKink - 1U));
  BOOST_TEST(fit.ndf == static_cast<std::int16_t>((2 * static_cast<int>(kLayers)) - 5 - 2));
  BOOST_TEST(fit.outwardChi2 == kinked.chi2);
  // Not zero: both passes linearise about one helix through every cluster,
  // which a 50 mrad kink leaves; a scatter of milliradians barely moves it.
  BOOST_TEST(fit.chi2 <= study::roadChi2Cut(iteration, kLayers, fit.kinks));

  // A straight road uses none.
  const auto straight = kinkedChain(kKink, 0.);
  const auto clean = fitMeasurements(kPropagator, trackingGeometry(), kNoMaterial, straight.data(), static_cast<int>(kLayers), budget);
  BOOST_TEST(clean.outwardKinks == 0U);
  BOOST_TEST(clean.kinks == 0U);
}

BOOST_AUTO_TEST_CASE(KinksLeaveTheFitADegreeOfFreedom)
{
  // n - 3 kinks at most: 2 n - 5 - 2 k stays at least one.
  const KinkBudget two{.max = 2};
  BOOST_TEST(detail::kinkBudgetFor(two, 3).max == 0);
  BOOST_TEST(detail::kinkBudgetFor(two, 4).max == 1);
  BOOST_TEST(detail::kinkBudgetFor(two, 5).max == 2);
  BOOST_TEST(detail::kinkBudgetFor(two, 7).max == 2);
  BOOST_TEST(detail::kinkBudgetFor(KinkBudget{}, 7).max == 0);

  // The second and third clusters thrown far to either side: each inner step
  // of the outward pass is a kink while the budget lasts.
  const auto iteration = parameters();
  auto chain = kinkedChain(kLayers - 1U, 0.);
  chain[1].loc0 += 40.F * std::sqrt(chain[1].varLoc0);
  chain[2].loc0 -= 40.F * std::sqrt(chain[2].varLoc0);
  const KinkBudget budget{.max = 2, .gate = iteration.kinkChi2Min};
  const auto five = fitOutward(kPropagator, trackingGeometry(), kNoMaterial, chain.data(), 5, kNoCut.data(), budget);
  BOOST_TEST(five.kinks == 2U);
  const auto four = fitOutward(kPropagator, trackingGeometry(), kNoMaterial, chain.data(), 4, kNoCut.data(), budget);
  BOOST_TEST(four.kinks == 1U);
  const auto fit = fitMeasurements(kPropagator, trackingGeometry(), kNoMaterial, chain.data(), 4, budget);
  BOOST_REQUIRE((fit.status == PropagationStatus::kOk));
  BOOST_TEST(fit.outwardKinks <= 1U);
  BOOST_TEST(fit.kinks <= 1U);
  BOOST_TEST(fit.ndf >= 1);
}

BOOST_AUTO_TEST_CASE(AKinkCannotAbsorbAWrongLastCluster)
{
  const auto iteration = parameters();
  auto chain = kinkedChain(kLayers - 1U, 0.);
  chain[kLayers - 1U].loc0 += 30.F * std::sqrt(chain[kLayers - 1U].varLoc0);
  const auto fit = fitOutward(kPropagator, trackingGeometry(), kNoMaterial, chain.data(), static_cast<int>(kLayers), iteration.roadChi2MaxByNdf, kinkBudget(iteration));
  BOOST_TEST((fit.status == RoadFitStatus::kRejected));
  BOOST_TEST(fit.kinks == 0U);
}

BOOST_AUTO_TEST_CASE(TheLastStepTakesTheAllowanceNotAFreeDirection)
{
  // A scatter at the second-to-last cluster shows on the last step alone.
  // The refit's budget never takes it there: rejected. The road fit's takes
  // the iteration's kink allowance there (against the road-end fakes): the
  // step again with its scattering scaled, its increment
  // within the gate. A kink inside the allowance passes and frees no degree of
  // freedom; one far beyond it -- a wrong last cluster -- does not, where a
  // freed direction would take any.
  // The scattering is taken at the filtered state's momentum, the chain's
  // 1.5 GeV, not the whole-chain circle's 1.29 the kink itself bends
  // 0.014 rad is inside the allowance there.
  constexpr UniformMaterial material{.x2X0PerCm = 0.01F, .rho = 0.F};
  auto iteration = parameters();
  const auto road = [&](const std::array<Measurement, kLayers>& chain) { return takesLastStepAllowance(iteration) ? fitRoad<true>(iteration, kPropagator, trackingGeometry(), material, chain.data(), static_cast<int>(kLayers)) : fitRoad<false>(iteration, kPropagator, trackingGeometry(), material, chain.data(), static_cast<int>(kLayers)); };
  const auto inside = kinkedChain(kLayers - 2U, 0.014);
  const auto beyond = kinkedChain(kLayers - 2U, 0.08);
  const auto refitBudget = fitOutward(kPropagator, trackingGeometry(), material, inside.data(), static_cast<int>(kLayers), iteration.roadChi2MaxByNdf, kinkBudget(iteration));
  BOOST_TEST((refitBudget.status == RoadFitStatus::kRejected));
  BOOST_TEST((road(inside).status == RoadFitStatus::kRejected));
  iteration.kinkScale2 = 25.F;
  const auto on = road(inside);
  BOOST_TEST((on.status == RoadFitStatus::kAccepted));
  BOOST_TEST(on.kinks == 1U);
  BOOST_TEST(on.chi2 <= study::roadChi2Cut(iteration, kLayers, 0U));
  BOOST_TEST((road(beyond).status == RoadFitStatus::kRejected));
  // The refit's passes keep the rule: a kink at the innermost cluster going
  // inward would leave the track's direction there unchecked.
  const auto refit = fitMeasurements(kPropagator, trackingGeometry(), material, inside.data(), static_cast<int>(kLayers), kinkBudget(iteration));
  BOOST_REQUIRE((refit.status == PropagationStatus::kOk));
  BOOST_TEST(refit.outwardKinks == 0U);
}

namespace
{
struct RoadsRun {
  Controller tracker;
  Frame::TimeSlice slice;
  Controller::Result result;
  RoadCount count;

  explicit RoadsRun(const Fixture& fixture)
    : tracker(iterationSettings(), Fixture::kBeamConstraint, scatteringTable<Det>(std::array<double, kLayers>{}), std::make_shared<BoundedMemoryResource>(), 2U, 2),
      slice(*fixture.frame, fixture.geometry.allBlocks(), 0)
  {
    BOOST_REQUIRE(tracker.beginTimeframe(*fixture.frame) == StageStatus::kComplete);
    result = test::runSlice(tracker, 0U, *fixture.frame, slice, {});
    BOOST_REQUIRE(result.status == StageStatus::kComplete);
    count.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    const auto iteration = parameters();
    BOOST_REQUIRE(fixedPointRoads(count, SerialForEach{}, iteration.linkSources, slice.core.data(), roadStartRule(iteration), *result.tracklets, *result.triplets, *result.neighbours) == StageStatus::kComplete);
  }
};
} // namespace

BOOST_AUTO_TEST_CASE(TheStageEmitsOneRecordPerRoad)
{
  constexpr std::size_t kTracks = 6U;
  const Fixture fixture{trackHits(kTracks)};
  const auto iteration = parameters();
  const RoadsRun roads{fixture};
  const auto& run = roads.result;
  const auto& count = roads.count;
  const SerialForEach serial{};
  BOOST_REQUIRE(count.roadCount() == kTracks);

  const Stage stage{iteration, *fixture.frame, trackingGeometry(), kNoMaterial, kPropagator, *run.tracklets, *run.triplets, count.walk(*run.triplets, *run.neighbours)};
  RecordOutput records;
  records.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  records.prepare(count.roadCount());
  BOOST_REQUIRE(runSerial(records, stage) == StageStatus::kComplete);
  BOOST_REQUIRE(records.size() == kTracks);

  // Each road is one track, its clusters' rows in order, innermost first.
  const auto view = records.view();
  std::vector<bool> found(kTracks, false);
  for (std::uint32_t record = 0; record < records.size(); ++record) {
    const auto& clusters = view.clusters[record];
    BOOST_TEST(view.count[record] == kLayers);
    BOOST_TEST(view.chi2[record] < 1.F);
    const auto track = fixture.frame->measurement(clusters.index[0]).clusterIndex / kLayers;
    for (std::size_t layer = 0; layer < kMaxClusters; ++layer) {
      BOOST_TEST(clusters.index[layer] == (layer < kLayers ? fixture.row((kLayers * track) + layer) : ClusterList<kMaxClusters>::kNoCluster));
    }
    found[track] = true;
  }
  BOOST_TEST(std::ranges::all_of(found, [](bool f) { return f; }));

  for (const int threads : {1, 2, 4}) {
    RecordOutput parallel;
    const auto parallelPool = std::make_shared<BoundedMemoryResource>();
    parallel.setMemoryPool(parallelPool);
    parallel.prepare(count.roadCount());
    detail::TaskArena cpuArena{threads};
    CpuStage<SelectionRecordOf<kMaxClusters>::template RecordT> cpu{cpuArena};
    cpu.setMemoryPool(parallelPool);
    BOOST_REQUIRE(cpu.run(parallel, 1U, stage) == StageStatus::kComplete);
    BOOST_REQUIRE(parallel.size() == records.size());
    const auto other = parallel.view();
    BOOST_TEST(std::memcmp(other.clusters, view.clusters, records.size() * sizeof(*view.clusters)) == 0);
    BOOST_TEST(std::memcmp(other.count, view.count, records.size() * sizeof(*view.count)) == 0);
    BOOST_TEST(std::memcmp(other.chi2, view.chi2, records.size() * sizeof(*view.chi2)) == 0);
  }
}

BOOST_AUTO_TEST_CASE(TheStageDropsRoadsBelowTheFloorOfTheirClusterCount)
{
  // The six tracks at 1.5, 1.75, 2, 2.25, 2.5 and 1.5 GeV, all on every
  // layer: a 2.1 GeV floor for that count keeps the two above it, and a floor
  // for one cluster fewer keeps them all.
  constexpr std::size_t kTracks = 6U;
  const Fixture fixture{trackHits(kTracks)};
  const RoadsRun roads{fixture};
  const auto& run = roads.result;
  const auto& count = roads.count;
  const SerialForEach serial{};
  BOOST_REQUIRE(count.roadCount() == kTracks);

  const auto tracksKept = [&](double floorAll, double floorOneShort) {
    auto settings = iterationSettings();
    // No follow: the road is the final track, and the floors act on it.
    settings.follow = false;
    settings.minPtGeVByClusters[kLayers] = floorAll;
    settings.minPtGeVByClusters[kLayers - 1U] = floorOneShort;
    const auto iteration = iterationParameters<Det>(settings, Fixture::kBeamConstraint, scatteringTable<Det>(std::array<double, kLayers>{}));
    const Stage stage{iteration, *fixture.frame, trackingGeometry(), kNoMaterial, kPropagator, *run.tracklets, *run.triplets, count.walk(*run.triplets, *run.neighbours)};
    RecordOutput records;
    records.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    records.prepare(count.roadCount());
    BOOST_REQUIRE(runSerial(records, stage) == StageStatus::kComplete);
    std::vector<std::size_t> tracks;
    tracks.reserve(records.size());
    for (std::uint32_t record = 0; record < records.size(); ++record) {
      tracks.push_back(fixture.frame->measurement(records.view().clusters[record].index[0]).clusterIndex / kLayers);
    }
    std::ranges::sort(tracks);
    return tracks;
  };
  BOOST_TEST(tracksKept(2.1, 0.) == (std::vector<std::size_t>{3U, 4U}), boost::test_tools::per_element());
  BOOST_TEST(tracksKept(0., 100.) == (std::vector<std::size_t>{0U, 1U, 2U, 3U, 4U, 5U}), boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(ARefusedRoadIsFittedAgainWithoutItsInnermostClusters)
{
  // The six tracks on every layer, refused at their full length by a floor
  // no track clears: with innerTrim each
  // road's record is its longest suffix the floors admit, never shorter
  // than minClusters; without it, or with it off, there is none.
  constexpr std::size_t kTracks = 6U;
  const Fixture fixture{trackHits(kTracks)};
  const RoadsRun roads{fixture};
  const auto& run = roads.result;
  const auto& count = roads.count;
  const SerialForEach serial{};
  BOOST_REQUIRE(count.roadCount() == kTracks);

  // The records, each its cluster count times 10 plus its first layer, by
  // \a Trim's stage.
  const auto recordsOf = [&]<bool Trim>(int innerTrim, int minClusters, std::size_t refusedBelow) {
    auto settings = iterationSettings();
    settings.innerTrim = innerTrim;
    settings.minClusters = minClusters;
    settings.follow = false;
    for (std::size_t clusters = refusedBelow; clusters <= kLayers; ++clusters) {
      settings.minPtGeVByClusters[clusters] = 100.;
    }
    const auto iteration = iterationParameters<Det>(settings, Fixture::kBeamConstraint, scatteringTable<Det>(std::array<double, kLayers>{}));
    const RoadFitStage<Det, UniformMaterial, false, false, false, Trim> stage{iteration, *fixture.frame, trackingGeometry(), kNoMaterial, kPropagator, *run.tracklets, *run.triplets, count.walk(*run.triplets, *run.neighbours)};
    RecordOutput records;
    records.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    records.prepare(count.roadCount());
    BOOST_REQUIRE(runSerial(records, stage) == StageStatus::kComplete);
    std::vector<std::uint32_t> kept;
    for (std::uint32_t record = 0; record < records.size(); ++record) {
      const auto innermost = fixture.frame->measurement(records.view().clusters[record].index[0]).clusterIndex;
      kept.push_back((10U * records.view().count[record]) + static_cast<std::uint32_t>(innermost % kLayers));
      // The clusters are the road's outermost, in order.
      const auto track = innermost / kLayers;
      for (std::uint32_t position = 0; position < records.view().count[record]; ++position) {
        BOOST_TEST(records.view().clusters[record].index[position] == fixture.row((kLayers * track) + (innermost % kLayers) + position));
      }
    }
    return kept;
  };
  const std::vector<std::uint32_t> none;
  const std::vector<std::uint32_t> seven(kTracks, 71U);
  const std::vector<std::uint32_t> five(kTracks, 53U);
  BOOST_TEST(recordsOf.template operator()<false>(3, 4, kLayers) == none, boost::test_tools::per_element());
  BOOST_TEST(recordsOf.template operator()<true>(0, 4, kLayers) == none, boost::test_tools::per_element());
  BOOST_TEST(recordsOf.template operator()<true>(1, 4, kLayers) == seven, boost::test_tools::per_element());
  BOOST_TEST(recordsOf.template operator()<true>(3, 4, kLayers - 2U) == five, boost::test_tools::per_element());
  BOOST_TEST(recordsOf.template operator()<true>(2, 4, kLayers - 2U) == none, boost::test_tools::per_element()); // two dropped are not enough
  BOOST_TEST(recordsOf.template operator()<true>(3, static_cast<int>(kLayers) - 1, kLayers - 2U) == none,
             boost::test_tools::per_element()); // nor below minClusters
  // A road the stage takes whole is taken whole.
  BOOST_TEST(recordsOf.template operator()<true>(3, 4, kLayers + 1U) == std::vector<std::uint32_t>(kTracks, 80U), boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(SelectedTracksAreRefittedWithBothPasses)
{
  // Hits to tracks: tracklets, triplets, neighbours, levels, roads, their
  // fit, the selection, the refit.
  constexpr std::size_t kTracks = 6U;
  const Fixture fixture{trackHits(kTracks)};
  const auto iteration = parameters();
  const RoadsRun roads{fixture};
  const auto& run = roads.result;
  const auto& count = roads.count;
  const SerialForEach serial{};
  const Stage stage{iteration, *fixture.frame, trackingGeometry(), kNoMaterial, kPropagator, *run.tracklets, *run.triplets, count.walk(*run.triplets, *run.neighbours)};
  RecordOutput records;
  records.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  records.prepare(count.roadCount());
  BOOST_REQUIRE(runSerial(records, stage) == StageStatus::kComplete);
  GlobalSelection<kMaxClusters> selection;
  selection.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  BOOST_REQUIRE(selection.append(records, SerialForEach{}) == StageStatus::kComplete);
  std::vector<std::uint8_t> mask(fixture.frame->size(), kFree);
  BOOST_REQUIRE(selection.select(serial, [](std::span<std::uint32_t> indices, const auto& less) { std::ranges::sort(indices, less); }, 0U, mask) == StageStatus::kComplete);
  BOOST_REQUIRE(selection.selected().size() == kTracks);

  TrackRefit<Det> refit;
  refit.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  BOOST_REQUIRE(refit.run(serial, selection, iteration, *fixture.frame, trackingGeometry(), kNoMaterial, kPropagator, 0U) == StageStatus::kComplete);
  BOOST_REQUIRE(refit.size() == kTracks);
  const auto tracks = refit.tracks().view();
  for (std::uint32_t track = 0; track < kTracks; ++track) {
    const auto& rows = tracks.clusters[track];
    const auto index = fixture.frame->measurement(rows.index[0]).clusterIndex / kLayers;
    const auto h = helix(index);
    BOOST_TEST_CONTEXT("track " << index)
    {
      // Both passes, as fitMeasurements gives them on the same rows.
      BOOST_REQUIRE(tracks.status[track] == PropagationStatus::kOk);
      BOOST_TEST(tracks.count[track] == kLayers);
      std::array<Measurement, kLayers> copied{};
      const RowChain<Det> chain{.frame = fixture.frame.get(), .rows = rows.index};
      for (std::size_t i = 0; i < kLayers; ++i) {
        copied[i] = chain[static_cast<int>(i)];
      }
      const auto full = fitMeasurements(kPropagator, trackingGeometry(), kNoMaterial, copied.data(), static_cast<int>(kLayers));
      BOOST_TEST(tracks.chi2[track] == full.chi2);

      BOOST_TEST(tracks.ndf[track] == (2 * static_cast<int>(kLayers)) - 5);
      // On the innermost surface, and the generating helix.
      const auto& state = tracks.state[track];
      const auto first = fixture.frame->measurement(rows.index[0]);
      BOOST_TEST(state.surface == first.surface);
      BOOST_TEST(std::abs(state.loc0() - first.loc0) < 1.e-4F);
      BOOST_TEST(std::abs(state.loc1() - first.loc1) < 1.e-4F);
      BOOST_TEST(state.theta() == static_cast<Scalar>(std::atan2(1., h.slope)), boost::test_tools::tolerance(1.e-4F));
      const double momentum = h.ptGeV * std::hypot(1., h.slope);
      BOOST_TEST(std::abs(state.qOverP()) == static_cast<Scalar>(1. / momentum), boost::test_tools::tolerance(1.e-3F));
    }
  }

  for (const int threads : {1, 2, 4}) {
    detail::TaskArena cpuArena{threads};
    CpuStage<TrackOf<kMaxClusters>::template RecordT> cpu{cpuArena};
    cpu.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    TrackRefit<Det> parallel;
    parallel.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    BOOST_REQUIRE(parallel.run([&](std::uint32_t elements, const auto& body) { cpu.forEachElement(elements, 1U, body); }, selection, iteration, *fixture.frame, trackingGeometry(), kNoMaterial, kPropagator, 0U) == StageStatus::kComplete);
    // Bit for bit.
    const auto other = parallel.tracks().view();
    const auto bits = [](Scalar value) { return std::bit_cast<std::uint32_t>(value); };
    for (std::uint32_t track = 0; track < kTracks; ++track) {
      const auto& a = other.state[track];
      const auto& b = tracks.state[track];
      BOOST_TEST(std::ranges::equal(a.par, b.par, {}, bits, bits));
      BOOST_TEST(std::ranges::equal(a.cov, b.cov, {}, bits, bits));
      BOOST_TEST(other.timePS[track] == tracks.timePS[track]);
      BOOST_TEST(other.widthPS[track] == tracks.widthPS[track]);
      BOOST_TEST(a.surface == b.surface);
      BOOST_TEST(bits(other.chi2[track]) == bits(tracks.chi2[track]));
      BOOST_TEST(std::ranges::equal(other.clusters[track].index, tracks.clusters[track].index));
    }
  }
}

} // namespace o2::gtrk

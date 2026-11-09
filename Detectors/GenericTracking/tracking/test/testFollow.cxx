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

/// \file testFollow.cxx
/// \brief The follow: prediction, proposal and
///        resolution of clusters beyond a track, the outward follow
///        and the sensors by place its hole
///        test reads.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Follow
#include <boost/test/unit_test.hpp>
#include <numbers>

#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Base/LayerLinks.h"
#include "ToyBarrel4Disc8.h"
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/SensorCells.h"
#include "GenericTracking/Tracking/SlicePlan.h"

#include "ToyBarrel8Fixture.h"
#include "SerialReferences.h"
#include "UniformMaterial.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <span>
#include <utility>
#include <vector>

namespace o2::gtrk
{

namespace
{

template <DetectorLayout Det>
IterationParameters<Det> parameters()
{
  auto settings = toybarrel8::iterationSettings();
  return iterationParameters<Det>(settings, VertexConstraint{}, scatteringTable<Det>(uniformLayerMaterial<Det>(0.)));
}

Scalar phiDifference(Scalar a, Scalar b) { return o2::gpu::CAMath::Abs(o2::math_utils::detail::toPMPi<Scalar>(a - b)); }

constexpr std::size_t kMaxClusters = maxChainLayers<toybarrel8::Det>;
using ToyTracks = BoundedColumns<TrackOf<kMaxClusters>::template RecordT>;
using ToyProposal = FollowProposal<kMaxClusters>;

/// One track of \a fixture's rows \a ordered, innermost first, fitted.
ToyTracks fittedRows(const toybarrel8::Fixture& fixture, const std::vector<std::uint32_t>& ordered, const std::shared_ptr<BoundedMemoryResource>& pool)
{
  ClusterList<kMaxClusters> rows;
  std::ranges::fill(rows.index, ClusterList<kMaxClusters>::kNoCluster);
  std::ranges::copy(ordered, rows.index);
  const auto count = ordered.size();
  const UniformMaterial material;
  const auto fit = fitMeasurements(Propagator{}, toybarrel8::trackingGeometry(), material, RowChain<toybarrel8::Det>{.frame = fixture.frame.get(), .rows = rows.index}, static_cast<int>(count));
  BOOST_REQUIRE(fit.status == PropagationStatus::kOk);
  ToyTracks tracks;
  tracks.setMemoryResource(pool.get());
  tracks.resize(1U);
  auto view = tracks.view();
  view.state[0] = fit.state;
  view.outwardState[0] = fit.outwardState;
  view.chi2[0] = fit.chi2;
  view.ndf[0] = fit.ndf;
  view.status[0] = fit.status;
  view.clusters[0] = rows;
  view.count[0] = static_cast<std::uint8_t>(count);
  view.iteration[0] = 0U;
  view.addedOutward[0] = 0U;
  view.addedInward[0] = 0U;
  return tracks;
}

ToyTracks fittedTrack(const toybarrel8::Fixture& fixture, std::size_t firstLayer, std::size_t lastLayer, const std::shared_ptr<BoundedMemoryResource>& pool)
{
  std::vector<std::uint32_t> rows;
  for (std::size_t layer = firstLayer; layer <= lastLayer; ++layer) {
    rows.push_back(fixture.row(layer));
  }
  return fittedRows(fixture, rows, pool);
}

std::array<std::uint32_t, toybarrel8::kLayers> widthBounds(const toybarrel8::Fixture& fixture) { return frameBounds(*fixture.frame).widthPS; }

/// The toy's sensors by place, built once.
const SensorCells<toybarrel8::Det>& toyCells()
{
  static const SensorCells<toybarrel8::Det> cells(toybarrel8::trackingGeometry());
  return cells;
}

/// The material the kink tests step through: enough that the kink
/// allowance's scaled scattering opens the step.
constexpr UniformMaterial kKinkMaterial{.x2X0PerCm = 0.002F, .rho = 0.F};

/// The kinks the budget test turns the track by beyond layer 5 and inside
/// layer 2 [rad]: each in the middle of the range the nominal step refuses
/// and the allowance admits.
constexpr double kKinkOut = 0.0175;
constexpr double kKinkIn = 0.0225;
/// The kink allowance of the budget test [sigma]: both kinks fit well
/// inside it, neither in the nominal step.
constexpr double kBudgetKinkSigma = 30.;

/// An iteration that follows, \a holes allowed, a kink allowance of
/// \a kinkSigma with \a maxKinks, the beam constraint with \a beamLine, the
/// end decided again with \a followEnd.
IterationParameters<toybarrel8::Det> followParameters(int holes = 0, double kinkSigma = 1., int maxKinks = 1, bool beamLine = false, bool followEnd = false)
{
  auto settings = toybarrel8::iterationSettings();
  settings.follow = true;
  settings.followEnd = followEnd;
  settings.followHoles = holes;
  settings.kinkSigma = kinkSigma;
  settings.maxKinks = maxKinks;
  settings.beamConstraint = beamLine;
  return iterationParameters<toybarrel8::Det>(settings, VertexConstraint{}, scatteringTable<toybarrel8::Det>(uniformLayerMaterial<toybarrel8::Det>(0.)));
}

/// Track \a track of \a tracks followed (`FollowStage`, the end decided
/// again with \a DecidesEnd), its kink steps into \a steps where given,
/// inward too with \a inward; serially, or on \a threads in chunks of
/// \a chunkSize; with \a Greedy, the smallest chi2 instead of the lookahead; with
/// \a EndDensity, the ends decided on the last row's density; with
/// \a BeamEnd, the inward end against the beam line.
template <bool BeamLine = false, bool DecidesEnd = false, bool Greedy = false, bool EndDensity = false, bool BeamEnd = false>
ToyProposal followTrack(const toybarrel8::Fixture& fixture, const ToyTracks& tracks, const std::vector<std::uint8_t>& used, const IterationParameters<toybarrel8::Det>& parameters, const UniformMaterial& material = UniformMaterial{}, follow::KinkSteps* steps = nullptr, bool inward = true, int threads = 0, std::uint32_t chunkSize = 1U)
{
  const auto widths = widthBounds(fixture);
  std::vector<ToyProposal> proposals(tracks.size());
  std::vector<follow::KinkSteps> kinks(tracks.size());
  // The end decided again starts from the fit without it, a pass of its own.
  std::vector<BoundTrackState> starts(tracks.size());
  if constexpr (DecidesEnd) {
    const follow::EndFit<toybarrel8::Det, UniformMaterial> fit{
      .parameters = &parameters,
      .frame = fixture.frame.get(),
      .geometry = &toybarrel8::trackingGeometry(),
      .material = &material,
      .propagator = Propagator{},
      .tracks = tracks.view(),
      .starts = starts.data(),
      .trackOffset = 0U,
    };
    for (std::uint32_t track = 0; track < tracks.size(); ++track) {
      fit(track);
    }
  }
  const FollowStage<toybarrel8::Det, UniformMaterial, BeamLine, DecidesEnd, Greedy, EndDensity, BeamEnd> stage{&parameters, fixture.frame.get(), &toybarrel8::trackingGeometry(), &material, Propagator{}, tracks.view(), used.data(), widths.data(), toyCells().view(), proposals.data(), kinks.data(), inward, 0U, DecidesEnd ? starts.data() : nullptr};
  if (threads == 0) {
    for (std::uint32_t track = 0; track < tracks.size(); ++track) {
      stage(track);
    }
  } else {
    detail::TaskArena cpuArena{threads};
    CpuStage<TrackletT> cpu{cpuArena};
    cpu.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    cpu.forEachElement(static_cast<std::uint32_t>(tracks.size()), chunkSize, stage);
  }
  if (steps != nullptr) {
    *steps = kinks.front();
  }
  return proposals.front();
}

bool sameProposal(const ToyProposal& a, const ToyProposal& b) { return a.outwardCount == b.outwardCount && a.inwardCount == b.inwardCount && a.replacedEnd == b.replacedEnd && std::ranges::equal(a.outward, b.outward) && std::ranges::equal(a.inward, b.inward); }

/// The follow's passes over \a tracks (`Follow::run`), serially.
StageStatus followAll(const toybarrel8::Fixture& fixture, const IterationParameters<toybarrel8::Det>& parameters, ToyTracks& tracks, std::vector<std::uint8_t>& used, Follow<toybarrel8::Det>& follow, const UniformMaterial& material = UniformMaterial{}) { return follow.run(SerialForEach{}, parameters, *fixture.frame, toybarrel8::trackingGeometry(), material, Propagator{}, tracks, 0U, true, false, used, widthBounds(fixture), toyCells().view()); }

template <DetectorLayout Det>
SensorID crossedSensor(const TrackingGeometry<Det>& geometry, const Propagator& propagator, const BoundTrackState& state, LayerID target, PropagationDirection direction)
{
  const auto& from = geometry.surface(state.surface);
  const SensorID begin = firstSensorOfLayer<Det>(target);
  const SensorID end = firstSensorOfLayer<Det>(static_cast<LayerID>(target + 1U));
  for (SensorID sensor = begin; sensor < end; ++sensor) {
    BoundTrackState probe = state;
    if (propagate(propagator, probe, from, geometry.surface(sensor), direction).status == PropagationStatus::kOk) {
      return sensor;
    }
  }
  return kInvalidSensorID;
}

template <DetectorLayout Det>
void checkAgainstSensor(const TrackingGeometry<Det>& geometry, const IterationParameters<Det>& iteration, const Propagator& propagator, BoundTrackState state, LayerID target, PropagationDirection direction)
{
  const auto source = geometry.surface(state.surface).layer;
  const auto link = direction == PropagationDirection::kForward ? linkOf<Det>(source, target) : linkOf<Det>(target, source);
  BOOST_REQUIRE(link < nLinks<Det>);
  const auto prediction = predictLayerCrossing(propagator, iteration, geometry.surface(state.surface), state, target, direction);
  BOOST_REQUIRE(prediction.status == PropagationStatus::kOk);
  const auto sensor = crossedSensor(geometry, propagator, state, target, direction);
  BOOST_REQUIRE(sensor != kInvalidSensorID);
  const auto propagated = propagate(propagator, state, geometry.surface(state.surface), geometry.surface(sensor), direction);
  BOOST_REQUIRE(propagated.status == PropagationStatus::kOk);
  Scalar global[3]{};
  surfaceToGlobal(geometry.surface(sensor), state.loc0(), state.loc1(), global);
  const auto phi = azimuth(global[0], global[1]);
  const auto second = secondCoordinate<Det>(target, o2::gpu::CAMath::Hypot(global[0], global[1]), global[2]);
  BOOST_TEST(phiDifference(prediction.phi, phi) <= iteration.windows[link].phiHalfWidth);
  BOOST_TEST(o2::gpu::CAMath::Abs(prediction.second - second) <= iteration.windows[link].secondHalfWidth);
}

/// \a hits with the ones inside layer \a layer turned there by \a angle in
/// the transverse plane: a kink at that layer.
std::vector<toybarrel8::Hit> kinkedAt(std::vector<toybarrel8::Hit> hits, LayerID layer, double angle)
{
  const double at = layerDistance<toybarrel8::Det>(layer);
  for (auto& hit : hits) {
    if (hit.point.layer < layer) {
      const double shift = angle * (at - static_cast<double>(hit.point.r));
      hit.point = toybarrel8::moved(hit.point, 0., shift);
      hit.loc0 += static_cast<Scalar>(shift);
    }
  }
  return hits;
}

/// \a hits with the ones outside layer \a layer turned there by \a angle.
std::vector<toybarrel8::Hit> kinkedBeyond(std::vector<toybarrel8::Hit> hits, LayerID layer, double angle)
{
  const double at = layerDistance<toybarrel8::Det>(layer);
  for (auto& hit : hits) {
    if (hit.point.layer > layer) {
      const double shift = angle * (static_cast<double>(hit.point.r) - at);
      hit.point = toybarrel8::moved(hit.point, 0., shift);
      hit.loc0 += static_cast<Scalar>(shift);
    }
  }
  return hits;
}

/// The first step of a follow from \a state toward \a row, as the walk
/// makes it (`updateCandidate`, nominal): its chi2 and sqrt|V| [cm^2].
std::pair<double, double> firstStep(const toybarrel8::Fixture& fixture, const BoundTrackState& state, std::uint32_t row, PropagationDirection direction, const UniformMaterial& material = UniformMaterial{})
{
  BoundTrackState updated;
  CovarianceFactors factors;
  Scalar chi2 = 0.F;
  Scalar determinant = 0.F;
  BOOST_REQUIRE((follow::updateCandidate<false, true>(Propagator{}, toybarrel8::trackingGeometry(), material, measurementAt(*fixture.frame, row), direction, state, factorCovariance(state.cov), updated, factors, chi2, 1.F, &determinant)));
  return {static_cast<double>(chi2), std::sqrt(static_cast<double>(determinant))};
}

/// Whether a query at \a surface's centre, with no reach, visits it.
template <DetectorLayout Det>
bool visitedAtItsCentre(const SensorCells<Det>& cells, const TrackingSurface& surface)
{
  const bool barrel = Detector<Det>::isBarrelLayer(surface.layer);
  Scalar centre[3]{};
  surfaceToGlobal(surface, 0.F, 0.F, centre);
  const auto phi = o2::math_utils::detail::to02Pi<Scalar>(std::atan2(centre[1], centre[0]));
  const Scalar second = barrel ? centre[2] : std::hypot(centre[0], centre[1]);
  return cells.view().anySensor(surface.layer, phi, 0.F, second, 0.F, [&](std::uint32_t id) { return id == surface.id; });
}

} // namespace

BOOST_AUTO_TEST_CASE(ToyBarrelPredictionMatchesTheHelixAndTheCrossedSensor)
{
  using namespace toybarrel8;
  constexpr double kStartPhi = 0.37;
  constexpr double kSlope = 0.04;
  constexpr double kPt = 1.7;
  const auto hits = sensorHelixHits(kStartPhi, kSlope, kPt);
  BOOST_REQUIRE(hits.size() == kLayers);
  const auto& geometry = trackingGeometry();
  const auto iteration = parameters<Det>();
  const Propagator propagator;
  const auto theta = static_cast<Scalar>(std::atan2(1., kSlope));
  const Scalar sinTheta = std::sin(theta);
  const double curvature = 1. / curvatureRadius(kPt, 20.);
  const double arc = static_cast<double>(hits[0].point.z) / kSlope;

  BoundTrackState state;
  state.par[kLoc0] = hits[0].loc0;
  state.par[kLoc1] = hits[0].loc1;
  state.par[kPhi] = static_cast<Scalar>(kStartPhi + (curvature * arc));
  state.par[kTheta] = theta;
  state.par[kQOverP] = -sinTheta / static_cast<Scalar>(kPt);
  state.surface = hits[0].surface;
  const auto prediction = predictLayerCrossing(propagator, iteration, geometry.surface(state.surface), state, 1, PropagationDirection::kForward);
  BOOST_REQUIRE(prediction.status == PropagationStatus::kOk);
  const auto nominal = helixAt(kStartPhi, kSlope, kPt, 1);
  BOOST_TEST(phiDifference(prediction.phi, nominal.phi) < 2.e-6F);
  BOOST_TEST(o2::gpu::CAMath::Abs(prediction.second - nominal.z) < 2.e-5F);
  checkAgainstSensor(geometry, iteration, propagator, state, 1, PropagationDirection::kForward);
}

BOOST_AUTO_TEST_CASE(BarrelAndDiscPredictionsMatchCrossedSensors)
{
  using Mixed = detectors::ToyBarrel4Disc8;
  const auto detector = toyGeometry<Mixed>();
  const TrackingGeometry<Mixed> geometry{detector};
  const auto iteration = parameters<Mixed>();
  const Propagator propagator;

  BoundTrackState barrel;
  barrel.surface = firstSensorOfLayer<Mixed>(0);
  const auto& barrelSurface = geometry.surface(barrel.surface);
  barrel.par[kPhi] = o2::gpu::CAMath::ATan2(barrelSurface.origin[1], barrelSurface.origin[0]);
  barrel.par[kTheta] = 1.4F;
  barrel.par[kQOverP] = 0.2F;
  checkAgainstSensor(geometry, iteration, propagator, barrel, 1, PropagationDirection::kForward);

  constexpr LayerID kDisc = 8;
  constexpr LayerID kNextDisc = 9;
  BoundTrackState disc;
  disc.surface = firstSensorOfLayer<Mixed>(kDisc);
  const auto& discSurface = geometry.surface(disc.surface);
  const Scalar radius = o2::gpu::CAMath::Hypot(discSurface.origin[0], discSurface.origin[1]);
  disc.par[kPhi] = o2::gpu::CAMath::ATan2(discSurface.origin[1], discSurface.origin[0]);
  disc.par[kTheta] = o2::gpu::CAMath::ATan2(radius, discSurface.origin[2]);
  disc.par[kQOverP] = 0.2F;
  checkAgainstSensor(geometry, iteration, propagator, disc, kNextDisc, PropagationDirection::kForward);
}

/// A track at eta ~ 0 reaches a forward disc only after looping round many
/// times. That crossing is not the next layer of its arc, and its covariance
/// is past what float can factorise, so the prediction refuses it; a track
/// that reaches the same disc within half a turn is predicted.
BOOST_AUTO_TEST_CASE(PredictionRefusesACrossingPastHalfATurn)
{
  using Mixed = detectors::ToyBarrel4Disc8;
  const auto detector = toyGeometry<Mixed>();
  const TrackingGeometry<Mixed> geometry{detector};
  const auto iteration = parameters<Mixed>();
  const Propagator propagator;

  constexpr LayerID kBarrel = 0;
  LayerID disc = kInvalidLayerID;
  for (std::size_t link = 0; link < nLinks<Mixed> && disc == kInvalidLayerID; ++link) {
    if (Mixed::layerLinks[link].from == kBarrel && !Detector<Mixed>::isBarrelLayer(Mixed::layerLinks[link].to)) {
      disc = Mixed::layerLinks[link].to;
    }
  }
  BOOST_REQUIRE(disc != kInvalidLayerID);
  const auto& discSurface = geometry.surface(firstSensorOfLayer<Mixed>(disc));
  const bool forwardZ = discSurface.origin[2] > 0.F;

  BoundTrackState state;
  state.surface = firstSensorOfLayer<Mixed>(kBarrel);
  const auto& start = geometry.surface(state.surface);
  state.par[kPhi] = o2::gpu::CAMath::ATan2(start.origin[1], start.origin[0]);
  // cot(theta) = 0.0103 and p = 0.5 GeV: a track that loops many times
  // before it reaches the disc.
  constexpr Scalar kNearlyTransverse = 1.5605F;
  state.par[kTheta] = forwardZ ? kNearlyTransverse : std::numbers::pi_v<Scalar> - kNearlyTransverse;
  state.par[kQOverP] = 2.F;
  const auto looper = predictLayerCrossing(propagator, iteration, start, state, disc, PropagationDirection::kForward);
  BOOST_TEST((looper.status == PropagationStatus::kTooManyTurns));

  state.par[kTheta] = o2::gpu::CAMath::ATan2(o2::gpu::CAMath::Hypot(discSurface.origin[0], discSurface.origin[1]), discSurface.origin[2]);
  const auto direct = predictLayerCrossing(propagator, iteration, start, state, disc, PropagationDirection::kForward);
  BOOST_TEST((direct.status == PropagationStatus::kOk));
}

/// A predicted covariance that is not positive definite can give a negative
/// or NaN chi2, or a NaN posterior. The update is refused rather than handed
/// back with a chi2 that can pass the cut, or win it when negative.
BOOST_AUTO_TEST_CASE(UpdateRefusesAnIndefiniteCovariance)
{
  using Mixed = detectors::ToyBarrel4Disc8;
  const auto detector = toyGeometry<Mixed>();
  const TrackingGeometry<Mixed> geometry{detector};
  const Propagator propagator;
  const UniformMaterial material;

  BoundTrackState state;
  state.surface = firstSensorOfLayer<Mixed>(0);
  const auto& start = geometry.surface(state.surface);
  state.par[kPhi] = o2::gpu::CAMath::ATan2(start.origin[1], start.origin[0]);
  state.par[kTheta] = 1.4F;
  state.par[kQOverP] = 0.2F;
  for (int i = 0; i < kNBoundParams; ++i) {
    state.variance(i) = 1.e-4F;
  }
  const auto sensor = crossedSensor(geometry, propagator, state, 1, PropagationDirection::kForward);
  BOOST_REQUIRE(sensor != kInvalidSensorID);
  BoundTrackState atTarget = state;
  BOOST_REQUIRE(propagate(propagator, atTarget, start, geometry.surface(sensor), PropagationDirection::kForward).status == PropagationStatus::kOk);
  const Measurement measurement{
    .loc0 = atTarget.loc0(),
    .loc1 = atTarget.loc1(),
    .timePS = test::toyTimeBox(test::kToyTrackTimeNS).start,
    .widthPS = test::toyTimeBox(test::kToyTrackTimeNS).width,
    .varLoc0 = 1.e-6F,
    .covLoc01 = 0.F,
    .varLoc1 = 1.e-6F,
    .clusterIndex = 0U,
    .surface = sensor,
  };

  BoundTrackState output;
  CovarianceFactors outputFactors;
  Scalar chi2 = 0.F;
  BOOST_REQUIRE(follow::updateCandidate<false>(propagator, geometry, material, measurement, PropagationDirection::kForward, state, factorCovariance(state.cov), output, outputFactors, chi2, 1.F));
  BOOST_TEST(chi2 >= 0.F);

  // The source may also be the destination buffer. Material must still
  // start at the original source point, before publishing the crossing.
  const UniformMaterial filled{.x2X0PerCm = 1.e-3F, .rho = 0.1F};
  Propagator massive = propagator;
  massive.massGeV = 0.13957F;
  BoundTrackState separate;
  CovarianceFactors separateFactors;
  Scalar separateChi2 = 0.F;
  BOOST_REQUIRE(follow::updateCandidate<false>(massive, geometry, filled, measurement, PropagationDirection::kForward, state, factorCovariance(state.cov), separate, separateFactors, separateChi2, 1.F));
  auto aliased = state;
  auto aliasedFactors = factorCovariance(state.cov);
  Scalar aliasedChi2 = 0.F;
  BOOST_REQUIRE(follow::updateCandidate<false>(massive, geometry, filled, measurement, PropagationDirection::kForward, aliased, aliasedFactors, aliased, aliasedFactors, aliasedChi2, 1.F));
  // Bit identity is intentional, including floating-point encodings; these records have no padding.
  // NOLINTNEXTLINE(bugprone-suspicious-memory-comparison)
  BOOST_TEST(std::memcmp(&aliased, &separate, sizeof(BoundTrackState)) == 0);
  // Bit identity is intentional, including floating-point encodings; these records have no padding.
  // NOLINTNEXTLINE(bugprone-suspicious-memory-comparison)
  BOOST_TEST(std::memcmp(&aliasedFactors, &separateFactors, sizeof(CovarianceFactors)) == 0);
  BOOST_TEST(exactlyEqual(aliasedChi2, separateChi2));

  // A negative loc0 variance, as a covariance past what float can hold ends
  // up with: its factors are not those of a covariance, and the step refuses
  // them before a negative chi2 could win the comparison.
  state.variance(kLoc0) = -1.e-2F;
  Measurement displaced = measurement;
  displaced.loc0 += 0.01F;
  BOOST_TEST(!follow::updateCandidate<false>(propagator, geometry, material, displaced, PropagationDirection::kForward, state, factorCovariance(state.cov), output, outputFactors, chi2, 1.F));
}

/// Navigation, both ways: a track on layers 0-4 takes 5-7 outward, one on
/// 3-7 takes 2-0 inward, one on 2-5 both; the same on any threads; a track
/// on every layer is not followed.
BOOST_AUTO_TEST_CASE(AFollowTakesTheRowsBeyondTheTrackBothWays)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(0.42, 0.025, 1.8);
  BOOST_REQUIRE(hits.size() == kLayers);
  const Fixture fixture{hits};
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto parameters = followParameters();

  const auto outerTracks = fittedTrack(fixture, 0U, 4U, pool);
  const auto outward = followTrack(fixture, outerTracks, used, parameters);
  BOOST_REQUIRE(outward.outwardCount == 3U);
  BOOST_TEST(outward.inwardCount == 0U);
  for (std::size_t i = 0; i < 3U; ++i) {
    BOOST_TEST(outward.outward[i] == fixture.row(5U + i));
  }
  const auto innerTracks = fittedTrack(fixture, 3U, 7U, pool);
  const auto inward = followTrack(fixture, innerTracks, used, parameters);
  BOOST_TEST(inward.outwardCount == 0U);
  BOOST_REQUIRE(inward.inwardCount == 3U);
  for (std::size_t i = 0; i < 3U; ++i) {
    BOOST_TEST(inward.inward[i] == fixture.row(2U - i));
  }
  const auto middleTracks = fittedTrack(fixture, 2U, 5U, pool);
  const auto both = followTrack(fixture, middleTracks, used, parameters);
  BOOST_REQUIRE(both.outwardCount == 2U);
  BOOST_REQUIRE(both.inwardCount == 2U);
  BOOST_TEST(both.outward[0] == fixture.row(6U));
  BOOST_TEST(both.outward[1] == fixture.row(7U));
  BOOST_TEST(both.inward[0] == fixture.row(1U));
  BOOST_TEST(both.inward[1] == fixture.row(0U));

  for (const int threads : {1, 2, 4}) {
    for (const std::uint32_t chunk : {1U, 2U, 7U}) {
      BOOST_TEST(sameProposal(followTrack(fixture, middleTracks, used, parameters, UniformMaterial{}, nullptr, true, threads, chunk), both));
    }
  }
  const auto whole = fittedTrack(fixture, 0U, kLayers - 1U, pool);
  const auto none = followTrack(fixture, whole, used, parameters);
  BOOST_TEST(none.outwardCount + none.inwardCount == 0U);
}

/// Inward stays off where the iteration shares innermost rows.
BOOST_AUTO_TEST_CASE(InwardStaysOffWithSharedInnermostRows)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(0.42, 0.025, 1.8);
  const Fixture fixture{hits};
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, 2U, 5U, pool);
  const auto proposal = followTrack(fixture, tracks, used, followParameters(), UniformMaterial{}, nullptr, false);
  BOOST_TEST(proposal.outwardCount == 2U);
  BOOST_TEST(proposal.inwardCount == 0U);
}

/// Inward from the beam line's update: the
/// same rows on a track from the origin.
BOOST_AUTO_TEST_CASE(InwardFromTheBeamLineTakesTheSameRows)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(1.1, -0.02, 2.2);
  const Fixture fixture{hits};
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, 3U, 7U, pool);
  const auto plain = followTrack(fixture, tracks, used, followParameters());
  const auto beam = followTrack<true>(fixture, tracks, used, followParameters(0, 1., 1, true));
  BOOST_REQUIRE(plain.inwardCount == 3U);
  BOOST_TEST(sameProposal(beam, plain));
}

/// BeamEnd: a track on
/// layers 3-7 from the origin walks inward to rows 2, 1 and 0 with the
/// variant as without. Its row on layer 0 moved along r-phi: where the walk
/// without the variant still takes it, the variant keeps it while the step's
/// prediction updated with the beam line holds it under the gate and refuses
/// it beyond, the walk's other rows as before; the refusal comes before the
/// step's own.
BOOST_AUTO_TEST_CASE(AnInwardLastRowOffTheBeamLineIsRefused)
{
  using namespace toybarrel8;
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto parameters = followParameters(0, 1., 1, true);
  const auto clean = sensorHelixHits(1.1, -0.02, 2.2);
  const auto follow = [&](const std::vector<Hit>& hits) {
    const Fixture fixture{hits};
    const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
    const auto tracks = fittedTrack(fixture, 3U, 7U, pool);
    return std::pair{followTrack<true>(fixture, tracks, used, parameters), followTrack<true, false, false, false, true>(fixture, tracks, used, parameters)};
  };
  const auto [plain, beamEnd] = follow(clean);
  BOOST_REQUIRE(plain.inwardCount == 3U);
  BOOST_TEST(sameProposal(beamEnd, plain));
  bool keptSeen = false;
  bool refusedSeen = false;
  bool keptAfterRefusal = false;
  for (int step = 1; step <= 40; ++step) {
    const double offset = 0.0005 * step;
    auto hits = clean;
    hits[0].point = moved(hits[0].point, 0., offset);
    hits[0].loc0 += static_cast<Scalar>(offset);
    const auto [without, with] = follow(hits);
    if (without.inwardCount != 3U) {
      continue;
    }
    if (with.inwardCount == 2U) {
      refusedSeen = true;
      BOOST_TEST(with.inward[0] == without.inward[0]);
      BOOST_TEST(with.inward[1] == without.inward[1]);
    } else {
      BOOST_TEST(sameProposal(with, without));
      keptSeen = true;
      keptAfterRefusal = keptAfterRefusal || refusedSeen;
    }
  }
  BOOST_TEST(keptSeen);
  BOOST_TEST(refusedSeen);
  BOOST_TEST(!keptAfterRefusal);
}

BOOST_AUTO_TEST_CASE(AFollowRejectsUsedNoiseAndDisplacedRows)
{
  using namespace toybarrel8;
  const auto clean = sensorHelixHits(1.1, -0.02, 2.2);
  BOOST_REQUIRE(clean.size() == kLayers);
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto parameters = followParameters();

  // Row 5 taken: a hole on layer 5 ends the follow.
  const Fixture fixture{clean};
  const auto tracks = fittedTrack(fixture, 0U, 4U, pool);
  std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  used[fixture.row(5U)] = kExclusive;
  BOOST_TEST(followTrack(fixture, tracks, used, parameters).outwardCount == 0U);

  auto withNoise = clean;
  Hit noise = clean[5];
  noise.point = moved(noise.point, 0.02, 0.02);
  noise.loc0 += 0.02F;
  noise.loc1 += 0.02F;
  withNoise.push_back(noise);
  const Fixture noiseFixture{withNoise};
  const auto noiseTracks = fittedTrack(noiseFixture, 0U, 4U, pool);
  const std::vector<std::uint8_t> noiseUsed(noiseFixture.frame->size(), kFree);
  const auto noiseProposal = followTrack(noiseFixture, noiseTracks, noiseUsed, parameters);
  BOOST_REQUIRE(noiseProposal.outwardCount == 3U);
  BOOST_TEST(noiseProposal.outward[0] == noiseFixture.row(5U));

  // Row 5 displaced: no row on layer 5 within the gate, a hole.
  auto displaced = clean;
  displaced[5].point = moved(displaced[5].point, 0.02, 0.02);
  displaced[5].loc0 += 0.02F;
  displaced[5].loc1 += 0.02F;
  const Fixture displacedFixture{displaced};
  const auto displacedTracks = fittedTrack(displacedFixture, 0U, 4U, pool);
  const std::vector<std::uint8_t> displacedUsed(displacedFixture.frame->size(), kFree);
  BOOST_TEST(followTrack(displacedFixture, displacedTracks, displacedUsed, parameters).outwardCount == 0U);
}

BOOST_AUTO_TEST_CASE(AFollowSkipsATrackWithoutFittedEndStates)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(1.1, -0.02, 2.2);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  auto tracks = fittedTrack(fixture, 2U, 5U, pool);
  tracks.view().status[0] = PropagationStatus::kNoIntersection;
  tracks.view().state[0].surface = kInvalidSensorID;
  tracks.view().outwardState[0].surface = kInvalidSensorID;
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto proposal = followTrack(fixture, tracks, used, followParameters());
  BOOST_TEST(proposal.outwardCount == 0U);
  BOOST_TEST(proposal.inwardCount == 0U);
}

/// A hole: a layer whose sensor the track crosses without a free row. Past
/// the budget the direction ends; within it, it goes on -- in the toy no
/// link leads past a layer, so it ends there all the same, its holes spent.
BOOST_AUTO_TEST_CASE(AHoleEndsTheDirectionPastItsBudget)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(1.1, -0.02, 2.2);
  const Fixture fixture{hits};
  std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  used[fixture.row(5U)] = kExclusive;
  used[fixture.row(1U)] = kExclusive;
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, 2U, 4U, pool);
  for (const int holes : {0, 1}) {
    const auto proposal = followTrack(fixture, tracks, used, followParameters(holes));
    BOOST_TEST(proposal.outwardCount == 0U);
    BOOST_TEST(proposal.inwardCount == 0U);
  }
}

/// The unconfirmed drop: row 6 taken, the follow takes row 5,
/// then crosses layer 6's sensor without a row, and drops row 5, which
/// nothing confirms; with every row free it reaches the last layer, beyond
/// which nothing is crossed, and keeps all three.
BOOST_AUTO_TEST_CASE(AnUnconfirmedLastRowIsDroppedAndOneWithNothingBeyondKept)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(0.42, 0.025, 1.8);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, 0U, 4U, pool);
  std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  used[fixture.row(6U)] = kExclusive;
  BOOST_TEST(followTrack(fixture, tracks, used, followParameters()).outwardCount == 0U);
  const std::vector<std::uint8_t> free(fixture.frame->size(), kFree);
  BOOST_TEST(followTrack(fixture, tracks, free, followParameters()).outwardCount == 3U);
}

/// The unconfirmed drop inward where the step that found nothing saw a row near: row 0 moved so that it lies within the step's
/// window but above the gate. The inward walk takes row 1, refuses the moved
/// row 0 and crosses layer 0's sensor without a row, and keeps row 1 all the
/// same; the moved row is not taken. Moved out of the window, nothing is
/// near and row 1 goes as before. Outward the same evidence does not keep a
/// row: a track on layers 0-4 with row 6 moved loses row 5 either way. The
/// offsets are found by the step's own outcomes.
BOOST_AUTO_TEST_CASE(ARowRefusedAboveTheGateKeepsTheRowBeforeInwardOnly)
{
  using namespace toybarrel8;
  struct Seen : follow::NoObserver {
    std::vector<std::tuple<LayerID, std::uint32_t, StepOutcome>>* rows = nullptr;
    void candidate(LayerID layer, std::uint32_t row, StepOutcome outcome, Scalar /*chi2*/) const { rows->emplace_back(layer, row, outcome); }
  };
  const auto hits = sensorHelixHits(0.42, 0.025, 1.8);
  const auto parameters = followParameters();
  const UniformMaterial material;
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto outcomeOf = [](const std::vector<std::tuple<LayerID, std::uint32_t, StepOutcome>>& seen, std::uint32_t row) {
    for (const auto& [layer, at, outcome] : seen) {
      if (at == row) {
        return std::optional{outcome};
      }
    }
    return std::optional<StepOutcome>{};
  };
  for (const bool inwardCase : {true, false}) {
    bool gatedSeen = false;
    bool outsideSeen = false;
    for (int step = 0; step < 38; ++step) {
      const double offset = 0.001 + (0.0005 * step);
      auto movedHits = hits;
      const std::size_t movedLayer = inwardCase ? 0U : 6U;
      movedHits[movedLayer].point = moved(movedHits[movedLayer].point, offset, offset);
      movedHits[movedLayer].loc0 += static_cast<Scalar>(offset);
      movedHits[movedLayer].loc1 += static_cast<Scalar>(offset);
      const Fixture fixture{movedHits};
      const auto tracks = inwardCase ? fittedTrack(fixture, 2U, 6U, pool) : fittedTrack(fixture, 0U, 4U, pool);
      std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
      if (inwardCase) {
        used[fixture.row(7U)] = kExclusive; // nothing outward
      }
      const auto widths = widthBounds(fixture);
      std::vector<ToyProposal> proposals(1U);
      std::vector<follow::KinkSteps> kinks(1U);
      const FollowStage<Det, UniformMaterial> stage{&parameters, fixture.frame.get(), &trackingGeometry(), &material, Propagator{}, tracks.view(), used.data(), widths.data(), toyCells().view(), proposals.data(), kinks.data(), true, 0U, nullptr};
      std::vector<std::tuple<LayerID, std::uint32_t, StepOutcome>> seen;
      const Seen observer{{}, &seen};
      const auto view = tracks.view();
      stage.follow(view.state[0], view.outwardState[0], view.clusters[0], view.count[0], proposals[0], kinks[0], observer);
      const auto outcome = outcomeOf(seen, fixture.row(movedLayer));
      BOOST_REQUIRE(outcome.has_value());
      if (*outcome == StepOutcome::kChi2) { // NOLINT(bugprone-unchecked-optional-access): BOOST_REQUIRE checks it
        gatedSeen = true;
        if (inwardCase) {
          BOOST_TEST(proposals[0].inwardCount == 1U);
          BOOST_TEST(proposals[0].inward[0] == fixture.row(1U));
        } else {
          BOOST_TEST(proposals[0].outwardCount == 0U);
        }
      } else if (*outcome == StepOutcome::kOutsideWindow) {
        outsideSeen = true;
        BOOST_TEST((inwardCase ? proposals[0].inwardCount : proposals[0].outwardCount) == 0U);
      }
    }
    BOOST_TEST(gatedSeen);
    BOOST_TEST(outsideSeen);
  }
}

/// EndDensity's layout-end refusal: a track
/// on layers 1-7 has room for one row, so its inward walk takes row 0 and
/// leaves through the capacity, at layer 0, which no link leads inward
/// from: the row goes where its update's sqrt|V| is above the key, stays
/// below it, and stays without the variant or with the key 0. Outward, at
/// the last layer, the key never acts; inward, on a walk that stops at
/// layer 0 for want of targets, it does.
BOOST_AUTO_TEST_CASE(AtTheLayoutsEndAWidePredictionLosesItsLastRow)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(0.42, 0.025, 1.8);
  const Fixture fixture{hits};
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, 1U, 7U, pool);
  const auto [chi2, sqrtDet] = firstStep(fixture, tracks.view().state[0], fixture.row(0U), PropagationDirection::kBackward);
  auto parameters = followParameters();
  BOOST_REQUIRE(followTrack(fixture, tracks, used, parameters).inwardCount == 1U);
  parameters.followEndMaxSqrtDet = static_cast<Scalar>(0.99 * sqrtDet);
  BOOST_TEST(followTrack(fixture, tracks, used, parameters).inwardCount == 1U);
  BOOST_TEST((followTrack<false, false, false, true>(fixture, tracks, used, parameters).inwardCount == 0U));
  parameters.followEndMaxSqrtDet = static_cast<Scalar>(1.01 * sqrtDet);
  BOOST_TEST((followTrack<false, false, false, true>(fixture, tracks, used, parameters).inwardCount == 1U));
  parameters.followEndMaxSqrtDet = 0.F;
  BOOST_TEST((followTrack<false, false, false, true>(fixture, tracks, used, parameters).inwardCount == 1U));

  parameters.followEndMaxSqrtDet = 1e-12F;
  const auto outward = followTrack<false, false, false, true>(fixture, fittedTrack(fixture, 0U, 4U, pool), used, parameters);
  BOOST_TEST(outward.outwardCount == 3U);
  std::vector<std::uint8_t> outerTaken(fixture.frame->size(), kFree);
  outerTaken[fixture.row(7U)] = kExclusive;
  const auto inward = followTrack<false, false, false, true>(fixture, fittedTrack(fixture, 2U, 6U, pool), outerTaken, parameters);
  BOOST_REQUIRE(inward.inwardCount == 1U);
  BOOST_TEST(inward.inward[0] == fixture.row(1U));
}

/// The layout-end refusal and a kink step's witness: a track turned at layer 2 takes row 1
/// by the kink allowance, standing on row 0; refusing row 0 at the layout's
/// end leaves the kink step without its witness, and `KeepWitnesses` takes
/// row 1 with it.
BOOST_AUTO_TEST_CASE(ARefusedLastRowTakesTheKinkItWitnessed)
{
  using namespace toybarrel8;
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto kinked = kinkedAt(sensorHelixHits(0.42, 0.025, 1.8), 2U, 0.012);
  const Fixture fixture{kinked};
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto tracks = fittedTrack(fixture, 2U, 7U, pool);
  auto parameters = followParameters(0, 10., 1);
  follow::KinkSteps steps;
  const auto kept = followTrack<false, false, false, true>(fixture, tracks, used, parameters, kKinkMaterial, &steps);
  BOOST_REQUIRE(kept.inwardCount == 2U);
  BOOST_REQUIRE(steps.inward == 1U);
  BOOST_REQUIRE(steps.inwardWitnessed == 1U);
  parameters.followEndMaxSqrtDet = 1e-12F;
  auto refused = followTrack<false, false, false, true>(fixture, tracks, used, parameters, kKinkMaterial, &steps);
  BOOST_REQUIRE(refused.inwardCount == 1U);
  BOOST_TEST(refused.inward[0] == fixture.row(1U));
  follow::KeepWitnesses<kMaxClusters>{.proposals = &refused, .steps = &steps}(0U);
  BOOST_TEST(refused.inwardCount == 0U);
}

/// EndDensity's density keep: row 6 taken, the outward walk takes row 5 and
/// crosses layer 6's sensor without a row, and the unconfirmed drop removes row 5 -- unless its
/// density g = exp(-chi2 / 2) / (2 pi sqrt|V|) is above the key: kept just
/// below g, dropped just above it, without the variant or with the key 0.
/// Inward the key never keeps: row 0 taken, the walk drops row 1 as before.
BOOST_AUTO_TEST_CASE(ADenseLastRowStaysOutwardWhereTheUnconfirmedDropWouldRemoveIt)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(0.42, 0.025, 1.8);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, 0U, 4U, pool);
  std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  used[fixture.row(6U)] = kExclusive;
  const auto [chi2, sqrtDet] = firstStep(fixture, tracks.view().outwardState[0], fixture.row(5U), PropagationDirection::kForward);
  const double density = std::exp(-0.5 * chi2) / (2. * std::numbers::pi * sqrtDet);
  auto parameters = followParameters();
  parameters.followKeepDensity = static_cast<Scalar>(0.99 * density);
  BOOST_TEST(followTrack(fixture, tracks, used, parameters).outwardCount == 0U);
  const auto kept = followTrack<false, false, false, true>(fixture, tracks, used, parameters);
  BOOST_REQUIRE(kept.outwardCount == 1U);
  BOOST_TEST(kept.outward[0] == fixture.row(5U));
  parameters.followKeepDensity = static_cast<Scalar>(1.01 * density);
  BOOST_TEST((followTrack<false, false, false, true>(fixture, tracks, used, parameters).outwardCount == 0U));
  parameters.followKeepDensity = 0.F;
  BOOST_TEST((followTrack<false, false, false, true>(fixture, tracks, used, parameters).outwardCount == 0U));

  parameters.followKeepDensity = 1e-30F;
  std::vector<std::uint8_t> innerTaken(fixture.frame->size(), kFree);
  innerTaken[fixture.row(0U)] = kExclusive;
  innerTaken[fixture.row(7U)] = kExclusive;
  BOOST_TEST((followTrack<false, false, false, true>(fixture, fittedTrack(fixture, 2U, 6U, pool), innerTaken, parameters).inwardCount == 0U));
}

/// The lookahead on a clean track: a noise row beside the true one on layer 5, the
/// next layer confirms the true one.
BOOST_AUTO_TEST_CASE(TheConfirmationTakesTheTrueRowsOnACleanTrack)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(1.1, -0.02, 2.2);
  auto withNoise = hits;
  Hit noise = hits[5];
  noise.point = moved(noise.point, 0.002, 0.002);
  noise.loc0 += 0.002F;
  noise.loc1 += 0.002F;
  withNoise.push_back(noise);
  const Fixture fixture{withNoise};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, 0U, 4U, pool);
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto proposal = followTrack(fixture, tracks, used, followParameters());
  BOOST_REQUIRE(proposal.outwardCount == 3U);
  for (std::size_t i = 0; i < 3U; ++i) {
    BOOST_TEST(proposal.outward[i] == fixture.row(5U + i));
  }
}

/// The lookahead where the greedy choice is wrong: the particle kinks at layer 5, its
/// rows 5-7 shifted in r phi by 12, 24 and 36 um; a noise row where its row
/// 5 would have been, nearer the prediction. The next layer confirms the
/// particle's row and refuses the noise row.
BOOST_AUTO_TEST_CASE(TheConfirmationTakesTheRowTheNextLayerConfirms)
{
  using namespace toybarrel8;
  auto hits = sensorHelixHits(1.1, -0.02, 2.2);
  const Hit noise = hits[5];
  hits = kinkedBeyond(std::move(hits), 4U, 0.006);
  hits.push_back(noise);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, 0U, 4U, pool);
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto proposal = followTrack(fixture, tracks, used, followParameters());
  BOOST_REQUIRE(proposal.outwardCount == 3U);
  for (std::size_t i = 0; i < 3U; ++i) {
    BOOST_TEST(proposal.outward[i] == fixture.row(5U + i));
  }
  // The greedy choice takes the noise
  // row, nearer the prediction, where the lookahead takes the confirmed one; the next
  // layer refuses it, so nothing stands on it and the unconfirmed drop removes it: the
  // direction ends with no row.
  const auto greedy = followTrack<false, false, true>(fixture, tracks, used, followParameters());
  BOOST_TEST(greedy.outwardCount == 0U);
}

/// The lookahead where nothing can confirm, on the last layer outward: two rows pass the
/// step, and it takes the smaller chi2, the true row, whichever of them comes
/// first in the frame.
BOOST_AUTO_TEST_CASE(UnconfirmedRowsRankByTheirChi2)
{
  using namespace toybarrel8;
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto hits = sensorHelixHits(0.42, 0.025, 1.8);
  bool decoyFirst = false;
  for (const bool decoyEarlier : {false, true}) {
    for (const double offset : {-0.0015, 0.0015}) {
      auto withDecoy = hits;
      Hit decoy = hits[7];
      decoy.point = moved(decoy.point, offset, offset);
      decoy.loc0 += static_cast<Scalar>(offset);
      decoy.loc1 += static_cast<Scalar>(offset);
      withDecoy.push_back(decoy);
      if (decoyEarlier) {
        std::swap(withDecoy[7], withDecoy[kLayers]);
      }
      const std::size_t trueCluster = decoyEarlier ? kLayers : 7U;
      const std::size_t decoyCluster = decoyEarlier ? 7U : kLayers;
      const Fixture fixture{withDecoy};
      decoyFirst = decoyFirst || fixture.row(decoyCluster) < fixture.row(trueCluster);
      const auto tracks = fittedTrack(fixture, 0U, 6U, pool);
      const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
      const auto proposal = followTrack(fixture, tracks, used, followParameters());
      BOOST_REQUIRE(proposal.outwardCount == 1U);
      BOOST_TEST(proposal.outward[0] == fixture.row(trueCluster));
    }
  }
  BOOST_TEST(decoyFirst);
}

/// The kink pass (the recovery pass): a track turned
/// at layer 3 by more than the nominal step admits is followed inward by the
/// kink allowance, its row on layer 2 standing on layer 1's, which confirms
/// it; without the allowance, or with the budget spent, no kink step;
/// without the rows inside, nothing confirms it and it takes nothing.
BOOST_AUTO_TEST_CASE(AKinkStepNeedsItsWitness)
{
  using namespace toybarrel8;
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto kinked = kinkedAt(sensorHelixHits(0.42, 0.025, 1.8), 3U, 0.012);
  const Fixture fixture{kinked};
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto tracks = fittedTrack(fixture, 3U, 7U, pool);
  follow::KinkSteps steps;
  BOOST_TEST(followTrack(fixture, tracks, used, followParameters(0, 1., 1), kKinkMaterial, &steps).inwardCount == 0U);
  const auto recovered = followTrack(fixture, tracks, used, followParameters(0, 10., 1), kKinkMaterial, &steps);
  BOOST_REQUIRE(recovered.inwardCount == 3U);
  for (std::size_t i = 0; i < 3U; ++i) {
    BOOST_TEST(recovered.inward[i] == fixture.row(2U - i));
  }
  BOOST_TEST(steps.inward == 1U);
  BOOST_TEST(steps.outward == 0U);
  BOOST_TEST(steps.inwardWitnessed == 1U);
  BOOST_TEST(followTrack(fixture, tracks, used, followParameters(0, 10., 0), kKinkMaterial, &steps).inwardCount == 0U);

  const std::vector<Hit> outer(kinked.begin() + 2, kinked.end());
  const Fixture unwitnessed{outer};
  const std::vector<std::uint8_t> outerUsed(unwitnessed.frame->size(), kFree);
  const auto outerTracks = fittedTrack(unwitnessed, 1U, 5U, pool);
  BOOST_TEST(followTrack(unwitnessed, outerTracks, outerUsed, followParameters(0, 10., 1), kKinkMaterial, &steps).inwardCount == 0U);
  BOOST_TEST(steps.inward == 0U);
}

/// The nominal pass first: on a clean track with the kink allowance, every
/// row is taken by the nominal step, and no kink is charged.
BOOST_AUTO_TEST_CASE(TheNominalPassComesFirst)
{
  using namespace toybarrel8;
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const Fixture fixture{sensorHelixHits(0.42, 0.025, 1.8)};
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto tracks = fittedTrack(fixture, 2U, 5U, pool);
  follow::KinkSteps steps;
  const auto proposal = followTrack(fixture, tracks, used, followParameters(0, 10., 1), kKinkMaterial, &steps);
  BOOST_TEST(proposal.outwardCount == 2U);
  BOOST_TEST(proposal.inwardCount == 2U);
  BOOST_TEST(steps.outward + steps.inward == 0U);
  BOOST_TEST(steps.outwardWitnessed + steps.inwardWitnessed == 0U);
}

/// One kink budget for both directions: a track turned beyond layer 5 and
/// inside layer 2 takes the outward kink step first; with one kink, the
/// inward one is not made; with two, both are.
BOOST_AUTO_TEST_CASE(TheDirectionsShareTheKinkBudget)
{
  using namespace toybarrel8;
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto hits = kinkedAt(kinkedBeyond(sensorHelixHits(0.42, 0.025, 1.8), 5U, kKinkOut), 2U, kKinkIn);
  const Fixture fixture{hits};
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto tracks = fittedTrack(fixture, 2U, 5U, pool);
  follow::KinkSteps steps;
  // Without the allowance neither side is reached.
  const auto none = followTrack(fixture, tracks, used, followParameters(0, 1., 2), kKinkMaterial, &steps);
  BOOST_TEST(none.outwardCount + none.inwardCount == 0U);
  const auto one = followTrack(fixture, tracks, used, followParameters(0, kBudgetKinkSigma, 1), kKinkMaterial, &steps);
  BOOST_TEST(steps.outward == 1U);
  BOOST_TEST(steps.inward == 0U);
  BOOST_TEST(one.outwardCount == 2U);
  BOOST_TEST(one.inwardCount == 0U);
  const auto two = followTrack(fixture, tracks, used, followParameters(0, kBudgetKinkSigma, 2), kKinkMaterial, &steps);
  BOOST_TEST(steps.outward == 1U);
  BOOST_TEST(steps.inward == 1U);
  BOOST_TEST(two.outwardCount == 2U);
  BOOST_TEST(two.inwardCount == 2U);
}

/// One time box for both directions: the road's rows at 30 ns, its outer
/// rows at 55 ns and its inner ones at 10 ns. Each fits the road's box
/// alone; once the outer rows have narrowed it, the inner ones no longer do.
BOOST_AUTO_TEST_CASE(TheDirectionsShareTheTimeBox)
{
  using namespace toybarrel8;
  auto hits = sensorHelixHits(0.42, 0.025, 1.8);
  for (std::size_t layer = 0; layer < kLayers; ++layer) {
    hits[layer].timeNS = 30.;
    if (layer < 3U) {
      hits[layer].timeNS = 10.;
    } else if (layer > 5U) {
      hits[layer].timeNS = 55.;
    }
  }
  const Fixture fixture{hits, 80.};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, 3U, 5U, pool);
  const std::vector<std::uint8_t> free(fixture.frame->size(), kFree);
  const auto both = followTrack(fixture, tracks, free, followParameters());
  BOOST_TEST(both.outwardCount == 2U);
  BOOST_TEST(both.inwardCount == 0U);
  std::vector<std::uint8_t> outerTaken(fixture.frame->size(), kFree);
  outerTaken[fixture.row(6U)] = kExclusive;
  outerTaken[fixture.row(7U)] = kExclusive;
  const auto inward = followTrack(fixture, tracks, outerTaken, followParameters());
  BOOST_TEST(inward.outwardCount == 0U);
  BOOST_TEST(inward.inwardCount == 3U);
}

/// A kink step keeps its witness through the bidding (the recovery pass): a
/// direction cut by the bidding on the witness drops
/// the kink step with it, and the rows from it on.
BOOST_AUTO_TEST_CASE(AKinkStepFallsWithItsWitness)
{
  BOOST_TEST(follow::keptWithWitnesses(3U, 0b001U) == 3U);
  BOOST_TEST(follow::keptWithWitnesses(1U, 0b001U) == 0U);
  BOOST_TEST(follow::keptWithWitnesses(2U, 0b011U) == 0U);
  BOOST_TEST(follow::keptWithWitnesses(3U, 0b011U) == 3U);
  BOOST_TEST(follow::keptWithWitnesses(2U, 0b100U) == 2U);
  BOOST_TEST(follow::keptWithWitnesses(3U, 0b100U) == 2U);
  BOOST_TEST(follow::keptWithWitnesses(0U, 0b111U) == 0U);

  ToyProposal proposal;
  proposal.inward[0] = 7U;
  proposal.inward[1] = 8U;
  proposal.inward[2] = 9U;
  proposal.inwardCount = 1U; // the bidding cut it on its second row
  proposal.outward[0] = 3U;
  proposal.outwardCount = 1U;
  const follow::KinkSteps steps{.outwardWitnessed = 0U, .inwardWitnessed = 0b001U, .outward = 0U, .inward = 1U};
  follow::KeepWitnesses<kMaxClusters>{.proposals = &proposal, .steps = &steps}(0U);
  BOOST_TEST(proposal.inwardCount == 0U);
  BOOST_TEST(proposal.outwardCount == 1U);
}

BOOST_AUTO_TEST_CASE(BiddingKeepsTheBetterTrackAndTruncatesTheOther)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(2.2, 0.015, 2.4);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  auto tracks = fittedTrack(fixture, 0U, 4U, pool);
  tracks.resize(2U);
  tracks.write(1U, tracks.read(0U));
  std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  Follow<Det> follow;
  follow.setMemoryPool(pool);
  BOOST_REQUIRE(followAll(fixture, followParameters(), tracks, used, follow) == StageStatus::kComplete);
  const auto proposals = follow.proposals();
  BOOST_REQUIRE(proposals.size() == 2U);
  BOOST_TEST(proposals[0].outwardCount == 3U);
  BOOST_TEST(proposals[1].outwardCount == 0U);
  const auto view = tracks.view();
  BOOST_TEST(view.count[0] == kLayers);
  BOOST_TEST(view.addedOutward[0] == 3U);
  BOOST_TEST(view.addedInward[0] == 0U);
  BOOST_TEST(view.count[1] == 5U);
  BOOST_TEST(view.addedOutward[1] == 0U);
  for (std::size_t layer = 5U; layer < kLayers; ++layer) {
    BOOST_TEST(used[fixture.row(layer)] == kExclusive);
  }
}

BOOST_AUTO_TEST_CASE(ARejectedRefitLeavesItsRowsFree)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(2.8, -0.01, 1.6);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  auto tracks = fittedTrack(fixture, 0U, 4U, pool);
  std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  auto parameters = followParameters();
  std::ranges::fill(parameters.roadChi2MaxByNdf, -1.F);
  Follow<Det> follow;
  follow.setMemoryPool(pool);
  BOOST_REQUIRE(followAll(fixture, parameters, tracks, used, follow) == StageStatus::kComplete);
  BOOST_REQUIRE(follow.proposals().size() == 1U);
  BOOST_TEST(follow.proposals()[0].outwardCount == 0U);
  BOOST_TEST(tracks.view().count[0] == 5U);
  for (std::size_t layer = 5U; layer < kLayers; ++layer) {
    BOOST_TEST(used[fixture.row(layer)] == kFree);
  }
}

/// The followed prefix: a refit the cut
/// refuses for the whole extension keeps its longest accepted prefix with
/// `followPrefix`, and its rows only are marked; without it, none.
BOOST_AUTO_TEST_CASE(ARefusedRefitKeepsItsLongestAcceptedPrefix)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(2.8, -0.01, 1.6);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  for (const bool prefix : {false, true}) {
    auto tracks = fittedTrack(fixture, 0U, 4U, pool);
    std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
    auto parameters = followParameters();
    parameters.followPrefix = prefix;
    // The track with all three outward rows refused, 2 * 8 - 5 degrees of
    // freedom; with two of them, and fewer, accepted.
    parameters.roadChi2MaxByNdf[(2U * kLayers) - static_cast<std::size_t>(kNBoundParams)] = -1.F;
    Follow<Det> follow;
    follow.setMemoryPool(pool);
    BOOST_REQUIRE(followAll(fixture, parameters, tracks, used, follow) == StageStatus::kComplete);
    BOOST_REQUIRE(follow.proposals().size() == 1U);
    BOOST_TEST_CONTEXT("followPrefix " << prefix)
    {
      const auto view = std::as_const(tracks).view();
      BOOST_TEST(follow.proposals()[0].outwardCount == (prefix ? 2U : 0U));
      BOOST_TEST(view.count[0] == (prefix ? 7U : 5U));
      BOOST_TEST(view.addedOutward[0] == (prefix ? 2U : 0U));
      for (std::size_t layer = 5U; layer < kLayers; ++layer) {
        BOOST_TEST(used[fixture.row(layer)] == (prefix && layer < 7U ? kExclusive : kFree));
      }
    }
  }
}

/// The unconfirmed drop's row is not marked: it stays free after the passes.
BOOST_AUTO_TEST_CASE(AnUnconfirmedRowStaysFreeAfterThePasses)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(0.42, 0.025, 1.8);
  const Fixture fixture{hits};
  std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  for (std::size_t layer = 0; layer <= 4U; ++layer) {
    used[fixture.row(layer)] = kExclusive;
  }
  used[fixture.row(6U)] = kExclusive;
  const auto pool = std::make_shared<BoundedMemoryResource>();
  auto tracks = fittedTrack(fixture, 0U, 4U, pool);
  Follow<Det> follow;
  follow.setMemoryPool(pool);
  BOOST_REQUIRE(followAll(fixture, followParameters(), tracks, used, follow) == StageStatus::kComplete);
  BOOST_TEST(std::as_const(tracks).view().count[0] == 5U);
  BOOST_TEST(used[fixture.row(5U)] == kFree);
}

/// The trim of an ambiguous last row: a second row within the gate on the outermost layer, 5 um
/// off the particle's, makes the follow's last row one of two passing. With
/// a `followTrimMargin` above the gate the row is dropped and left free, the
/// track refitted one row shorter; it stays with the margin 0 (off), with a
/// margin below the second row's chi2 less its own, without the second
/// row, and where the length after the follow would drop the track one row
/// shorter.
BOOST_AUTO_TEST_CASE(AnAmbiguousLastRowIsTrimmedWhereTheTrackKeepsItsCuts)
{
  using namespace toybarrel8;
  const auto clean = sensorHelixHits(0.42, 0.025, 1.8);
  auto ambiguous = clean;
  Hit other = clean[kLayers - 1U];
  other.point = moved(other.point, 0.0005, 0.0005);
  other.loc0 += 0.0005F;
  other.loc1 += 0.0005F;
  ambiguous.push_back(other);
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto run = [&](const std::vector<Hit>& hits, float margin, std::uint32_t minAfter) {
    const Fixture fixture{hits};
    auto tracks = fittedTrack(fixture, 0U, 4U, pool);
    std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
    auto parameters = followParameters();
    parameters.followTrimMargin = margin;
    parameters.minClustersAfterFollow = minAfter;
    Follow<Det> follow;
    follow.setMemoryPool(pool);
    BOOST_REQUIRE(followAll(fixture, parameters, tracks, used, follow) == StageStatus::kComplete);
    BOOST_REQUIRE(follow.proposals().size() == 1U);
    const auto view = std::as_const(tracks).view();
    return std::tuple{static_cast<std::uint32_t>(view.count[0]), static_cast<std::uint32_t>(follow.proposals()[0].outwardCount), static_cast<std::uint32_t>(view.addedOutward[0]), used[fixture.row(kLayers - 1U)]};
  };
  const auto [count, outward, added, last] = run(ambiguous, 100.F, 0U);
  BOOST_TEST(count == kLayers - 1U);
  BOOST_TEST(outward == 2U);
  BOOST_TEST(added == 2U);
  BOOST_TEST(last == kFree);
  using Case = std::tuple<std::vector<Hit>, float, std::uint32_t>;
  for (const auto& [hits, margin, minAfter] : {Case{ambiguous, 0.F, 0U}, Case{ambiguous, 1e-6F, 0U}, Case{clean, 100.F, 0U}, Case{ambiguous, 100.F, static_cast<std::uint32_t>(kLayers)}}) {
    BOOST_TEST_CONTEXT("margin " << margin << ", rows " << hits.size() << ", length after the follow " << minAfter)
    {
      const auto [kept, keptOutward, keptAdded, keptLast] = run(hits, margin, minAfter);
      BOOST_TEST(kept == kLayers);
      BOOST_TEST(keptOutward == 3U);
      BOOST_TEST(keptLast == kExclusive);
    }
  }
}

/// The passes add the rows in order -- inward ones first, then the road,
/// then the outward ones -- count them on the track, after what an earlier
/// follow added, and mark them.
BOOST_AUTO_TEST_CASE(TheFollowsPassesAddItsRowsToTheTrackAndMarkThem)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(0.42, 0.025, 1.8);
  const Fixture fixture{hits};
  std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto pool = std::make_shared<BoundedMemoryResource>();
  auto tracks = fittedTrack(fixture, 2U, 5U, pool);
  for (std::size_t layer = 2; layer <= 5U; ++layer) {
    used[fixture.row(layer)] = kExclusive;
  }
  tracks.view().addedOutward[0] = 1U; // an earlier follow's row, kept in the count
  Follow<Det> follow;
  follow.setMemoryPool(pool);
  BOOST_REQUIRE(followAll(fixture, followParameters(), tracks, used, follow) == StageStatus::kComplete);
  const auto view = std::as_const(tracks).view();
  BOOST_TEST(view.count[0] == kLayers);
  BOOST_TEST(view.addedOutward[0] == 3U);
  BOOST_TEST(view.addedInward[0] == 2U);
  for (std::size_t layer = 0; layer < kLayers; ++layer) {
    BOOST_TEST(view.clusters[0].index[layer] == fixture.row(layer));
    BOOST_TEST(used[fixture.row(layer)] == kExclusive);
  }
}

/// The road's end decided again: the particle kinks at
/// layer 4, its rows 4-7 shifted in r phi by 12, 24, 36 and 48 um; the road
/// ends on a noise row where its row 4 would have been. Refitted without it,
/// the follow takes the particle's row 4, which layer 5 confirms, in its
/// place, and goes on outward.
BOOST_AUTO_TEST_CASE(TheFollowReplacesARoadEndTheNextLayerRefuses)
{
  using namespace toybarrel8;
  auto hits = sensorHelixHits(1.1, -0.02, 2.2);
  const Hit noise = hits[4];
  hits = kinkedBeyond(std::move(hits), 3U, 0.006);
  hits.push_back(noise);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto noiseRow = fixture.row(kLayers);
  auto tracks = fittedRows(fixture, {fixture.row(0U), fixture.row(1U), fixture.row(2U), fixture.row(3U), noiseRow}, pool);
  std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  for (std::size_t layer = 0; layer <= 3U; ++layer) {
    used[fixture.row(layer)] = kExclusive;
  }
  used[noiseRow] = kExclusive;

  const auto plain = followTrack(fixture, tracks, used, followParameters());
  BOOST_TEST(plain.replacedEnd == ClusterList<kMaxClusters>::kNoCluster);
  const auto parameters = followParameters(0, 1., 1, false, true);
  const auto proposal = followTrack<false, true>(fixture, tracks, used, parameters);
  BOOST_TEST(proposal.replacedEnd == noiseRow);
  BOOST_REQUIRE(proposal.outwardCount == 4U);
  for (std::size_t i = 0; i < 4U; ++i) {
    BOOST_TEST(proposal.outward[i] == fixture.row(4U + i));
  }

  Follow<Det> follow;
  follow.setMemoryPool(pool);
  BOOST_REQUIRE(followAll(fixture, parameters, tracks, used, follow) == StageStatus::kComplete);
  const auto view = std::as_const(tracks).view();
  BOOST_REQUIRE(view.count[0] == kLayers);
  BOOST_TEST(view.addedOutward[0] == 4U);
  for (std::size_t layer = 0; layer < kLayers; ++layer) {
    BOOST_TEST(view.clusters[0].index[layer] == fixture.row(layer));
    BOOST_TEST(used[fixture.row(layer)] == kExclusive);
  }
  BOOST_TEST(used[noiseRow] == kFree);
}

/// On a clean road the follow takes its end again: the proposal is the
/// plain follow's, the end kept.
BOOST_AUTO_TEST_CASE(TheFollowKeepsACleanRoadEnd)
{
  using namespace toybarrel8;
  const auto hits = sensorHelixHits(0.42, 0.025, 1.8);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, 0U, 4U, pool);
  std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  for (std::size_t layer = 0; layer <= 4U; ++layer) {
    used[fixture.row(layer)] = kExclusive;
  }
  const auto plain = followTrack(fixture, tracks, used, followParameters());
  const auto again = followTrack<false, true>(fixture, tracks, used, followParameters(0, 1., 1, false, true));
  BOOST_TEST(again.replacedEnd == ClusterList<kMaxClusters>::kNoCluster);
  BOOST_TEST(sameProposal(again, plain));
}

/// Without the follow, the end alone is decided again: the particle's row
/// replaces the noise row, nothing beyond it is taken, nothing inward.
BOOST_AUTO_TEST_CASE(WithoutTheFollowOnlyTheEndIsDecidedAgain)
{
  using namespace toybarrel8;
  auto hits = sensorHelixHits(1.1, -0.02, 2.2);
  const Hit noise = hits[4];
  hits = kinkedBeyond(std::move(hits), 3U, 0.006);
  hits.push_back(noise);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto noiseRow = fixture.row(kLayers);
  const auto tracks = fittedRows(fixture, {fixture.row(1U), fixture.row(2U), fixture.row(3U), noiseRow}, pool);
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  auto parameters = followParameters(0, 1., 1, false, true);
  parameters.follow = false;
  const auto proposal = followTrack<false, true>(fixture, tracks, used, parameters);
  BOOST_TEST(proposal.replacedEnd == noiseRow);
  BOOST_REQUIRE(proposal.outwardCount == 1U);
  BOOST_TEST(proposal.outward[0] == fixture.row(4U));
  BOOST_TEST(proposal.inwardCount == 0U);
}

BOOST_AUTO_TEST_CASE(EverySensorIsListedAtItsCentre)
{
  const auto& toy = toybarrel8::trackingGeometry();
  const SensorCells<toybarrel8::Det> toyCells(toy);
  for (const auto& surface : toy.surfaces()) {
    BOOST_TEST(visitedAtItsCentre(toyCells, surface));
  }
  using Mixed = detectors::ToyBarrel4Disc8;
  const auto detector = toyGeometry<Mixed>();
  const TrackingGeometry<Mixed> mixed{detector};
  const SensorCells<Mixed> mixedCells(mixed);
  std::size_t missed = 0U;
  for (const auto& surface : mixed.surfaces()) {
    missed += visitedAtItsCentre(mixedCells, surface) ? 0U : 1U;
  }
  BOOST_TEST(missed == 0U);
}

BOOST_AUTO_TEST_CASE(AQueryAtAFootprintsOuterEdgeVisitsItsSensor)
{
  // Each layer's outermost footprints in the second coordinate: a query at
  // their edge, with no reach, visits them; one beyond, none.
  const auto& toy = toybarrel8::trackingGeometry();
  const SensorCells<toybarrel8::Det> cells(toy);
  for (const auto& surface : toy.surfaces()) {
    const auto footprint = SensorCells<toybarrel8::Det>::footprint(surface, true);
    const auto phi = static_cast<Scalar>(footprint.phiCentre < 0. ? footprint.phiCentre + (2 * std::numbers::pi_v<double>) : footprint.phiCentre);
    for (const double edge : {footprint.secondLow, footprint.secondHigh}) {
      BOOST_TEST(cells.view().anySensor(surface.layer, phi, 0.F, static_cast<Scalar>(edge), 0.F, [&](std::uint32_t id) { return id == surface.id; }));
    }
  }
  const auto& layer = cells.layers()[0];
  const Scalar beyond = layer.secondLow + (static_cast<Scalar>(layer.nSecond) * layer.secondBin) + 1.F;
  BOOST_TEST(!cells.view().anySensor(0U, 0.F, std::numbers::pi_v<Scalar>, beyond, 0.F, [](std::uint32_t /*id*/) { return true; }));
}

} // namespace o2::gtrk

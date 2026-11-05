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

/// \file testTrackFitter.cxx
/// \brief Seeding, the whole-chain seed and bidirectional fitting.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Track Fitter
#include <boost/test/unit_test.hpp>
#include <numbers>

#include "GridTestSupport.h"
#include "Oracles.h"
#include "GenericTracking/Base/Geometry.h"
#include "ToyBarrel4Disc8.h"
#include "ToyBarrel3.h"
#include "GenericTracking/Tracking/BeamConstraint.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/TrackFitter.h"

#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/TrackUtils.h"
#include "UniformMaterial.h"

#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

namespace o2::gtrk
{
namespace
{

struct PullMoments {
  double mean = 0.;
  double width = 0.;
};

PullMoments moments(std::span<const double> values)
{
  const auto n = static_cast<double>(values.size());
  double mean = 0.;
  for (const double value : values) {
    mean += value / n;
  }
  double variance = 0.;
  for (const double value : values) {
    variance += (value - mean) * (value - mean);
  }
  return {.mean = mean, .width = std::sqrt(variance / (n - 1.))};
}

using Toy = detectors::ToyBarrel3;

void collectPulls(int trial, const BoundTrackState& state, const BoundTrackState& truth, std::array<std::vector<double>, kNBoundParams>& pull)
{
  for (int i = 0; i < kNBoundParams; ++i) {
    double diff = static_cast<double>(state.par[i]) - static_cast<double>(truth.par[i]);
    if (i == kPhi) {
      math_utils::detail::bringToPMPi<double>(diff);
    }
    const auto variance = static_cast<double>(state.variance(i));
    BOOST_TEST_CONTEXT("trial " << trial << " parameter " << i << " variance " << variance) { BOOST_REQUIRE(variance > 0.); }
    pull[static_cast<std::size_t>(i)].push_back(diff / std::sqrt(variance));
  }
}

Measurement measurementFromState(const BoundTrackState& state, std::uint32_t index)
{
  return Measurement{.loc0 = state.loc0(),
                     .loc1 = state.loc1(),
                     .timePS = test::toyTimeBox(test::kToyTrackTimeNS).start,
                     .widthPS = test::toyTimeBox(test::kToyTrackTimeNS).width,
                     .varLoc0 = 2.5e-6F,
                     .covLoc01 = 0.F,
                     .varLoc1 = 2.5e-6F,
                     .clusterIndex = index,
                     .surface = state.surface};
}

struct SyntheticChain {
  Geometry<Toy> geometry = toyGeometry<Toy>();
  TrackingGeometry<Toy> tracking{geometry};
  Propagator propagator;
  std::array<BoundTrackState, 3> truth;
  std::array<Measurement, 3> measurements;

  SyntheticChain()
  {
    truth[0].par[kLoc0] = 0.F;
    truth[0].par[kLoc1] = -0.2F;
    truth[0].par[kPhi] = 0.F;
    truth[0].par[kTheta] = 1.5F;
    truth[0].par[kQOverP] = 0.5F;
    truth[0].surface = firstSensorOfLayer<Toy>(0);

    for (std::size_t layer = 1; layer < truth.size(); ++layer) {
      truth[layer] = truth[layer - 1];
      const auto& from = tracking.surface(truth[layer - 1].surface);
      const auto toID = firstSensorOfLayer<Toy>(static_cast<LayerID>(layer));
      const auto& to = tracking.surface(toID);
      const auto result = propagate(propagator, truth[layer], from, to, PropagationDirection::kForward);
      BOOST_REQUIRE(result.status == PropagationStatus::kOk);
    }
    for (std::size_t i = 0; i < measurements.size(); ++i) {
      measurements[i] = measurementFromState(truth[i], static_cast<std::uint32_t>(100 + i));
    }
  }
};

double transverseRadius2At(const Propagator& prop, const BoundTrackState& state, const TrackingSurface& surface, double path)
{
  Scalar global[3]{};
  surfaceToGlobal(surface, state.loc0(), state.loc1(), global);
  const double sinTheta = std::sin(static_cast<double>(state.theta()));
  const double transverse = path * sinTheta;
  const double curvature = prop.curvature(static_cast<double>(state.qOverP()), sinTheta);
  const double turn = curvature * transverse;
  const auto phi = static_cast<double>(state.phi());
  const double dx = transverse * (std::cos(phi) * sinc(turn) - std::sin(phi) * versine(turn));
  const double dy = transverse * (std::sin(phi) * sinc(turn) + std::cos(phi) * versine(turn));
  const double x = static_cast<double>(global[0]) + dx;
  const double y = static_cast<double>(global[1]) + dy;
  return x * x + y * y;
}

} // namespace

BOOST_AUTO_TEST_CASE(ThreePointSeedRecoversTheGeneratingHelix)
{
  const SyntheticChain chain;
  const auto seed = makeThreePointSeed<Scalar>(chain.propagator, chain.measurements[0], chain.measurements[1], chain.measurements[2],
                                               chain.tracking.surface(chain.measurements[0].surface), chain.tracking.surface(chain.measurements[1].surface),
                                               chain.tracking.surface(chain.measurements[2].surface));
  const auto exact = makeThreePointSeed<ExactScalar>(
    chain.propagator, chain.measurements[0], chain.measurements[1], chain.measurements[2], chain.tracking.surface(chain.measurements[0].surface),
    chain.tracking.surface(chain.measurements[1].surface), chain.tracking.surface(chain.measurements[2].surface));

  BOOST_TEST(std::abs(seed.loc0() - chain.truth[0].loc0()) < 1e-6F);
  BOOST_TEST(std::abs(seed.loc1() - chain.truth[0].loc1()) < 1e-6F);
  BOOST_TEST(std::abs(seed.phi() - chain.truth[0].phi()) < 1e-5F);
  BOOST_TEST(std::abs(seed.theta() - chain.truth[0].theta()) < 1e-5F);
  BOOST_TEST(std::abs(seed.qOverP() - chain.truth[0].qOverP()) < 1e-3F);
  BOOST_TEST(seed.surface == chain.truth[0].surface);

  BOOST_TEST(std::abs(seed.loc0() - exact.loc0()) < 1e-6F);
  BOOST_TEST(std::abs(seed.loc1() - exact.loc1()) < 1e-6F);
  BOOST_TEST(std::abs(seed.phi() - exact.phi()) < 2e-6F);
  BOOST_TEST(std::abs(seed.theta() - exact.theta()) < 2e-6F);
  BOOST_TEST(std::abs(seed.qOverP() - exact.qOverP()) < 2e-4F);
}

/// At fixed non-zero chord the asin argument saturates as |curvature| tends
/// to infinity. Both arithmetic types must reach the finite half-turn limit;
/// this also checks that the clamp handles infinity before asin sees it.
BOOST_AUTO_TEST_CASE(TransverseArcHasFiniteInfiniteCurvatureLimit)
{
  constexpr Scalar kChord = 3.F;
  const auto scalar = detail::transverseArc<Scalar>(std::numeric_limits<Scalar>::infinity(), kChord, 0.F);
  const auto exact = detail::transverseArc<ExactScalar>(std::numeric_limits<ExactScalar>::infinity(), static_cast<ExactScalar>(kChord), 0.);
  const auto negative = detail::transverseArc<Scalar>(-std::numeric_limits<Scalar>::infinity(), kChord, 0.F);
  const auto expected = 0.5 * std::acos(-1.) * static_cast<double>(kChord);

  BOOST_TEST(std::isfinite(scalar));
  BOOST_TEST(std::isfinite(exact));
  BOOST_TEST(std::isfinite(negative));
  BOOST_TEST(std::abs(static_cast<double>(scalar) - expected) < 1e-6);
  BOOST_TEST(std::abs(exact - expected) < 1e-12);
  BOOST_TEST(std::abs(static_cast<double>(negative) - expected) < 1e-6);
  BOOST_TEST(std::abs(static_cast<double>(scalar) - exact) < 1e-6);
}

/// Shrinking a similar three-point configuration sends its curvature to
/// infinity. The calculation must stay finite at every finite scale, and the
/// Scalar result must follow the ExactScalar oracle rather than losing the
/// triangle area or overflowing an intermediate chord product.
BOOST_AUTO_TEST_CASE(ThreePointSeedRemainsFiniteAsCurvatureDiverges)
{
  const Propagator prop;
  Scalar previousMagnitude = 0.F;
  for (const Scalar scale : {1e-2F, 1e-5F, 1e-8F}) {
    const TrackingSurface firstSurface;
    TrackingSurface secondSurface;
    TrackingSurface thirdSurface;
    secondSurface.origin[0] = scale;
    thirdSurface.origin[0] = scale;
    thirdSurface.origin[1] = scale;

    Measurement first;
    const Measurement second;
    const Measurement third;
    first.widthPS = 1'000U;

    const auto scalar = makeThreePointSeed<Scalar>(prop, first, second, third, firstSurface, secondSurface, thirdSurface);
    const auto exact = makeThreePointSeed<ExactScalar>(prop, first, second, third, firstSurface, secondSurface, thirdSurface);
    const Scalar magnitude = std::abs(scalar.qOverP());
    BOOST_TEST_CONTEXT("scale " << scale)
    {
      BOOST_TEST(std::isfinite(scalar.qOverP()));
      BOOST_TEST(std::isfinite(exact.qOverP()));
      BOOST_TEST(magnitude > previousMagnitude);
      BOOST_TEST(std::abs((scalar.qOverP() - exact.qOverP()) / exact.qOverP()) < 2e-5F);
    }
    previousMagnitude = magnitude;
  }
}

BOOST_AUTO_TEST_CASE(IndependentPassesFitTheSameNoiselessChain)
{
  const SyntheticChain chain;
  const UniformMaterial material;
  const auto fit = fitMeasurements(chain.propagator, chain.tracking, material, chain.measurements.data(), static_cast<int>(chain.measurements.size()));

  BOOST_TEST(static_cast<int>(fit.status) == static_cast<int>(PropagationStatus::kOk));
  BOOST_TEST(fit.ndf == 1);
  BOOST_TEST(fit.outwardState.surface == chain.truth[2].surface);
  BOOST_TEST(fit.state.surface == chain.truth[0].surface);
  BOOST_TEST(std::abs(fit.state.phi() - chain.truth[0].phi()) < 1e-5F);
  BOOST_TEST(std::abs(fit.state.theta() - chain.truth[0].theta()) < 1e-5F);
  BOOST_TEST(std::abs(fit.state.qOverP() - chain.truth[0].qOverP()) < 1e-3F);
  BOOST_TEST(std::abs(fit.chi2) < 1e-4F);
  BOOST_TEST(std::abs(fit.outwardChi2) < 1e-4F);
}

/// A primary from a vertex in the diamond's tail crosses the far end of barrel
/// layer 0 at |eta| ~ 3.7, then the forward discs. Its
/// layer-0 incidence, |n . d| ~ 0.048, is an ordinary track's, and both passes
/// fit it. A threshold of 0.05 would pass the outward pass, which starts on
/// layer 0 and so never checks it, and fail the inward pass, arriving there
/// last.
BOOST_AUTO_TEST_CASE(LayerZeroAtShallowIncidenceFitsBothWays)
{
  using Mixed = detectors::ToyBarrel4Disc8;
  const auto detector = toyGeometry<Mixed>();
  const TrackingGeometry<Mixed> tracking{detector};
  const Propagator propagator;
  const UniformMaterial material;

  SensorID start = firstSensorOfLayer<Mixed>(0);
  for (SensorID sensor = start; sensor < firstSensorOfLayer<Mixed>(1); ++sensor) {
    if (tracking.surface(sensor).origin[2] > tracking.surface(start).origin[2]) {
      start = sensor;
    }
  }
  const auto& first = tracking.surface(start);
  std::array<BoundTrackState, 5> truth;
  truth[0].par[kPhi] = std::atan2(first.origin[1], first.origin[0]);
  truth[0].par[kTheta] = 0.048F;
  truth[0].par[kQOverP] = 0.05F;
  truth[0].surface = start;
  const double cosIncidence = std::sin(0.048) * (std::cos(static_cast<double>(truth[0].phi())) * static_cast<double>(first.normal[0]) +
                                                 std::sin(static_cast<double>(truth[0].phi())) * static_cast<double>(first.normal[1]));
  BOOST_TEST_REQUIRE(std::abs(cosIncidence) < 0.05);
  BOOST_TEST_REQUIRE(std::abs(cosIncidence) > static_cast<double>(propagator.minAbsCosIncidence));

  constexpr std::array<LayerID, 4> kDiscs{8, 9, 10, 11};
  for (std::size_t i = 0; i < kDiscs.size(); ++i) {
    const auto& from = tracking.surface(truth[i].surface);
    for (SensorID sensor = firstSensorOfLayer<Mixed>(kDiscs[i]); sensor < firstSensorOfLayer<Mixed>(static_cast<LayerID>(kDiscs[i] + 1U)); ++sensor) {
      BoundTrackState probe = truth[i];
      if (propagate(propagator, probe, from, tracking.surface(sensor), PropagationDirection::kForward).status == PropagationStatus::kOk) {
        truth[i + 1U] = probe;
        break;
      }
    }
    BOOST_TEST_REQUIRE(truth[i + 1U].surface != truth[i].surface);
  }
  std::array<Measurement, 5> measurements;
  for (std::size_t i = 0; i < measurements.size(); ++i) {
    measurements[i] = measurementFromState(truth[i], static_cast<std::uint32_t>(i));
  }

  const auto fit = fitMeasurements(propagator, tracking, material, measurements.data(), static_cast<int>(measurements.size()));
  BOOST_TEST(static_cast<int>(fit.status) == static_cast<int>(PropagationStatus::kOk));
  BOOST_TEST(fit.state.surface == start);
  BOOST_TEST(std::abs(fit.state.theta() - truth[0].theta()) < 1e-5F);

  const Propagator former{.minAbsCosIncidence = 0.05F};
  const auto failed = fitMeasurements(former, tracking, material, measurements.data(), static_cast<int>(measurements.size()));
  BOOST_TEST(static_cast<int>(failed.status) == static_cast<int>(PropagationStatus::kGrazing));
  BOOST_TEST(static_cast<int>(failed.failedPass) == static_cast<int>(FitPass::kInward));
  BOOST_TEST(failed.failedMeasurement == 0);
}

BOOST_AUTO_TEST_CASE(SmearedChainsGiveUnbiasedParametersAndUnitPulls)
{
  // A fit can be self-consistent and still wrong in
  // exactly one way that matters here: if the inward pass inherited the
  // outward pass's covariance it would consume every measurement twice, the
  // estimate would be too confident by root two, and the pull widths would
  // come out near 0.707 while every other test in this file still passed.
  // That number is asserted against explicitly.
  const SyntheticChain chain;
  const UniformMaterial material;

  constexpr int kTrials = 3000;
  const Scalar sigmaLoc = 1.5811388e-3F; // sqrt(2.5e-6) cm

  std::mt19937 rng(20260919);
  std::normal_distribution<double> gauss(0., 1.);

  std::array<std::vector<double>, kNBoundParams> pull;

  for (int trial = 0; trial < kTrials; ++trial) {
    auto smeared = chain.measurements;
    for (auto& m : smeared) {
      m.loc0 += static_cast<Scalar>(gauss(rng) * static_cast<double>(sigmaLoc));
      m.loc1 += static_cast<Scalar>(gauss(rng) * static_cast<double>(sigmaLoc));
    }

    const auto fit = fitMeasurements(chain.propagator, chain.tracking, material, smeared.data(), static_cast<int>(smeared.size()));
    BOOST_REQUIRE(fit.status == PropagationStatus::kOk);

    collectPulls(trial, fit.state, chain.truth[0], pull);
  }

  const char* names[kNBoundParams] = {"loc0", "loc1", "phi", "theta", "qOverP"};
  for (int i = 0; i < kNBoundParams; ++i) {
    const auto& p = pull[static_cast<std::size_t>(i)];
    const auto n = static_cast<double>(p.size());

    const auto [meanPull, width] = moments(p);

    BOOST_TEST_CONTEXT("parameter " << names[i] << " pull mean " << meanPull << " width " << width)
    {
      // The mean of kTrials unit-width pulls has standard error 1/sqrt(N),
      // which is 0.018 here; four of those is a comfortable band that still
      // catches any real bias.
      BOOST_TEST(std::abs(meanPull) < 4. / std::sqrt(n));
      // The width itself is known to 1/sqrt(2N) = 1.3 %, so 10 % is generous
      // -- and 0.707 is 29 % low, far outside it.
      BOOST_TEST(std::abs(width - 1.) < 0.10);
    }
  }
}

/// A chain over the first \a N barrel layers of ToyBarrel4Disc8: tilted,
/// overlapping sensors and 2N spatial constraints against five parameters, so
/// the fit has genuine redundancy. ToyBarrel3 cannot show this -- three points
/// determine a helix exactly, and the accumulated chi2 is zero however badly a
/// point is moved. With N = 4 the chain runs from 3 to 16 cm.
template <std::size_t N>
struct BarrelChain {
  using Mixed = detectors::ToyBarrel4Disc8;
  Geometry<Mixed> geometry = toyGeometry<Mixed>();
  TrackingGeometry<Mixed> tracking{geometry};
  Propagator propagator;
  std::array<BoundTrackState, N> truth;
  std::array<Measurement, N> measurements;
  bool built = false;

  /// The sensor of \a layer that this track actually crosses. A barrel layer
  /// overlaps in azimuth by about a fifth, so "the sensor at this phi" is not
  /// a well-defined question -- the first one whose bounds contain the
  /// intersection will do.
  [[nodiscard]] SensorID crossedSensor(const BoundTrackState& state, const TrackingSurface& from, LayerID layer) const
  {
    const SensorID begin = firstSensorOfLayer<Mixed>(layer);
    const SensorID end = firstSensorOfLayer<Mixed>(static_cast<LayerID>(layer + 1));
    for (SensorID id = begin; id < end; ++id) {
      BoundTrackState probe = state;
      const auto result = propagate(propagator, probe, from, tracking.surface(id), PropagationDirection::kForward);
      if (result.status == PropagationStatus::kOk) {
        return id;
      }
    }
    return kInvalidSensorID;
  }

  explicit BarrelChain(Scalar qOverP = 0.5F) { generate(qOverP); }

  /// (Re)build the chain. With \a x2X0PerCm the track scatters as the fit
  /// models it (`StepMaterial`): each step's material in two thin scatterers,
  /// the half behind the sensor just left kicked on leaving it, the half in
  /// front of the sensor reached kicked on arrival, both with the Highland
  /// theta0 of half the step's chord. truth[k] is the state at the mid-plane
  /// of sensor k: after the kick in front of it, before the one behind it,
  /// which is what both passes estimate there. With \a rho [g/cm^3] the
  /// track loses the fit's mean energy over each step's chord on reaching the
  /// next sensor, without straggling.
  void generate(Scalar qOverP, Scalar x2X0PerCm = 0.F, std::mt19937* rng = nullptr, Scalar rho = 0.F)
  {
    built = false;
    truth[0] = BoundTrackState{};
    truth[0].par[kPhi] = 0.12F; // along the stave tilt, so the track leaves the sensor
    truth[0].par[kTheta] = 1.4F;
    truth[0].par[kQOverP] = qOverP;
    truth[0].surface = firstSensorOfLayer<Mixed>(0);

    std::normal_distribution<double> gauss(0., 1.);
    const auto kick = [&](BoundTrackState& state, Scalar theta0) {
      const auto sinTheta = std::sin(state.theta());
      state.par[kTheta] += static_cast<Scalar>(gauss(*rng)) * theta0;
      state.par[kPhi] += static_cast<Scalar>(gauss(*rng)) * theta0 / sinTheta;
    };
    for (std::size_t layer = 1; layer < truth.size(); ++layer) {
      const auto& from = tracking.surface(truth[layer - 1].surface);
      BoundTrackState outgoing = truth[layer - 1];
      Scalar theta0Half = 0.F;
      if (x2X0PerCm > 0.F && rng != nullptr) {
        const auto unkicked = crossedSensor(outgoing, from, static_cast<LayerID>(layer));
        if (unkicked == kInvalidSensorID) {
          return;
        }
        BoundTrackState probe = outgoing;
        (void)propagate(propagator, probe, from, tracking.surface(unkicked), PropagationDirection::kForward);
        Scalar a[3]{};
        Scalar b[3]{};
        surfaceToGlobal(from, outgoing.loc0(), outgoing.loc1(), a);
        surfaceToGlobal(tracking.surface(unkicked), probe.loc0(), probe.loc1(), b);
        const Scalar chord = std::hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
        theta0Half = highlandTheta0(0.5F * x2X0PerCm * chord, betaMomentum(1.F / std::abs(outgoing.qOverP()), propagator.massGeV));
        kick(outgoing, theta0Half); // behind the sensor just left
      }
      const auto id = crossedSensor(outgoing, from, static_cast<LayerID>(layer));
      if (id == kInvalidSensorID) {
        return;
      }
      truth[layer] = outgoing;
      const auto result = propagate(propagator, truth[layer], from, tracking.surface(id), PropagationDirection::kForward);
      if (result.status != PropagationStatus::kOk) {
        return;
      }
      if (rho > 0.F) {
        Scalar a[3]{};
        Scalar b[3]{};
        surfaceToGlobal(from, truth[layer - 1].loc0(), truth[layer - 1].loc1(), a);
        surfaceToGlobal(tracking.surface(id), truth[layer].loc0(), truth[layer].loc1(), b);
        Scalar variance = 0.F;
        if (!detail::applyEnergyLoss(propagator, rho * std::hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]), PropagationDirection::kForward, truth[layer],
                                     variance)) {
          return;
        }
      }
      if (theta0Half > 0.F) {
        kick(truth[layer], theta0Half); // in front of the sensor reached: the direction at its mid-plane
      }
    }
    for (std::size_t i = 0; i < measurements.size(); ++i) {
      measurements[i] = measurementFromState(truth[i], static_cast<std::uint32_t>(200 + i));
    }
    built = true;
  }
};

/// The chi2 is the positions': two per measurement against the five spatial
/// parameters, and a measurement's time, its box, does not enter it.
BOOST_AUTO_TEST_CASE(TheChiSquaredIsSpatialAlone)
{
  const BarrelChain<4> chain;
  BOOST_REQUIRE(chain.built);
  const UniformMaterial material;

  const auto clean = fitMeasurements(chain.propagator, chain.tracking, material, chain.measurements.data(), static_cast<int>(chain.measurements.size()));
  BOOST_REQUIRE(clean.status == PropagationStatus::kOk);
  BOOST_TEST(clean.ndf == 3);
  BOOST_TEST(std::abs(clean.chi2) < 1e-2F);

  auto spatialOutlier = chain.measurements;
  spatialOutlier[1].loc0 += 0.02F; // about thirteen sigma
  const auto spatial = fitMeasurements(chain.propagator, chain.tracking, material, spatialOutlier.data(), static_cast<int>(spatialOutlier.size()));
  BOOST_REQUIRE(spatial.status == PropagationStatus::kOk);
  BOOST_TEST(spatial.chi2 > 10.F);

  auto timeOutlier = chain.measurements;
  timeOutlier[2].timePS += 30'000U;
  const auto timed = fitMeasurements(chain.propagator, chain.tracking, material, timeOutlier.data(), static_cast<int>(timeOutlier.size()));
  BOOST_REQUIRE(timed.status == PropagationStatus::kOk);
  BOOST_TEST(timed.chi2 == clean.chi2);
}

/// The whole-chain seed recovers the helix that made the chain, over the
/// whole barrel, at both charges and from 0.4 to 20 GeV/c.
/// The positions it fits are float, which is what bounds the tolerances.
BOOST_AUTO_TEST_CASE(ChainSeedRecoversTheGeneratingHelix)
{
  for (const Scalar qOverP : {2.5F, -2.5F, 0.5F, -0.1F, 0.05F}) {
    const BarrelChain<4> chain(qOverP);
    BOOST_TEST_CONTEXT("q/p " << qOverP)
    {
      BOOST_REQUIRE(chain.built);
      const auto scalarResult = makeChainSeed<Scalar>(chain.propagator, chain.tracking, chain.measurements.data(), static_cast<int>(chain.measurements.size()));
      const auto exactResult =
        makeChainSeed<ExactScalar>(chain.propagator, chain.tracking, chain.measurements.data(), static_cast<int>(chain.measurements.size()));
      const auto& reference = scalarResult.state;
      const auto& exact = exactResult.state;
      const auto& truth = chain.truth[0];
      BOOST_TEST(reference.surface == truth.surface);
      BOOST_TEST(std::abs(reference.loc0() - truth.loc0()) < 1e-4F);
      BOOST_TEST(std::abs(reference.loc1() - truth.loc1()) < 1e-4F);
      BOOST_TEST(std::abs(reference.phi() - truth.phi()) < 1e-5F);
      BOOST_TEST(std::abs(reference.theta() - truth.theta()) < 1e-5F);
      BOOST_TEST(std::abs(reference.qOverP() / truth.qOverP() - 1.F) < 1e-3F);

      BOOST_TEST(std::abs(reference.loc0() - exact.loc0()) < 5e-5F);
      BOOST_TEST(std::abs(reference.loc1() - exact.loc1()) < 5e-5F);
      BOOST_TEST(std::abs(reference.phi() - exact.phi()) < 5e-6F);
      BOOST_TEST(std::abs(reference.theta() - exact.theta()) < 5e-6F);
      BOOST_TEST(std::abs(reference.qOverP() - exact.qOverP()) < 5e-4F);
      BOOST_TEST(std::abs(scalarResult.qOverPSigma - exactResult.qOverPSigma) < 1e-3F);
    }
  }
}

/// Karimaki's fit stays regular as the curvature vanishes; a straight track is
/// not a special case, and must not come out as a NaN or a spurious bend.
BOOST_AUTO_TEST_CASE(ChainSeedOfAStraightTrackHasNoCurvature)
{
  const BarrelChain<4> chain(0.F);
  BOOST_REQUIRE(chain.built);
  const auto reference = makeChainSeed<Scalar>(chain.propagator, chain.tracking, chain.measurements.data(), static_cast<int>(chain.measurements.size())).state;
  const auto exact = makeChainSeed<ExactScalar>(chain.propagator, chain.tracking, chain.measurements.data(), static_cast<int>(chain.measurements.size())).state;
  BOOST_TEST(std::isfinite(reference.qOverP()));
  BOOST_TEST(std::abs(reference.qOverP()) < 1e-3F);
  BOOST_TEST(std::abs(reference.phi() - chain.truth[0].phi()) < 1e-5F);
  BOOST_TEST(std::abs(reference.theta() - chain.truth[0].theta()) < 1e-5F);
  BOOST_TEST(std::abs(reference.loc0() - exact.loc0()) < 5e-5F);
  BOOST_TEST(std::abs(reference.loc1() - exact.loc1()) < 5e-5F);
  BOOST_TEST(std::abs(reference.phi() - exact.phi()) < 5e-6F);
  BOOST_TEST(std::abs(reference.theta() - exact.theta()) < 5e-6F);
  BOOST_TEST(std::abs(reference.qOverP() - exact.qOverP()) < 5e-4F);
}

/// Why the seed takes the whole chain: a 0.4 GeV/c track crosses every barrel
/// layer; seeded from the first three clusters alone, a 16 um measurement
/// error becomes a curvature that misses the outer sensors or linearises the
/// fit about the wrong helix. Every chain must fit, and the pulls must be
/// unit.
///
/// With material, and scattering in the truth to match it. Without material
/// the points fix q/p to 1e-4, which no prior both non-informative and
/// representable in a float covariance can reach: see kChainQOverPSpread. A
/// real track crosses about 1.5e-3 X0 per cm here, what the material LUT
/// reports for the barrel, and is scattering-limited long before that.
///
/// Both passes are checked: the inward state against the mid-plane of the
/// first sensor and the outward one against the mid-plane of the last, with
/// unit pulls, and each pass's spatial chi2 per degree of freedom at one.
/// The outward pass is the road fit's, whose chi2 decides which roads are
/// kept; its chi2 is one only with the process noise placed as the material
/// is.
BOOST_AUTO_TEST_CASE(LowMomentumChainsFitWithUnitPulls)
{
  constexpr Scalar kQOverP = 2.5F;
  constexpr Scalar kX2X0PerCm = 1.5e-3F;
  constexpr std::size_t kLayers = 4;
  BarrelChain<kLayers> chain(kQOverP);
  BOOST_REQUIRE(chain.built);
  const UniformMaterial material{.x2X0PerCm = kX2X0PerCm};

  constexpr int kTrials = 1000;
  const Scalar sigmaLoc = 1.5811388e-3F; // sqrt(2.5e-6) cm

  std::mt19937 rng(20260923);
  std::normal_distribution<double> gauss(0., 1.);
  std::array<std::vector<double>, kNBoundParams> inwardPull;
  std::array<std::vector<double>, kNBoundParams> outwardPull;
  double inwardChi2PerDof = 0.;
  double outwardChi2PerDof = 0.;

  for (int trial = 0; trial < kTrials; ++trial) {
    // A kink can send the track through the 0.5 mm gap between two modules;
    // draw again rather than fit a chain with a missing layer.
    for (int attempt = 0; attempt < 10; ++attempt) {
      chain.generate(kQOverP, kX2X0PerCm, &rng);
      if (chain.built) {
        break;
      }
    }
    BOOST_REQUIRE(chain.built);
    auto smeared = chain.measurements;
    for (auto& m : smeared) {
      m.loc0 += static_cast<Scalar>(gauss(rng) * static_cast<double>(sigmaLoc));
      m.loc1 += static_cast<Scalar>(gauss(rng) * static_cast<double>(sigmaLoc));
    }
    const auto fit = fitMeasurements(chain.propagator, chain.tracking, material, smeared.data(), static_cast<int>(smeared.size()));
    BOOST_TEST_CONTEXT("trial " << trial << ", propagation status " << static_cast<int>(fit.status)) { BOOST_REQUIRE(fit.status == PropagationStatus::kOk); }
    collectPulls(trial, fit.state, chain.truth[0], inwardPull);
    collectPulls(trial, fit.outwardState, chain.truth[kLayers - 1], outwardPull);
    inwardChi2PerDof += static_cast<double>(fit.chi2) / fit.ndf / kTrials;
    outwardChi2PerDof += static_cast<double>(fit.outwardChi2) / fit.ndf / kTrials;
  }

  const char* names[kNBoundParams] = {"loc0", "loc1", "phi", "theta", "qOverP"};
  for (const auto* pass : {&inwardPull, &outwardPull}) {
    for (int i = 0; i < kNBoundParams; ++i) {
      const auto& p = (*pass)[static_cast<std::size_t>(i)];
      const auto n = static_cast<double>(p.size());
      const auto [mean, width] = moments(p);
      BOOST_TEST_CONTEXT((pass == &inwardPull ? "inward" : "outward") << " parameter " << names[i] << " pull mean " << mean << " width " << width)
      {
        BOOST_TEST(std::abs(mean) < 4. / std::sqrt(n));
        BOOST_TEST(std::abs(width - 1.) < 0.10);
      }
    }
  }
  BOOST_TEST_CONTEXT("inward chi2 per spatial dof " << inwardChi2PerDof << ", outward " << outwardChi2PerDof)
  {
    BOOST_TEST(std::abs(inwardChi2PerDof - 1.) < 0.15);
    BOOST_TEST(std::abs(outwardChi2PerDof - 1.) < 0.15);
  }
}

/// The same chains losing energy: the chain seed is a circle through the
/// clusters without energy loss, so a pass that carried it, unfitted, to its
/// first update would charge the seed's miss there against a 1 mm position
/// prior. Each pass starts on its first measurement, the inward one from the
/// outward pass's end state, and both chi2 per dof are one. The position and
/// direction pulls are unit; q/p's is printed, not held: this truth loses
/// energy without straggling, which the fit's q/p noise includes, and is not
/// the fit's model of where the loss is taken.
BOOST_AUTO_TEST_CASE(LowMomentumChainsLosingEnergyFitWithUnitPulls)
{
  constexpr Scalar kQOverP = 2.5F;
  constexpr Scalar kX2X0PerCm = 1.5e-3F;
  constexpr Scalar kRho = 0.3F; // about 2.5 % of the momentum over the chain
  constexpr std::size_t kLayers = 4;
  BarrelChain<kLayers> chain(kQOverP);
  BOOST_REQUIRE(chain.built);
  const UniformMaterial material{.x2X0PerCm = kX2X0PerCm, .rho = kRho};

  constexpr int kTrials = 1000;
  const Scalar sigmaLoc = 1.5811388e-3F;
  std::mt19937 rng(20261021);
  std::normal_distribution<double> gauss(0., 1.);
  std::array<std::vector<double>, kNBoundParams> inwardPull;
  double inwardChi2PerDof = 0.;
  double outwardChi2PerDof = 0.;
  double lostShare = 0.;
  for (int trial = 0; trial < kTrials; ++trial) {
    for (int attempt = 0; attempt < 10; ++attempt) {
      chain.generate(kQOverP, kX2X0PerCm, &rng, kRho);
      if (chain.built) {
        break;
      }
    }
    BOOST_REQUIRE(chain.built);
    lostShare += (1. - static_cast<double>(kQOverP / chain.truth[kLayers - 1].qOverP())) / kTrials;
    auto smeared = chain.measurements;
    for (auto& m : smeared) {
      m.loc0 += static_cast<Scalar>(gauss(rng) * static_cast<double>(sigmaLoc));
      m.loc1 += static_cast<Scalar>(gauss(rng) * static_cast<double>(sigmaLoc));
    }
    const auto fit = fitMeasurements(chain.propagator, chain.tracking, material, smeared.data(), static_cast<int>(smeared.size()));
    BOOST_TEST_CONTEXT("trial " << trial << ", propagation status " << static_cast<int>(fit.status)) { BOOST_REQUIRE(fit.status == PropagationStatus::kOk); }
    collectPulls(trial, fit.state, chain.truth[0], inwardPull);
    inwardChi2PerDof += static_cast<double>(fit.chi2) / fit.ndf / kTrials;
    outwardChi2PerDof += static_cast<double>(fit.outwardChi2) / fit.ndf / kTrials;
  }
  BOOST_TEST_MESSAGE("momentum lost over the chain " << lostShare << ", chi2 per dof inward " << inwardChi2PerDof << " outward " << outwardChi2PerDof);
  BOOST_TEST(lostShare > 0.01);
  for (int i = 0; i < kNBoundParams; ++i) {
    const auto& p = inwardPull[static_cast<std::size_t>(i)];
    const bool held = i != o2::gtrk::kQOverP;
    const auto n = static_cast<double>(p.size());
    const auto [mean, width] = moments(p);
    BOOST_TEST_MESSAGE("inward parameter " << i << " pull mean " << mean << " width " << width);
    BOOST_TEST_CONTEXT("inward parameter " << i << " pull mean " << mean << " width " << width)
    {
      BOOST_TEST((!held || std::abs(mean) < 4. / std::sqrt(n)));
      BOOST_TEST((!held || std::abs(width - 1.) < 0.10));
    }
  }
  BOOST_TEST_CONTEXT("inward chi2 per spatial dof " << inwardChi2PerDof << ", outward " << outwardChi2PerDof)
  {
    BOOST_TEST(std::abs(inwardChi2PerDof - 1.) < 0.15);
    BOOST_TEST(std::abs(outwardChi2PerDof - 1.) < 0.15);
  }
}

BOOST_AUTO_TEST_CASE(PerigeePhiIsTheDirectionAtTheTransverseClosestApproach)
{
  const SyntheticChain chain;
  const auto& state = chain.truth[0];
  const auto& surface = chain.tracking.surface(state.surface);
  // The closest approach by search: the smallest transverse radius on a
  // coarse scan, refined by ternary search.
  double best = 0.;
  for (double path = -100.; path <= 100.; path += 0.1) {
    if (transverseRadius2At(chain.propagator, state, surface, path) < transverseRadius2At(chain.propagator, state, surface, best)) {
      best = path;
    }
  }
  double low = best - 0.1;
  double high = best + 0.1;
  for (int i = 0; i < 200; ++i) {
    const double left = low + ((high - low) / 3.);
    const double right = high - ((high - low) / 3.);
    if (transverseRadius2At(chain.propagator, state, surface, left) < transverseRadius2At(chain.propagator, state, surface, right)) {
      high = right;
    } else {
      low = left;
    }
  }
  const double path = 0.5 * (low + high);
  const double sinTheta = std::sin(static_cast<double>(state.theta()));
  const double turn = chain.propagator.curvature(static_cast<double>(state.qOverP()), sinTheta) * path * sinTheta;
  const double expected = static_cast<double>(state.phi()) + turn;
  const double difference = std::remainder(static_cast<double>(perigeePhi(chain.propagator, state, surface)) - expected, 2. * std::numbers::pi);
  BOOST_TEST(std::abs(difference) < 1e-5);
}

BOOST_AUTO_TEST_CASE(StraightTrackPerigeeHasNoCurvatureSpecialCase)
{
  TrackingSurface surface;
  surface.origin[0] = 2.F;
  surface.origin[1] = 1.F;
  surface.origin[2] = 3.F;
  BoundTrackState state;
  state.par[kPhi] = 0.F;
  state.par[kTheta] = std::numbers::pi_v<Scalar> / 2;
  state.par[kQOverP] = 0.F;

  const Propagator propagator;
  BOOST_TEST(std::abs(perigeePhi(propagator, state, surface)) < 1e-6F);
}

BOOST_AUTO_TEST_CASE(TheEnergyLossIsO2s)
{
  // The Bethe-Bloch, from slow protons to ultra-relativistic electrons.
  for (const float betaGamma : {0.3F, 0.7F, 1.F, 3.5F, 10.F, 100.F, 1000.F, 1e4F}) {
    BOOST_TEST_CONTEXT("beta gamma " << betaGamma)
    {
      BOOST_TEST(betheBlochSolid(betaGamma) == o2::track::BetheBlochSolidOpt<float>(betaGamma), boost::test_tools::tolerance(1e-6F));
    }
  }
  // A step's loss against O2's own correction of a pion: outward loses
  // (x rho < 0), inward regains, the thick step in sub-steps; q/p scales by
  // O2's factor and its variance grows by O2's straggling.
  constexpr float kMass = 0.13957F;
  for (const float momentum : {0.12F, 0.3F, 2.F}) {
    for (const float xRho : {-0.08F, 0.08F, -1.5F}) {
      BOOST_TEST_CONTEXT("p " << momentum << " GeV, x rho " << xRho)
      {
        // At tan(lambda) = 0, O2's q/pT is our q/p.
        o2::track::TrackParCov o2Track(0.F, 0.F, {0.F, 0.F, 0.F, 0.F, 1.F / momentum},
                                       {1e-4F, 0.F, 1e-4F, 0.F, 0.F, 1e-4F, 0.F, 0.F, 0.F, 1e-4F, 0.F, 0.F, 0.F, 0.F, 1e-4F});
        o2Track.setPID(o2::track::PID::Pion);
        const float variance = o2Track.getSigma1Pt2();
        BOOST_REQUIRE(o2Track.correctForMaterial(0.F, xRho));
        const auto loss = energyLoss(1.F / momentum, kMass, xRho);
        BOOST_REQUIRE(loss.ok);
        BOOST_TEST(loss.momentumScale / momentum == o2Track.getQ2Pt(), boost::test_tools::tolerance(2e-6F));
        BOOST_TEST(variance + loss.varQOverP == o2Track.getSigma1Pt2(), boost::test_tools::tolerance(2e-6F));
        BOOST_TEST((xRho < 0.F ? loss.momentumScale > 1.F : loss.momentumScale < 1.F));
      }
    }
  }
  // No material: nothing; stopped: not ok.
  BOOST_TEST(energyLoss(10.F, kMass, 0.F).momentumScale == 1.F);
  BOOST_TEST(!energyLoss(1.F / 0.05F, kMass, -5.F).ok);
  // beta p: below p.
}

} // namespace o2::gtrk

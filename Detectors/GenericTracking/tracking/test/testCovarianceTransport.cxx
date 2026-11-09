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

/// \file testCovarianceTransport.cxx
/// \brief The Jacobian, the covariance transport and the scattering term.
///
/// The Monte Carlo case is the one that matters: a Jacobian can be
/// self-consistently wrong, but a sample covariance built by actually
/// propagating twenty thousand perturbed states cannot agree with it by
/// accident.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK CovarianceTransport
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/CovarianceTransport.h"
#include "GenericTracking/Tracking/Material.h"
#include "Oracles.h"

#include <array>
#include <cmath>
#include <cstring>
#include <random>
#include <vector>

namespace o2::gtrk
{

static TrackingSurface barrelPlane(double x, SensorID id)
{
  TrackingSurface s{};
  s.origin[0] = static_cast<Scalar>(x);
  s.normal[0] = 1.F;
  s.axis0[1] = 1.F; // loc0 around the beam
  s.axis1[2] = 1.F; // loc1 along it
  s.bounds = SurfaceBounds{.halfX = 1e6F, .halfY = 1e6F};
  s.id = id;
  return s;
}

static TrackingSurface discPlane(double z, SensorID id)
{
  TrackingSurface s{};
  s.origin[2] = static_cast<Scalar>(z);
  s.normal[2] = 1.F;
  s.axis0[0] = 1.F;
  s.axis1[1] = 1.F;
  s.bounds = SurfaceBounds{.halfX = 1e6F, .halfY = 1e6F};
  s.id = id;
  return s;
}

static BoundTrackState makeState(SensorID id, double phi, double theta, double qOverP)
{
  BoundTrackState st;
  st.par[kPhi] = static_cast<Scalar>(phi);
  st.par[kTheta] = static_cast<Scalar>(theta);
  st.par[kQOverP] = static_cast<Scalar>(qOverP);
  st.surface = id;
  return st;
}

BOOST_AUTO_TEST_CASE(JacobianPublishesItsCrossingWithoutChangingTheResult)
{
  const Propagator prop;
  const auto from = barrelPlane(3.4, 1);
  std::mt19937 random{4703U};
  std::uniform_real_distribution<double> phi{-1.2, 1.2};
  std::uniform_real_distribution<double> theta{0.15, 2.99};
  std::uniform_real_distribution<double> qOverP{-20., 20.};
  unsigned successful = 0U;
  unsigned failed = 0U;
  for (unsigned trial = 0U; trial < 1000U; ++trial) {
    const auto to = trial % 2U == 0U ? barrelPlane(17.2, 2) : discPlane(trial % 4U == 1U ? 29. : -29., 2);
    const auto seed = makeState(1, phi(random), theta(random), qOverP(random));
    SolenoidJacobian expectedJacobian;
    const auto expectedStatus = jacobianAnalytic(prop, seed, from, to, PropagationDirection::kForward, expectedJacobian);
    auto expectedState = seed;
    const auto expectedResult = propagate(prop, expectedState, from, to, PropagationDirection::kForward);
    auto actualState = seed; // In-place aliasing is used by filterStep.
    SolenoidJacobian actualJacobian;
    PropagationResult actualResult;
    const auto actualStatus = jacobianAnalytic(prop, actualState, from, to, PropagationDirection::kForward, actualJacobian, &actualState, &actualResult);
    BOOST_TEST(static_cast<unsigned>(actualStatus) == static_cast<unsigned>(expectedStatus));
    if (actualStatus == PropagationStatus::kOk) {
      // Bit identity is intentional, including floating-point encodings; these records have no padding.
      // NOLINTNEXTLINE(bugprone-suspicious-memory-comparison)
      BOOST_TEST(std::memcmp(&actualState, &expectedState, sizeof(BoundTrackState)) == 0);
      // Bit identity is intentional, including floating-point encodings; these records have no padding.
      // NOLINTNEXTLINE(bugprone-suspicious-memory-comparison)
      BOOST_TEST(std::memcmp(&actualJacobian, &expectedJacobian, sizeof(SolenoidJacobian)) == 0);
      BOOST_TEST(static_cast<unsigned>(actualResult.status) == static_cast<unsigned>(expectedResult.status));
      BOOST_TEST(exactlyEqual(actualResult.pathLength, expectedResult.pathLength));
      ++successful;
    } else {
      // Bit identity is intentional, including floating-point encodings; these records have no padding.
      // NOLINTNEXTLINE(bugprone-suspicious-memory-comparison)
      BOOST_TEST(std::memcmp(&actualState, &seed, sizeof(BoundTrackState)) == 0);
      BOOST_TEST(static_cast<unsigned>(actualResult.status) == static_cast<unsigned>(PropagationStatus::kNotFinite));
      ++failed;
    }
  }
  BOOST_TEST(successful > 100U);
  BOOST_TEST(failed > 10U);
}

/// Propagate many samples drawn from a diagonal source covariance and compare
/// the sample covariance at the target with F P F^T.
static void checkMonteCarloClosure(const TrackingSurface& from, const TrackingSurface& to, const BoundTrackState& seed, double tolerance, bool analytic = false)
{
  const Propagator prop;
  constexpr int kSamples = 20000;

  const Scalar sigma[kNBoundParams] = {5e-3F, 5e-3F, 1e-3F, 1e-3F, 1e-3F};

  BoundTrackState st = seed;
  for (int i = 0; i < kNBoundParams; ++i) {
    st.variance(i) = sigma[i] * sigma[i];
  }

  Scalar step[kNBoundParams] = {};
  defaultDifferenceSteps(step);
  BoundTrackState transported = st;
  if (analytic) {
    SolenoidJacobian jac;
    BOOST_TEST((jacobianAnalytic(prop, st, from, to, PropagationDirection::kForward, jac) == PropagationStatus::kOk));
    auto factors = factorCovariance(st.cov);
    BOOST_REQUIRE(transportFactors(jac, StepNoise{}, factors));
    covarianceFromFactors(factors, transported.cov);
  } else {
    Scalar jac[kJacobianSize] = {};
    BOOST_TEST((jacobianFiniteDifference(prop, st, from, to, PropagationDirection::kForward, step, jac) == PropagationStatus::kOk));
    transportCovarianceDense(jac, transported.cov);
  }
  BOOST_TEST((propagate(prop, transported, from, to, PropagationDirection::kForward).status == PropagationStatus::kOk));

  std::mt19937 rng(20260919);
  std::normal_distribution<double> gauss(0., 1.);

  std::vector<std::array<double, kNBoundParams>> out;
  out.reserve(kSamples);
  for (int n = 0; n < kSamples; ++n) {
    BoundTrackState s = seed;
    for (int i = 0; i < kNBoundParams; ++i) {
      s.par[i] += static_cast<Scalar>(gauss(rng) * static_cast<double>(sigma[i]));
    }
    if (propagate(prop, s, from, to, PropagationDirection::kForward).status != PropagationStatus::kOk) {
      continue;
    }
    std::array<double, kNBoundParams> row{};
    for (int i = 0; i < kNBoundParams; ++i) {
      row[static_cast<std::size_t>(i)] = static_cast<double>(s.par[i]);
    }
    out.push_back(row);
  }
  BOOST_TEST(out.size() > kSamples / 2);

  const auto n = static_cast<double>(out.size());
  std::array<double, kNBoundParams> mean{};
  for (const auto& row : out) {
    for (int i = 0; i < kNBoundParams; ++i) {
      mean[static_cast<std::size_t>(i)] += row[static_cast<std::size_t>(i)] / n;
    }
  }

  for (int i = 0; i < kNBoundParams; ++i) {
    for (int j = 0; j <= i; ++j) {
      double acc = 0.;
      for (const auto& row : out) {
        acc += (row[static_cast<std::size_t>(i)] - mean[static_cast<std::size_t>(i)]) * (row[static_cast<std::size_t>(j)] - mean[static_cast<std::size_t>(j)]);
      }
      const double sample = acc / (n - 1.);
      const auto predicted = static_cast<double>(transported.covAt(i, j));
      const double scale = std::sqrt(static_cast<double>(transported.variance(i)) * static_cast<double>(transported.variance(j)));
      // Compared against the geometric mean of the two variances, so an
      // off-diagonal element is judged on the same footing as a diagonal one.
      BOOST_TEST(std::abs(sample - predicted) < tolerance * scale);
      // With a correct Jacobian the residual is a per cent or two: the Monte
      // Carlo's own 1/sqrt(2N), about 0.5 %, and the second-order terms a
      // linearised transport does not carry and the sample covariance does.
      // A 10 % error in a dominant Jacobian entry shows up well above 5 %, so
      // the window is comfortably above the floor and well below a real fault.
      // Do not "fix" a residual of a per cent or two here.
    }
  }
}

BOOST_AUTO_TEST_CASE(TransportedCovarianceMatchesAMonteCarloInTheBarrel) { checkMonteCarloClosure(barrelPlane(3.4, 1), barrelPlane(17.2, 2), makeState(1, 0.1, 1.2, 1. / 2.), 0.05); }

BOOST_AUTO_TEST_CASE(TransportedCovarianceMatchesAMonteCarloAtHighEta)
{
  // sin(theta) = 0.05 puts the track on the outermost disc at r = 7.6 cm, the
  // inner ring. This case matters because the curvature carries a
  // 1/sin(theta), so theta and q/p are amplified twentyfold and strongly
  // correlated exactly here -- and nowhere near mid-rapidity, where a
  // barrel-only test would never see it.
  checkMonteCarloClosure(discPlane(0., 1), discPlane(152., 2), makeState(1, 0.3, 0.05, 1. / 5.), 0.05);
}

BOOST_AUTO_TEST_CASE(AnalyticJacobianAgreesWithTheReferenceAcrossTheAcceptance)
{
  // 50 MeV to 100 GeV, both charges, theta from 0.05 (the outermost disc's
  // inner ring) to 2.4, barrel-to-barrel and through the endcap.
  //
  // The tolerance is two-sided on purpose. A pure relative test is meaningless
  // on entries near zero, where the reference is noise: a central difference
  // of a float output carries about 1e-6 cm / 2h, which is 5e-4 at the phi
  // step at normal incidence, and 1 / |n . d| times that at the target -- the
  // float rounding of the start moves the intersection that much more along a
  // plane the track crosses shallowly (0.012 at theta = 0.05 on a barrel
  // plane). A pure absolute test is meaningless on the large entries, which
  // reach 150 at high eta. So: an absolute floor for what the reference cannot
  // resolve, and a relative term for what it can.
  const Propagator prop;
  const TrackingSurface pairs[3][2] = {{barrelPlane(3.4, 1), barrelPlane(7.0, 2)}, {barrelPlane(3.4, 1), barrelPlane(17.2, 2)}, {discPlane(0., 1), discPlane(152., 2)}};

  int checked = 0;
  for (const auto& pair : pairs) {
    for (const double pT : {0.05, 0.1, 0.5, 1.0, 10.0, 100.0}) {
      for (const double theta : {0.05, 0.3, 0.9, 1.5708, 2.4}) {
        for (const double charge : {1., -1.}) {
          for (const double phi : {0.0, 0.7, -2.1}) {
            const BoundTrackState st = makeState(1, phi, theta, charge * std::sin(theta) / pT);
            Scalar step[kNBoundParams] = {};
            defaultDifferenceSteps(step);

            Scalar ref[kJacobianSize] = {};
            SolenoidJacobian sol;
            SolenoidJacobian exactSol;
            Scalar ana[kJacobianSize] = {};
            Scalar exact[kJacobianSize] = {};
            if (jacobianFiniteDifference(prop, st, pair[0], pair[1], PropagationDirection::kForward, step, ref) != PropagationStatus::kOk) {
              continue; // unreachable target: a fact about the track
            }
            BoundTrackState arrived = st;
            BOOST_REQUIRE(propagate(prop, arrived, pair[0], pair[1], PropagationDirection::kForward).status == PropagationStatus::kOk);
            const auto arrivedTheta = static_cast<double>(arrived.theta());
            const auto arrivedPhi = static_cast<double>(arrived.phi());
            const double incidence = std::abs((std::sin(arrivedTheta) * ((std::cos(arrivedPhi) * static_cast<double>(pair[1].normal[0])) + (std::sin(arrivedPhi) * static_cast<double>(pair[1].normal[1])))) + (std::cos(arrivedTheta) * static_cast<double>(pair[1].normal[2])));
            const double floor = 1e-3 / incidence;
            BOOST_TEST((jacobianAnalytic<Scalar>(prop, st, pair[0], pair[1], PropagationDirection::kForward, sol) == PropagationStatus::kOk));
            BOOST_TEST((jacobianAnalytic<ExactScalar>(prop, st, pair[0], pair[1], PropagationDirection::kForward, exactSol) == PropagationStatus::kOk));
            expandJacobian(sol, ana);
            expandJacobian(exactSol, exact);

            for (int i = 0; i < kJacobianSize; ++i) {
              const double diff = std::abs(static_cast<double>(ana[i]) - static_cast<double>(ref[i]));
              const double precisionDiff = std::abs(static_cast<double>(ana[i]) - static_cast<double>(exact[i]));
              BOOST_TEST_CONTEXT("pT " << pT << ", theta " << theta << ", charge " << charge << ", phi " << phi << ", element " << i)
              {
                BOOST_TEST(diff < floor + 5e-3 * std::abs(static_cast<double>(ref[i])));
                BOOST_TEST(precisionDiff < 2e-5 + 5e-5 * std::abs(static_cast<double>(exact[i])));
              }
            }
            ++checked;
          }
        }
      }
    }
  }
  // If a change makes most targets unreachable the loop above would pass
  // vacuously, so pin how much of it actually ran.
  BOOST_TEST(checked > 250);
}

BOOST_AUTO_TEST_CASE(ScalarJacobianMatchesExactNearTheGrazingLimit)
{
  const Propagator prop;
  const auto from = barrelPlane(0., 1);
  const auto to = barrelPlane(1., 2);

  for (const double incidence : {0.0101, 0.011, 0.02, 0.05, 0.06}) {
    const auto phi = std::acos(incidence);
    const BoundTrackState st = makeState(1, phi, 0.5 * std::acos(-1.), 0.);
    SolenoidJacobian scalarSol;
    SolenoidJacobian exactSol;
    const auto scalarStatus = jacobianAnalytic<Scalar>(prop, st, from, to, PropagationDirection::kForward, scalarSol);
    const auto exactStatus = jacobianAnalytic<ExactScalar>(prop, st, from, to, PropagationDirection::kForward, exactSol);

    BOOST_TEST_CONTEXT("incidence " << incidence)
    {
      BOOST_REQUIRE(scalarStatus == PropagationStatus::kOk);
      BOOST_REQUIRE(exactStatus == scalarStatus);
      Scalar scalar[kJacobianSize] = {};
      Scalar exact[kJacobianSize] = {};
      expandJacobian(scalarSol, scalar);
      expandJacobian(exactSol, exact);
      for (int i = 0; i < kJacobianSize; ++i) {
        const auto expected = static_cast<double>(exact[i]);
        const double diff = std::abs(static_cast<double>(scalar[i]) - expected);
        BOOST_TEST(diff < 2e-4 + 1e-4 * std::abs(expected));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(TheAnalyticJacobianAlsoClosesAgainstTheMonteCarlo)
{
  // The reference and the analytic form can only agree with each other by
  // accident if they share a mistake; agreeing with propagated samples as well
  // rules that out.
  checkMonteCarloClosure(barrelPlane(3.4, 1), barrelPlane(17.2, 2), makeState(1, 0.1, 1.2, 1. / 2.), 0.05, true);
  checkMonteCarloClosure(discPlane(0., 1), discPlane(152., 2), makeState(1, 0.3, 0.05, 1. / 5.), 0.05, true);
}

BOOST_AUTO_TEST_CASE(HighlandReproducesASensorByHand)
{
  // 125 um of silicon, X0 = 9.37 cm, so x/X0 = 1.334e-3 and Highland gives
  // 0.372 mrad at 1 GeV. The whole material model rests on this number, so it
  // is checked rather than remembered.
  const Scalar x2X0 = 0.0125F / 9.37F;
  BOOST_TEST(std::abs(static_cast<double>(x2X0) - 1.334e-3) < 1e-5);

  const Scalar t0 = highlandTheta0(x2X0, 1.F);
  BOOST_TEST(std::abs(static_cast<double>(t0) - 3.718e-4) < 1e-6);

  // 1/p scaling, and two sensors -- the overlap case -- 46 % wider.
  BOOST_TEST(std::abs(static_cast<double>(highlandTheta0(x2X0, 10.F)) - 3.718e-5) < 1e-7);
  const auto ratio = static_cast<double>(highlandTheta0(2.F * x2X0, 1.F)) / static_cast<double>(t0);
  BOOST_TEST(std::abs(ratio - 1.46) < 0.01);

  BOOST_TEST(std::abs(static_cast<double>(highlandTheta0(0.F, 1.F))) < 1e-12);
}

BOOST_AUTO_TEST_CASE(ScatteringProjectsOntoTheTwoAngles)
{
  const auto noise = scatteringNoise(1e-3F, 0.5F);
  const double sine = std::sin(0.5);
  BOOST_TEST(std::abs(static_cast<double>(noise.theta) - 1e-6) < 1e-12);
  BOOST_TEST(std::abs(static_cast<double>(noise.phi) - 1e-6 / (sine * sine)) < 1e-12);
}

/// The factored prediction is the double one, F (P + Q_before) F^T +
/// Q_after, to float precision on every element, over the planes and momenta
/// of the structured transport's test, for full covariances and every kind of
/// noise, straggling included.
BOOST_AUTO_TEST_CASE(FactoredPredictionMatchesTheDoubleOne)
{
  const Propagator prop;
  const TrackingSurface pairs[3][2] = {{barrelPlane(3.4, 1), barrelPlane(7.0, 2)}, {barrelPlane(3.4, 1), barrelPlane(17.2, 2)}, {discPlane(0., 1), discPlane(152., 2)}};
  std::mt19937 rng(11);
  std::uniform_real_distribution<double> unit(-1., 1.);
  std::uniform_real_distribution<double> positive(0.1, 1.);
  int checked = 0;
  double worst = 0.;
  for (const auto& pair : pairs) {
    for (const double pT : {0.1, 1.0, 50.0}) {
      for (const double theta : {0.05, 0.9}) {
        const BoundTrackState seed = makeState(1, 0.4, theta, std::sin(theta) / pT);
        SolenoidJacobian jac;
        BOOST_REQUIRE(jacobianAnalytic(prop, seed, pair[0], pair[1], PropagationDirection::kForward, jac) == PropagationStatus::kOk);
        for (int trial = 0; trial < 20; ++trial) {
          // P = L L^T with L lower triangular: positive definite, every
          // element correlated, the scales of a fit's (cm, rad, 1/GeV).
          constexpr std::array<double, kNBoundParams> kScale{1e-2, 1e-2, 1e-3, 1e-3, 1e-2};
          double l[kNBoundParams][kNBoundParams]{};
          for (int i = 0; i < kNBoundParams; ++i) {
            for (int j = 0; j < i; ++j) {
              l[i][j] = 0.7 * unit(rng) * kScale[static_cast<std::size_t>(i)];
            }
            l[i][i] = positive(rng) * kScale[static_cast<std::size_t>(i)];
          }
          Scalar cov[kBoundCovSize]{};
          for (int i = 0; i < kNBoundParams; ++i) {
            for (int j = 0; j <= i; ++j) {
              double sum = 0.;
              for (int k = 0; k <= j; ++k) {
                sum += l[i][k] * l[j][k];
              }
              cov[covIndex(i, j)] = static_cast<Scalar>(sum);
            }
          }
          // A third of the trials without noise, the rest with every kind.
          const auto pick = [&](double scale) { return trial % 3 == 0 ? 0.F : static_cast<Scalar>(positive(rng) * scale); };
          const StepNoise noise{
            .phiBefore = pick(1e-6),
            .thetaBefore = pick(1e-6),
            .phiAfter = pick(1e-6),
            .thetaAfter = pick(1e-6),
            .qOverPAfter = pick(1e-5),
          };

          auto factors = factorCovariance(cov);
          BOOST_REQUIRE(transportFactors(jac, noise, factors));
          Scalar fromFactors[kBoundCovSize]{};
          covarianceFromFactors(factors, fromFactors);
          ExactScalar exact[kBoundCovSize]{};
          predictCovarianceExact(jac, noise, cov, exact);
          for (int i = 0; i < kNBoundParams; ++i) {
            for (int j = 0; j <= i; ++j) {
              const double scale = std::sqrt(exact[covIndex(i, i)] * exact[covIndex(j, j)]);
              const double error = std::abs(static_cast<double>(fromFactors[covIndex(i, j)]) - exact[covIndex(i, j)]) / scale;
              worst = std::max(worst, error);
              BOOST_TEST_CONTEXT("pT " << pT << " theta " << theta << " trial " << trial << " element " << i << j) { BOOST_TEST(error < 1e-5); }
            }
          }
          ++checked;
        }
      }
    }
  }
  BOOST_TEST_MESSAGE("factored prediction against double: " << checked << " predictions, worst element " << worst << " of its sigmas' product");
  BOOST_TEST(checked == 360);
}

/// Why the fit carries factors:
/// a state fixed by one measurement to a 2.9 um pixel, its direction still
/// the diffuse 0.1 rad, carried from one VD layer to the next at eta 3, where
/// loc1 (z) moves by r / sin^2(theta) = 71 cm per radian of theta. The
/// predicted loc1 variance is 50 cm^2; what the next update needs is loc1's
/// variance given the direction and q/p -- the pixel's and the scattering's,
/// about 1e-8 of it at 50 GeV. As factors it is D_1, and it comes out to
/// float's precision on itself. The float P of the same step, factored, has
/// lost it: rounding 50 cm^2 to float costs more than the whole of it.
BOOST_AUTO_TEST_CASE(FactoredPredictionKeepsWhatAFloatCovarianceLoses)
{
  SolenoidJacobian jac;
  jac.loc0[kLoc0] = 1.F;
  jac.loc1[kLoc1] = 1.F;
  jac.loc1[kTheta] = -71.F;
  jac.phi[kPhi] = 1.F;
  Scalar cov[kBoundCovSize]{};
  cov[covIndex(kLoc0, kLoc0)] = 8.33e-8F; // (10 um)^2 / 12
  cov[covIndex(kLoc1, kLoc1)] = 8.33e-8F;
  cov[covIndex(kPhi, kPhi)] = 1e-2F;
  cov[covIndex(kTheta, kTheta)] = 1e-2F;
  cov[covIndex(kQOverP, kQOverP)] = 1e-2F;
  // theta0 = 1.2e-5 each half, 50 GeV through 1e-3 X0/cm over 7 cm.
  const Scalar t2 = 1.2e-5F * 1.2e-5F;
  const Scalar s2 = 0.0993F * 0.0993F; // sin^2 theta at eta 3
  const StepNoise noise{.phiBefore = t2 / s2, .thetaBefore = t2, .phiAfter = t2 / s2, .thetaAfter = t2};

  ExactScalar exact[kBoundCovSize]{};
  predictCovarianceExact(jac, noise, cov, exact);
  ExactScalar exactUpper[kBoundUpperSize]{};
  ExactScalar exactDiagonal[kNBoundParams]{};
  detail::factorUD(exact, exactUpper, exactDiagonal);
  const double conditional = exactDiagonal[kLoc1];
  BOOST_TEST(conditional > 0.);
  BOOST_TEST(conditional / exact[covIndex(kLoc1, kLoc1)] < 1e-7); // beyond a float's 6e-8 of P_11 only by a little

  auto factors = factorCovariance(cov);
  BOOST_REQUIRE(transportFactors(jac, noise, factors));
  BOOST_TEST(std::abs(static_cast<double>(factors.diagonal[kLoc1]) / conditional - 1.) < 1e-5);

  BoundTrackState asState;
  for (int i = 0; i < kBoundCovSize; ++i) {
    asState.cov[i] = cov[i];
  }
  asState.par[kTheta] = 0.0995F;
  asState.variance(kPhi) += noise.phiBefore;
  asState.variance(kTheta) += noise.thetaBefore;
  Scalar dense[kJacobianSize];
  expandJacobian(jac, dense);
  transportCovarianceDense(dense, asState.cov);
  asState.variance(kPhi) += noise.phiAfter;
  asState.variance(kTheta) += noise.thetaAfter;
  const auto rounded = factorCovariance(asState.cov);
  BOOST_TEST_MESSAGE("loc1 given the rest: double " << conditional << ", factors " << factors.diagonal[kLoc1] << ", factored float P " << rounded.diagonal[kLoc1]);
  BOOST_TEST(std::abs(static_cast<double>(rounded.diagonal[kLoc1]) / conditional - 1.) > 0.5);
}

} // namespace o2::gtrk

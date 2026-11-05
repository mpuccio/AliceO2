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

/// \file testMeasurementUpdate.cxx
/// \brief Closure tests for the spatial Kalman update.

#define BOOST_TEST_MODULE Test GTRK Measurement Update
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <numbers>

#include "GenericTracking/Tracking/MeasurementUpdate.h"
#include "Oracles.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <random>

namespace o2::gtrk
{
namespace
{

BoundTrackState makePrior()
{
  BoundTrackState state;
  state.par[kPhi] = 0.2F;
  state.par[kTheta] = 1.1F;
  state.par[kQOverP] = 0.5F;
  state.variance(kLoc0) = 4.F;
  state.variance(kLoc1) = 9.F;
  state.variance(kPhi) = 0.04F;
  state.variance(kTheta) = 0.09F;
  state.variance(kQOverP) = 0.16F;
  return state;
}

Measurement makeMeasurement()
{
  Measurement measurement;
  measurement.loc0 = 1.F;
  measurement.loc1 = -2.F;
  measurement.varLoc0 = 1.F;
  measurement.varLoc1 = 4.F;
  return measurement;
}

} // namespace

BOOST_AUTO_TEST_CASE(UpdateMatchesTheScalarDiagonalAnswer)
{
  auto state = makePrior();
  const auto result = updateWithMeasurement(state, makeMeasurement());

  BOOST_TEST(std::abs(static_cast<double>(state.loc0()) - 0.8) < 1e-6);
  BOOST_TEST(std::abs(static_cast<double>(state.loc1()) + 18. / 13.) < 1e-6);
  BOOST_TEST(std::abs(static_cast<double>(state.variance(kLoc0)) - 0.8) < 1e-6);
  BOOST_TEST(std::abs(static_cast<double>(state.variance(kLoc1)) - 36. / 13.) < 1e-6);

  BOOST_TEST(std::abs(static_cast<double>(result.chi2) - (1. / 5. + 4. / 13.)) < 1e-6);

  BOOST_TEST(state.phi() == 0.2F);
  BOOST_TEST(state.theta() == 1.1F);
  BOOST_TEST(state.qOverP() == 0.5F);
}

BOOST_AUTO_TEST_CASE(CorrelationsUpdateAllFiveParameters)
{
  auto state = makePrior();
  state.covAt(kLoc1, kLoc0) = 0.5F;
  state.covAt(kPhi, kLoc0) = 0.1F;
  state.covAt(kTheta, kLoc1) = -0.2F;
  state.covAt(kQOverP, kLoc0) = 0.3F;

  auto measurement = makeMeasurement();
  measurement.covLoc01 = -0.25F;

  const auto before = state;
  const auto result = updateWithMeasurement(state, measurement);

  BOOST_TEST(state.phi() != before.phi());
  BOOST_TEST(state.theta() != before.theta());
  BOOST_TEST(state.qOverP() != before.qOverP());
  BOOST_TEST(std::abs(static_cast<double>(state.loc0()) - 0.752646775746) < 1e-6);
  BOOST_TEST(std::abs(static_cast<double>(state.loc1()) + 1.316650625602) < 1e-6);
  BOOST_TEST(std::abs(static_cast<double>(result.chi2) - 0.523580365736) < 1e-6);

  for (int i = 0; i < kNBoundParams; ++i) {
    BOOST_TEST(state.variance(i) > 0.F);
    for (int j = 0; j < kNBoundParams; ++j) {
      BOOST_TEST(state.covAt(i, j) == state.covAt(j, i));
    }
  }
}

BOOST_AUTO_TEST_CASE(FactorizedScalarUpdateMatchesExactOracle)
{
  // The spatial prior is millions of times wider than the measurement, the
  // cancellation regime in which a direct float covariance update fails.
  // Build a correlated positive-definite covariance as L L^T so every
  // state parameter participates without relying on hand-tuned entries.
  constexpr Scalar lower[kNBoundParams][kNBoundParams] = {
    {3.F, 0.F, 0.F, 0.F, 0.F},
    {0.4F, 2.F, 0.F, 0.F, 0.F},
    {0.05F, -0.02F, 0.1F, 0.F, 0.F},
    {-0.03F, 0.04F, 0.01F, 0.12F, 0.F},
    {0.8F, -0.5F, 0.2F, -0.1F, 4.F},
  };

  BoundTrackState scalar;
  scalar.par[kLoc0] = 0.2F;
  scalar.par[kLoc1] = -0.4F;
  scalar.par[kPhi] = 0.3F;
  scalar.par[kTheta] = 1.2F;
  scalar.par[kQOverP] = 2.5F;
  for (int row = 0; row < kNBoundParams; ++row) {
    for (int column = 0; column <= row; ++column) {
      Scalar value = 0.F;
      for (int k = 0; k < kNBoundParams; ++k) {
        value += lower[row][k] * lower[column][k];
      }
      scalar.covAt(row, column) = value;
    }
  }
  auto exact = scalar;

  Measurement measurement;
  measurement.loc0 = 0.202F;
  measurement.loc1 = -0.397F;
  measurement.varLoc0 = 2.5e-6F;
  measurement.varLoc1 = 1.5e-6F;
  measurement.covLoc01 = 5.e-7F;

  const auto scalarResult = updateWithMeasurement(scalar, measurement);
  const auto exactResult = updateWithMeasurementExact(exact, measurement);

  for (int i = 0; i < kNBoundParams; ++i) {
    const Scalar scale = std::max(1.F, std::abs(exact.par[i]));
    BOOST_TEST_CONTEXT("parameter " << i) { BOOST_TEST(std::abs(scalar.par[i] - exact.par[i]) < 2e-5F * scale); }
  }

  for (int row = 0; row < kNBoundParams; ++row) {
    for (int column = 0; column <= row; ++column) {
      const Scalar expected = exact.covAt(row, column);
      const Scalar covarianceScale = std::sqrt(exact.variance(row) * exact.variance(column));
      const Scalar tolerance = row == column ? 2e-5F * std::abs(expected) : 1e-4F * covarianceScale;
      BOOST_TEST_CONTEXT("covariance (" << row << ", " << column << ")") { BOOST_TEST(std::abs(scalar.covAt(row, column) - expected) < tolerance); }
    }
  }

  BOOST_TEST(std::abs(scalarResult.chi2 - exactResult.chi2) < 2e-5F);
  BOOST_TEST(scalar.variance(kLoc0) > 0.F);
  BOOST_TEST(scalar.variance(kLoc1) > 0.F);
  BOOST_TEST(scalar.variance(kQOverP) > 0.F);
}

/// The innovation covariance's determinant (the follow's EndDensity cuts on
/// it) from the two sequential
/// innovations, against the determinant of V = H P H^T + R in double for the
/// covariance P = U D U^T the factors hold -- the follow carries the factors,
/// not P. Over the follow's range: position sigmas of 10 um to 1 mm,
/// correlated with each other up to 0.99 and with the direction and q/p;
/// measurement sigmas of 5 to 300 um, correlated up to 0.95. The float
/// product is the determinant algebraically, not exactly: within a relative
/// 1e-5 of the double one here, so a cut on it decides as the double one
/// except that close to the cut. Priors far more degenerate lose more, about
/// float's epsilon over their conditioning.
BOOST_AUTO_TEST_CASE(TheInnovationsMultiplyToTheJointDeterminant)
{
  std::mt19937 random(0x4454524bU);
  std::uniform_real_distribution<double> unit(0., 1.);
  std::normal_distribution<double> normal;
  double worst = 0.;
  for (int sample = 0; sample < 200000; ++sample) {
    // The prior as L L^T in double, factored in double, the factors in float.
    const double sigma0 = std::pow(10., -3. + (2. * unit(random)));
    const double sigma1 = std::pow(10., -3. + (2. * unit(random)));
    const double rho = (unit(random) < 0.5 ? -1. : 1.) * (1. - std::pow(10., -2. * unit(random)));
    double lower[kNBoundParams][kNBoundParams]{};
    lower[kLoc0][kLoc0] = sigma0;
    lower[kLoc1][kLoc0] = sigma1 * rho;
    lower[kLoc1][kLoc1] = sigma1 * std::sqrt(1. - (rho * rho));
    for (int row = kPhi; row < kNBoundParams; ++row) {
      const double scale = std::pow(10., -4. + (3. * unit(random)));
      for (int column = 0; column < row; ++column) {
        lower[row][column] = 0.5 * scale * normal(random);
      }
      lower[row][row] = scale;
    }
    double packed[kBoundCovSize]{};
    for (int row = 0; row < kNBoundParams; ++row) {
      for (int column = 0; column <= row; ++column) {
        double value = 0.;
        for (int k = 0; k < kNBoundParams; ++k) {
          value += lower[row][k] * lower[column][k];
        }
        packed[covIndex(row, column)] = value;
      }
    }
    double upper[kBoundUpperSize]{};
    double diagonal[kNBoundParams]{};
    detail::factorUD(packed, upper, diagonal);
    CovarianceFactors factors;
    for (std::size_t i = 0; i < std::size(upper); ++i) {
      factors.upper[i] = static_cast<Scalar>(upper[i]);
    }
    for (int i = 0; i < kNBoundParams; ++i) {
      factors.diagonal[i] = static_cast<Scalar>(diagonal[i]);
    }
    // The oracle: P from the float factors, in double.
    const auto covariance = [&](int i, int j) {
      double sum = 0.;
      for (int k = std::max(i, j); k < kNBoundParams; ++k) {
        const double ui = k == i ? 1. : static_cast<double>(factors.upper[detail::strictUpperIndex(i, k)]);
        const double uj = k == j ? 1. : static_cast<double>(factors.upper[detail::strictUpperIndex(j, k)]);
        sum += ui * static_cast<double>(factors.diagonal[k]) * uj;
      }
      return sum;
    };
    Measurement measurement;
    const double m0 = std::pow(10., -3.3 + (1.8 * unit(random)));
    const double m1 = std::pow(10., -3.3 + (1.8 * unit(random)));
    measurement.varLoc0 = static_cast<Scalar>(m0 * m0);
    measurement.varLoc1 = static_cast<Scalar>(m1 * m1);
    measurement.covLoc01 = static_cast<Scalar>(0.95 * ((2. * unit(random)) - 1.) * m0 * m1);
    const double v00 = covariance(kLoc0, kLoc0) + static_cast<double>(measurement.varLoc0);
    const double v11 = covariance(kLoc1, kLoc1) + static_cast<double>(measurement.varLoc1);
    const double v01 = covariance(kLoc0, kLoc1) + static_cast<double>(measurement.covLoc01);
    const double expected = (v00 * v11) - (v01 * v01);
    BoundTrackState state;
    state.par[kTheta] = 1.2F;
    state.par[kQOverP] = 1.5F;
    Scalar determinant = 0.F;
    const auto result = updateWithMeasurement<true>(state, factors, measurement, &determinant);
    BOOST_REQUIRE(updateUsable(result, factors));
    worst = std::max(worst, std::abs(static_cast<double>(determinant) - expected) / expected);
  }
  BOOST_TEST_MESSAGE("determinant: worst relative difference " << worst);
  BOOST_TEST(worst < 1e-5);
}

BOOST_AUTO_TEST_CASE(DirectionIsCanonicalAfterTheUpdate)
{
  auto state = makePrior();
  state.par[kPhi] = std::numbers::pi_v<Scalar> - 0.05F;
  state.par[kTheta] = 0.02F;
  state.covAt(kTheta, kLoc0) = 0.2F;
  state.covAt(kPhi, kLoc0) = 0.2F;

  auto measurement = makeMeasurement();
  measurement.loc0 = -1.F;
  measurement.loc1 = 0.F;

  const auto result = updateWithMeasurement(state, measurement);

  BOOST_TEST(result.chi2 > 0.F);
  BOOST_TEST(state.phi() >= -std::numbers::pi_v<Scalar>);
  BOOST_TEST(state.phi() <= std::numbers::pi_v<Scalar>);
  BOOST_TEST(state.theta() >= 0.F);
  BOOST_TEST(state.theta() <= std::numbers::pi_v<Scalar>);
  BOOST_TEST(state.covAt(kTheta, kLoc0) < 0.F);
}

BOOST_AUTO_TEST_CASE(SyntheticResidualsAndPosteriorPullsClose)
{
  constexpr int kSamples = 50000;
  constexpr double kPriorSigmaLoc0 = 2.;
  constexpr double kPriorSigmaLoc1 = 3.;
  constexpr double kMeasurementSigmaLoc0 = 1.;
  constexpr double kMeasurementSigmaLoc1 = 2.;

  std::mt19937 random(0x4754524bU);
  std::normal_distribution<double> normal;

  std::array<double, 4> sum{};
  std::array<double, 4> sum2{};
  for (int sample = 0; sample < kSamples; ++sample) {
    const double truth0 = kPriorSigmaLoc0 * normal(random);
    const double truth1 = kPriorSigmaLoc1 * normal(random);

    auto state = makePrior();
    auto measurement = makeMeasurement();
    measurement.loc0 = static_cast<Scalar>(truth0 + kMeasurementSigmaLoc0 * normal(random));
    measurement.loc1 = static_cast<Scalar>(truth1 + kMeasurementSigmaLoc1 * normal(random));

    const auto result = updateWithMeasurement(state, measurement);
    const std::array<double, 4> values{
      static_cast<double>(result.residualLoc0) / std::sqrt(kPriorSigmaLoc0 * kPriorSigmaLoc0 + kMeasurementSigmaLoc0 * kMeasurementSigmaLoc0),
      static_cast<double>(result.residualLoc1) / std::sqrt(kPriorSigmaLoc1 * kPriorSigmaLoc1 + kMeasurementSigmaLoc1 * kMeasurementSigmaLoc1),
      (static_cast<double>(state.loc0()) - truth0) / std::sqrt(static_cast<double>(state.variance(kLoc0))),
      (static_cast<double>(state.loc1()) - truth1) / std::sqrt(static_cast<double>(state.variance(kLoc1))),
    };

    for (std::size_t i = 0; i < values.size(); ++i) {
      sum[i] += values[i];
      sum2[i] += values[i] * values[i];
    }
  }

  for (std::size_t i = 0; i < sum.size(); ++i) {
    const double mean = sum[i] / kSamples;
    const double width = std::sqrt(sum2[i] / kSamples - mean * mean);
    BOOST_TEST(std::abs(mean) < 0.015);
    BOOST_TEST(std::abs(width - 1.) < 0.015);
  }
}

/// A predicted covariance that is not positive definite gives a negative
/// spatial chi2, which the road fit's cut would pass and the quality order
/// would rank first; `updateUsable` refuses it, and accepts a sound update.
BOOST_AUTO_TEST_CASE(UpdateUsableRefusesAnIndefinitePrior)
{
  const auto measurement = makeMeasurement();

  BoundTrackState sound = makePrior();
  auto soundFactors = factorCovariance(sound.cov);
  const auto good = updateWithMeasurement(sound, soundFactors, measurement);
  BOOST_TEST(updateUsable(good, soundFactors));

  BoundTrackState indefinite = makePrior();
  // The loc0/loc1 block with a correlation beyond one: the innovation's
  // spatial determinant, (4 + 1)(9 + 4) - 10^2, is negative, so its marginal
  // chi2 is too.
  indefinite.covAt(kLoc0, kLoc1) = 10.F;
  auto indefiniteFactors = factorCovariance(indefinite.cov);
  const auto bad = updateWithMeasurement(indefinite, indefiniteFactors, measurement);
  BOOST_TEST(!bad.positiveDefinite);
  BOOST_TEST(!(bad.chi2 >= 0.F));
  BOOST_TEST(!updateUsable(bad, indefiniteFactors));
  BOOST_TEST(good.positiveDefinite);

  // A posterior variance that is not positive is refused on its own.
  BoundTrackState collapsed = makePrior();
  auto collapsedFactors = factorCovariance(collapsed.cov);
  const auto fine = updateWithMeasurement(collapsed, collapsedFactors, measurement);
  collapsedFactors.diagonal[kQOverP] = 0.F;
  BOOST_TEST(!updateUsable(fine, collapsedFactors));
}

} // namespace o2::gtrk

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

/// \file PullStatistics.h
/// \brief Pull collection, widths, reports and closure checks for sensor chains.
#ifndef GTRK_TEST_PULLSTATISTICS_H_
#define GTRK_TEST_PULLSTATISTICS_H_
#include <boost/test/unit_test.hpp>
#include "GenericTracking/Tracking/TrackFitter.h"
#include "MathUtils/detail/trigonometric.h"
#include <fmt/format.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <string>
#include <vector>
namespace o2::gtrk::test::pulls
{
/// One pull distribution: its plain width (the sample standard deviation),
/// its robust width (1.4826 times the median absolute deviation), and the
/// fraction beyond 5 -- a tail shows in the first and the last, not the second.
struct PullWidth {
  double plain = 0.;
  double robust = 0.;
  double beyond5 = 0.;
};

inline PullWidth widthOf(std::vector<double> pulls)
{
  PullWidth w;
  const auto n = static_cast<double>(pulls.size());
  double mean = 0.;
  for (const double p : pulls) {
    mean += p / n;
  }
  double var = 0.;
  for (const double p : pulls) {
    var += (p - mean) * (p - mean);
    w.beyond5 += std::abs(p) > 5. ? 1. / n : 0.;
  }
  w.plain = std::sqrt(var / (n - 1.));
  const auto middle = pulls.begin() + static_cast<std::ptrdiff_t>(pulls.size() / 2);
  std::ranges::nth_element(pulls, middle);
  const double median = *middle;
  for (double& p : pulls) {
    p = std::abs(p - median);
  }
  std::ranges::nth_element(pulls, middle);
  w.robust = 1.4826 * *middle;
  return w;
}

/// The pulls of the fitted state against truth[0] for loc0, loc1, phi, theta
/// and q/p.
struct PullPoint {
  double ptGeV = 0.;
  double eta = 0.;
  int fits = 0;
  int failed = 0;
  int lostVariance = 0; ///< fits with a variance that is not positive
  std::array<PullWidth, 5> width{};
};

constexpr std::array<const char*, 5> kNames{"loc0", "loc1", "phi", "theta", "q/p"};
/// Collect fitted-state pulls; detector-specific generation and smearing
/// live in the trial callback so that each scan keeps its random sequence.
template <typename Chain, typename Trial>
PullPoint collectPulls(const Chain& chain, double ptGeV, double eta, int trials, const Trial& trialFit)
{
  std::array<std::vector<double>, 5> pulls;
  PullPoint point{.ptGeV = ptGeV, .eta = eta};
  for (int trial = 0; trial < trials; ++trial) {
    const auto fit = trialFit(trial);
    if (fit.status != PropagationStatus::kOk) {
      ++point.failed;
      continue;
    }
    ++point.fits;
    bool positive = true;
    for (int i = 0; i < kNBoundParams; ++i) {
      positive = positive && fit.state.variance(i) > 0.F;
    }
    if (!positive) {
      ++point.lostVariance;
      continue;
    }
    for (int i = 0; i < kNBoundParams; ++i) {
      double diff = static_cast<double>(fit.state.par[i]) - static_cast<double>(chain.truth[0].par[i]);
      if (i == kPhi) {
        math_utils::detail::bringToPMPi<double>(diff);
      }
      pulls[static_cast<std::size_t>(i)].push_back(diff / std::sqrt(static_cast<double>(fit.state.variance(i))));
    }
  }
  for (std::size_t i = 0; i < pulls.size(); ++i) {
    point.width[i] = widthOf(pulls[i]);
  }
  return point;
}

inline std::string report(const std::string& title, std::span<const PullPoint> points, bool tails = true)
{
  std::string out = fmt::format("{}: plain / robust pull width{}\n", title, tails ? ", % beyond 5" : "");
  out += fmt::format("{:>5} {:>5} {:>5} {:>4}", "pT", "eta", "fits", "lost");
  for (const auto* name : kNames) {
    out += fmt::format(" {:>{}}", name, tails ? 19 : 13);
  }
  out += '\n';
  for (const auto& p : points) {
    out += fmt::format("{:5.2f} {:5.1f} {:5d} {:4d}", p.ptGeV, p.eta, p.fits, p.lostVariance);
    for (const auto& w : p.width) {
      out += fmt::format("   {:5.2f} / {:4.2f}", w.plain, w.robust);
      if (tails) {
        out += fmt::format(" {:4.1f}", 100. * w.beyond5);
      }
    }
    out += '\n';
  }
  return out;
}

inline void checkUnit(std::span<const PullPoint> points, double tolerance)
{
  for (const auto& p : points) {
    BOOST_TEST_CONTEXT("pT " << p.ptGeV << " eta " << p.eta)
    {
      BOOST_TEST(p.failed == 0);
      BOOST_TEST(p.lostVariance == 0);
      for (std::size_t i = 0; i < kNames.size(); ++i) {
        BOOST_TEST_CONTEXT(kNames[i] << " plain " << p.width[i].plain << " robust " << p.width[i].robust)
        {
          BOOST_TEST(std::abs(p.width[i].plain - 1.) < tolerance);
          BOOST_TEST(std::abs(p.width[i].robust - 1.) < tolerance);
        }
      }
    }
  }
}
} // namespace o2::gtrk::test::pulls
#endif // GTRK_TEST_PULLSTATISTICS_H_

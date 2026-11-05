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

#include "GenericTracking/Study/TrackPulls.h"
#include "GenericTracking/Study/Statistics.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <span>
#include <utility>

namespace o2::gtrk::study
{

namespace
{

/// 1.4826 times the median absolute deviation: a Gaussian's sigma, blind to
/// a tail. Reorders \a values.
double robustSigma(std::span<double> values)
{
  const double centre = median(values);
  std::vector<double> deviations;
  deviations.reserve(values.size());
  for (const double value : values) {
    deviations.push_back(std::abs(value - centre));
  }
  return 1.4826 * median(deviations);
}

/// Their copies: the medians reorder them.
std::string summary(std::vector<double> residuals, std::vector<double> pulls, double scale)
{
  if (residuals.empty()) {
    return "     --";
  }
  std::string text = std::format("{:7} residual median {:8.2f} sigma {:8.2f}", residuals.size(), scale * median(residuals), scale * robustSigma(residuals));
  if (!pulls.empty()) {
    const auto tail = std::ranges::count_if(pulls, [](double pull) { return std::abs(pull) > 5.; });
    text += std::format("; pull median {:6.3f} sigma {:6.3f} beyond 5 {:5.2f} %", median(pulls), robustSigma(pulls),
                        percent(static_cast<std::uint64_t>(tail), pulls.size()));
  }
  return text;
}

} // namespace

PullTable::PullTable(std::string name, std::vector<PullParameter> parameters)
  : mName(std::move(name)), mParameters(std::move(parameters)), mSamples(mParameters.size())
{
}

void PullTable::add(double pt, std::span<const double> residuals, std::span<const double> sigmas)
{
  if (pt < kPtEdges.front() || residuals.size() != mParameters.size() || sigmas.size() != mParameters.size()) {
    return;
  }
  const auto bin = std::min<std::size_t>(static_cast<std::size_t>(std::ranges::upper_bound(kPtEdges, pt) - kPtEdges.begin()) - 1U, kPtBins - 1U);
  ++mTracks;
  for (std::size_t parameter = 0; parameter < mParameters.size(); ++parameter) {
    auto& samples = mSamples[parameter][bin];
    samples.residuals.push_back(residuals[parameter]);
    if (sigmas[parameter] > 0.) {
      samples.pulls.push_back(residuals[parameter] / sigmas[parameter]);
    }
  }
}

std::vector<std::string> PullTable::report() const
{
  std::vector<std::string> lines;
  lines.push_back(std::format("  {}: {} good tracks", mName, mTracks));
  for (std::size_t parameter = 0; parameter < mParameters.size(); ++parameter) {
    const auto& p = mParameters[parameter];
    std::vector<double> allResiduals;
    std::vector<double> allPulls;
    for (std::size_t bin = 0; bin < kPtBins; ++bin) {
      const auto& samples = mSamples[parameter][bin];
      allResiduals.insert(allResiduals.end(), samples.residuals.begin(), samples.residuals.end());
      allPulls.insert(allPulls.end(), samples.pulls.begin(), samples.pulls.end());
      lines.push_back(std::format("    {:>8} [{}] pT [{:g}, {:g}): {}", p.name, p.unit, kPtEdges[bin], kPtEdges[bin + 1U],
                                  summary(samples.residuals, samples.pulls, p.scale)));
    }
    lines.push_back(std::format("    {:>8} [{}] all pT:      {}", p.name, p.unit, summary(allResiduals, allPulls, p.scale)));
  }
  return lines;
}

} // namespace o2::gtrk::study

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

/// \file TrackPulls.h
/// \brief Residuals and pulls of fitted track parameters against MC truth,
///        by parameter and pT.
///
/// Robust statistics throughout: a pull's width is 1.4826 times its median
/// absolute deviation, so a tail of wrong assignments does not hide the core,
/// and the tail is reported on its own as the fraction beyond 5 sigma.

#ifndef GTRK_STUDY_TRACKPULLS_H_
#define GTRK_STUDY_TRACKPULLS_H_

#include <array>
#include <cstddef>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace o2::gtrk::study
{

/// One fitted parameter as a table shows it.
struct PullParameter {
  std::string name;
  std::string unit;  ///< of the residual as reported
  double scale = 1.; ///< from the parameter's own unit to `unit`
};

class PullTable
{
 public:
  static constexpr std::array<double, 4> kPtEdges{0.1, 0.3, 1., std::numeric_limits<double>::infinity()};
  static constexpr std::size_t kPtBins = kPtEdges.size() - 1U;

  PullTable(std::string name, std::vector<PullParameter> parameters);

  /// One track: each parameter's residual (fitted less true) and its fitted
  /// sigma; a sigma of zero records the residual without a pull.
  void add(double pt, std::span<const double> residuals, std::span<const double> sigmas);

  /// One line per parameter and pT bin, and one for all pT.
  [[nodiscard]] std::vector<std::string> report() const;

 private:
  struct Samples {
    std::vector<double> residuals;
    std::vector<double> pulls;
  };

  std::string mName;
  std::vector<PullParameter> mParameters;
  std::vector<std::array<Samples, kPtBins>> mSamples; ///< per parameter, per pT bin
  std::size_t mTracks = 0U;
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_TRACKPULLS_H_

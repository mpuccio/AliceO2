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

/// \file Statistics.h
/// \brief The median, the quantile and the percentage the study's reports share.

#ifndef GTRK_STUDY_STATISTICS_H_
#define GTRK_STUDY_STATISTICS_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>

namespace o2::gtrk::study
{

/// The median of \a values, the upper of the two middle ones for an even
/// count, 0 for none. Reorders \a values (nth_element); pass a copy to keep
/// their order.
[[nodiscard]] inline double median(std::span<double> values)
{
  if (values.empty()) {
    return 0.;
  }
  const auto middle = values.begin() + static_cast<std::ptrdiff_t>(values.size() / 2U);
  std::ranges::nth_element(values, middle);
  return *middle;
}

/// The element at rank floor(\a q (n - 1)) of \a values (n of them, not
/// empty: the caller says what an empty set gives), \a q in [0, 1]: 0 the
/// smallest, 1 the largest. Reorders \a values (nth_element).
template <typename T>
[[nodiscard]] T quantile(std::span<T> values, double q)
{
  const auto at = values.begin() + static_cast<std::ptrdiff_t>(q * static_cast<double>(values.size() - 1U));
  std::ranges::nth_element(values, at);
  return *at;
}

/// \a part as a percentage of \a whole, 0 for none.
[[nodiscard]] inline double percent(std::uint64_t part, std::uint64_t whole)
{
  return whole == 0U ? 0. : 100. * static_cast<double>(part) / static_cast<double>(whole);
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_STATISTICS_H_

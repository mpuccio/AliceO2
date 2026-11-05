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

#include "GenericTracking/Base/CapacityEstimator.h"

#include <algorithm>
#include <cmath>

namespace o2::gtrk
{

namespace
{

struct Config {
  float alpha{0.2F}; ///< the ratio's exponential smoothing
  float marginInit{1.30F};
  float marginMin{1.10F};
  float marginMax{4.00F};
  float marginUp{1.50F};        ///< the most a timeframe that grew raises the margin by
  float marginDown{0.98F};      ///< what a timeframe without growth lowers it by
  std::size_t floorSlots{1024}; ///< per thread: before anything is learned, and never below
};
constexpr Config kConfig{};

} // namespace

std::size_t CapacityEstimator::perThread(std::size_t site, std::size_t rows) const
{
  const auto& entry = mEntries[site];
  if (!entry.learned || rows == 0U) {
    return kConfig.floorSlots;
  }
  const double raw = static_cast<double>(entry.ratio) * static_cast<double>(rows) * static_cast<double>(entry.margin);
  return std::max(kConfig.floorSlots, static_cast<std::size_t>(std::ceil(raw)));
}

void CapacityEstimator::learn(std::size_t site, std::size_t largestUse, std::size_t granted, std::size_t rows, bool grew)
{
  if (rows == 0U) {
    return;
  }
  auto& entry = mEntries[site];
  const auto sample = static_cast<float>(static_cast<double>(largestUse) / static_cast<double>(rows));
  const bool first = !entry.learned;
  entry.ratio = first ? sample : (kConfig.alpha * sample) + ((1.F - kConfig.alpha) * entry.ratio);
  if (first) {
    entry.margin = kConfig.marginInit;
  } else if (grew) {
    // By the shortfall: a little at least, marginUp at most.
    const float shortfall = granted != 0U ? static_cast<float>(static_cast<double>(largestUse) / static_cast<double>(granted)) : kConfig.marginUp;
    entry.margin = std::min(kConfig.marginMax, entry.margin * std::clamp(shortfall, 1.02F, kConfig.marginUp));
  } else {
    entry.margin = std::max(kConfig.marginMin, entry.margin * kConfig.marginDown);
  }
  entry.learned = true;
}

} // namespace o2::gtrk

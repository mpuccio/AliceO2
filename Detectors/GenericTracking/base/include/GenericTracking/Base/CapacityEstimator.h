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

/// \file CapacityEstimator.h
/// \brief Per stage, the staging each thread's buffer is given at a
///        timeframe's start, learned from the timeframes before.

#ifndef GTRK_BASE_CAPACITYESTIMATOR_H_
#define GTRK_BASE_CAPACITYESTIMATOR_H_

#include <array>
#include <cstddef>

namespace o2::gtrk
{

/// What each stage's per-thread staging buffer is given at a timeframe's
/// start: the largest use any
/// of its buffers had over a timeframe, per measurement row, smoothed over the
/// timeframes before, times a margin that rises when a buffer grew during a
/// timeframe and falls when none did. Within a timeframe the buffers grow as
/// they must and are not predicted again. One entry per stage, a site
/// index below `kSites` (the tracker's `StagingSite`), read and written by
/// the tracker's thread alone.
///
/// ITSMFTTracking's CapacityEstimator is the same model, keyed, locked and
/// transactional, in a library this one cannot link.
class CapacityEstimator
{
 public:
  static constexpr std::size_t kSites = 4U;

  /// The records to stage per thread for \a site in a timeframe of \a rows
  /// measurement rows.
  [[nodiscard]] std::size_t perThread(std::size_t site, std::size_t rows) const;
  /// Learn from a timeframe of \a rows rows: \a largestUse, the most records
  /// one of \a site's buffers held at the end of a run; \a granted, the
  /// records per thread it was given; \a grew, whether a buffer grew past
  /// them.
  void learn(std::size_t site, std::size_t largestUse, std::size_t granted, std::size_t rows, bool grew);

 private:
  struct Entry {
    float ratio{0.F}; ///< records per measurement row
    float margin{0.F};
    bool learned{false};
  };
  std::array<Entry, kSites> mEntries{};
};

} // namespace o2::gtrk

#endif /* GTRK_BASE_CAPACITYESTIMATOR_H_ */

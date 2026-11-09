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

#ifndef GTRK_TEST_STAGEOUTPUTCOMPARISON_H_
#define GTRK_TEST_STAGEOUTPUTCOMPARISON_H_

#include <cstddef>
#include <cstring>

namespace o2::gtrk
{

/// Compare the element-to-output mapping shared by every StageOutput.
template <typename Output>
[[nodiscard]] bool sameStageShape(const Output& lhs, const Output& rhs)
{
  return lhs.size() == rhs.size() && lhs.counts().size() == rhs.counts().size() && lhs.offsets().size() == rhs.offsets().size() && std::memcmp(lhs.counts().data(), rhs.counts().data(), lhs.counts().size_bytes()) == 0 && std::memcmp(lhs.offsets().data(), rhs.offsets().data(), lhs.offsets().size_bytes()) == 0;
}

/// Compare counts, offsets and every schema column byte for byte.
template <typename Output>
[[nodiscard]] bool sameStageOutput(const Output& lhs, const Output& rhs)
{
  if (!sameStageShape(lhs, rhs)) {
    return false;
  }

  const auto left = lhs.view();
  const auto right = rhs.view();
  bool same = true;
  left.apply([&](const auto&... leftColumn) {
    right.apply([&](const auto&... rightColumn) {
      static_assert(sizeof...(leftColumn) == sizeof...(rightColumn));
      same = ((std::memcmp(leftColumn, rightColumn, static_cast<std::size_t>(lhs.size()) * sizeof(*leftColumn)) == 0) && ...);
    });
  });
  return same;
}

} // namespace o2::gtrk

#endif // GTRK_TEST_STAGEOUTPUTCOMPARISON_H_

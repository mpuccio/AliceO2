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

/// \file Lanes.h
/// \brief A lane's share of a walk that several lanes make together, and
///        the lanes that make one in step.

#ifndef GTRK_BASE_LANES_H_
#define GTRK_BASE_LANES_H_

#include "GPUCommonDef.h"

#include <cstdint>

namespace o2::gtrk
{

/// A share of one walk's rows: lane `lane` of `count` visits the `lane`-th row
/// of every range and each `count`-th after it, so the lanes together visit
/// every row once. One lane of one visits them all in order: the walk of a
/// single thread.
struct Lanes {
  std::uint32_t lane = 0U;
  std::uint32_t count = 1U;
};

/// The lanes of a group walk in step: each takes
/// its `Lanes` share of a walk's rows, and every lane makes every step, so
/// that the group can decide together after each. A group type gives its
/// width `kLanes`, the lane's index `lane()`, `ballot(set)` -- bit `i` for
/// lane `i` -- and `broadcast(value, from)`, lane `from`'s value of a
/// trivially copyable type on every lane. Every lane of a group calls every
/// collective, in the same order.
///
/// This one is a single lane: the host's, and a one-thread kernel's. Its
/// `lane()` is 0 and `ballot(set)` is `set`; `broadcast` only groups of
/// several lanes need. The devices' groups of several lanes are `WarpGroup`
/// (LaneGroup.cuh).
struct SoloGroup {
  static constexpr std::uint32_t kLanes = 1U;
  // Members, not static: the interface is the device groups', whose lanes are
  // their state; a static one would leave `group.lane()` a static call through
  // an instance.
  // NOLINTBEGIN(readability-convert-member-functions-to-static)
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint32_t lane() const { return 0U; }
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint64_t ballot(bool set) const { return set ? 1U : 0U; }
  // NOLINTEND(readability-convert-member-functions-to-static)
};

} // namespace o2::gtrk

#endif // GTRK_BASE_LANES_H_

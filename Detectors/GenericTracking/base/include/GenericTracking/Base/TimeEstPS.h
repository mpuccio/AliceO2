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

/// \file TimeEstPS.h
/// \brief A measurement's time as a box in picoseconds: O2's `its::TimeEstBC`
///        with picoseconds for bunch crossings.
///
/// A box is `[start, start + width)`, both unsigned 32-bit picoseconds from
/// the timeframe's start. A 32-orbit timeframe is 2.85e9 ps of the 4.29e9 the
/// type holds. Nothing in a timeframe is before its start, so nothing stored
/// is negative.
///
/// A bunch crossing is 10^10 / 400789 ps = 24950.78 ps (O2's
/// `LHCBunchSpacingNS`), not a whole number of picoseconds. A box made from BC
/// takes its start and its length each rounded down (`psOfBC`): every frame of
/// one length is equally wide, two boxes overlap exactly when their BC
/// intervals do -- the rounding, under a picosecond, cannot undo an overlap of
/// a BC nor make one -- and two frames that touch in BC leave at most a
/// picosecond between them. A point time in such a gap lies in neither box:
/// a vertex's containment is quantised to the picosecond, not exact.

#ifndef GTRK_BASE_TIMEESTPS_H_
#define GTRK_BASE_TIMEESTPS_H_

#include "GenericTracking/Base/Definitions.h"

#include <cstdint>

namespace o2::gtrk
{

/// A bunch crossing is `kPSPerBCNumerator / kPSPerBCDenominator` ps, exactly
/// (O2's LHCBunchSpacingNS, a double, is not).
inline constexpr std::int64_t kPSPerBCNumerator = 10'000'000'000;
inline constexpr std::int64_t kPSPerBCDenominator = 400'789;

/// \a bc bunch crossings in picoseconds, rounded down, \a bc >= 0.
GPUhdi() constexpr std::int64_t psOfBC(std::int64_t bc)
{
  return (bc * kPSPerBCNumerator) / kPSPerBCDenominator;
}

/// The bunch crossings whose `psOfBC` is \a ps, \a ps >= 0: exact for every
/// value `psOfBC` returns, which lies less than a picosecond -- far less than
/// a BC -- below the exact time.
GPUhdi() constexpr std::int64_t bcOfPS(std::int64_t ps)
{
  return ((ps * kPSPerBCDenominator) + kPSPerBCNumerator - 1) / kPSPerBCNumerator;
}

/// A time box `[start, start + width)` in picoseconds, O2's `its::TimeEstBC`.
struct TimeEstPS {
  std::uint32_t start = 0U;
  std::uint32_t width = 0U;

  [[nodiscard]] GPUhdi() constexpr std::uint32_t lower() const { return start; }

  /// The box's end, saturating at the type's range as O2's.
  [[nodiscard]] GPUhdi() constexpr std::uint32_t upper() const
  {
    constexpr std::uint32_t kMax = std::numeric_limits<std::uint32_t>::max();
    return start > kMax - width ? kMax : start + width;
  }

  /// The boxes overlap.
  [[nodiscard]] GPUhdi() constexpr bool isCompatible(const TimeEstPS& other) const { return upper() > other.lower() && other.upper() > lower(); }

  /// Narrows to the intersection with \a other, which must overlap.
  GPUhdi() constexpr TimeEstPS& operator+=(const TimeEstPS& other)
  {
    const std::uint32_t lo = lower() > other.lower() ? lower() : other.lower();
    const std::uint32_t hi = upper() < other.upper() ? upper() : other.upper();
    start = lo;
    width = hi - lo;
    return *this;
  }
};

/// The box of a readout frame `[startBC, endBC)` widened by \a marginBC at
/// both ends and clamped at the timeframe's start in BC, as O2's
/// `getROFTimeBounds(rof, true)`; its start and its length rounded down.
GPUhdi() constexpr TimeEstPS frameBox(std::int64_t startBC, std::int64_t endBC, std::int64_t marginBC)
{
  const std::int64_t lowerBC = startBC > marginBC ? startBC - marginBC : 0;
  const std::int64_t upperBC = endBC + marginBC;
  return {.start = static_cast<std::uint32_t>(psOfBC(lowerBC)), .width = static_cast<std::uint32_t>(psOfBC(upperBC - lowerBC))};
}

} // namespace o2::gtrk

#endif // GTRK_BASE_TIMEESTPS_H_

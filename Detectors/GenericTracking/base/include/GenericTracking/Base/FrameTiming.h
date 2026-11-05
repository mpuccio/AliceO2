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

/// \file FrameTiming.h
/// \brief A layer's readout frames, which its time boxes, the grid's time
///        blocks and the slices' reach follow.

#ifndef GTRK_BASE_FRAMETIMING_H_
#define GTRK_BASE_FRAMETIMING_H_

#include "GenericTracking/Base/TimeEstPS.h"

#include <cstdint>

namespace o2::gtrk
{

/// Frame `k` of a layer starts at `firstStartBC + k * lengthBC` BC from the
/// timeframe's first orbit: O2's `ROFLookupTables` layer timing
/// (`getROFStartInBC`), the delay and the bias together -- for the ITS and
/// the MFT from the detector's `DPLAlpideParam` in the CCDB, for ALICE 3's
/// TRK from `TRKAlmiraParam`. The frames' records do not hold the bias:
/// frame `k`'s record is at `k * lengthBC`, in the MC and in data alike, and
/// O2's digitiser puts a collision at BC `t` into the frame of `t - bias`.
/// Every layer has frames: the grid refuses a length that is not positive.
/// O2's `ROFTimingLayer` (ITSMFTTracking) holds the same frames, in its ITS/MFT
/// library; this one also serves ALICE 3's TRK.
struct FrameTiming {
  std::int32_t lengthBC = 0;
  std::int32_t firstStartBC = 0;

  /// Whether \a box is exactly one of these frames, as `frameBox` makes it
  /// with no margin: its start a frame's start and its width the length, each
  /// rounded down from BC. A box widened by a margin is not.
  [[nodiscard]] constexpr bool isFrame(const TimeEstPS& box) const
  {
    const std::int64_t startBC = bcOfPS(box.start);
    return startBC >= firstStartBC && (startBC - firstStartBC) % lengthBC == 0 && box.start == psOfBC(startBC) && box.width == psOfBC(lengthBC);
  }

  friend constexpr bool operator==(const FrameTiming&, const FrameTiming&) = default;
};

} // namespace o2::gtrk

#endif // GTRK_BASE_FRAMETIMING_H_

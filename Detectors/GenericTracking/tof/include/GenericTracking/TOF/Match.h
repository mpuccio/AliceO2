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

#ifndef GTRK_TOF_MATCH_H_
#define GTRK_TOF_MATCH_H_

#include <array>
#include <cstdint>
#include <limits>

namespace o2::gtrk::tof
{
/// Fine time in ns from the beginning of an integer, absolute BC. Normalized
/// to [0, LHCBunchSpacingNS). Never flatten an absolute timestamp into float.
struct Timestamp {
  std::uint64_t bc = 0;
  float fineNS = 0.F;
};

/// Signed time difference, subtracting integer BCs before converting to double.
double differenceNS(Timestamp a, Timestamp b);
Timestamp shifted(Timestamp time, double deltaNS);

inline constexpr std::size_t kHypotheses = 5;
inline constexpr std::array<int, kHypotheses> kPDG{11, 13, 211, 321, 2212};
inline constexpr float kUnset = std::numeric_limits<float>::quiet_NaN();

struct TimeHypothesis {
  Timestamp time;               ///< time at the track's innermost fitted point
  float sigmaNS = kUnset;       ///< measurement and propagated-track uncertainty
  float innerFlightNS = kUnset; ///< flight from that point to inner TOF
  float outerFlightNS = kUnset; ///< flight from that point to outer TOF
  float pairChi2 = kUnset;      ///< inner/outer timing compatibility
  bool valid = false;
};

/// One record per input track, in input order. Indices address the original
/// TF3/CLUSTERS array of this timeframe; -1 means unmatched. Time has the same
/// spatial reference as TrackRecord<TRKFT3>::state, including for secondaries.
struct Match {
  std::int32_t inner = -1;
  std::int32_t outer = -1;
  float innerSpatialChi2 = kUnset;
  float outerSpatialChi2 = kUnset;
  std::uint32_t innerCandidates = 0;
  std::uint32_t outerCandidates = 0;
  std::array<TimeHypothesis, kHypotheses> hypotheses{};
};
} // namespace o2::gtrk::tof
#endif

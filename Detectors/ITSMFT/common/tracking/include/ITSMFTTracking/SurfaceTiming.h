// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACETIMING_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACETIMING_H_

#include <cstdint>
#include <limits>
#include <type_traits>

#include "GPUCommonDef.h"

#ifndef GPUCA_GPUCODE
#include "CommonDataFormat/InteractionRecord.h"
#endif

namespace o2::itsmft::tracking
{

// TimeFrame-relative bunch-crossing coordinate. Host timing remains signed
// 64-bit; a future 32-bit device rebase must be an explicit, checked step
// (Architecture.md 7.2), not implemented here.
using TFBC = int64_t;

// Half-open [begin, end) TF-relative BC interval for one source ROF.
// sourceROF is source-local: cross-source compatibility must use interval
// intersection (see intersects()), never equal ROF ordinals.
struct ROFIntervalBC {
  TFBC begin{0};
  TFBC end{0};
  uint32_t sourceROF{std::numeric_limits<uint32_t>::max()};
  uint32_t flags{0};

  GPUhdi() constexpr bool isValid() const noexcept { return begin <= end; }
  GPUhdi() constexpr TFBC length() const noexcept { return end - begin; }
};

static_assert(std::is_standard_layout_v<ROFIntervalBC>);
static_assert(std::is_trivially_copyable_v<ROFIntervalBC>);

// Per-source readout timing configuration. Field names and signs mirror
// o2::its::LayerTiming (ITStracking/ROFLookupTables.h): the ROF start is the
// source ROF's own real InteractionRecord (so both continuous and triggered
// readout are supported without a rofLength*ordinal formula) plus rofDelay
// plus rofBias; rofAddTimeErr is the additional uncertainty margin, applied
// only on demand via widen() -- never folded into the stored interval, so
// that uncertainty handling stays an explicit choice of the compatibility
// policy (Architecture.md 7.2), matching LayerTiming::getROFTimeBounds(bool).
struct ROFTimingConfig {
  TFBC rofLength{0};
  TFBC rofDelay{0};
  TFBC rofBias{0};
  TFBC rofAddTimeErr{0};
};

static_assert(std::is_standard_layout_v<ROFTimingConfig>);
static_assert(std::is_trivially_copyable_v<ROFTimingConfig>);

enum class TimingBuildError : uint8_t {
  None,
  InvalidROFLength,
  Overflow
};

struct ROFIntervalBuildResult {
  ROFIntervalBC interval{};
  TimingBuildError error{TimingBuildError::None};

  constexpr bool ok() const noexcept { return error == TimingBuildError::None; }
};

#ifndef GPUCA_GPUCODE

namespace detail
{
inline bool checkedAddBC(TFBC a, TFBC b, TFBC& out) noexcept
{
  if (b >= 0) {
    if (a > std::numeric_limits<TFBC>::max() - b) {
      return false;
    }
  } else {
    if (a < std::numeric_limits<TFBC>::min() - b) {
      return false;
    }
  }
  out = a + b;
  return true;
}
} // namespace detail

// origin is the single explicit InteractionRecord chosen for the loaded
// frame; rofIR is the source ROF's own real interaction record (e.g.
// ROFRecord::getBCData()). Continuous and triggered readout are both
// supported because rofIR is taken as-is per ROF rather than re-derived from
// a rofLength * ordinal formula.
inline ROFIntervalBuildResult computeROFIntervalBC(const o2::InteractionRecord& rofIR,
                                                   const o2::InteractionRecord& origin,
                                                   const ROFTimingConfig& cfg,
                                                   uint32_t sourceROF) noexcept
{
  if (cfg.rofLength <= 0) {
    return {{}, TimingBuildError::InvalidROFLength};
  }
  const TFBC anchor = static_cast<TFBC>(rofIR.differenceInBC(origin));
  TFBC begin{0};
  TFBC withDelay{0};
  TFBC end{0};
  if (!detail::checkedAddBC(anchor, cfg.rofDelay, withDelay) ||
      !detail::checkedAddBC(withDelay, cfg.rofBias, begin) ||
      !detail::checkedAddBC(begin, cfg.rofLength, end)) {
    return {{}, TimingBuildError::Overflow};
  }
  return {ROFIntervalBC{begin, end, sourceROF, 0}, TimingBuildError::None};
}

// Widen an interval symmetrically by an explicit margin (e.g. rofAddTimeErr)
// on demand; uncertainty is never folded in automatically.
inline ROFIntervalBC widen(const ROFIntervalBC& interval, TFBC margin) noexcept
{
  return ROFIntervalBC{interval.begin - margin, interval.end + margin, interval.sourceROF, interval.flags};
}

#endif // GPUCA_GPUCODE

// Cross-source compatibility uses interval intersection, never equal ROF
// ordinals. Half-open intervals: adjacent (touching) intervals do not
// intersect.
GPUhdi() constexpr bool intersects(const ROFIntervalBC& a, const ROFIntervalBC& b) noexcept
{
  return a.begin < b.end && b.begin < a.end;
}

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACETIMING_H_ */

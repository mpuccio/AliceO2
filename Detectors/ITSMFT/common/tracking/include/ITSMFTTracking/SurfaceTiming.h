// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACETIMING_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACETIMING_H_

#include <cstddef>
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

  // The default-constructed interval (begin==end==0, sourceROF==invalid) is
  // deliberately the invalid sentinel. A genuinely stored ROF interval must
  // carry a real source ROF ordinal and a strictly positive half-open
  // extent; a zero-length or reversed interval is not valid even if its
  // sourceROF looks real.
  GPUhdi() constexpr bool isValid() const noexcept
  {
    return sourceROF != std::numeric_limits<uint32_t>::max() && begin < end;
  }
  // Signed `end - begin` can overflow TFBC for a valid interval spanning
  // more than INT64_MAX BC (e.g. begin close to INT64_MIN, end close to
  // INT64_MAX). Computing the subtraction in uint64_t is well-defined
  // (unsigned arithmetic never overflows) and reproduces the exact
  // mathematical distance for every begin <= end pair representable in
  // TFBC, since that distance is always within [0, UINT64_MAX].
  GPUhdi() constexpr uint64_t length() const noexcept
  {
    return static_cast<uint64_t>(end) - static_cast<uint64_t>(begin);
  }
};

static_assert(std::is_standard_layout_v<ROFIntervalBC>);
static_assert(std::is_trivially_copyable_v<ROFIntervalBC>);

// A track/cell/road-level TF-relative BC time estimate, expressed as the
// same half-open [begin, end) convention as ROFIntervalBC, but without a
// source-ROF identity: this represents a merged estimate potentially
// combining several sources/ROFs, not one source's own single ROF. Standard-
// layout and trivially-copyable (unlike o2::its::TimeEstBC, whose own
// TimeStampWithError/TimeStamp base hierarchy declares non-static data
// members at more than one level and is therefore not standard-layout):
// this is the common, detector-neutral type CommonTrack
// (ITSMFTTracking/CommonTrack.h) uses for its own timestamp field.
struct CommonTrackTimestamp {
  TFBC begin{0};
  TFBC end{0};

  GPUhdi() constexpr bool isValid() const noexcept { return begin < end; }
  // Half-open interval intersection test, matching ROFIntervalBC's
  // intersects() semantics exactly: adjacent (touching) intervals do not
  // intersect, and an invalid interval never intersects anything.
  GPUhdi() constexpr bool isCompatible(const CommonTrackTimestamp& other) const noexcept
  {
    return isValid() && other.isValid() && begin < other.end && other.begin < end;
  }
};

static_assert(std::is_standard_layout_v<CommonTrackTimestamp>);
static_assert(std::is_trivially_copyable_v<CommonTrackTimestamp>);
static_assert(sizeof(CommonTrackTimestamp) == 16);
static_assert(alignof(CommonTrackTimestamp) == alignof(TFBC));
static_assert(offsetof(CommonTrackTimestamp, begin) == 0);
static_assert(offsetof(CommonTrackTimestamp, end) == 8);

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
  InvalidSourceROF,
  Overflow
};

struct ROFIntervalBuildResult {
  ROFIntervalBC interval{};
  TimingBuildError error{TimingBuildError::None};

  // A default-constructed (never-assigned) result must not read as ok(): its
  // `interval` is the invalid sentinel even though `error` defaults to None.
  // Every genuine success path builds an interval with a real sourceROF and
  // begin < end, so requiring interval.isValid() only excludes the
  // uninitialized/default case, never a real success.
  constexpr bool ok() const noexcept { return error == TimingBuildError::None && interval.isValid(); }
};

enum class WidenError : uint8_t {
  None,
  InvalidInterval,
  InvalidMargin,
  LowerBoundOverflow,
  UpperBoundOverflow
};

struct WidenResult {
  ROFIntervalBC interval{};
  WidenError error{WidenError::None};

  // Same self-consistency invariant as ROFIntervalBuildResult::ok(): a
  // default-constructed result must not read as ok() while its `interval` is
  // still the invalid sentinel.
  constexpr bool ok() const noexcept { return error == WidenError::None && interval.isValid(); }
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
  if (sourceROF == std::numeric_limits<uint32_t>::max()) {
    return {{}, TimingBuildError::InvalidSourceROF};
  }
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

// Widen an interval symmetrically by an explicit, non-negative margin (e.g.
// rofAddTimeErr) on demand; uncertainty is never folded in automatically.
// Every failure mode is reported rather than silently wrapping: an invalid
// input interval, a negative margin, or overflow in either bound.
inline WidenResult widen(const ROFIntervalBC& interval, TFBC margin) noexcept
{
  if (!interval.isValid()) {
    return {{}, WidenError::InvalidInterval};
  }
  if (margin < 0) {
    return {{}, WidenError::InvalidMargin};
  }
  TFBC newBegin{0};
  if (!detail::checkedAddBC(interval.begin, -margin, newBegin)) {
    return {{}, WidenError::LowerBoundOverflow};
  }
  TFBC newEnd{0};
  if (!detail::checkedAddBC(interval.end, margin, newEnd)) {
    return {{}, WidenError::UpperBoundOverflow};
  }
  return {ROFIntervalBC{newBegin, newEnd, interval.sourceROF, interval.flags}, WidenError::None};
}

#endif // GPUCA_GPUCODE

// Cross-source compatibility uses interval intersection, never equal ROF
// ordinals. Half-open intervals: adjacent (touching) intervals do not
// intersect. An invalid interval never intersects anything.
GPUhdi() constexpr bool intersects(const ROFIntervalBC& a, const ROFIntervalBC& b) noexcept
{
  return a.isValid() && b.isValid() && a.begin < b.end && b.begin < a.end;
}

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACETIMING_H_ */

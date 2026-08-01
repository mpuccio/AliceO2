// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.

#ifndef ALICEO2_ITSMFT_TRACKING_CLOCKTIMINGPUBLICATIONVIEW_H_
#define ALICEO2_ITSMFT_TRACKING_CLOCKTIMINGPUBLICATIONVIEW_H_

#ifndef GPUCA_GPUCODE

#include <cstdint>
#include <limits>
#include <optional>

#include "ITSMFTTracking/SurfaceTiming.h"
#include "ITStracking/ROFLookupTables.h"

namespace o2::itsmft::tracking
{

// Host-only immutable output boundary around the established clock-layer
// implementation. It deliberately delegates symmetry, clamp and ROF lookup
// to LayerTiming/TimeEstBC; no second timing formula belongs in adapters.
class ClockTimingPublicationView
{
 public:
  explicit ClockTimingPublicationView(const o2::its::LayerTiming& clock) : mClock{clock} {}

  std::optional<o2::its::TimeEstBC> makeTimeEstBC(const CommonTrackTimestamp& timestamp) const noexcept
  {
    if (!timestamp.isValid() || timestamp.begin < 0 || timestamp.end < 0 ||
        timestamp.begin > std::numeric_limits<uint32_t>::max() || timestamp.end > std::numeric_limits<uint32_t>::max()) {
      return std::nullopt;
    }
    const auto width = static_cast<uint64_t>(timestamp.end) - static_cast<uint64_t>(timestamp.begin);
    if (width > std::numeric_limits<uint16_t>::max()) {
      return std::nullopt;
    }
    return o2::its::TimeEstBC{static_cast<uint32_t>(timestamp.begin), static_cast<uint16_t>(width)};
  }

  std::optional<o2::its::TimeStamp> makeOutputTimestamp(const CommonTrackTimestamp& timestamp) const noexcept
  {
    const auto asymmetric = makeTimeEstBC(timestamp);
    if (!asymmetric) {
      return std::nullopt;
    }
    auto symmetric = asymmetric->makeSymmetrical();
    const float clamp = mClock.mROFLength * 0.5f;
    if (symmetric.getTimeStampError() > clamp) {
      symmetric.setTimeStampError(clamp);
    }
    return symmetric;
  }

  int getROF(const o2::its::TimeStamp& timestamp) const noexcept { return mClock.getROF(timestamp); }
  uint32_t getROFCount() const noexcept { return mClock.mNROFsTF; }
  const o2::its::LayerTiming& getLegacyClockLayer() const noexcept { return mClock; }

 private:
  o2::its::LayerTiming mClock;
};

} // namespace o2::itsmft::tracking

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_CLOCKTIMINGPUBLICATIONVIEW_H_

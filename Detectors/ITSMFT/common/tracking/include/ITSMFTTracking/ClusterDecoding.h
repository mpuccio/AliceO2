// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_CLUSTERDECODING_H_
#define ALICEO2_ITSMFT_TRACKING_CLUSTERDECODING_H_

#include <cstddef>
#include <cstdint>

#include <gsl/gsl>

#include "DataFormatsITSMFT/ClusterPattern.h"

namespace o2::itsmft::tracking
{

// Typed failures at the host compact-cluster decoding boundary. These values
// are intentionally independent of MultiSourceLoadError: a loader maps them
// while adding source/ROF/external-cluster context.
enum class ClusterDecodeError : uint8_t {
  None,
  MissingDictionary,
  TruncatedExplicitPattern,
  MalformedExplicitPattern,
  InvalidPatternId,
  InvalidSensor,
  InvalidLayer,
  InvalidLayerMapping,
  GeometryUnavailable,
  OtherMalformedInput
};

// Host-only cursor for the source-local explicit-pattern byte stream. The
// cursor owns no storage and always retains both the current position and the
// end through a span plus an offset. It validates the ClusterPattern encoding
// (row byte, column byte, then ceil(row*column/8) bitmap bytes) before calling
// ClusterPattern's iterator constructor, which is itself unbounded.
class BoundedPatternCursor
{
 public:
  explicit BoundedPatternCursor(gsl::span<const unsigned char> bytes) noexcept : mBytes(bytes) {}

  size_t consumed() const noexcept { return mPosition; }
  size_t remaining() const noexcept { return mBytes.size() - mPosition; }
  bool empty() const noexcept { return remaining() == 0; }

  ClusterDecodeError acquirePattern(o2::itsmft::ClusterPattern& pattern) noexcept
  {
    const auto available = remaining();
    if (available < 2) {
      return ClusterDecodeError::TruncatedExplicitPattern;
    }

    const auto rowSpan = mBytes[mPosition];
    const auto columnSpan = mBytes[mPosition + 1];
    if (rowSpan == 0 || columnSpan == 0 ||
        rowSpan > o2::itsmft::ClusterPattern::MaxRowSpan ||
        columnSpan > o2::itsmft::ClusterPattern::MaxColSpan) {
      return ClusterDecodeError::MalformedExplicitPattern;
    }

    const size_t nBits = static_cast<size_t>(rowSpan) * columnSpan;
    const size_t payloadBytes = (nBits + 7) / 8;
    const size_t encodedBytes = 2 + payloadBytes;
    if (available < encodedBytes) {
      return ClusterDecodeError::TruncatedExplicitPattern;
    }

    auto iterator = mBytes.begin() + mPosition;
    o2::itsmft::ClusterPattern decoded{iterator};
    if (decoded.getNPixels() == 0) {
      return ClusterDecodeError::MalformedExplicitPattern;
    }
    pattern = decoded;
    mPosition += encodedBytes;
    return ClusterDecodeError::None;
  }

 private:
  gsl::span<const unsigned char> mBytes{};
  size_t mPosition{0};
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_CLUSTERDECODING_H_ */

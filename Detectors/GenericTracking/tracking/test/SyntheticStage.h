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

/// \file SyntheticStage.h
/// \brief The synthetic stage of testStageOutput as a host-device functor, so
///        the serial reference on the host and GpuStage on a device run one
///        and the same per-element function.

#ifndef GTRK_TEST_SYNTHETICSTAGE_H_
#define GTRK_TEST_SYNTHETICSTAGE_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/Tracklet.h"

#include <cstdint>

namespace o2::gtrk::synthetic
{

/// How many records each element emits.
enum class Kind : std::uint8_t {
  kSample,    ///< 1 2 0 2 0: the canonical-output example of testStageOutput
  kUneven,    ///< six for the first 40 elements, then none: dense then empty
  kIrregular, ///< (7 e) mod 11
  kOverflow,  ///< one element emits one record more than an ordinal can number
};

struct Count {
  Kind kind = Kind::kIrregular;

  [[nodiscard]] GPUhdi() std::uint32_t operator()(std::uint32_t element) const
  {
    switch (kind) {
      case Kind::kSample:
        if (element == 0U) {
          return 1U;
        }
        return element == 1U || element == 3U ? 2U : 0U;
      case Kind::kUneven:
        return element < 40U ? 6U : 0U;
      case Kind::kIrregular:
        return (element * 7U) % 11U;
      case Kind::kOverflow:
        return element == 1U ? static_cast<std::uint32_t>(StageOutput<TrackletT>::kMaxPerElement) + 1U : 1U;
    }
    return 0U;
  }
};

/// Record k of element e: both in every column, so a record placed anywhere
/// but offset[e] + k is caught.
[[nodiscard]] GPUhdi() Tracklet syntheticRecord(std::uint32_t element, std::uint32_t k)
{
  return Tracklet{
    .source = element,
    .target = (element * 1000U) + k,
    .link = static_cast<std::uint16_t>(element % 5U),
    .lambda = static_cast<Scalar>(element) + (static_cast<Scalar>(k) * 0.125F),
    .phi = -static_cast<Scalar>(element) - (static_cast<Scalar>(k) * 0.25F),
  };
}

/// testStageOutput's Generator.
template <typename Counts>
struct CountGenerator {
  Counts count;
  template <typename Emitter>
  GPUhdi() void operator()(std::uint32_t element, Emitter& emitter) const
  {
    const std::uint32_t records = count(element);
    for (std::uint32_t k = 0; k < records; ++k) {
      emitter.emit(syntheticRecord(element, k));
    }
  }
};

using Generator = CountGenerator<Count>;

/// The flat form's synthetic test: a
/// candidate row is accepted unless it is 1 modulo \a divisor; 0 accepts all.
struct Accept {
  std::uint32_t divisor = 3U;
  [[nodiscard]] GPUhdi() bool operator()(std::uint32_t row) const { return divisor == 0U || row % divisor != 1U; }
};

/// The serial reference of the flat form: element e's candidates are the
/// rows e * 1000 + k, k below `count(e)`, in that order; the accepted ones
/// are emitted as `syntheticRecord(e, k)`.
struct FlatGenerator {
  Count count;
  Accept accept;
  template <typename Emitter>
  GPUhdi() void operator()(std::uint32_t element, Emitter& emitter) const
  {
    const std::uint32_t candidates = count(element);
    for (std::uint32_t k = 0; k < candidates; ++k) {
      if (accept((element * 1000U) + k)) {
        emitter.emit(syntheticRecord(element, k));
      }
    }
  }
};

} // namespace o2::gtrk::synthetic

#endif // GTRK_TEST_SYNTHETICSTAGE_H_

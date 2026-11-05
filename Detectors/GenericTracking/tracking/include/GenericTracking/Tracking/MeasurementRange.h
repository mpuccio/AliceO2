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

/// \file MeasurementRange.h
/// \brief A half-open range of measurement rows.

#ifndef GTRK_TRACKING_MEASUREMENTRANGE_H_
#define GTRK_TRACKING_MEASUREMENTRANGE_H_

#include "GenericTracking/Base/Definitions.h"

#include <cstdint>

namespace o2::gtrk
{

struct MeasurementRange {
  std::uint32_t begin = 0;
  std::uint32_t end = 0;

  [[nodiscard]] GPUhdi() GPUAlwaysInline() constexpr std::uint32_t size() const { return end - begin; }
  [[nodiscard]] GPUhdi() GPUAlwaysInline() constexpr bool empty() const { return begin == end; }
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_MEASUREMENTRANGE_H_

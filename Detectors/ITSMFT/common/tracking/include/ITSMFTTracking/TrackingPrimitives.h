// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKINGPRIMITIVES_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKINGPRIMITIVES_H_

#include "DataFormatsITS/TimeEstBC.h"
#include "GPUCommonDef.h"
#include "ITStracking/Constants.h"
#include "ITStracking/MathUtils.h"

#include <type_traits>

namespace o2::itsmft::tracking
{

// Event-owned locator sorted into the per-surface spatial index. Physics and
// identity remain authoritative in GlobalMeasurement; this type contains only
// the values needed by the hot candidate-search loops.
struct MeasurementLocator {
  GPUhdDefault() MeasurementLocator() = default;
  GPUhd() MeasurementLocator(float x, float y, float z, int measurementIndex)
    : xCoordinate{x}, yCoordinate{y}, zCoordinate{z}, phi{o2::its::math_utils::computeNormalizedPhi(x, y)}, radius{o2::its::math_utils::hypot(x, y)}, clusterId{measurementIndex}, indexTableBinIndex{0}
  {
  }

  float xCoordinate{0.f};
  float yCoordinate{0.f};
  float zCoordinate{0.f};
  float phi{0.f};
  float radius{0.f};
  int clusterId{o2::its::constants::UnusedIndex};
  int indexTableBinIndex{o2::its::constants::UnusedIndex};
};

// Per-iteration connection between two sorted measurement locators.
struct Tracklet {
  GPUhdDefault() Tracklet() = default;
  GPUhd() Tracklet(int first, int second, float tanL, float azimuth, const o2::its::TimeEstBC& time)
    : firstClusterIndex{first}, secondClusterIndex{second}, tanLambda{tanL}, phi{azimuth}, mTime{time}
  {
  }

  GPUhd() bool operator<(const Tracklet& other) const noexcept
  {
    return firstClusterIndex != other.firstClusterIndex ? firstClusterIndex < other.firstClusterIndex
                                                        : secondClusterIndex < other.secondClusterIndex;
  }
  GPUhd() bool operator==(const Tracklet& other) const noexcept
  {
    return firstClusterIndex == other.firstClusterIndex && secondClusterIndex == other.secondClusterIndex;
  }
  GPUhd() bool isCompatible(const Tracklet& other) const { return mTime.isCompatible(other.mTime); }
  GPUhd() auto& getTimeStamp() noexcept { return mTime; }
  GPUhd() const auto& getTimeStamp() const noexcept { return mTime; }

  int firstClusterIndex{o2::its::constants::UnusedIndex};
  int secondClusterIndex{o2::its::constants::UnusedIndex};
  float tanLambda{o2::its::constants::UnsetValue};
  float phi{o2::its::constants::UnsetValue};
  o2::its::TimeEstBC mTime;
};

static_assert(std::is_trivially_copyable_v<MeasurementLocator>);
static_assert(std::is_trivially_copyable_v<Tracklet>);

} // namespace o2::itsmft::tracking

#endif // ALICEO2_ITSMFT_TRACKING_TRACKINGPRIMITIVES_H_

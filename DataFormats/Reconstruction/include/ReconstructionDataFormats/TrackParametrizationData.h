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

/// @file TrackParametrizationData.h
/// @brief Common storage for track parameterizations.

#ifndef INCLUDE_RECONSTRUCTIONDATAFORMATS_TRACKPARAMETRIZATIONDATA_H_
#define INCLUDE_RECONSTRUCTIONDATAFORMATS_TRACKPARAMETRIZATIONDATA_H_

#include "GPUCommonDef.h"
#include "GPUCommonRtypes.h"

#include <array>

namespace o2::track
{

template <typename value_T, int nParams>
class TrackParametrizationData
{
 public:
  using value_t = value_T;

  GPUdDefault() TrackParametrizationData() = default;
  GPUd() explicit TrackParametrizationData(value_t s) : mX{s} {}
  GPUdDefault() TrackParametrizationData(const TrackParametrizationData&) = default;
  GPUdDefault() TrackParametrizationData(TrackParametrizationData&&) = default;
  GPUhdDefault() TrackParametrizationData& operator=(const TrackParametrizationData&) = default;
  GPUhdDefault() TrackParametrizationData& operator=(TrackParametrizationData&&) = default;
  GPUdDefault() ~TrackParametrizationData() = default;

 protected:
  value_t mX = 0;          ///< Intrinsic coordinate of the track parameterization
  value_t mP[nParams] = {}; ///< Local track parameters

  ClassDefNV(TrackParametrizationData, 1);
};

template <typename value_T, int nCov>
class TrackCovarianceData
{
 public:
  using value_t = value_T;

  GPUdDefault() TrackCovarianceData() = default;
  GPUdDefault() TrackCovarianceData(const TrackCovarianceData&) = default;
  GPUdDefault() TrackCovarianceData(TrackCovarianceData&&) = default;
  GPUhdDefault() TrackCovarianceData& operator=(const TrackCovarianceData&) = default;
  GPUhdDefault() TrackCovarianceData& operator=(TrackCovarianceData&&) = default;
  GPUdDefault() ~TrackCovarianceData() = default;

 protected:
  std::array<value_t, nCov> mC{}; ///< Packed covariance matrix

  ClassDefNV(TrackCovarianceData, 1);
};

} // namespace o2::track

#endif

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
///
/// \file TrackingParamRef.h
/// \brief Resolve CA tracker configurable params per detector without duplicating ITS registration
///

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKINGPARAMREF_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKINGPARAMREF_H_

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/Constants.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITStracking/TrackingConfigParam.h"

namespace o2::itsmft::tracking
{

/// MFT uses o2::itsmft::TrackerParamConfig; ITS production params stay in O2::ITStracking.
template <o2::detectors::DetID::ID DetId>
struct TrackerParamRef;

template <>
struct TrackerParamRef<o2::detectors::DetID::MFT> {
  using Type = o2::itsmft::TrackerParamConfig<o2::detectors::DetID::MFT>;
  static const Type& get() { return Type::Instance(); }
  static constexpr int nLayers() { return Type::getNLayers(); }
};

template <>
struct TrackerParamRef<o2::detectors::DetID::ITS> {
  using Type = o2::its::TrackerParamConfig;
  static const Type& get() { return Type::Instance(); }
  static constexpr int nLayers() { return constants::ITSNLayers; }
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKINGPARAMREF_H_ */

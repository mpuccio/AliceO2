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
/// \file MFTForwardRefit.h
/// \brief MFT-specific forward Kalman refit for CA track seeds
///

#ifndef ALICEO2_ITSMFT_TRACKING_MFTFORWARDREFIT_H_
#define ALICEO2_ITSMFT_TRACKING_MFTFORWARDREFIT_H_

#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/Constants.h"
#include "ITSMFTTracking/MFTCATrack.h"
#include "ITSMFTTracking/TimeFrame.h"

namespace o2::itsmft::tracking
{

bool refitTrackFwd(const TrackSeed<constants::MFTNLayers>& seed,
                   MFTCATrack& track,
                   const TimeFrame<constants::MFTNLayers>& tf,
                   const TrackingParameters& params,
                   float bz);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_MFTFORWARDREFIT_H_ */

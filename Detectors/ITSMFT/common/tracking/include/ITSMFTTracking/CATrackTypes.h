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
/// \file CATrackTypes.h
/// \brief Per-detector CA track storage type
///

#ifndef ALICEO2_ITSMFT_TRACKING_CATRACKTYPES_H_
#define ALICEO2_ITSMFT_TRACKING_CATRACKTYPES_H_

#include "DataFormatsITS/TrackITS.h"
#include "ITSMFTTracking/Constants.h"

namespace o2::itsmft::tracking
{

class MFTCATrack;

template <int NLayers>
struct CATrackTypeHelper {
  using type = o2::its::TrackITSExt;
};

template <int NLayers>
using CATrackType = typename CATrackTypeHelper<NLayers>::type;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_CATRACKTYPES_H_ */

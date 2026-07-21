// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACESTATEOPERATIONRESULT_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACESTATEOPERATIONRESULT_H_

#include <cstdint>

namespace o2::itsmft::tracking
{

enum class OperationFailureReason : uint8_t {
  SourceFamilyMismatch = 0,
  NonFiniteInput = 1,
  NonFiniteOutput = 2,
  InvalidCovariance = 3,
  UnreachableTarget = 4,
  PropagationFailure = 5,
  MaterialFailure = 6,
  PredictedChi2Failure = 7,
  UpdateFailure = 8,
  RotationFailure = 9,
  AlphaMismatch = 10,
  ReferenceCoordinateMismatch = 11
};

static_assert(sizeof(OperationFailureReason) == sizeof(uint8_t));

} // namespace o2::itsmft::tracking

#endif // ALICEO2_ITSMFT_TRACKING_SURFACESTATEOPERATIONRESULT_H_

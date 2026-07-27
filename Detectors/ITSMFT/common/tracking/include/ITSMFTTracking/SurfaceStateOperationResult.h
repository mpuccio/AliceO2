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
  ReferenceCoordinateMismatch = 11,
  // forward::buildSeed's own established strict-boundary rejections
  // (the retained detail::mftFwdFitCellClusters / buildCellSeed<DiskDisk>
  // initializer treats an insufficient inner/outer z-ordering margin, or any
  // inner-middle/inner-outer separation below its 1e-6f minimum, as a hard
  // rejection before any direction estimate is computed -- not a NaN/Inf
  // artifact of the arithmetic itself). Distinct from NonFiniteInput (raw
  // inputs are finite; the *geometry* is degenerate) and from
  // NonFiniteOutput (no output was even attempted yet). barrel::buildSeed
  // never raises this: its retained oracle (o2::its::track::buildTrackSeed)
  // has no analogous early rejection, only NonFiniteOutput-detected numeric
  // fallbacks.
  SeedGeometryDegenerate = 12,
  // barrel::/forward::buildSeed(SeedAnchor, ...): the anchor argument is
  // not one of the enum's locked values (SeedAnchor::Inner/Outer).
  // Distinct from NonFiniteInput -- the raw measurement/bz/absCharge
  // inputs may be perfectly well-formed; it is the anchor selector itself
  // that is invalid -- so an out-of-range SeedAnchor is never
  // misclassified as a numeric-input problem.
  InvalidSeedAnchor = 13
};

static_assert(sizeof(OperationFailureReason) == sizeof(uint8_t));

} // namespace o2::itsmft::tracking

#endif // ALICEO2_ITSMFT_TRACKING_SURFACESTATEOPERATIONRESULT_H_

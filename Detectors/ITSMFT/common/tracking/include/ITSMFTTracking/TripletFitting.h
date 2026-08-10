// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRIPLETFITTING_H_
#define ALICEO2_ITSMFT_TRACKING_TRIPLETFITTING_H_

#include <array>

#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"

namespace o2::itsmft::tracking
{

struct SymmetricCovariance3D {
  double xx{0.};
  double xy{0.};
  double xz{0.};
  double yy{0.};
  double yz{0.};
  double zz{0.};
};

struct TripletFitObservation {
  std::array<double, 3> position{};
  SymmetricCovariance3D covariance{};
};

struct TripletFitProcessNoise {
  // Variance of the physical space-angle kick at the middle observation.
  double angularVariance{0.};
};

struct LocalTripletFitResult {
  double curvature{0.};
  double curvatureVariance{0.};
  double chi2{0.};
  double referenceTransverseCurvature{0.};
  double referenceSinTheta{0.};
  double phiTilde{0.};
  double thetaTilde{0.};
  double rhoPhi{0.};
  double rhoTheta{0.};
  double gammaThetaTheta{0.};
  double gammaThetaPhi{0.};
  double gammaPhiPhi{0.};
};

bool makeCylinderTripletFitObservation(const SurfaceDescriptor& surface,
                                       const SurfaceMeasurement& measurement,
                                       TripletFitObservation& observation) noexcept;

bool makeDiskTripletFitObservation(const SurfaceDescriptor& surface,
                                   const SurfaceMeasurement& measurement,
                                   TripletFitObservation& observation) noexcept;

bool makeTripletFitObservation(const SurfaceDescriptor& surface,
                               const SurfaceMeasurement& measurement,
                               TripletFitObservation& observation) noexcept;

bool fitLocalTripletUniformSolenoid(
  const std::array<TripletFitObservation, 3>& observations,
  const TripletFitProcessNoise& processNoise,
  LocalTripletFitResult& result) noexcept;

double fittedTripletTransverseMomentum(const LocalTripletFitResult& result,
                                       double bz) noexcept;

double fittedTripletTransverseMomentum(const LocalTripletFitResult& result,
                                       double bz, uint8_t absCharge) noexcept;

} // namespace o2::itsmft::tracking

#endif // ALICEO2_ITSMFT_TRACKING_TRIPLETFITTING_H_

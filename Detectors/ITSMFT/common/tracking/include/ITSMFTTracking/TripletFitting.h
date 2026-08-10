// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRIPLETFITTING_H_
#define ALICEO2_ITSMFT_TRACKING_TRIPLETFITTING_H_

#include <array>
#include <cstdint>

#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ReconstructionDataFormats/PID.h"

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

struct TripletFitMaterial {
  std::array<double, 3> unitNormal{};
  NominalSurfaceMaterial nominal{};
};

enum class MaterialTripletFitStatus : uint8_t {
  Success,
  InvalidInput,
  LocalFitFailure,
  MomentumUnresolved,
  MomentumBelowMaterialRange,
  ExcessiveScattering,
  StoppedInMaterial,
  MaterialEvaluationFailure,
  NoConvergence
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

struct MaterialAwareTripletFitResult {
  LocalTripletFitResult local{};
  double momentum{0.};
  double transverseMomentum{0.};
  double angularVariance{0.};
  double incidenceCosine{0.};
  uint8_t iterations{0};
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

bool makeCylinderTripletFitMaterial(const SurfaceDescriptor& surface,
                                    const SurfaceMeasurement& measurement,
                                    TripletFitMaterial& material) noexcept;

bool makeDiskTripletFitMaterial(const SurfaceDescriptor& surface,
                                const SurfaceMeasurement& measurement,
                                TripletFitMaterial& material) noexcept;

bool makeTripletFitMaterial(const SurfaceDescriptor& surface,
                            const SurfaceMeasurement& measurement,
                            TripletFitMaterial& material) noexcept;

bool fitLocalTripletUniformSolenoid(
  const std::array<TripletFitObservation, 3>& observations,
  const TripletFitProcessNoise& processNoise,
  LocalTripletFitResult& result) noexcept;

double fittedTripletTransverseMomentum(const LocalTripletFitResult& result,
                                       double bz) noexcept;

double fittedTripletTransverseMomentum(const LocalTripletFitResult& result,
                                       double bz, uint8_t absCharge) noexcept;

MaterialTripletFitStatus fitLocalTripletWithMaterial(
  const std::array<TripletFitObservation, 3>& observations,
  const TripletFitMaterial& middleMaterial,
  double bz, uint8_t absCharge, o2::track::PID pid,
  MaterialAwareTripletFitResult& result) noexcept;

} // namespace o2::itsmft::tracking

#endif // ALICEO2_ITSMFT_TRACKING_TRIPLETFITTING_H_

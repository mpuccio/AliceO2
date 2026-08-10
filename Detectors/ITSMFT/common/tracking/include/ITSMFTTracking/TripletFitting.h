// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRIPLETFITTING_H_
#define ALICEO2_ITSMFT_TRACKING_TRIPLETFITTING_H_

#include <array>
#include <type_traits>

#include "GPUCommonDef.h"
#include "GPUCommonMath.h"
#include "ITSMFTTracking/GlobalMeasurement.h"

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

struct TripletKinkVector {
  float theta{0.f};
  float phi{0.f};
};

// Two rows of the hit-coordinate Jacobian H for one of the three observations.
struct TripletHitJacobian {
  std::array<float, 3> theta{};
  std::array<float, 3> phi{};
};

// Linearized local-triplet factor from Eq. (19) of the General Triplet Track
// Fit. H is evaluated at kappaRef = -Psi_phi / rho_phi. Hit slot i is bound to
// CellSeed::getClusterReference(i); measurement and MS covariances stay outside
// this factor and are assembled when adjacent triplets are compared.
struct TripletFitFactor {
  TripletKinkVector psi{};
  TripletKinkVector rho{};
  std::array<TripletHitJacobian, 3> h{};

  GPUhdi() bool isValid() const noexcept
  {
    if (!o2::gpu::GPUCommonMath::Finite(psi.theta) ||
        !o2::gpu::GPUCommonMath::Finite(psi.phi) ||
        !o2::gpu::GPUCommonMath::Finite(rho.theta) ||
        !o2::gpu::GPUCommonMath::Finite(rho.phi) || rho.phi == 0.f) {
      return false;
    }
    for (const auto& hit : h) {
      for (int coordinate = 0; coordinate < 3; ++coordinate) {
        if (!o2::gpu::GPUCommonMath::Finite(hit.theta[coordinate]) ||
            !o2::gpu::GPUCommonMath::Finite(hit.phi[coordinate])) {
          return false;
        }
      }
    }
    return true;
  }
};

static_assert(std::is_standard_layout_v<TripletFitFactor>);
static_assert(std::is_trivially_copyable_v<TripletFitFactor>);
static_assert(sizeof(TripletFitFactor) == 88);

struct AdjacentTripletFitResult {
  double curvature{0.};
  double curvatureVariance{0.};
  double chi2{0.};
};

bool makeTripletFitObservation(const GlobalMeasurement& measurement,
                               TripletFitObservation& observation) noexcept;

bool makeTripletFitFactor(
  const std::array<TripletFitObservation, 3>& observations,
  TripletFitFactor& factor) noexcept;

// Closed-form minimization of Eq. (19) for adjacent triplets over one common
// curvature. observations are the four unique ordered hits; angularVariance
// contains the physical space-angle MS variance for each triplet.
bool fitAdjacentTripletFactors(
  const TripletFitFactor& firstFactor,
  const TripletFitFactor& secondFactor,
  const std::array<TripletFitObservation, 4>& observations,
  const std::array<double, 2>& angularVariance,
  AdjacentTripletFitResult& result) noexcept;

} // namespace o2::itsmft::tracking

#endif // ALICEO2_ITSMFT_TRACKING_TRIPLETFITTING_H_

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

#ifndef ALICEO2_ITSMFT_TRACKING_MATERIALPHYSICS_H_
#define ALICEO2_ITSMFT_TRACKING_MATERIALPHYSICS_H_

#include <cstdint>
#include "GPUCommonDef.h"
#include <cmath>
#include "ReconstructionDataFormats/TrackParametrization.h"
#include "ReconstructionDataFormats/TrackUtils.h"

#include "ReconstructionDataFormats/PID.h"

// Shared CPU/GPU material calculation.

namespace o2::itsmft::tracking::material
{

// Material traversal direction relative to the particle momentum, independent
// of any propagation or covariance sign convention in the caller.
enum class MaterialTraversalDirection : uint8_t {
  AlongMomentum = 0,
  OppositeMomentum = 1
};

// Unsigned, path-integrated material budget. Both fields are non-negative;
// direction is supplied separately.
struct IntegratedMaterialBudget {
  float xOverX0;             ///< thickness in units of radiation length
  float arealDensityGPerCm2; ///< crossed length*density, g/cm^2
};

// Scalar material-physics kernel for charged particles.
// pid supplies the mass; absCharge supplies |q| for energy-loss and
// scattering scale factors. absCharge must be nonzero and need not equal
// PID::getCharge().
//
// Validation precedence (first failure wins): invalid direction, negative
// material, non-positive momentum, invalid PID, then a charged massless PID.
// The PID range is checked before accessing its mass.
//
// For charged massive states, non-positive beta^2 is rejected before either
// material-effect calculation.
//
// For charged massive states, momentumGeV is the caller-selected physical
// momentum; no covariance projection is performed. Energy loss uses the same
// capped-substep Bethe-Bloch algorithm as
// o2::track::TrackParametrizationWithError::correctForMaterial(). The
// requested substep count is
//   1 + floor(|dE_full| / eKin * o2::track::ELoss2EKinThreshInv)
// with a range-bounded float-to-int conversion, capped at
// o2::track::MaxELossIter (50). All arealDensityGPerCm2 is processed; only the
// granularity changes. Bethe-Bloch is recomputed from the current momentum
// at each substep.
// MaterialTraversalDirection::AlongMomentum subtracts energy per substep;
// OppositeMomentum adds it. Reject particles stopped in material or ending
// with momentum below 0.01 GeV/c.
//
// highlandTheta2Rad2 and relativeInverseMomentumVariance use the simplified
// O2 Highland variance (no logarithmic correction) and pre-material momentum,
// energy, and beta. Both scale with absCharge^2. highlandTheta2Rad2 > pi^2
// is rejected.
//
// This kernel does not construct track states, detector geometry, or
// ITS/MFT/topology/propagation objects.
// Output values are committed only on success.
GPUdi() bool calculateMaterialPhysics(
  float momentumGeV,
  o2::track::PID pid,
  uint8_t absCharge,
  MaterialTraversalDirection direction,
  IntegratedMaterialBudget material,
  float& momentumAfterGeV,
  float& highlandTheta2Rad2,
  float& relativeInverseMomentumVariance) noexcept;



namespace materialphysics_detail
{
constexpr float kHighlandConst2 = 0.0136f * 0.0136f;
constexpr float kStragglingConst = 0.0007f;
constexpr float kMinMomentumGeV = 0.01f;

// Compute the capped substep count without an out-of-range float-to-int
// conversion.
GPUdi() uint8_t classifySubsteps(float fullStepEnergyLossGeV, float kineticEnergyGeV) noexcept
{
  const float ratio = std::fabs(fullStepEnergyLossGeV) / kineticEnergyGeV * o2::track::ELoss2EKinThreshInv;
  if (ratio >= static_cast<float>(o2::track::MaxELossIter)) {
    return static_cast<uint8_t>(o2::track::MaxELossIter);
  }
  // Keep the conversion in range even when ratio is unordered. Subsequent
  // arithmetic remains responsible for propagating invalid inputs.
  const float boundedRatio = ratio < static_cast<float>(o2::track::MaxELossIter) ? ratio : 0.f;
  const int requested = 1 + static_cast<int>(boundedRatio);
  return static_cast<uint8_t>(requested);
}

} // namespace materialphysics_detail

GPUdi() bool calculateMaterialPhysics(
  float momentumGeV,
  o2::track::PID pid,
  uint8_t absCharge,
  MaterialTraversalDirection direction,
  IntegratedMaterialBudget material,
  float& momentumAfterGeV,
  float& outHighlandTheta2Rad2,
  float& outRelativeInverseMomentumVariance) noexcept
{
  using namespace materialphysics_detail;
  if (direction != MaterialTraversalDirection::AlongMomentum && direction != MaterialTraversalDirection::OppositeMomentum) {
    return false;
  }
  if (material.xOverX0 < 0.f || material.arealDensityGPerCm2 < 0.f) {
    return false;
  }
  if (momentumGeV <= 0.f) {
    return false;
  }
  if (pid.getID() >= o2::track::PID::NIDsTot) {
    return false;
  }
  const float mass = pid.getMass();
  if (mass == 0.f) {
    return false;
  }

  const float q2 = static_cast<float>(absCharge) * static_cast<float>(absCharge);
  const float p0 = momentumGeV;
  const float p0Squared = p0 * p0;
  const float e0 = std::sqrt(p0Squared + mass * mass);
  const float beta2 = p0Squared / (e0 * e0);
  if (beta2 <= 0.f) {
    return false;
  }

  float e = e0;
  float p = p0;

  if (material.arealDensityGPerCm2 > 0.f) {
    const float ekin = e0 - mass;
    const float bg0 = p0 / mass;
    const float dedx0 = o2::track::BetheBlochSolidOpt<float>(bg0) * q2;
    const float fullStepEnergyLoss = dedx0 * material.arealDensityGPerCm2;

    const uint8_t substeps = materialphysics_detail::classifySubsteps(fullStepEnergyLoss, ekin);

    const float arealDensityStep = material.arealDensityGPerCm2 / static_cast<float>(substeps);
    for (uint8_t i = 0; i < substeps; ++i) {
      const float bg = p / mass;
      const float dedx = o2::track::BetheBlochSolidOpt<float>(bg) * q2;
      const float dE = dedx * arealDensityStep;
      e = (direction == MaterialTraversalDirection::AlongMomentum) ? (e - dE) : (e + dE);
      if (e <= mass) {
        return false;
      }
      p = std::sqrt(e * e - mass * mass);
    }
  }

  if (p < kMinMomentumGeV) {
    return false;
  }
  const float signedEnergyChangeGeV = e - e0;

  float highlandTheta2Rad2 = 0.f;
  if (material.xOverX0 > 0.f) {
    highlandTheta2Rad2 = kHighlandConst2 / (beta2 * p0 * p0) * material.xOverX0 * q2;
    if (highlandTheta2Rad2 > o2::constants::math::PI * o2::constants::math::PI) {
      return false;
    }
  }

  float relativeInverseMomentumVariance = 0.f;
  if (signedEnergyChangeGeV != 0.f) {
    relativeInverseMomentumVariance = kStragglingConst * kStragglingConst * std::fabs(signedEnergyChangeGeV) * e0 * e0 / (p0 * p0 * p0 * p0);
  }

  momentumAfterGeV = p;
  outHighlandTheta2Rad2 = highlandTheta2Rad2;
  outRelativeInverseMomentumVariance = relativeInverseMomentumVariance;
  return true;
}


} // namespace o2::itsmft::tracking::material

#endif // ALICEO2_ITSMFT_TRACKING_MATERIALPHYSICS_H_

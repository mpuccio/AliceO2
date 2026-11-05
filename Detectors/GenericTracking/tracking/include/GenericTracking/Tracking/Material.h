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

/// \file Material.h
/// \brief What a propagation step crossed, and the scattering it causes.

#ifndef GTRK_TRACKING_MATERIAL_H_
#define GTRK_TRACKING_MATERIAL_H_

#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Trig.h"

#include "GPUCommonMath.h"

namespace o2::gtrk
{

/// beta p [GeV] at momentum \a momentumGeV for the mass hypothesis
/// \a massGeV (`Propagator::massGeV`, positive): what Highland's formula
/// divides by.
GPUhdi() GPUAlwaysInline() Scalar betaMomentum(Scalar momentumGeV, Scalar massGeV)
{
  return momentumGeV * momentumGeV / o2::gpu::GPUCommonMath::Sqrt((momentumGeV * momentumGeV) + (massGeV * massGeV));
}

/// Highland's scattering angle for a step, in radians.
///
/// theta0 = (13.6 MeV / beta c p) z sqrt(x/X0) [1 + 0.038 ln(x/X0)], with z = 1
/// and \a betaMomentumGeV = beta p from the mass hypothesis (`betaMomentum`).
/// 13.6 MeV is 0.0136 in the GeV this state is parameterised in. O2's
/// `correctForMaterial` leaves out the logarithmic term.
///
/// The logarithmic correction turns negative below x/X0 = exp(-1/0.038), about
/// 4e-12, which is far outside anything a detector produces but is guarded
/// because a negative sigma is worse than a zero one.
GPUhdi() GPUAlwaysInline() Scalar highlandTheta0(Scalar x2X0, Scalar betaMomentumGeV)
{
  if (x2X0 <= 0.F || betaMomentumGeV <= 0.F) {
    return 0.F;
  }
  const Scalar corr = 1.F + 0.038F * logOf(x2X0);
  if (corr <= 0.F) {
    return 0.F;
  }
  return (0.0136F / betaMomentumGeV) * o2::gpu::GPUCommonMath::Sqrt(x2X0) * corr;
}

/// O2's mean Bethe-Bloch for silicon, `o2::track::BetheBlochSolidOpt`, at
/// \a betaGamma, in GeV / (g/cm^2): its constants and its arithmetic, for
/// host and device alike, where O2's is device-only in a GPU translation
/// unit. testTrackFitter.cxx's `TheEnergyLossIsO2s` holds it to O2's.
GPUhdi() GPUAlwaysInline() Scalar betheBlochSolid(Scalar betaGamma)
{
  constexpr Scalar kK = 0.307075e-3F; // [GeV cm^2 / g]
  // O2's 2.303, not ln 10: the same numbers as O2.
  constexpr Scalar kP1 = 0.20F * 2.303F; // NOLINT(modernize-use-std-numbers)
  constexpr Scalar kP2 = 3.00F * 2.303F; // NOLINT(modernize-use-std-numbers)
  constexpr Scalar kMeanZA = 0.49848F;
  constexpr Scalar kLhwI = -1.7175226F;         // log(28.816e-9 sqrt(rho Z/A) / I)
  constexpr Scalar kLog2MeToMeanI = 8.6839805F; // log(2 m_e / I)
  const Scalar bg2 = betaGamma * betaGamma;
  const Scalar beta2 = bg2 / (1.F + bg2);
  Scalar d2 = 0.F;
  const Scalar x = o2::gpu::GPUCommonMath::Log(betaGamma);
  if (x > kP2) {
    d2 = kLhwI - 0.5F + x;
  } else if (x > kP1) {
    const Scalar r = (kP2 - x) / (kP2 - kP1);
    d2 = kLhwI - 0.5F + x + ((0.5F - kLhwI - kP1) * r * r * r);
  }
  const Scalar dedx = kK * kMeanZA / beta2 * (kLog2MeToMeanI + x + x - beta2 - d2);
  return dedx > 0.F ? dedx : 0.F;
}

/// One step's energy loss, as O2 corrects for it
/// (`TrackParametrizationWithError::correctForMaterial`, with a linearisation
/// reference).
struct EnergyLoss {
  Scalar momentumScale = 1.F; ///< p before over p after: the factor on q/p
  Scalar varQOverP = 0.F;     ///< the straggling, added to q/p's variance
  bool ok = true;             ///< false: stopped, or below 10 MeV after it
};

/// The mean energy lost -- or, against the momentum, regained -- over
/// \a xRho [g/cm^2] by a particle of mass \a massGeV at \a qOverP: O2's
/// Bethe-Bloch for silicon (`betheBlochSolid`) [GeV / (g/cm^2)], in sub-steps of at most 2.5 % of the kinetic energy, and
/// O2's straggling, 0.0007 sqrt(|dE|) [GeV] on the energy. \a xRho is negative
/// along the momentum. A MatLUT carries no material's I or Z/A, which is why
/// silicon's stand in for all of it, as in O2.
GPUhdi() EnergyLoss energyLoss(Scalar qOverP, Scalar massGeV, Scalar xRho)
{
  EnergyLoss result;
  if (!(o2::gpu::GPUCommonMath::Abs(xRho) > 0.F)) {
    return result;
  }
  constexpr Scalar kMinMomentum = 0.01F;          // O2's kMinP
  constexpr Scalar kStepsPerKineticEnergy = 40.F; // O2's ELoss2EKinThreshInv
  constexpr int kMaxSteps = 50;                   // O2's MaxELossIter
  const Scalar mass2 = massGeV * massGeV;
  const Scalar massInv = 1.F / massGeV;
  const Scalar p0 = 1.F / o2::gpu::GPUCommonMath::Abs(qOverP);
  const Scalar e0 = o2::gpu::GPUCommonMath::Sqrt((p0 * p0) + mass2);
  Scalar p = p0;
  Scalar e = e0;
  Scalar dE = betheBlochSolid(p * massInv) * xRho;
  int steps = 1 + static_cast<int>(o2::gpu::GPUCommonMath::Abs(dE) / (e - massGeV) * kStepsPerKineticEnergy);
  steps = steps > kMaxSteps ? kMaxSteps : steps;
  if (steps > 1) {
    dE /= static_cast<Scalar>(steps);
    xRho /= static_cast<Scalar>(steps);
  }
  while (steps-- > 0) {
    e += dE;
    if (!(e > massGeV)) {
      result.ok = false;
      return result;
    }
    p = o2::gpu::GPUCommonMath::Sqrt((e * e) - mass2);
    if (steps > 0) {
      dE = betheBlochSolid(p * massInv) * xRho;
    }
  }
  if (p < kMinMomentum) {
    result.ok = false;
    return result;
  }
  result.momentumScale = p0 / p;
  const Scalar sigmaQOverP = 0.0007F * o2::gpu::GPUCommonMath::Sqrt(o2::gpu::GPUCommonMath::Abs(e - e0)) * e0 / (p0 * p0) * o2::gpu::GPUCommonMath::Abs(qOverP);
  result.varQOverP = sigmaQOverP * sigmaQOverP;
  return result;
}

/// The variances a step's scattering adds to phi and theta, for a covariance
/// kept as factors (StepNoise): zero for no material.
struct ScatteringNoise {
  Scalar phi = 0.F;
  Scalar theta = 0.F;
};

GPUhdi() GPUAlwaysInline() ScatteringNoise scatteringNoise(Scalar theta0, Scalar theta)
{
  if (theta0 <= 0.F) {
    return {};
  }
  const Scalar t2 = theta0 * theta0;
  const auto sinTheta = sinOf(theta);
  return {.phi = t2 / (sinTheta * sinTheta), .theta = t2};
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_MATERIAL_H_

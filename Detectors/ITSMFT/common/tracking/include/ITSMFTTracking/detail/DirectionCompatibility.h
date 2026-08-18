// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License Version 3, copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_DIRECTIONCOMPATIBILITY_H_
#define ALICEO2_ITSMFT_TRACKING_DIRECTIONCOMPATIBILITY_H_

#include <array>

#ifndef GPUCA_GPUCODE
#include "ITSMFTTracking/GlobalMeasurement.h"
#endif

namespace o2::itsmft::tracking
{

#ifndef GPUCA_GPUCODE

struct DirectionObservation {
  double r{0.};
  double z{0.};
  double varianceR{0.};
  double covarianceRZ{0.};
  double varianceZ{0.};
  double radialUnitX{0.};
  double radialUnitY{0.};
};

struct TransverseDirectionObservation {
  double x{0.};
  double y{0.};
  double varianceX{0.};
  double covarianceXY{0.};
  double varianceY{0.};
};

struct DirectionProcessNoise {
  // Variance of a thin angular kick at the middle observation.
  double angularVariance{0.};
};

struct CellDirectionCompatibility {
  double residual{0.};
  double variance{0.};
  double chi2{0.};
};

struct TransverseDirectionCompatibility {
  double deltaPhi{0.};
  double maximumBending{0.};
  double variance{0.};
  double chi2{0.};
};

bool makeTransverseDirectionObservation(const GlobalMeasurement& measurement,
                                        TransverseDirectionObservation& observation) noexcept;

bool trackletDirectionsAreTransverselyCompatible(
  const std::array<TransverseDirectionObservation, 3>& observations,
  float firstPhi, float secondPhi,
  const DirectionProcessNoise& processNoise,
  float bz, float trackletMinPt, float nSigmaCut,
  TransverseDirectionCompatibility& compatibility) noexcept;

bool makeDirectionObservation(const GlobalMeasurement& measurement,
                              DirectionObservation& observation) noexcept;

bool cellDirectionsAreCompatible(const std::array<DirectionObservation, 3>& observations,
                                 const DirectionProcessNoise& processNoise,
                                 float beamPositionVariance,
                                 float nSigmaCut,
                                 CellDirectionCompatibility& compatibility) noexcept;

inline bool cellDirectionsAreCompatible(const std::array<DirectionObservation, 3>& observations,
                                        const DirectionProcessNoise& processNoise,
                                        float nSigmaCut,
                                        CellDirectionCompatibility& compatibility) noexcept
{
  return cellDirectionsAreCompatible(observations, processNoise, 0.f, nSigmaCut, compatibility);
}

#endif

} // namespace o2::itsmft::tracking

#endif

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

/// \file UniformMaterial.h
/// \brief The same material everywhere: the tests' material model.
///
/// The fixture the Jacobian and the covariance transport are validated
/// against: finite differences and a Monte Carlo perturbation both need a
/// budget that can be computed by hand, which a binary LUT cannot give.

#ifndef GTRK_TEST_UNIFORMMATERIAL_H_
#define GTRK_TEST_UNIFORMMATERIAL_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"

#include "DetectorsBase/MatCell.h"
#include "GPUCommonMath.h"

#include <array>
#include <cstddef>

namespace o2::gtrk
{

struct UniformMaterial {
  Scalar x2X0PerCm = 0.F; ///< radiation lengths per cm
  Scalar rho = 0.F;       ///< density [g/cm^3]

  /// GPUdi like MatLayerCylSet::getMatBudget, whose interface this mirrors:
  /// MatBudget's constructor is device-only in a GPU translation unit.
  [[nodiscard]] GPUdi() o2::base::MatBudget getMatBudget(float x0, float y0, float z0, float x1, float y1, float z1) const
  {
    const float dx = x1 - x0;
    const float dy = y1 - y0;
    const float dz = z1 - z0;
    const float len = o2::gpu::GPUCommonMath::Sqrt(dx * dx + dy * dy + dz * dz);
    o2::base::MatBudget budget;
    budget.meanX2X0 = x2X0PerCm * len;
    budget.meanRho = rho;
    budget.length = len;
    return budget;
  }
};

/// The same `x/X0` for every layer, as `layerMaterial` gives a `LayerBudget`
/// from a LUT.
template <DetectorLayout Det>
[[nodiscard]] inline std::array<double, static_cast<std::size_t>(Detector<Det>::nLayers)> uniformLayerMaterial(double x2X0)
{
  std::array<double, static_cast<std::size_t>(Detector<Det>::nLayers)> budget;
  budget.fill(x2X0);
  return budget;
}

} // namespace o2::gtrk

#endif // GTRK_TEST_UNIFORMMATERIAL_H_

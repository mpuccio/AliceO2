// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License Version 3, copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_CELLFINDING_H_
#define ALICEO2_ITSMFT_TRACKING_CELLFINDING_H_

#include <array>
#include <limits>

#ifndef GPUCA_GPUCODE
#include <cmath>

#include <gsl/span>

#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/SurfaceCatalogView.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceStateOperationResult.h"
#include "ITSMFTTracking/detail/TrackingKernelParameters.h"
#endif

namespace o2::itsmft::tracking
{

#ifndef GPUCA_GPUCODE

inline bool isRecognizedMatCorrType(o2::base::PropagatorF::MatCorrType corrType) noexcept
{
  return corrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE ||
         corrType == o2::base::PropagatorF::MatCorrType::USEMatCorrTGeo ||
         corrType == o2::base::PropagatorF::MatCorrType::USEMatCorrLUT;
}

struct AttachHitConfigView {
  gsl::span<const NominalSurfaceMaterial> layerMaterial;
  o2::base::PropagatorF::MatCorrType corrType{o2::base::PropagatorF::MatCorrType::USEMatCorrNONE};

  bool isValid(size_t expectedLayers) const noexcept
  {
    if (layerMaterial.size() < expectedLayers || !isRecognizedMatCorrType(corrType)) {
      return false;
    }
    for (size_t layer = 0; layer < expectedLayers; ++layer) {
      const auto& material = layerMaterial[layer];
      if (!isFiniteParam(material.xOverX0) || material.xOverX0 < 0.f ||
          !isFiniteParam(material.arealDensityGPerCm2) || material.arealDensityGPerCm2 < 0.f) {
        return false;
      }
    }
    return true;
  }
};

inline AttachHitConfigView bindAttachHitConfig(gsl::span<const NominalSurfaceMaterial> layerMaterial,
                                               const TrackingParameters& params) noexcept
{
  return {layerMaterial, params.CorrType};
}

enum class MaterialCorrectionModeSupport : uint8_t {
  Supported,
  Unsupported,
  InvalidMode,
  InvalidSurfaceKind
};

inline MaterialCorrectionModeSupport materialCorrectionModeSupport(SurfaceKind kind,
                                                                   o2::base::PropagatorF::MatCorrType corrType) noexcept
{
  if (!isRecognizedMatCorrType(corrType)) {
    return MaterialCorrectionModeSupport::InvalidMode;
  }
  if (kind != SurfaceKind::Cylinder && kind != SurfaceKind::Disk) {
    return MaterialCorrectionModeSupport::InvalidSurfaceKind;
  }
  if (kind == SurfaceKind::Cylinder && corrType != o2::base::PropagatorF::MatCorrType::USEMatCorrNONE) {
    return MaterialCorrectionModeSupport::Unsupported;
  }
  return MaterialCorrectionModeSupport::Supported;
}

struct LayerGeometryConfigView {
  gsl::span<const float> layerRadii;
  gsl::span<const NominalSurfaceMaterial> layerMaterial;

  bool isValid(size_t expectedLayers) const noexcept
  {
    return layerRadii.size() >= expectedLayers && layerMaterial.size() >= expectedLayers;
  }
};

inline LayerGeometryConfigView bindLayerGeometryConfig(const TrackingParameters& params,
                                                       const AttachHitConfigView& attachHitConfig) noexcept
{
  return {gsl::span<const float>{params.LayerRadii.data(), params.LayerRadii.size()}, attachHitConfig.layerMaterial};
}

bool buildCylinderCellSeed(const SurfaceMeasurement& measurementInner,
                           const SurfaceMeasurement& measurementMiddle,
                           const SurfaceMeasurement& measurementOuter,
                           const std::array<NominalSurfaceMaterial, 3>& material,
                           float bz, uint8_t absCharge, o2::track::PID pid,
                           SurfaceKinematicState& outState, float& chi2,
                           const TrackingKernelParameters& params,
                           OperationFailureReason& reason) noexcept;

bool buildDiskCellSeed(const SurfaceMeasurement& measurementInner,
                       const SurfaceMeasurement& measurementMiddle,
                       const SurfaceMeasurement& measurementOuter,
                       const std::array<NominalSurfaceMaterial, 3>& material,
                       float bz, uint8_t absCharge, o2::track::PID pid,
                       SurfaceKinematicState& outState, float& chi2,
                       const TrackingKernelParameters& params,
                       OperationFailureReason& reason) noexcept;

bool attachCylinderHit(SurfaceKinematicState& state, const SurfaceMeasurement& measurement,
                       const NominalSurfaceMaterial& material, float bz, float& chi2,
                       const TrackingKernelParameters& params,
                       OperationFailureReason& reason) noexcept;

bool attachDiskHit(SurfaceKinematicState& state, const SurfaceMeasurement& measurement,
                   const NominalSurfaceMaterial& material, float bz, float& chi2,
                   const TrackingKernelParameters& params,
                   OperationFailureReason& reason) noexcept;

bool cellsCylinderAreCompatible(const SurfaceKinematicState& current,
                                const SurfaceKinematicState& next,
                                int currentSecondClusterIndex, int nextFirstClusterIndex,
                                float bz, const TrackingKernelParameters& params) noexcept;

bool cellsDiskAreCompatible(const SurfaceKinematicState& current,
                            const SurfaceKinematicState& next,
                            int currentSecondClusterIndex, int nextFirstClusterIndex,
                            float bz, const TrackingKernelParameters& params) noexcept;

inline bool passesCylinderCellRoadPrecut(const GlobalPoint3F&, const GlobalPoint3F&, const GlobalPoint3F&,
                                         int, int, int, gsl::span<const float>,
                                         const TrackingKernelParameters&) noexcept
{
  return true;
}

inline bool passesDiskCellRoadPrecut(const GlobalPoint3F& pointInner, const GlobalPoint3F& pointMiddle,
                                     const GlobalPoint3F& pointOuter, int layerInner, int layerMiddle,
                                     int layerOuter, gsl::span<const float> perLayerReferenceZ,
                                     const TrackingKernelParameters& params) noexcept
{
  const auto distanceToSeedLineSquared = [](const GlobalPoint3F& from, const GlobalPoint3F& to,
                                            const GlobalPoint3F& point) -> float {
    const float dxSeed = to.x - from.x;
    const float dySeed = to.y - from.y;
    const float dzSeed = to.z - from.z;
    if (std::abs(dzSeed) < 1.e-9f) {
      return std::numeric_limits<float>::max();
    }
    const float invdzSeed = (point.z - from.z) / dzSeed;
    const float dx = point.x - (from.x + dxSeed * invdzSeed);
    const float dy = point.y - (from.y + dySeed * invdzSeed);
    return dx * dx + dy * dy;
  };
  const auto conicalRoadR2Scale = [](float zFrom, float zTo) -> float {
    if (std::abs(zFrom) < 1.e-6f) {
      return 1.f;
    }
    const float dCone = 1.f + (zTo - zFrom) / zFrom;
    return dCone * dCone;
  };
  const float zInner = perLayerReferenceZ[layerInner];
  const float zMiddle = perLayerReferenceZ[layerMiddle];
  const float zOuter = perLayerReferenceZ[layerOuter];
  const float r2Cut = params.cellRoadRCut * params.cellRoadRCut;
  return distanceToSeedLineSquared(pointInner, pointOuter, pointMiddle) < r2Cut * conicalRoadR2Scale(zInner, zMiddle) &&
         distanceToSeedLineSquared(pointInner, pointMiddle, pointOuter) < r2Cut * conicalRoadR2Scale(zInner, zOuter) &&
         distanceToSeedLineSquared(pointMiddle, pointOuter, pointInner) < r2Cut * conicalRoadR2Scale(zMiddle, zInner);
}

#endif

} // namespace o2::itsmft::tracking

#endif

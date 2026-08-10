// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License Version 3, copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_CANDIDATEFINDING_H_
#define ALICEO2_ITSMFT_TRACKING_CANDIDATEFINDING_H_

#include <array>
#include <variant>
#ifndef GPUCA_GPUCODE
#include <gsl/span>

#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceStateOperationResult.h"
#include "ITSMFTTracking/detail/TrackingKernelParameters.h"
#endif

#ifndef GPUCA_GPUCODE
namespace o2::dataformats
{
template <typename Stamp>
class Vertex;
}
namespace o2::its
{
class TimeEstBC;
using Vertex = o2::dataformats::Vertex<TimeEstBC>;
struct Cluster;
} // namespace o2::its
#endif

namespace o2::itsmft::tracking
{

#ifndef GPUCA_GPUCODE

struct CylinderTrackletProjectionState {
  int fromLayer;
  int toLayer;
  float meanDeltaR;
  float targetMinR;
  float targetMaxR;
  float sourcePositionResolution;
  float transitionMSAngle;
  float transitionPhiCut;
};

struct DiskTrackletProjectionState {
  int fromLayer;
  int toLayer;
  float fromZ;
  float toZ;
  float meanDeltaZ;
  float sourceReferenceRadius;
  float transitionMSAngle;
  float transitionBendingAngle;
};

struct CylinderTrackletSearchWindow {
  int4 bins;
  float tanLambda;
  float sigmaZ;
  float phiCut;
  float nSigmaCut;

  bool acceptCandidate(const SurfaceMeasurement& sourceMeasurement,
                       const o2::its::Cluster& sourceLocator,
                       const SurfaceMeasurement& targetMeasurement,
                       const o2::its::Cluster& targetLocator,
                       float& tanLambdaOut) const;
};

struct DiskTrackletSearchWindow {
  int4 bins;
  float xProj;
  float yProj;
  float sigmaX;
  float sigmaY;
  float nSigmaCut;

  bool acceptCandidate(const SurfaceMeasurement& sourceMeasurement,
                       const SurfaceMeasurement& targetMeasurement,
                       float& tanLambdaOut) const;
};

using TrackletProjectionState = std::variant<CylinderTrackletProjectionState, DiskTrackletProjectionState>;
using TrackletSearchWindow = std::variant<CylinderTrackletSearchWindow, DiskTrackletSearchWindow>;

bool bindTrackletProjectionState(SurfaceKind kind, int fromLayer, int toLayer,
                                 gsl::span<const float> layerRadii,
                                 gsl::span<const float> diskReferenceZ,
                                 float targetMinR, float targetMaxR,
                                 float sourcePositionResolution,
                                 float transitionMSAngle, float transitionBendingAngle,
                                 TrackletProjectionState& out) noexcept;

bool projectTrackletSearchWindow(const SurfaceMeasurement& sourceMeasurement,
                                 const o2::its::Cluster& sourceLocator,
                                 const o2::its::Vertex& vertex,
                                 const TrackletProjectionState& transitionState,
                                 float bz, const o2::itsmft::IndexTableUtilsCore& indexUtils,
                                 const TrackingKernelParameters& params,
                                 TrackletSearchWindow& out);

int4 trackletSearchBins(const TrackletSearchWindow& window) noexcept;
int trackletSearchRowCount(const TrackletSearchWindow& window,
                           const o2::itsmft::IndexTableUtilsCore& indexUtils) noexcept;
int trackletSearchRowBin(const TrackletSearchWindow& window, int offset,
                         const o2::itsmft::IndexTableUtilsCore& indexUtils) noexcept;

bool acceptTrackletCandidate(const TrackletSearchWindow& window,
                             const SurfaceMeasurement& sourceMeasurement,
                             const o2::its::Cluster& sourceLocator,
                             const SurfaceMeasurement& targetMeasurement,
                             const o2::its::Cluster& targetLocator,
                             float& tanLambdaOut) noexcept;

bool projectCylinderSearchWindow(const SurfaceMeasurement& sourceMeasurement,
                                 const o2::its::Cluster& sourceLocator,
                                 const o2::its::Vertex& vertex,
                                 const CylinderTrackletProjectionState& transitionState,
                                 float bz, const o2::itsmft::IndexTableUtilsCore& indexUtils,
                                 const TrackingKernelParameters& params,
                                 CylinderTrackletSearchWindow& out);

bool projectDiskSearchWindow(const SurfaceMeasurement& sourceMeasurement,
                             const o2::its::Cluster& sourceLocator,
                             const o2::its::Vertex& vertex,
                             const DiskTrackletProjectionState& transitionState,
                             float bz, const o2::itsmft::IndexTableUtilsCore& indexUtils,
                             const TrackingKernelParameters& params,
                             DiskTrackletSearchWindow& out);

struct CylinderLayerScatteringInputs {
  float layerxX0;
};
struct DiskLayerScatteringInputs {
  float layerxX0;
  float layerRadius;
  float referenceCoordinate;
};

float cylinderLayerMultipleScatteringAngle(const CylinderLayerScatteringInputs& inputs, float trackletMinPt);
float diskLayerMultipleScatteringAngle(const DiskLayerScatteringInputs& inputs, float trackletMinPt);

struct DiskReferenceCoordinateView {
  gsl::span<const float> perLayerReferenceZ;
  bool isValid(size_t expectedLayers) const noexcept { return perLayerReferenceZ.size() >= expectedLayers; }
};

DiskReferenceCoordinateView bindLegacyMFTReferenceCoordinates() noexcept;

float clampCylinderTransitionCurvature(float oneOverR, float r2) noexcept;
float clampDiskTransitionCurvature(float oneOverR, float r2) noexcept;

struct TransitionScatteringBendingPrep {
  float msAngle;
  float phiCut;
};

TransitionScatteringBendingPrep prepareTransitionScatteringAndBending(
  gsl::span<const float> perLayerMSAngle, int fromLayer, int toLayer,
  float r1, float r2, float clampedOneOverR, float res1, float res2) noexcept;

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

struct DirectionObservation {
  double r{0.};
  double z{0.};
  double varianceR{0.};
  double covarianceRZ{0.};
  double varianceZ{0.};
};

struct DirectionProcessNoise {
  // Variance of an equivalent thin angular kick at the middle observation.
  double angularVariance{0.};
};

struct CellDirectionCompatibility {
  double residual{0.};
  double variance{0.};
  double chi2{0.};
};

bool makeCylinderDirectionObservation(const SurfaceDescriptor& surface,
                                      const SurfaceMeasurement& measurement,
                                      DirectionObservation& observation) noexcept;

bool makeDiskDirectionObservation(const SurfaceDescriptor& surface,
                                  const SurfaceMeasurement& measurement,
                                  DirectionObservation& observation) noexcept;

bool makeDirectionObservation(const SurfaceDescriptor& surface,
                              const SurfaceMeasurement& measurement,
                              DirectionObservation& observation) noexcept;

bool cellDirectionsAreCompatible(const std::array<DirectionObservation, 3>& observations,
                                 const DirectionProcessNoise& processNoise,
                                 float nSigmaCut,
                                 CellDirectionCompatibility& compatibility) noexcept;

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

bool buildCellSeed(SurfaceKind kind,
                   const SurfaceMeasurement& measurementInner,
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

#endif

} // namespace o2::itsmft::tracking

#endif

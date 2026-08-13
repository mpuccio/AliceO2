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

#include "ITSMFTTracking/GlobalMeasurement.h"
#include "ITSMFTTracking/IndexTableUtils.h"
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
  int toLayer;
  float meanDeltaR;
  float targetMinR;
  float targetMaxR;
  float sourcePositionResolution;
  float transitionMSAngle;
  float transitionPhiCut;
};

struct DiskTrackletProjectionState {
  int toLayer;
  float fromZ;
  float toZ;
  float meanDeltaZ;
  float sourceReferenceRadius;
  float transitionMSAngle;
  float transitionPhiCut;
};

struct CylinderTrackletSearchWindow {
  int4 bins;
  float tanLambda;
  float sigmaZ;
  float phiCut;
  float nSigmaCut;

  bool acceptCandidate(const GlobalMeasurement& sourceMeasurement,
                       const o2::its::Cluster& sourceLocator,
                       const GlobalMeasurement& targetMeasurement,
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

  bool acceptCandidate(const GlobalMeasurement& sourceMeasurement,
                       const GlobalMeasurement& targetMeasurement,
                       float& tanLambdaOut) const;
};

using TrackletProjectionState = std::variant<CylinderTrackletProjectionState, DiskTrackletProjectionState>;
using TrackletSearchWindow = std::variant<CylinderTrackletSearchWindow, DiskTrackletSearchWindow>;

bool bindTrackletProjectionState(SurfaceKind kind, int fromLayer, int toLayer,
                                 gsl::span<const float> layerRadii,
                                 gsl::span<const float> diskReferenceZ,
                                 float targetMinR, float targetMaxR,
                                 float sourcePositionResolution,
  float transitionMSAngle, float transitionPhiCut,
                                 TrackletProjectionState& out) noexcept;

bool projectTrackletSearchWindow(const GlobalMeasurement& sourceMeasurement,
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
                             const GlobalMeasurement& sourceMeasurement,
                             const o2::its::Cluster& sourceLocator,
                             const GlobalMeasurement& targetMeasurement,
                             const o2::its::Cluster& targetLocator,
                             float& tanLambdaOut) noexcept;

bool projectCylinderSearchWindow(const GlobalMeasurement& sourceMeasurement,
                                 const o2::its::Cluster& sourceLocator,
                                 const o2::its::Vertex& vertex,
                                 const CylinderTrackletProjectionState& transitionState,
                                 float bz, const o2::itsmft::IndexTableUtilsCore& indexUtils,
                                 const TrackingKernelParameters& params,
                                 CylinderTrackletSearchWindow& out);

bool projectDiskSearchWindow(const GlobalMeasurement& sourceMeasurement,
                             const o2::its::Cluster& sourceLocator,
                             const o2::its::Vertex& vertex,
                             const DiskTrackletProjectionState& transitionState,
                             float bz, const o2::itsmft::IndexTableUtilsCore& indexUtils,
                             const TrackingKernelParameters& params,
                             DiskTrackletSearchWindow& out);

bool buildCylinderCellSeed(const GlobalMeasurement& globalInner,
                           const GlobalMeasurement& globalMiddle,
                           const SurfaceMeasurement& measurementInner,
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
                   const GlobalMeasurement& globalInner,
                   const GlobalMeasurement& globalMiddle,
                   const GlobalMeasurement& globalOuter,
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

#endif

} // namespace o2::itsmft::tracking

#endif

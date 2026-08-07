// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License Version 3, copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKLETFINDING_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKLETFINDING_H_

#include <array>

#include "ITSMFTTracking/detail/TrackingKernelParameters.h"

#ifndef GPUCA_GPUCODE
#include <gsl/span>

#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"

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
  float meanDeltaZ;
  float nSigmaCut;

  bool acceptCandidate(const SurfaceMeasurement& sourceMeasurement,
                       const SurfaceMeasurement& targetMeasurement,
                       float& tanLambdaOut) const;
};

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

#endif

} // namespace o2::itsmft::tracking

#endif

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License Version 3, copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_CANDIDATEFINDING_H_
#define ALICEO2_ITSMFT_TRACKING_CANDIDATEFINDING_H_

#include <array>
#ifndef GPUCA_GPUCODE
#include <gsl/span>

#include "ITSMFTTracking/GlobalMeasurement.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/TrackingPrimitives.h"
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
} // namespace o2::its
#endif

namespace o2::itsmft::tracking
{

#ifndef GPUCA_GPUCODE

struct TrackletProjectionCache {
  int fromLayer;
  int toLayer;
  float fromRadius;
  float toRadius;
  float targetMinR;
  float targetMaxR;
  float targetMinZ;
  float targetMaxZ;
  float sourcePositionResolution;
  float edgeMSAngle;
  float edgePhiCut;
};

struct TrackletSearchWindow {
  int4 bins;
  float prediction[2];
  float variance[3];
};

bool bindTrackletProjectionCache(int fromLayer, int toLayer,
                                 gsl::span<const float> layerRadii,
                                 float targetMinR, float targetMaxR, float targetMinZ, float targetMaxZ,
                                 float sourcePositionResolution,
                                 float edgeMSAngle, float edgePhiCut,
                                 TrackletProjectionCache& out) noexcept;

bool projectTrackletSearchWindow(const GlobalMeasurement& sourceMeasurement,
                                 const o2::its::Vertex& vertex,
                                 float beamX, float beamY,
                                 float beamPositionVariance,
                                 SurfaceKind kind,
                                 const TrackletProjectionCache& edgeCache,
                                 float bz, const o2::itsmft::IndexTableUtilsCore& indexUtils,
                                 const TrackingKernelParameters& params,
                                 TrackletSearchWindow& out);

inline bool projectTrackletSearchWindow(const GlobalMeasurement& sourceMeasurement,
                                        const o2::its::Vertex& vertex,
                                        float beamX, float beamY,
                                        SurfaceKind kind,
                                        const TrackletProjectionCache& edgeCache,
                                        float bz, const o2::itsmft::IndexTableUtilsCore& indexUtils,
                                        const TrackingKernelParameters& params,
                                        TrackletSearchWindow& out)
{
  return projectTrackletSearchWindow(sourceMeasurement, vertex, beamX, beamY, 0.f, kind,
                                     edgeCache, bz, indexUtils, params, out);
}

bool acceptTrackletCandidate(const TrackletSearchWindow& window,
                             const GlobalMeasurement& sourceMeasurement,
                             const GlobalMeasurement& targetMeasurement,
                             SurfaceKind kind, float nSigmaCut,
                             float& tanLambdaOut) noexcept;

#endif

} // namespace o2::itsmft::tracking

#endif

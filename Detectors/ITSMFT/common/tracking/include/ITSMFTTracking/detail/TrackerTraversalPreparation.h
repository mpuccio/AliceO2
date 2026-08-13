// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License Version 3, copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKERTRAVERSALPREPARATION_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKERTRAVERSALPREPARATION_H_

#ifndef GPUCA_GPUCODE
#include <gsl/span>

#endif

namespace o2::itsmft::tracking
{

#ifndef GPUCA_GPUCODE

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

float clampTransitionCurvature(float oneOverR, float outerRadius) noexcept;

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

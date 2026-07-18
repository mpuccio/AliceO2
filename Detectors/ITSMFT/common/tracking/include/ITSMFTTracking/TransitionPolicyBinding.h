// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYBINDING_H_
#define ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYBINDING_H_

#include "ITSMFTTracking/TransitionPolicyState.h"

// Host-only: TrackingParameters (Configuration.h) owns std::vector members
// and is not itself device-compatible; only the bound Params struct crosses
// into device-facing policy operations. This binding step has no device
// counterpart and must never be compiled for device code.
#ifndef GPUCA_GPUCODE

#include "ITSMFTTracking/Configuration.h"
#include <gsl/span>

namespace o2::itsmft::tracking
{

/// Host view of the existing per-surface material and propagation-correction
/// configuration consumed by attachHit<Tag>. The view borrows one iteration's
/// TrackingParameters storage and is bound/validated with the typed family
/// Params before traversal starts.
struct AttachHitPolicyConfigView {
  gsl::span<const float> layerxX0;
  o2::base::PropagatorF::MatCorrType corrType{o2::base::PropagatorF::MatCorrType::USEMatCorrNONE};

  bool isValid(size_t expectedLayers) const noexcept
  {
    if (layerxX0.size() < expectedLayers ||
        (corrType != o2::base::PropagatorF::MatCorrType::USEMatCorrNONE &&
         corrType != o2::base::PropagatorF::MatCorrType::USEMatCorrTGeo &&
         corrType != o2::base::PropagatorF::MatCorrType::USEMatCorrLUT)) {
      return false;
    }
    for (size_t layer = 0; layer < expectedLayers; ++layer) {
      if (!isFiniteParam(layerxX0[layer]) || layerxX0[layer] < 0.f) {
        return false;
      }
    }
    return true;
  }
};

inline AttachHitPolicyConfigView bindAttachHitPolicyConfig(const TrackingParameters& params) noexcept
{
  return {gsl::span<const float>{params.LayerxX0.data(), params.LayerxX0.size()}, params.CorrType};
}

/// Binds one iteration's legacy TrackingParameters into a typed,
/// bounds-checkable policy parameter block (TransitionPolicyState.h). Callers
/// must invoke this once per iteration, outside any candidate/neighbour/road
/// loop (D007 / Architecture.md 10.1) -- never per-candidate.
template <TransitionPolicyTag Tag>
typename TransitionPolicyTraits<Tag>::Params bindTransitionPolicyParams(const TrackingParameters& params) noexcept;

template <>
inline CylinderCylinderPolicyParams bindTransitionPolicyParams<TransitionPolicyTag::CylinderCylinder>(const TrackingParameters& params) noexcept
{
  CylinderCylinderPolicyParams out;
  out.trackletMinPt = params.TrackletMinPt;
  out.cellDeltaTanLambdaSigma = params.CellDeltaTanLambdaSigma;
  out.nSigmaCut = params.NSigmaCut;
  out.maxChi2ClusterAttachment = params.MaxChi2ClusterAttachment;
  out.maxChi2NDF = params.MaxChi2NDF;
  return out;
}

template <>
inline DiskDiskPolicyParams bindTransitionPolicyParams<TransitionPolicyTag::DiskDisk>(const TrackingParameters& params) noexcept
{
  DiskDiskPolicyParams out;
  out.trackletMinPt = params.TrackletMinPt;
  out.cellDeltaTanLambdaSigma = params.CellDeltaTanLambdaSigma;
  out.cellRoadRCut = params.CellRoadRCut;
  out.trackletMinAbsX = params.TrackletMinAbsX;
  out.nSigmaCut = params.NSigmaCut;
  out.maxChi2ClusterAttachment = params.MaxChi2ClusterAttachment;
  out.maxChi2NDF = params.MaxChi2NDF;
  return out;
}

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYBINDING_H_ */

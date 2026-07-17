// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYBINDING_H_
#define ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYBINDING_H_

#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/TransitionPolicyState.h"

namespace o2::itsmft::tracking
{

/// Host-only: binds one iteration's legacy TrackingParameters into a typed,
/// bounds-checkable policy parameter block (TransitionPolicyState.h). Callers
/// must invoke this once per iteration, outside any candidate/neighbour/road
/// loop (D007 / Architecture.md 10.1) -- never per-candidate. TrackingParameters
/// owns std::vector members and is not itself device-compatible; only the
/// bound Params struct crosses into policy operations.
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

#endif /* ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYBINDING_H_ */

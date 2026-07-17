// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file TransitionPolicyOperations.cxx
/// \brief Out-of-line cellsAreCompatible<Tag> specializations (TransitionPolicyOperations.h)
///
/// Only this translation unit may include MFTFwdTrackHelpers.h on behalf of
/// the D007 policy operation boundary, so the common public policy header
/// stays free of MFT-specific constants/TimeFrame/MFTCATrack dependencies.

#include "ITSMFTTracking/TransitionPolicyOperations.h"

#include "ITSMFTTracking/MFTFwdTrackHelpers.h"

namespace o2::itsmft::tracking
{

template <>
bool cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(
  const CellSeedTpl<o2::track::TrackParCovF>& currentCell,
  const CellSeedTpl<o2::track::TrackParCovF>& nextCell,
  float bz,
  const CylinderCylinderPolicyParams& params)
{
  auto propagated = nextCell;
  if (!propagated.rotate(currentCell.getAlpha()) || !propagated.propagateTo(currentCell.getX(), bz)) {
    return false;
  }
  return currentCell.getPredictedChi2(propagated) <= params.maxChi2ClusterAttachment;
}

template <>
bool cellsAreCompatible<TransitionPolicyTag::DiskDisk>(
  const CellSeedTpl<o2::track::TrackParCovFwd>& currentCell,
  const CellSeedTpl<o2::track::TrackParCovFwd>& nextCell,
  float bz,
  const DiskDiskPolicyParams& params)
{
  return detail::mftFwdCellsAreCompatible(currentCell, nextCell, bz, params.maxChi2ClusterAttachment);
}

} // namespace o2::itsmft::tracking

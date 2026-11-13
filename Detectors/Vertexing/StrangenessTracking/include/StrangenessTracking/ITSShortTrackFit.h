// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_STRANGENESS_TRACKING_ITSSHORTTRACKFIT_H
#define O2_STRANGENESS_TRACKING_ITSSHORTTRACKFIT_H

#include "DetectorsBase/Propagator.h"
#include "ReconstructionDataFormats/BaseCluster.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/Track.h"
#include <array>

namespace o2::strangeness_tracking
{
// Input slots are physical layers 0, 1, 2. The ITS geometric seed is only an
// initializer: broad covariance precedes exactly three ITS measurement updates.
// Output is the PV-independent hit fit at L0. The PV is used only for pointing
// compatibility on a temporary DCA-propagated clone, with no PV update. Optional
// outerOut is the hit-fit clone transported to L2 without further measurements.
// Hit-only chi2/NDF uses NDF=1 and excludes the separate PV-pointing score.
// Outputs (including optional pointing DCA residuals) change only on success.
bool fitITSShortTrack(
  const std::array<o2::BaseCluster<float>, 3>& trackingClusters,
  const o2::dataformats::PrimaryVertex& pv, float bz,
  o2::base::Propagator::MatCorrType corr,
  const std::array<float, 3>& layerxX0, float minPt, float maxPointingChi2,
  float maxChi2NDF, o2::track::TrackParCovF& out, float& hitChi2,
  float& pointingChi2, int& ndf, std::array<float, 2>* pointingDCA = nullptr,
  o2::track::TrackParCovF* outerOut = nullptr);
} // namespace o2::strangeness_tracking

#endif

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITS_COMMON_TRACKING_PARAMETERS_H_
#define ALICEO2_ITS_COMMON_TRACKING_PARAMETERS_H_

#include <vector>

#include "ITSMFTTracking/Configuration.h"
#include "ITStracking/Configuration.h"

namespace o2::its::commontracking
{

/// Losslessly translates one already-constructed legacy ITS parameter set.
///
/// This function stores inactive/seeding masks and follower settings, but the
/// common tracker does not execute those semantics yet. Likewise,
/// FataliseUponFailure, benchmark/memory diagnostics, and complete UPC
/// vertex/mask orchestration are not provided by this compile-only onboarding
/// slice. A future opt-in workflow must reject requested unsupported behavior.
o2::itsmft::TrackingParameters translateTrackingParameters(const o2::its::TrackingParameters& legacy);

/// Calls the production legacy ITS constructor once and translates its result.
/// No mode normalization, override replay, field rescaling, or fallback occurs.
std::vector<o2::itsmft::TrackingParameters> makeTrackingParameters(o2::its::TrackingMode::Type mode);

} // namespace o2::its::commontracking

#endif // ALICEO2_ITS_COMMON_TRACKING_PARAMETERS_H_

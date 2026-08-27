// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file ConfigPreflight.h
/// \brief Driver-level configuration and vertex-constraint preflight for the
///        ITS common-CA tracker workflow.
///
/// These three checks are meant to be called, in order, from
/// defineDataProcessing() -- i.e. before any DataProcessorSpec/device is
/// constructed -- so a rejected configuration never reaches a running
/// device:
///
///   1. applyConfigKeyValuesOrFatal(rawConfigKeyValues)
///   2. requireSyncTrackingModeOrFatal(trackingMode)
///   3. requireDiamondVertexConstraintOrFatal()  (call AFTER step 1, since it
///      reads ITSCommonCATrackerParam::Instance() post-update)
///
/// These checks encode this workflow's policy: it has no per-event vertexing
/// capability, while the shared tracking configuration remains generic.

#ifndef ALICEO2_ITS_CA_WORKFLOW_CONFIGPREFLIGHT_H_
#define ALICEO2_ITS_CA_WORKFLOW_CONFIGPREFLIGHT_H_

#include <string>

#include "ITSMFTTracking/Configuration.h"

namespace o2::its::ca
{

/// Rejects a raw --configKeyValues string carrying an ITSCATrackerParam.*
/// override before applying the accepted string to ConfigurableParam.
void applyConfigKeyValuesOrFatal(const std::string& configKeyValues);

/// Fatals unless mode == Sync, naming the rejected mode explicitly, before
/// device construction.
void requireSyncTrackingModeOrFatal(o2::itsmft::TrackingMode::Type mode);

/// Fatals unless ITSCommonCATrackerParam::useDiamond is set. Call after
/// applyConfigKeyValuesOrFatal() so command-line overrides are observed.
void requireDiamondVertexConstraintOrFatal();

} // namespace o2::its::ca

#endif // ALICEO2_ITS_CA_WORKFLOW_CONFIGPREFLIGHT_H_

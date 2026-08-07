// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
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
/// \brief Driver-level configuration/vertex-constraint preflight for the
///        opt-in ITS common-CA tracker workflow (Gate 3 workflow-onboarding
///        Slice 2).
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
/// None of these belong in the shared common-tracking library: they encode
/// this *workflow's* policy (no real per-event vertexing capability yet),
/// not a general property of the standalone workflow or of
/// TrackingMode::getTrackingParameters() -- both of which continue to accept
/// useDiamond=false for ITS Sync (see workflow-onboarding Slice 1's own
/// tests), since a future caller with real vertexing may not need this
/// restriction.

#ifndef ALICEO2_ITS_CA_WORKFLOW_CONFIGPREFLIGHT_H_
#define ALICEO2_ITS_CA_WORKFLOW_CONFIGPREFLIGHT_H_

#include <string>

#include "ITSMFTTracking/Configuration.h"

namespace o2::its::ca
{

/// Rejects a raw --configKeyValues string carrying any legacy
/// "ITSCATrackerParam.*" override (see
/// o2::itsmft::tracking::checkITSCommonCAConfigKeyValues()) with a fatal
/// diagnostic naming "ITSCommonCATrackerParam" as the supported namespace,
/// *before* ever calling o2::conf::ConfigurableParam::updateFromString().
/// Only reaches updateFromString() -- applying it verbatim -- once the
/// preflight has accepted the string.
void applyConfigKeyValuesOrFatal(const std::string& configKeyValues);

/// Fatals unless mode == Sync, naming the rejected mode explicitly. Every
/// other o2::itsmft::TrackingMode::Type value (Off, Unset, Async, Cosmics)
/// is rejected -- none is silently mapped onto Sync. This is deliberately
/// redundant with (and runs strictly before) the identical fail-closed gate
/// already inside o2::itsmft::TrackingMode::getTrackingParameters(ITS, ...)
/// (workflow-onboarding Slice 1): that gate only fires once the device is
/// already running (CATrackerDPL::initialiseTracking(), at the first
/// TimeFrame), too late to satisfy "reject before device construction".
void requireSyncTrackingModeOrFatal(o2::itsmft::TrackingMode::Type mode);

/// Fatals unless o2::itsmft::ITSCommonCATrackerParam::Instance().useDiamond
/// is set. This is the one Sync vertex/beam-constraint mode the opt-in ITS
/// common-CA workflow can reproduce faithfully: the common tracker (no real
/// per-event vertexing capability for ITS) would otherwise silently run
/// tracklet/cell finding against an always-empty per-ROF primary-vertex
/// table (o2::itsmft::tracking::TrackerTraits<7>::computeLayerTracklets()'s
/// non-diamond branch), producing spuriously empty/degenerate tracking
/// instead of a loud failure. Must be called after
/// applyConfigKeyValuesOrFatal(), so it observes any --configKeyValues
/// override of ITSCommonCATrackerParam.useDiamond.
void requireDiamondVertexConstraintOrFatal();

} // namespace o2::its::ca

#endif // ALICEO2_ITS_CA_WORKFLOW_CONFIGPREFLIGHT_H_

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file ConfigPreflight.h
/// \brief Driver-level configuration preflight for the opt-in combined
///        ITS+MFT common-CA tracker workflow (Gate 4 C4).
///
/// Deliberately its own, independent implementation -- not a call into
/// ITSCAWorkflow/ConfigPreflight.h -- so this library never links
/// O2::ITSCAWorkflow merely to reuse its Sync/useDiamond checks (that would
/// be an unnecessary workflow-to-workflow dependency, conflicting with the
/// isolation rationale this whole library exists for). The *observable*
/// fatal conditions and messages match ITS's own opt-in route where the
/// same condition applies (Sync-only, ITSCommonCATrackerParam.useDiamond),
/// since the combined workflow's application plan inherits exactly the same
/// single-shared-layout, no-real-vertexer constraint ITS's own opt-in workflow
/// already documents.
///
/// Called, in order, from defineDataProcessing() -- before any
/// DataProcessorSpec/device is constructed -- so a rejected configuration
/// never reaches a running device:
///
///   1. requireCombinedTrackingEnabledOrFatal()
///   2. requireSyncTrackingModeOrFatal(trackingMode)
///   3. requireDiamondVertexConstraintOrFatal()
///   4. requireNoMFTIRFrameConfigOrFatal()
///
/// A fifth check, the MFT continuous-readout GRP fact, cannot be evaluated
/// here: it is CCDB-sourced and only available once the device is running
/// (see CombinedCATrackerSpec.h's own doc for where that runs instead).

#ifndef ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_CONFIGPREFLIGHT_H_
#define ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_CONFIGPREFLIGHT_H_

#include "ITSMFTTracking/Configuration.h"

namespace o2::itsmft::combined
{

/// Fatals unless o2::itsmft::ITSMFTCombinedCATrackerParam::Instance().enabled
/// is set. This opt-in combined workflow is a separate executable from
/// every ITS-only/MFT-only workflow already, but this flag adds a second,
/// config-level guard: a pipeline template that accidentally invokes this
/// binary without explicitly setting `enabled=true` fatals immediately
/// rather than silently running a combined ITS+MFT tracking pass no one
/// asked for.
void requireCombinedTrackingEnabledOrFatal();

/// Fatals unless mode == Sync, naming the rejected mode explicitly.
/// The combined workflow's application plan accepts exactly one
/// TrackingParameters iteration per detector (the shape both
/// TrackingMode::getTrackingParameters(ITS, Sync) and (MFT, Sync) produce);
/// every other mode either produces a different iteration count or fatals
/// on its own deeper inside TrackingMode::getTrackingParameters(), too late
/// to satisfy "reject before device construction".
void requireSyncTrackingModeOrFatal(o2::itsmft::TrackingMode::Type mode);

/// Fatals unless o2::itsmft::ITSCommonCATrackerParam::Instance().useDiamond
/// is set -- identical condition and rationale to
/// o2::its::ca::requireDiamondVertexConstraintOrFatal(): the combined
/// workflow's ITS leg has no real per-event vertexing capability either,
/// so a non-diamond configuration would otherwise silently run
/// tracklet/cell finding against an always-empty per-ROF primary vertex
/// table instead of failing loudly.
void requireDiamondVertexConstraintOrFatal();

/// Fatals if o2::mft::MFTTrackingParam::Instance().irFramesOnly is set.
/// The combined ITS+MFT common-CA tracking flow has no IR-frame parameter and
/// its ROF mask is always "every ROF enabled" -- unlike the single-detector
/// MFT opt-in workflow's own configureROFMask(), it never applies IR-frame
/// or multiplicity-cut filtering. This workflow deliberately declares no
/// IRFramesITS InputSpec at all (CombinedCATrackerSpec.cxx); this check is
/// the config-level half of refusing to silently run unmasked triggered
/// MFT data -- see requireContinuousMFTReadoutOrFatal() (called at runtime,
/// once GRPECS is available) for the other half.
void requireNoMFTIRFrameConfigOrFatal();

/// Fatals unless `continuousReadout` is true. Pure function over the
/// already-fetched
/// o2::base::GRPGeomHelper::instance().getGRPECS()->isDetContinuousReadOut(
/// o2::detectors::DetID::MFT) fact, so it is independently testable without
/// a real GRPECS CCDB object. Called once, from
/// CombinedCATrackerDPL::updateTimeDependentParams(), the first time GRPECS
/// is available -- necessarily at runtime, not from defineDataProcessing(),
/// since this fact is CCDB-sourced -- but strictly before the first
/// process() call, so no triggered MFT TF is ever processed. This is the
/// runtime half of refusing to silently run unmasked triggered MFT data;
/// see requireNoMFTIRFrameConfigOrFatal() for the config-level half.
void requireContinuousMFTReadoutOrFatal(bool continuousReadout);

} // namespace o2::itsmft::combined

#endif // ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_CONFIGPREFLIGHT_H_

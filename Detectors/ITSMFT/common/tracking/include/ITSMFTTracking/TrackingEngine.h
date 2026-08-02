// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// ADR 0007 (Detectors/ITSMFT/common/tracking/doc/decisions/
// 0007-generic-tracking-engine-boundary.md) decisions 3, 6:
// GenericTrackingEngineMigration.md M1's engine contract. One concrete
// TrackingEngine executes one event over an explicitly ordered
// TrackingParticipant (TrackingParticipant.h) schedule and applies the
// whole-event all-or-nothing failure contract Gate 4's
// CombinedTimeFrameCoordinator::process() already proved for the fixed
// ITS+MFT case, generalized here to an arbitrary participant collection.
// M1 is a contract/seam slice only: this engine is exercised in this slice
// by the focused contract test's fake participants, not yet by any
// production caller -- M2 wraps today's Tracker<NLayers>/TrackerTraits
// <NLayers>/LegacyTrackerScratch<NLayers> composition in a concrete
// TrackingParticipant and has CombinedTimeFrameCoordinator delegate to this
// engine without amending this header.

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKINGENGINE_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKINGENGINE_H_

#ifndef GPUCA_GPUCODE

#include <cstddef>
#include <vector>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "ITSMFTTracking/ParticipantId.h"
#include "ITSMFTTracking/TrackingParticipant.h"

namespace o2::itsmft::tracking
{

struct ParticipantEventResult {
  ParticipantId participant{};
  ParticipantOutcome outcome{ParticipantOutcome::Structural};
  std::size_t trackCount{0};
};

// Whole-event result (ADR 0007 decision 3's all-or-nothing contract,
// generalized from CombinedTimeFrameCoordinator::CombinedTrackingResult's
// fixed nITSTracks/nMFTTracks pair to an arbitrary ordered participant
// collection). `outcome` classifies the event as a whole: Success only if
// every scheduled participant's load() and track() both returned Success;
// otherwise the classification of whichever call first did not (see
// TrackingEngine::executeEvent()). `participants` carries one entry per
// scheduled participant, in schedule order, populated only on a fully
// successful event -- left empty on any failure, matching
// CombinedTrackingResult's existing zeroed-count convention.
struct EventResult {
  ParticipantOutcome outcome{ParticipantOutcome::Structural};
  std::vector<ParticipantEventResult> participants;
};

// ADR 0007 decision 3: the one concrete tracking-engine class. No
// interface/executor split -- a second public "TrackingExecutor"
// abstraction must not be introduced alongside this one.
class TrackingEngine
{
 public:
  // Executes one event over `schedule` against the one shared `frame`.
  // `schedule`'s order is the caller-supplied execution order (decision 6)
  // -- read as given, never sorted, deduplicated, or otherwise reordered by
  // this call, and never inferred from any static per-detector catalog or
  // from participant construction order. Every participant is loaded then
  // tracked, in that order; on any non-Success outcome, or any exception
  // thrown by either call, execution stops there, every participant in
  // `schedule` -- not only the one that failed -- has eventReset() called
  // on it exactly once, and this returns a zero-participant EventResult
  // carrying that failure's classification.
  EventResult executeEvent(TimeFrame& frame, gsl::span<TrackingParticipant* const> schedule, const o2::InteractionRecord& origin);
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE
#endif /* ALICEO2_ITSMFT_TRACKING_TRACKINGENGINE_H_ */

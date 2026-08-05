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
// TrackingEngine executes one already atomically loaded event over an
// explicitly ordered TrackingParticipant (TrackingParticipant.h) schedule
// and applies the whole-event all-or-nothing failure contract Gate 4's
// combined ITS+MFT tracking flow already proved for the fixed ITS+MFT case,
// generalized here to an arbitrary participant collection.
//
// TrackingEngine never loads anything itself. Atomic multi-source loading
// (every source's decoding and legacy backfill staged and committed
// together, or the live TimeFrame left completely untouched) is a distinct,
// stronger contract than executeEvent()'s own per-participant track() loop
// could honor if it also drove loading -- see
// TrackingParticipant.h's file-level comment. Today that atomic load is
// MultiSourceTimeFrameLoader::loadITSAndMFT(); M2 generalizes it into a
// participant-count-agnostic event loader and calls executeEvent() only
// after that loader reports success. This header does not change when that
// happens.
//
// M1 is a contract/seam slice only: this engine is exercised in this slice
// by the focused contract test's fake participants, not yet by any
// production caller -- M2 wraps today's Tracker<NLayers>/TrackerTraits
// <NLayers>/SurfaceTrackingScratch composition in a concrete
// TrackingParticipant; M3 has the combined DPL task
// (CombinedCATrackerSpec.cxx) compose loader + engine + the two concrete
// participants directly, all without amending this header.

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKINGENGINE_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKINGENGINE_H_

#ifndef GPUCA_GPUCODE

#include <cstddef>
#include <vector>

#include <gsl/gsl>

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
// generalized from the combined ITS+MFT tracking flow's own fixed
// nITSTracks/nMFTTracks pair to an arbitrary ordered participant
// collection). `outcome` classifies the event as a whole: Success only if
// every scheduled participant's track() returned Success; otherwise the
// classification of whichever call first did not (see
// TrackingEngine::executeEvent()). `participants` carries one entry per
// scheduled participant, in schedule order, populated only on a fully
// successful event -- left empty on any failure, matching that same
// zeroed-count convention.
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
  // Executes one already atomically loaded event: every scheduled
  // participant's decoding and legacy backfill must already be committed
  // into `frame` before this is called -- executeEvent() itself performs no
  // loading (see the file-level comment above). `schedule`'s order is the
  // caller-supplied execution order (decision 6) -- read as given, never
  // sorted, deduplicated, or otherwise reordered by this call, and never
  // inferred from any static per-detector catalog or from participant
  // construction order.
  //
  // Every participant's track() runs in that order. On any non-Success
  // outcome, or any exception thrown by track(), execution stops there,
  // resetEvent(frame, schedule) is called exactly once (resetting every
  // participant in `schedule` -- not only the one that failed -- and then
  // wiping `frame`'s shared storage), and this returns a zero-participant
  // EventResult carrying that failure's classification.
  EventResult executeEvent(TimeFrame& frame, gsl::span<TrackingParticipant* const> schedule);

  // Resets every participant in `schedule` (each participant's own
  // eventReset(), which by TrackingParticipant's contract touches only that
  // participant's own scratch/compatibility state) and then wipes `frame`'s
  // shared storage exactly once (TimeFrame::wipe()) -- always in that
  // order, so every participant has released its own per-event references
  // before the shared storage they pointed into is cleared.
  //
  // Public so a caller whose atomic event load itself failed -- before
  // executeEvent() was ever called, so no TrackingEngine call has touched
  // `frame` yet -- can reach the same all-participant/shared-frame reset
  // contract executeEvent() applies internally on a tracking failure,
  // without duplicating it.
  void resetEvent(TimeFrame& frame, gsl::span<TrackingParticipant* const> schedule) const noexcept;
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE
#endif /* ALICEO2_ITSMFT_TRACKING_TRACKINGENGINE_H_ */

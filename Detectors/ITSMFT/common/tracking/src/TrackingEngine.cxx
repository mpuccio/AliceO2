// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/TrackingEngine.h"

#include <exception>

namespace o2::itsmft::tracking
{

namespace
{
// `frame` is forwarded through unread by this function itself -- every
// participant's own eventReset() is what actually touches it. Taking it
// here (rather than each caller passing it individually) keeps
// executeEvent()'s two failure-handling branches identical.
void resetAll(gsl::span<TrackingParticipant* const> schedule, TimeFrame& frame) noexcept
{
  for (auto* participant : schedule) {
    participant->eventReset(frame);
  }
}
} // namespace

EventResult TrackingEngine::executeEvent(TimeFrame& frame, gsl::span<TrackingParticipant* const> schedule, const o2::InteractionRecord& origin)
{
  EventResult result;
  try {
    result.participants.reserve(schedule.size());
    for (auto* participant : schedule) {
      const auto loadResult = participant->load(frame, origin);
      if (loadResult.outcome != ParticipantOutcome::Success) {
        resetAll(schedule, frame);
        result.outcome = loadResult.outcome;
        result.participants.clear();
        return result;
      }
      const auto trackingResult = participant->track(frame);
      if (trackingResult.outcome != ParticipantOutcome::Success) {
        resetAll(schedule, frame);
        result.outcome = trackingResult.outcome;
        result.participants.clear();
        return result;
      }
      result.participants.push_back({participant->id(), trackingResult.outcome, trackingResult.trackCount});
    }
  } catch (const std::exception&) {
    resetAll(schedule, frame);
    result.outcome = ParticipantOutcome::Structural;
    result.participants.clear();
    return result;
  } catch (...) {
    resetAll(schedule, frame);
    result.outcome = ParticipantOutcome::Structural;
    result.participants.clear();
    return result;
  }

  result.outcome = ParticipantOutcome::Success;
  return result;
}

} // namespace o2::itsmft::tracking

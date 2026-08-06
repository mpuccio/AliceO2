// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/TrackingEngine.h"

#include <exception>

#include "ITSMFTTracking/TimeFrame.h"

namespace o2::itsmft::tracking
{

void TrackingEngine::resetEvent(TimeFrame& frame, gsl::span<TrackingParticipant* const> schedule) const noexcept
{
  // Participants release adapter-local references before the frame clears
  // the generic event storage they borrow.
  for (auto* participant : schedule) {
    participant->eventReset(frame);
  }
  frame.resetEvent();
}

EventResult TrackingEngine::executeEvent(TimeFrame& frame, gsl::span<TrackingParticipant* const> schedule)
{
  EventResult result;
  const auto resetCount = frame.getEventResetCount();
  const auto resetIfNeeded = [&]() noexcept {
    if (frame.getEventResetCount() == resetCount) {
      resetEvent(frame, schedule);
    }
  };
  try {
    result.participants.reserve(schedule.size());
    for (auto* participant : schedule) {
      const auto trackingResult = participant->track(frame);
      if (trackingResult.outcome != ParticipantOutcome::Success) {
        resetIfNeeded();
        result.outcome = trackingResult.outcome;
        result.participants.clear();
        return result;
      }
      result.participants.push_back({participant->id(), trackingResult.outcome, trackingResult.trackCount});
    }
  } catch (const std::exception&) {
    resetIfNeeded();
    result.outcome = ParticipantOutcome::Structural;
    result.participants.clear();
    return result;
  } catch (...) {
    resetIfNeeded();
    result.outcome = ParticipantOutcome::Structural;
    result.participants.clear();
    return result;
  }

  result.outcome = ParticipantOutcome::Success;
  return result;
}

} // namespace o2::itsmft::tracking

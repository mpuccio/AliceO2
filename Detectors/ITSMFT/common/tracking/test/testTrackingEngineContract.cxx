// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// GenericTrackingEngineMigration.md M1 focused contract test (ADR 0007
// decisions 3, 5, 6), corrected after the initial M1 slice: TrackingEngine
// no longer drives per-participant loading (that would violate the
// existing atomic multi-source loading contract -- see
// TrackingParticipant.h's file-level comment); it only runs an already
// atomically loaded event's track() calls, in schedule order. Covers:
//  - ParticipantId is opaque: no implicit conversion to/from a raw integer,
//    distinct from every other Identifier16 instantiation, and this whole
//    file -- which constructs, compares, and schedules ParticipantIds --
//    never includes DetectorsCommonDataFormats/DetID.h or any other
//    detector header;
//  - TrackingEngine::executeEvent() runs track() for a preloaded event in
//    exactly the caller-supplied schedule order, not participant id order,
//    declaration order, or any other inferred order;
//  - the whole-event all-or-nothing failure contract: a RecoverableDropped
//    or Structural track() outcome (via a typed return or a thrown
//    exception) stops execution there, resets every scheduled participant
//    -- including ones never reached -- wipes the shared TimeFrame exactly
//    once, and zeroes the result;
//  - TrackingEngine::resetEvent() -- the operation a caller whose atomic
//    event load itself failed uses directly, without ever calling
//    executeEvent() -- applies that exact same all-participant/shared-frame
//    contract;
//  - ParticipantOutcome/ParticipantTrackingResult/ParticipantEventResult/
//    EventResult/ParticipantPublicationExport are all usable, by this same
//    test, without any detector header.
//
// This file's own include list is deliberately restricted to
// ITSMFTTracking/TrackingEngine.h (which itself pulls in
// TrackingParticipant.h and ParticipantId.h), ITSMFTTracking/TimeFrame.h
// and ITSMFTTracking/CommonTrack.h (needed only to construct a real
// TimeFrame instance and to populate its shared storage directly, standing
// in for a real atomic load -- ADR 0007 decision 1 classifies both as
// permanent/detector-neutral core, not migration artifacts this slice must
// avoid), and ordinary standard/Boost.Test headers -- no Tracker.h,
// SurfaceTrackingScratch.h, TrackerTraits.h, TransitionPolicy*, DPL/workflow,
// or output-writer header. testTrackingEngineDependencyBoundary.cxx
// additionally scans the production headers themselves for those exact
// tokens.

#define BOOST_TEST_MODULE ITSMFT TrackingEngineContract
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <cstdint>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "ITSMFTTracking/CommonTrack.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingEngine.h"

using namespace o2::itsmft::tracking;

namespace
{

// Deliberately populates `frame`'s shared storage with self-consistent
// content, standing in for a real atomic load, so tests can observe that a
// whole-event failure wipes it. Content mirrors testCommonTrack.cxx's own
// populateCommonResults() fixture pattern.
void populateSharedFrame(TimeFrame& frame)
{
  frame.getTrackClusterIndices().push_back(TrackClusterReference{SurfaceId{0}, SurfaceMeasurementIndex{0}});
  CommonTrack track{};
  track.firstClusterRef = 0;
  track.clusterRefEnd = 1;
  track.hitSurfaces.set(SurfaceId{0});
  frame.getCommonTracks().push_back(track);
}

bool sharedFrameIsWiped(const TimeFrame& frame) noexcept
{
  return frame.getCommonTracks().empty() && frame.getTrackClusterIndices().empty();
}

// A fake, detector-free TrackingParticipant. Records every track()/
// eventReset() call, in call order, into a shared log so tests can assert
// the exact cross-participant execution order the engine used -- not just
// each result's final content.
class FakeParticipant final : public TrackingParticipant
{
 public:
  FakeParticipant(ParticipantId id, std::vector<std::string>& log,
                  ParticipantOutcome trackOutcome = ParticipantOutcome::Success,
                  std::size_t trackCount = 0)
    : mId{id}, mLog{log}, mTrackOutcome{trackOutcome}, mTrackCount{trackCount}
  {
  }

  ParticipantId id() const noexcept override { return mId; }
  gsl::span<const SurfaceId> ownedSurfaces() const noexcept override { return mSurfaces; }

  ParticipantTrackingResult track(TimeFrame&) override
  {
    ++mTrackCalls;
    mLog.get().push_back("track:" + std::to_string(mId.value()));
    return {mTrackOutcome, mTrackOutcome == ParticipantOutcome::Success ? mTrackCount : 0};
  }

  void eventReset(TimeFrame&) noexcept override
  {
    ++mResetCalls;
    mLog.get().push_back("reset:" + std::to_string(mId.value()));
  }

  std::optional<ParticipantPublicationExport> publicationExport() const override
  {
    if (mTrackOutcome != ParticipantOutcome::Success) {
      return std::nullopt;
    }
    return ParticipantPublicationExport{mId, gsl::span<const SurfaceId>{mSurfaces}};
  }

  int trackCalls() const noexcept { return mTrackCalls; }
  int resetCalls() const noexcept { return mResetCalls; }

 private:
  ParticipantId mId;
  std::reference_wrapper<std::vector<std::string>> mLog;
  std::vector<SurfaceId> mSurfaces{};
  ParticipantOutcome mTrackOutcome;
  std::size_t mTrackCount;
  int mTrackCalls{0};
  int mResetCalls{0};
};

class ThrowingParticipant final : public TrackingParticipant
{
 public:
  ThrowingParticipant(ParticipantId id, std::vector<std::string>& log) : mId{id}, mLog{log} {}

  ParticipantId id() const noexcept override { return mId; }
  gsl::span<const SurfaceId> ownedSurfaces() const noexcept override { return {}; }

  ParticipantTrackingResult track(TimeFrame&) override
  {
    mLog.get().push_back("track:" + std::to_string(mId.value()));
    throw std::runtime_error("ThrowingParticipant::track() structural failure");
  }

  void eventReset(TimeFrame&) noexcept override
  {
    ++mResetCalls;
    mLog.get().push_back("reset:" + std::to_string(mId.value()));
  }

  std::optional<ParticipantPublicationExport> publicationExport() const override { return std::nullopt; }

  int resetCalls() const noexcept { return mResetCalls; }

 private:
  ParticipantId mId;
  std::reference_wrapper<std::vector<std::string>> mLog;
  int mResetCalls{0};
};

} // namespace

// --- ParticipantId opacity -------------------------------------------------

BOOST_AUTO_TEST_CASE(ParticipantIdIsOpaque)
{
  static_assert(!std::is_convertible_v<int, ParticipantId>,
                "ParticipantId must not be implicitly constructible from a raw integer");
  static_assert(!std::is_convertible_v<ParticipantId, int>,
                "ParticipantId must not implicitly decay to a raw integer");
  static_assert(std::is_same_v<ParticipantId::ValueType, uint16_t>);

  const ParticipantId a{3};
  const ParticipantId b{3};
  const ParticipantId c{4};
  BOOST_CHECK(a == b);
  BOOST_CHECK(a != c);
  BOOST_CHECK(a < c);
  BOOST_CHECK(a.isValid());
  BOOST_CHECK_EQUAL(a.value(), 3u);
  BOOST_CHECK(!ParticipantId::invalid().isValid());
  BOOST_CHECK(!ParticipantId{}.isValid()); // default-constructed == invalid, same as every other Identifier16
}

// --- Explicit schedule order -----------------------------------------------

BOOST_AUTO_TEST_CASE(ScheduleOrderIsExplicitNotInferredOnAPreloadedEvent)
{
  std::vector<std::string> log;
  // Deliberately descending/unsorted ids, and a schedule order that matches
  // neither id order nor declaration order below.
  FakeParticipant high{ParticipantId{9}, log, ParticipantOutcome::Success, 5};
  FakeParticipant low{ParticipantId{1}, log, ParticipantOutcome::Success, 2};
  FakeParticipant mid{ParticipantId{4}, log, ParticipantOutcome::Success, 7};

  const std::array<TrackingParticipant*, 3> schedule{&mid, &high, &low};

  TimeFrame frame; // stands in for an already atomically loaded event
  TrackingEngine engine;
  const auto result = engine.executeEvent(frame, gsl::span<TrackingParticipant* const>{schedule.data(), schedule.size()});

  BOOST_CHECK(result.outcome == ParticipantOutcome::Success);
  BOOST_REQUIRE_EQUAL(result.participants.size(), 3u);
  BOOST_CHECK(result.participants[0].participant == mid.id());
  BOOST_CHECK(result.participants[1].participant == high.id());
  BOOST_CHECK(result.participants[2].participant == low.id());
  BOOST_CHECK_EQUAL(result.participants[0].trackCount, 7u);
  BOOST_CHECK_EQUAL(result.participants[1].trackCount, 5u);
  BOOST_CHECK_EQUAL(result.participants[2].trackCount, 2u);

  const std::vector<std::string> expectedLog{"track:4", "track:9", "track:1"};
  BOOST_CHECK_EQUAL_COLLECTIONS(log.begin(), log.end(), expectedLog.begin(), expectedLog.end());
}

// --- Whole-event all-or-nothing failure contract ---------------------------

BOOST_AUTO_TEST_CASE(RecoverableDroppedTrackingResultResetsEveryScheduledParticipantAndWipesFrame)
{
  std::vector<std::string> log;
  FakeParticipant first{ParticipantId{0}, log};
  FakeParticipant second{ParticipantId{1}, log, ParticipantOutcome::RecoverableDropped};
  FakeParticipant third{ParticipantId{2}, log};
  const std::array<TrackingParticipant*, 3> schedule{&first, &second, &third};

  TimeFrame frame;
  populateSharedFrame(frame);
  BOOST_REQUIRE(!sharedFrameIsWiped(frame));

  TrackingEngine engine;
  const auto result = engine.executeEvent(frame, gsl::span<TrackingParticipant* const>{schedule.data(), schedule.size()});

  BOOST_CHECK(result.outcome == ParticipantOutcome::RecoverableDropped);
  BOOST_CHECK(result.participants.empty());
  BOOST_CHECK(sharedFrameIsWiped(frame));

  BOOST_CHECK_EQUAL(first.trackCalls(), 1);
  BOOST_CHECK_EQUAL(second.trackCalls(), 1);
  // `third` is never reached -- execution stops at `second` -- yet it is
  // still reset, proving the contract resets every scheduled participant,
  // not only the ones actually executed.
  BOOST_CHECK_EQUAL(third.trackCalls(), 0);

  BOOST_CHECK_EQUAL(first.resetCalls(), 1);
  BOOST_CHECK_EQUAL(second.resetCalls(), 1);
  BOOST_CHECK_EQUAL(third.resetCalls(), 1);
  BOOST_CHECK_EQUAL(frame.getEventResetCount(), 1u);
}

BOOST_AUTO_TEST_CASE(StructuralTrackOutcomeResetsEveryScheduledParticipantAndWipesFrame)
{
  std::vector<std::string> log;
  FakeParticipant first{ParticipantId{0}, log, ParticipantOutcome::Structural};
  FakeParticipant second{ParticipantId{1}, log};
  const std::array<TrackingParticipant*, 2> schedule{&first, &second};

  TimeFrame frame;
  populateSharedFrame(frame);

  TrackingEngine engine;
  const auto result = engine.executeEvent(frame, gsl::span<TrackingParticipant* const>{schedule.data(), schedule.size()});

  BOOST_CHECK(result.outcome == ParticipantOutcome::Structural);
  BOOST_CHECK(result.participants.empty());
  BOOST_CHECK(sharedFrameIsWiped(frame));
  BOOST_CHECK_EQUAL(second.trackCalls(), 0); // never reached
  BOOST_CHECK_EQUAL(first.resetCalls(), 1);
  BOOST_CHECK_EQUAL(second.resetCalls(), 1);
  BOOST_CHECK_EQUAL(frame.getEventResetCount(), 1u);
}

BOOST_AUTO_TEST_CASE(ExceptionFromAParticipantIsStructuralAndResetsEveryScheduledParticipantAndWipesFrame)
{
  std::vector<std::string> log;
  FakeParticipant first{ParticipantId{0}, log};
  ThrowingParticipant second{ParticipantId{1}, log};
  FakeParticipant third{ParticipantId{2}, log};
  const std::array<TrackingParticipant*, 3> schedule{&first, &second, &third};

  TimeFrame frame;
  populateSharedFrame(frame);

  TrackingEngine engine;
  const auto result = engine.executeEvent(frame, gsl::span<TrackingParticipant* const>{schedule.data(), schedule.size()});

  BOOST_CHECK(result.outcome == ParticipantOutcome::Structural);
  BOOST_CHECK(result.participants.empty());
  BOOST_CHECK(sharedFrameIsWiped(frame));
  BOOST_CHECK_EQUAL(third.trackCalls(), 0); // never reached

  BOOST_CHECK_EQUAL(first.resetCalls(), 1);
  BOOST_CHECK_EQUAL(second.resetCalls(), 1);
  BOOST_CHECK_EQUAL(third.resetCalls(), 1);
  BOOST_CHECK_EQUAL(frame.getEventResetCount(), 1u);
}

// --- Direct resetEvent() (used after an atomic load failure) ---------------

BOOST_AUTO_TEST_CASE(DirectResetEventAppliesTheSameAllParticipantAndSharedFrameContract)
{
  std::vector<std::string> log;
  FakeParticipant first{ParticipantId{0}, log};
  FakeParticipant second{ParticipantId{1}, log};
  const std::array<TrackingParticipant*, 2> schedule{&first, &second};

  TimeFrame frame;
  populateSharedFrame(frame);
  BOOST_REQUIRE(!sharedFrameIsWiped(frame));

  // No executeEvent() call at all: this is the path a caller whose own
  // atomic event load already failed uses directly.
  TrackingEngine engine;
  engine.resetEvent(frame, gsl::span<TrackingParticipant* const>{schedule.data(), schedule.size()});

  BOOST_CHECK(sharedFrameIsWiped(frame));
  BOOST_CHECK_EQUAL(first.trackCalls(), 0);
  BOOST_CHECK_EQUAL(second.trackCalls(), 0);
  BOOST_CHECK_EQUAL(first.resetCalls(), 1);
  BOOST_CHECK_EQUAL(second.resetCalls(), 1);
  BOOST_CHECK_EQUAL(frame.getEventResetCount(), 1u);

  const std::vector<std::string> expectedLog{"reset:0", "reset:1"};
  BOOST_CHECK_EQUAL_COLLECTIONS(log.begin(), log.end(), expectedLog.begin(), expectedLog.end());
}

// --- Result/outcome representations are usable without detector headers ---
// (structural claim: this whole file, up to and including this test case,
// never includes DetectorsCommonDataFormats/DetID.h or any workflow/writer
// header -- see the file-level comment above.)

BOOST_AUTO_TEST_CASE(ResultAndOutcomeTypesAreUsableWithoutDetectorHeaders)
{
  static_assert(std::is_same_v<decltype(ParticipantTrackingResult::outcome), ParticipantOutcome>);
  static_assert(std::is_same_v<decltype(ParticipantTrackingResult::trackCount), std::size_t>);
  static_assert(std::is_same_v<decltype(ParticipantEventResult::participant), ParticipantId>);
  static_assert(std::is_same_v<decltype(ParticipantEventResult::outcome), ParticipantOutcome>);
  static_assert(std::is_same_v<decltype(EventResult::outcome), ParticipantOutcome>);
  static_assert(std::is_same_v<decltype(EventResult::participants), std::vector<ParticipantEventResult>>);
  static_assert(std::is_same_v<decltype(ParticipantPublicationExport::id), ParticipantId>);

  const ParticipantTrackingResult trackingResult{ParticipantOutcome::Success, 42};
  BOOST_CHECK(trackingResult.outcome == ParticipantOutcome::Success);
  BOOST_CHECK_EQUAL(trackingResult.trackCount, 42u);

  const EventResult eventResult{};
  BOOST_CHECK(eventResult.outcome == ParticipantOutcome::Structural); // default
  BOOST_CHECK(eventResult.participants.empty());

  const ParticipantPublicationExport publicationExport{ParticipantId{7}, {}};
  BOOST_CHECK(publicationExport.id == ParticipantId{7});
  BOOST_CHECK(publicationExport.orderedSurfaces.empty());
}

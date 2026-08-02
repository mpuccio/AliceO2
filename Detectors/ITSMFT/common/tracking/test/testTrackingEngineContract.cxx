// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// GenericTrackingEngineMigration.md M1 focused contract test (ADR 0007
// decisions 3, 5, 6). Covers:
//  - ParticipantId is opaque: no implicit conversion to/from a raw integer,
//    distinct from every other Identifier16 instantiation, and this whole
//    file -- which constructs, compares, and schedules ParticipantIds --
//    never includes DetectorsCommonDataFormats/DetID.h or any other
//    detector header;
//  - TrackingEngine::executeEvent() executes a schedule in exactly the
//    caller-supplied order, not participant id order, declaration order, or
//    any other inferred order;
//  - the whole-event all-or-nothing failure contract: a RecoverableDropped
//    or Structural outcome from any scheduled participant (via a typed
//    return or a thrown exception) stops execution there, resets every
//    scheduled participant -- including ones never reached -- and zeroes
//    the result;
//  - ParticipantOutcome/ParticipantLoadResult/ParticipantTrackingResult/
//    ParticipantEventResult/EventResult/ParticipantPublicationExport are
//    all usable, by this same test, without any detector header.
//
// This file's own include list is deliberately restricted to
// ITSMFTTracking/TrackingEngine.h (which itself pulls in
// TrackingParticipant.h and ParticipantId.h), ITSMFTTracking/TimeFrame.h
// (needed only to construct a real TimeFrame instance to pass by
// reference -- ADR 0007 decision 1 classifies TimeFrame itself as
// permanent/detector-neutral core, not a migration artifact this slice
// must avoid), and ordinary standard/Boost.Test headers -- no
// CATracker.h, LegacyTrackerScratch.h, TrackerTraits.h, TransitionPolicy*,
// DPL/workflow, or output-writer header. testTrackingEngineDependencyBoundary
// .cxx additionally scans the production headers themselves for those exact
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

#include "CommonDataFormat/InteractionRecord.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingEngine.h"

using namespace o2::itsmft::tracking;

namespace
{

// A fake, detector-free TrackingParticipant. Records every load()/track()/
// eventReset() call, in call order, into a shared log so tests can assert
// the exact cross-participant execution order the engine used -- not just
// each result's final content.
class FakeParticipant final : public TrackingParticipant
{
 public:
  FakeParticipant(ParticipantId id, std::vector<std::string>& log,
                  ParticipantOutcome loadOutcome = ParticipantOutcome::Success,
                  ParticipantOutcome trackOutcome = ParticipantOutcome::Success,
                  std::size_t trackCount = 0)
    : mId{id}, mLog{log}, mLoadOutcome{loadOutcome}, mTrackOutcome{trackOutcome}, mTrackCount{trackCount}
  {
  }

  ParticipantId id() const noexcept override { return mId; }
  gsl::span<const SurfaceId> ownedSurfaces() const noexcept override { return mSurfaces; }

  ParticipantLoadResult load(TimeFrame&, const o2::InteractionRecord&) override
  {
    ++mLoadCalls;
    mLog.get().push_back("load:" + std::to_string(mId.value()));
    return {mLoadOutcome};
  }

  ParticipantTrackingResult track(TimeFrame&) override
  {
    BOOST_REQUIRE_MESSAGE(mLoadOutcome == ParticipantOutcome::Success,
                          "track() must never be called after this participant's own non-Success load()");
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

  int loadCalls() const noexcept { return mLoadCalls; }
  int trackCalls() const noexcept { return mTrackCalls; }
  int resetCalls() const noexcept { return mResetCalls; }

 private:
  ParticipantId mId;
  std::reference_wrapper<std::vector<std::string>> mLog;
  std::vector<SurfaceId> mSurfaces{};
  ParticipantOutcome mLoadOutcome;
  ParticipantOutcome mTrackOutcome;
  std::size_t mTrackCount;
  int mLoadCalls{0};
  int mTrackCalls{0};
  int mResetCalls{0};
};

class ThrowingParticipant final : public TrackingParticipant
{
 public:
  ThrowingParticipant(ParticipantId id, std::vector<std::string>& log) : mId{id}, mLog{log} {}

  ParticipantId id() const noexcept override { return mId; }
  gsl::span<const SurfaceId> ownedSurfaces() const noexcept override { return {}; }

  ParticipantLoadResult load(TimeFrame&, const o2::InteractionRecord&) override
  {
    mLog.get().push_back("load:" + std::to_string(mId.value()));
    throw std::runtime_error("ThrowingParticipant::load() structural failure");
  }

  ParticipantTrackingResult track(TimeFrame&) override
  {
    BOOST_FAIL("track() must never be called on a participant whose load() threw");
    return {};
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

BOOST_AUTO_TEST_CASE(ScheduleOrderIsExplicitNotInferred)
{
  std::vector<std::string> log;
  // Deliberately descending/unsorted ids, and a schedule order that matches
  // neither id order nor declaration order below.
  FakeParticipant high{ParticipantId{9}, log, ParticipantOutcome::Success, ParticipantOutcome::Success, 5};
  FakeParticipant low{ParticipantId{1}, log, ParticipantOutcome::Success, ParticipantOutcome::Success, 2};
  FakeParticipant mid{ParticipantId{4}, log, ParticipantOutcome::Success, ParticipantOutcome::Success, 7};

  const std::array<TrackingParticipant*, 3> schedule{&mid, &high, &low};

  TimeFrame frame;
  TrackingEngine engine;
  const auto result = engine.executeEvent(frame, gsl::span<TrackingParticipant* const>{schedule.data(), schedule.size()}, o2::InteractionRecord{0, 0});

  BOOST_CHECK(result.outcome == ParticipantOutcome::Success);
  BOOST_REQUIRE_EQUAL(result.participants.size(), 3u);
  BOOST_CHECK(result.participants[0].participant == mid.id());
  BOOST_CHECK(result.participants[1].participant == high.id());
  BOOST_CHECK(result.participants[2].participant == low.id());
  BOOST_CHECK_EQUAL(result.participants[0].trackCount, 7u);
  BOOST_CHECK_EQUAL(result.participants[1].trackCount, 5u);
  BOOST_CHECK_EQUAL(result.participants[2].trackCount, 2u);

  const std::vector<std::string> expectedLog{"load:4", "track:4", "load:9", "track:9", "load:1", "track:1"};
  BOOST_CHECK_EQUAL_COLLECTIONS(log.begin(), log.end(), expectedLog.begin(), expectedLog.end());
}

// --- Whole-event all-or-nothing failure contract ---------------------------

BOOST_AUTO_TEST_CASE(RecoverableDroppedTrackingResultResetsEveryScheduledParticipant)
{
  std::vector<std::string> log;
  FakeParticipant first{ParticipantId{0}, log};
  FakeParticipant second{ParticipantId{1}, log, ParticipantOutcome::Success, ParticipantOutcome::RecoverableDropped};
  FakeParticipant third{ParticipantId{2}, log};
  const std::array<TrackingParticipant*, 3> schedule{&first, &second, &third};

  TimeFrame frame;
  TrackingEngine engine;
  const auto result = engine.executeEvent(frame, gsl::span<TrackingParticipant* const>{schedule.data(), schedule.size()}, o2::InteractionRecord{0, 0});

  BOOST_CHECK(result.outcome == ParticipantOutcome::RecoverableDropped);
  BOOST_CHECK(result.participants.empty());

  BOOST_CHECK_EQUAL(first.loadCalls(), 1);
  BOOST_CHECK_EQUAL(first.trackCalls(), 1);
  BOOST_CHECK_EQUAL(second.loadCalls(), 1);
  BOOST_CHECK_EQUAL(second.trackCalls(), 1);
  // `third` is never reached -- execution stops at `second` -- yet it is
  // still reset, proving the contract resets every scheduled participant,
  // not only the ones actually executed.
  BOOST_CHECK_EQUAL(third.loadCalls(), 0);
  BOOST_CHECK_EQUAL(third.trackCalls(), 0);

  BOOST_CHECK_EQUAL(first.resetCalls(), 1);
  BOOST_CHECK_EQUAL(second.resetCalls(), 1);
  BOOST_CHECK_EQUAL(third.resetCalls(), 1);
}

BOOST_AUTO_TEST_CASE(StructuralLoadOutcomeResetsEveryScheduledParticipant)
{
  std::vector<std::string> log;
  FakeParticipant first{ParticipantId{0}, log, ParticipantOutcome::Structural};
  FakeParticipant second{ParticipantId{1}, log};
  const std::array<TrackingParticipant*, 2> schedule{&first, &second};

  TimeFrame frame;
  TrackingEngine engine;
  const auto result = engine.executeEvent(frame, gsl::span<TrackingParticipant* const>{schedule.data(), schedule.size()}, o2::InteractionRecord{0, 0});

  BOOST_CHECK(result.outcome == ParticipantOutcome::Structural);
  BOOST_CHECK(result.participants.empty());
  BOOST_CHECK_EQUAL(first.trackCalls(), 0); // load() itself already failed
  BOOST_CHECK_EQUAL(second.loadCalls(), 0); // never reached
  BOOST_CHECK_EQUAL(first.resetCalls(), 1);
  BOOST_CHECK_EQUAL(second.resetCalls(), 1);
}

BOOST_AUTO_TEST_CASE(ExceptionFromAParticipantIsStructuralAndResetsEveryScheduledParticipant)
{
  std::vector<std::string> log;
  FakeParticipant first{ParticipantId{0}, log};
  ThrowingParticipant second{ParticipantId{1}, log};
  FakeParticipant third{ParticipantId{2}, log};
  const std::array<TrackingParticipant*, 3> schedule{&first, &second, &third};

  TimeFrame frame;
  TrackingEngine engine;
  const auto result = engine.executeEvent(frame, gsl::span<TrackingParticipant* const>{schedule.data(), schedule.size()}, o2::InteractionRecord{0, 0});

  BOOST_CHECK(result.outcome == ParticipantOutcome::Structural);
  BOOST_CHECK(result.participants.empty());
  BOOST_CHECK_EQUAL(third.loadCalls(), 0); // never reached

  BOOST_CHECK_EQUAL(first.resetCalls(), 1);
  BOOST_CHECK_EQUAL(second.resetCalls(), 1);
  BOOST_CHECK_EQUAL(third.resetCalls(), 1);
}

// --- Result/outcome representations are usable without detector headers ---
// (structural claim: this whole file, up to and including this test case,
// never includes DetectorsCommonDataFormats/DetID.h or any workflow/writer
// header -- see the file-level comment above.)

BOOST_AUTO_TEST_CASE(ResultAndOutcomeTypesAreUsableWithoutDetectorHeaders)
{
  static_assert(std::is_same_v<decltype(ParticipantLoadResult::outcome), ParticipantOutcome>);
  static_assert(std::is_same_v<decltype(ParticipantTrackingResult::outcome), ParticipantOutcome>);
  static_assert(std::is_same_v<decltype(ParticipantTrackingResult::trackCount), std::size_t>);
  static_assert(std::is_same_v<decltype(ParticipantEventResult::participant), ParticipantId>);
  static_assert(std::is_same_v<decltype(ParticipantEventResult::outcome), ParticipantOutcome>);
  static_assert(std::is_same_v<decltype(EventResult::outcome), ParticipantOutcome>);
  static_assert(std::is_same_v<decltype(EventResult::participants), std::vector<ParticipantEventResult>>);
  static_assert(std::is_same_v<decltype(ParticipantPublicationExport::id), ParticipantId>);

  const ParticipantLoadResult loadResult{};
  BOOST_CHECK(loadResult.outcome == ParticipantOutcome::Structural); // default

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

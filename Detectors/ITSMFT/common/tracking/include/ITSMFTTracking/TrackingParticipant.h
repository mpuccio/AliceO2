// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// ADR 0007 (Detectors/ITSMFT/common/tracking/doc/decisions/
// 0007-generic-tracking-engine-boundary.md) decision 5: the detector-leg
// boundary. GenericTrackingEngineMigration.md M1 commits this contract
// additively; M2 implements a concrete legacy-CA participant against it
// without amending it. Its minimal interface -- identity, owned surfaces,
// load, track, event reset, publication export -- must not mention ITS,
// MFT, NLayers, LegacyTrackerScratch, source 0/1, or any current output
// type (decision 5); every one of those stays inside concrete participant
// implementations and their own adapters.

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKINGPARTICIPANT_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKINGPARTICIPANT_H_

#ifndef GPUCA_GPUCODE

#include <cstddef>
#include <cstdint>
#include <optional>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "ITSMFTTracking/ParticipantId.h"
#include "ITSMFTTracking/SurfaceId.h"

namespace o2::itsmft::tracking
{

// ADR 0007 decision 1: the permanent, detector-neutral event owner. Only
// forward-declared here -- TrackingParticipant/TrackingEngine's public
// contract stays usable, and its own dependency-boundary test compilable,
// without pulling in TimeFrame.h's own ITStracking-namespaced
// infrastructure includes (BoundedAllocator.h, Cluster.h). Concrete
// TrackingParticipant implementations and TrackingEngine.cxx include the
// full ITSMFTTracking/TimeFrame.h where the complete type is actually
// needed.
struct TimeFrame;

// Reuses TrackingOutcome's (CATracker.h) three-value vocabulary --
// Success; RecoverableDropped, a per-TF data failure a participant's own
// policy opted to drop; Structural, everything else, including any
// exception (see TrackingEngine::executeEvent()) -- but redeclares it here
// rather than including either CATracker.h or
// CombinedTimeFrameCoordinator.h: both transitively pull in
// LegacyTrackerScratch.h/TrackerTraits.h/TimeFrame.h/
// DetectorTraversalBinding.h, exactly the headers this public contract must
// stay free of.
enum class ParticipantOutcome : uint8_t {
  Success,
  RecoverableDropped,
  Structural
};

struct ParticipantLoadResult {
  ParticipantOutcome outcome{ParticipantOutcome::Structural};
};

struct ParticipantTrackingResult {
  ParticipantOutcome outcome{ParticipantOutcome::Structural};
  std::size_t trackCount{0};
};

// Generic per-participant publication-export boundary: identity plus the
// ordered surface span the engine's schedule already carries. Deliberately
// excludes detector identity, cluster-source id, and any clock/timing type
// -- TrackingInterface.h's CommonTrackPublicationExport (which carries
// o2::detectors::DetID::ID and ClockTimingPublicationView) remains the
// existing, untouched adapter bridge; a detector adapter builds that richer
// export from this plus its own concrete TrackingParticipant.
struct ParticipantPublicationExport {
  ParticipantId id{};
  gsl::span<const SurfaceId> orderedSurfaces;
};

// The detector-leg boundary (ADR 0007 decision 5). A pure interface: the
// engine's schedule (TrackingEngine::executeEvent()) is a runtime-ordered,
// heterogeneous collection of participants -- an ITS leg and an MFT leg
// today, a future ALICE 3 leg tomorrow -- built from plan data the engine
// itself does not choose (decision 6). That heterogeneity is decided at run
// time, not compile time, which is exactly the case dynamic polymorphism
// exists for: a compile-time mechanism (templates, a variant/visitor closed
// over today's two detectors) cannot represent "one more concrete kind
// added later without recompiling the engine or widening a closed set".
// This is the only abstraction this slice introduces beyond plain data, and
// it is deliberately minimal: six operations, no shared state, no
// registry, no factory.
class TrackingParticipant
{
 public:
  virtual ~TrackingParticipant() = default;

  virtual ParticipantId id() const noexcept = 0;
  virtual gsl::span<const SurfaceId> ownedSurfaces() const noexcept = 0;

  // Loads this participant's own per-event input into `frame` for the
  // event at `origin`. How the participant obtains that input -- which
  // cluster source(s), which decoder -- is entirely the concrete adapter's
  // own business, never a fixed source index or DPL input shape mentioned
  // by this contract.
  virtual ParticipantLoadResult load(TimeFrame& frame, const o2::InteractionRecord& origin) = 0;

  // Tracks this participant's own owned surfaces against `frame`'s current
  // normalized content, appending accepted results to `frame`'s shared
  // CommonTrack/TrackClusterReference storage.
  virtual ParticipantTrackingResult track(TimeFrame& frame) = 0;

  // Clears this participant's own per-event state. Never responsible for
  // any other participant's state, or for any of `frame`'s shared storage
  // beyond this participant's own contribution -- applying the whole-event
  // all-or-nothing contract across an entire schedule is
  // TrackingEngine::executeEvent()'s job, not this method's.
  virtual void eventReset(TimeFrame& frame) noexcept = 0;

  // Engaged only after this participant's own successful track(); disengaged
  // by eventReset() and by any state prior to a first successful track().
  virtual std::optional<ParticipantPublicationExport> publicationExport() const = 0;
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE
#endif /* ALICEO2_ITSMFT_TRACKING_TRACKINGPARTICIPANT_H_ */

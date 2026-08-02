// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// ADR 0007 (Detectors/ITSMFT/common/tracking/doc/decisions/
// 0007-generic-tracking-engine-boundary.md) decisions 3, 5, 6:
// GenericTrackingEngineMigration.md M1's opaque TrackingEngine-schedule
// identity.

#ifndef ALICEO2_ITSMFT_TRACKING_PARTICIPANTID_H_
#define ALICEO2_ITSMFT_TRACKING_PARTICIPANTID_H_

#include <cstdint>
#include <type_traits>

#include "ITSMFTTracking/SurfaceId.h"

namespace o2::itsmft::tracking
{

struct ParticipantIdTag;

// Reuses the same tag-dispatched detail::Identifier16<Tag> the common core
// already uses for SurfaceId/TransitionId/CellTopologyId/ClusterSourceId
// (SurfaceId.h) rather than introducing a second bespoke opaque-id
// mechanism. Deliberately not o2::detectors::DetID::ID: DetID names a
// fixed, closed detector enumeration the generic core must never depend on
// (ADR 0007 decisions 2, 5); ParticipantId instead names a position in one
// event's TrackingEngine schedule, assigned by whatever builds that
// schedule (a future adapter/plan, not this contract, see decision 6). Two
// participants for the same detector, or participants for two different
// detectors, are simply two different ParticipantId values -- the type
// itself carries no detector semantics, and (like every Identifier16
// instantiation) has no implicit conversion to or from a raw integer or any
// other id type.
using ParticipantId = detail::Identifier16<ParticipantIdTag>;

static_assert(std::is_standard_layout_v<ParticipantId> && std::is_trivially_copyable_v<ParticipantId>);
static_assert(sizeof(ParticipantId) == sizeof(uint16_t));

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_PARTICIPANTID_H_ */

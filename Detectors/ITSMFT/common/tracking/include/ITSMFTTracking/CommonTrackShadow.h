// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_COMMONTRACKSHADOW_H_
#define ALICEO2_ITSMFT_TRACKING_COMMONTRACKSHADOW_H_

#ifndef GPUCA_GPUCODE

#include <cstddef>
#include <vector>

#include "ITSMFTTracking/CommonTrack.h"
#include "ITSMFTTracking/TimeFrame.h"

namespace o2::itsmft::tracking
{

// A fully local, validated candidate for one serial owner-thread publication.
// No TimeFrame storage is touched until publishCommonTrackShadow() below.
struct CommonTrackShadowRecord {
  CommonTrack track{};
  std::vector<TrackClusterReference> references;
};

enum class CommonTrackShadowPublishStep : uint8_t {
  BeforeReferenceReserve,
  BeforeTrackReserve,
  BeforeReferences,
  BeforeTrack
};

inline bool validateCommonTrackShadowRecord(const TimeFrame& frame, const CommonTrackShadowRecord& record) noexcept
{
  if (!record.track.innerState.hasRecognizedFamily() || !record.track.outerState.hasRecognizedFamily() ||
      record.track.innerState.family != record.track.outerState.family || !record.track.timestamp.isValid() ||
      record.references.empty()) {
    return false;
  }
  SurfaceMask surfaces;
  for (const auto& reference : record.references) {
    const auto* measurement = frame.getNormalizedFrame().getMeasurement(reference.surface, reference.index);
    if (measurement == nullptr || measurement->surface != reference.surface) {
      return false;
    }
    surfaces.set(reference.surface);
  }
  return surfaces == record.track.hitSurfaces;
}

// The caller is the serial acceptTracks() owner thread. The optional callback
// is deliberately a test seam: production passes a no-op and never needs a
// lock. Both collections retain their original content and sizes if reserve,
// validation, or either append throws; CommonTrack is appended last so an
// incomplete range can never become visible as a track.
template <typename PublishHook>
bool publishCommonTrackShadow(TimeFrame& frame, const CommonTrackShadowRecord& record, PublishHook&& hook)
{
  if (!validateCommonTrackShadowRecord(frame, record)) {
    return false;
  }

  auto& tracks = frame.getCommonTracks();
  auto& references = frame.getTrackClusterIndices();
  const auto oldTrackSize = tracks.size();
  const auto oldReferenceSize = references.size();

  try {
    // Both fallible allocations precede every append. A reserve failure leaves
    // the collections' logical contents unchanged.
    hook(CommonTrackShadowPublishStep::BeforeReferenceReserve);
    references.reserve(oldReferenceSize + record.references.size());
    hook(CommonTrackShadowPublishStep::BeforeTrackReserve);
    tracks.reserve(oldTrackSize + 1);
    hook(CommonTrackShadowPublishStep::BeforeReferences);
    for (const auto& reference : record.references) {
      references.push_back(reference);
    }
    CommonTrack committed = record.track;
    committed.firstClusterRef = static_cast<uint32_t>(oldReferenceSize);
    committed.clusterRefEnd = static_cast<uint32_t>(references.size());
    hook(CommonTrackShadowPublishStep::BeforeTrack);
    tracks.push_back(committed);
  } catch (...) {
    references.resize(oldReferenceSize);
    tracks.resize(oldTrackSize);
    throw;
  }
  return true;
}

inline bool publishCommonTrackShadow(TimeFrame& frame, const CommonTrackShadowRecord& record)
{
  return publishCommonTrackShadow(frame, record, [](CommonTrackShadowPublishStep) {});
}

} // namespace o2::itsmft::tracking

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_COMMONTRACKSHADOW_H_

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_ACCEPTEDTRACKSHADOWPUBLISHER_H_
#define ALICEO2_ITSMFT_TRACKING_ACCEPTEDTRACKSHADOWPUBLISHER_H_

#ifndef GPUCA_GPUCODE

#include <optional>

#include "ITSMFTTracking/CommonTrackShadow.h"
#include "ITSMFTTracking/MFTCATrack.h"
#include "ITSMFTTracking/MFTPublicationCompatibility.h"

namespace o2::itsmft::tracking
{

// Detector-neutral default: ITS publishes only its CommonTrack shadow and
// never allocates, populates, or queries MFT compatibility state.
template <int NLayers>
class AcceptedTrackShadowPublisher
{
 public:
  template <typename Track>
  std::optional<uint32_t> publish(TimeFrame& frame, const CommonTrackShadowRecord& record, const Track&) const
  {
    return publishCommonTrackShadow(frame, record);
  }

  void adoptMFTPublicationCompatibility(MFTPublicationCompatibility*) noexcept {}
};

// This specialization is the narrow MFT-only compatibility hook. It is
// invoked only by serial final acceptTracks(), after refit has completed.
template <>
class AcceptedTrackShadowPublisher<o2::mft::constants::mft::LayersNumber>
{
 public:
  std::optional<uint32_t> publish(TimeFrame& frame, const CommonTrackShadowRecord& record, const MFTCATrack& track) const
  {
    if (mSidecar == nullptr) {
      return std::nullopt;
    }
    MFTPublicationCompatibilityTransaction compatibility{*mSidecar, track.getTrack().getInvQPtSeed(), track.getTrack().getChi2QPtSeed()};
    return publishCommonTrackShadow(frame, record, compatibility, [](CommonTrackShadowPublishStep) {});
  }

  void adoptMFTPublicationCompatibility(MFTPublicationCompatibility* sidecar) noexcept { mSidecar = sidecar; }

 private:
  MFTPublicationCompatibility* mSidecar = nullptr;
};

} // namespace o2::itsmft::tracking

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_ACCEPTEDTRACKSHADOWPUBLISHER_H_

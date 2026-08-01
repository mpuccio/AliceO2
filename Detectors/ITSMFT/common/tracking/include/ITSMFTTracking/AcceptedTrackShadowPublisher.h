// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_ACCEPTEDTRACKSHADOWPUBLISHER_H_
#define ALICEO2_ITSMFT_TRACKING_ACCEPTEDTRACKSHADOWPUBLISHER_H_

#ifndef GPUCA_GPUCODE

#include <optional>

#include "ITSMFTTracking/CommonTrackShadow.h"
#include "ITSMFTTracking/ITSSharedClusterCompatibility.h"
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
  void adoptITSSharedClusterCompatibility(ITSSharedClusterCompatibility*) noexcept {}

  template <typename Tracks>
  bool sealITSSharedClusterCompatibility(const Tracks&) const noexcept
  {
    return true;
  }
};

// Narrow ITS-only compatibility hook. Pending associations are appended with
// the CommonTrack transaction at final acceptance; final shared status is
// sealed from the same pre-sort legacy slots after markTracks().
template <>
class AcceptedTrackShadowPublisher<ITSNLayers>
{
 public:
  template <typename Track>
  std::optional<uint32_t> publish(TimeFrame& frame, const CommonTrackShadowRecord& record, const Track&) const
  {
    if (mSidecar == nullptr) {
      // Standalone tracker/unit-test use without an interface-owned bridge
      // remains a pure CommonTrack shadow path. Production ITS interfaces
      // always adopt the bridge before tracking starts.
      return publishCommonTrackShadow(frame, record);
    }
    ITSSharedClusterCompatibilityTransaction compatibility{*mSidecar};
    return publishCommonTrackShadow(frame, record, compatibility, [](CommonTrackShadowPublishStep) {});
  }

  void adoptMFTPublicationCompatibility(MFTPublicationCompatibility*) noexcept {}
  void adoptITSSharedClusterCompatibility(ITSSharedClusterCompatibility* sidecar) noexcept { mSidecar = sidecar; }

  template <typename Tracks>
  bool sealITSSharedClusterCompatibility(const Tracks& tracks) const
  {
    return mSidecar == nullptr || mSidecar->sealFromMarkedTracks(tracks);
  }

 private:
  ITSSharedClusterCompatibility* mSidecar = nullptr;
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
    MFTPublicationCompatibilityTransaction compatibility{*mSidecar, track.getTrack().getInvQPtSeed(), track.getTrack().getChi2QPtSeed(),
                                                         track.getSeedPattern(), static_cast<float>(track.getTrack().getOutParam().getTrackChi2())};
    return publishCommonTrackShadow(frame, record, compatibility, [](CommonTrackShadowPublishStep) {});
  }

  void adoptMFTPublicationCompatibility(MFTPublicationCompatibility* sidecar) noexcept { mSidecar = sidecar; }
  void adoptITSSharedClusterCompatibility(ITSSharedClusterCompatibility*) noexcept {}

  template <typename Tracks>
  bool sealITSSharedClusterCompatibility(const Tracks&) const noexcept
  {
    return true;
  }

 private:
  MFTPublicationCompatibility* mSidecar = nullptr;
};

} // namespace o2::itsmft::tracking

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_ACCEPTEDTRACKSHADOWPUBLISHER_H_

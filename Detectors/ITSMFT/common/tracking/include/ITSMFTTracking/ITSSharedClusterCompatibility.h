// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_ITSSHAREDCLUSTERCOMPATIBILITY_H_
#define ALICEO2_ITSMFT_TRACKING_ITSSHAREDCLUSTERCOMPATIBILITY_H_

#ifndef GPUCA_GPUCODE

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace o2::itsmft::tracking
{

// Temporary ITS output-compatibility data. The pending sequence is an
// explicit, pre-sort bridge from a legacy accepted-track slot to its global
// CommonTrack index. It is never exposed to an output adapter. After the
// final serial markTracks() pass, entries() exposes only the sealed sparse
// global-index sequence carrying the final legacy shared-cluster bit.
struct ITSSharedClusterCompatibilityEntry {
  uint32_t commonTrackIndex{};
  bool hasSharedClusters{};
};

enum class ITSSharedClusterCompatibilitySealStep : uint8_t {
  BeforeReserve
};

class ITSSharedClusterCompatibility
{
 public:
  const std::vector<ITSSharedClusterCompatibilityEntry>& entries() const noexcept { return mEntries; }
  bool isSealed() const noexcept { return mSealed; }
  size_t pendingSize() const noexcept { return mPendingIndices.size(); }

  void clear() noexcept
  {
    mPendingIndices.clear();
    mEntries.clear();
    mSealed = false;
  }

  // Called only after the final serial markTracks() pass and before legacy
  // rectify/sort. `tracks` is deliberately supplied here, while the explicit
  // pending indices are still aligned with its un-reordered accepted slots.
  template <typename Tracks, typename Hook>
  bool sealFromMarkedTracks(const Tracks& tracks, Hook&& hook)
  {
    if (mSealed || tracks.size() != mPendingIndices.size()) {
      return false;
    }
    std::vector<ITSSharedClusterCompatibilityEntry> staged;
    hook(ITSSharedClusterCompatibilitySealStep::BeforeReserve);
    staged.reserve(mPendingIndices.size());
    uint32_t previous = 0;
    bool havePrevious = false;
    for (size_t i = 0; i < mPendingIndices.size(); ++i) {
      const auto index = mPendingIndices[i];
      if ((havePrevious && previous >= index)) {
        return false;
      }
      staged.push_back({index, tracks[i].hasSharedClusters()});
      previous = index;
      havePrevious = true;
    }
    mEntries.swap(staged);
    mSealed = true;
    return true;
  }

  template <typename Tracks>
  bool sealFromMarkedTracks(const Tracks& tracks)
  {
    return sealFromMarkedTracks(tracks, [](ITSSharedClusterCompatibilitySealStep) {});
  }

 private:
  friend class ITSSharedClusterCompatibilityTransaction;
  std::vector<uint32_t> mPendingIndices;
  std::vector<ITSSharedClusterCompatibilityEntry> mEntries;
  bool mSealed = false;
};

// Local participant in publishCommonTrackShadow(). It records the explicit
// pre-sort association in the same transaction as the new CommonTrack.
class ITSSharedClusterCompatibilityTransaction
{
 public:
  explicit ITSSharedClusterCompatibilityTransaction(ITSSharedClusterCompatibility& sidecar)
    : mSidecar{sidecar}, mOldSize{sidecar.mPendingIndices.size()}
  {
  }

  bool validate(uint32_t commonTrackIndex) const noexcept
  {
    return !mSidecar.mSealed && (mSidecar.mPendingIndices.empty() || mSidecar.mPendingIndices.back() < commonTrackIndex);
  }
  void reserve() { mSidecar.mPendingIndices.reserve(mOldSize + 1); }
  void append(uint32_t commonTrackIndex)
  {
    if (!validate(commonTrackIndex)) {
      throw std::logic_error{"invalid ITS shared-cluster compatibility index"};
    }
    mSidecar.mPendingIndices.push_back(commonTrackIndex);
  }
  void rollback() noexcept { mSidecar.mPendingIndices.resize(mOldSize); }

 private:
  ITSSharedClusterCompatibility& mSidecar;
  size_t mOldSize;
};

} // namespace o2::itsmft::tracking

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_ITSSHAREDCLUSTERCOMPATIBILITY_H_

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_MFTPUBLICATIONCOMPATIBILITY_H_
#define ALICEO2_ITSMFT_TRACKING_MFTPUBLICATIONCOMPATIBILITY_H_

#ifndef GPUCA_GPUCODE

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace o2::itsmft::tracking
{

// MFT output-compatibility sidecar kept outside CommonTrack and TimeFrame.
// Indices address the shared CommonTrack collection, not an MFT-local slot.
struct MFTPublicationCompatibilityEntry {
  uint32_t commonTrackIndex{};
  double invQPtSeed{};
  double chi2QPtSeed{};
  float outParamChi2{};
  uint16_t seedPattern{};
};

class MFTPublicationCompatibility
{
 public:
  const std::vector<MFTPublicationCompatibilityEntry>& entries() const noexcept { return mEntries; }
  void clear() noexcept { mEntries.clear(); }

  // Materialize compatibility fields from detector-neutral accepted results.
  template <typename Results>
  bool replaceFromAcceptedResults(const Results& results)
  {
    std::vector<MFTPublicationCompatibilityEntry> staged;
    staged.reserve(results.size());
    uint32_t previous = 0;
    bool havePrevious = false;
    for (const auto& result : results) {
      const auto index = result.commonTrackIndex;
      if ((havePrevious && previous >= index) || index == std::numeric_limits<uint32_t>::max()) {
        return false;
      }
      staged.push_back({index,
                        static_cast<double>(result.seed.getQOverPt()),
                        0.,
                        0.f,
                        result.seed.getHitLayerMask().value()});
      previous = index;
      havePrevious = true;
    }
    mEntries.swap(staged);
    return true;
  }

  const MFTPublicationCompatibilityEntry* find(uint32_t commonTrackIndex, size_t commonTrackCount) const noexcept
  {
    if (commonTrackIndex >= commonTrackCount) {
      return nullptr;
    }
    const auto it = std::lower_bound(mEntries.begin(), mEntries.end(), commonTrackIndex,
                                     [](const auto& entry, uint32_t index) { return entry.commonTrackIndex < index; });
    return it != mEntries.end() && it->commonTrackIndex == commonTrackIndex ? &*it : nullptr;
  }

 private:
  friend class MFTPublicationCompatibilityTransaction;
  std::vector<MFTPublicationCompatibilityEntry> mEntries;
};

// Transactional participant for one compatibility entry. Invalid ordering
// fails before append, and rollback restores the prior sequence size.
class MFTPublicationCompatibilityTransaction
{
 public:
  MFTPublicationCompatibilityTransaction(MFTPublicationCompatibility& sidecar, double invQPtSeed, double chi2QPtSeed, uint16_t seedPattern,
                                         float outParamChi2 = 0.f)
    : mSidecar{sidecar}, mEntry{0u, invQPtSeed, chi2QPtSeed, outParamChi2, seedPattern}, mOldSize{sidecar.mEntries.size()}
  {
  }

  bool validate(uint32_t commonTrackIndex) noexcept
  {
    mEntry.commonTrackIndex = commonTrackIndex;
    return mSidecar.mEntries.empty() || mSidecar.mEntries.back().commonTrackIndex < commonTrackIndex;
  }

  void reserve() { mSidecar.mEntries.reserve(mOldSize + 1); }
  void append(uint32_t commonTrackIndex)
  {
    if (commonTrackIndex != mEntry.commonTrackIndex ||
        (!mSidecar.mEntries.empty() && mSidecar.mEntries.back().commonTrackIndex >= commonTrackIndex)) {
      throw std::logic_error{"invalid MFT publication compatibility index"};
    }
    mSidecar.mEntries.push_back(mEntry);
  }
  void rollback() noexcept { mSidecar.mEntries.resize(mOldSize); }

 private:
  MFTPublicationCompatibility& mSidecar;
  MFTPublicationCompatibilityEntry mEntry;
  size_t mOldSize;
};

} // namespace o2::itsmft::tracking

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_MFTPUBLICATIONCOMPATIBILITY_H_

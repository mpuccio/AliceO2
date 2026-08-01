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
#include <stdexcept>
#include <vector>

namespace o2::itsmft::tracking
{

// Temporary MFT output-compatibility data. This deliberately remains outside
// CommonTrack and TimeFrame: it is only needed until the TrackMFT publication
// schema can deliberately evolve. Indices address the one shared TimeFrame
// CommonTrack collection, not an MFT-local position, so ITS/MFT entries may
// interleave in a later combined owner.
struct MFTPublicationCompatibilityEntry {
  uint32_t commonTrackIndex{};
  double invQPtSeed{};
  double chi2QPtSeed{};
};

class MFTPublicationCompatibility
{
 public:
  const std::vector<MFTPublicationCompatibilityEntry>& entries() const noexcept { return mEntries; }
  void clear() noexcept { mEntries.clear(); }

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

// Local participant in publishCommonTrackShadow(). The caller stages one
// entry before touching any owner; duplicate/non-monotonic keys fail before
// reserve, and rollback restores the original sparse sequence size.
class MFTPublicationCompatibilityTransaction
{
 public:
  MFTPublicationCompatibilityTransaction(MFTPublicationCompatibility& sidecar, double invQPtSeed, double chi2QPtSeed)
    : mSidecar{sidecar}, mEntry{0u, invQPtSeed, chi2QPtSeed}, mOldSize{sidecar.mEntries.size()}
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

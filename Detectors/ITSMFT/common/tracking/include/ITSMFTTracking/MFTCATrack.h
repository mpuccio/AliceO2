// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file MFTCATrack.h
/// \brief MFT CA track with native forward parameters and layer-indexed cluster refs
///

#ifndef ALICEO2_ITSMFT_TRACKING_MFTCATRACK_H_
#define ALICEO2_ITSMFT_TRACKING_MFTCATRACK_H_

#include <array>
#include <cstdint>

#include "DataFormatsITS/TimeEstBC.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "ITSMFTTracking/Cell.h"
#include "ITStracking/Constants.h"
#include "MFTTracking/Constants.h"

namespace o2::itsmft::tracking
{

/// MFT CA output track: forward fit parameters plus ITS-style cluster indexing for CA bookkeeping.
class MFTCATrack
{
  enum UserBits {
    kSharedClusters = 1 << 28
  };

 public:
  static constexpr int MaxClusters = o2::its::TrackITSExt::MaxClusters;

  MFTCATrack()
  {
    mIndex.fill(o2::its::constants::UnusedIndex);
  }

  const o2::mft::TrackMFT& getTrack() const { return mTrack; }
  o2::mft::TrackMFT& getTrack() { return mTrack; }

  float getChi2() const { return static_cast<float>(mTrack.getTrackChi2()); }
  float getPhi() const { return static_cast<float>(mTrack.getPhi()); }
  float getEta() const { return static_cast<float>(mTrack.getEta()); }
  double getCharge() const { return mTrack.getCharge(); }

  int getNumberOfClusters() const
  {
    int n{0};
    for (int layer = 0; layer < MaxClusters; ++layer) {
      if (hasHitOnLayer(layer)) {
        ++n;
      }
    }
    return n;
  }

  uint32_t getPattern() const { return mPattern; }
  void setPattern(uint32_t p) { mPattern = p; }
  bool hasHitOnLayer(uint32_t layer) const { return mPattern & (0x1u << layer); }

  uint32_t getFirstClusterLayer() const
  {
    uint32_t s{0};
    while (!(mPattern & (1u << s)) && s < MaxClusters) {
      ++s;
    }
    return s;
  }

  uint32_t getLastClusterLayer() const
  {
    uint32_t r{0}, v{mPattern & ((1u << MaxClusters) - 1)};
    while (v >>= 1) {
      ++r;
    }
    return r;
  }

  int getClusterIndex(int layer) const { return mIndex[layer]; }

  int getFirstLayerClusterIndex() const { return getClusterIndex(getFirstClusterLayer()); }

  void setClusterIndex(int layer, int idx)
  {
    mIndex[layer] = idx;
    if (idx != o2::its::constants::UnusedIndex) {
      mPattern |= 0x1u << layer;
    }
  }

  void setExternalClusterIndex(int layer, int idx, bool newCluster = false)
  {
    if (newCluster) {
      mPattern |= 0x1u << layer;
    }
    mIndex[layer] = idx;
  }

  void setClusterSize(int layer, int size)
  {
    if (layer >= 10) {
      return;
    }
    if (size > 63) {
      size = 63;
    }
    mClusterSizes &= ~(0x3fULL << (layer * 6));
    mClusterSizes |= (static_cast<uint64_t>(size) << (layer * 6));
  }

  int getClusterSize(int layer) const
  {
    if (layer >= 10) {
      return 0;
    }
    return (mClusterSizes >> (layer * 6)) & 0x3f;
  }

  auto& getTimeStamp() { return mTime; }
  const auto& getTimeStamp() const { return mTime; }

  uint16_t getSeedPattern() const { return mSeedPattern; }
  void setSeedPattern(uint16_t pattern) { mSeedPattern = pattern; }

  void setSharedClusters(bool toggle = true)
  {
    mClusterSizes = toggle ? (mClusterSizes | kSharedClusters) : (mClusterSizes & ~kSharedClusters);
  }

  bool hasSharedClusters() const { return mClusterSizes & kSharedClusters; }

 private:
  o2::mft::TrackMFT mTrack;
  std::array<int, MaxClusters> mIndex = {};
  uint32_t mPattern = 0;
  o2::its::TimeStamp mTime;
  uint16_t mSeedPattern{0};
  uint64_t mClusterSizes = 0;
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_MFTCATRACK_H_ */

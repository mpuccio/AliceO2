// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_DETECTORPUBLICATIONADAPTER_H_
#define ALICEO2_ITSMFT_TRACKING_DETECTORPUBLICATIONADAPTER_H_

#ifndef GPUCA_GPUCODE

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <limits>
#include <vector>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/detail/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/detail/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/TrackingOperationAdapter.h"
#include "ITStracking/MathUtils.h"
#include "MFTTracking/Constants.h"

namespace o2::itsmft::tracking
{

// The generic tracker publishes CommonTrack results. This adapter owns only
// detector compatibility sidecars; typed accepted-track vectors stay outside.
template <int NLayers>
class DetectorPublicationAdapter
{
 public:
  void adoptMFTPublicationCompatibility(MFTPublicationCompatibility*) noexcept {}
  void adoptITSSharedClusterCompatibility(ITSSharedClusterCompatibility*) noexcept {}
  ITSSharedClusterCompatibility* getITSSharedClusterCompatibility() const noexcept { return nullptr; }
  MFTPublicationCompatibility* getMFTPublicationCompatibility() const noexcept { return nullptr; }

  bool completeAccepted(gsl::span<const TrackingCandidate>,
                        const TrackingParameters&,
                        const SurfaceTrackingScratch&,
                        bool) const noexcept
  {
    return true;
  }
  void reset() noexcept {}
};

template <>
class DetectorPublicationAdapter<ITSNLayers>
{
 public:
  void adoptMFTPublicationCompatibility(MFTPublicationCompatibility*) noexcept {}
  void adoptITSSharedClusterCompatibility(ITSSharedClusterCompatibility* sidecar) noexcept { mSidecar = sidecar; }
  ITSSharedClusterCompatibility* getITSSharedClusterCompatibility() const noexcept { return mSidecar; }
  MFTPublicationCompatibility* getMFTPublicationCompatibility() const noexcept { return nullptr; }

  bool completeAccepted(gsl::span<const TrackingCandidate> candidates,
                        const TrackingParameters& params,
                        const SurfaceTrackingScratch& scratch,
                        bool final)
  {
    if (mSidecar == nullptr) {
      return true;
    }
    if (!stageSharedClusterFlags(candidates, params, scratch)) {
      return false;
    }
    return !final || mSidecar->replaceFromAcceptedResults(candidates, mSharedClusterFlags);
  }

  void reset() noexcept
  {
    mSharedClusterFlags.clear();
    if (mSidecar != nullptr) {
      mSidecar->clear();
    }
  }

 private:
  bool stageSharedClusterFlags(gsl::span<const TrackingCandidate> candidates,
                               const TrackingParameters& params,
                               const SurfaceTrackingScratch& scratch)
  {
    uint32_t maxIndex = 0;
    for (const auto& candidate : candidates) {
      if (candidate.commonTrackIndex == std::numeric_limits<uint32_t>::max()) {
        return false;
      }
      maxIndex = std::max(maxIndex, candidate.commonTrackIndex);
    }
    if (mSharedClusterFlags.size() <= maxIndex) {
      mSharedClusterFlags.resize(static_cast<size_t>(maxIndex) + 1, 0);
    }
    if (!params.AllowSharingFirstCluster) {
      return true;
    }
    for (size_t first = 0; first < candidates.size(); ++first) {
      const auto& firstCandidate = candidates[first];
      const int firstLayer = firstCandidate.getFirstClusterLayer();
      if (firstLayer < 0) {
        continue;
      }
      const int firstCluster = firstCandidate.getClusterIndex(firstLayer);
      for (size_t second = first + 1; second < candidates.size(); ++second) {
        const auto& secondCandidate = candidates[second];
        const int secondLayer = secondCandidate.getFirstClusterLayer();
        if (secondLayer != firstLayer || secondCandidate.getClusterIndex(secondLayer) != firstCluster) {
          continue;
        }
        if (scratch.getClusterROF(firstLayer, firstCluster) != scratch.getClusterROF(secondLayer, secondCandidate.getClusterIndex(secondLayer))) {
          continue;
        }
        if (!o2::its::math_utils::isPhiDifferenceBelow(firstCandidate.phi, secondCandidate.phi, params.SharedClusterMaxDeltaPhi)) {
          continue;
        }
        if (std::abs(firstCandidate.eta - secondCandidate.eta) > params.SharedClusterMaxDeltaEta) {
          continue;
        }
        if (params.SharedClusterOppositeSign && firstCandidate.charge == secondCandidate.charge) {
          continue;
        }
        mSharedClusterFlags[firstCandidate.commonTrackIndex] = 1;
        mSharedClusterFlags[secondCandidate.commonTrackIndex] = 1;
      }
    }
    return true;
  }

  ITSSharedClusterCompatibility* mSidecar = nullptr;
  std::vector<uint8_t> mSharedClusterFlags;
};

template <>
class DetectorPublicationAdapter<o2::mft::constants::mft::LayersNumber>
{
 public:
  void adoptMFTPublicationCompatibility(MFTPublicationCompatibility* sidecar) noexcept { mSidecar = sidecar; }
  void adoptITSSharedClusterCompatibility(ITSSharedClusterCompatibility*) noexcept {}
  ITSSharedClusterCompatibility* getITSSharedClusterCompatibility() const noexcept { return nullptr; }
  MFTPublicationCompatibility* getMFTPublicationCompatibility() const noexcept { return mSidecar; }

  bool completeAccepted(gsl::span<const TrackingCandidate> candidates,
                        const TrackingParameters&,
                        const SurfaceTrackingScratch&,
                        bool) const
  {
    if (mSidecar == nullptr) {
      return true;
    }
    // Preserve initialized/default compatibility fields and the undefined
    // invQPtSeed writer byte required by the legacy output contract.
    return mSidecar->replaceFromAcceptedResults(candidates);
  }

  void reset() noexcept
  {
    if (mSidecar != nullptr) {
      mSidecar->clear();
    }
  }

 private:
  MFTPublicationCompatibility* mSidecar = nullptr;
};

} // namespace o2::itsmft::tracking

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_DETECTORPUBLICATIONADAPTER_H_

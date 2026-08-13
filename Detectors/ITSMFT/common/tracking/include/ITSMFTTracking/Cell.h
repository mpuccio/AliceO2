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
/// \file Cell.h
/// \brief CA cell/track seed types with hole-layer support (ITS PR #15390)
///

#ifndef ALICEO2_ITSMFT_TRACKING_INCLUDE_CACELL_H_
#define ALICEO2_ITSMFT_TRACKING_INCLUDE_CACELL_H_

#include <array>
#include <cstdint>

#include "DataFormatsITS/TimeEstBC.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IdTypes.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMask.h"
#include "ITSMFTTracking/TripletFitting.h"
#include "ITStracking/Constants.h"
#include "GPUCommonDef.h"

namespace o2::itsmft::tracking
{

struct CellNeighbour {
  int cellTopology{-1};
  int cell{-1};
  int nextCellTopology{-1};
  int nextCell{-1};
  int level{-1};
};

struct CellClusterReference {
  int surfacePosition{o2::its::constants::UnusedIndex};
  int clusterIndex{o2::its::constants::UnusedIndex};
};

/// Common non-`SurfaceKind`-templated CA cell representation.
class CellSeed final
{
 public:
  GPUhdDefault() CellSeed() = default;
  GPUhd() CellSeed(int innerL, int cl0, int cl1, int cl2, int trkl0, int trkl1, const SurfaceKinematicState& state, float chi2, const o2::its::TimeEstBC& time)
    : CellSeed(LayerMask(innerL, innerL + 1, innerL + 2), cl0, cl1, cl2, trkl0, trkl1, state, chi2, time)
  {
  }
  GPUhd() CellSeed(LayerMask hitLayerMask, int cl0, int cl1, int cl2, int trkl0, int trkl1, const SurfaceKinematicState& state, float chi2, const o2::its::TimeEstBC& time)
    : mState(state), mChi2(chi2), mLevel(1), mTime(time)
  {
    setHitLayerMask(hitLayerMask);
    auto& clusters = mClusters;
    clusters[0] = cl0;
    clusters[1] = cl1;
    clusters[2] = cl2;
    setFirstTrackletIndex(trkl0);
    setSecondTrackletIndex(trkl1);
  }
  GPUhdDefault() CellSeed(const CellSeed&) = default;
  GPUhdDefault() ~CellSeed() = default;
  GPUhdDefault() CellSeed(CellSeed&&) = default;
  GPUhdDefault() CellSeed& operator=(const CellSeed&) = default;
  GPUhdDefault() CellSeed& operator=(CellSeed&&) = default;

  GPUhd() LayerMask getHitLayerMask() const { return LayerMask{mHitLayerMask}; }
  GPUhd() void setHitLayerMask(LayerMask mask) { mHitLayerMask = mask.value(); }
  GPUhd() int getInnerLayer() const { return getHitLayerMask().first(); }
  GPUhd() int getFirstTrackletIndex() const { return mTracklets[0]; }
  GPUhd() void setFirstTrackletIndex(int trkl) { mTracklets[0] = trkl; }
  GPUhd() int getSecondTrackletIndex() const { return mTracklets[1]; }
  GPUhd() void setSecondTrackletIndex(int trkl) { mTracklets[1] = trkl; }
  GPUhd() float getChi2() const { return mChi2; }
  GPUhd() void setChi2(float chi2) { mChi2 = chi2; }
  GPUhd() int getLevel() const { return mLevel; }
  GPUhd() void setLevel(int level) { mLevel = level; }
  GPUhd() int* getLevelPtr() { return &mLevel; }
  GPUhd() auto& getTimeStamp() noexcept { return mTime; }
  GPUhd() const auto& getTimeStamp() const noexcept { return mTime; }
  GPUhd() SurfaceKinematicState& state() noexcept { return mState; }
  GPUhd() const SurfaceKinematicState& state() const noexcept { return mState; }
  GPUhd() float getQOverPt() const noexcept { return mState.parameters[4]; }
  GPUhd() int getFirstClusterIndex() const { return mClusters[0]; }
  GPUhd() int getSecondClusterIndex() const { return mClusters[1]; }
  GPUhd() int getThirdClusterIndex() const { return mClusters[2]; }
  GPUhd() auto& getClusters() { return mClusters; }
  GPUhd() const auto& getClusters() const { return mClusters; }
  GPUhd() TripletFitFactor& tripletFactor() noexcept { return mTripletFactor; }
  GPUhd() const TripletFitFactor& tripletFactor() const noexcept { return mTripletFactor; }
  GPUhd() CellClusterReference getClusterReference(int requestedSlot) const noexcept
  {
    if (requestedSlot < 0 || requestedSlot >= o2::its::constants::ClustersPerCell) {
      return {};
    }
    const auto mask = getHitLayerMask();
    int slot = 0;
    for (int position = 0; position < 32; ++position) {
      if (mask.has(position) && slot++ == requestedSlot) {
        return {position, mClusters[requestedSlot]};
      }
    }
    return {};
  }
  GPUhd() int getCluster(int layer) const
  {
    const int slot = getHitLayerMask().slot(layer);
    return (slot >= 0 && slot < o2::its::constants::ClustersPerCell) ? mClusters[slot] : o2::its::constants::UnusedIndex;
  }

 private:
  SurfaceKinematicState mState{};
  uint32_t mHitLayerMask{0};
  float mChi2{o2::its::constants::UnsetValue};
  int mLevel{o2::its::constants::UnusedIndex};
  std::array<int, 2> mTracklets = o2::its::constants::helpers::initArray<int, 2, o2::its::constants::UnusedIndex>();
  std::array<int, o2::its::constants::ClustersPerCell> mClusters =
    o2::its::constants::helpers::initArray<int, o2::its::constants::ClustersPerCell, o2::its::constants::UnusedIndex>();
  o2::its::TimeEstBC mTime;
  TripletFitFactor mTripletFactor{};
};

static_assert(std::is_trivially_copyable_v<CellSeed>);

/// GPU-portable, non-templated whole-track seed with one cluster slot per
/// adopted-plan position. Fixed MaxLayoutSurfaces capacity is required for
/// device use, where heap allocation is unavailable.
///
/// TrackSeed uses SurfaceMask rather than CellSeed's positional LayerMask:
/// each set bit is a position in its fixed array, not a global SurfaceId.
///
/// This fixed-capacity value is the sole common-CA whole-track seed
/// representation.
class TrackSeed final
{
 public:
  static constexpr int MaxSurfaces = static_cast<int>(MaxLayoutSurfaces);

  GPUhdDefault() TrackSeed() = default;
  GPUhdDefault() TrackSeed(const TrackSeed&) = default;
  GPUhdDefault() ~TrackSeed() = default;
  GPUhdDefault() TrackSeed(TrackSeed&&) = default;
  GPUhdDefault() TrackSeed& operator=(const TrackSeed&) = default;
  GPUhdDefault() TrackSeed& operator=(TrackSeed&&) = default;

  // CellSeed's hit mask is positional in the same fixed-capacity domain.
  GPUhd() explicit TrackSeed(const CellSeed& cs)
    : mState(cs.state()), mChi2(cs.getChi2()), mLevel(cs.getLevel()), mTracklets{cs.getFirstTrackletIndex(), cs.getSecondTrackletIndex()}, mTime(cs.getTimeStamp())
  {
    const auto hitMask = cs.getHitLayerMask();
    int slot = 0;
    for (int position = 0; position < MaxSurfaces; ++position) {
      if (hitMask.has(position)) {
        mClusters[position] = cs.getClusters()[slot++];
        mSurfaceMask.set(SurfaceId{static_cast<uint16_t>(position)});
      }
    }
  }

  GPUhd() SurfaceMask getSurfaceMask() const noexcept { return mSurfaceMask; }
  GPUhd() void setSurfaceMask(SurfaceMask mask) noexcept { mSurfaceMask = mask; }
  GPUhd() int getActiveSurfaceCount() const noexcept { return mSurfaceMask.count(); }
  GPUhd() int getInnerLayer() const noexcept { return mSurfaceMask.first(); }
  GPUhd() bool hasCluster(int position) const noexcept
  {
    return position >= 0 && position < MaxSurfaces && mSurfaceMask.has(SurfaceId{static_cast<uint16_t>(position)});
  }

  // Bounds-checked: an out-of-[0, MaxSurfaces) position safely
  // returns UnusedIndex instead of indexing out of bounds.
  GPUhd() int getCluster(int position) const noexcept
  {
    return (position >= 0 && position < MaxSurfaces) ? mClusters[position] : o2::its::constants::UnusedIndex;
  }

  // Keep a LayerMask view for the positional hole/acceptance code while
  // SurfaceMask remains the authoritative fixed-capacity representation.
  GPUhd() LayerMask getHitLayerMask() const noexcept
  {
    return LayerMask{mSurfaceMask.value()};
  }
  GPUhd() void setCluster(int position, int clusterIndex) noexcept
  {
    if (position >= 0 && position < MaxSurfaces) {
      mClusters[position] = clusterIndex;
    }
  }

  GPUhd() int getFirstClusterIndex() const noexcept { return getClusterBySlot(0); }
  GPUhd() int getSecondClusterIndex() const noexcept { return getClusterBySlot(1); }
  GPUhd() int getThirdClusterIndex() const noexcept { return getClusterBySlot(2); }

  GPUhd() auto& getClusters() noexcept { return mClusters; }
  GPUhd() const auto& getClusters() const noexcept { return mClusters; }

  GPUhd() int getFirstTrackletIndex() const noexcept { return mTracklets[0]; }
  GPUhd() void setFirstTrackletIndex(int trkl) noexcept { mTracklets[0] = trkl; }
  GPUhd() int getSecondTrackletIndex() const noexcept { return mTracklets[1]; }
  GPUhd() void setSecondTrackletIndex(int trkl) noexcept { mTracklets[1] = trkl; }

  GPUhd() float getChi2() const noexcept { return mChi2; }
  GPUhd() void setChi2(float chi2) noexcept { mChi2 = chi2; }
  GPUhd() int getLevel() const noexcept { return mLevel; }
  GPUhd() void setLevel(int level) noexcept { mLevel = level; }

  GPUhd() auto& getTimeStamp() noexcept { return mTime; }
  GPUhd() const auto& getTimeStamp() const noexcept { return mTime; }

  GPUhd() SurfaceKinematicState& state() noexcept { return mState; }
  GPUhd() const SurfaceKinematicState& state() const noexcept { return mState; }
  // Raw signed q/pT in slot 4 for cylinder and disk states; never squared.
  GPUhd() float getQOverPt() const noexcept { return mState.parameters[4]; }

 private:
  GPUhd() int getClusterBySlot(int requestedSlot) const noexcept
  {
    int slot = 0;
    for (int position = 0; position < MaxSurfaces; ++position) {
      if (hasCluster(position)) {
        if (slot++ == requestedSlot) {
          return mClusters[position];
        }
      }
    }
    return o2::its::constants::UnusedIndex;
  }

  SurfaceKinematicState mState{};
  SurfaceMask mSurfaceMask{};
  float mChi2{o2::its::constants::UnsetValue};
  int mLevel{o2::its::constants::UnusedIndex};
  std::array<int, 2> mTracklets = o2::its::constants::helpers::initArray<int, 2, o2::its::constants::UnusedIndex>();
  std::array<int, MaxSurfaces> mClusters = o2::its::constants::helpers::initArray<int, MaxSurfaces, o2::its::constants::UnusedIndex>();
  o2::its::TimeEstBC mTime;
};

// TrackSeed crosses the host/device boundary by value. TimeEstBC prevents a
// standard-layout assertion; trivially copyable is the required property.
static_assert(std::is_trivially_copyable_v<TrackSeed>);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_INCLUDE_CACELL_H_ */

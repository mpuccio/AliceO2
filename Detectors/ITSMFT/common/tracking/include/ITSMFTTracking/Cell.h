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
#include "ITSMFTTracking/LayerMask.h"
#include "ITSMFTTracking/SurfaceId.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMask.h"
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

/// Stage-B activation (Architecture.md Sec 10/11): shared Cell/TrackSeed
/// metadata, composing a named SurfaceKinematicState rather than inheriting a
/// detector-selected track parametrization. Neither TrackParCovF nor
/// TrackParCovFwd is a base of this type or of anything derived from it.
template <int NClusters>
class SeedMetadataBase
{
 public:
  GPUhd() LayerMask getHitLayerMask() const { return LayerMask{mHitLayerMask}; }
  GPUhd() void setHitLayerMask(LayerMask mask) { mHitLayerMask = mask.value(); }
  GPUhd() int getInnerLayer() const { return getHitLayerMask().first(); }
  GPUhd() int getFirstTrackletIndex() const { return mTracklets[0]; };
  GPUhd() void setFirstTrackletIndex(int trkl) { mTracklets[0] = trkl; };
  GPUhd() int getSecondTrackletIndex() const { return mTracklets[1]; };
  GPUhd() void setSecondTrackletIndex(int trkl) { mTracklets[1] = trkl; };
  GPUhd() float getChi2() const { return mChi2; };
  GPUhd() void setChi2(float chi2) { mChi2 = chi2; };
  GPUhd() int getLevel() const { return mLevel; };
  GPUhd() void setLevel(int level) { mLevel = level; };
  GPUhd() int* getLevelPtr() { return &mLevel; }
  GPUhd() auto& getTimeStamp() noexcept { return mTime; }
  GPUhd() const auto& getTimeStamp() const noexcept { return mTime; }

  GPUhd() SurfaceKinematicState& state() noexcept { return mState; }
  GPUhd() const SurfaceKinematicState& state() const noexcept { return mState; }

  /// Raw signed q/pT, common to both families (no barrel/forward branch):
  /// slot 4 is the raw signed q/pT parameter for both Barrel and Forward
  /// Stage-B state conventions -- never squared. The former common
  /// getQ2Pt() accessor squared this value (correct convention for neither
  /// family) and has been removed; do not reintroduce it.
  GPUhd() float getQOverPt() const noexcept { return mState.parameters[4]; }

 protected:
  GPUhdDefault() SeedMetadataBase() = default;
  GPUhdDefault() SeedMetadataBase(const SeedMetadataBase&) = default;
  GPUhdDefault() ~SeedMetadataBase() = default;
  GPUhdDefault() SeedMetadataBase(SeedMetadataBase&&) = default;
  GPUhdDefault() SeedMetadataBase& operator=(const SeedMetadataBase&) = default;
  GPUhdDefault() SeedMetadataBase& operator=(SeedMetadataBase&&) = default;
  GPUhd() SeedMetadataBase(const SurfaceKinematicState& state, float chi2, int level, const o2::its::TimeEstBC& time)
    : mState(state), mChi2(chi2), mLevel(level), mTime(time)
  {
  }
  GPUhd() auto& clustersRaw() { return mClusters; }
  GPUhd() const auto& clustersRaw() const { return mClusters; }

 private:
  SurfaceKinematicState mState{};
  uint16_t mHitLayerMask{0};
  float mChi2{o2::its::constants::UnsetValue};
  int mLevel{o2::its::constants::UnusedIndex};
  std::array<int, 2> mTracklets = o2::its::constants::helpers::initArray<int, 2, o2::its::constants::UnusedIndex>();
  std::array<int, NClusters> mClusters = o2::its::constants::helpers::initArray<int, NClusters, o2::its::constants::UnusedIndex>();
  o2::its::TimeEstBC mTime;
};

/// Common (non-family-templated) CA cell representation.
class CellSeed final : public SeedMetadataBase<o2::its::constants::ClustersPerCell>
{
  using Base = SeedMetadataBase<o2::its::constants::ClustersPerCell>;

 public:
  GPUhdDefault() CellSeed() = default;
  GPUhd() CellSeed(int innerL, int cl0, int cl1, int cl2, int trkl0, int trkl1, const SurfaceKinematicState& state, float chi2, const o2::its::TimeEstBC& time)
    : CellSeed(LayerMask(innerL, innerL + 1, innerL + 2), cl0, cl1, cl2, trkl0, trkl1, state, chi2, time)
  {
  }
  GPUhd() CellSeed(LayerMask hitLayerMask, int cl0, int cl1, int cl2, int trkl0, int trkl1, const SurfaceKinematicState& state, float chi2, const o2::its::TimeEstBC& time)
    : Base(state, chi2, 1, time)
  {
    this->setHitLayerMask(hitLayerMask);
    auto& clusters = this->clustersRaw();
    clusters[0] = cl0;
    clusters[1] = cl1;
    clusters[2] = cl2;
    this->setFirstTrackletIndex(trkl0);
    this->setSecondTrackletIndex(trkl1);
  }
  GPUhdDefault() CellSeed(const CellSeed&) = default;
  GPUhdDefault() ~CellSeed() = default;
  GPUhdDefault() CellSeed(CellSeed&&) = default;
  GPUhdDefault() CellSeed& operator=(const CellSeed&) = default;
  GPUhdDefault() CellSeed& operator=(CellSeed&&) = default;

  GPUhd() int getFirstClusterIndex() const { return this->clustersRaw()[0]; };
  GPUhd() int getSecondClusterIndex() const { return this->clustersRaw()[1]; };
  GPUhd() int getThirdClusterIndex() const { return this->clustersRaw()[2]; };
  GPUhd() auto& getClusters() { return this->clustersRaw(); }
  GPUhd() const auto& getClusters() const { return this->clustersRaw(); }
  GPUhd() int getCluster(int layer) const
  {
    const int slot = this->getHitLayerMask().slot(layer);
    return (slot >= 0 && slot < o2::its::constants::ClustersPerCell) ? this->clustersRaw()[slot] : o2::its::constants::UnusedIndex;
  }
};

/// M6c (doc/design/0002-m6-generic-workspace-migration.md Sec 4.2, 9):
/// GPU-portable, non-templated whole-track seed -- the generic successor to
/// the former layer-count-specific seed. Its per-surface cluster-index array is indexed
/// positionally, one slot per position in the adopted plan's ownedSurfaces()
/// order, exactly like the former seed's per-layer indexing -- but
/// fixed at MaxLayoutSurfaces (SurfaceId.h) capacity rather than templated
/// on NLayers, since this type must remain usable on device (GPUhd()),
/// where heap allocation is unavailable. MaxLayoutSurfaces already bounds
/// every owned-surface position in this library (SurfaceMask,
/// SurfaceGraphBuilder, so no graph this library
/// can validly build can ever exceed this capacity -- no new bound is
/// invented here.
///
/// Deliberately does not derive from SeedMetadataBase<N> (unlike CellSeed):
/// SeedMetadataBase's own hit-mask field is a
/// 16-bit LayerMask, which cannot mark all MaxLayoutSurfaces=32 positions
/// active. Reusing it here would either silently truncate at 16 active
/// positions or require widening SeedMetadataBase itself -- a shared base
/// CellSeed still uses in production. TrackSeed instead duplicates
/// SeedMetadataBase's small metadata surface
/// directly and uses SurfaceMask (already 32-bit, and already reused
/// positionally elsewhere in this library -- see SurfaceMask.h's own
/// positionalSurfaceMask()) as its active-surface mask: each set bit is a
/// *position* in this seed's own fixed array, never a numeric comparison
/// against a real global SurfaceId.
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

  // CellSeed's own hit mask is a 16-bit LayerMask (SeedMetadataBase<N>'s
  // stored width, independent of N): every set bit already lies in
  // [0, 15], well within MaxSurfaces, so this conversion never needs to
  // know NLayers -- this fixed-capacity conversion uses the same positional
  // mapping as the former layer-count-specific constructor, generalized from
  // an NLayers-wide loop to a fixed 16-wide one (CellSeed's mask can never set
  // a bit beyond position 15 in the first place).
  GPUhd() explicit TrackSeed(const CellSeed& cs)
    : mState(cs.state()), mChi2(cs.getChi2()), mLevel(cs.getLevel()), mTracklets{cs.getFirstTrackletIndex(), cs.getSecondTrackletIndex()}, mTime(cs.getTimeStamp())
  {
    const auto hitMask = cs.getHitLayerMask();
    int slot = 0;
    for (int position = 0; position < 16; ++position) {
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

  // The active common-CA layer range is at most ten positions today. Keep a
  // LayerMask view for the hole/acceptance code while SurfaceMask remains the
  // authoritative fixed-capacity representation.
  GPUhd() LayerMask getHitLayerMask() const noexcept
  {
    return LayerMask{static_cast<uint16_t>(mSurfaceMask.value())};
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
  // See SeedMetadataBase<N>::getQOverPt()'s own doc: raw signed q/pT, never
  // squared, common to both Barrel and Forward Stage-B state conventions.
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

// TrackSeed is a GPU value type (GPUhd() throughout, exactly like CellSeed/
// CellSeed above), so the applicable property is copyability by
// value across the host/device boundary, not standard-layout: compiled and
// checked here, TrackSeed is *not* standard-layout, because its embedded
// o2::its::TimeEstBC (o2::dataformats::TimeStampWithError<uint32_t,
// uint16_t> deriving from TimeStamp<uint32_t>) has non-static data members
// declared in more than one class of its own hierarchy -- a property of
// TimeEstBC itself, unrelated to anything TrackSeed adds. This is not a
// TrackSeed-specific problem: CellSeed embeds the
// exact same mTime member via SeedMetadataBase and carry no
// standard-layout/trivially-copyable static_assert either. is_trivially_copyable
// is the property this codebase's own device-value-type convention actually
// needs (byte-for-byte copyable, no user-defined copy/move/destructor
// logic) and the one TrackSeed does satisfy, checked below.
static_assert(std::is_trivially_copyable_v<TrackSeed>);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_INCLUDE_CACELL_H_ */

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

#include "DetectorsCommonDataFormats/DetID.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITS/TimeEstBC.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/LayerMask.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/TransitionPolicy.h"
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

/// Compatibility alias: resolves to the same common CellSeed type regardless
/// of NLayers (kept only so existing ITSNLayers/MFTNLayers call sites do not
/// need to change to the bare name).
template <int NLayers>
using CellSeedN = CellSeed;

/// TrackSeed remains templated only by NLayers: its cluster-index array width
/// is a temporary legacy boundary (raw int indices, LayerMask), not a
/// detector-selected state family.
template <int NLayers>
class TrackSeedTpl final : public SeedMetadataBase<NLayers>
{
  using Base = SeedMetadataBase<NLayers>;

 public:
  GPUhdDefault() TrackSeedTpl() = default;
  GPUhd() TrackSeedTpl(const CellSeed& cs)
    : Base(cs.state(), cs.getChi2(), cs.getLevel(), cs.getTimeStamp())
  {
    this->setHitLayerMask(cs.getHitLayerMask());
    this->setFirstTrackletIndex(cs.getFirstTrackletIndex());
    this->setSecondTrackletIndex(cs.getSecondTrackletIndex());
    auto& clusters = this->clustersRaw();
    int slot = 0;
    const auto hitMask = cs.getHitLayerMask();
    for (int layer = 0; layer < NLayers; ++layer) {
      if (hitMask.has(layer)) {
        clusters[layer] = cs.getClusters()[slot++];
      }
    }
  }
  GPUhdDefault() TrackSeedTpl(const TrackSeedTpl&) = default;
  GPUhdDefault() ~TrackSeedTpl() = default;
  GPUhdDefault() TrackSeedTpl(TrackSeedTpl&&) = default;
  GPUhdDefault() TrackSeedTpl& operator=(const TrackSeedTpl&) = default;
  GPUhdDefault() TrackSeedTpl& operator=(TrackSeedTpl&&) = default;

  GPUhd() int getFirstClusterIndex() const { return getClusterBySlot(0); }
  GPUhd() int getSecondClusterIndex() const { return getClusterBySlot(1); }
  GPUhd() int getThirdClusterIndex() const { return getClusterBySlot(2); }
  GPUhd() auto& getClusters() { return this->clustersRaw(); }
  GPUhd() const auto& getClusters() const { return this->clustersRaw(); }
  GPUhd() int getCluster(int layer) const { return this->clustersRaw()[layer]; }

 private:
  GPUhd() int getClusterBySlot(int requestedSlot) const
  {
    int slot = 0;
    const auto hitMask = this->getHitLayerMask();
    for (int layer = 0; layer < NLayers; ++layer) {
      if (hitMask.has(layer)) {
        if (slot++ == requestedSlot) {
          return this->clustersRaw()[layer];
        }
      }
    }
    return o2::its::constants::UnusedIndex;
  }
};

/// Compatibility alias: TrackSeedN<NLayers> resolves to TrackSeedTpl<NLayers>,
/// differing between families only by cluster-array width.
template <int NLayers>
using TrackSeedN = TrackSeedTpl<NLayers>;

template <int NLayers>
struct CATrackTypeHelper {
  using type = o2::its::TrackITSExt;
};

template <int NLayers>
using CATrackType = typename CATrackTypeHelper<NLayers>::type;

/// Temporary NLayers -> StateFamily compatibility boundary (Architecture.md
/// Sec 10.1): the common Cell/TrackSeed representation no longer encodes
/// family in its C++ type, so orchestration boundaries that need to validate
/// "this TrackerTraits<NLayers> instantiation may process this
/// TransitionPolicyTraits<Tag>::Family" compare this against Traits::Family
/// directly, instead of comparing Cell/TrackSeed types. ITSNLayers maps to
/// Barrel, MFTNLayers (o2::mft::constants::mft::LayersNumber) maps to
/// Forward. This inference is temporary and must not be expanded into a new
/// durable detector abstraction.
template <int NLayers>
GPUhdi() constexpr StateFamily stateFamilyFromNLayers() noexcept
{
  return detIdFromNLayers<NLayers>() == o2::detectors::DetID::MFT ? StateFamily::Forward : StateFamily::Barrel;
}

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_INCLUDE_CACELL_H_ */

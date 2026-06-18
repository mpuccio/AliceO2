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
#include <type_traits>

#include "DetectorsCommonDataFormats/DetID.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITS/TimeEstBC.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/LayerMask.h"
#include "ITStracking/Constants.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/TrackFwd.h"
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

template <int NClusters, typename TrackParT>
class SeedBase : public TrackParT
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

  /// Road-length filter: barrel q/pT² for ITS, (invQPt)² for forward MFT seeds.
  GPUhd() float getQ2Pt() const
  {
    if constexpr (std::is_same_v<TrackParT, o2::track::TrackParCovFwd>) {
      const float invQPt = static_cast<float>(TrackParT::getInvQPt());
      return invQPt * invQPt;
    } else {
      return TrackParT::getQ2Pt();
    }
  }

 protected:
  GPUhdDefault() SeedBase() = default;
  GPUhdDefault() SeedBase(const SeedBase&) = default;
  GPUhdDefault() ~SeedBase() = default;
  GPUhdDefault() SeedBase(SeedBase&&) = default;
  GPUhdDefault() SeedBase& operator=(const SeedBase&) = default;
  GPUhdDefault() SeedBase& operator=(SeedBase&&) = default;
  GPUhd() SeedBase(const TrackParT& tpc, float chi2, int level, const o2::its::TimeEstBC& time)
    : TrackParT(tpc), mChi2(chi2), mLevel(level), mTime(time)
  {
  }
  GPUhd() auto& clustersRaw() { return mClusters; }
  GPUhd() const auto& clustersRaw() const { return mClusters; }

 private:
  uint16_t mHitLayerMask{0};
  float mChi2{o2::its::constants::UnsetValue};
  int mLevel{o2::its::constants::UnusedIndex};
  std::array<int, 2> mTracklets = o2::its::constants::helpers::initArray<int, 2, o2::its::constants::UnusedIndex>();
  std::array<int, NClusters> mClusters = o2::its::constants::helpers::initArray<int, NClusters, o2::its::constants::UnusedIndex>();
  o2::its::TimeEstBC mTime;
};

template <typename TrackParT = o2::track::TrackParCovF>
class CellSeedTpl final : public SeedBase<o2::its::constants::ClustersPerCell, TrackParT>
{
  using Base = SeedBase<o2::its::constants::ClustersPerCell, TrackParT>;

 public:
  GPUhdDefault() CellSeedTpl() = default;
  GPUhd() CellSeedTpl(int innerL, int cl0, int cl1, int cl2, int trkl0, int trkl1, const TrackParT& tpc, float chi2, const o2::its::TimeEstBC& time)
    : CellSeedTpl(LayerMask(innerL, innerL + 1, innerL + 2), cl0, cl1, cl2, trkl0, trkl1, tpc, chi2, time)
  {
  }
  GPUhd() CellSeedTpl(LayerMask hitLayerMask, int cl0, int cl1, int cl2, int trkl0, int trkl1, const TrackParT& tpc, float chi2, const o2::its::TimeEstBC& time)
    : Base(tpc, chi2, 1, time)
  {
    this->setHitLayerMask(hitLayerMask);
    auto& clusters = this->clustersRaw();
    clusters[0] = cl0;
    clusters[1] = cl1;
    clusters[2] = cl2;
    this->setFirstTrackletIndex(trkl0);
    this->setSecondTrackletIndex(trkl1);
  }
  GPUhdDefault() CellSeedTpl(const CellSeedTpl&) = default;
  GPUhdDefault() ~CellSeedTpl() = default;
  GPUhdDefault() CellSeedTpl(CellSeedTpl&&) = default;
  GPUhdDefault() CellSeedTpl& operator=(const CellSeedTpl&) = default;
  GPUhdDefault() CellSeedTpl& operator=(CellSeedTpl&&) = default;

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

/// ITS default: barrel track parameters in cells.
using CellSeed = CellSeedTpl<o2::track::TrackParCovF>;

template <int NLayers, typename TrackParT = o2::track::TrackParCovF>
class TrackSeedTpl final : public SeedBase<NLayers, TrackParT>
{
  using Base = SeedBase<NLayers, TrackParT>;

 public:
  GPUhdDefault() TrackSeedTpl() = default;
  GPUhd() TrackSeedTpl(const CellSeedTpl<TrackParT>& cs)
    : Base(static_cast<const TrackParT&>(cs), cs.getChi2(), cs.getLevel(), cs.getTimeStamp())
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

template <int NLayers>
using TrackSeed = TrackSeedTpl<NLayers, o2::track::TrackParCovF>;

template <int NLayers>
struct CATrackTypeHelper {
  using type = o2::its::TrackITSExt;
};

template <int NLayers>
using CATrackType = typename CATrackTypeHelper<NLayers>::type;

/// Per-detector track parametrization stored in CA cells and extended seeds.
template <int NLayers>
struct CASeedTrackPar {
  static constexpr o2::detectors::DetID::ID DetId = detIdFromNLayers<NLayers>();
  using type = std::conditional_t<DetId == o2::detectors::DetID::MFT, o2::track::TrackParCovFwd, o2::track::TrackParCovF>;
};

template <int NLayers>
using CellSeedN = CellSeedTpl<typename CASeedTrackPar<NLayers>::type>;

template <int NLayers>
using TrackSeedN = TrackSeedTpl<NLayers, typename CASeedTrackPar<NLayers>::type>;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_INCLUDE_CACELL_H_ */

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

#include <algorithm>
#include <format>
#include <limits>

#include "ITSMFTTracking/Configuration.h"

namespace o2::itsmft
{

std::string TrackingParameters::asString() const
{
  std::string str = std::format("NColB:{} NRowB:{} PerVtx:{} DropFail:{} TtklMinPt:{:.2f} MinCl:{}", ColBins, RowBins, PerPrimaryVertexProcessing, DropTFUponFailure, TrackletMinPt, MinTrackLength);
  auto isSet = [](auto e) { return e >= 0; };
  auto isAnySet = [&isSet](auto v) { return !v.empty() && std::any_of(v.begin(), v.end(), isSet); };
  bool first = true;
  for (int il = NLayers; il >= MinTrackLength; il--) {
    int slot = NLayers - il;
    if (slot < (int)MinPt.size() && MinPt[slot] > 0) {
      if (first) {
        first = false;
        str += " MinPt: ";
      }
      str += std::format("L{}:{:.2f} ", il, MinPt[slot]);
    }
  }
  if (isAnySet(SystError2Row) || isAnySet(SystError2Col)) {
    str += " SystErrRow/Col:";
    for (size_t i = 0; i < SystError2Row.size(); i++) {
      str += std::format("{:.2e}/{:.2e} ", SystError2Row[i], SystError2Col[i]);
    }
  }
  if (isAnySet(AddTimeError)) {
    str += " AddTimeError:";
    for (unsigned int i : AddTimeError) {
      str += std::format("{} ", i);
    }
  }
  if (SharedMaxClusters) {
    str += std::format(" ShaMaxCls:{} ", SharedMaxClusters);
  }
  if (AllowSharingFirstCluster) {
    str += std::format(" ShaClsDPhi:{} ShaClsDEta:{} ShaClsSign:{}", SharedClusterMaxDeltaPhi, SharedClusterMaxDeltaEta, SharedClusterOppositeSign);
  }
  if (MaxHoles) {
    str += std::format(" MaxHoles:{} HoleMask:{:016b}", MaxHoles, HoleLayerMask);
  }
  if (std::numeric_limits<size_t>::max() != MaxMemory) {
    str += std::format(" MemLimit {:.2f} GB", double(MaxMemory) / (1024.f * 1024.f * 1024.f));
  }
  return str;
}

std::string VertexingParameters::asString() const
{
  std::string str = std::format("NColB:{} NRowB:{} MinVtxCont:{} SupLowMultDebris:{} MaxTrkltCls:{} ZCut:{} PhCut:{} PairCut:{} ClCut:{} SeedRad:{}x{}",
                                ColBins, RowBins, clusterContributorsCut, suppressLowMultDebris, maxTrackletsPerCluster, zCut, phiCut, pairCut, clusterCut, seedMemberRadiusTime, seedMemberRadiusZ);
  if (std::numeric_limits<size_t>::max() != MaxMemory) {
    str += std::format(" MemLimit {:.2f} GB", double(MaxMemory) / (1024.f * 1024.f * 1024.f));
  }
  return str;
}

} // namespace o2::itsmft

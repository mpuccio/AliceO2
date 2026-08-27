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
#include <array>
#include <cctype>
#include <cmath>
#include <format>
#include <limits>
#include <string_view>
#include <vector>

#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITSMFTTracking/Constants.h"
#include "MFTTracking/Constants.h"

namespace
{
constexpr bool iequals(std::string_view a, std::string_view b)
{
  return std::equal(a.begin(), a.end(), b.begin(), b.end(),
                    [](char x, char y) { return std::tolower(x) == std::tolower(y); });
}
} // namespace

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
    str += std::format(" MaxHoles:{}", MaxHoles);
  }
  if (!InactiveLayerMask.empty()) {
    str += std::format(" InactiveMask:{}", InactiveLayerMask.asString());
  }
  if (!SeedingLayers.empty()) {
    str += std::format(" SeedingLayers:{}", SeedingLayers.asString());
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

void resetDetectorDefaults(TrackingParameters& p, detectors::DetID::ID detId)
{
  if (detId == detectors::DetID::ITS) {
    p = TrackingParameters{};
    p.MinPt.assign(TrackerParamConfig<detectors::DetID::ITS>::MaxTrackLength - TrackerParamConfig<detectors::DetID::ITS>::MinTrackLength + 1, 0.f);
    return;
  }

  if (detId == detectors::DetID::MFT) {
    namespace mftc = o2::mft::constants;
    namespace mft = mftc::mft;
    constexpr int nLayers = o2::mft::constants::mft::LayersNumber;

    p = TrackingParameters{};
    p.NLayers = nLayers;
    p.LayerZ.clear();
    p.LayerZ.reserve(nLayers);
    for (float z : mft::LayerZCoordinate()) {
      p.LayerZ.push_back(std::abs(z));
    }
    p.LayerColHalfExtent.assign(mftc::index_table::RMax.begin(), mftc::index_table::RMax.end());
    p.IndexRowMin = -20.f;
    p.IndexRowMax = 20.f;
    p.LayerRadii.resize(nLayers);
    for (int i{0}; i < nLayers; ++i) {
      p.LayerRadii[i] = 0.5f * (mftc::index_table::RMin[i] + mftc::index_table::RMax[i]);
    }
    p.LayerxX0.assign(tracking::kNominalMFTLayerX0.begin(), tracking::kNominalMFTLayerX0.end());
    p.LayerResolution.assign(nLayers, mft::Resolution);
    p.SystError2Row.assign(nLayers, 0.f);
    p.SystError2Col.assign(nLayers, 0.f);
    p.AddTimeError.assign(nLayers, 0u);
    p.ColBins = 64;
    p.RowBins = 128;
    p.UseDiamond = true;
    p.PerPrimaryVertexProcessing = false;
    p.StartLayerMask = (1u << nLayers) - 1u;
    p.MinPt.assign(TrackerParamConfig<detectors::DetID::MFT>::MaxTrackLength - TrackerParamConfig<detectors::DetID::MFT>::MinTrackLength + 1, 0.f);
    return;
  }

  LOGP(fatal, "Unsupported detector id {} in resetDetectorDefaults", static_cast<int>(detId));
}

namespace TrackingMode
{

Type fromString(std::string_view str)
{
  constexpr std::array smodes = {
    std::pair{"sync", Sync},
    std::pair{"async", Async},
    std::pair{"cosmics", Cosmics},
    std::pair{"unset", Unset},
    std::pair{"off", Off}};

  const auto it = std::find_if(smodes.begin(), smodes.end(), [&str](const auto& pair) {
    return iequals(str, pair.first);
  });
  if (it == smodes.end()) {
    LOGP(fatal, "Unrecognized CA tracking mode '{}'", str);
  }
  return it->second;
}

std::string toString(Type mode)
{
  switch (mode) {
    case Sync:
      return "sync";
    case Async:
      return "async";
    case Cosmics:
      return "cosmics";
    case Unset:
      return "unset";
    case Off:
      return "off";
  }
  LOGP(fatal, "Unrecognized CA tracking mode {}", static_cast<int>(mode));
  return "";
}

std::vector<TrackingParameters> getTrackingParameters(detectors::DetID::ID detId, Type mode)
{
  if (detId == detectors::DetID::ITS) {
    // Only Sync is implemented for ITS; reject all other modes explicitly.
    if (mode != Sync) {
      LOGP(fatal, "ITS common-CA tracking mode '{}' is not supported yet; only 'sync' is implemented", toString(mode));
    }

    const auto& tc = ITSCommonCATrackerParam::Instance();
    std::vector<TrackingParameters> trackParams(1);
    auto& p = trackParams[0];
    resetDetectorDefaults(p, detectors::DetID::ITS);
    p.MinTrackLength = TrackerParamConfig<detectors::DetID::ITS>::MinTrackLength;
    p.DropTFUponFailure = tc.dropTFUponFailure;
    p.PrintMemory = tc.printMemory;
    p.MaxMemory = tc.maxMemory;
    p.SaveTimeBenchmarks = tc.saveTimeBenchmarks;
    p.UseDiamond = tc.useDiamond;
    for (int iD = 0; iD < 3; ++iD) {
      p.Diamond[iD] = tc.diamondPos[iD];
    }
    p.PVres = tc.pvRes > 0 ? tc.pvRes : p.PVres;
    return trackParams;
  }
  if (detId != detectors::DetID::MFT) {
    LOGP(fatal, "Unsupported detector id {} in getTrackingParameters", static_cast<int>(detId));
  }

  const auto& tc = TrackerParamConfig<detectors::DetID::MFT>::Instance();
  std::vector<TrackingParameters> trackParams;

  if (mode == Off) {
    return trackParams;
  }
  if (mode == Unset) {
    LOGP(fatal, "CA tracking mode is unset; set --tracking-mode or {}.trackingMode", TrackerParamConfig<detectors::DetID::MFT>::getParamName());
  }

  if (mode == Async) {
    trackParams.resize(3);
    for (auto& param : trackParams) {
      resetDetectorDefaults(param, detId);
    }

    trackParams[1].TrackletMinPt = 0.15f;
    trackParams[2].TrackletMinPt = 0.08f;

    trackParams[0].MinPt[0] = 1.f / 12.f; // 10 clusters
    trackParams[1].MinPt[0] = 1.f / 12.f;

    trackParams[2].MinTrackLength = TrackerParamConfig<detectors::DetID::MFT>::MinTrackLength;
    trackParams[2].MinPt[0] = 1.f / 12.f; // 10 clusters
    trackParams[2].MinPt[1] = 1.f / 8.f;  // 9 clusters
    trackParams[2].MinPt[2] = 1.f / 5.f;  // 8 clusters
    trackParams[2].MinPt[3] = 1.f / 3.f;  // 7 clusters
    trackParams[2].MinPt[4] = 1.f / 2.f;  // 6 clusters
    trackParams[2].MinPt[5] = 1.f / 1.f;  // 5 clusters

    for (int ip = 0; ip < static_cast<int>(trackParams.size()); ip++) {
      auto& param = trackParams[ip];
      if (ip < o2::its::constants::MaxIter) {
        if (tc.startLayerMask[ip] > 0) {
          param.StartLayerMask = tc.startLayerMask[ip];
        }
        if (tc.minTrackLgtIter[ip] > 0) {
          param.MinTrackLength = tc.minTrackLgtIter[ip];
        }
        for (int ilg = tc.MaxTrackLength; ilg >= tc.MinTrackLength; ilg--) {
          const int lslot0 = tc.MaxTrackLength - ilg;
          const int lslot = lslot0 + ip * (tc.MaxTrackLength - tc.MinTrackLength + 1);
          if (tc.minPtIterLgt[lslot] > 0.f) {
            param.MinPt[lslot0] = tc.minPtIterLgt[lslot];
          }
        }
      }
    }
  } else if (mode == Sync) {
    trackParams.resize(1);
    resetDetectorDefaults(trackParams[0], detId);
    trackParams[0].MinTrackLength = TrackerParamConfig<detectors::DetID::MFT>::MinTrackLength;
  } else if (mode == Cosmics) {
    trackParams.resize(1);
    resetDetectorDefaults(trackParams[0], detId);
    trackParams[0].MinTrackLength = TrackerParamConfig<detectors::DetID::MFT>::MinTrackLength;
    trackParams[0].ColBins = 32;
    trackParams[0].RowBins = 64;
    trackParams[0].PVres = 1.e5f;
    trackParams[0].MaxChi2ClusterAttachment = 60.f;
    trackParams[0].MaxChi2NDF = 40.f;
  } else {
    LOGP(fatal, "Unsupported CA tracking mode {}", toString(mode));
  }

  for (auto& param : trackParams) {
    param.PassFlags.reset();
  }
  if (!trackParams.empty()) {
    trackParams[0].PassFlags.set(IterationStep::FirstPass, IterationStep::RebuildClusterLUT);
  }

  const float bFactor = std::abs(o2::base::Propagator::Instance()->getNominalBz()) / 5.0066791f;
  const float bFactorTracklets = bFactor < 0.01f ? 1.f : bFactor;

  for (auto& p : trackParams) {
    p.TrackletMinPt *= bFactorTracklets;
    for (int ilg = tc.MaxTrackLength; ilg >= tc.MinTrackLength; ilg--) {
      const int lslot = tc.MaxTrackLength - ilg;
      if (lslot < static_cast<int>(p.MinPt.size())) {
        p.MinPt[lslot] *= bFactor;
      }
    }

    p.ReseedIfShorter = tc.reseedIfShorter;
    p.RepeatRefitOut = tc.repeatRefitOut;
    p.ShiftRefToCluster = tc.shiftRefToCluster;
    p.CreateArtefactLabels = tc.createArtefactLabels;
    p.PrintMemory = tc.printMemory;
    p.MaxMemory = tc.maxMemory;
    p.DropTFUponFailure = tc.dropTFUponFailure;
    p.SaveTimeBenchmarks = tc.saveTimeBenchmarks;
    p.FataliseUponFailure = tc.fataliseUponFailure;
    p.AllowSharingFirstCluster = tc.allowSharingFirstCluster;
    p.SharedClusterMaxDeltaPhi = tc.sharedClusterMaxDeltaPhi;
    p.SharedClusterMaxDeltaEta = tc.sharedClusterMaxDeltaEta;
    p.SharedClusterOppositeSign = tc.sharedClusterOppositeSign;
    p.PerPrimaryVertexProcessing = tc.perPrimaryVertexProcessing;

    const auto iter = &p - trackParams.data();
    if (iter < o2::its::constants::MaxIter) {
      p.MaxHoles = tc.maxHolesIter[iter];
    }

    if (tc.useMatCorrTGeo) {
      p.CorrType = o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrTGeo;
    } else if (tc.useFastMaterial) {
      p.CorrType = o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrNONE;
    } else {
      p.CorrType = o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrLUT;
    }

    for (int i{0}; i < TrackerParamConfig<detectors::DetID::MFT>::getNLayers(); ++i) {
      p.SystError2Row[i] = tc.sysErr2Row[i] > 0 ? tc.sysErr2Row[i] : p.SystError2Row[i];
      p.SystError2Col[i] = tc.sysErr2Col[i] > 0 ? tc.sysErr2Col[i] : p.SystError2Col[i];
      p.AddTimeError[i] = tc.addTimeError[i];
    }

    p.MaxChi2ClusterAttachment = tc.maxChi2ClusterAttachment > 0 ? tc.maxChi2ClusterAttachment : p.MaxChi2ClusterAttachment;
    p.MaxChi2NDF = tc.maxChi2NDF > 0 ? tc.maxChi2NDF : p.MaxChi2NDF;
    p.ColBins = tc.LUTbinsU > 0 ? tc.LUTbinsU : p.ColBins;
    p.RowBins = tc.LUTbinsV > 0 ? tc.LUTbinsV : p.RowBins;
    p.PVres = tc.pvRes > 0 ? tc.pvRes : p.PVres;
    p.NSigmaCut *= tc.nSigmaCut > 0 ? tc.nSigmaCut : 1.f;
    p.TrackletMinPt *= tc.minPt > 0 ? tc.minPt : 1.f;
    for (int iD{0}; iD < 3; ++iD) {
      p.Diamond[iD] = tc.diamondPos[iD];
    }
    if (detId == detectors::DetID::MFT) {
      p.UseDiamond = true;
      p.PerPrimaryVertexProcessing = false;
    } else {
      p.UseDiamond = tc.useDiamond;
    }
  }

  if (trackParams.size() > static_cast<size_t>(tc.nIterations)) {
    trackParams.resize(tc.nIterations);
  }

  return trackParams;
}

} // namespace TrackingMode
} // namespace o2::itsmft

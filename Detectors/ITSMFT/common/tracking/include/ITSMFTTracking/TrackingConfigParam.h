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

#ifndef ALICEO2_ITSMFT_TRACKING_CONFIG_PARAM_H_
#define ALICEO2_ITSMFT_TRACKING_CONFIG_PARAM_H_

#include <array>
#include <limits>
#include <string_view>

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include "DetectorsCommonDataFormats/DetID.h"

namespace o2::itsmft::tracking
{
/// ITS CA layer count.
constexpr int ITSNLayers = 7;
/// MFT CA half-disk layer count.
constexpr int MFTNLayers = 10;
/// Maximum CA iterations.
constexpr int MaxIter = 4;
inline constexpr std::array<float, ITSNLayers> kITSLookupZHalfExtent{
  16.333f + 1.f, 16.333f + 1.f, 16.333f + 1.f,
  42.140f + 1.f, 42.140f + 1.f, 73.745f + 1.f, 73.745f + 1.f};
} // namespace o2::itsmft::tracking

namespace o2::itsmft
{

/// ITS vertexer settings for the existing configuration boundary.
struct VertexerParamConfig : public o2::conf::ConfigurableParamHelper<VertexerParamConfig> {
  bool saveTimeBenchmarks = false; // dump metrics to a file

  int nIterations = 1;         // Number of vertexing passes.
  int vertPerRofThreshold = 0; // Vertices per ROF that trigger a second iteration.

  // Geometrical cuts for Pb-Pb tracklet selection.
  float zCut = 0.002f;
  float phiCut = 0.005f;
  float pairCut = 0.017321f;
  float clusterCut = 0.170048f;
  float coarseZWindow = 0.055458f;
  float seedDedupZCut = 0.116685f;
  float refitDedupZCut = 0.039855f;
  float duplicateZCut = 0.200097f;
  float finalSelectionZCut = 0.034535f;
  float duplicateDistance2Cut = 0.005117f;
  float tanLambdaCut = 0.002f; // tanLambda = deltaZ/deltaR
  float nSigmaCut = 0.0164651f;
  float maxZPositionAllowed = 25.f; // 4x sZ of the beam

  // Artefact selection.
  int clusterContributorsCut = 3; // Minimum contributors for an accepted final vertex.
  int suppressLowMultDebris = 16; // Suppress lower-multiplicity vertices after finding one in a ROF.
  int seedMemberRadiusTime = 0;
  int seedMemberRadiusZ = 2;
  int maxTrackletsPerCluster = 100;
  int phiSpan = -1;
  int zSpan = -1;
  int ZBins = 1;     // Number of z bins in the z-phi index table.
  int PhiBins = 128; // Number of phi bins in the z-phi index table.

  bool useTruthSeeding{false}; // Replace seed vertices with MC truth.

  int nThreads = 1;
  bool printMemory = false;
  size_t maxMemory = std::numeric_limits<size_t>::max();
  bool dropTFUponFailure = false;

  O2ParamDef(VertexerParamConfig, "ITSVertexerParam");
};

/// Minimal configuration for opt-in ITS common-CA tracking.
/// It is distinct from the dormant TrackerParamConfig<DetID::ITS> and does
/// not use the registered name "ITSCATrackerParam", which belongs to the
/// legacy o2::its::TrackerParamConfig.
/// It exposes only fields consumed by
/// TrackingMode::getTrackingParameters(DetID::ITS, Sync); legacy-only knobs
/// are intentionally absent. Defaults match TrackingParameters' Sync
/// baseline, so leaving this configuration unchanged preserves the detector
/// defaults.
///
/// diamondPos, pvRes, and useDiamond define the static vertex/beam constraint
/// consumed by the shared TrackerTraits.
struct ITSCommonCATrackerParam : public o2::conf::ConfigurableParamHelper<ITSCommonCATrackerParam> {
  bool dropTFUponFailure = false;
  bool printMemory = false;
  size_t maxMemory = std::numeric_limits<size_t>::max();
  bool saveTimeBenchmarks = false;
  bool useDiamond = false;
  float diamondPos[3] = {0.f, 0.f, 0.f}; // Diamond vertex position when useDiamond is set.
  float pvRes = -1.f;                    // Diamond-vertex PV resolution; <=0 keeps the default.
  uint16_t holeLayerMask = 0;            // Detector layers that may be absent from accepted tracks.

  /// Number of tbb::task_arena threads for the ITS common-CA tracker.
  /// This dedicated field is separate from the legacy ITS configuration.
  /// Must be > 0; validated where consumed because ConfigurableParam
  /// structs cannot reject construction.
  int nThreads = 1;

  O2ParamDef(ITSCommonCATrackerParam, "ITSCommonCATrackerParam");
};

template <int N>
struct TrackerParamConfig : public o2::conf::ConfigurableParamHelper<TrackerParamConfig<N>> {
  static constexpr std::string_view getParamName()
  {
    return N == o2::detectors::DetID::ITS ? TrackerParamName[0] : TrackerParamName[1];
  }

  static constexpr int MinTrackLength = 4;
  static constexpr int MaxTrackLength = N == o2::detectors::DetID::ITS ? o2::itsmft::tracking::ITSNLayers
                                                                       : o2::itsmft::tracking::MFTNLayers;

  static constexpr int getNLayers()
  {
    return N == o2::detectors::DetID::ITS ? o2::itsmft::tracking::ITSNLayers : o2::itsmft::tracking::MFTNLayers;
  }

  bool useMatCorrTGeo = false;                                                                    // Use full geometry for material correction; default is the LUT.
  bool useFastMaterial = false;                                                                   // Use a faster material approximation in fits.
  int addTimeError[getNLayers()] = {0};                                                           // Tracking window width in BC.
  int minTrackLgtIter[o2::itsmft::tracking::MaxIter] = {};                                        // Minimum track length per iteration; <=0 uses code defaults.
  uint8_t startLayerMask[o2::itsmft::tracking::MaxIter] = {};                                     // Start-layer mask per iteration.
  int maxHolesIter[o2::itsmft::tracking::MaxIter] = {};                                           // Maximum missing internal layers per iteration.
  uint16_t holeLayerMask = 0;                                                                     // Detector layers that may be absent from accepted tracks.
  float minPtIterLgt[o2::itsmft::tracking::MaxIter * (MaxTrackLength - MinTrackLength + 1)] = {}; // Minimum pT by track length; <=0 uses code defaults.
  float sysErr2Row[getNLayers()] = {0};                                                           // Systematic error squared along local X per layer.
  float sysErr2Col[getNLayers()] = {0};                                                           // Systematic error squared along local Z per layer.
  float maxChi2ClusterAttachment = -1.f;
  float maxChi2NDF = -1.f;
  float nSigmaCut = -1.f;
  float deltaTanLres = -1.f;
  float minPt = -1.f;
  float pvRes = -1.f;
  int LUTbinsU = N == o2::detectors::DetID::MFT ? 64 : -1;                              // LUT bins along the first coordinate (ITS: phi, MFT: global x).
  int LUTbinsV = N == o2::detectors::DetID::MFT ? 128 : -1;                             // LUT bins along the second coordinate (ITS: z, MFT: global y).
  float diamondPos[3] = {0.f, 0.f, 0.f};                                                // Diamond vertex for MFT seeds, or ITS when useDiamond.
  bool useDiamond = N == o2::detectors::DetID::MFT;                                     // MFT always uses diamond; ITS opts in.
  bool perPrimaryVertexProcessing = false;                                              // Track separately for each vertex hypothesis.
  bool saveTimeBenchmarks = false;                                                      // Dump metrics to a file.
  bool overrideBeamEstimation = false;                                                  // ITS only: use the meanVertex CCDB beam seed.
  int trackingMode = -1;                                                                // -1: unset (use --tracking-mode); 0: sync, 1: async, 2: cosmics.
  bool doUPCIteration = false;                                                          // Add an iteration for UPC events on tagged vertices; use nIterations=2.
  int nIterations = N == o2::detectors::DetID::MFT ? 1 : o2::itsmft::tracking::MaxIter; // Number of iterations.
  int reseedIfShorter = 6;                                                              // Reseed short final-refit tracks with a circle.
  bool shiftRefToCluster{true};                                                         // Shift the linearization reference to the cluster after update.
  bool repeatRefitOut{false};                                                           // Repeat outward refit using the inward refit as a seed.
  bool createArtefactLabels{false};                                                     // Create labels for artefacts on the fly.

  int nThreads = 1;
  bool printMemory = false;
  size_t maxMemory = std::numeric_limits<size_t>::max();
  bool dropTFUponFailure = false;
  bool fataliseUponFailure = true; // Granular fatalisation control in async mode.

  // Selection of tracks sharing clusters.
  bool allowSharingFirstCluster = false;  // Allow sharing the first cluster.
  float sharedClusterMaxDeltaPhi = 0.05f; // Maximum delta phi at the cluster.
  float sharedClusterMaxDeltaEta = 0.03f; // Maximum delta eta at the cluster.
  bool sharedClusterOppositeSign = false; // Require opposite-sign tracklets.

  O2ParamDef(TrackerParamConfig, getParamName().data());

 private:
  static constexpr std::string_view TrackerParamName[2] = {"ITSCATrackerParam", "MFTCATrackerParam"};

  static_assert(N == o2::detectors::DetID::ITS || N == o2::detectors::DetID::MFT, "only DetID::ITS or DetID::MFT are allowed");
};

template <int N>
TrackerParamConfig<N> TrackerParamConfig<N>::sInstance;

} // namespace o2::itsmft

namespace framework
{
template <typename T>
struct is_messageable;
template <>
struct is_messageable<o2::itsmft::VertexerParamConfig> : std::true_type {
};
template <>
struct is_messageable<o2::itsmft::TrackerParamConfig<o2::detectors::DetID::ITS>> : std::true_type {
};
template <>
struct is_messageable<o2::itsmft::TrackerParamConfig<o2::detectors::DetID::MFT>> : std::true_type {
};
} // namespace framework

#endif /* ALICEO2_ITSMFT_TRACKING_CONFIG_PARAM_H_ */

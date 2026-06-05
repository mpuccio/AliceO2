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

#include <limits>
#include <string_view>

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/Constants.h"

namespace o2::itsmft
{

namespace tracking_constants
{
constexpr int MaxIter = 4;
} // namespace tracking_constants

/// ITS vertexer settings (not used for MFT)
struct VertexerParamConfig : public o2::conf::ConfigurableParamHelper<VertexerParamConfig> {
  bool saveTimeBenchmarks = false; // dump metrics on file

  int nIterations = 1;         // Number of vertexing passes to perform.
  int vertPerRofThreshold = 0; // Maximum number of vertices per ROF to trigger second a iteration.

  // geometrical cuts for tracklet selection for Pb-Pb
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

  // Artefacts selections
  int clusterContributorsCut = 3; // minimum number of contributors for an accepted final vertex
  int suppressLowMultDebris = 16; // suppress all vertices below this threshold if a vertex was already found in a rof
  int seedMemberRadiusTime = 0;
  int seedMemberRadiusZ = 2;
  int maxTrackletsPerCluster = 100;
  int phiSpan = -1;
  int zSpan = -1;
  int ZBins = 1;     // z-phi index table configutation: number of z bins
  int PhiBins = 128; // z-phi index table configutation: number of phi bins

  bool useTruthSeeding{false}; // overwrite seeding vertices with MC truth

  int nThreads = 1;
  bool printMemory = false;
  size_t maxMemory = std::numeric_limits<size_t>::max();
  bool dropTFUponFailure = false;

  O2ParamDef(VertexerParamConfig, "ITSVertexerParam");
};

template <int N>
struct TrackerParamConfig : public o2::conf::ConfigurableParamHelper<TrackerParamConfig<N>> {
  static constexpr std::string_view getParamName()
  {
    return N == o2::detectors::DetID::ITS ? TrackerParamName[0] : TrackerParamName[1];
  }

  static constexpr int MinTrackLength = N == o2::detectors::DetID::ITS ? 4 : 5;
  static constexpr int MaxTrackLength = N == o2::detectors::DetID::ITS ? o2::itsmft::tracking::constants::ITSNLayers
                                                                        : o2::itsmft::tracking::constants::MFTNLayers;

  static constexpr int getNLayers()
  {
    return N == o2::detectors::DetID::ITS ? o2::itsmft::tracking::constants::ITSNLayers : o2::itsmft::tracking::constants::MFTNLayers;
  }

  bool useMatCorrTGeo = false;                                                         // use full geometry to corect for material budget accounting in the fits. Default is to use the material budget LUT.
  bool useFastMaterial = false;                                                        // use faster material approximation for material budget accounting in the fits.
  int addTimeError[getNLayers()] = {0};                                                // configure the width of the window in BC to be considered for the tracking.
  int minTrackLgtIter[tracking_constants::MaxIter] = {};                               // minimum track length at each iteration, used only if >0, otherwise use code defaults
  uint8_t startLayerMask[tracking_constants::MaxIter] = {};                            // mask of start layer for this iteration (if >0)
  int maxHolesIter[tracking_constants::MaxIter] = {};                                  // maximum number of missing internal layers allowed in the CA topology for each iteration
  uint16_t holeLayerMaskIter[tracking_constants::MaxIter] = {};                          // layers that may be skipped by the CA topology for each iteration
  float minPtIterLgt[tracking_constants::MaxIter * (MaxTrackLength - MinTrackLength + 1)] = {}; // min.pT for given track length at this iteration, used only if >0, otherwise use code defaults
  float sysErr2Row[getNLayers()] = {0}; // systematic error^2 along ALPIDE rows (local X) per layer
  float sysErr2Col[getNLayers()] = {0}; // systematic error^2 along ALPIDE columns (local Z) per layer
  float maxChi2ClusterAttachment = -1.f;
  float maxChi2NDF = -1.f;
  float nSigmaCut = -1.f;
  float deltaTanLres = -1.f;
  float minPt = -1.f;
  float pvRes = -1.f;
  int LUTbinsU = -1; // number of LUT bins along the first coordinate (ITS: phi, MFT: global x)
  int LUTbinsV = -1; // number of LUT bins along the second coordinate (ITS: z, MFT: global y)
  float diamondPos[3] = {0.f, 0.f, 0.f};   // override the position of the vertex
  bool useDiamond = false;                 // enable overriding the vertex position
  bool perPrimaryVertexProcessing = false; // perform the full tracking considering the vertex hypotheses one at the time.
  bool saveTimeBenchmarks = false;         // dump metrics on file
  bool overrideBeamEstimation = false;     // use beam position from meanVertex CCDB object
  int trackingMode = -1;                   // -1: unset, 0=sync, 1=async, 2=cosmics used by gpuwf only
  bool doUPCIteration = false;             // Perform an additional iteration for UPC events on tagged vertices. You want to combine this config with VertexerParamConfig.nIterations=2
  int nIterations = tracking_constants::MaxIter;                                       // overwrite the number of iterations
  int reseedIfShorter = 6;                 // for the final refit reseed the track with circle if they are shorter than this value
  bool shiftRefToCluster{true};            // TrackFit: after update shift the linearization reference to cluster
  bool repeatRefitOut{false};              // repeat outward refit using inward refit as a seed
  bool createArtefactLabels{false};        // create on-the-fly labels for the artefacts

  int nThreads = 1;
  bool printMemory = false;
  size_t maxMemory = std::numeric_limits<size_t>::max();
  bool dropTFUponFailure = false;
  bool fataliseUponFailure = true;       // granular management of the fatalisation in async mode

  // Selections on tracks sharing clusters
  bool allowSharingFirstCluster = false;  // allow first cluster sharing among tracks
  float sharedClusterMaxDeltaPhi = 0.05f; // Maximum allowed delta phi at the cluster position
  float sharedClusterMaxDeltaEta = 0.03f; // Maximum allowed delta eta at the cluster position
  bool sharedClusterOppositeSign = false; // Require opposite sign of the tracklets

  O2ParamDef(TrackerParamConfig, getParamName().data());

 private:
  static constexpr std::string_view TrackerParamName[2] = {"ITSCATrackerParam", "MFTCATrackerParam"};
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

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

namespace o2::itsmft::tracking
{
/// ITS CA layer count (matches legacy o2::its::TrackingParameters::NLayers default).
constexpr int ITSNLayers = 7;
/// MFT CA half-disk layer count (same as o2::mft::constants::mft::LayersNumber).
constexpr int MFTNLayers = 10;
/// Max CA iterations (same as o2::its::constants::MaxIter).
constexpr int MaxIter = 4;
} // namespace o2::itsmft::tracking

namespace o2::itsmft
{

/// ITS vertexer settings (header only for now; not registered or used while ITS
/// tracking stays on O2::ITStracking)
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

/// Dedicated, minimal configuration for the opt-in ITS common-CA tracking
/// path (workflow-onboarding Slice 1: Sync mode only). Deliberately a
/// distinct, non-templated type: it neither instantiates nor renames the
/// still-dormant TrackerParamConfig<DetID::ITS> below, and its registered
/// name is not "ITSCATrackerParam" -- that name already belongs to the
/// unrelated, frozen legacy o2::its::TrackerParamConfig in O2::ITStracking,
/// which stays transitively linked into any executable using this library
/// (O2::ITSMFTTracking -> O2::ITStracking, for the frozen final-refit
/// boundary) regardless of which type is used here.
/// Exposes only fields actually consumed by
/// TrackingMode::getTrackingParameters(DetID::ITS, Sync); legacy-only knobs
/// (TrackFollower*, UPC orchestration, vertexing, FataliseUponFailure) are
/// intentionally absent, not merely unread. Defaults match
/// TrackingParameters' own struct defaults (the documented Sync baseline),
/// so an unmodified ITSCommonCATrackerParam changes nothing relative to
/// resetDetectorDefaults(..., DetID::ITS) alone.
///
/// workflow-onboarding Slice 2: diamondPos/pvRes added alongside the
/// pre-existing useDiamond, mirroring the legacy
/// o2::its::TrackerParamConfig fields of the same name/semantics. These
/// three fields together are the one Sync vertex/beam-constraint mode the
/// opt-in ITS common-CA workflow can reproduce faithfully without a real
/// per-event vertexer (see ITSCAWorkflow/ConfigPreflight.h): a static,
/// config-supplied diamond vertex consumed identically by the shared
/// TrackerTraits (Detectors/ITSMFT/common/tracking/src/TrackerTraits.cxx)
/// regardless of detector.
struct ITSCommonCATrackerParam : public o2::conf::ConfigurableParamHelper<ITSCommonCATrackerParam> {
  bool dropTFUponFailure = false;
  bool printMemory = false;
  size_t maxMemory = std::numeric_limits<size_t>::max();
  bool saveTimeBenchmarks = false;
  bool useDiamond = false;
  float diamondPos[3] = {0.f, 0.f, 0.f}; // diamond vertex position, consumed only when useDiamond is set
  float pvRes = -1.f;                    // PV resolution override for the diamond vertex; <=0 keeps the struct default

  /// Number of tbb::task_arena threads for the ITS common-CA tracker
  /// (ITSMFTTrackingInterface<7>::initialiseTracker()). Deliberately this
  /// struct's own field, not TrackerParamRef<ITS>::get().nThreads (the
  /// frozen legacy o2::its::TrackerParamConfig, registered "ITSCATrackerParam"
  /// -- see this struct's own doc comment above on why the two are kept
  /// distinct): the common ITS interface must read this dedicated value.
  /// Must be > 0; validated where consumed, not here (ConfigurableParam
  /// structs cannot fail construction).
  int nThreads = 1;

  O2ParamDef(ITSCommonCATrackerParam, "ITSCommonCATrackerParam");
};

struct ITSMFTCommonCAOutputParam : public o2::conf::ConfigurableParamHelper<ITSMFTCommonCAOutputParam> {
  bool useCommonTrackOutput = false;
  O2ParamDef(ITSMFTCommonCAOutputParam, "ITSMFTCommonCAOutputParam");
};

inline bool isCommonTrackOutputEnabled()
{
  return ITSMFTCommonCAOutputParam::Instance().useCommonTrackOutput;
}

/// Gate 4 C4: default-false, ROOT-visible opt-in gate for the combined
/// ITS+MFT DPL workflow (Detectors/ITSMFT/common/workflow-combined-ca/,
/// o2-itsmft-combined-ca-tracker-workflow). A distinct struct/registered
/// name from ITSMFTCommonCAOutputParam (which only ever toggles output
/// *format* on the single-detector opt-in workflows): this flag instead
/// gates whether the combined executable's own defineDataProcessing()
/// refuses to run at all. Checked before any DataProcessorSpec is
/// constructed -- see the combined workflow's own ConfigPreflight -- so a
/// pipeline template that accidentally invokes this binary without
/// explicitly setting `enabled=true` fatals immediately rather than
/// silently replacing o2-its-ca-tracker-workflow/o2-mft-ca-tracker-workflow
/// or either legacy workflow.
struct ITSMFTCombinedCATrackerParam : public o2::conf::ConfigurableParamHelper<ITSMFTCombinedCATrackerParam> {
  bool enabled = false;
  O2ParamDef(ITSMFTCombinedCATrackerParam, "ITSMFTCombinedCATrackerParam");
};

template <int N>
struct TrackerParamConfig : public o2::conf::ConfigurableParamHelper<TrackerParamConfig<N>> {
  static constexpr std::string_view getParamName()
  {
    return N == o2::detectors::DetID::ITS ? TrackerParamName[0] : TrackerParamName[1];
  }

  static constexpr int MinTrackLength = N == o2::detectors::DetID::ITS ? 4 : 5;
  static constexpr int MaxTrackLength = N == o2::detectors::DetID::ITS ? o2::itsmft::tracking::ITSNLayers
                                                                        : o2::itsmft::tracking::MFTNLayers;

  static constexpr int getNLayers()
  {
    return N == o2::detectors::DetID::ITS ? o2::itsmft::tracking::ITSNLayers : o2::itsmft::tracking::MFTNLayers;
  }

  bool useMatCorrTGeo = false;                                                         // use full geometry to corect for material budget accounting in the fits. Default is to use the material budget LUT.
  bool useFastMaterial = false;                                                        // use faster material approximation for material budget accounting in the fits.
  int addTimeError[getNLayers()] = {0};                                                // configure the width of the window in BC to be considered for the tracking.
  int minTrackLgtIter[o2::itsmft::tracking::MaxIter] = {};                               // minimum track length at each iteration, used only if >0, otherwise use code defaults
  uint8_t startLayerMask[o2::itsmft::tracking::MaxIter] = {};                            // mask of start layer for this iteration (if >0)
  int maxHolesIter[o2::itsmft::tracking::MaxIter] = {};                                  // maximum number of missing internal layers allowed in the CA topology for each iteration
  uint16_t holeLayerMaskIter[o2::itsmft::tracking::MaxIter] = {};                          // layers that may be skipped by the CA topology for each iteration
  float minPtIterLgt[o2::itsmft::tracking::MaxIter * (MaxTrackLength - MinTrackLength + 1)] = {}; // min.pT for given track length at this iteration, used only if >0, otherwise use code defaults
  float sysErr2Row[getNLayers()] = {0}; // systematic error^2 along ALPIDE rows (local X) per layer
  float sysErr2Col[getNLayers()] = {0}; // systematic error^2 along ALPIDE columns (local Z) per layer
  float maxChi2ClusterAttachment = -1.f;
  float maxChi2NDF = -1.f;
  float nSigmaCut = -1.f;
  float deltaTanLres = -1.f;
  float minPt = -1.f;
  float pvRes = -1.f;
  int LUTbinsU = N == o2::detectors::DetID::MFT ? 64 : -1; // number of LUT bins along the first coordinate (ITS: phi, MFT: global x)
  int LUTbinsV = N == o2::detectors::DetID::MFT ? 128 : -1; // number of LUT bins along the second coordinate (ITS: z, MFT: global y)
  float diamondPos[3] = {0.f, 0.f, 0.f};   // diamond vertex position (MFT CA tracklet seed; ITS when useDiamond)
  bool useDiamond = N == o2::detectors::DetID::MFT; // MFT CA: always diamond; ITS: opt-in via param
  bool perPrimaryVertexProcessing = false; // perform the full tracking considering the vertex hypotheses one at the time.
  bool saveTimeBenchmarks = false;         // dump metrics on file
  bool overrideBeamEstimation = false;     // ITS only: meanVertex CCDB beam seed (MFT CA always uses diamond)
  int trackingMode = -1; // -1: unset (use --tracking-mode), 0=sync, 1=async, 2=cosmics
  bool doUPCIteration = false;             // Perform an additional iteration for UPC events on tagged vertices. You want to combine this config with VertexerParamConfig.nIterations=2
  int nIterations = N == o2::detectors::DetID::MFT ? 1 : o2::itsmft::tracking::MaxIter; // overwrite the number of iterations
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
  float cellRoadRCut = -1.f; // MFT: max distance to seed line; <=0 uses default (0.05 cm)
  float trackletMinAbsX = -1.f; // MFT: min |x| (cm) for tracklet seeds and accepted tracks; <=0 uses code default

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

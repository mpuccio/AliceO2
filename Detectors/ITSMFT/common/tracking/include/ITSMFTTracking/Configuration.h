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
/// \file Configuration.h
/// \brief Shared CA tracking configuration for ITS and MFT
///

#ifndef ALICEO2_ITSMFT_TRACKING_CONFIGURATION_H_
#define ALICEO2_ITSMFT_TRACKING_CONFIGURATION_H_

#include <cstdint>
#ifndef GPUCA_GPUCODE_DEVICE
#include <limits>
#include <string>
#include <string_view>
#include <vector>
#endif

#include "CommonUtils/EnumFlags.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/Constants.h"
#include "ITSMFTTracking/LayerMask.h"

namespace o2::itsmft
{

inline constexpr int ClustersPerCell = 3;

// Steering of dedicated steps in an iteration
enum class IterationStep : uint8_t {
  FirstPass = 0,
  RebuildClusterLUT,
  UseUPCMask,
  SelectUPCVertices,
  ResetVertices,
  SkipROFsAboveThreshold,
  MarkVerticesAsUPC,
};
using IterationSteps = o2::utils::EnumFlags<IterationStep>;

struct TrackingParameters {
  int CellMinimumLevel() const noexcept
  {
    const int minClusters = MinTrackLength - (MaxHoles > 0 ? MaxHoles : 0);
    const int effectiveMinClusters = minClusters > ClustersPerCell ? minClusters : ClustersPerCell;
    return effectiveMinClusters - ClustersPerCell + 1;
  }
  int NeighboursPerRoad() const noexcept { return NLayers - 3; }
  int CellsPerRoad() const noexcept { return NLayers - 2; }
  int TrackletsPerRoad() const noexcept { return NLayers - 1; }
  std::string asString() const;

  IterationSteps PassFlags{IterationStep::FirstPass, IterationStep::RebuildClusterLUT};
  int NLayers = tracking::constants::ITSNLayers;
  std::vector<uint32_t> AddTimeError = {0, 0, 0, 0, 0, 0, 0};
  std::vector<float> LayerZ = {16.333f + 1, 16.333f + 1, 16.333f + 1, 42.140f + 1, 42.140f + 1, 73.745f + 1, 73.745f + 1};
  std::vector<float> LayerColHalfExtent{}; // index-table column half extent (ITS: z, MFT: global x); falls back to LayerZ
  float IndexRowMin{0.f};                  // index-table row origin (MFT: global y min; unused for ITS phi-z)
  float IndexRowMax{0.f};                  // index-table row span end (MFT: global y max; 0 => TwoPI for ITS)
  std::vector<float> LayerRadii = {2.33959f, 3.14076f, 3.91924f, 19.6213f, 24.5597f, 34.388f, 39.3329f};
  std::vector<float> LayerxX0 = {5.e-3f, 5.e-3f, 5.e-3f, 1.e-2f, 1.e-2f, 1.e-2f, 1.e-2f};
  std::vector<float> LayerResolution = {5.e-4f, 5.e-4f, 5.e-4f, 5.e-4f, 5.e-4f, 5.e-4f, 5.e-4f};
  std::vector<float> SystError2Row = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f}; // systematic error^2 along local row (ALPIDE X) per layer
  std::vector<float> SystError2Col = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f}; // systematic error^2 along local column (ALPIDE Z) per layer
  int ColBins{256}; // ITS: ZBins
  int RowBins{128}; // ITS: PhiBins
  bool UseDiamond = false;
  float Diamond[3] = {0.f, 0.f, 0.f};
  float DiamondCov[6] = {25.e-6f, 0.f, 0.f, 25.e-6f, 0.f, 36.f};

  /// General parameters
  int MinTrackLength = 7;
  int MaxHoles = 0;
  tracking::LayerMask HoleLayerMask = 0;
  float NSigmaCut = 5;
  float PVres = 1.e-2f;
  /// Trackleting cuts
  float TrackletMinPt = 0.3f;
  /// Cell finding cuts
  float CellDeltaTanLambdaSigma = 0.007f;
  float CellRoadRCut = 0.05f; // MFT: max distance to seed line (classic ROADclsRCut)
  float TrackletMinAbsX = 0.f; // MFT: reject clusters/tracks with |x| below this (cm); 0 = disabled
  /// Fitter parameters
  o2::base::PropagatorImpl<float>::MatCorrType CorrType = o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrNONE;
  float MaxChi2ClusterAttachment = 60.f;
  float MaxChi2NDF = 30.f;
  int ReseedIfShorter = 6; // reseed for the final fit track with the length shorter than this
  std::vector<float> MinPt = {0.f, 0.f, 0.f, 0.f};
  tracking::LayerMask StartLayerMask = 0x7F;
  bool RepeatRefitOut = false;   // repeat outward refit using inward refit as a seed
  bool ShiftRefToCluster = true; // TrackFit: after update shift the linearization reference to cluster
  bool PerPrimaryVertexProcessing = false;
  bool SaveTimeBenchmarks = false;
  bool DoUPCIteration = false;
  bool FataliseUponFailure = true;
  bool CreateArtefactLabels{false};
  bool PrintMemory = false; // print allocator usage in epilog report
  size_t MaxMemory = std::numeric_limits<size_t>::max();
  bool DropTFUponFailure = false;

  // Selections on tracks sharing clusters
  bool AllowSharingFirstCluster = false;
  float SharedClusterMaxDeltaPhi = 0.05f; // For tracks sharing clusters, maximum allowed delta phi at the cluster position
  float SharedClusterMaxDeltaEta = 0.03f; // For tracks sharing clusters, maximum allowed delta eta at the cluster position
  bool SharedClusterOppositeSign = false; // For tracks sharing clusters, require opposite sign of the tracklets
  int SharedMaxClusters = 0;              // Maximal allowed shared clusters (excluding first cluster)
};

/// Reset tracking parameters to detector geometry defaults (ITS: struct defaults; MFT: MFTTracking/Constants.h).
void resetDetectorDefaults(TrackingParameters& params, o2::detectors::DetID::ID detId);

namespace TrackingMode
{
enum Type : int8_t {
  Unset = -1,
  Sync = 0,
  Async = 1,
  Cosmics = 2,
  Off = 3,
};

Type fromString(std::string_view str);
std::string toString(Type mode);
std::vector<TrackingParameters> getTrackingParameters(o2::detectors::DetID::ID detId, Type mode);

} // namespace TrackingMode

struct VertexingParameters {
  std::string asString() const;

  IterationSteps PassFlags{IterationStep::FirstPass, IterationStep::ResetVertices};
  std::vector<float> LayerZ = {16.333f + 1, 16.333f + 1, 16.333f + 1, 42.140f + 1, 42.140f + 1, 73.745f + 1, 73.745f + 1};
  std::vector<float> LayerRadii = {2.33959f, 3.14076f, 3.91924f, 19.6213f, 24.5597f, 34.388f, 39.3329f};
  int vertPerRofThreshold = 0; // Maximum number of vertices per ROF to trigger second a round
  int ColBins = 1;
  int RowBins = 128;
  float zCut = -1.f;
  float phiCut = -1.f;
  float pairCut = -1.f;
  float clusterCut = -1.f;
  float coarseZWindow = -1.f;
  float seedDedupZCut = -1.f;
  float refitDedupZCut = -1.f;
  float duplicateZCut = -1.f;
  float finalSelectionZCut = -1.f;
  float duplicateDistance2Cut = -1.f;
  float tanLambdaCut = -1.f;
  float NSigmaCut = -1;
  float maxZPositionAllowed = -1.f;
  int clusterContributorsCut = -1;
  int suppressLowMultDebris = -1;
  int seedMemberRadiusTime = -1;
  int seedMemberRadiusZ = -1;
  int maxTrackletsPerCluster = -1;
  int phiSpan = -1;
  int zSpan = -1;
  bool SaveTimeBenchmarks = false;

  bool useTruthSeeding = false; // overwrite found vertices with MC events

  int nThreads = 1;
  bool PrintMemory = false; // print allocator usage in epilog report
  size_t MaxMemory = std::numeric_limits<size_t>::max();
  bool DropTFUponFailure = false;
};

} // namespace o2::itsmft

#endif /* ALICEO2_ITSMFT_TRACKING_CONFIGURATION_H_ */

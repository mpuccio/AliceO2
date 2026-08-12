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

#ifndef GPUCA_GPUCODE
#include <cmath>
#include <gsl/span>
#include "ITSMFTTracking/SurfaceDescriptor.h"
#endif

#ifndef GPUCA_GPUCODE_DEVICE
#include <limits>
#include <string>
#include <string_view>
#include <vector>
#endif

#include "CommonUtils/EnumFlags.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/SurfaceMask.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITStracking/TrackingConfigParam.h"

namespace o2::itsmft
{

inline constexpr int ClustersPerCell = 3;

// Steering of dedicated steps in an iteration
enum class IterationStep : uint16_t {
  FirstPass = 0,
  RebuildClusterLUT = 1,
  UseUPCMask = 2,
  SelectUPCVertices = 3,
  ResetVertices = 4,
  SkipROFsAboveThreshold = 5,
  MarkVerticesAsUPC = 6,
  TrackFollowerTop = 7,
  TrackFollowerBot = 8,
};
using IterationSteps = o2::utils::EnumFlags<IterationStep>;

static_assert(sizeof(IterationStep) == sizeof(uint16_t));
static_assert(sizeof(IterationSteps) == sizeof(uint16_t));
static_assert(static_cast<uint16_t>(IterationStep::FirstPass) == 0);
static_assert(static_cast<uint16_t>(IterationStep::RebuildClusterLUT) == 1);
static_assert(static_cast<uint16_t>(IterationStep::UseUPCMask) == 2);
static_assert(static_cast<uint16_t>(IterationStep::SelectUPCVertices) == 3);
static_assert(static_cast<uint16_t>(IterationStep::ResetVertices) == 4);
static_assert(static_cast<uint16_t>(IterationStep::SkipROFsAboveThreshold) == 5);
static_assert(static_cast<uint16_t>(IterationStep::MarkVerticesAsUPC) == 6);
static_assert(static_cast<uint16_t>(IterationStep::TrackFollowerTop) == 7);
static_assert(static_cast<uint16_t>(IterationStep::TrackFollowerBot) == 8);

struct TrackingParameters {
  int CellMinimumLevel() const noexcept
  {
    const int minClusters = MinTrackLength - (MaxHoles > 0 ? MaxHoles : 0);
    const int effectiveMinClusters = minClusters > ClustersPerCell ? minClusters : ClustersPerCell;
    return effectiveMinClusters - ClustersPerCell + 1;
  }
  // Adapter/frozen-ITS compatibility accessors.  The common runtime tracker
  // derives its road start level from SurfaceTrackingScratch; these remain
  // only for the frozen ITStracking and GPU-facing consumers of this shared
  // configuration record.
  int NeighboursPerRoad() const noexcept { return NLayers - 3; }
  int CellsPerRoad() const noexcept { return NLayers - 2; }
  int TrackletsPerRoad() const noexcept { return NLayers - 1; }
  std::string asString() const;

  IterationSteps PassFlags{IterationStep::FirstPass, IterationStep::RebuildClusterLUT};
  int NLayers = tracking::ITSNLayers;
  std::vector<uint32_t> AddTimeError = {0, 0, 0, 0, 0, 0, 0};
  std::vector<float> LayerZ = {16.333f + 1, 16.333f + 1, 16.333f + 1, 42.140f + 1, 42.140f + 1, 73.745f + 1, 73.745f + 1};
  std::vector<float> LayerColHalfExtent{}; // index-table column half extent (ITS: z, MFT: global x); falls back to LayerZ
  float IndexRowMin{0.f};                  // index-table row origin (MFT: global y min; unused for ITS phi-z)
  float IndexRowMax{0.f};                  // index-table row span end (MFT: global y max; 0 => TwoPI for ITS)
  std::vector<float> LayerRadii = {2.33959f, 3.14076f, 3.91924f, 19.6213f, 24.5597f, 34.388f, 39.3329f};
  std::vector<float> LayerxX0{tracking::kNominalITSLayerX0.begin(), tracking::kNominalITSLayerX0.end()};
  std::vector<float> LayerResolution = {5.e-4f, 5.e-4f, 5.e-4f, 5.e-4f, 5.e-4f, 5.e-4f, 5.e-4f};
  std::vector<float> SystError2Row = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f}; // systematic error^2 along local row (ALPIDE X) per layer
  std::vector<float> SystError2Col = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f}; // systematic error^2 along local column (ALPIDE Z) per layer
  int ColBins{256};                                                       // ITS: ZBins
  int RowBins{128};                                                       // ITS: PhiBins
  bool UseDiamond = false;
  float Diamond[3] = {0.f, 0.f, 0.f};
  float DiamondCov[6] = {25.e-6f, 0.f, 0.f, 25.e-6f, 0.f, 36.f};

  /// General parameters
  int MinTrackLength = 7;
  int MaxHoles = 0;
  tracking::LayerMask HoleLayerMask = 0;
  // Reserved compatibility storage; adapters reject non-empty values because
  // the common CA does not consume these declarations.
  tracking::LayerMask InactiveLayerMask = 0;
  tracking::LayerMask SeedingLayers = 0;
  float NSigmaCut = 5;
  float PVres = 1.e-2f;
  /// Trackleting cuts
  float TrackletMinPt = 0.3f;
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
  // Reserved compatibility storage; top/bottom follower execution is not part
  // of the common tracker.
  float TrackFollowerNSigmaCutZ = 1.f;
  float TrackFollowerNSigmaCutPhi = 1.f;
  int TrackFollowerMaxHypotheses = 1;
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

#ifndef GPUCA_GPUCODE

inline bool isRecognizedMatCorrType(o2::base::PropagatorF::MatCorrType corrType) noexcept
{
  return corrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE ||
         corrType == o2::base::PropagatorF::MatCorrType::USEMatCorrTGeo ||
         corrType == o2::base::PropagatorF::MatCorrType::USEMatCorrLUT;
}

struct AttachHitConfigView {
  gsl::span<const tracking::NominalSurfaceMaterial> layerMaterial;
  o2::base::PropagatorF::MatCorrType corrType{o2::base::PropagatorF::MatCorrType::USEMatCorrNONE};

  bool isValid(size_t expectedLayers) const noexcept
  {
    if (layerMaterial.size() < expectedLayers || !isRecognizedMatCorrType(corrType)) {
      return false;
    }
    for (size_t layer = 0; layer < expectedLayers; ++layer) {
      const auto& material = layerMaterial[layer];
      if (!std::isfinite(material.xOverX0) || material.xOverX0 < 0.f ||
          !std::isfinite(material.arealDensityGPerCm2) || material.arealDensityGPerCm2 < 0.f) {
        return false;
      }
    }
    return true;
  }
};

inline AttachHitConfigView bindAttachHitConfig(gsl::span<const tracking::NominalSurfaceMaterial> layerMaterial,
                                               const TrackingParameters& params) noexcept
{
  return {layerMaterial, params.CorrType};
}

#endif

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

namespace o2::itsmft::tracking
{

/// MFT uses o2::itsmft::TrackerParamConfig; ITS production params stay in O2::ITStracking.
template <o2::detectors::DetID::ID DetId>
struct TrackerParamRef;

template <>
struct TrackerParamRef<o2::detectors::DetID::MFT> {
  using Type = o2::itsmft::TrackerParamConfig<o2::detectors::DetID::MFT>;
  static const Type& get() { return Type::Instance(); }
  static constexpr int nLayers() { return Type::getNLayers(); }
};

template <>
struct TrackerParamRef<o2::detectors::DetID::ITS> {
  using Type = o2::its::TrackerParamConfig;
  static const Type& get() { return Type::Instance(); }
  static constexpr int nLayers() { return ITSNLayers; }
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_CONFIGURATION_H_ */

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSCommonTracking/CommonTrackingParameters.h"

#include <algorithm>
#include <iterator>

namespace o2::its::commontracking
{
namespace
{
o2::itsmft::IterationSteps translateIterationSteps(const o2::its::IterationSteps& legacy)
{
  o2::itsmft::IterationSteps translated;
  translated.reset();

  const auto copyByName = [&legacy, &translated](o2::its::IterationStep legacyStep, o2::itsmft::IterationStep commonStep) {
    if (legacy[legacyStep]) {
      translated.set(commonStep);
    }
  };
  copyByName(o2::its::IterationStep::FirstPass, o2::itsmft::IterationStep::FirstPass);
  copyByName(o2::its::IterationStep::RebuildClusterLUT, o2::itsmft::IterationStep::RebuildClusterLUT);
  copyByName(o2::its::IterationStep::UseUPCMask, o2::itsmft::IterationStep::UseUPCMask);
  copyByName(o2::its::IterationStep::SelectUPCVertices, o2::itsmft::IterationStep::SelectUPCVertices);
  copyByName(o2::its::IterationStep::ResetVertices, o2::itsmft::IterationStep::ResetVertices);
  copyByName(o2::its::IterationStep::SkipROFsAboveThreshold, o2::itsmft::IterationStep::SkipROFsAboveThreshold);
  copyByName(o2::its::IterationStep::MarkVerticesAsUPC, o2::itsmft::IterationStep::MarkVerticesAsUPC);
  copyByName(o2::its::IterationStep::TrackFollowerTop, o2::itsmft::IterationStep::TrackFollowerTop);
  copyByName(o2::its::IterationStep::TrackFollowerBot, o2::itsmft::IterationStep::TrackFollowerBot);

  return translated;
}
} // namespace

o2::itsmft::TrackingParameters translateTrackingParameters(const o2::its::TrackingParameters& legacy)
{
  o2::itsmft::TrackingParameters translated;
  translated.PassFlags = translateIterationSteps(legacy.PassFlags);
  translated.NLayers = legacy.NLayers;
  translated.AddTimeError = legacy.AddTimeError;
  translated.LayerZ = legacy.LayerZ;
  translated.LayerRadii = legacy.LayerRadii;
  translated.LayerxX0 = legacy.LayerxX0;
  translated.LayerResolution = legacy.LayerResolution;
  translated.SystError2Row = legacy.SystErrorY2;
  translated.SystError2Col = legacy.SystErrorZ2;
  translated.ColBins = legacy.ZBins;
  translated.RowBins = legacy.PhiBins;
  translated.UseDiamond = legacy.UseDiamond;
  std::copy(std::begin(legacy.Diamond), std::end(legacy.Diamond), std::begin(translated.Diamond));
  std::copy(std::begin(legacy.DiamondCov), std::end(legacy.DiamondCov), std::begin(translated.DiamondCov));

  translated.MinTrackLength = legacy.MinTrackLength;
  translated.MaxHoles = legacy.MaxHoles;
  translated.HoleLayerMask = legacy.HoleLayerMask.value();
  translated.InactiveLayerMask = legacy.InactiveLayerMask.value();
  translated.SeedingLayers = legacy.SeedingLayers.value();
  translated.NSigmaCut = legacy.NSigmaCut;
  translated.PVres = legacy.PVres;
  translated.TrackletMinPt = legacy.TrackletMinPt;
  translated.CorrType = legacy.CorrType;
  translated.MaxChi2ClusterAttachment = legacy.MaxChi2ClusterAttachment;
  translated.MaxChi2NDF = legacy.MaxChi2NDF;
  translated.ReseedIfShorter = legacy.ReseedIfShorter;
  translated.MinPt = legacy.MinPt;
  translated.StartLayerMask = legacy.StartLayerMask.value();
  translated.RepeatRefitOut = legacy.RepeatRefitOut;
  translated.ShiftRefToCluster = legacy.ShiftRefToCluster;
  translated.PerPrimaryVertexProcessing = legacy.PerPrimaryVertexProcessing;
  translated.SaveTimeBenchmarks = legacy.SaveTimeBenchmarks;
  translated.DoUPCIteration = legacy.DoUPCIteration;
  translated.FataliseUponFailure = legacy.FataliseUponFailure;
  translated.CreateArtefactLabels = legacy.CreateArtefactLabels;
  translated.TrackFollowerNSigmaCutZ = legacy.TrackFollowerNSigmaCutZ;
  translated.TrackFollowerNSigmaCutPhi = legacy.TrackFollowerNSigmaCutPhi;
  translated.TrackFollowerMaxHypotheses = legacy.TrackFollowerMaxHypotheses;
  translated.PrintMemory = legacy.PrintMemory;
  translated.MaxMemory = legacy.MaxMemory;
  translated.DropTFUponFailure = legacy.DropTFUponFailure;

  translated.AllowSharingFirstCluster = legacy.AllowSharingFirstCluster;
  translated.SharedClusterMaxDeltaPhi = legacy.SharedClusterMaxDeltaPhi;
  translated.SharedClusterMaxDeltaEta = legacy.SharedClusterMaxDeltaEta;
  translated.SharedClusterOppositeSign = legacy.SharedClusterOppositeSign;
  translated.SharedMaxClusters = legacy.SharedMaxClusters;
  return translated;
}

std::vector<o2::itsmft::TrackingParameters> makeTrackingParameters(o2::its::TrackingMode::Type mode)
{
  const auto legacy = o2::its::TrackingMode::getTrackingParameters(mode);
  std::vector<o2::itsmft::TrackingParameters> translated;
  translated.reserve(legacy.size());
  std::transform(legacy.begin(), legacy.end(), std::back_inserter(translated), translateTrackingParameters);
  return translated;
}

} // namespace o2::its::commontracking

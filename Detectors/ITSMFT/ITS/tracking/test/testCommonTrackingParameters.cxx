// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITS common tracking parameter adapter
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <limits>
#include <string>
#include <system_error>
#include <vector>

#include <sys/resource.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "DetectorsBase/Propagator.h"
#include "ITSCommonTracking/CommonTrackingParameters.h"
#include "ITStracking/TrackingConfigParam.h"

namespace
{
using LegacyParameters = o2::its::TrackingParameters;
using CommonParameters = o2::itsmft::TrackingParameters;

o2::its::TrackerParamConfig& trackerConfig()
{
  return const_cast<o2::its::TrackerParamConfig&>(o2::its::TrackerParamConfig::Instance());
}

class GlobalTrackingStateGuard
{
 public:
  GlobalTrackingStateGuard()
    : mNominalBz{o2::base::Propagator::Instance(true)->getNominalBz()}
  {
    auto& tc = trackerConfig();
#define PRESERVE_CONFIG_FIELD(field) preserve(tc.field)
    PRESERVE_CONFIG_FIELD(useMatCorrTGeo);
    PRESERVE_CONFIG_FIELD(useFastMaterial);
    preserveArray(tc.addTimeError);
    preserveArray(tc.minTrackLgtIter);
    preserveArray(tc.startLayerMask);
    preserveArray(tc.maxHolesIter);
    preserveArray(tc.minPtIterLgt);
    preserveArray(tc.sysErrY2);
    preserveArray(tc.sysErrZ2);
    PRESERVE_CONFIG_FIELD(maxChi2ClusterAttachment);
    PRESERVE_CONFIG_FIELD(maxChi2NDF);
    PRESERVE_CONFIG_FIELD(nSigmaCut);
    PRESERVE_CONFIG_FIELD(deltaTanLres);
    PRESERVE_CONFIG_FIELD(minPt);
    PRESERVE_CONFIG_FIELD(pvRes);
    PRESERVE_CONFIG_FIELD(LUTbinsPhi);
    PRESERVE_CONFIG_FIELD(LUTbinsZ);
    preserveArray(tc.diamondPos);
    PRESERVE_CONFIG_FIELD(useDiamond);
    PRESERVE_CONFIG_FIELD(perPrimaryVertexProcessing);
    PRESERVE_CONFIG_FIELD(saveTimeBenchmarks);
    PRESERVE_CONFIG_FIELD(overrideBeamEstimation);
    PRESERVE_CONFIG_FIELD(trackingMode);
    PRESERVE_CONFIG_FIELD(doUPCIteration);
    PRESERVE_CONFIG_FIELD(nIterations);
    PRESERVE_CONFIG_FIELD(reseedIfShorter);
    PRESERVE_CONFIG_FIELD(shiftRefToCluster);
    PRESERVE_CONFIG_FIELD(repeatRefitOut);
    PRESERVE_CONFIG_FIELD(createArtefactLabels);
    preserveArray(tc.trackFollowerTop);
    preserveArray(tc.trackFollowerBot);
    PRESERVE_CONFIG_FIELD(trackFollowerNSigmaCutZ);
    PRESERVE_CONFIG_FIELD(trackFollowerNSigmaCutPhi);
    PRESERVE_CONFIG_FIELD(trackFollowerMaxHypotheses);
    PRESERVE_CONFIG_FIELD(nThreads);
    PRESERVE_CONFIG_FIELD(printMemory);
    PRESERVE_CONFIG_FIELD(maxMemory);
    PRESERVE_CONFIG_FIELD(dropTFUponFailure);
    PRESERVE_CONFIG_FIELD(fataliseUponFailure);
    PRESERVE_CONFIG_FIELD(allowSharingFirstCluster);
    PRESERVE_CONFIG_FIELD(sharedClusterMaxDeltaPhi);
    PRESERVE_CONFIG_FIELD(sharedClusterMaxDeltaEta);
    PRESERVE_CONFIG_FIELD(sharedClusterOppositeSign);
#undef PRESERVE_CONFIG_FIELD
  }

  ~GlobalTrackingStateGuard() noexcept
  {
    for (auto restore = mRestorers.rbegin(); restore != mRestorers.rend(); ++restore) {
      (*restore)();
    }
    o2::base::Propagator::Instance()->setNominalBz(mNominalBz);
  }

  GlobalTrackingStateGuard(const GlobalTrackingStateGuard&) = delete;
  GlobalTrackingStateGuard& operator=(const GlobalTrackingStateGuard&) = delete;

 private:
  template <typename T>
  void preserve(T& value)
  {
    const T original = value;
    mRestorers.emplace_back([&value, original] { value = original; });
  }

  template <typename T, size_t N>
  void preserveArray(T (&values)[N])
  {
    for (auto& value : values) {
      preserve(value);
    }
  }

  std::vector<std::function<void()>> mRestorers;
  float mNominalBz;
};

template <typename LegacyVector, typename CommonVector>
void checkVector(const LegacyVector& legacy, const CommonVector& common)
{
  BOOST_REQUIRE_EQUAL(legacy.size(), common.size());
  BOOST_CHECK_EQUAL_COLLECTIONS(legacy.begin(), legacy.end(), common.begin(), common.end());
}

void checkFlags(const LegacyParameters& legacy, const CommonParameters& common)
{
  const auto check = [&legacy, &common](o2::its::IterationStep legacyStep, o2::itsmft::IterationStep commonStep) {
    BOOST_CHECK_EQUAL(legacy.PassFlags[legacyStep], common.PassFlags[commonStep]);
  };
  check(o2::its::IterationStep::FirstPass, o2::itsmft::IterationStep::FirstPass);
  check(o2::its::IterationStep::RebuildClusterLUT, o2::itsmft::IterationStep::RebuildClusterLUT);
  check(o2::its::IterationStep::UseUPCMask, o2::itsmft::IterationStep::UseUPCMask);
  check(o2::its::IterationStep::SelectUPCVertices, o2::itsmft::IterationStep::SelectUPCVertices);
  check(o2::its::IterationStep::ResetVertices, o2::itsmft::IterationStep::ResetVertices);
  check(o2::its::IterationStep::SkipROFsAboveThreshold, o2::itsmft::IterationStep::SkipROFsAboveThreshold);
  check(o2::its::IterationStep::MarkVerticesAsUPC, o2::itsmft::IterationStep::MarkVerticesAsUPC);
  check(o2::its::IterationStep::TrackFollowerTop, o2::itsmft::IterationStep::TrackFollowerTop);
  check(o2::its::IterationStep::TrackFollowerBot, o2::itsmft::IterationStep::TrackFollowerBot);
}

void checkParity(const LegacyParameters& legacy, const CommonParameters& common)
{
  checkFlags(legacy, common);
  BOOST_CHECK_EQUAL(legacy.NLayers, common.NLayers);
  checkVector(legacy.AddTimeError, common.AddTimeError);
  checkVector(legacy.LayerZ, common.LayerZ);
  checkVector(legacy.LayerRadii, common.LayerRadii);
  checkVector(legacy.LayerxX0, common.LayerxX0);
  checkVector(legacy.LayerResolution, common.LayerResolution);
  checkVector(legacy.SystErrorY2, common.SystError2Row);
  checkVector(legacy.SystErrorZ2, common.SystError2Col);
  BOOST_CHECK_EQUAL(legacy.ZBins, common.ColBins);
  BOOST_CHECK_EQUAL(legacy.PhiBins, common.RowBins);
  BOOST_CHECK_EQUAL(legacy.UseDiamond, common.UseDiamond);
  BOOST_CHECK_EQUAL_COLLECTIONS(std::begin(legacy.Diamond), std::end(legacy.Diamond), std::begin(common.Diamond), std::end(common.Diamond));
  BOOST_CHECK_EQUAL_COLLECTIONS(std::begin(legacy.DiamondCov), std::end(legacy.DiamondCov), std::begin(common.DiamondCov), std::end(common.DiamondCov));

  BOOST_CHECK_EQUAL(legacy.MinTrackLength, common.MinTrackLength);
  BOOST_CHECK_EQUAL(legacy.MaxHoles, common.MaxHoles);
  BOOST_CHECK_EQUAL(legacy.InactiveLayerMask.value(), common.InactiveLayerMask.value());
  BOOST_CHECK_EQUAL(legacy.SeedingLayers.value(), common.SeedingLayers.value());
  BOOST_CHECK_EQUAL(legacy.getActiveLayerMask().value(), common.getActiveLayerMask().value());
  BOOST_CHECK_EQUAL(legacy.getSeedingLayerMask().value(), common.getSeedingLayerMask().value());
  const auto configuredLayerSpan = o2::itsmft::tracking::LayerMask::span(0, legacy.NLayers - 1);
  BOOST_CHECK_EQUAL(legacy.getNonSeedingLayerMask().value() & configuredLayerSpan.value(),
                    common.getNonSeedingLayerMask().value());
  BOOST_CHECK_EQUAL(legacy.getNSeedingLayers(), common.getNSeedingLayers());
  BOOST_CHECK_EQUAL(legacy.getMinSeedingClusters(), common.getMinSeedingClusters());
  BOOST_CHECK_EQUAL(legacy.CellMinimumLevel(), common.CellMinimumLevel());
  BOOST_CHECK_EQUAL(legacy.NeighboursPerRoad(), common.NeighboursPerRoad());
  BOOST_CHECK_EQUAL(legacy.CellsPerRoad(), common.CellsPerRoad());
  BOOST_CHECK_EQUAL(legacy.TrackletsPerRoad(), common.TrackletsPerRoad());
  BOOST_CHECK_EQUAL(legacy.NSigmaCut, common.NSigmaCut);
  BOOST_CHECK_EQUAL(legacy.PVres, common.PVres);
  BOOST_CHECK_EQUAL(legacy.TrackletMinPt, common.TrackletMinPt);
  BOOST_CHECK(legacy.CorrType == common.CorrType);
  BOOST_CHECK_EQUAL(legacy.MaxChi2ClusterAttachment, common.MaxChi2ClusterAttachment);
  BOOST_CHECK_EQUAL(legacy.MaxChi2NDF, common.MaxChi2NDF);
  BOOST_CHECK_EQUAL(legacy.ReseedIfShorter, common.ReseedIfShorter);
  checkVector(legacy.MinPt, common.MinPt);
  BOOST_CHECK_EQUAL(legacy.StartLayerMask.value(), common.StartLayerMask.value());
  BOOST_CHECK_EQUAL(legacy.RepeatRefitOut, common.RepeatRefitOut);
  BOOST_CHECK_EQUAL(legacy.ShiftRefToCluster, common.ShiftRefToCluster);
  BOOST_CHECK_EQUAL(legacy.PerPrimaryVertexProcessing, common.PerPrimaryVertexProcessing);
  BOOST_CHECK_EQUAL(legacy.SaveTimeBenchmarks, common.SaveTimeBenchmarks);
  BOOST_CHECK_EQUAL(legacy.DoUPCIteration, common.DoUPCIteration);
  BOOST_CHECK_EQUAL(legacy.FataliseUponFailure, common.FataliseUponFailure);
  BOOST_CHECK_EQUAL(legacy.CreateArtefactLabels, common.CreateArtefactLabels);
  BOOST_CHECK_EQUAL(legacy.TrackFollowerNSigmaCutZ, common.TrackFollowerNSigmaCutZ);
  BOOST_CHECK_EQUAL(legacy.TrackFollowerNSigmaCutPhi, common.TrackFollowerNSigmaCutPhi);
  BOOST_CHECK_EQUAL(legacy.TrackFollowerMaxHypotheses, common.TrackFollowerMaxHypotheses);
  BOOST_CHECK_EQUAL(legacy.PrintMemory, common.PrintMemory);
  BOOST_CHECK_EQUAL(legacy.MaxMemory, common.MaxMemory);
  BOOST_CHECK_EQUAL(legacy.DropTFUponFailure, common.DropTFUponFailure);
  BOOST_CHECK_EQUAL(legacy.AllowSharingFirstCluster, common.AllowSharingFirstCluster);
  BOOST_CHECK_EQUAL(legacy.SharedClusterMaxDeltaPhi, common.SharedClusterMaxDeltaPhi);
  BOOST_CHECK_EQUAL(legacy.SharedClusterMaxDeltaEta, common.SharedClusterMaxDeltaEta);
  BOOST_CHECK_EQUAL(legacy.SharedClusterOppositeSign, common.SharedClusterOppositeSign);
  BOOST_CHECK_EQUAL(legacy.SharedMaxClusters, common.SharedMaxClusters);

  BOOST_CHECK(common.LayerColHalfExtent.empty());
  BOOST_CHECK_EQUAL(common.IndexRowMin, 0.f);
  BOOST_CHECK_EQUAL(common.IndexRowMax, 0.f);
}

void checkModeParity(o2::its::TrackingMode::Type mode)
{
  const auto legacy = o2::its::TrackingMode::getTrackingParameters(mode);
  const auto common = o2::its::commontracking::makeTrackingParameters(mode);
  BOOST_REQUIRE_EQUAL(legacy.size(), common.size());
  for (size_t i = 0; i < legacy.size(); ++i) {
    checkParity(legacy[i], common[i]);
  }
}

void resetGlobals(float bz = 5.0066791f)
{
  auto& tc = trackerConfig();
  tc.useMatCorrTGeo = false;
  tc.useFastMaterial = false;
  std::fill(std::begin(tc.addTimeError), std::end(tc.addTimeError), 0);
  std::fill(std::begin(tc.minTrackLgtIter), std::end(tc.minTrackLgtIter), 0);
  std::fill(std::begin(tc.startLayerMask), std::end(tc.startLayerMask), 0);
  std::fill(std::begin(tc.maxHolesIter), std::end(tc.maxHolesIter), 0);
  std::fill(std::begin(tc.minPtIterLgt), std::end(tc.minPtIterLgt), 0.f);
  std::fill(std::begin(tc.sysErrY2), std::end(tc.sysErrY2), 0.f);
  std::fill(std::begin(tc.sysErrZ2), std::end(tc.sysErrZ2), 0.f);
  tc.maxChi2ClusterAttachment = -1.f;
  tc.maxChi2NDF = -1.f;
  tc.nSigmaCut = -1.f;
  tc.deltaTanLres = -1.f;
  tc.minPt = -1.f;
  tc.pvRes = -1.f;
  tc.LUTbinsPhi = -1;
  tc.LUTbinsZ = -1;
  std::fill(std::begin(tc.diamondPos), std::end(tc.diamondPos), 0.f);
  tc.useDiamond = false;
  tc.perPrimaryVertexProcessing = false;
  tc.saveTimeBenchmarks = false;
  tc.overrideBeamEstimation = false;
  tc.trackingMode = -1;
  tc.doUPCIteration = false;
  tc.nIterations = o2::its::constants::MaxIter;
  tc.reseedIfShorter = 6;
  tc.shiftRefToCluster = true;
  tc.repeatRefitOut = false;
  tc.createArtefactLabels = false;
  std::fill(std::begin(tc.trackFollowerTop), std::end(tc.trackFollowerTop), false);
  std::fill(std::begin(tc.trackFollowerBot), std::end(tc.trackFollowerBot), false);
  tc.trackFollowerNSigmaCutZ = 1.f;
  tc.trackFollowerNSigmaCutPhi = 1.f;
  tc.trackFollowerMaxHypotheses = 1;
  tc.nThreads = 1;
  tc.printMemory = false;
  tc.maxMemory = std::numeric_limits<size_t>::max();
  tc.dropTFUponFailure = false;
  tc.fataliseUponFailure = true;
  tc.allowSharingFirstCluster = false;
  tc.sharedClusterMaxDeltaPhi = 0.05f;
  tc.sharedClusterMaxDeltaEta = 0.03f;
  tc.sharedClusterOppositeSign = false;
  o2::base::Propagator::Instance()->setNominalBz(bz);
}

void setSentinelOverrides()
{
  auto& tc = trackerConfig();
  tc.useMatCorrTGeo = true;
  for (int layer = 0; layer < 7; ++layer) {
    tc.addTimeError[layer] = 11 + layer;
    tc.sysErrY2[layer] = 0.001f * (layer + 1);
    tc.sysErrZ2[layer] = 0.002f * (layer + 1);
  }
  for (int iteration = 0; iteration < o2::its::constants::MaxIter; ++iteration) {
    tc.minTrackLgtIter[iteration] = 4 + iteration;
    tc.startLayerMask[iteration] = static_cast<uint8_t>(1u << iteration);
    tc.maxHolesIter[iteration] = iteration;
    tc.trackFollowerTop[iteration] = (iteration % 2) == 0;
    tc.trackFollowerBot[iteration] = (iteration % 2) != 0;
    for (int slot = 0; slot < 4; ++slot) {
      tc.minPtIterLgt[iteration * 4 + slot] = 0.01f * (1 + iteration * 4 + slot);
    }
  }
  tc.maxChi2ClusterAttachment = 17.f;
  tc.maxChi2NDF = 19.f;
  tc.nSigmaCut = 1.3f;
  tc.deltaTanLres = 1.7f;
  tc.minPt = 1.9f;
  tc.pvRes = 0.023f;
  tc.LUTbinsPhi = 37;
  tc.LUTbinsZ = 73;
  tc.diamondPos[0] = 0.11f;
  tc.diamondPos[1] = -0.22f;
  tc.diamondPos[2] = 0.33f;
  tc.useDiamond = true;
  tc.perPrimaryVertexProcessing = true;
  tc.saveTimeBenchmarks = true;
  tc.doUPCIteration = true;
  tc.nIterations = 4;
  tc.reseedIfShorter = 5;
  tc.shiftRefToCluster = false;
  tc.repeatRefitOut = true;
  tc.createArtefactLabels = true;
  tc.trackFollowerNSigmaCutZ = 2.3f;
  tc.trackFollowerNSigmaCutPhi = 2.9f;
  tc.trackFollowerMaxHypotheses = 7;
  tc.printMemory = true;
  tc.maxMemory = 123456789u;
  tc.dropTFUponFailure = true;
  tc.fataliseUponFailure = false;
  tc.allowSharingFirstCluster = true;
  tc.sharedClusterMaxDeltaPhi = 0.071f;
  tc.sharedClusterMaxDeltaEta = 0.083f;
  tc.sharedClusterOppositeSign = true;
}

class TemporaryArtifactDirectory
{
 public:
  TemporaryArtifactDirectory()
  {
    std::array<char, 64> pattern{};
    const std::string prefix = "/private/tmp/o2-its-common-parameters-XXXXXX";
    std::copy(prefix.begin(), prefix.end(), pattern.begin());
    auto* created = mkdtemp(pattern.data());
    BOOST_REQUIRE(created != nullptr);
    mPath = created;
  }

  ~TemporaryArtifactDirectory()
  {
    std::error_code error;
    std::filesystem::remove_all(mPath, error);
  }

  const std::filesystem::path& path() const { return mPath; }

 private:
  std::filesystem::path mPath;
};

void checkInvalidMode(o2::its::TrackingMode::Type mode)
{
  TemporaryArtifactDirectory artifacts;
  const auto testExecutable = std::filesystem::canonical(boost::unit_test::framework::master_test_suite().argv[0]);
  const auto helper = testExecutable.parent_path().parent_path() / "bin" / "o2-its-tracking-invalid-common-tracking-mode";
  BOOST_REQUIRE(std::filesystem::is_regular_file(helper));
  const auto modeArgument = std::to_string(static_cast<int>(mode));
  const pid_t child = fork();
  BOOST_REQUIRE(child >= 0);
  if (child == 0) {
    struct rlimit noCore{
      0, 0};
    setrlimit(RLIMIT_CORE, &noCore);
    if (chdir(artifacts.path().c_str()) != 0) {
      _exit(120);
    }
    execl(helper.c_str(), helper.filename().c_str(), modeArgument.c_str(), nullptr);
    _exit(121);
  }

  int status = 0;
  BOOST_REQUIRE_EQUAL(waitpid(child, &status, 0), child);
  BOOST_CHECK(WIFSIGNALED(status) || (WIFEXITED(status) && WEXITSTATUS(status) != 0));
  // FairRoot may create a core_dump_* artifact even with RLIMIT_CORE=0. It is
  // confined to this fresh directory and removed by its RAII owner.
}
} // namespace

BOOST_AUTO_TEST_CASE(ArbitraryLegacyObjectIsTranslatedFieldByField)
{
  LegacyParameters legacy;
  legacy.PassFlags.set(o2::its::IterationStep::FirstPass);
  legacy.PassFlags.set(o2::its::IterationStep::UseUPCMask);
  legacy.PassFlags.set(o2::its::IterationStep::TrackFollowerTop);
  legacy.PassFlags.set(o2::its::IterationStep::TrackFollowerBot);
  legacy.NLayers = 6;
  legacy.AddTimeError = {1, 2, 3, 4, 5, 6, 7};
  legacy.LayerZ = {1.f, 2.f, 3.f};
  legacy.LayerRadii = {4.f, 5.f, 6.f};
  legacy.LayerxX0 = {0.01f, 0.02f, 0.03f};
  legacy.LayerResolution = {0.04f, 0.05f, 0.06f};
  legacy.SystErrorY2 = {0.07f, 0.08f, 0.09f};
  legacy.SystErrorZ2 = {0.10f, 0.11f, 0.12f};
  legacy.ZBins = 41;
  legacy.PhiBins = 43;
  legacy.UseDiamond = true;
  const std::array diamond{1.f, 2.f, 3.f};
  std::copy(diamond.begin(), diamond.end(), std::begin(legacy.Diamond));
  const std::array covariance{4.f, 5.f, 6.f, 7.f, 8.f, 9.f};
  std::copy(covariance.begin(), covariance.end(), std::begin(legacy.DiamondCov));
  legacy.MinTrackLength = 5;
  legacy.MaxHoles = 2;
  legacy.InactiveLayerMask = 0x24;
  legacy.SeedingLayers = 0x5a;
  legacy.NSigmaCut = 3.1f;
  legacy.PVres = 3.2f;
  legacy.TrackletMinPt = 3.3f;
  legacy.CorrType = o2::base::PropagatorF::MatCorrType::USEMatCorrTGeo;
  legacy.MaxChi2ClusterAttachment = 3.5f;
  legacy.MaxChi2NDF = 3.6f;
  legacy.ReseedIfShorter = 4;
  legacy.MinPt = {3.7f, 3.8f, 3.9f, 4.f};
  legacy.StartLayerMask = 0x45;
  legacy.RepeatRefitOut = true;
  legacy.ShiftRefToCluster = false;
  legacy.PerPrimaryVertexProcessing = true;
  legacy.SaveTimeBenchmarks = true;
  legacy.DoUPCIteration = true;
  legacy.FataliseUponFailure = false;
  legacy.CreateArtefactLabels = true;
  legacy.TrackFollowerNSigmaCutZ = 4.1f;
  legacy.TrackFollowerNSigmaCutPhi = 4.2f;
  legacy.TrackFollowerMaxHypotheses = 8;
  legacy.PrintMemory = true;
  legacy.MaxMemory = 987654u;
  legacy.DropTFUponFailure = true;
  legacy.AllowSharingFirstCluster = true;
  legacy.SharedClusterMaxDeltaPhi = 4.3f;
  legacy.SharedClusterMaxDeltaEta = 4.4f;
  legacy.SharedClusterOppositeSign = true;
  legacy.SharedMaxClusters = 2;

  checkParity(legacy, o2::its::commontracking::translateTrackingParameters(legacy));
}

BOOST_AUTO_TEST_CASE(ProductionModesAndUPCHaveFieldParity)
{
  GlobalTrackingStateGuard restore;
  resetGlobals();
  checkModeParity(o2::its::TrackingMode::Sync);
  checkModeParity(o2::its::TrackingMode::Async);
  checkModeParity(o2::its::TrackingMode::Cosmics);

  trackerConfig().doUPCIteration = true;
  checkModeParity(o2::its::TrackingMode::Async);
  BOOST_CHECK_EQUAL(o2::its::commontracking::makeTrackingParameters(o2::its::TrackingMode::Async).size(), 4u);
}

BOOST_AUTO_TEST_CASE(DistinctSentinelOverridesPreserveConstructionSemantics)
{
  GlobalTrackingStateGuard restore;
  resetGlobals();
  setSentinelOverrides();
  checkModeParity(o2::its::TrackingMode::Sync);
  checkModeParity(o2::its::TrackingMode::Async);
  checkModeParity(o2::its::TrackingMode::Cosmics);

  const auto async = o2::its::commontracking::makeTrackingParameters(o2::its::TrackingMode::Async);
  BOOST_REQUIRE_EQUAL(async.size(), 4u);
  for (size_t iteration = 0; iteration < async.size(); ++iteration) {
    BOOST_CHECK_EQUAL(async[iteration].StartLayerMask.value(), uint16_t{1} << iteration);
    BOOST_CHECK_EQUAL(async[iteration].MaxHoles, iteration);
  }
  BOOST_CHECK(async[0].PassFlags[o2::itsmft::IterationStep::TrackFollowerTop]);
  BOOST_CHECK(async[1].PassFlags[o2::itsmft::IterationStep::TrackFollowerBot]);

  const auto sync = o2::its::commontracking::makeTrackingParameters(o2::its::TrackingMode::Sync);
  BOOST_REQUIRE_EQUAL(sync.size(), 1u);
  BOOST_CHECK_EQUAL(sync[0].StartLayerMask.value(), 0x7fu);
  BOOST_CHECK_EQUAL(sync[0].MinTrackLength, 4);
  BOOST_CHECK_EQUAL(sync[0].MaxHoles, 0);
}

BOOST_AUTO_TEST_CASE(MaterialModesAndIterationTruncationHaveParity)
{
  GlobalTrackingStateGuard restore;
  resetGlobals();
  auto& tc = trackerConfig();

  tc.useMatCorrTGeo = false;
  tc.useFastMaterial = false;
  checkModeParity(o2::its::TrackingMode::Async);
  BOOST_CHECK(o2::its::commontracking::makeTrackingParameters(o2::its::TrackingMode::Async)[0].CorrType == o2::base::PropagatorF::MatCorrType::USEMatCorrLUT);

  tc.useFastMaterial = true;
  checkModeParity(o2::its::TrackingMode::Async);
  BOOST_CHECK(o2::its::commontracking::makeTrackingParameters(o2::its::TrackingMode::Async)[0].CorrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE);

  tc.useMatCorrTGeo = true;
  checkModeParity(o2::its::TrackingMode::Async);
  BOOST_CHECK(o2::its::commontracking::makeTrackingParameters(o2::its::TrackingMode::Async)[0].CorrType == o2::base::PropagatorF::MatCorrType::USEMatCorrTGeo);

  tc.doUPCIteration = true;
  tc.nIterations = 2;
  checkModeParity(o2::its::TrackingMode::Async);
  BOOST_CHECK_EQUAL(o2::its::commontracking::makeTrackingParameters(o2::its::TrackingMode::Async).size(), 2u);
}

BOOST_AUTO_TEST_CASE(MagneticFieldScalingIsCopiedWithoutRescaling)
{
  GlobalTrackingStateGuard restore;
  for (const float bz : {0.f, 2.50333955f, 5.0066791f, -5.0066791f}) {
    resetGlobals(bz);
    checkModeParity(o2::its::TrackingMode::Async);
    const auto translated = o2::its::commontracking::makeTrackingParameters(o2::its::TrackingMode::Async);
    const float factor = std::abs(bz) / 5.0066791f;
    const float trackletFactor = factor < 0.01f ? 1.f : factor;
    BOOST_CHECK_CLOSE(translated[0].TrackletMinPt, 0.3f * trackletFactor, 1.e-4f);
    BOOST_CHECK_CLOSE(translated[1].TrackletMinPt, 0.2f * trackletFactor, 1.e-4f);
    BOOST_CHECK_SMALL(translated[0].MinPt[0] - (1.f / 12.f) * factor, 1.e-7f);
  }
}

BOOST_AUTO_TEST_CASE(InvalidModesPreserveLegacyFatalBehaviorOutsideWorktree)
{
  GlobalTrackingStateGuard restore;
  resetGlobals();
  checkInvalidMode(o2::its::TrackingMode::Unset);
  checkInvalidMode(o2::its::TrackingMode::Off);
  checkInvalidMode(static_cast<o2::its::TrackingMode::Type>(99));
}

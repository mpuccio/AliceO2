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

#include "GenericTracking/TOF/Workflow.h"
#include "GenericTracking/TOF/Clusters.h"
#include "GenericTracking/TOF/Transport.h"
#include "GenericTracking/Workflow/ConditionsTRK.h"
#include "GenericTracking/Workflow/GeometryProvider.h"
#include "GenericTracking/Workflow/TrackerSpec.h"

#include "CommonConstants/LHCConstants.h"
#include "DPLUtils/MakeRootTreeWriterSpec.h"
#include "DataFormatsIOTOF/Cluster.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/Task.h"
#include "Framework/TimingInfo.h"
#include "IOTOFBase/GeometryTGeo.h"

#include <TFile.h>
#include <TList.h>
#include <TStreamerInfo.h>
#include <TTree.h>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>

namespace o2::gtrk::tof
{
using Det = detectors::TRKFT3;
using namespace o2::framework;

namespace
{
class MatcherTask : public Task
{
 public:
  explicit MatcherTask(std::shared_ptr<workflow::Conditions<Det, true>> conditions) : mConditions(std::move(conditions)) {}
  void init(InitContext& ic) final
  {
    mParameters.maxSpatialChi2 = ic.options().get<double>("tof-max-spatial-chi2");
    mParameters.maxTimeChi2 = ic.options().get<double>("tof-max-time-chi2");
    mParameters.timeMarginNS = ic.options().get<double>("tof-time-margin-ns");
    mParameters.maxStepCM = ic.options().get<double>("tof-max-step-cm");
    mParameters.maxPathCM = ic.options().get<double>("tof-max-path-cm");
    mSigmaNS = static_cast<float>(ic.options().get<double>("tof-time-resolution-ps") * 1e-3);
    if (!(mParameters.maxStepCM > 0.) || !(mParameters.maxPathCM > 0.) || !(mSigmaNS > 0.F)) {
      throw std::invalid_argument("TOF step, path and timing resolution must be positive");
    }
    mConditions->init(ic);
  }
  void finaliseCCDB(ConcreteDataMatcher& matcher, void* object) final { mConditions->finaliseCCDB(matcher, object); }
  void run(ProcessingContext& pc) final
  {
    mConditions->checkUpdates(pc);
    if (!mTracking) {
      mGeometry = workflow::makeGeometry<Det>(mConditions->measurementDepthsCM());
      mTracking = std::make_unique<TrackingGeometry<Det>>(*mGeometry);
      mTOF = o2::iotof::GeometryTGeo::Instance();
      mTOF->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
    }
    const auto input = pc.inputs().get<gsl::span<workflow::TrackRecord<Det>>>("tracks");
    const auto clusters = pc.inputs().get<gsl::span<o2::iotof::Cluster>>("tofClusters");
    const auto hits = convertClusters({clusters.data(), clusters.size()}, *mTOF, mSigmaNS);
    const auto* material = mConditions->materialLut();
    for (const auto& hit : hits) {
      if (std::hypot(hit.surface.origin[0], hit.surface.origin[1]) > material->getRMax() || std::abs(hit.surface.origin[2]) > material->getZMax()) {
        throw std::runtime_error("TOF sensor lies beyond TRK/Param/MatLUT: extend the LUT in radius and z to cover the TOF before matching");
      }
    }
    const auto orbit = pc.services().get<TimingInfo>().firstTForbit;
    const Timestamp start{static_cast<std::uint64_t>(orbit) * o2::constants::lhc::LHCMaxBunches, 0.F};
    std::vector<Track> tracks(input.size());
    for (std::size_t i = 0; i < input.size(); ++i) {
      const auto& record = input[i];
      if (record.status != PropagationStatus::kOk || record.state.surface >= Det::nSensors || record.outwardState.surface >= Det::nSensors) {
        continue;
      }
      tracks[i] = {.state = record.state, .outwardState = record.outwardState, .surface = mTracking->planeOf(record.state), .outwardSurface = mTracking->planeOf(record.outwardState), .boxStart = shifted(start, record.timePS * 1e-3), .boxWidthNS = record.widthPS * 1e-3, .valid = true};
    }
    Propagator propagator{.bzKiloGauss = static_cast<float>(workflow::nominalBzKiloGauss()),
                          .massGeV = static_cast<float>(recoParam<Det>().massHypothesisGeV)};
    const auto matches = match(tracks, hits, [&](const Track& track, const Hit& hit) { return predict(track, hit, propagator, *material, mParameters); }, mParameters);
    pc.outputs().snapshot(Output{workflow::gDataOrigin, "TRKTOFMATCH", 0}, matches);
    const auto outer = std::count_if(matches.begin(), matches.end(), [](const Match& m) { return m.outer >= 0; });
    const auto inner = std::count_if(matches.begin(), matches.end(), [](const Match& m) { return m.inner >= 0; });
    LOGP(info, "TOF: {} tracks, {} usable clusters ({} skipped), {} outer and {} inner matches", input.size(), hits.size(), clusters.size() - hits.size(), outer, inner);
  }

 private:
  std::shared_ptr<workflow::Conditions<Det, true>> mConditions;
  std::unique_ptr<const Geometry<Det>> mGeometry;
  std::unique_ptr<TrackingGeometry<Det>> mTracking;
  o2::iotof::GeometryTGeo* mTOF = nullptr;
  Parameters mParameters;
  float mSigmaNS = 0.020F;
};

class ClusterReader : public Task
{
 public:
  void init(InitContext& ic) final
  {
    const auto name = ic.options().get<std::string>("tof-cluster-infile");
    mFile.reset(TFile::Open(name.c_str()));
    if (!mFile || mFile->IsZombie()) {
      throw std::runtime_error("cannot open TOF clusters: " + name);
    }
    const auto* info = dynamic_cast<TStreamerInfo*>(mFile->GetStreamerInfoList()->FindObject("o2::iotof::Cluster"));
    if (!info || info->GetClassVersion() < 3) {
      throw std::runtime_error("TOF clusters lack BC/fine timing; re-cluster the digits with the timing bridge");
    }
    mTree = mFile->Get<TTree>("o2sim");
    if (!mTree || mTree->SetBranchAddress("TF3Cluster", &mPointer) < 0) {
      throw std::runtime_error("missing o2sim/TF3Cluster in " + name);
    }
  }
  void run(ProcessingContext& pc) final
  {
    if (mEntry < mTree->GetEntries()) {
      if (mTree->GetEntry(mEntry++) <= 0) {
        throw std::runtime_error("failed to read a TOF cluster entry");
      }
      pc.outputs().snapshot(Output{"TF3", "CLUSTERS", 0}, mClusters);
    }
    if (mEntry >= mTree->GetEntries()) {
      pc.services().get<ControlService>().endOfStream();
      pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
    }
  }

 private:
  std::unique_ptr<TFile> mFile;
  TTree* mTree = nullptr;
  std::vector<o2::iotof::Cluster> mClusters;
  std::vector<o2::iotof::Cluster>* mPointer = &mClusters;
  Long64_t mEntry = 0;
};
} // namespace

DataProcessorSpec getMatcherSpec()
{
  std::vector<InputSpec> inputs{{"tracks", workflow::gDataOrigin, workflow::description<Det>(workflow::kTracks), 0, Lifetime::Timeframe},
                                {"tofClusters", "TF3", "CLUSTERS", 0, Lifetime::Timeframe}};
  auto conditions = std::make_shared<workflow::Conditions<Det, true>>(inputs);
  return {"gtrk-tof-matcher", inputs, {{workflow::gDataOrigin, "TRKTOFMATCH", 0, Lifetime::Timeframe}}, AlgorithmSpec{adaptFromTask<MatcherTask>(conditions)}, Options{{"tof-max-spatial-chi2", VariantType::Double, 25., {"maximum spatial matching chi2"}}, {"tof-max-time-chi2", VariantType::Double, 25., {"maximum inner/outer time chi2"}}, {"tof-time-margin-ns", VariantType::Double, 0.1, {"coarse tracking interval margin in ns"}}, {"tof-time-resolution-ps", VariantType::Double, 20., {"single-cluster time sigma in ps"}}, {"tof-max-step-cm", VariantType::Double, 2., {"maximum flight integration step in cm"}}, {"tof-max-path-cm", VariantType::Double, 1500., {"maximum flight path in cm"}}}};
}

DataProcessorSpec getClusterReaderSpec()
{
  return {"gtrk-tof-cluster-reader", {}, {{"TF3", "CLUSTERS", 0, Lifetime::Timeframe}}, AlgorithmSpec{adaptFromTask<ClusterReader>()}, Options{{"tof-cluster-infile", VariantType::String, "tf3clusters.root", {"BC/fine-time TF3 cluster file"}}}};
}

DataProcessorSpec getMatchWriterSpec()
{
  using Writer = MakeRootTreeWriterSpec;
  return Writer("gtrk-tof-match-writer", "o2_gtrk_tof.root", Writer::TreeAttributes{"o2sim", "GTRK TOF associations and track times"},
                Writer::BranchDefinition<std::vector<Match>>{InputSpec{"tofMatches", workflow::gDataOrigin, "TRKTOFMATCH", 0}, "TRKTOFMatch"})();
}
} // namespace o2::gtrk::tof

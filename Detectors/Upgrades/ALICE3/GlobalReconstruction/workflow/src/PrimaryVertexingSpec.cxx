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

/// @file  PrimaryVertexingSpec.cxx

#include "ALICE3GlobalReconstructionWorkflow/PrimaryVertexingSpec.h"

#include "CommonConstants/LHCConstants.h"
#include "CommonDataFormat/BunchFilling.h"
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/NameConf.h"
#include "DataFormatsITS/TrackITS.h"
#include "DetectorsVertexing/PVertexer.h"
#include "DetectorsVertexing/PVertexerHelpers.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/DeviceSpec.h"
#include "Framework/Task.h"
#include "MathUtils/Utils.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCEventLabel.h"

#include <TStopwatch.h>

#include <cmath>
#include <vector>

using namespace o2::framework;

namespace o2::trk
{

namespace
{
using GTrackID = o2::dataformats::GlobalTrackID;
using PVertex = o2::dataformats::PrimaryVertex;
using V2TRef = o2::dataformats::VtxTrackRef;
using VtxTrackIndex = o2::dataformats::VtxTrackIndex;

class PrimaryVertexingSpec final : public Task
{
 public:
  PrimaryVertexingSpec(bool useMC, bool skip) : mUseMC(useMC), mSkip(skip) {}

  void init(InitContext& ic) final
  {
    mTimer.Stop();
    mTimer.Reset();
    mVertexer.setPoolDumpDirectory(ic.options().get<std::string>("pool-dumps-directory"));
    mVertexer.setTrackSources(GTrackID::getSourceMask(GTrackID::ITS));
    mVertexer.setITSROFrameLength(ic.options().get<float>("alice3-pv-rof-length-bc") * o2::constants::lhc::LHCBunchSpacingMUS);
    o2::BunchFilling bunchFilling;
    for (int bc = 0; bc < o2::constants::lhc::LHCMaxBunches; ++bc) {
      bunchFilling.setBC(bc);
    }
    mVertexer.setBunchFilling(bunchFilling);
    mVertexer.init();
    mVertexer.setBz(ic.options().get<float>("alice3-pv-bz"));
  }

  void run(ProcessingContext& pc) final
  {
    const double timeCPU0 = mTimer.CpuTime();
    const double timeReal0 = mTimer.RealTime();
    mTimer.Start(false);

    std::vector<PVertex> vertices;
    std::vector<VtxTrackIndex> vertexTrackIDs;
    std::vector<V2TRef> v2tRefs;
    std::vector<o2::MCEventLabel> vertexLabels;

    if (!mSkip) {
      const auto inputTracks = pc.inputs().get<gsl::span<o2::its::TrackITS>>("tracks");
      gsl::span<const o2::MCCompLabel> inputLabels;
      if (mUseMC) {
        inputLabels = pc.inputs().get<gsl::span<o2::MCCompLabel>>("labels");
      }

      std::vector<o2::vertexing::TrackWithTimeStamp> tracks;
      std::vector<GTrackID> gids;
      std::vector<o2::MCCompLabel> trackLabels;
      tracks.reserve(inputTracks.size());
      gids.reserve(inputTracks.size());
      if (mUseMC) {
        trackLabels.reserve(inputTracks.size());
      }

      const auto& params = o2::vertexing::PVertexerParams::Instance();
      for (size_t iTrack = 0; iTrack < inputTracks.size(); ++iTrack) {
        const auto& track = inputTracks[iTrack];
        if (track.getX() > params.trackMaxX) {
          continue;
        }
        if (o2::math_utils::numberOfBitsSet(track.getPattern() & 0x7) < params.minIBHits) {
          continue;
        }
        const auto& ts = track.getTimeStamp();
        const float timeMUS = static_cast<float>(ts.getTimeStamp()) * o2::constants::lhc::LHCBunchSpacingMUS;
        const float timeErrorMUS = static_cast<float>(ts.getTimeStampError()) * o2::constants::lhc::LHCBunchSpacingMUS;
        if (timeErrorMUS >= params.maxTimeErrorMUS) {
          continue;
        }

        tracks.emplace_back(o2::vertexing::TrackWithTimeStamp{track, {timeMUS, timeErrorMUS}});
        gids.emplace_back(static_cast<int>(iTrack), GTrackID::ITS);
        if (mUseMC) {
          trackLabels.emplace_back(iTrack < inputLabels.size() ? inputLabels[iTrack] : o2::MCCompLabel{});
        }
      }

      std::vector<o2::vertexing::InteractionCandidate> interactionCandidates;
      mVertexer.process(tracks, gids, interactionCandidates, vertices, vertexTrackIDs, v2tRefs, trackLabels, vertexLabels);
    }

    pc.outputs().snapshot(Output{"GLO", "PVTX", 0}, vertices);
    pc.outputs().snapshot(Output{"GLO", "PVTX_CONTID", 0}, vertexTrackIDs);
    pc.outputs().snapshot(Output{"GLO", "PVTX_CONTIDREFS", 0}, v2tRefs);
    if (mUseMC) {
      pc.outputs().snapshot(Output{"GLO", "PVTX_MCTR", 0}, vertexLabels);
    }

    mTimer.Stop();
    LOGP(info, "ALICE3 primary vertexing found {} PVs from TRK tracks, timing CPU/Real: {:.3f}/{:.3f} s",
         vertices.size(), mTimer.CpuTime() - timeCPU0, mTimer.RealTime() - timeReal0);

    static bool first = true;
    if (first) {
      first = false;
      if (pc.services().get<const o2::framework::DeviceSpec>().inputTimesliceId == 0) {
        o2::conf::ConfigurableParam::write(o2::base::NameConf::getConfigOutputFileName(pc.services().get<const o2::framework::DeviceSpec>().name, o2::vertexing::PVertexerParams::Instance().getName()), o2::vertexing::PVertexerParams::Instance().getName());
      }
    }
  }

  void endOfStream(EndOfStreamContext& ec) final
  {
    mVertexer.end();
    LOGF(info, "ALICE3 primary vertexing total timing: Cpu: %.3e Real: %.3e s in %d slots",
         mTimer.CpuTime(), mTimer.RealTime(), mTimer.Counter() - 1);
  }

 private:
  bool mUseMC{false};
  bool mSkip{false};
  o2::vertexing::PVertexer mVertexer;
  TStopwatch mTimer;
};

} // namespace

DataProcessorSpec getPrimaryVertexingSpec(bool useMC, bool skip)
{
  std::vector<InputSpec> inputs;
  inputs.emplace_back("tracks", "TRK", "TRACKS", 0, Lifetime::Timeframe);
  if (useMC) {
    inputs.emplace_back("labels", "TRK", "TRACKSMCTR", 0, Lifetime::Timeframe);
  }

  std::vector<OutputSpec> outputs;
  outputs.emplace_back("GLO", "PVTX", 0, Lifetime::Timeframe);
  outputs.emplace_back("GLO", "PVTX_CONTID", 0, Lifetime::Timeframe);
  outputs.emplace_back("GLO", "PVTX_CONTIDREFS", 0, Lifetime::Timeframe);
  if (useMC) {
    outputs.emplace_back("GLO", "PVTX_MCTR", 0, Lifetime::Timeframe);
  }

  return DataProcessorSpec{
    "alice3-primary-vertexing",
    inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<PrimaryVertexingSpec>(useMC, skip)},
    Options{{"pool-dumps-directory", VariantType::String, "", {"Destination directory for the tracks pool dumps"}},
            {"alice3-pv-bz", VariantType::Float, 5.f, {"Nominal ALICE3 Bz field in kG for primary vertex refits"}},
            {"alice3-pv-rof-length-bc", VariantType::Float, 1.f, {"Effective TRK readout length in BC used by pvertexer time clustering"}}}};
}

} // namespace o2::trk

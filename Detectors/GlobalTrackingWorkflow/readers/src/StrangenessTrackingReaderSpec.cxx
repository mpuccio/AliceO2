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

#include "GlobalTrackingWorkflowReaders/StrangenessTrackingReaderSpec.h"
#include <vector>
#include "Framework/ControlService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/Logger.h"
#include "CommonUtils/NameConf.h"
#include "CommonDataFormat/TimeStamp.h"
#include "ReconstructionDataFormats/StrangeTrack.h"
#include "ReconstructionDataFormats/StrangenessTrackingStats.h"
#include "DataFormatsITS/TrackITS.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "TFile.h"
#include "TTree.h"

using namespace o2::framework;
namespace o2::strangeness_tracking
{
class StrangenessTrackingReader : public Task
{
 public:
  StrangenessTrackingReader(bool mc, bool extensions, bool shortTracks) : mUseMC(mc), mExtensions(extensions), mShortTracks(shortTracks) {}
  void init(InitContext& ic) final
  {
    mFileName = o2::utils::Str::concat_string(o2::utils::Str::rectifyDirectory(ic.options().get<std::string>("input-dir")), ic.options().get<std::string>("strange-tracks-infile"));
    mFile.reset(TFile::Open(mFileName.c_str()));
    if (!mFile || mFile->IsZombie()) {
      LOG(fatal) << "Cannot open strange-track input " << mFileName;
    }
    mTree.reset(static_cast<TTree*>(mFile->Get("o2sim")));
    if (!mTree) {
      LOG(fatal) << "Missing o2sim tree in " << mFileName;
    }
    // Legacy mode does not require any of the new branches.
    mTree->SetBranchAddress("StrangeTracks", &mTracksPtr);
    if (mUseMC) {
      mTree->SetBranchAddress("StrangeTrackMCLab", &mLabelsPtr);
    }
    auto bind = [&](const char* name, auto address) {
      if (!mTree->GetBranch(name) || mTree->SetBranchAddress(name, address) < 0) {
        LOG(fatal) << "Missing or incompatible enabled strangeness-tracking branch " << name << " in " << mFileName;
      }
    };
    if (mExtensions) {
      bind("ExtendedStrangeTracks", &mExtendedPtr);
      bind("ExtendedStrangeTrackTimes", &mTimesPtr);
      if (mUseMC) {
        bind("ExtendedStrangeTrackMCLab", &mExtendedLabelsPtr);
      }
    }
    if (mShortTracks) {
      bind("ITSShortTracks", &mShortPtr);
      bind("ITSShortTrackPVRefs", &mShortPVPtr);
      bind("ITSShortTrackClusRefs", &mShortClusPtr);
      if (mUseMC) {
        bind("ITSShortTrackMCLab", &mShortLabelsPtr);
      }
    }
    if (mExtensions || mShortTracks) {
      bind("StrangenessTrackingStats", &mStatsPtr);
    }
  }
  void run(ProcessingContext& pc) final
  {
    const auto entry = mTree->GetReadEntry() + 1;
    const bool noEntry = entry >= mTree->GetEntries();
    if (!noEntry) {
      const auto bytes = mTree->GetEntry(entry);
      if ((mExtensions || mShortTracks) && bytes <= 0) {
        LOG(fatal) << "Cannot read strange-track entry " << entry << " in " << mFileName;
      }
    }
    if (mExtensions && (mTimes.size() != mExtended.size() || (mUseMC && mExtendedLabels.size() != mExtended.size()))) {
      LOG(fatal) << "Extended strange-track/time/MC cardinality mismatch at entry " << entry;
    }
    if (mShortTracks) {
      if (mShortPV.size() != mShort.size() || mShortClus.size() != 3 * mShort.size() || (mUseMC && mShortLabels.size() != mShort.size())) {
        LOG(fatal) << "ITS short-track/PV/cluster/MC cardinality mismatch at entry " << entry;
      }
      for (size_t i = 0; i < mShort.size(); ++i) {
        if (mShort[i].getNClusters() != 3 || mShort[i].getFirstClusterEntry() != 3 * i) {
          LOG(fatal) << "Invalid ITS short-track cluster reference range at entry " << entry;
        }
      }
    }
    pc.outputs().snapshot(Output{"GLO", "STRANGETRACKS", 0}, mTracks);
    if (mUseMC) {
      pc.outputs().snapshot(Output{"GLO", "STRANGETRACKS_MC", 0}, mLabels);
    }
    if (mExtensions) {
      pc.outputs().snapshot(Output{"GLO", "EXTSTRTRACKS", 0}, mExtended);
      pc.outputs().snapshot(Output{"GLO", "EXTSTRTIMES", 0}, mTimes);
      if (mUseMC) {
        pc.outputs().snapshot(Output{"GLO", "EXTSTRTRACKS_MC", 0}, mExtendedLabels);
      }
    }
    if (mShortTracks) {
      pc.outputs().snapshot(Output{"GLO", "ITSSHORTTRACKS", 0}, mShort);
      pc.outputs().snapshot(Output{"GLO", "ITSSHORTPV", 0}, mShortPV);
      pc.outputs().snapshot(Output{"GLO", "ITSSHORTCLREF", 0}, mShortClus);
      if (mUseMC) {
        pc.outputs().snapshot(Output{"GLO", "ITSSHORTMC", 0}, mShortLabels);
      }
    }
    if (mExtensions || mShortTracks) {
      pc.outputs().snapshot(Output{"GLO", "STRTRACKSTATS", 0}, mStats);
    }
    if (noEntry || mTree->GetReadEntry() + 1 >= mTree->GetEntries()) {
      pc.services().get<ControlService>().endOfStream();
      pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
    }
  }

 private:
  bool mUseMC = false, mExtensions = false, mShortTracks = false;
  std::vector<o2::dataformats::StrangeTrack> mTracks, *mTracksPtr = &mTracks, mExtended, *mExtendedPtr = &mExtended;
  std::vector<o2::MCCompLabel> mLabels, *mLabelsPtr = &mLabels, mExtendedLabels, *mExtendedLabelsPtr = &mExtendedLabels, mShortLabels, *mShortLabelsPtr = &mShortLabels;
  std::vector<o2::dataformats::TimeStampWithError<float, float>> mTimes, *mTimesPtr = &mTimes;
  std::vector<o2::its::TrackITS> mShort, *mShortPtr = &mShort;
  std::vector<int> mShortPV, *mShortPVPtr = &mShortPV, mShortClus, *mShortClusPtr = &mShortClus;
  o2::dataformats::StrangenessTrackingStats mStats, *mStatsPtr = &mStats;
  std::unique_ptr<TFile> mFile;
  std::unique_ptr<TTree> mTree;
  std::string mFileName;
};

DataProcessorSpec getStrangenessTrackingReaderSpec(bool useMC, bool enableExtensions, bool enableShortTracks)
{
  std::vector<OutputSpec> outputs{{"GLO", "STRANGETRACKS", 0, Lifetime::Timeframe}};
  if (useMC) {
    outputs.emplace_back("GLO", "STRANGETRACKS_MC", 0, Lifetime::Timeframe);
  }
  if (enableExtensions) {
    outputs.emplace_back("GLO", "EXTSTRTRACKS", 0, Lifetime::Timeframe);
    outputs.emplace_back("GLO", "EXTSTRTIMES", 0, Lifetime::Timeframe);
    if (useMC) {
      outputs.emplace_back("GLO", "EXTSTRTRACKS_MC", 0, Lifetime::Timeframe);
    }
  }
  if (enableShortTracks) {
    outputs.emplace_back("GLO", "ITSSHORTTRACKS", 0, Lifetime::Timeframe);
    outputs.emplace_back("GLO", "ITSSHORTPV", 0, Lifetime::Timeframe);
    outputs.emplace_back("GLO", "ITSSHORTCLREF", 0, Lifetime::Timeframe);
    if (useMC) {
      outputs.emplace_back("GLO", "ITSSHORTMC", 0, Lifetime::Timeframe);
    }
  }
  if (enableExtensions || enableShortTracks) {
    outputs.emplace_back("GLO", "STRTRACKSTATS", 0, Lifetime::Timeframe);
  }
  return {"strangeness-tracking-reader", Inputs{}, outputs,
          AlgorithmSpec{adaptFromTask<StrangenessTrackingReader>(useMC, enableExtensions, enableShortTracks)},
          Options{{"strange-tracks-infile", VariantType::String, "o2_strange_tracks.root", {"Name of the input strange tracks file"}},
                  {"input-dir", VariantType::String, "none", {"Input directory"}}}};
}
} // namespace o2::strangeness_tracking

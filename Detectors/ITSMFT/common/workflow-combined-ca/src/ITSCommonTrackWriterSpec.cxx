// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTCombinedCAWorkflow/ITSCommonTrackWriterSpec.h"

#include <vector>

#include "DPLUtils/MakeRootTreeWriterSpec.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::framework;

namespace o2::itsmft::combined
{

template <typename T>
using BranchDefinition = MakeRootTreeWriterSpec::BranchDefinition<T>;
using LabelsType = std::vector<o2::MCCompLabel>;

DataProcessorSpec getITSCommonTrackWriterSpec(bool useMC)
{
  auto tracksSize = std::make_shared<size_t>(0);
  auto tracksSizeGetter = [tracksSize](std::vector<o2::its::TrackITS> const& tracks) {
    *tracksSize = tracks.size();
  };
  auto logger = [tracksSize](std::vector<o2::itsmft::ROFRecord> const& rofs) {
    LOG(info) << "ITSMFTCombinedCATrackWriter pulled " << *tracksSize << " ITS tracks, in " << rofs.size() << " RO frames";
  };
  return MakeRootTreeWriterSpec("itsmft-combined-ca-its-track-writer",
                                "o2trac_its_ca.root",
                                MakeRootTreeWriterSpec::TreeAttributes{"o2sim", "Tree with ITS common-CA tracks"},
                                BranchDefinition<std::vector<o2::its::TrackITS>>{InputSpec{"tracks", "ITS", "TRACKS", 0},
                                                                                 "ITSTrack",
                                                                                 tracksSizeGetter},
                                BranchDefinition<std::vector<int>>{InputSpec{"trackClIdx", "ITS", "TRACKCLSID", 0},
                                                                   "ITSTrackClusIdx"},
                                BranchDefinition<std::vector<o2::itsmft::ROFRecord>>{InputSpec{"ROframes", "ITS", "ITSTrackROF", 0},
                                                                                     "ITSTracksROF",
                                                                                     logger},
                                BranchDefinition<LabelsType>{InputSpec{"labels", "ITS", "TRACKSMCTR", 0},
                                                             "ITSTrackMCTruth",
                                                             (useMC ? 1 : 0),
                                                             ""})();
}

} // namespace o2::itsmft::combined

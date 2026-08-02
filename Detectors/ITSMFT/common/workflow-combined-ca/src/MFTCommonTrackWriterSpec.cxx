// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTCombinedCAWorkflow/MFTCommonTrackWriterSpec.h"

#include <vector>

#include "DPLUtils/MakeRootTreeWriterSpec.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::framework;

namespace o2::itsmft::combined
{

template <typename T>
using BranchDefinition = MakeRootTreeWriterSpec::BranchDefinition<T>;
using LabelsType = std::vector<o2::MCCompLabel>;

DataProcessorSpec getMFTCommonTrackWriterSpec(bool useMC)
{
  auto tracksSize = std::make_shared<int>(0);
  auto tracksSizeGetter = [tracksSize](std::vector<o2::mft::TrackMFT> const& tracks) {
    *tracksSize = tracks.size();
  };
  auto logger = [tracksSize](std::vector<o2::itsmft::ROFRecord> const& rofs) {
    LOG(info) << "ITSMFTCombinedCATrackWriter pulled " << *tracksSize << " MFT tracks, in " << rofs.size() << " RO frames";
  };
  return MakeRootTreeWriterSpec("itsmft-combined-ca-mft-track-writer",
                                "mfttracks.root",
                                MakeRootTreeWriterSpec::TreeAttributes{"o2sim", "Tree with MFT tracks"},
                                BranchDefinition<std::vector<o2::mft::TrackMFT>>{InputSpec{"tracks", "MFT", "TRACKS", 0},
                                                                                 "MFTTrack",
                                                                                 tracksSizeGetter},
                                BranchDefinition<std::vector<int>>{InputSpec{"trackClIdx", "MFT", "TRACKCLSID", 0},
                                                                   "MFTTrackClusIdx"},
                                BranchDefinition<std::vector<uint16_t>>{InputSpec{"trackSeedPat", "MFT", "TRACKSEEDPAT", 0},
                                                                        "MFTTrackSeedPattern",
                                                                        1,
                                                                        ""},
                                BranchDefinition<std::vector<o2::itsmft::ROFRecord>>{InputSpec{"ROframes", "MFT", "MFTTrackROF", 0},
                                                                                     "MFTTracksROF",
                                                                                     logger},
                                BranchDefinition<LabelsType>{InputSpec{"labels", "MFT", "TRACKSMCTR", 0},
                                                             "MFTTrackMCTruth",
                                                             (useMC ? 1 : 0),
                                                             ""})();
}

} // namespace o2::itsmft::combined

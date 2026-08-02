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

/// @file   MFTCATrackWriterSpec.h
///
/// Gate 4 C4: relocated unmodified from
/// Detectors/ITSMFT/MFT/workflow/include/MFTWorkflow/TrackWriterSpec.h (same
/// namespace/signature, same o2::mft::getTrackWriterSpec symbol) into this
/// neutral, common location so both the legacy o2-mft-reco-workflow / opt-in
/// o2-mft-ca-tracker-workflow (via O2::MFTWorkflow) and the opt-in combined
/// ITS+MFT workflow (via O2::ITSMFTCombinedCAWorkflow) share exactly one
/// compiled writer-spec implementation instead of each linking or
/// duplicating it separately.

#ifndef O2_ITSMFT_CAWRITER_MFTCATRACKWRITERSPEC_H_
#define O2_ITSMFT_CAWRITER_MFTCATRACKWRITERSPEC_H_

#include "Framework/DataProcessorSpec.h"

namespace o2::mft
{

/// create a processor spec
/// write MFT tracks a root file
o2::framework::DataProcessorSpec getTrackWriterSpec(bool useMC, bool useCA = false);

} // namespace o2::mft

#endif // O2_ITSMFT_CAWRITER_MFTCATRACKWRITERSPEC_H_

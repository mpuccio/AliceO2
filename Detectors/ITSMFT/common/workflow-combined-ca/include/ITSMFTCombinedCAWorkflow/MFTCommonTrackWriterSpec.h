// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file MFTCommonTrackWriterSpec.h
/// \brief Gate 4 C4: local MFT CommonTrack ROOT writer for the combined
///        workflow. Same DPLUtils::MakeRootTreeWriterSpec construction,
///        default output file ("mfttracks.root"), tree ("o2sim"), and
///        branch names as MFTWorkflow/TrackWriterSpec.cxx's own
///        getTrackWriterSpec(useMC, /*useCA=*/true) -- a deliberate,
///        independent reimplementation so this library never links
///        O2::MFTWorkflow (see this directory's CMakeLists.txt).

#ifndef ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_MFTCOMMONTRACKWRITERSPEC_H_
#define ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_MFTCOMMONTRACKWRITERSPEC_H_

#include "Framework/DataProcessorSpec.h"

namespace o2::itsmft::combined
{

o2::framework::DataProcessorSpec getMFTCommonTrackWriterSpec(bool useMC);

} // namespace o2::itsmft::combined

#endif // ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_MFTCOMMONTRACKWRITERSPEC_H_

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
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
/// \file ITSCATrackWriterSpec.h
/// \brief Vertex-free ITS common-CA track writer. Writes a distinct file
///        (o2trac_its_ca.root) with no vertex branches.
///
/// Gate 4 C4: relocated unmodified from
/// Detectors/ITSMFT/ITS/workflow-ca/include/ITSCAWorkflow/TrackWriterSpec.h
/// (same namespace/signature, same o2::its::ca::getTrackWriterSpec symbol)
/// into this neutral, common location so o2-its-ca-tracker-workflow (via
/// O2::ITSCAWorkflow) and the opt-in combined ITS+MFT workflow (via
/// O2::ITSMFTCombinedCAWorkflow) share exactly one compiled writer-spec
/// implementation instead of each linking or duplicating it separately.

#ifndef O2_ITSMFT_CAWRITER_ITSCATRACKWRITERSPEC_H_
#define O2_ITSMFT_CAWRITER_ITSCATRACKWRITERSPEC_H_

#include "Framework/DataProcessorSpec.h"

namespace o2::its::ca
{

/// create a processor spec: write ITS common-CA tracks to a dedicated root
/// file (o2trac_its_ca.root), distinct from the legacy o2trac_its.root. No
/// vertex branch is ever declared.
o2::framework::DataProcessorSpec getTrackWriterSpec(bool useMC);

} // namespace o2::its::ca

#endif // O2_ITSMFT_CAWRITER_ITSCATRACKWRITERSPEC_H_

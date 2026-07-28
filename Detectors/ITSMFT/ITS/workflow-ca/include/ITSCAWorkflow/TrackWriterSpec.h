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
/// \file TrackWriterSpec.h
/// \brief Vertex-free ITS common-CA track writer (Gate 3 workflow-onboarding
///        Slice 2), modeled on MFTWorkflow/TrackWriterSpec.h. Writes a
///        distinct file (o2trac_its_ca.root) with no vertex branches.

#ifndef O2_ITS_CA_WORKFLOW_TRACKWRITERSPEC_H_
#define O2_ITS_CA_WORKFLOW_TRACKWRITERSPEC_H_

#include "Framework/DataProcessorSpec.h"

namespace o2::its::ca
{

/// create a processor spec: write ITS common-CA tracks to a dedicated root
/// file (o2trac_its_ca.root), distinct from the legacy o2trac_its.root. No
/// vertex branch is ever declared.
o2::framework::DataProcessorSpec getTrackWriterSpec(bool useMC);

} // namespace o2::its::ca

#endif // O2_ITS_CA_WORKFLOW_TRACKWRITERSPEC_H_

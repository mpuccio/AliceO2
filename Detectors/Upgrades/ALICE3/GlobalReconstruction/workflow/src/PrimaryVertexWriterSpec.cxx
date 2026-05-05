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

/// @file   PrimaryVertexWriterSpec.cxx

#include "ALICE3GlobalReconstructionWorkflow/PrimaryVertexWriterSpec.h"

#include "CommonDataFormat/RangeReference.h"
#include "DetectorsVertexing/PVertexerHelpers.h"
#include "DPLUtils/MakeRootTreeWriterSpec.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "SimulationDataFormat/MCEventLabel.h"

#include <vector>

using namespace o2::framework;

namespace o2::trk
{

template <typename T>
using BranchDefinition = MakeRootTreeWriterSpec::BranchDefinition<T>;

DataProcessorSpec getPrimaryVertexWriterSpec(bool useMC)
{
  auto logger = [](std::vector<o2::dataformats::PrimaryVertex> const& vertices) {
    LOG(info) << "ALICE3PrimaryVertexWriter pulled " << vertices.size() << " vertices";
  };

  return MakeRootTreeWriterSpec("alice3-primary-vertex-writer",
                                "o2_primary_vertex.root",
                                MakeRootTreeWriterSpec::TreeAttributes{"o2sim", "Tree with ALICE3 primary vertices"},
                                BranchDefinition<std::vector<o2::dataformats::PrimaryVertex>>{InputSpec{"vertices", "GLO", "PVTX", 0}, "PrimaryVertex", logger},
                                BranchDefinition<std::vector<o2::dataformats::VtxTrackRef>>{InputSpec{"v2tref", "GLO", "PVTX_CONTIDREFS", 0}, "PV2TrackRefs"},
                                BranchDefinition<std::vector<o2::dataformats::VtxTrackIndex>>{InputSpec{"vttrackID", "GLO", "PVTX_CONTID", 0}, "PVTrackIndices"},
                                BranchDefinition<std::vector<o2::MCEventLabel>>{InputSpec{"labels", "GLO", "PVTX_MCTR", 0}, "PVMCTruth", (useMC ? 1 : 0), ""})();
}

} // namespace o2::trk

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

/// \file ClusterConverterSpec.h
/// \brief A detector's clusters, as O2's cluster readers publish them, into
///        GTRK/<DET>CLUSTERS (DataSpecs.h).

#ifndef GTRK_WORKFLOW_CLUSTERCONVERTERSPEC_H_
#define GTRK_WORKFLOW_CLUSTERCONVERTERSPEC_H_

#include "GenericTracking/Base/LayerSpec.h"

#include "Framework/DataProcessorSpec.h"

namespace o2::gtrk::workflow
{

/// Defined in the detector's own unit (ITS.cxx, MFT.cxx, TRK.cxx). The ITS's
/// and MFT's staggered input has one stream per layer, concatenated in layer order.
template <DetectorLayout Det>
o2::framework::DataProcessorSpec getClusterConverterSpec(bool useMC, bool staggered = false);

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_CLUSTERCONVERTERSPEC_H_

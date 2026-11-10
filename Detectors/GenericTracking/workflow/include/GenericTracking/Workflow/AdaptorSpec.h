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

/// \file AdaptorSpec.h
/// \brief The O2-format adaptor of the ITS or the MFT, gtrk-<det>-adaptor: the
///        tracking device's products (TrackerSpec.h) as O2's tracker of the
///        detector publishes its own, for O2's consumers unchanged.

#ifndef GTRK_WORKFLOW_ADAPTORSPEC_H_
#define GTRK_WORKFLOW_ADAPTORSPEC_H_

#include "Framework/DataProcessorSpec.h"
#include "GenericTracking/Base/LayerSpec.h"

namespace o2::gtrk::workflow
{

/// The adaptor of \a Det (the ITS or the MFT); with \a useMC the track labels
/// from the clusters' (and, for the ITS, the vertices' outputs, without truth).
/// With staggered input, subscribe to each layer and publish composed O2 cluster
/// references; GTRK measurement indices remain dense in layer-stream order.
template <DetectorLayout Det>
o2::framework::DataProcessorSpec getAdaptorSpec(bool useMC, bool staggered = false);

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_ADAPTORSPEC_H_

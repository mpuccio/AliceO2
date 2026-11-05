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

/// \file MCVertexReaderSpec.h
/// \brief Publishes one smeared MC vertex seed per collision in a timeframe.

#ifndef GTRK_STUDY_MCVERTEXREADERSPEC_H_
#define GTRK_STUDY_MCVERTEXREADERSPEC_H_

#include "GenericTracking/Base/LayerSpec.h"

#include "Framework/DataProcessorSpec.h"

namespace o2::gtrk::study
{

/// \a Det's MC vertices as the native seeds (DataSpecs.h's `kVertices`),
/// smeared by GTRKMCVertex's widths (MCVertexParam.h).
template <DetectorLayout Det>
o2::framework::DataProcessorSpec getMCVertexReaderSpec();

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_MCVERTEXREADERSPEC_H_

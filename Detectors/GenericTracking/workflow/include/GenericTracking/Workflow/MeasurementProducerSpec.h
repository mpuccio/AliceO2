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

/// \file MeasurementProducerSpec.h
/// \brief Clusters in; the grid-ordered measurement columns, their prefix
///        offsets and the grid's layout with the beam line
///        (`MeasurementGeometry`) out. One device per detector, its own
///        instantiation.

#ifndef GTRK_WORKFLOW_MEASUREMENTPRODUCERSPEC_H_
#define GTRK_WORKFLOW_MEASUREMENTPRODUCERSPEC_H_

#include "GenericTracking/Base/LayerSpec.h"

#include "Framework/DataProcessorSpec.h"

namespace o2::gtrk::workflow
{

template <DetectorLayout Det>
o2::framework::DataProcessorSpec getMeasurementProducerSpec();

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_MEASUREMENTPRODUCERSPEC_H_

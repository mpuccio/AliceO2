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

/// \file MFT.cxx
/// \brief The MFT's own: its chips placed, its cluster converter and its
///        measurement producer.

#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Workflow/ClusterConverterSpec.h"
#include "GenericTracking/Workflow/Conditions.h"
#include "GenericTracking/Workflow/MeasurementProducerSpec.h"

#include "ClusterConverterITSMFT.h"
#include "MeasurementProducer.h"

#include "MFTBase/GeometryTGeo.h"

#include <vector>

namespace o2::gtrk
{

template <>
std::vector<Sensor> alice::placeChips<detectors::MFT>()
{
  return alice::placeAlpides(*o2::mft::GeometryTGeo::Instance());
}

template <>
o2::framework::DataProcessorSpec workflow::getClusterConverterSpec<detectors::MFT>(bool useMC, bool staggered)
{
  return workflow::getITSMFTClusterConverterSpec<detectors::MFT>(useMC, staggered);
}

template o2::framework::DataProcessorSpec workflow::getMeasurementProducerSpec<detectors::MFT>();

} // namespace o2::gtrk

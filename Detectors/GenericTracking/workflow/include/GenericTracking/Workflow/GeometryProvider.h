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

/// \file GeometryProvider.h
/// \brief A detector's geometry for a workflow device or a tool, placed from
///        the detector's GeometryTGeo.
///
/// The layout is placed from O2's GeometryTGeo, aligned as it is, and
/// checked against the layout's sensor counts (GenericTracking/Alice/Placement.h); its
/// layer envelopes are measured from the placements and logged. A workflow
/// device has the GeometryTGeo from the CCDB (Conditions.h) and places
/// on its first timeframe; a tool loads an o2sim geometry file first. The
/// layers' measurement depths are the caller's: the TRK's topology
/// dictionaries' in a device (Conditions.h), the dump's in the bench.
/// The holder keeps the result for its lifetime.

#ifndef GTRK_WORKFLOW_GEOMETRYPROVIDER_H_
#define GTRK_WORKFLOW_GEOMETRYPROVIDER_H_

#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Base/LayerEnvelopes.h"

#include <fairlogger/Logger.h>

#include <memory>
#include <span>

namespace o2::gtrk::workflow
{

/// The geometry from the detector's GeometryTGeo as it is: adopted from the
/// CCDB, or built from a loaded geometry file; \a depthsCM, each layer's
/// measurement depth, or none for 0 on every layer.
template <DetectorLayout Det>
[[nodiscard]] std::unique_ptr<const Geometry<Det>> makeGeometry(std::span<const double> depthsCM = {})
{
  const auto sensors = alice::placeChecked<Det>();
  auto geometry = std::make_unique<const Geometry<Det>>(sensors, depthsCM);
  LOG(info) << Det::name << " placed from its GeometryTGeo, layer envelopes measured:\n"
            << describeEnvelopes(PlacedEnvelopes<Det>::all(), Detector<Det>::nBarrelLayers);
  return geometry;
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_GEOMETRYPROVIDER_H_

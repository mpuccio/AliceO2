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

/// \file Placement.h
/// \brief A detector's chips as our sensors, from O2's GeometryTGeo: the frame
///        every detector shares; each detector's `placeChips` is its own unit's.
///
/// Host-only, in the workflow library: the tracking library does not depend
/// on it.
///
/// **The frame.** O2's L2G matrix maps a chip's local frame to the global one:
/// x across the rows, y along the normal, z along the columns, with the
/// origin already moved from the physical sensor's centre into the effective
/// sensitive layer (`GeometryTGeo::extractMatrixSensor`). Ours has y along the
/// columns and z along the normal, so a sensor's axes are ALICE's (x, z, -y),
/// which is right-handed. On the ITS that normal points away from the beam;
/// on the MFT it points along +z on odd faces and -z on even ones.
///
/// **The origin.** In the plane, the active matrix is not centred on the chip:
/// the readout edge is wider than the opposite one, which moves the matrix
/// centre along local x by `(PassiveEdgeReadOut - PassiveEdgeTop) / 2`,
/// 0.58 mm. A sensor's origin is put there, so that the bounds of
/// `detectors::alpide::kSensor` are symmetric about it, and a cluster's loc0 is
/// ALICE's local x less `kActiveCentreX`.
///
/// **Chip IDs** are O2's. The ITS numbers its chips layer by layer; the MFT
/// does not, and Geometry<Det> for a placed layout does not require it.

#ifndef GTRK_ALICE_PLACEMENT_H_
#define GTRK_ALICE_PLACEMENT_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/Sensor.h"

#include "MathUtils/Cartesian.h"
#include "MathUtils/Utils.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace o2::gtrk::alice
{

/// Where the active matrix's centre sits in ALICE's local x [cm].
extern const double kActiveCentreX;

/// Our flat sensor on layer \a layer from its L2G matrix: axes ALICE's
/// (x, z, -y), origin moved by \a activeCentreX along ALICE's x onto the
/// active matrix's centre (0 where the sensor volume is the active matrix).
[[nodiscard]] Sensor flatSensor(const o2::math_utils::Transform3D& l2g, int layer, double activeCentreX);

/// Our curved sensor on layer \a layer from its L2G matrix: O2's frame whole,
/// the cylinder's, z along its axis and the sensor's arc centred on +x
/// (SensorSpec::unrolledToFrame).
[[nodiscard]] Sensor curvedSensor(const o2::math_utils::Transform3D& l2g, int layer);

/// Every chip of \a geometry as a flat sensor with the ALPIDE's active
/// matrix (the ITS's and MFT's), indexed by O2's chip ID, `layer`
/// `GeometryTGeo::getLayer`, `ring` 0.
template <typename GeometryTGeo>
[[nodiscard]] std::vector<Sensor> placeAlpides(GeometryTGeo& geometry)
{
  geometry.fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  std::vector<Sensor> sensors;
  sensors.reserve(static_cast<std::size_t>(geometry.getNumberOfChips()));
  for (int chip = 0; chip < geometry.getNumberOfChips(); ++chip) {
    sensors.push_back(flatSensor(geometry.getMatrixL2G(chip), geometry.getLayer(chip), kActiveCentreX));
  }
  return sensors;
}

/// Every chip of \a Det as our Sensor, indexed by our sensor ID, from the
/// detector's GeometryTGeo: the one adopted from the CCDB, or one built from
/// the geometry in gGeoManager. Defined in the detector's own unit (ITS.cxx,
/// MFT.cxx, TRK.cxx).
template <DetectorLayout Det>
[[nodiscard]] std::vector<Sensor> placeChips();

/// The placed geometry of \a Det from the geometry in gGeoManager, checked
/// against the layout's sensor counts, in total and per layer. Throws
/// std::runtime_error on a disagreement: that is a wrong file, or an O2
/// geometry whose chip count changed. Where the chips sit is not checked
/// against anything compiled; the Geometry built from these sensors measures
/// its own layer envelopes.
template <DetectorLayout Det>
[[nodiscard]] std::vector<Sensor> placeChecked()
{
  auto sensors = placeChips<Det>();
  if (sensors.size() != static_cast<std::size_t>(Det::nSensors)) {
    throw std::runtime_error(std::string{Det::name} + ": " + std::to_string(sensors.size()) + " chips in the geometry, " + std::to_string(Det::nSensors) +
                             " in the layout");
  }
  std::vector<int> counts(static_cast<std::size_t>(Detector<Det>::nLayers), 0);
  for (const auto& sensor : sensors) {
    if (sensor.layer >= counts.size()) {
      throw std::runtime_error(std::string{Det::name} + ": a chip on layer " + std::to_string(sensor.layer) + ", beyond the layout's " +
                               std::to_string(Detector<Det>::nLayers) + " layers");
    }
    ++counts[sensor.layer];
  }
  for (std::size_t layer = 0; layer < counts.size(); ++layer) {
    const int stated = Detector<Det>::sensorCount(static_cast<LayerID>(layer));
    if (counts[layer] != stated) {
      throw std::runtime_error(std::string{Det::name} + ": layer " + std::to_string(layer) + " has " + std::to_string(counts[layer]) +
                               " chips in the geometry, " + std::to_string(stated) + " in the layout");
    }
  }
  return sensors;
}

} // namespace o2::gtrk::alice

#endif // GTRK_ALICE_PLACEMENT_H_

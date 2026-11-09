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

/// \file Material.h
/// \brief The ALICE 3 TRK's material as placed in gGeoManager, for the
///        material LUT (o2-gtrk-build-matlut) and the imported layouts
///        (o2-gtrk-import-geometry): the chips' metal stacks, the VD's side
///        walls and the cylinders around the beam line.

#ifndef GTRK_TOOLS_MATERIAL_H_
#define GTRK_TOOLS_MATERIAL_H_

#include "GenericTracking/Base/LayerEnvelopes.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/Sensor.h"

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

class TGeoNavigator;

namespace o2::gtrk::tools
{

/// The threads of a fill from the geometry: \a threads, or where it is
/// negative NTHREADS_MATBUD as O2's own fill reads it
/// (`MatLayerCylSet::getNThreadsFromEnv`); at least one.
[[nodiscard]] int fillThreads(int threads);

/// \a work(i, navigator) for every i below \a count, shared out over
/// \a threads threads (at least one), each taking the next item with a TGeo
/// navigator of its own; on one thread in order, on gGeoManager's current
/// navigator. Each item writes its own result, so the split does not move
/// it.
void forEachOnNavigators(std::size_t count, int threads, const std::function<void(std::size_t, TGeoNavigator*)>& work);

/// The metal stack O2 builds on a chip's sensor, by face [cm].
struct ChipStack {
  double inner = 0.; ///< on the face towards the beam
  double outer = 0.; ///< on the face away from it
};

/// The TRK's material on each barrel sensor's silicon, per barrel layer: the
/// metal stack O2 builds on one face of a sensor
/// (`GeometryTGeo::getTRKMetalStackPattern`) -- the outer on the VD and most
/// staves, the inner on the OT's first layer -- as thick as the simulation
/// made it (a VD chip is its configured x/X0 of silicon), so read from the
/// geometry in gGeoManager: along each of \a sensors' normals from its
/// centre, both ways, through the stack's volumes, the layer's thickest.
/// Throws std::runtime_error without a geometry or where a sensor's centre is
/// not in a sensor volume.
[[nodiscard]] std::vector<ChipStack> trkStackThickness(std::span<const Sensor> sensors, std::span<const SensorSpec> specs, int nBarrelLayers);

/// \a envelopes, the silicon's, with each barrel layer's \a stack: the hull
/// of its sensors taken that much thicker on each face, measured as
/// `measureEnvelopes` measures (at kEnvelopeResolution); a layer's z extent
/// and position stay its silicon's, a disc layer whole. The material bands of
/// the LUT and of the passive scatterers end there, so that a chip's stack is
/// its layer's material.
[[nodiscard]] std::vector<LayerEnvelope> withStacks(std::span<const LayerEnvelope> envelopes, std::span<const Sensor> sensors, std::span<const SensorSpec> specs, std::span<const ChipStack> stack);

/// The VD petals' side walls: O2's `VD_SideWall`s, flat plates between the
/// petals, along r and z [cm].
struct SideWalls {
  double thickness = 0.; ///< across the plate, in r-phi; 0 without walls
  double rMin = 0.;      ///< the radii the plates span
  double rMax = 0.;
};

/// The TRK's VD side walls as placed in gGeoManager: the thinnest plate and
/// the radii all of them span, from each plate's box (its local y across,
/// its local x along r) and its global matrix. A plate thinner than a
/// material LUT's cell lands on every track through the cell. Throws
/// std::runtime_error without a geometry.
[[nodiscard]] SideWalls trkSideWalls();

/// The radii where solid material begins or ends on cylinders around the beam
/// line, as placed in gGeoManager, below \a rMax [cm]: the beam pipe, the
/// iris, the VD's walls and cold plates, the service and support shells.
/// Along radial lines from the beam line at every \a dz in |z| < \a zMax and
/// at 256 azimuths a golden angle apart, the radii where a line enters or
/// leaves solid material (a gas, below 0.01 g/cm^3 -- air is 0.0012 -- does not
/// count; touching solids are one), to 1 um. A radius is a face where more
/// than half of the azimuths at one z meet it: a cylinder, or an arc around
/// most of it, meets every line at its radius; a flat face's radius changes
/// with the azimuth and a radial wall has none. Sorted. A material LUT band
/// ending there keeps the cylinder's material at its radius. The z positions
/// are shared out over \a threads, each with its
/// own TGeo navigator. Throws std::runtime_error without a geometry.
[[nodiscard]] std::vector<double> beamLineCylinderFaces(double rMax, double zMax, double dz, int threads = 1);

} // namespace o2::gtrk::tools

#endif // GTRK_TOOLS_MATERIAL_H_

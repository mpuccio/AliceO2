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

/// \file MeasurementGridBuilder.h
/// \brief Writing the grid of MeasurementGrid.h: bucket assignment and the
///        counting sort that realises the row order.
///
/// This is the producer's half of the index and it is host-side by nature:
/// turning `(surface, loc0, loc1)` into a global position needs `Geometry<Det>`,
/// which is double and lives on the host. The tracker never does it -- it
/// receives rows already in grid order with that position shipped beside
/// them, which is the whole point of the producer's boundary.
///
/// It allocates nothing. The producer asks the grid geometry how many buckets
/// there are, allocates the prefix offsets directly in the outgoing message,
/// and hands them here as a span. A caller can either recompute each bucket in
/// the scatter pass or retain one bucket number per input row.

#ifndef GTRK_TRACKING_MEASUREMENTGRIDBUILDER_H_
#define GTRK_TRACKING_MEASUREMENTGRIDBUILDER_H_

#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <span>

namespace o2::gtrk
{

/// A measurement's global position as shipped -- Cartesian, with transverse
/// radius and azimuth cached -- and the layer and second coordinate the grid
/// bins it by.
struct GridPoint {
  LayerID layer = kInvalidLayerID;
  Scalar x = 0.F; ///< global position [cm]
  Scalar y = 0.F;
  Scalar z = 0.F;
  Scalar r = 0.F;      ///< transverse radius [cm]
  Scalar phi = 0.F;    ///< in [0, 2pi)
  Scalar second = 0.F; ///< z on a barrel layer, r on a disc [cm]
};

/// Where a cluster sits.
///
/// The sensor-local position goes through the sensor's own transform, so a
/// tilted stave and a staggered disc ring are handled by the placement rather
/// than by a nominal-radius approximation -- which matters here: on a tilted
/// stave a nominal-phi binning would put rows in the wrong cell near its edge.
///
/// The transform is double, as all geometry is; the position narrows to float
/// once, and radius and azimuth are derived from those floats, never from the
/// double values. They are shipped as they are computed here, so a row's
/// bucket and the exact cut applied to it read the same numbers.
///
/// \a loc0 and \a loc1 are the sensor's unrolled coordinates: on a curved
/// sensor (the ALICE 3 VD) the point is on the cylinder at the angle loc0 /
/// radius (`SensorSpec::unrolledToFrame`), on a flat one the local point;
/// either at the layer's measurement depth (`Geometry::measurementDepth`).
///
/// The beam line's \a beamX and \a beamY are subtracted in double before the
/// narrowing: x, y, r and phi are about the corrected beam line, as the pair
/// cuts expect them, while
/// z and the local position stay the detector's.
template <DetectorLayout Det>
[[nodiscard]] inline GridPoint gridPoint(const Geometry<Det>& geometry, SensorID chip, double loc0, double loc1, double beamX = 0., double beamY = 0.)
{
  double local[3] = {0., 0., 0.};
  geometry.specOf(chip).unrolledToFrame(loc0, loc1, geometry.measurementDepth(geometry.sensor(chip).layer), local);
  double global[3] = {0., 0., 0.};
  geometry.toGlobal(chip, local, global);

  const auto layer = geometry.sensor(chip).layer;
  const auto x = static_cast<Scalar>(global[0] - beamX);
  const auto y = static_cast<Scalar>(global[1] - beamY);
  const auto z = static_cast<Scalar>(global[2]);
  const auto r = o2::gpu::CAMath::Hypot(x, y);
  return {.layer = layer, .x = x, .y = y, .z = z, .r = r, .phi = azimuth(x, y), .second = secondCoordinate<Det>(layer, r, z)};
}

/// A measurement's variances along z, r and r-phi.
struct PointVariances {
  Scalar z = 0.F;
  Scalar r = 0.F;
  Scalar rPhi = 0.F;
};

/// The local covariance `(varLoc0, covLoc01, varLoc1)` projected through the
/// sensor's rotation onto the global z axis and the radial and azimuthal
/// directions at \a point, in double as all geometry is. The local axes are
/// the sensor's; on a curved sensor the first turns with the point, at the
/// angle \a loc0 / radius, the second is the cylinder's axis.
///
/// The directions are those of the corrected beam line: r and r-phi are
/// calibrated quantities, like the radius and azimuth they accompany.
template <DetectorLayout Det>
[[nodiscard]] inline PointVariances pointVariances(const Geometry<Det>& geometry, SensorID chip, double loc0, const GridPoint& point, Scalar varLoc0, Scalar covLoc01, Scalar varLoc1)
{
  const auto& rot = geometry.sensor(chip).rot;
  double u[3] = {rot[0], rot[3], rot[6]}; // local axis 0, in global coordinates
  double v[3] = {rot[1], rot[4], rot[7]}; // local axis 1
  if (const auto& spec = geometry.specOf(chip); spec.curved()) {
    // In the cylinder's frame: across (-sin a, cos a, 0), along (0, 0, 1).
    const double angle = loc0 / spec.radius;
    const double across[3] = {-std::sin(angle), std::cos(angle), 0.};
    for (int axis = 0; axis < 3; ++axis) {
      u[axis] = (rot[(3 * axis) + 0] * across[0]) + (rot[(3 * axis) + 1] * across[1]);
      v[axis] = rot[(3 * axis) + 2];
    }
  }
  const auto x = static_cast<double>(point.x);
  const auto y = static_cast<double>(point.y);
  const double r = std::hypot(x, y);
  const double radial[3] = {r > 0. ? x / r : 1., r > 0. ? y / r : 0., 0.};
  const double azimuthal[3] = {-radial[1], radial[0], 0.};
  const double axial[3] = {0., 0., 1.};
  const auto project = [&](const double (&n)[3]) {
    const double a = (u[0] * n[0]) + (u[1] * n[1]) + (u[2] * n[2]);
    const double b = (v[0] * n[0]) + (v[1] * n[1]) + (v[2] * n[2]);
    return static_cast<Scalar>((static_cast<double>(varLoc0) * a * a) + (2. * static_cast<double>(covLoc01) * a * b) + (static_cast<double>(varLoc1) * b * b));
  };
  return {.z = project(axial), .r = project(radial), .rPhi = project(azimuthal)};
}

/// The counting sort that puts rows in `(layer, block, phi, second, time)`
/// order, in place in the offset array itself.
///
/// The offsets hold the per-bucket counts, then the prefix sum turns them into
/// starts, then `place` consumes them so that each entry ends up holding its
/// bucket's *end*, and `finish` rotates the array one step right to turn ends
/// back into starts. That last step is what makes the array a prefix sum again,
/// and it is exact rather than approximately right: after it, `offsets[b]` is
/// the first row of bucket b and `offsets[bucketCount]` is the row count.
///
/// The box-start key needs no work of its own: a caller that counts and
/// places each bucket's rows in box-start order gets each bucket out in that
/// order, a counting sort being stable.
class MeasurementGridBuilder
{
 public:
  /// \a offsets must be `GridGeometry::bucketCount() + 1` long.
  explicit MeasurementGridBuilder(std::span<std::uint32_t> offsets) : mOffsets(offsets)
  {
    assert(!offsets.empty());
    for (auto& entry : mOffsets) {
      entry = 0;
    }
  }

  /// Pass one, once per measurement.
  void count(std::uint32_t bucket)
  {
    assert(bucket + 1 < mOffsets.size());
    ++mOffsets[bucket];
  }

  /// Between the passes: counts become starts.
  void prefixSum()
  {
    std::uint32_t running = 0;
    for (auto& entry : mOffsets) {
      const auto here = entry;
      entry = running;
      running += here;
    }
  }

  /// Pass two, once per measurement and in the input's own order: the row this
  /// measurement occupies in the reordered columns.
  [[nodiscard]] std::uint32_t place(std::uint32_t bucket)
  {
    assert(bucket + 1 < mOffsets.size());
    return mOffsets[bucket]++;
  }

  /// After the second pass: ends become starts again.
  void finish()
  {
    for (std::size_t bucket = mOffsets.size() - 1; bucket > 0; --bucket) {
      mOffsets[bucket] = mOffsets[bucket - 1];
    }
    mOffsets[0] = 0;
  }

 private:
  std::span<std::uint32_t> mOffsets;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_MEASUREMENTGRIDBUILDER_H_

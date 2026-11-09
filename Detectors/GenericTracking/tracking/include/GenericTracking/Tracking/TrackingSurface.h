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

/// \file TrackingSurface.h
/// \brief The sensor planes as the tracking sees them, and the conversion
///        between a bound position and the global frame.
///
/// This is the device side of the host-to-tracking boundary: the double
/// Geometry<Det> narrows here, once, into fixed-size float records;
/// Geometry<Det> is the host side.
///
/// A surface's three axes are the **columns** of Sensor::rot, which is stored
/// row-major and maps local to global.  Writing the rows instead is a
/// transpose that a same-surface round trip cannot see, because it would use
/// the same wrong matrix in both directions.  testTrackingSurface.cxx
/// therefore checks against Geometry<Det>::toGlobal, an independent
/// implementation, and not against this file's own inverse.

#ifndef GTRK_TRACKING_TRACKINGSURFACE_H_
#define GTRK_TRACKING_TRACKINGSURFACE_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Trig.h"

#include "GPUCommonMath.h"

#include <array>
#include <cstdint>
#include <type_traits>

namespace o2::gtrk
{

/// The active area of a sensor, a rectangle in its own local frame.
struct SurfaceBounds {
  Scalar halfX = 0.F; ///< half width along loc0
  Scalar halfY = 0.F; ///< half length along loc1

  /// Written so that a NaN coordinate is outside.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool isInside(Scalar loc0, Scalar loc1) const { return loc1 >= -halfY && loc1 <= halfY && loc0 >= -halfX && loc0 <= halfX; }
};

/// One sensor plane.
///
/// Held as a record rather than as columns of a struct-of-arrays because a
/// propagation reads a whole surface at once -- origin, both in-plane axes and
/// the normal -- so the fifteen floats travel together.
struct TrackingSurface {
  Scalar origin[3]{}; ///< sensor centre in the global frame [cm]
  Scalar axis0[3]{};  ///< global direction of local loc0
  Scalar axis1[3]{};  ///< global direction of local loc1
  Scalar normal[3]{}; ///< global plane normal, axis0 x axis1

  SurfaceBounds bounds{};

  /// The surface's own chip id.  A propagation is handed the surface and not a
  /// geometry to look it up in, so without this the resulting state could not
  /// say which plane it is bound to.
  SensorID id = kInvalidSensorID;

  LayerID layer = kInvalidLayerID;

  /// A curved sensor's own entry is its cylinder, not a plane: the origin a
  /// point of its axis, axis0 and axis1 the sensor's local x and y, the normal
  /// its axis, the bounds unrolled, and this the radius its unrolled x is
  /// measured at (SensorSpec::radius). 0 on a flat sensor, and on the planes
  /// tangentPlane() makes. [cm]
  Scalar radius = 0.F;
  /// On a plane tangentPlane() made: its angle in the sensor, which a state
  /// bound to it carries (BoundTrackState::surfacePhi). 0 otherwise. [rad]
  Scalar surfacePhi = 0.F;
  /// A curved sensor's measurement depth (Geometry::measurementDepth): its
  /// tangent planes touch at radius + depth, while its unrolled x stays
  /// measured at the radius. 0 on a flat sensor, whose origin already lies at
  /// its depth, and on the planes tangentPlane() makes. [cm]
  Scalar depth = 0.F;
};

/// The plane touching the curved sensor \a cylinder at its local angle
/// \a angle, at radius + depth: normal
/// radial, axis1 along the cylinder's axis, axis0 = axis1 x normal along
/// increasing angle. Its origin is placed so that the touching point is at
/// loc0 = radius * angle, loc1 = 0 at the axis point: a point of the sensor near there has its unrolled
/// coordinates as its plane coordinates, and a measurement keeps its own. The
/// bounds are the sensor's unrolled ones.
GPUhdi() GPUAlwaysInline() TrackingSurface tangentPlane(const TrackingSurface& cylinder, Scalar angle)
{
  Scalar sinAngle = 0.F;
  Scalar cosAngle = 0.F;
  sincosOf(angle, sinAngle, cosAngle);
  TrackingSurface plane = cylinder;
  const Scalar arc = cylinder.radius * angle;
  for (int i = 0; i < 3; ++i) {
    const Scalar radial = cosAngle * cylinder.axis0[i] + sinAngle * cylinder.axis1[i];
    const Scalar along = -sinAngle * cylinder.axis0[i] + cosAngle * cylinder.axis1[i];
    plane.origin[i] = cylinder.origin[i] + (cylinder.radius + cylinder.depth) * radial - arc * along;
    plane.axis0[i] = along;
    plane.axis1[i] = cylinder.normal[i];
    plane.normal[i] = radial;
  }
  plane.radius = 0.F;
  plane.depth = 0.F;
  plane.surfacePhi = angle;
  return plane;
}

/// global = origin + loc0 * axis0 + loc1 * axis1
GPUhdi() GPUAlwaysInline() void surfaceToGlobal(const TrackingSurface& s, Scalar loc0, Scalar loc1, Scalar g[3])
{
  for (int i = 0; i < 3; ++i) {
    g[i] = s.origin[i] + loc0 * s.axis0[i] + loc1 * s.axis1[i];
  }
}

/// The inverse, including the out-of-plane component.
///
/// l[2] is returned rather than discarded so that a caller can check that a
/// point it believes is on the surface actually is, instead of assuming it.
GPUhdi() GPUAlwaysInline() void surfaceToLocal(const TrackingSurface& s, const Scalar g[3], Scalar l[3])
{
  const Scalar d0 = g[0] - s.origin[0];
  const Scalar d1 = g[1] - s.origin[1];
  const Scalar d2 = g[2] - s.origin[2];
  l[0] = d0 * s.axis0[0] + d1 * s.axis0[1] + d2 * s.axis0[2];
  l[1] = d0 * s.axis1[0] + d1 * s.axis1[1] + d2 * s.axis1[2];
  l[2] = d0 * s.normal[0] + d1 * s.normal[1] + d2 * s.normal[2];
}

/// \a n dot \a d, the cosine of incidence on a surface.  Its magnitude is what
/// the grazing check thresholds on.
GPUhdi() GPUAlwaysInline() Scalar cosIncidence(const TrackingSurface& s, const Scalar d[3]) { return s.normal[0] * d[0] + s.normal[1] * d[1] + s.normal[2] * d[2]; }

/// Every sensor plane of one detector, as the tracking sees them.
///
/// Fixed extent, so the whole thing is one contiguous block sized at compile
/// time and uploads as a single buffer.  Built from Geometry<Det>, which is
/// double and stays on the host; this is where that narrows.
template <DetectorLayout Det>
class TrackingGeometry
{
 public:
  explicit TrackingGeometry(const Geometry<Det>& geo) { build(geo); }

  /// A sensor's own entry: its plane, or a curved sensor's cylinder
  /// (TrackingSurface::radius). Its layer is the sensor's either way.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() const TrackingSurface& surface(SensorID id) const { return mSurfaces[id]; }

  /// What planeOf returns: the table entry itself where every sensor is flat,
  /// a plane built on the spot where some are curved.
  using Plane = std::conditional_t<Detector<Det>::hasCurvedSensors, TrackingSurface, const TrackingSurface&>;

  /// The plane a state is bound to: its sensor's, or on a curved sensor the
  /// tangent plane at the state's angle.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Plane planeOf(const BoundTrackState& state) const { return plane(state.surface, state.surfacePhi); }

  /// The plane a measurement is on: its sensor's, or on a curved sensor the
  /// tangent plane where it sits, at the angle its unrolled loc0 gives.
  template <typename M>
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Plane planeOf(const M& measurement) const
  {
    if constexpr (Detector<Det>::hasCurvedSensors) {
      const TrackingSurface& s = mSurfaces[measurement.surface];
      return s.radius > 0.F ? tangentPlane(s, measurement.loc0 / s.radius) : s;
    } else {
      return mSurfaces[measurement.surface];
    }
  }

  /// Sensor \a id's plane, or a curved sensor's tangent plane at \a angle.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Plane plane(SensorID id, [[maybe_unused]] Scalar angle) const
  {
    if constexpr (Detector<Det>::hasCurvedSensors) {
      const TrackingSurface& s = mSurfaces[id];
      return s.radius > 0.F ? tangentPlane(s, angle) : s;
    } else {
      return mSurfaces[id];
    }
  }

  [[nodiscard]] const std::array<TrackingSurface, Det::nSensors>& surfaces() const { return mSurfaces; }

 private:
  void build(const Geometry<Det>& geo);

  std::array<TrackingSurface, Det::nSensors> mSurfaces{};
};

template <DetectorLayout Det>
void TrackingGeometry<Det>::build(const Geometry<Det>& geo)
{
  for (SensorID id = 0; id < static_cast<SensorID>(Det::nSensors); ++id) {
    const Sensor& src = geo.sensor(id);
    const SensorSpec& spec = geo.specOf(id);
    TrackingSurface& dst = mSurfaces[id];
    // A flat sensor's plane at its layer's measurement depth, along its normal.
    const double depth = geo.measurementDepth(src.layer);
    const double shift = spec.curved() ? 0. : depth;

    for (int i = 0; i < 3; ++i) {
      dst.origin[i] = static_cast<Scalar>(src.trans[i] + (shift * src.rot[3 * i + 2]));
      // Columns, not rows: rot is row-major local -> global, so column j is
      // the global direction of local axis j.
      dst.axis0[i] = static_cast<Scalar>(src.rot[3 * i + 0]);
      dst.axis1[i] = static_cast<Scalar>(src.rot[3 * i + 1]);
      dst.normal[i] = static_cast<Scalar>(src.rot[3 * i + 2]);
    }

    dst.bounds.halfX = static_cast<Scalar>(spec.halfX);
    dst.bounds.halfY = static_cast<Scalar>(spec.halfY);
    dst.radius = spec.curved() ? static_cast<Scalar>(spec.radius) : 0.F;
    dst.depth = spec.curved() ? static_cast<Scalar>(depth) : 0.F;

    dst.id = id;
    dst.layer = src.layer;
  }
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRACKINGSURFACE_H_

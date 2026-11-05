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

/// \file ShapeCut.h
/// \brief A candidate's line against its measurements' cluster shapes.
///
/// A particle crossing a sensor leaves a cluster as long along the columns as
/// the sensor's depth times its local slope there: the slope of its line
/// along the columns over its length along the sensor's normal. The
/// measurement producer turns each cluster's spans into the interval of local
/// slopes a particle above its detector's calibration floor leaves such a
/// cluster at (`MeasurementT::slopeLow`, `slopeHigh`; ClusterShapeBand.h):
/// empty where none does (low infinite, high 0), [0, `kSlopeUnbounded`)
/// without a shape to go by. An iteration
/// with `shapeTest` holds every pair and triplet to it at each of its
/// measurements.
///
/// The slope is the barrel's: the normal radial at the measurement, the
/// columns along z. A calibrated layer is a barrel layer (the producer refuses
/// others); every other measurement's interval is [0, `kSlopeUnbounded`): no
/// upper limit short of the largest float.

#ifndef GTRK_TRACKING_SHAPECUT_H_
#define GTRK_TRACKING_SHAPECUT_H_

#include "GenericTracking/Tracking/Definitions.h"

#include "GPUCommonDef.h"
#include "GPUCommonMath.h"

namespace o2::gtrk
{

/// The upper end of a measurement's slope interval without a shape: the
/// largest float, the sentinel for an interval without an upper limit.
inline constexpr Scalar kSlopeUnbounded = FLT_MAX;

/// Whether the line along (\a dx, \a dy, \a dz) crosses the measurement at
/// (\a x, \a y), transverse radius \a r, at a local slope in [\a low,
/// \a high): |dz| over the line's length along the radial normal there,
/// compared without the division (a line along the sensor, of no length
/// along the normal, passes none).
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool slopeAdmits(Scalar low, Scalar high, Scalar x, Scalar y, Scalar r, Scalar dx, Scalar dy, Scalar dz)
{
  const Scalar across = o2::gpu::CAMath::Abs((dx * x) + (dy * y)); // the length along the normal, times r
  const Scalar along = o2::gpu::CAMath::Abs(dz) * r;
  return along >= low * across && along < high * across;
}

/// \a slopeAdmits for measurement \a row of \a columns, the line from
/// (\a fromX, \a fromY, \a fromZ) to (\a toX, \a toY, \a toZ).
template <typename Columns, typename Row>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool shapeAdmits(const Columns& columns, Row row, Scalar fromX, Scalar fromY, Scalar fromZ, Scalar toX, Scalar toY,
                                                          Scalar toZ)
{
  return slopeAdmits(columns.slopeLow[row], columns.slopeHigh[row], columns.x[row], columns.y[row], columns.r[row], toX - fromX, toY - fromY, toZ - fromZ);
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_SHAPECUT_H_

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

/// \file Measurement.h
/// \brief What a cluster looks like once the fit can consume it, as one
///        record and as the SoA columns the tracker reads.
///
/// The record is here; filling one is the workflow's (MeasurementProducer.h):
/// the cluster record carries ClassDefNV and Rtypes.h, which nothing here may
/// include.
///
/// Record and columns are one declaration, `MeasurementT`. A producer
/// allocates one message part per column and one for the grid offsets; a
/// consumer binds a pointer per column and reads `view[i].timePS`. Nothing
/// here serialises or knows about a transport; only the owned columns of the
/// tests and the replay format allocate.
#ifndef GTRK_TRACKING_MEASUREMENT_H_
#define GTRK_TRACKING_MEASUREMENT_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/TimeEstPS.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/ShapeCut.h"

#include <memlayout/wrapper.h>

#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <vector>

namespace o2::gtrk
{

/// One cluster, in the frame the update works in, declared once for every
/// storage: `Measurement` below is the record, `MeasurementStore` and
/// `MeasurementView` are its columns. The defaults apply to the record only;
/// in any other storage a field defaults to an empty column.
///
/// The spatial part is the sensor-local position and its 2x2 covariance. The
/// cross term is part of the update's contract, as in O2's cluster update;
/// the producer writes 0.
///
/// Instantiated with C = const_reference, every member *is* a reference --
/// that is what an SoA element proxy is -- and `apply` is upstream's macro,
/// which takes a forwarding reference it does not forward.  Neither is a
/// defect here.  The macro also uses these camelBack member names as enum
/// constants.
// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
template <template <class> class C>
struct MeasurementT {
  MEMLAYOUT_APPLY(MeasurementT, loc0, loc1, timePS, widthPS, varLoc0, covLoc01, varLoc1, clusterIndex, surface, x, y, z, r, phi, varZ, varR, varRPhi, slopeLow,
                  slopeHigh)
  C<Scalar> loc0 = memlayout::default_value<C>(0.F); ///< sensor-local position [cm]
  C<Scalar> loc1 = memlayout::default_value<C>(0.F);
  /// The time box `[timePS, timePS + widthPS)` [ps from the timeframe's
  /// start]: the readout frame or tick, widened by its margin
  /// (TimeEstPS.h).
  C<std::uint32_t> timePS = memlayout::default_value<C>(std::uint32_t{0});
  C<std::uint32_t> widthPS = memlayout::default_value<C>(std::uint32_t{0});

  C<Scalar> varLoc0 = memlayout::default_value<C>(0.F);
  C<Scalar> covLoc01 = memlayout::default_value<C>(0.F);
  C<Scalar> varLoc1 = memlayout::default_value<C>(0.F);

  C<std::uint32_t> clusterIndex = memlayout::default_value<C>(std::uint32_t{0});
  C<SensorID> surface = memlayout::default_value<C>(kInvalidSensorID);

  /// Global position [cm], computed once by the producer from the sensor and
  /// the double geometry. Everything that compares positions across layers --
  /// the exact spatial cuts, tracklets, triplets -- reads these; the fit works
  /// in the local frame above.
  C<Scalar> x = memlayout::default_value<C>(0.F);
  C<Scalar> y = memlayout::default_value<C>(0.F);
  C<Scalar> z = memlayout::default_value<C>(0.F);
  /// The same position's transverse radius [cm] and azimuth in [0, 2pi),
  /// cached by the producer so that no consumer takes a square root or an
  /// arctangent per candidate.
  C<Scalar> r = memlayout::default_value<C>(0.F);
  C<Scalar> phi = memlayout::default_value<C>(0.F);
  /// The local covariance projected through the sensor's rotation onto z, r
  /// and r-phi of the corrected beam line [cm^2]: what the pair cuts read,
  /// computed once by the producer. The r-z
  /// covariance is zero for sensors whose plane contains z or lies across it.
  C<Scalar> varZ = memlayout::default_value<C>(0.F);
  C<Scalar> varR = memlayout::default_value<C>(0.F);
  C<Scalar> varRPhi = memlayout::default_value<C>(0.F);
  /// The local slopes [slopeLow, slopeHigh) at which a particle above its
  /// detector's calibration floor leaves this cluster's shape: its line's
  /// slope along the sensor's columns over its length along the normal
  /// (ShapeCut.h), made by the producer from its detector's calibrated band
  /// (ClusterShapeBand.h);
  /// empty (infinite, 0) where none does, [0, kSlopeUnbounded) without a
  /// shape to go by. An iteration with `shapeTest` reads them.
  C<Scalar> slopeLow = memlayout::default_value<C>(0.F);
  C<Scalar> slopeHigh = memlayout::default_value<C>(kSlopeUnbounded);
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)

/// The record: one measurement, as the fit and the QA hold it.
using Measurement = MeasurementT<memlayout::value>;

/// The columns, owned. Useful on the host; not what goes on the wire.
using MeasurementStore = memlayout::Wrapper<MeasurementT, std::vector, memlayout::Layout::soa>;

/// The columns, borrowed and read-only. `const_pointer` is the contract: a
/// consumer cannot write through this view, whatever the message is backed by.
using MeasurementView = memlayout::Wrapper<MeasurementT, memlayout::const_pointer, memlayout::Layout::soa>;

/// One message part per field in the MEMLAYOUT_APPLY order. The wrapper
/// generates these IDs from that same field list.
using MeasurementColumn = Measurement::Column;

inline constexpr std::uint32_t kNMeasurementColumns = static_cast<std::uint32_t>(MeasurementColumn::ColumnCount);
/// Prefix offsets for every `(layer, time block, phi cell, second cell)`
/// bucket, followed by the total measurement count.
inline constexpr std::uint32_t kColGridOffsets = kNMeasurementColumns;
inline constexpr std::uint32_t kNMeasurementParts = kNMeasurementColumns + 1;

/// Build one container per schema field. The factory receives the field type
/// and its wire-column ID, so producers and consumers cannot silently omit a
/// field when MeasurementT grows.
template <template <class> class Container, class Factory>
[[nodiscard]] auto makeMeasurementColumns(Factory factory)
{
  std::size_t column = 0U;
  return Measurement{}.apply([&]<class... T>(const T&...) {
    return memlayout::Wrapper<MeasurementT, Container, memlayout::Layout::soa>{
      MeasurementT<Container>{factory(std::type_identity<T>{}, static_cast<MeasurementColumn>(column++))...}};
  });
}

/// Resize every owned column together.
template <class Columns>
void resizeMeasurementColumns(Columns& columns, std::size_t rows)
{
  columns.apply([rows](auto&... column) { (column.resize(rows), ...); });
}

/// Bind any contiguous measurement columns as the tracker's read-only view.
template <class Columns>
[[nodiscard]] MeasurementView measurementView(const Columns& columns)
{
  return MeasurementView{columns.apply([](const auto&... column) { return MeasurementT<memlayout::const_pointer>{column.data()...}; })};
}

template <class Columns>
void writeMeasurement(Columns& columns, std::size_t row, const Measurement& measurement)
{
  auto target = columns[row];
  target.apply(measurement, [](auto& destination, const auto& source) -> auto& {
    destination = source;
    return destination;
  });
}

/// Bytes occupied by one row across all measurement columns, excluding grid
/// offsets and transport metadata.
inline constexpr std::size_t kMeasurementRowBytes = Measurement{}.apply([](const auto&... field) { return (sizeof(field) + ... + 0U); });

/// A row's time box.
template <typename Row>
[[nodiscard]] GPUhdi() GPUAlwaysInline() TimeEstPS timeBox(const Row& row)
{
  return {.start = row.timePS, .width = row.widthPS};
}

/// The same for row \a row of \a columns, reading only the two columns it holds.
template <typename Columns>
[[nodiscard]] GPUhdi() GPUAlwaysInline() TimeEstPS timeBoxAt(const Columns& columns, std::uint32_t row)
{
  return {.start = columns.timePS[row], .width = columns.widthPS[row]};
}

/// Whether a measurement with time box \a row can come from a collision in
/// \a vertex's box: the rule of ITSMFT's `isVertexCompatible`, the vertex's end inclusive so
/// that a vertex of width 0 -- a point -- at the box's start belongs to it.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool vertexTimeCompatible(const TimeEstPS& row, const TimeEstPS& vertex)
{
  return vertex.upper() >= row.lower() && vertex.lower() < row.upper();
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_MEASUREMENT_H_

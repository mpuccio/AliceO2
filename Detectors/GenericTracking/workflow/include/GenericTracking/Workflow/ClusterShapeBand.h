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

/// \file ClusterShapeBand.h
/// \brief A detector's calibrated cluster-shape band, a CCDB object
///        (`<DET>/Calib/ClusterShapeBand`): per cluster shape the slope
///        interval a measurement carries.
///
/// Made from the simulation: per
/// calibrated layer, the most rows a particle above the calibration's floor
/// leaves, and per column span the hull of the local slopes (ShapeCut.h) at
/// which such a particle leaves that many columns -- the calibration's band
/// inverted, one interval per shape. A cluster
/// without a pattern, or on a layer the band does not calibrate, gets every
/// slope; one of more rows than its layer's most, or of a span no slope
/// gives, none. The slope is the barrel's (the normal radial, the columns
/// along z), so a band calibrates barrel layers only.

#ifndef GTRK_WORKFLOW_CLUSTERSHAPEBAND_H_
#define GTRK_WORKFLOW_CLUSTERSHAPEBAND_H_

#include "GenericTracking/Tracking/ShapeCut.h"

#include <Rtypes.h>

#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace o2::gtrk::workflow
{

struct ClusterShapeBand {
  /// A measurement's slope interval [first, second).
  using Interval = std::pair<Scalar, Scalar>;
  static constexpr Interval kEvery{0.F, kSlopeUnbounded};
  static constexpr Interval kNone{std::numeric_limits<Scalar>::infinity(), 0.F};
  /// The column spans a cluster can have, 0 to 255.
  static constexpr std::size_t kSpans = 256;

  std::vector<int> maxRows;     ///< per layer from the first, its most rows; 0 on a layer the band does not calibrate, as every one after
  std::vector<float> slopeLow;  ///< per layer and column span, at layer * kSpans + span: the interval's start
  std::vector<float> slopeHigh; ///< and its end

  /// Whether it is a band of a layout of \a barrelLayers barrel layers: of
  /// its first layers at most, every layer's spans there.
  [[nodiscard]] bool fits(std::size_t barrelLayers) const { return maxRows.size() <= barrelLayers && slopeLow.size() == maxRows.size() * kSpans && slopeHigh.size() == maxRows.size() * kSpans; }

  /// The slope interval of a cluster of \a rowSpan rows and \a colSpan
  /// columns on layer \a layer.
  [[nodiscard]] Interval interval(std::size_t layer, std::uint8_t rowSpan, std::uint8_t colSpan) const
  {
    if (colSpan == 0U || layer >= maxRows.size() || maxRows[layer] == 0) {
      return kEvery;
    }
    if (rowSpan > maxRows[layer]) {
      return kNone;
    }
    const auto at = (layer * kSpans) + colSpan;
    return {slopeLow[at], slopeHigh[at]};
  }

  /// One line: the layers it calibrates, with their most rows.
  [[nodiscard]] std::string describe() const
  {
    std::string layers;
    for (std::size_t layer = 0; layer < maxRows.size(); ++layer) {
      if (maxRows[layer] != 0) {
        layers += std::format("{}layer {} ({} rows)", layers.empty() ? "" : ", ", layer, maxRows[layer]);
      }
    }
    return "cluster-shape band: " + (layers.empty() ? std::string{"no layer"} : layers);
  }

  ClassDefNV(ClusterShapeBand, 1);
};

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_CLUSTERSHAPEBAND_H_

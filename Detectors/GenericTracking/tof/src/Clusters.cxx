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

#include "GenericTracking/TOF/Clusters.h"
#include "DataFormatsIOTOF/Cluster.h"
#include "IOTOFBase/GeometryTGeo.h"
#include "IOTOFBase/IOTOFBaseParam.h"
#include "CommonConstants/LHCConstants.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace o2::gtrk::tof
{
std::vector<Hit> convertClusters(std::span<const o2::iotof::Cluster> clusters, const o2::iotof::GeometryTGeo& geometry, float sigmaTimeNS)
{
  if (!(sigmaTimeNS > 0.F) || !std::isfinite(sigmaTimeNS) || clusters.size() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    throw std::invalid_argument("invalid TOF timing resolution or cluster count");
  }
  std::vector<Hit> hits;
  hits.reserve(clusters.size());
  for (std::size_t index = 0; index < clusters.size(); ++index) {
    const auto& cluster = clusters[index];
    const auto chip = cluster.getChipID();
    if (!geometry.isValidIOTOFChipIndex(chip)) {
      throw std::invalid_argument("TOF matching currently supports segmented barrel chips only");
    }
    if (!(cluster.getTime() >= 0.F && cluster.getTime() < o2::constants::lhc::LHCBunchSpacingNS)) {
      throw std::invalid_argument("TOF cluster fine time must be relative to its BC");
    }
    const auto layer = geometry.getIOTOFLayer(chip);
    const auto& specs = o2::iotof::ChipSpecificsParam::Instance();
    const auto mask = cluster.getPattern();
    if (mask == o2::iotof::Cluster::InvalidPatternID || mask == 0 || cluster.getRowSpan() == 0 || cluster.getColSpan() == 0 || cluster.getRowSpan() * cluster.getColSpan() > 16) {
      continue;
    }
    if (cluster.getRow() + cluster.getRowSpan() > static_cast<unsigned>(specs.NRows) || cluster.getCol() + cluster.getColSpan() > static_cast<unsigned>(specs.NCols)) {
      throw std::invalid_argument("TOF cluster lies outside its chip's pixel matrix");
    }
    double x = 0., z = 0.;
    int pixels = 0;
    for (int bit = 0; bit < cluster.getRowSpan() * cluster.getColSpan(); ++bit) {
      if ((mask & (1U << bit)) == 0) {
        continue;
      }
      const auto point = geometry.detectorToLocal(cluster.getRow() + bit / cluster.getColSpan(), cluster.getCol() + bit % cluster.getColSpan(), chip);
      x += point.X();
      z += point.Z();
      ++pixels;
    }
    if (pixels == 0) {
      continue;
    }
    const auto first = geometry.detectorToLocal(0, 0, chip);
    const auto last = geometry.detectorToLocal(specs.NRows - 1, specs.NCols - 1, chip);
    const float cx = (first.X() + last.X()) / 2.F;
    const float cz = (first.Z() + last.Z()) / 2.F;
    const auto& matrix = geometry.getMatrixL2G(chip);
    using Point = o2::math_utils::Point3D<double>;
    const auto origin = matrix * Point{cx, 0.F, cz};
    const auto axis0 = matrix * Point{cx + 1.F, 0.F, cz};
    const auto axis1 = matrix * Point{cx, 0.F, cz + 1.F};
    Hit hit{.time = {cluster.getBc(), cluster.getTime()},
            .loc0 = static_cast<float>(x / pixels - cx),
            .loc1 = static_cast<float>(z / pixels - cz),
            .varLoc0 = specs.PitchRow * specs.PitchRow / 12.F,
            .varLoc1 = specs.PitchCol * specs.PitchCol / 12.F,
            .sigmaTimeNS = sigmaTimeNS,
            .index = static_cast<std::int32_t>(index),
            .layer = layer == 0 ? Layer::Inner : Layer::Outer};
    auto& surface = hit.surface;
    surface.origin[0] = origin.X();
    surface.origin[1] = origin.Y();
    surface.origin[2] = origin.Z();
    surface.axis0[0] = axis0.X() - origin.X();
    surface.axis0[1] = axis0.Y() - origin.Y();
    surface.axis0[2] = axis0.Z() - origin.Z();
    surface.axis1[0] = axis1.X() - origin.X();
    surface.axis1[1] = axis1.Y() - origin.Y();
    surface.axis1[2] = axis1.Z() - origin.Z();
    for (int i = 0; i < 3; ++i) {
      const int j = (i + 1) % 3, k = (i + 2) % 3;
      surface.normal[i] = surface.axis0[j] * surface.axis1[k] - surface.axis0[k] * surface.axis1[j];
    }
    surface.id = chip;
    surface.bounds = {.halfX = specs.ActiveMatrixSizeRows() / 2.F, .halfY = specs.ActiveMatrixSizeCols() / 2.F};
    hits.push_back(hit);
  }
  return hits;
}
} // namespace o2::gtrk::tof

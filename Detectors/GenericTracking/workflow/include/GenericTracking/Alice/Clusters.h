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

/// \file Clusters.h
/// \brief O2's compact clusters as our clusters: the steps every format takes,
///        and the ITS's and MFT's format, `o2::itsmft::CompClusterExt`.
///
/// Two steps, so that the second can be tested without a dictionary:
///
/// - `localClusters` gives each compact cluster the local position, errors
///   and pixel count of O2's ITS/MFT decode, `itsmft::ioutils::extractClusterData`:
///   the dictionary's position and `Err2X`/`Err2Z` for a common topology, the
///   pattern's centre of gravity with the dictionary's errors for a group, and
///   the pattern's centre of gravity with half-pitch errors for an invalid
///   pattern ID; and the pattern's row and column spans, the dictionary's or
///   the stream's. The closure in `o2-gtrk-cluster-closure` holds them to
///   O2's result.
/// - `toClusters` maps them into our sensor frame (Placement.h) and gives each
///   its readout frame as its time box (TimeEstPS.h), in O2's cluster order, so
///   that a cluster's index here is its index in O2's cluster arrays and the
///   MC labels carry over by index.
///
/// ALICE 3's format, `o2::trkft3::TRKCluster` and `FT3Cluster`, is TRK.h's;
/// it takes the same two steps.

#ifndef GTRK_ALICE_CLUSTERS_H_
#define GTRK_ALICE_CLUSTERS_H_

#include "GenericTracking/Base/TimeEstPS.h"
#include "GenericTracking/Clustering/Cluster.h"

#include "CommonConstants/LHCConstants.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace o2::itsmft
{
class CompClusterExt;
class ROFRecord;
class TopologyDictionary;
} // namespace o2::itsmft

namespace o2::gtrk::alice
{

/// A compact cluster's position and errors in ALICE's sensor-local frame.
struct LocalCluster {
  std::uint32_t chip = 0;
  float x = 0.F;       ///< ALICE local x, across the rows [cm]
  float z = 0.F;       ///< ALICE local z, along the columns [cm]
  float sigma2X = 0.F; ///< [cm^2]
  float sigma2Z = 0.F;
  std::uint32_t nPixels = 0;
  std::uint8_t rowSpan = 0; ///< the pattern's rows
  std::uint8_t colSpan = 0; ///< the pattern's columns
};

/// Every compact cluster's local position and errors, in order, reading the
/// patterns of the group and invalid ones as it goes. Throws
/// std::runtime_error for a pattern ID beyond \a dictionary; patterns that
/// run out before a cluster's fail gsl::span's bounds check.
[[nodiscard]] std::vector<LocalCluster> localClusters(std::span<const o2::itsmft::CompClusterExt> clusters, std::span<const unsigned char> patterns,
                                                      const o2::itsmft::TopologyDictionary& dictionary);

/// How long a detector's readout frames are.
struct ReadoutFrames {
  int lengthBC = 0;
};

/// Our clusters from the local ones, in the same order. \a tfFirstOrbit is
/// the timeframe's first orbit, which a cluster's time is measured from.
/// Throws std::runtime_error if a ROF does not cover its clusters exactly and
/// in order, or starts before the timeframe.
[[nodiscard]] std::vector<Cluster> toClusters(std::span<const LocalCluster> local, std::span<const o2::itsmft::ROFRecord> rofs, ReadoutFrames frames,
                                              std::uint32_t tfFirstOrbit);

/// The readout frames as \a rofs lay them out: O2 writes one record per frame
/// of the timeframe, empty or not, so their spacing is the frame length. That
/// is the length the data were taken with, whatever a configuration dump says:
/// in the MC the ITS's comes from CCDB (ITS/Config/AlpideParam, 198 BC) while
/// the workflows' .ini files print the parameter's default, 891. Throws
/// std::runtime_error unless there are two frames or more, evenly spaced.
[[nodiscard]] ReadoutFrames readoutFrames(std::span<const o2::itsmft::ROFRecord> rofs);

namespace detail
{

/// `toClusters` for any format's readout frames, our sensor's origin at
/// \a xOrigin in ALICE's local x.
template <typename ROFRecord>
[[nodiscard]] std::vector<Cluster> toClusters(std::span<const LocalCluster> local, std::span<const ROFRecord> rofs, ReadoutFrames frames,
                                              std::uint32_t tfFirstOrbit, double xOrigin)
{
  if (frames.lengthBC <= 0) {
    throw std::runtime_error("readout frame length not set");
  }
  std::vector<Cluster> result(local.size());
  std::size_t next = 0;
  for (const auto& rof : rofs) {
    const auto& ir = rof.getBCData();
    if (ir.orbit < tfFirstOrbit) {
      throw std::runtime_error(std::format("a readout frame starts in orbit {}, before the timeframe's {}", ir.orbit, tfFirstOrbit));
    }
    if (std::cmp_not_equal(rof.getFirstEntry(), next) || rof.getNEntries() < 0 || next + static_cast<std::size_t>(rof.getNEntries()) > local.size()) {
      throw std::runtime_error(std::format("readout frame at orbit {} bc {} does not continue the clusters at {}", ir.orbit, ir.bc, next));
    }
    const auto startBC = (static_cast<std::int64_t>(ir.orbit - tfFirstOrbit) * o2::constants::lhc::LHCMaxBunches) + ir.bc;
    const auto box = frameBox(startBC, startBC + frames.lengthBC, 0);
    for (int i = 0; i < rof.getNEntries(); ++i, ++next) {
      const auto& in = local[next];
      auto& out = result[next];
      out.chipID = in.chip;
      out.xLocal = static_cast<float>(static_cast<double>(in.x) - xOrigin);
      out.yLocal = in.z;
      out.sigmaXLocal = std::sqrt(in.sigma2X);
      out.sigmaYLocal = std::sqrt(in.sigma2Z);
      out.timePS = box.start;
      out.widthPS = box.width;
      out.rowSpan = in.rowSpan;
      out.colSpan = in.colSpan;
    }
  }
  if (next != local.size()) {
    throw std::runtime_error(std::format("the readout frames cover {} of {} clusters", next, local.size()));
  }
  return result;
}

/// `readoutFrames` for any format's readout frames.
template <typename ROFRecord>
[[nodiscard]] ReadoutFrames readoutFrames(std::span<const ROFRecord> rofs)
{
  if (rofs.size() < 2) {
    throw std::runtime_error("fewer than two readout frames: their length cannot be read off");
  }
  const auto spacing = [&rofs](std::size_t i) { return rofs[i + 1].getBCData().differenceInBC(rofs[i].getBCData()); };
  const auto length = spacing(0);
  for (std::size_t i = 1; i + 1 < rofs.size(); ++i) {
    if (spacing(i) != length) {
      throw std::runtime_error(std::format("readout frames {} and {} are {} BC apart, frames 0 and 1 {}", i, i + 1, spacing(i), length));
    }
  }
  if (length <= 0) {
    throw std::runtime_error("readout frames out of order");
  }
  return {.lengthBC = static_cast<int>(length)};
}

} // namespace detail

} // namespace o2::gtrk::alice

#endif // GTRK_ALICE_CLUSTERS_H_

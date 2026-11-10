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

#ifndef GTRK_WORKFLOW_ADAPTORHELPERS_H_
#define GTRK_WORKFLOW_ADAPTORHELPERS_H_

#include "ITSMFTInput.h"
#include "GenericTracking/Alice/Clusters.h"
#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Workflow/TrackLabel.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"

#include <algorithm>
#include <array>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace o2::gtrk::workflow::adaptor
{

/// A particle's label, the fake flag unless \a good; O2's noise label for no
/// particle.
[[nodiscard]] inline o2::MCCompLabel labelOf(std::uint64_t particle, bool good) { return particle == kNoTrackParticle ? o2::MCCompLabel(true) : labelOfParticle(particle, !good); }

/// A timeframe's clusters as the adaptor needs them: each one's readout frame
/// (its record in CLUSTERSROF) and that frame's start in BC from the
/// timeframe's start, without the delay and bias.
struct ClusterFrames {
  std::vector<std::uint32_t> record;  ///< per cluster
  std::vector<std::int64_t> recordBC; ///< per record
};

[[nodiscard]] inline ClusterFrames clusterFrames(std::span<const o2::itsmft::ROFRecord> rofs, std::size_t clusters, std::uint32_t tfFirstOrbit)
{
  ClusterFrames result;
  result.record.assign(clusters, std::numeric_limits<std::uint32_t>::max());
  result.recordBC.reserve(rofs.size());
  std::size_t next = 0U;
  for (std::size_t r = 0; r < rofs.size(); ++r) {
    if (rofs[r].getFirstEntry() < 0 || static_cast<std::size_t>(rofs[r].getFirstEntry()) != next ||
        rofs[r].getNEntries() < 0 || static_cast<std::size_t>(rofs[r].getNEntries()) > clusters - next) {
      throw std::runtime_error("cluster ROFs do not cover the input in order");
    }
    next += static_cast<std::size_t>(rofs[r].getNEntries());
    const auto& ir = rofs[r].getBCData();
    result.recordBC.push_back(((static_cast<std::int64_t>(ir.orbit) - tfFirstOrbit) * o2::constants::lhc::LHCMaxBunches) + ir.bc);
    for (int i = 0; i < rofs[r].getNEntries(); ++i) {
      result.record[static_cast<std::size_t>(rofs[r].getFirstEntry()) + static_cast<std::size_t>(i)] = static_cast<std::uint32_t>(r);
    }
  }
  if (next != clusters) {
    throw std::runtime_error("cluster ROFs do not cover the input");
  }
  return result;
}

/// A track's time as ITS takes it (`TrackerTraits`): the intersection of its
/// clusters' frames `[start + delay + bias, + length)` in BC from the
/// timeframe's start; where they do not intersect -- a track across frames,
/// which only a margin lets the tracker make -- the record's box, the frames
/// widened by the margin. [lower, upper) in BC.
struct BCRange {
  std::int64_t lower = 0;
  std::int64_t upper = 0;
};

template <typename Record, typename LayerOf>
[[nodiscard]] BCRange trackTime(const Record& record, std::span<const std::uint32_t> clusterIndex, const ClusterFrames& frames, std::span<const FrameTiming> timings, const LayerOf& layerOf)
{
  BCRange nominal{.lower = std::numeric_limits<std::int64_t>::min(), .upper = std::numeric_limits<std::int64_t>::max()};
  for (std::uint8_t i = 0; i < record.count; ++i) {
    const auto row = record.clusters.index[i];
    const auto& timing = timings[layerOf(row)];
    const auto start = frames.recordBC[frames.record[clusterIndex[row]]] + timing.firstStartBC;
    nominal.lower = std::max(nominal.lower, start);
    nominal.upper = std::min(nominal.upper, start + timing.lengthBC);
  }
  if (nominal.lower < nominal.upper) {
    return nominal;
  }
  return {.lower = bcOfPS(record.timePS), .upper = bcOfPS(static_cast<std::int64_t>(record.timePS) + record.widthPS)};
}

/// The earliest frame (its CLUSTERSROF record) among a track's clusters.
template <typename Record>
[[nodiscard]] std::uint32_t earliestRecord(const Record& record, std::span<const std::uint32_t> clusterIndex, const ClusterFrames& frames)
{
  std::uint32_t earliest = std::numeric_limits<std::uint32_t>::max();
  for (std::uint8_t i = 0; i < record.count; ++i) {
    earliest = std::min(earliest, frames.record[clusterIndex[record.clusters.index[i]]]);
  }
  return earliest;
}

/// Each cluster's pixel count, the converter's decoding (`localClusters`).
[[nodiscard]] inline std::vector<std::uint32_t> pixelCounts(std::span<const o2::itsmft::CompClusterExt> compact, std::span<const unsigned char> patterns, const o2::itsmft::TopologyDictionary& dictionary)
{
  const auto local = alice::localClusters(compact, patterns, dictionary);
  std::vector<std::uint32_t> result(local.size());
  std::ranges::transform(local, result.begin(), [](const auto& cluster) { return cluster.nPixels; });
  return result;
}

/// The track's particle by the shared rule (TrackLabel.h), its rows'
/// clusters' labels, innermost first; \a carried gets, per row, whether it
/// carries that particle.
template <typename Record>
[[nodiscard]] o2::MCCompLabel trackLabel(const Record& record, std::span<const std::uint32_t> clusterIndex, const ITSMFTInput& input, std::span<bool> carried)
{
  constexpr std::size_t kRows = std::extent_v<std::remove_cvref_t<decltype(record.clusters.index)>>;
  std::array<RowParticles, kRows> rows{};
  std::array<std::uint32_t, kRows> track{};
  for (std::uint8_t i = 0; i < record.count; ++i) {
    const auto all = input.labels(clusterIndex[record.clusters.index[i]]);
    rows[i] = rowParticles(std::span<const o2::MCCompLabel>{all.data(), all.size()});
    track[i] = i;
  }
  const std::span<const std::uint32_t> order{track.data(), record.count};
  const auto majority = majorityParticle(rows, order);
  const bool good = majority.particle != kNoTrackParticle && majority.carried == record.count;
  for (std::uint8_t i = 0; i < record.count; ++i) {
    carried[i] = majority.particle != kNoTrackParticle && rows[i].carries(majority.particle);
  }
  return labelOf(majority.particle, good);
}

/// Flatten ROF identities in the same layer order as the converter's clusters.
[[nodiscard]] inline ClusterFrames clusterFrames(const ITSMFTInput& input, std::uint32_t orbit)
{
  ClusterFrames result;
  result.record.reserve(input.size());
  for (const auto& stream : input.streams()) {
    auto frames = clusterFrames(stream.rofs, stream.compact.size(), orbit);
    const auto offset = static_cast<std::uint32_t>(result.recordBC.size());
    for (const auto record : frames.record) {
      result.record.push_back(offset + record);
    }
    result.recordBC.insert(result.recordBC.end(), frames.recordBC.begin(), frames.recordBC.end());
  }
  return result;
}

[[nodiscard]] inline std::vector<std::uint32_t> pixelCounts(const ITSMFTInput& input, const o2::itsmft::TopologyDictionary& dictionary)
{
  std::vector<std::uint32_t> result;
  result.reserve(input.size());
  for (const auto& stream : input.streams()) {
    if (!stream.compact.empty()) {
      const auto pixels = pixelCounts(stream.compact, stream.patterns, dictionary);
      result.insert(result.end(), pixels.begin(), pixels.end());
    }
  }
  return result;
}

/// O2's output clock is the layer with the most frames per timeframe. Equal
/// lengths retain the first layer, independently of the input's occupancy.
[[nodiscard]] inline std::size_t clockLayer(std::span<const FrameTiming> timings)
{
  if (timings.empty() || std::ranges::any_of(timings, [](const auto& timing) { return timing.lengthBC <= 0; })) {
    throw std::invalid_argument("output clock needs positive per-layer frame lengths");
  }
  return static_cast<std::size_t>(std::ranges::min_element(timings, {}, &FrameTiming::lengthBC) - timings.begin());
}

/// Frame records carry the delay, whereas frame intervals also carry the bias.
/// Normalise a delayed record across an orbit boundary.
template <typename Records>
inline void fillClockRecords(Records& result, std::size_t count, int lengthBC, int delayBC, std::uint32_t orbit)
{
  if (lengthBC <= 0 || delayBC < 0) {
    throw std::invalid_argument("invalid output clock record");
  }
  result.resize(count);
  for (std::size_t i = 0; i < count; ++i) {
    const auto bc = static_cast<std::int64_t>(i) * lengthBC + delayBC;
    auto& rof = result[i];
    rof.setBCData(o2::InteractionRecord(static_cast<std::uint16_t>(bc % o2::constants::lhc::LHCMaxBunches),
                                      orbit + static_cast<std::uint32_t>(bc / o2::constants::lhc::LHCMaxBunches)));
    rof.setROFrame(static_cast<std::uint32_t>(i));
    rof.setNEntries(0);
    rof.setFirstEntry(-1);
  }
}

} // namespace o2::gtrk::workflow::adaptor

#endif // GTRK_WORKFLOW_ADAPTORHELPERS_H_

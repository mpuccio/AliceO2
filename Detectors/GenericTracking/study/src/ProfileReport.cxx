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

#include "GenericTracking/Study/ProfileReport.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <format>
#include <ostream>
#include <string>
#include <string_view>
#include <vector>

namespace o2::gtrk::study
{

namespace
{

using enum profile::Phase;

/// The chunk duration below which \a fraction of the chunks lie, from the
/// histogram: the upper edge of the bin that crosses it [s].
[[nodiscard]] double chunkQuantile(const profile::PhaseStats& stats, double fraction)
{
  if (stats.chunks == 0U) {
    return 0.;
  }
  const auto target = fraction * static_cast<double>(stats.chunks);
  std::uint64_t below = 0U;
  for (std::size_t bin = 0; bin < profile::kChunkBins; ++bin) {
    below += stats.chunkHistogram[bin];
    if (static_cast<double>(below) >= target) {
      return std::ldexp(1., static_cast<int>(bin) + 1) * 1.e-9;
    }
  }
  return stats.maxChunkSeconds;
}

} // namespace

void report(std::ostream& out, const profile::Profile& profile, profile::Phase root, int threads)
{
  const auto& top = profile[root];
  if (top.calls == 0U) {
    return;
  }
  const auto parentOf = [](std::size_t index) { return info(static_cast<profile::Phase>(index)).parent; };
  // Per phase: its own serial time -- what neither a child, a parallel loop
  // nor an opaque parallel call covers -- and that summed over its subtree;
  // children follow their parents in the enum, so a reverse pass sums them.
  std::array<double, profile::kPhaseCount> childSeconds{};
  std::array<std::uint64_t, profile::kPhaseCount> childFaults{};
  std::array<std::uint64_t, profile::kPhaseCount> childBytes{};
  for (std::size_t index = 0; index < profile::kPhaseCount; ++index) {
    if (parentOf(index) != kNone) {
      const auto parent = static_cast<std::size_t>(parentOf(index));
      childSeconds[parent] += profile.phases[index].seconds;
      childFaults[parent] += profile.phases[index].minorFaults;
      childBytes[parent] += profile.phases[index].allocatedBytes;
    }
  }
  std::array<double, profile::kPhaseCount> ownSerial{};
  std::array<double, profile::kPhaseCount> subtreeSerial{};
  for (std::size_t index = 0; index < profile::kPhaseCount; ++index) {
    const auto& stats = profile.phases[index];
    ownSerial[index] = std::max(0., stats.seconds - childSeconds[index] - stats.regionSeconds - stats.opaqueSeconds);
    subtreeSerial[index] = ownSerial[index];
  }
  for (std::size_t index = profile::kPhaseCount; index-- > 0U;) {
    if (parentOf(index) != kNone) {
      subtreeSerial[static_cast<std::size_t>(parentOf(index))] += subtreeSerial[index];
    }
  }
  const auto depthOf = [&](std::size_t index) {
    int depth = 0;
    for (auto phase = static_cast<profile::Phase>(index); phase != root; phase = parentOf(static_cast<std::size_t>(phase))) {
      ++depth;
    }
    return depth;
  };

  out << std::format("{:<36} {:>8} {:>6} {:>6} {:>8} {:>11} {:>8} {:>5} {:>5} {:>7} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}\n", "phase", "wall s", "share", "calls", "serial s", "elements", "ns/elem", "util", "imbal", "tail s", "chunk50", "chunk99", "chunkmax", "faults", "alloc MB", "peak MB", "+peak MB");
  // The parallel columns describe the loops opened directly in a phase: on a
  // leaf's row, or on the `(other)` row under a parent.
  const auto parallelColumns = [](const profile::PhaseStats& stats) {
    if (stats.regions == 0U || stats.busySeconds <= 0. || stats.regionSeconds <= 0.) {
      return std::format("{:>11} {:>8} {:>5} {:>5} {:>7} {:>8} {:>8} {:>8}", "", "", "", "", "", "", "", "");
    }
    const double concurrency = stats.capacitySeconds / stats.regionSeconds;
    return std::format("{:>11} {:>8.1f} {:>5.2f} {:>5.2f} {:>7.4f} {:>8.1e} {:>8.1e} {:>8.1e}", stats.elements, stats.elements != 0U ? 1.e9 * stats.busySeconds / static_cast<double>(stats.elements) : 0., stats.busySeconds / stats.capacitySeconds, stats.maxSlotSeconds * concurrency / stats.busySeconds, stats.tailSeconds / concurrency, chunkQuantile(stats, 0.5), chunkQuantile(stats, 0.99), stats.maxChunkSeconds);
  };
  const auto row = [&](std::string_view label, int depth, double seconds, double parentSeconds, std::uint64_t calls, double serial, const std::string& parallel, std::uint64_t faults, std::uint64_t bytes, const std::string& peak) {
    const std::string indented = std::string(static_cast<std::size_t>(2 * depth), ' ') + std::string{label};
    out << std::format("{:<36} {:>8.4f} {:>5.1f}% {:>6} {:>8.4f} {} {:>8} {:>8.1f} {}\n", indented, seconds, parentSeconds > 0. ? 100. * seconds / parentSeconds : 100., calls, serial, parallel, faults, static_cast<double>(bytes) / 1.e6, peak);
  };
  // A peak is the phase's own; the `(other)` row has none.
  const auto peakColumn = [](const profile::PhaseStats& stats) { return std::format("{:>8.1f} {:>8.1f}", static_cast<double>(stats.peakBytes) / 1.e6, static_cast<double>(stats.peakGrowthBytes) / 1.e6); };

  double idle = 0.;
  double tails = 0.;
  double opaque = 0.;
  // The rows in tree order -- a phase, its `(other)` row, then its children
  // in enum order -- since a child need not follow its parent closely in the
  // enum (the paths' phases, the road fills): a depth-first walk, the
  // children pushed last first so the first comes out first.
  std::vector<std::size_t> pending{static_cast<std::size_t>(root)};
  while (!pending.empty()) {
    const auto index = pending.back();
    pending.pop_back();
    const auto& stats = profile.phases[index];
    if (stats.calls != 0U) {
      const double parentSeconds = static_cast<profile::Phase>(index) == root ? stats.seconds : profile[parentOf(index)].seconds;
      const bool parent = childSeconds[index] > 0.;
      row(info(static_cast<profile::Phase>(index)).name, depthOf(index), stats.seconds, parentSeconds, stats.calls, subtreeSerial[index], parent ? parallelColumns(profile::PhaseStats{}) : parallelColumns(stats), stats.minorFaults, stats.allocatedBytes, peakColumn(stats));
      if (parent) {
        row("(other)", depthOf(index) + 1, stats.seconds - childSeconds[index], stats.seconds, stats.calls, ownSerial[index], parallelColumns(stats), stats.minorFaults - std::min(stats.minorFaults, childFaults[index]), stats.allocatedBytes - std::min(stats.allocatedBytes, childBytes[index]), std::format("{:>8} {:>8}", "", ""));
      }
      idle += stats.capacitySeconds - stats.busySeconds;
      tails += stats.tailSeconds;
      opaque += stats.opaqueSeconds;
    }
    for (std::size_t child = profile::kPhaseCount; child-- > index + 1U;) {
      if (parentOf(child) == static_cast<profile::Phase>(index)) {
        pending.push_back(child);
      }
    }
  }
  const double wall = top.seconds;
  const double serial = subtreeSerial[static_cast<std::size_t>(root)];
  const auto perThread = [&](double threadSeconds) { return threads > 0 ? threadSeconds / threads : 0.; };
  out << std::format(
    "{}: {:.4f} s; serial {:.4f} s ({:.1f} %); opaque parallel {:.4f} s; parallel loops idle {:.4f} s of {} threads "
    "({:.1f} %), of which tails {:.4f} s\n",
    info(root).name, wall, serial, wall > 0. ? 100. * serial / wall : 0., opaque, perThread(idle), threads, wall > 0. ? 100. * perThread(idle) / wall : 0., perThread(tails));
}

void writeJson(std::ostream& out, const profile::Profile& profile)
{
  out << '{';
  bool first = true;
  for (std::size_t phase = 0; phase < profile::kPhaseCount; ++phase) {
    const auto& stats = profile.phases[phase];
    if (stats.calls == 0U) {
      continue;
    }
    out << (first ? "" : ",");
    first = false;
    out << std::format(R"("{}":{{"parent":"{}","seconds":{:.9f},"calls":{},"elements":{},"outputs":{},"regions":{},"chunks":{},)"
                       R"("regionSeconds":{:.9f},"capacitySeconds":{:.9f},"busySeconds":{:.9f},"maxSlotSeconds":{:.9f},"tailSeconds":{:.9f},)"
                       R"("maxChunkSeconds":{:.9f},"opaqueSeconds":{:.9f},"minorFaults":{},"dramBytes":{},"allocations":{},"allocatedBytes":{},"peakBytes":{},"peakGrowthBytes":{},"chunkHistogram":[)",
                       info(static_cast<profile::Phase>(phase)).name, info(static_cast<profile::Phase>(phase)).parent == kNone ? "" : info(info(static_cast<profile::Phase>(phase)).parent).name, stats.seconds, stats.calls, stats.elements, stats.outputs, stats.regions, stats.chunks, stats.regionSeconds, stats.capacitySeconds, stats.busySeconds, stats.maxSlotSeconds, stats.tailSeconds, stats.maxChunkSeconds, stats.opaqueSeconds, stats.minorFaults, stats.dramBytes, stats.allocations, stats.allocatedBytes, stats.peakBytes, stats.peakGrowthBytes);
    for (std::size_t bin = 0; bin < profile::kChunkBins; ++bin) {
      out << (bin == 0U ? "" : ",") << stats.chunkHistogram[bin];
    }
    out << "]}";
  }
  out << '}';
}

} // namespace o2::gtrk::study

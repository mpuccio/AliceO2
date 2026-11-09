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

#include "GenericTracking/Tracking/Profile.h"

#include "GenericTracking/Base/BoundedAllocator.h"

#if defined(__linux__)
#include <linux/perf_event.h>
#include <sys/syscall.h>
#endif
#include <sys/resource.h>
#include <unistd.h>

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

namespace o2::gtrk::profile
{

namespace
{

using enum Phase;

constexpr std::array<PhaseInfo, kPhaseCount> kInfo{{
  {.name = "run"},
  {.name = "plan"},
  {.name = "begin timeframe", .parent = kRun},
  {.name = "iteration", .parent = kRun},
  {.name = "slices", .parent = kIteration},
  {.name = "tracklets", .parent = kSlices},
  {.name = "tracklet setup", .parent = kTracklets},
  {.name = "tracklet reset", .parent = kTracklets},
  {.name = "tracklet generate", .parent = kTracklets},
  {.name = "tracklet scan", .parent = kTracklets},
  {.name = "tracklet place", .parent = kTracklets},
  {.name = "tracklet sort", .parent = kTracklets},
  {.name = "triplets", .parent = kSlices},
  {.name = "permutation", .parent = kTriplets},
  {.name = "triplet setup", .parent = kTriplets},
  {.name = "triplet reset", .parent = kTriplets},
  {.name = "triplet generate", .parent = kTriplets},
  {.name = "triplet scan", .parent = kTriplets},
  {.name = "triplet place", .parent = kTriplets},
  {.name = "neighbours", .parent = kSlices},
  {.name = "end index", .parent = kNeighbours},
  {.name = "neighbour setup", .parent = kNeighbours},
  {.name = "neighbour reset", .parent = kNeighbours},
  {.name = "neighbour generate", .parent = kNeighbours},
  {.name = "neighbour scan", .parent = kNeighbours},
  {.name = "neighbour place", .parent = kNeighbours},
  {.name = "roads", .parent = kIteration},
  {.name = "road count", .parent = kRoads},
  {.name = "road fills", .parent = kNeighbourSetup},
  {.name = "road starts", .parent = kRoadCount},
  {.name = "road fit setup", .parent = kRoads},
  {.name = "road fit reset", .parent = kRoads},
  {.name = "road fit generate", .parent = kRoads},
  {.name = "road fit scan", .parent = kRoads},
  {.name = "road fit place", .parent = kRoads},
  {.name = "selection append", .parent = kRoads},
  {.name = "selection", .parent = kIteration},
  {.name = "selection keys", .parent = kSelection},
  {.name = "selection order", .parent = kSelection},
  {.name = "selection deduplicate", .parent = kSelection},
  {.name = "selection rank", .parent = kSelection},
  {.name = "selection fill", .parent = kSelection},
  {.name = "selection slots", .parent = kSelection},
  {.name = "selection rounds", .parent = kSelection},
  {.name = "selection compact", .parent = kSelection},
  {.name = "refit", .parent = kIteration},
  {.name = "follow earlier", .parent = kIteration},
  {.name = "follow", .parent = kIteration},
  {.name = "follow propose", .parent = kFollow},
  {.name = "follow resolve", .parent = kFollow},
  {.name = "follow refit", .parent = kFollow},
  {.name = "follow mark", .parent = kFollow},
  {.name = "observers", .parent = kIteration},
  {.name = "tracklet paths", .parent = kTriplets},
  {.name = "triplet paths", .parent = kTriplets},
  {.name = "length cut", .parent = kIteration},
}};

/// Every parent precedes its children, so one pass in enum order is a
/// depth-first walk once children are grouped under their parent.
constexpr bool parentsFirst()
{
  for (std::size_t phase = 0; phase < kPhaseCount; ++phase) {
    if (kInfo[phase].parent != kNone && static_cast<std::size_t>(kInfo[phase].parent) >= phase) {
      return false;
    }
  }
  return true;
}
static_assert(parentsFirst());

// The profiled thread's context: what a zone deep inside a stage fills without
// every call in between carrying a profile.
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
thread_local Context* tCurrent = nullptr;

// The memory controllers' counters, once `countDram` opened them.
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::vector<int> gDramCounters;

[[nodiscard]] std::uint64_t dramBytes()
{
  std::uint64_t commands = 0U;
  for (const int counter : gDramCounters) {
    std::uint64_t value = 0U;
    if (read(counter, &value, sizeof(value)) == static_cast<ssize_t>(sizeof(value))) {
      commands += value;
    }
  }
  return 64U * commands;
}

[[nodiscard]] std::uint64_t minorFaults()
{
  rusage usage{};
  getrusage(RUSAGE_SELF, &usage);
  // glibc declares the field inside an anonymous union.
  return static_cast<std::uint64_t>(usage.ru_minflt); // NOLINT(cppcoreguidelines-pro-type-union-access)
}

[[nodiscard]] std::size_t chunkBin(std::int64_t nanoseconds)
{
  const auto value = static_cast<std::uint64_t>(std::max<std::int64_t>(nanoseconds, 1));
  return std::min<std::size_t>(static_cast<std::size_t>(std::bit_width(value)) - 1U, kChunkBins - 1U);
}

[[nodiscard]] double seconds(std::int64_t nanoseconds) { return static_cast<double>(nanoseconds) * 1.e-9; }

} // namespace

const PhaseInfo& info(Phase phase)
{
  static const PhaseInfo kNoPhase{.name = "none"};
  return phase < kCount ? kInfo[static_cast<std::size_t>(phase)] : kNoPhase;
}

PhaseStats& PhaseStats::operator+=(const PhaseStats& other)
{
  seconds += other.seconds;
  calls += other.calls;
  elements += other.elements;
  outputs += other.outputs;
  regions += other.regions;
  chunks += other.chunks;
  regionSeconds += other.regionSeconds;
  capacitySeconds += other.capacitySeconds;
  busySeconds += other.busySeconds;
  maxSlotSeconds += other.maxSlotSeconds;
  tailSeconds += other.tailSeconds;
  maxChunkSeconds = std::max(maxChunkSeconds, other.maxChunkSeconds);
  for (std::size_t bin = 0; bin < kChunkBins; ++bin) {
    chunkHistogram[bin] += other.chunkHistogram[bin];
  }
  opaqueSeconds += other.opaqueSeconds;
  minorFaults += other.minorFaults;
  dramBytes += other.dramBytes;
  allocations += other.allocations;
  allocatedBytes += other.allocatedBytes;
  peakBytes = std::max(peakBytes, other.peakBytes);
  peakGrowthBytes = std::max(peakGrowthBytes, other.peakGrowthBytes);
  return *this;
}

Profile& Profile::operator+=(const Profile& other)
{
  for (std::size_t phase = 0; phase < kPhaseCount; ++phase) {
    phases[phase] += other.phases[phase];
  }
  return *this;
}

bool countDram()
{
  if (!gDramCounters.empty()) {
    return true;
  }
  // One PMU per memory controller; its event 0x0a counts CAS commands, read
  // (mask 1) and write (mask 2) ones apart. System-wide, on the CPU the PMU names.
  for (int umc = 0;; ++umc) {
    const std::filesystem::path pmu{std::format("/sys/bus/event_source/devices/amd_umc_{}", umc)};
    if (!std::filesystem::exists(pmu)) {
      break;
    }
    std::uint32_t type = 0U;
    int cpu = 0;
    std::ifstream{pmu / "type"} >> type;
    std::ifstream{pmu / "cpumask"} >> cpu;
#if defined(__linux__)
    for (const std::uint64_t mask : {1U, 2U}) {
      perf_event_attr attr{};
      attr.size = sizeof(attr);
      attr.type = type;
      attr.config = 0x0aU | (mask << 8U);
      const auto counter = static_cast<int>(syscall(SYS_perf_event_open, &attr, -1, cpu, -1, 0UL)); // NOLINT(cppcoreguidelines-pro-type-vararg): perf_event_open has no libc wrapper
      if (counter < 0) {
        for (const int open : gDramCounters) {
          close(open);
        }
        gDramCounters.clear();
        return false;
      }
      gDramCounters.push_back(counter);
    }
#else
    // perf_event_open is Linux-only: no DRAM counters elsewhere.
    static_cast<void>(type);
    static_cast<void>(cpu);
    return false;
#endif
  }
  return !gDramCounters.empty();
}

std::int64_t now() noexcept { return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

Trace::Trace(std::size_t slots) : mSlots(slots) {}

void Trace::chunk(int slot, const Event& event)
{
  if (std::cmp_greater_equal(slot, 0) && std::cmp_less(slot, mSlots.size())) {
    mSlots[static_cast<std::size_t>(slot)].push_back(event);
  }
}

void Trace::clear()
{
  mZones.clear();
  for (auto& slot : mSlots) {
    slot.clear();
  }
}

void Trace::write(std::ostream& out) const
{
  std::int64_t origin = std::numeric_limits<std::int64_t>::max();
  for (const auto& event : mZones) {
    origin = std::min(origin, event.begin);
  }
  for (const auto& slot : mSlots) {
    for (const auto& event : slot) {
      origin = std::min(origin, event.begin);
    }
  }
  // Chrome's format is microseconds; tid 0 is the calling thread's zones,
  // tid s + 1 arena slot s.
  bool first = true;
  const auto emit = [&](const Event& event, std::size_t tid, std::string_view category) {
    out << (first ? "\n" : ",\n") << std::format(R"({{"name":"{}","cat":"{}","ph":"X","pid":1,"tid":{},"ts":{:.3f},"dur":{:.3f}}})", info(event.phase).name, category, tid, static_cast<double>(event.begin - origin) * 1.e-3, static_cast<double>(event.end - event.begin) * 1.e-3);
    first = false;
  };
  const auto name = [&](std::size_t tid, std::string_view label) {
    out << (first ? "\n" : ",\n") << std::format(R"({{"name":"thread_name","ph":"M","pid":1,"tid":{},"args":{{"name":"{}"}}}})", tid, label);
    first = false;
  };
  out << R"({"displayTimeUnit":"ms","traceEvents":[)";
  name(0U, "tracker");
  for (const auto& event : mZones) {
    emit(event, 0U, "zone");
  }
  for (std::size_t slot = 0; slot < mSlots.size(); ++slot) {
    if (mSlots[slot].empty()) {
      continue;
    }
    name(slot + 1U, std::format("slot {}", slot));
    for (const auto& event : mSlots[slot]) {
      emit(event, slot + 1U, "chunk");
    }
  }
  out << "\n]}\n";
}

Context* current() noexcept { return tCurrent; }

Activation::Activation(Profile& profile, BoundedMemoryResource* pool, Trace* trace)
  : mContext{.profile = &profile, .pool = pool, .trace = trace},
    mOuter{tCurrent}
{
  tCurrent = &mContext;
}

Activation::~Activation() { tCurrent = mOuter; }

Zone::Zone(Phase phase)
{
  if (phase >= kCount || tCurrent == nullptr) {
    return;
  }
  mContext = tCurrent;
  mPhase = phase;
  mOuter = mContext->phase;
  mContext->phase = phase;
  mFaults = minorFaults();
  mDram = dramBytes();
  if (mContext->pool != nullptr) {
    mAllocations = mContext->pool->getAllocationCount();
    mAllocatedBytes = mContext->pool->getAllocatedBytes();
    mOuterPeak = mContext->pool->resetPeakMemory();
    mEntryBytes = mContext->pool->getPeakMemory();
  }
  mBegin = now();
}

Zone::~Zone()
{
  if (mContext == nullptr) {
    return;
  }
  const auto end = now();
  auto& stats = (*mContext->profile)[mPhase];
  stats.seconds += seconds(end - mBegin);
  ++stats.calls;
  stats.minorFaults += minorFaults() - mFaults;
  stats.dramBytes += dramBytes() - mDram;
  if (mContext->pool != nullptr) {
    stats.allocations += mContext->pool->getAllocationCount() - mAllocations;
    stats.allocatedBytes += mContext->pool->getAllocatedBytes() - mAllocatedBytes;
    // The window's peak, then the outer window's back: the larger of the two.
    const std::uint64_t peak = mContext->pool->getPeakMemory();
    stats.peakBytes = std::max(stats.peakBytes, peak);
    stats.peakGrowthBytes = std::max(stats.peakGrowthBytes, peak - std::min(peak, mEntryBytes));
    mContext->pool->raisePeakMemory(mOuterPeak);
  }
  if (mContext->trace != nullptr) {
    mContext->trace->zone({.begin = mBegin, .end = end, .phase = mPhase});
  }
  mContext->phase = mOuter;
}

void count(Phase phase, std::uint64_t elements, std::uint64_t outputs)
{
  if (phase >= kCount || tCurrent == nullptr) {
    return;
  }
  auto& stats = (*tCurrent->profile)[phase];
  stats.elements += elements;
  stats.outputs += outputs;
}

void addOpaque(double seconds)
{
  if (tCurrent == nullptr || tCurrent->phase >= kCount) {
    return;
  }
  (*tCurrent->profile)[tCurrent->phase].opaqueSeconds += seconds;
}

void RegionSlots::resize(int slots) { mSlots.assign(static_cast<std::size_t>(std::max(slots, 0)), Slot{}); }

Region::Region(RegionSlots& slots, int concurrency, std::uint64_t elements)
{
  if (tCurrent == nullptr || tCurrent->phase >= kCount) {
    return;
  }
  mContext = tCurrent;
  mTrace = mContext->trace;
  mPhase = mContext->phase;
  mSlots = slots.data();
  mSlotCount = slots.size();
  mConcurrency = std::max(concurrency, 1);
  mElements = elements;
  mBegin = now();
}

void Region::end(int slot, std::int64_t begin) const noexcept
{
  if (mContext == nullptr || std::cmp_less(slot, 0) || std::cmp_greater_equal(slot, mSlotCount)) {
    return;
  }
  const auto end = now();
  const auto duration = end - begin;
  auto& record = mSlots[static_cast<std::size_t>(slot)];
  record.busy += duration;
  record.lastEnd = end;
  record.maxChunk = std::max(record.maxChunk, duration);
  ++record.chunks;
  ++record.histogram[chunkBin(duration)];
  if (mTrace != nullptr) {
    mTrace->chunk(slot, {.begin = begin, .end = end, .phase = mPhase});
  }
}

Region::~Region()
{
  if (mContext == nullptr) {
    return;
  }
  const auto end = now();
  auto& stats = (*mContext->profile)[mPhase];
  const double wall = seconds(end - mBegin);
  ++stats.regions;
  stats.elements += mElements;
  stats.regionSeconds += wall;
  stats.capacitySeconds += wall * mConcurrency;
  std::int64_t busiest = 0;
  for (std::size_t slot = 0; slot < mSlotCount; ++slot) {
    auto& record = mSlots[slot];
    if (record.chunks == 0U) {
      continue;
    }
    stats.chunks += record.chunks;
    stats.busySeconds += seconds(record.busy);
    stats.tailSeconds += seconds(end - record.lastEnd);
    stats.maxChunkSeconds = std::max(stats.maxChunkSeconds, seconds(record.maxChunk));
    for (std::size_t bin = 0; bin < kChunkBins; ++bin) {
      stats.chunkHistogram[bin] += record.histogram[bin];
    }
    busiest = std::max(busiest, record.busy);
    record = RegionSlots::Slot{};
  }
  stats.maxSlotSeconds += seconds(busiest);
}

} // namespace o2::gtrk::profile

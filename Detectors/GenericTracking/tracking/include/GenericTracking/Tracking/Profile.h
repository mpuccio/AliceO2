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

/// \file Profile.h
/// \brief Where the host tracking loop spends its time, phase by phase.
///
/// A `Profile` holds one
/// `PhaseStats` per `Phase`; the phases form a tree, and a phase's time
/// includes its children's, so a report can show what each child leaves
/// unattributed. Three things fill it:
///
/// - a `Zone` times the phase it names on the calling thread, with the minor
///   page faults, the pool allocations inside it and the pool's peak;
/// - a `Region` wraps one parallel loop over chunks: each chunk's duration on
///   its arena slot, folded when the loop ends into the phase open around it
///   as busy time, slot imbalance and tail;
/// - `count` adds elements or emitted records to a phase.
///
/// All three act only while an `Activation` is open on the calling thread, and
/// cost one thread-local read otherwise. Nothing here reads or writes what the
/// tracking computes: a profile changes timings, never output.

#ifndef GTRK_TRACKING_PROFILE_H_
#define GTRK_TRACKING_PROFILE_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string_view>
#include <vector>

namespace o2::gtrk
{
class BoundedMemoryResource;
} // namespace o2::gtrk

namespace o2::gtrk::profile
{

/// Every timed phase of `Tracker::run`, in report order. `kNone` names no
/// phase: a `Zone` of it does nothing.
enum class Phase : std::uint8_t {
  kRun,
  kPlan,
  kBegin,
  kIteration,
  kSlices,
  kTracklets,
  kTrackletSetup,
  kTrackletReset,
  kTrackletGenerate,
  kTrackletScan,
  kTrackletPlace,
  kTrackletSort,
  kTriplets,
  kPermutation,
  kTripletSetup,
  kTripletReset,
  kTripletGenerate,
  kTripletScan,
  kTripletPlace,
  kNeighbours,
  kEndIndex,
  kNeighbourSetup,
  kNeighbourReset,
  kNeighbourGenerate,
  kNeighbourScan,
  kNeighbourPlace,
  kRoads,
  kRoadCount,
  kRoadFills,
  kRoadStarts,
  kRoadFitSetup,
  kRoadFitReset,
  kRoadFitGenerate,
  kRoadFitScan,
  kRoadFitPlace,
  kSelectionAppend,
  kSelection,
  kSelectionKeys,
  kSelectionOrder,
  kSelectionDeduplicate,
  kSelectionRank,
  kSelectionFill,
  kSelectionSlots,
  kSelectionRounds,
  kSelectionCompact,
  kRefit,
  kFollowEarlier,
  kFollow,
  kFollowPropose,
  kFollowResolve,
  kFollowRefit,
  kFollowMark,
  kObservers,
  kTrackletPaths,
  kTripletPaths,
  kLengthCut,
  kCount,
  kNone
};

inline constexpr std::size_t kPhaseCount = static_cast<std::size_t>(Phase::kCount);

/// A phase's name in reports and traces, and the phase that contains it;
/// kNone for a root.
struct PhaseInfo {
  std::string_view name;
  Phase parent = Phase::kNone;
};

[[nodiscard]] const PhaseInfo& info(Phase phase);

/// Chunk durations in powers of two of nanoseconds: bin b holds [2^b, 2^(b+1)).
inline constexpr std::size_t kChunkBins = 40U;

/// What one phase did, summed over every time it was entered.
struct PhaseStats {
  double seconds = 0.;        ///< wall time on the calling thread, children included
  std::uint64_t calls = 0;    ///< times entered
  std::uint64_t elements = 0; ///< elements the parallel loops directly in it visited
  std::uint64_t outputs = 0;  ///< records it emitted, where it emits any
  // The parallel loops opened directly in the phase, not in a child.
  std::uint64_t regions = 0;
  std::uint64_t chunks = 0;
  double regionSeconds = 0.;   ///< their wall time
  double capacitySeconds = 0.; ///< their wall time times the arena's concurrency
  double busySeconds = 0.;     ///< the sum of their chunks' durations
  double maxSlotSeconds = 0.;  ///< per loop, the busiest slot's busy time, summed
  double tailSeconds = 0.;     ///< per loop and working slot, its last chunk's end to the loop's end
  double maxChunkSeconds = 0.;
  std::array<std::uint64_t, kChunkBins> chunkHistogram{};
  /// Parallel work inside a library call that is not cut into our chunks --
  /// the TBB sort. Neither serial nor chunk-attributed.
  double opaqueSeconds = 0.;
  std::uint64_t minorFaults = 0; ///< the process's, while the phase ran
  /// The memory controllers' DRAM traffic while the phase ran, reads and
  /// writes, children included: system-wide, so on a quiet machine
  /// (`countDram`); 0 without.
  std::uint64_t dramBytes = 0;
  std::uint64_t allocations = 0;    ///< from the tracker's pool
  std::uint64_t allocatedBytes = 0; ///< from the tracker's pool
  /// The most of the tracker's pool in use while the phase ran, over its
  /// calls: everything the pool holds, the phase's inputs and what outlives
  /// it included (O2's `stagepeak`). Zones on several threads at once share
  /// the pool's peak, so there a phase's may come out low; the run's is exact.
  std::uint64_t peakBytes = 0;
  /// The most that peak rose above what the pool held when the phase began,
  /// over its calls: what the phase itself added (O2's `stagealloc`).
  std::uint64_t peakGrowthBytes = 0;

  PhaseStats& operator+=(const PhaseStats& other);
};

/// One `PhaseStats` per phase.
struct Profile {
  std::array<PhaseStats, kPhaseCount> phases{};

  [[nodiscard]] PhaseStats& operator[](Phase phase) { return phases[static_cast<std::size_t>(phase)]; }
  [[nodiscard]] const PhaseStats& operator[](Phase phase) const { return phases[static_cast<std::size_t>(phase)]; }
  Profile& operator+=(const Profile& other);
};

/// Nanoseconds on the steady clock.
[[nodiscard]] std::int64_t now() noexcept;

/// Count each phase's DRAM traffic (`PhaseStats::dramBytes`) from now on:
/// the CAS commands of every memory controller, 64 bytes each, from the
/// kernel's AMD UMC counters (Zen 4 and later), read where a zone opens and
/// closes. False, and nothing counted, where they cannot be opened.
[[nodiscard]] bool countDram();

/// A timeline of a run, for Chrome's trace format: the zones on the
/// calling thread and every chunk on its arena slot.
class Trace
{
 public:
  struct Event {
    std::int64_t begin = 0;
    std::int64_t end = 0;
    Phase phase = Phase::kNone;
  };

  /// Room for arena slots [0, \a slots); a chunk on another slot is dropped.
  explicit Trace(std::size_t slots = 1024U);

  /// Events of slot \a slot; each slot is written by one thread at a time.
  void chunk(int slot, const Event& event);
  void zone(const Event& event) { mZones.push_back(event); }
  void clear();

  /// Write the events as a Chrome trace (chrome://tracing, ui.perfetto.dev).
  void write(std::ostream& out) const;

 private:
  std::vector<Event> mZones;
  std::vector<std::vector<Event>> mSlots;
};

/// What a profiled thread is filling.
struct Context {
  Profile* profile = nullptr;
  BoundedMemoryResource* pool = nullptr;
  Trace* trace = nullptr;
  Phase phase = Phase::kNone; ///< the innermost open zone's
};

/// The calling thread's context, or nullptr when it is not profiled.
[[nodiscard]] Context* current() noexcept;

/// Profile the calling thread into \a profile until destroyed; allocations
/// are read from \a pool, events go to \a trace. Nests: the outer activation
/// is restored at the end.
class Activation
{
 public:
  Activation(Profile& profile, BoundedMemoryResource* pool, Trace* trace = nullptr);
  ~Activation();
  Activation(const Activation&) = delete;
  Activation& operator=(const Activation&) = delete;
  Activation(Activation&&) = delete;
  Activation& operator=(Activation&&) = delete;

 private:
  Context mContext;
  Context* mOuter = nullptr;
};

/// Time \a phase from construction to destruction on the calling thread.
class Zone
{
 public:
  explicit Zone(Phase phase);
  ~Zone();
  Zone(const Zone&) = delete;
  Zone& operator=(const Zone&) = delete;
  Zone(Zone&&) = delete;
  Zone& operator=(Zone&&) = delete;

 private:
  Context* mContext = nullptr;
  Phase mPhase = Phase::kNone;
  Phase mOuter = Phase::kNone;
  std::int64_t mBegin = 0;
  std::uint64_t mFaults = 0;
  std::uint64_t mDram = 0;
  std::uint64_t mAllocations = 0;
  std::uint64_t mAllocatedBytes = 0;
  std::uint64_t mOuterPeak = 0;  ///< the pool's peak before the zone opened its window
  std::uint64_t mEntryBytes = 0; ///< the pool's use when it opened it
};

/// Add \a elements and \a outputs to \a phase of the calling thread's profile.
void count(Phase phase, std::uint64_t elements, std::uint64_t outputs);

/// Add \a seconds of opaque parallel work to the open phase.
void addOpaque(double seconds);

/// Per arena slot, what the chunks of one parallel loop did. Owned by
/// whoever runs the loops, sized once to its arena.
class RegionSlots
{
 public:
  struct alignas(64) Slot {
    std::int64_t busy = 0;
    std::int64_t lastEnd = 0;
    std::int64_t maxChunk = 0;
    std::uint64_t chunks = 0;
    std::array<std::uint32_t, kChunkBins> histogram{};
  };

  explicit RegionSlots(int slots = 0) { resize(slots); }
  void resize(int slots);
  [[nodiscard]] std::size_t size() const { return mSlots.size(); }
  [[nodiscard]] Slot* data() { return mSlots.data(); }

 private:
  std::vector<Slot> mSlots;
};

/// One parallel loop of \a elements elements on \a concurrency slots of
/// \a slots. Workers call `begin()` before a chunk and `end(slot, begin)`
/// after it; both do nothing when the loop is not profiled.
class Region
{
 public:
  Region(RegionSlots& slots, int concurrency, std::uint64_t elements);
  ~Region();
  Region(const Region&) = delete;
  Region& operator=(const Region&) = delete;
  Region(Region&&) = delete;
  Region& operator=(Region&&) = delete;

  [[nodiscard]] std::int64_t begin() const noexcept { return mContext != nullptr ? now() : 0; }
  void end(int slot, std::int64_t begin) const noexcept;

 private:
  Context* mContext = nullptr;
  Trace* mTrace = nullptr;
  Phase mPhase = Phase::kNone;
  RegionSlots::Slot* mSlots = nullptr;
  std::size_t mSlotCount = 0U;
  int mConcurrency = 1;
  std::uint64_t mElements = 0U;
  std::int64_t mBegin = 0;
};

/// Which phases a stage's reset, generate, scan and place are.
struct StagePhases {
  Phase reset = Phase::kNone;
  Phase generate = Phase::kNone;
  Phase scan = Phase::kNone;
  Phase place = Phase::kNone;
};

} // namespace o2::gtrk::profile

#endif // GTRK_TRACKING_PROFILE_H_

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

/// \file SerialReferences.h
/// \brief The serial references the tests hold the backends to: a stage run,
///        the tracklet sort, the selection's greedy walk.
///
/// A stage run (`SerialStage`, `runSerial`): elements in order, each one's records staged by (element, ordinal) as a
/// device backend stages them, then the counts scanned and the records placed
/// -- the backends' steps, one at a time, with StageOutput's building blocks.
/// Its failures come in the order the backends report them: an element past
/// `kMaxPerElement` records, then staging too small for the run (the staged
/// count is then the exact size a retry needs), then the scan's.

#ifndef GTRK_TEST_SERIALREFERENCES_H_
#define GTRK_TEST_SERIALREFERENCES_H_

#include "GenericTracking/Tracking/Selection.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/TrackletSort.h"
#include "GenericTracking/Tracking/Triplet.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace o2::gtrk
{

/// Clear \a output's counts, offsets and result, keeping every allocation.
template <typename Output>
void resetOutput(Output& output) noexcept
{
  std::ranges::fill(output.counts(), typename Output::Index{0});
  std::ranges::fill(output.offsets(), typename Output::Index{0});
  output.clearResult();
}

/// Exclusive scan of \a output's counts into its offsets, and the output
/// sized to the total, on the calling thread: the backends' scan, serially.
template <typename Output>
StageStatus serialScan(Output& output)
{
  using Index = typename Output::Index;
  const auto counts = output.counts();
  const auto offsets = output.offsets();
  std::uint64_t total = 0;
  for (std::size_t element = 0; element < counts.size(); ++element) {
    if (counts[element] > Output::kMaxPerElement) {
      return StageStatus::kOrdinalOverflow;
    }
    offsets[element] = static_cast<Index>(total);
    total += counts[element];
    if (total > std::numeric_limits<Index>::max()) {
      return StageStatus::kIndexOverflow;
    }
  }
  offsets[counts.size()] = static_cast<Index>(total);
  output.resizeOutput(static_cast<Index>(total));
  return StageStatus::kComplete;
}

/// Every element in order, on the calling thread: the `forEach` of the host
/// algorithms and the stages' builds, serially.
struct SerialForEach {
  template <typename Body>
  void operator()(std::uint32_t count, const Body& body) const
  {
    for (std::uint32_t i = 0; i < count; ++i) {
      body(i);
    }
  }
};

template <template <template <class> class> class S>
class SerialStage
{
 public:
  using Output = StageOutput<S>;
  using Index = typename Output::Index;
  using Record = typename Output::Record;

  /// What each element's generation hands its records to.
  class Emitter
  {
   public:
    void emit(const Record& record)
    {
      mStage->stage(mElement, mCount, record);
      ++mCount;
    }
    [[nodiscard]] Index count() const { return mCount; }

   private:
    friend class SerialStage;
    Emitter(SerialStage& stage, Index element) : mStage(&stage), mElement(element) {}

    SerialStage* mStage = nullptr;
    Index mElement = 0;
    Index mCount = 0;
  };

  /// Staging for \a capacity records; by default for any run.
  explicit SerialStage(std::size_t capacity = std::numeric_limits<std::size_t>::max()) : mCapacity(capacity) {}

  /// \a generate(element, emitter) for every element of \a output, in order,
  /// into \a output.
  template <typename Generate>
  StageStatus run(Output& output, Generate generate)
  {
    resetOutput(output);
    mStaged.clear();
    mStagedCount = 0;
    mOrdinalOverflow = false;
    const auto counts = output.counts();
    for (Index element = 0; element < output.elementCount(); ++element) {
      Emitter emitter{*this, element};
      generate(element, emitter);
      counts[element] = emitter.count();
    }
    if (mOrdinalOverflow) {
      return StageStatus::kOrdinalOverflow;
    }
    if (mStagedCount > mCapacity) {
      return StageStatus::kStagingOverflow;
    }
    return guardOutOfMemory([&] {
      const auto status = serialScan(output);
      if (status != StageStatus::kComplete) {
        return status;
      }
      const auto offsets = output.offsets();
      for (const auto& staged : mStaged) {
        output.output().write(std::size_t{offsets[staged.element]} + staged.ordinal, staged.record);
      }
      return StageStatus::kComplete;
    });
  }

  /// The records the last run emitted, past the capacity too.
  [[nodiscard]] std::uint64_t stagedCount() const { return mStagedCount; }

 private:
  struct Staged {
    Record record;
    Index element = 0;
    Index ordinal = 0;
  };

  void stage(Index element, Index ordinal, const Record& record)
  {
    const std::uint64_t slot = mStagedCount++;
    if (ordinal >= Output::kMaxPerElement) {
      mOrdinalOverflow = true;
      return;
    }
    if (slot >= mCapacity) {
      return;
    }
    mStaged.push_back({.record = record, .element = element, .ordinal = ordinal});
  }

  std::size_t mCapacity;
  std::vector<Staged> mStaged;
  std::uint64_t mStagedCount = 0;
  bool mOrdinalOverflow = false;
};

/// One serial run of \a generate into \a output, with staging for any count.
template <template <template <class> class> class S, typename Generate>
StageStatus runSerial(StageOutput<S>& output, Generate generate)
{
  return SerialStage<S>{}.run(output, generate);
}

/// The tracklet sort's serial reference for one range: an insertion sort over
/// the columns, quadratic in the range's length.
template <typename View>
void sortTrackletSegment(const View& tracklets, std::uint32_t begin, std::uint32_t end)
{
  for (auto at = begin + 1U; at < end; ++at) {
    const auto source = tracklets.source[at];
    const auto target = tracklets.target[at];
    const auto link = tracklets.link[at];
    const auto lambda = tracklets.lambda[at];
    const auto phi = tracklets.phi[at];
    auto slot = at;
    while (slot > begin && trackletBefore(link, lambda, target, tracklets.link[slot - 1U], tracklets.lambda[slot - 1U], tracklets.target[slot - 1U])) {
      tracklets.source[slot] = tracklets.source[slot - 1U];
      tracklets.target[slot] = tracklets.target[slot - 1U];
      tracklets.link[slot] = tracklets.link[slot - 1U];
      tracklets.lambda[slot] = tracklets.lambda[slot - 1U];
      tracklets.phi[slot] = tracklets.phi[slot - 1U];
      --slot;
    }
    tracklets.source[slot] = source;
    tracklets.target[slot] = target;
    tracklets.link[slot] = link;
    tracklets.lambda[slot] = lambda;
    tracklets.phi[slot] = phi;
  }
}

/// The tracklet sort's serial reference: every element's range of
/// \a output, in turn.
inline void sortTrackletSegments(StageOutput<TrackletT>& output)
{
  const auto view = output.view();
  const auto offsets = output.offsets();
  for (StageOutput<TrackletT>::Index element = 0; element < output.elementCount(); ++element) {
    sortTrackletSegment(view, offsets[element], offsets[element + 1U]);
  }
}

/// The triplets of each permuted element of placed \a triplets, as
/// `TripletPaths::build` reads them (the tracker hands it the staged ones).
inline auto placedTriplets(const StageOutput<TripletT>& triplets)
{
  return [&triplets](StageOutput<TripletT>::Index begin, StageOutput<TripletT>::Index end, const auto& body) {
    const auto view = triplets.view();
    for (auto element = begin; element < end; ++element) {
      for (auto triplet = triplets.offsets()[element]; triplet < triplets.offsets()[element + 1U]; ++triplet) {
        body(element, view.first[triplet], view.second[triplet]);
      }
    }
  };
}

/// The sequential greedy walk the rounds must match, over the
/// records in the order \a ranked gives, best first: accept each that does
/// not conflict with \a mask, and mark it.
template <std::size_t K, typename View>
[[nodiscard]] std::vector<std::uint32_t> selectSequentially(const View& records, std::span<const std::uint32_t> ranked, std::uint32_t sharedInnermost,
                                                            std::span<std::uint8_t> mask)
{
  std::vector<std::uint32_t> selected;
  for (const auto record : ranked) {
    if (selection::conflicts(records.clusters[record], records.count[record], mask.data(), sharedInnermost)) {
      continue;
    }
    selected.push_back(record);
    for (std::uint32_t position = 0; position < records.count[record]; ++position) {
      mask[records.clusters[record].index[position]] = selection::shareable(position, sharedInnermost) ? kShared : kExclusive;
    }
  }
  return selected;
}

} // namespace o2::gtrk

#endif // GTRK_TEST_SERIALREFERENCES_H_

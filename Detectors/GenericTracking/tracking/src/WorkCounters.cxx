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

#include "GenericTracking/Tracking/WorkCounters.h"

#ifdef GTRK_WORK_COUNTERS

#include <algorithm>
#include <mutex>
#include <vector>

namespace o2::gtrk
{

namespace
{

/// Every live thread's block, and what the threads that ended had counted.
struct ThreadBlock;
struct Registry {
  std::mutex mutex;
  std::vector<const ThreadBlock*> live;
  WorkCounts ended{};
  LinkWorkCounts endedLinks;
};

/// Never destroyed: a worker thread may end, and add its counts, after the
/// statics are gone.
Registry& registry()
{
  static auto* instance = new Registry; // NOLINT(cppcoreguidelines-owning-memory)
  return *instance;
}

void add(WorkCounts& sum, const WorkCounts& counts)
{
  for (std::size_t counter = 0; counter < kWorkCounterCount; ++counter) {
    sum[counter] += counts[counter];
  }
}

void add(LinkWorkCounts& sum, const LinkWorkCounts& counts)
{
  if (sum.size() < counts.size()) {
    sum.resize(counts.size());
  }
  for (std::size_t link = 0; link < counts.size(); ++link) {
    for (std::size_t counter = 0; counter < kLinkCounterCount; ++counter) {
      sum[link][counter] += counts[link][counter];
    }
  }
}

/// One thread's block, in the registry for as long as the thread lives.
struct ThreadBlock {
  WorkCounts counts{};
  LinkWorkCounts links;

  ThreadBlock()
  {
    const std::lock_guard lock{registry().mutex};
    registry().live.push_back(this);
  }
  ThreadBlock(const ThreadBlock&) = delete;
  ThreadBlock& operator=(const ThreadBlock&) = delete;
  ThreadBlock(ThreadBlock&&) = delete;
  ThreadBlock& operator=(ThreadBlock&&) = delete;
  ~ThreadBlock()
  {
    const std::lock_guard lock{registry().mutex};
    add(registry().ended, counts);
    add(registry().endedLinks, links);
    std::erase(registry().live, this);
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
thread_local ThreadBlock tBlock;

} // namespace

WorkCounts& detail::threadWorkCounts() { return tBlock.counts; }
LinkWorkCounts& detail::threadLinkWorkCounts() { return tBlock.links; }

WorkCounts workCounts()
{
  const std::lock_guard lock{registry().mutex};
  WorkCounts sum = registry().ended;
  for (const auto* block : registry().live) {
    add(sum, block->counts);
  }
  return sum;
}

LinkWorkCounts linkWorkCounts()
{
  const std::lock_guard lock{registry().mutex};
  LinkWorkCounts sum = registry().endedLinks;
  for (const auto* block : registry().live) {
    add(sum, block->links);
  }
  return sum;
}

} // namespace o2::gtrk

#endif

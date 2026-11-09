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

#include "GenericTracking/Tracking/Parallel.h"

#include "GenericTracking/Tracking/Profile.h"

#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_sort.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/task_scheduler_observer.h>

#if defined(__linux__)
#include <pthread.h>
#include <sched.h>
#endif

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <memory>
#include <vector>

namespace o2::gtrk::detail
{

namespace
{

#if defined(__linux__)
/// One thread per CPU: every thread entering the arena is bound
/// to the CPU of its arena slot, slot i to the i-th CPU of the process's
/// affinity mask, so no two of the arena's threads share a CPU. Left to the
/// kernel, a thread woken for a parallel loop is placed on a busy CPU of its
/// cache domain often enough to stall a chunk for tens of milliseconds. Only
/// when the arena has as many threads as the mask has CPUs,
/// i.e. the job owns its CPUs; with a wider mask other jobs share them, and
/// fixed CPUs would collide with theirs. A thread that is not the arena's
/// worker -- the caller of `parallelFor` -- gets its own mask back when it
/// leaves; a worker stays bound until it enters again, maybe on another slot.
class OneThreadPerCpu final : public tbb::task_scheduler_observer
{
 public:
  OneThreadPerCpu(tbb::task_arena& arena, std::vector<int> cpus) : tbb::task_scheduler_observer{arena}, mCpus{std::move(cpus)} { observe(true); }
  ~OneThreadPerCpu() override { observe(false); }
  OneThreadPerCpu(const OneThreadPerCpu&) = delete;
  OneThreadPerCpu& operator=(const OneThreadPerCpu&) = delete;
  OneThreadPerCpu(OneThreadPerCpu&&) = delete;
  OneThreadPerCpu& operator=(OneThreadPerCpu&&) = delete;

  /// The CPUs of the calling thread's affinity mask, in order.
  [[nodiscard]] static std::vector<int> allowedCpus()
  {
    cpu_set_t mask;
    CPU_ZERO(&mask);
    std::vector<int> cpus;
    if (sched_getaffinity(0, sizeof(mask), &mask) == 0) {
      for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
        if (CPU_ISSET(static_cast<std::size_t>(cpu), &mask)) {
          cpus.push_back(cpu);
        }
      }
    }
    return cpus;
  }

  void on_scheduler_entry(bool worker) override
  {
    const int slot = tbb::this_task_arena::current_thread_index();
    if (slot < 0 || static_cast<std::size_t>(slot) >= mCpus.size()) {
      return;
    }
    if (!worker && tCaller.depth++ == 0) {
      tCaller.saved = pthread_getaffinity_np(pthread_self(), sizeof(tCaller.mask), &tCaller.mask) == 0;
    }
    cpu_set_t one;
    CPU_ZERO(&one);
    CPU_SET(static_cast<std::size_t>(mCpus[static_cast<std::size_t>(slot)]), &one);
    pthread_setaffinity_np(pthread_self(), sizeof(one), &one);
  }

  void on_scheduler_exit(bool worker) override
  {
    if (worker || tCaller.depth == 0) {
      return;
    }
    if (--tCaller.depth == 0 && tCaller.saved) {
      pthread_setaffinity_np(pthread_self(), sizeof(tCaller.mask), &tCaller.mask);
    }
  }

 private:
  /// A calling thread's mask before it entered the first of its arenas.
  struct Caller {
    int depth = 0;
    bool saved = false;
    cpu_set_t mask{};
  };
  static thread_local Caller tCaller;

  std::vector<int> mCpus;
};

thread_local OneThreadPerCpu::Caller OneThreadPerCpu::tCaller;
#endif

} // namespace

struct TaskArena::Impl {
  explicit Impl(int maxConcurrency) : arena{maxConcurrency > 0 ? maxConcurrency : tbb::task_arena::automatic}
  {
#if defined(__linux__)
    auto cpus = OneThreadPerCpu::allowedCpus();
    arena.initialize();
    if (!cpus.empty() && static_cast<std::size_t>(arena.max_concurrency()) == cpus.size()) {
      pinning = std::make_unique<OneThreadPerCpu>(arena, std::move(cpus));
    }
#endif
  }

  tbb::task_arena arena;
#if defined(__linux__)
  std::unique_ptr<OneThreadPerCpu> pinning; ///< with as many threads as CPUs
#endif
};

TaskArena::TaskArena(int maxConcurrency) : mImpl{std::make_unique<Impl>(maxConcurrency)} {}

TaskArena::~TaskArena() = default;

int TaskArena::maxConcurrency() const noexcept { return mImpl->arena.max_concurrency(); }

void parallelSortIndices(TaskArena& arena, size_t count, void* context, ParallelSortLess less, uint32_t* first)
{
  const auto begin = profile::now();
  arena.mImpl->arena.execute([=] { tbb::parallel_sort(first, first + count, [context, less](uint32_t lhs, uint32_t rhs) { return less(context, lhs, rhs); }); });
  // Parallel, but not in chunks a profile can see.
  profile::addOpaque(static_cast<double>(profile::now() - begin) * 1.e-9);
}

void parallelFor(TaskArena& arena, size_t begin, size_t end, size_t grainSize, void* context, ParallelForBody body)
{
  if (begin >= end) {
    return;
  }
  // One task per arena slot, each claiming the next grain of [begin, end) from one counter until
  // none is left. The calling thread runs one of the tasks, so it works to the end. Not
  // tbb::parallel_for over a blocked_range: there the calling thread, an external waiter, stops
  // taking chunks after its first failed steal and sleeps.
  arena.mImpl->arena.execute([=, &arena] {
    std::atomic<size_t> next{begin};
    const int slots = arena.mImpl->arena.max_concurrency();
    tbb::parallel_for(
      0, slots, 1,
      [&](int /*slot*/) {
        for (;;) {
          const size_t first = next.fetch_add(grainSize, std::memory_order_relaxed);
          if (first >= end) {
            return;
          }
          body(context, first, std::min(first + grainSize, end));
        }
      },
      tbb::simple_partitioner{});
  });
}

int currentArenaThreadIndex() noexcept { return tbb::this_task_arena::current_thread_index(); }

} // namespace o2::gtrk::detail

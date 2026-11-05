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

/// \file Parallel.h
/// \brief Host parallelism without TBB in any header.
///
/// A parallel loop and a parallel sort, both in a task arena, implemented with
/// TBB in Parallel.cxx. VertexStage.h includes this header and is compiled
/// into the CUDA and HIP units, so TBB here would reach nvcc and hipcc.

#ifndef GTRK_TRACKING_PARALLEL_H_
#define GTRK_TRACKING_PARALLEL_H_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>

namespace o2::gtrk::detail
{

using ParallelForBody = void (*)(void*, size_t, size_t);
using ParallelSortLess = bool (*)(void*, size_t, size_t);

/// A TBB task arena kept behind the tracking library boundary.
class TaskArena
{
 public:
  explicit TaskArena(int maxConcurrency = 0);
  ~TaskArena();
  TaskArena(const TaskArena&) = delete;
  TaskArena& operator=(const TaskArena&) = delete;
  TaskArena(TaskArena&&) = delete;
  TaskArena& operator=(TaskArena&&) = delete;

  [[nodiscard]] int maxConcurrency() const noexcept;

 private:
  struct Impl;
  std::unique_ptr<Impl> mImpl;

  friend void parallelFor(TaskArena& arena, size_t begin, size_t end, size_t grainSize, void* context, ParallelForBody body);
  friend void parallelSortIndices(TaskArena& arena, size_t count, void* context, ParallelSortLess less, uint32_t* first);
};

/// The type-erased body of `parallelSort`.
void parallelSortIndices(TaskArena& arena, size_t count, void* context, ParallelSortLess less, uint32_t* first);

void parallelFor(TaskArena& arena, size_t begin, size_t end, size_t grainSize, void* context, ParallelForBody body);
[[nodiscard]] int currentArenaThreadIndex() noexcept;

/// Sort \a indices in place by \a less on \a arena's threads (a sort outside
/// it would take every core TBB has, not the tracker's threads). The sort is
/// unstable, so \a less has to be a strict total order -- no two distinct
/// indices equivalent -- for the result to be a serial sort's; every caller
/// ends its comparison in an index tie-break.
template <typename Less>
void parallelSort(TaskArena& arena, std::span<std::uint32_t> indices, Less less)
{
  parallelSortIndices(
    arena, indices.size(), static_cast<void*>(&less),
    [](void* context, std::size_t lhs, std::size_t rhs) {
      return (*static_cast<Less*>(context))(static_cast<std::uint32_t>(lhs), static_cast<std::uint32_t>(rhs));
    },
    indices.data());
}

/// Run the loop inside an arena, which limits its concurrency independently
/// of the process-wide TBB scheduler setting; \a grainSize > 0.
template <typename Function>
void parallelFor(TaskArena& arena, std::size_t begin, std::size_t end, std::size_t grainSize, Function function)
{
  parallelFor(arena, begin, end, grainSize, static_cast<void*>(&function),
              [](void* context, std::size_t first, std::size_t last) { (*static_cast<Function*>(context))(first, last); });
}

} // namespace o2::gtrk::detail

#endif // GTRK_TRACKING_PARALLEL_H_

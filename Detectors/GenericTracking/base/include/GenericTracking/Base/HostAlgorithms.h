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

/// \file HostAlgorithms.h
/// \brief An ordered compaction and an exclusive scan on the host threads.
///
/// Both cut [0, n) into chunks of consecutive indices, work on the chunks in
/// parallel, and join them with a scan of the chunks' totals on one thread:
/// one total per chunk. A chunk works through its indices in order, so the
/// result is the serial loop's at any thread count and any schedule.
///
/// \a forEach(count, body) runs body(i) for every i in [0, count), serially or
/// in parallel -- the tracker passes a backend's `forEachElement`. The chunks
/// are its elements.
///
/// Not the standard parallel algorithms: they run outside the tracker's arena,
/// allocate outside the bounded pool, and report no 32-bit overflow.

#ifndef GTRK_BASE_HOSTALGORITHMS_H_
#define GTRK_BASE_HOSTALGORITHMS_H_

#include "GenericTracking/Base/BoundedAllocator.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory_resource>
#include <span>

namespace o2::gtrk::host_algorithms
{

using Index = std::uint32_t;

/// How [0, \a n) is cut: chunks of \a chunkSize indices, the last one shorter.
struct Chunking {
  Index n = 0U;
  Index size = 1U;

  Chunking(Index items, Index chunkSize) : n(items), size(std::max<Index>(1U, chunkSize)) {}
  [[nodiscard]] Index count() const { return n == 0U ? 0U : 1U + ((n - 1U) / size); }
  [[nodiscard]] Index begin(Index chunk) const { return chunk * size; }
  [[nodiscard]] Index end(Index chunk) const { return static_cast<Index>(std::min<std::uint64_t>(n, (std::uint64_t{chunk} + 1U) * size)); }
};

/// The indices in [0, \a n) for which \a keep(i) holds, in increasing order,
/// into \a out, resized to their number. \a keep is called twice per index.
template <typename ForEach, typename Keep>
void compactIndices(const ForEach& forEach, Index n, Index chunkSize, const Keep& keep, BoundedVector<Index>& out)
{
  const Chunking chunks{n, chunkSize};
  BoundedVector<Index> offsets(std::size_t{chunks.count()} + 1U, Index{0}, out.get_allocator());
  forEach(chunks.count(), [&](Index chunk) {
    Index count = 0U;
    for (Index i = chunks.begin(chunk); i < chunks.end(chunk); ++i) {
      count += keep(i) ? 1U : 0U;
    }
    offsets[chunk + 1U] = count;
  });
  for (Index chunk = 0; chunk < chunks.count(); ++chunk) {
    offsets[chunk + 1U] += offsets[chunk];
  }
  resizeBoundedVector(out, offsets[chunks.count()]);
  forEach(chunks.count(), [&](Index chunk) {
    Index at = offsets[chunk];
    for (Index i = chunks.begin(chunk); i < chunks.end(chunk); ++i) {
      if (keep(i)) {
        out[at++] = i;
      }
    }
  });
}

/// `offsets[i]` = the sum of \a value(j) over j < i, for i in [0, \a n];
/// \a offsets holds n + 1 entries. False, with \a offsets unspecified, when the
/// total does not fit an Index -- as the serial loop's running check, since the
/// sums only grow. \a value(i) must not read \a offsets.
template <typename ForEach, typename Value>
[[nodiscard]] bool exclusiveScan(const ForEach& forEach, Index n, Index chunkSize, const Value& value, std::span<Index> offsets, std::pmr::memory_resource* resource)
{
  const Chunking chunks{n, chunkSize};
  BoundedVector<std::uint64_t> totals(std::size_t{chunks.count()} + 1U, std::uint64_t{0}, std::pmr::polymorphic_allocator<std::uint64_t>{resource});
  forEach(chunks.count(), [&](Index chunk) {
    std::uint64_t sum = 0U;
    for (Index i = chunks.begin(chunk); i < chunks.end(chunk); ++i) {
      sum += value(i);
    }
    totals[chunk + 1U] = sum;
  });
  for (Index chunk = 0; chunk < chunks.count(); ++chunk) {
    totals[chunk + 1U] += totals[chunk];
  }
  const std::uint64_t total = totals[chunks.count()];
  if (total > std::numeric_limits<Index>::max()) {
    return false;
  }
  forEach(chunks.count(), [&](Index chunk) {
    auto running = static_cast<Index>(totals[chunk]);
    for (Index i = chunks.begin(chunk); i < chunks.end(chunk); ++i) {
      offsets[i] = running;
      running += value(i);
    }
  });
  offsets[n] = static_cast<Index>(total);
  return true;
}

} // namespace o2::gtrk::host_algorithms

#endif // GTRK_BASE_HOSTALGORITHMS_H_

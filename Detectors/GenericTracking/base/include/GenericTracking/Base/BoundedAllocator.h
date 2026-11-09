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
///
/// \file BoundedAllocator.h
/// \brief A std::pmr::memory_resource that throws instead of exceeding a
///        budget, and the helpers that release and rebuild vectors on it.
///

#ifndef GTRK_BASE_BOUNDEDALLOCATOR_H_
#define GTRK_BASE_BOUNDEDALLOCATOR_H_

#include <array>
#include <atomic>
#include <limits>
#include <memory>
#include <memory_resource>
#include <new>
#include <string>
#include <utility>
#include <vector>

namespace o2::gtrk
{

class BoundedMemoryResource final : public std::pmr::memory_resource
{
 public:
  class MemoryLimitExceeded final : public std::bad_alloc
  {
   public:
    MemoryLimitExceeded(size_t attempted, size_t used, size_t max);
    [[nodiscard]] const char* what() const noexcept final;

   private:
    std::string mMsg;
  };

  /// At most \a maxBytes in use at once: an allocation past it throws
  /// MemoryLimitExceeded. Without \a upstream the shared caching pool serves.
  explicit BoundedMemoryResource(size_t maxBytes = std::numeric_limits<size_t>::max(), std::pmr::memory_resource* upstream = nullptr);

  /// As above, owning its \a upstream.
  explicit BoundedMemoryResource(std::unique_ptr<std::pmr::memory_resource> upstream, size_t maxBytes = std::numeric_limits<size_t>::max());

  [[nodiscard]] size_t getUsedMemory() const noexcept;
  [[nodiscard]] size_t getMaxMemory() const noexcept;
  [[nodiscard]] size_t getThrowCount() const noexcept;
  [[nodiscard]] size_t getPeakMemory() const noexcept;
  [[nodiscard]] size_t getPeakMemoryDelta() const noexcept;

  /// Open a window on the peak: peak and baseline become what is in use now,
  /// and the peak before is returned, for raisePeakMemory to restore when the
  /// window closes (nested stages).
  size_t resetPeakMemory() noexcept;
  /// The peak becomes at least \a peak.
  void raisePeakMemory(size_t peak) noexcept;
  /// Allocations served, and their bytes, since construction.
  [[nodiscard]] size_t getAllocationCount() const noexcept { return mAllocations.load(std::memory_order_relaxed); }
  [[nodiscard]] size_t getAllocatedBytes() const noexcept { return mAllocatedBytes.load(std::memory_order_relaxed); }
  void setMaxMemory(size_t max);

#if !defined(__HIPCC__) && !defined(__CUDACC__)
  [[nodiscard]] std::string asString() const;
  void print() const;
#endif

 private:
  static std::pmr::memory_resource* cachingUpstream();

  void* do_allocate(size_t bytes, size_t alignment) final;
  void do_deallocate(void* p, size_t bytes, size_t alignment) final;
  [[nodiscard]] bool do_is_equal(const std::pmr::memory_resource& other) const noexcept final;

  std::atomic<size_t> mMaxMemory{std::numeric_limits<size_t>::max()};
  std::atomic<size_t> mCountThrow{0};
  std::atomic<size_t> mUsedMemory{0};
  std::atomic<size_t> mAllocations{0};
  std::atomic<size_t> mAllocatedBytes{0};
  std::atomic<size_t> mPeakUsedMemory{0};
  std::atomic<size_t> mPeakBaselineMemory{0};
  std::unique_ptr<std::pmr::memory_resource> mOwnedUpstream;
  std::pmr::memory_resource* mUpstream{nullptr};
};

template <typename T>
using BoundedVector = std::pmr::vector<T>;

template <typename T>
using bounded_vector = BoundedVector<T>; // NOLINT(readability-identifier-naming): ITS's spelling

template <typename T>
inline void deepVectorClear(std::vector<T>& vec)
{
  std::vector<T>().swap(vec);
}

template <typename T>
inline void deepVectorClear(bounded_vector<T>& vec, std::pmr::memory_resource* mr = nullptr)
{
  std::pmr::memory_resource* tmr = (mr != nullptr) ? mr : vec.get_allocator().resource();
  std::destroy_at(&vec);
  new (&vec) bounded_vector<T>(std::pmr::polymorphic_allocator<T>{tmr});
}

template <typename T>
inline void deepVectorClear(std::vector<bounded_vector<T>>& vec, std::pmr::memory_resource* mr = nullptr)
{
  for (auto& v : vec) {
    deepVectorClear(v, mr);
  }
}

template <typename T, size_t S>
inline void deepVectorClear(std::array<bounded_vector<T>, S>& arr, std::pmr::memory_resource* mr = nullptr)
{
  for (size_t i{0}; i < S; ++i) {
    deepVectorClear(arr[i], mr);
  }
}

template <typename T>
inline void clearResizeBoundedVector(bounded_vector<T>& vec, size_t sz, std::pmr::memory_resource* mr = nullptr, const T& def = T())
{
  // Rebuild a valid empty vector before allocation or element construction can
  // throw. The old allocation is released first, as required by bounded pools.
  deepVectorClear(vec, mr);
  vec.assign(sz, def);
}

// \a vec at \a sz elements: in its storage where that holds them (the old
// elements kept, the added ones value-initialised; the caller rewrites the
// array whole before reading it), otherwise
// rebuilt as clearResizeBoundedVector rebuilds it, the old storage released
// first. A non-null different resource requires a rebuild.
template <typename T>
inline void resizeBoundedVector(BoundedVector<T>& vec, size_t sz, std::pmr::memory_resource* mr = nullptr)
{
  if (sz > vec.capacity() || (mr != nullptr && mr != vec.get_allocator().resource())) {
    clearResizeBoundedVector(vec, sz, mr);
    return;
  }
  vec.resize(sz);
}

// Release several bounded vectors onto one resource.
template <typename... T>
inline void deepVectorClear(std::pmr::memory_resource* mr, BoundedVector<T>&... vecs)
{
  (deepVectorClear(vecs, mr), ...);
}

template <typename T>
inline void clearResizeBoundedVector(std::vector<bounded_vector<T>>& vec, size_t size, std::pmr::memory_resource* mr)
{
  vec.clear();
  vec.reserve(size);
  for (size_t i = 0; i < size; ++i) {
    vec.emplace_back(std::pmr::polymorphic_allocator<bounded_vector<T>>{mr});
  }
}

} // namespace o2::gtrk

#endif /* GTRK_BASE_BOUNDEDALLOCATOR_H_ */

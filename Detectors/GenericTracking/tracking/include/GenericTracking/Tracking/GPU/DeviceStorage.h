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

/// \file DeviceStorage.h
/// \brief Typed, non-owning descriptors over DeviceArena storage.

#ifndef GTRK_TRACKING_GPU_DEVICESTORAGE_H_
#define GTRK_TRACKING_GPU_DEVICESTORAGE_H_

#include "GenericTracking/Base/Columns.h"
#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"

#include <memlayout/wrapper.h>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <new>
#include <vector>
#include <memory>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <unistd.h>

namespace o2::gtrk
{

namespace device_storage_detail
{
template <typename T>
[[nodiscard]] DeviceAllocation allocateArray(DeviceArena& arena, std::size_t count, ArenaEnd end)
{
  static_assert(std::is_trivially_copyable_v<T>);
  if (count > std::numeric_limits<std::size_t>::max() / sizeof(T)) {
    return arena.allocate(std::numeric_limits<std::size_t>::max(), 256U, end);
  }
  return arena.allocate(count * sizeof(T), 256U, end);
}

/// \a copy(to, from) for each column of \a from and the same column of \a to,
/// both column views of one record type.
template <typename From, typename To, typename Copy>
void forEachColumnPair(const From& from, const To& to, const Copy& copy)
{
  from.apply([&](const auto&... source) { to.apply([&](auto... target) { (copy(target, source), ...); }); });
}

/// The arena of a descriptor, \a what naming it in the error; with \a tag,
/// also checked live (std::logic_error for storage of a rewound scope).
[[nodiscard]] inline DeviceArena& requireAllocated(DeviceArena* arena, const char* what, const ArenaTag* tag = nullptr)
{
  if (arena == nullptr) {
    throw std::logic_error(std::string{what} + " has not been allocated");
  }
  if (tag != nullptr) {
    arena->requireLive(*tag);
  }
  return *arena;
}
} // namespace device_storage_detail

template <typename T>
class DeviceSpan
{
 public:
  static_assert(std::is_trivially_copyable_v<T>);

  [[nodiscard]] bool allocate(DeviceArena& arena, std::size_t count, ArenaEnd end = ArenaEnd::kBase)
  {
    const DeviceAllocation allocation = device_storage_detail::allocateArray<T>(arena, count, end);
    if (!allocation) {
      return false;
    }
    mArena = &arena;
    mData = static_cast<T*>(allocation.data);
    mSize = count;
    mTag = allocation.tag;
    return true;
  }

  [[nodiscard]] T* data()
  {
    requireLive();
    return mData;
  }
  [[nodiscard]] const T* data() const
  {
    requireLive();
    return mData;
  }
  [[nodiscard]] std::size_t size() const
  {
    requireLive();
    return mSize;
  }
  [[nodiscard]] DeviceArena& arena() const { return device_storage_detail::requireAllocated(mArena, "device span"); }

 private:
  void requireLive() const { static_cast<void>(device_storage_detail::requireAllocated(mArena, "device span", &mTag)); }

  DeviceArena* mArena = nullptr;
  T* mData = nullptr;
  std::size_t mSize = 0;
  ArenaTag mTag;
};

template <typename T>
class DeviceObject
{
 public:
  static_assert(std::is_trivially_copyable_v<T>);

  [[nodiscard]] bool allocate(DeviceArena& arena, ArenaEnd end = ArenaEnd::kBase) { return mStorage.allocate(arena, 1U, end); }
  [[nodiscard]] T* data() { return mStorage.data(); }
  [[nodiscard]] const T* data() const { return mStorage.data(); }
  [[nodiscard]] DeviceArena& arena() const { return mStorage.arena(); }

 private:
  DeviceSpan<T> mStorage;
};

template <template <template <class> class> class S>
class DeviceColumns
{
 public:
  using View = memlayout::Wrapper<S, memlayout::pointer, memlayout::Layout::soa>;
  using ConstView = memlayout::Wrapper<S, memlayout::const_pointer, memlayout::Layout::soa>;
  using Pointers = S<memlayout::pointer>;

  /// Allocate all columns or none. Successful allocations inherit the parent
  /// scope's lifetime; a failure rewinds the temporary scope.
  [[nodiscard]] bool allocate(DeviceArena& arena, std::size_t rows, ArenaEnd end = ArenaEnd::kBase)
  {
    auto transaction = arena.scope(end);
    Pointers pointers{};
    bool complete = true;
    pointers.apply([&](auto&... pointer) {
      const auto allocateOne = [&](auto& out) {
        if (!complete) {
          return;
        }
        using Element = std::remove_pointer_t<std::remove_reference_t<decltype(out)>>;
        const DeviceAllocation allocation = device_storage_detail::allocateArray<Element>(arena, rows, end);
        if (!allocation) {
          complete = false;
          return;
        }
        out = static_cast<Element*>(allocation.data);
      };
      (allocateOne(pointer), ...);
    });
    if (!complete) {
      return false;
    }
    transaction.commit();
    mArena = &arena;
    mPointers = pointers;
    mRows = rows;
    mTag = arena.currentTag(end);
    return true;
  }

  [[nodiscard]] std::size_t size() const
  {
    requireLive();
    return mRows;
  }
  [[nodiscard]] View view()
  {
    requireLive();
    return View{mPointers};
  }
  [[nodiscard]] ConstView view() const
  {
    requireLive();
    return ConstView{mPointers};
  }
  [[nodiscard]] DeviceArena& arena() const { return device_storage_detail::requireAllocated(mArena, "device columns"); }

 private:
  void requireLive() const { static_cast<void>(device_storage_detail::requireAllocated(mArena, "device columns", &mTag)); }

  DeviceArena* mArena = nullptr;
  Pointers mPointers{};
  std::size_t mRows = 0;
  ArenaTag mTag;
};

template <typename T>
void upload(const T* source, DeviceSpan<T>& destination, std::size_t count, DeviceStream& stream)
{
  if (count > destination.size()) {
    throw std::out_of_range("device span upload exceeds its destination");
  }
  if (count == 0U) {
    return;
  }
  destination.arena().runtime().copyToDevice(destination.data(), source, count * sizeof(T), stream);
}

template <typename T>
void download(const DeviceSpan<T>& source, T* destination, std::size_t count, DeviceStream& stream)
{
  if (count > source.size()) {
    throw std::out_of_range("device span download exceeds its source");
  }
  if (count == 0U) {
    return;
  }
  source.arena().runtime().copyToHost(destination, source.data(), count * sizeof(T), stream);
}

/// The same into \a destination, sized to \a count.
template <typename T>
void download(const DeviceSpan<T>& source, std::vector<T>& destination, std::size_t count, DeviceStream& stream)
{
  destination.resize(count);
  download(source, destination.data(), count, stream);
}

template <typename T>
void uploadObject(const T& source, DeviceObject<T>& destination, DeviceStream& stream)
{
  destination.arena().runtime().copyToDevice(destination.data(), &source, sizeof(T), stream);
}

/// Upload \a rows rows of any host column view -- a DPL message's columns
/// included, which are not a BoundedColumns. The record type cannot be deduced
/// through the nested ConstView, so callers name it: `upload<MeasurementT>(...)`.
template <template <template <class> class> class S>
void upload(typename DeviceColumns<S>::ConstView source, std::size_t rows, DeviceColumns<S>& destination, DeviceStream& stream)
{
  if (rows > destination.size()) {
    throw std::out_of_range("device column upload exceeds its destination");
  }
  if (rows == 0U) {
    return;
  }
  auto& runtime = destination.arena().runtime();
  device_storage_detail::forEachColumnPair(source, destination.view(), [&](auto* to, const auto* from) { runtime.copyToDevice(to, from, rows * sizeof(*to), stream); });
}

/// Download \a rows rows of \a source into \a destination from row \a
/// destinationOffset, which must hold them.
template <template <template <class> class> class S>
void download(const DeviceColumns<S>& source, BoundedColumns<S>& destination, std::size_t destinationOffset, std::size_t rows, DeviceStream& stream)
{
  if (rows > source.size() || destinationOffset > destination.size() || rows > destination.size() - destinationOffset) {
    throw std::out_of_range("device column download exceeds its source or its destination");
  }
  if (rows == 0U) {
    return;
  }
  auto& runtime = source.arena().runtime();
  device_storage_detail::forEachColumnPair(source.view(), destination.view(), [&](auto* to, const auto* from) { runtime.copyToHost(to + destinationOffset, from, rows * sizeof(*to), stream); });
}

template <template <template <class> class> class S>
void download(const DeviceColumns<S>& source, BoundedColumns<S>& destination, std::size_t rows, DeviceStream& stream)
{
  if (rows > source.size()) {
    throw std::out_of_range("device column download exceeds its source");
  }
  destination.resize(rows);
  download(source, destination, 0U, rows, stream);
}

namespace device_storage_detail
{
/// The copies of \a rows rows of \a source, from row \a sourceOffset, into
/// \a destination at row \a destinationOffset, a copy per column, appended to
/// \a copies.
template <template <template <class> class> class S>
void columnCopies(const DeviceColumns<S>& source, std::size_t sourceOffset, DeviceColumns<S>& destination, std::size_t destinationOffset, std::size_t rows, std::vector<DeviceCopy>& copies)
{
  if (sourceOffset > source.size() || rows > source.size() - sourceOffset || destinationOffset > destination.size() || rows > destination.size() - destinationOffset) {
    throw std::out_of_range("device column copy exceeds its source or its destination");
  }
  if (rows == 0U) {
    return;
  }
  forEachColumnPair(source.view(), destination.view(), [&](auto* to, const auto* from) { copies.push_back({.to = to + destinationOffset, .from = from + sourceOffset, .bytes = rows * sizeof(*to)}); });
}
} // namespace device_storage_detail

/// Copy \a rows rows of \a source, from row \a sourceOffset, into
/// \a destination at row \a destinationOffset, on the device: all columns
/// by one kernel, queued on \a stream.
template <template <template <class> class> class S>
void copyColumns(const DeviceColumns<S>& source, std::size_t sourceOffset, DeviceColumns<S>& destination, std::size_t destinationOffset, std::size_t rows, DeviceStream& stream)
{
  std::vector<DeviceCopy> copies;
  device_storage_detail::columnCopies(source, sourceOffset, destination, destinationOffset, rows, copies);
  if (!copies.empty()) {
    destination.arena().runtime().copyEachOnDevice(copies, stream);
  }
}

/// Records of a column set kept as chunks, each exactly sized, in the order
/// they came, and gathered into one column set when they are all there.
/// The chunks live in the arena scope that was
/// current on their end when they were appended; gathering them after it
/// closed throws, as any use of stale device storage does.
template <template <template <class> class> class S>
class DeviceColumnChunks
{
 public:
  /// Copy \a rows rows of \a source, from row \a sourceOffset, into a new
  /// chunk on \a end of \a arena. False, and nothing kept, if the arena refuses.
  [[nodiscard]] bool append(DeviceArena& arena, ArenaEnd end, const DeviceColumns<S>& source, std::size_t sourceOffset, std::size_t rows, DeviceStream& stream)
  {
    if (rows == 0U) {
      return true;
    }
    DeviceColumns<S> chunk;
    if (!chunk.allocate(arena, rows, end)) {
      return false;
    }
    copyColumns(source, sourceOffset, chunk, 0U, rows, stream);
    mChunks.push_back(chunk);
    mRows += rows;
    return true;
  }

  /// \a columns, all their rows, as the next chunk, where they were allocated.
  void adopt(const DeviceColumns<S>& columns)
  {
    mChunks.push_back(columns);
    mRows += columns.size();
  }

  /// Copy every chunk, in the order appended, into \a destination from row
  /// \a destinationOffset: all chunks' columns by one kernel.
  void gather(DeviceColumns<S>& destination, std::size_t destinationOffset, DeviceStream& stream) const
  {
    std::vector<DeviceCopy> copies;
    std::size_t at = destinationOffset;
    for (const auto& chunk : mChunks) {
      device_storage_detail::columnCopies(chunk, 0U, destination, at, chunk.size(), copies);
      at += chunk.size();
    }
    if (!copies.empty()) {
      destination.arena().runtime().copyEachOnDevice(copies, stream);
    }
  }

  /// Download every chunk, in the order appended, into \a destination from
  /// row \a destinationOffset, which must hold them: no gather on the device.
  void download(BoundedColumns<S>& destination, std::size_t destinationOffset, DeviceStream& stream) const
  {
    std::size_t at = destinationOffset;
    for (const auto& chunk : mChunks) {
      gtrk::download(chunk, destination, at, chunk.size(), stream);
      at += chunk.size();
    }
  }

  /// Forget the chunks; their storage goes with its scope.
  void clear() noexcept
  {
    mChunks.clear();
    mRows = 0U;
  }

  [[nodiscard]] std::size_t rows() const noexcept { return mRows; }

 private:
  std::vector<DeviceColumns<S>> mChunks;
  std::size_t mRows = 0U;
};

/// Host memory in whole pages of its own, page-aligned: a block a device can
/// page-lock (DeviceRuntime::registerHost) without touching a page that
/// anything else uses.
template <typename T>
struct PageAllocator {
  using value_type = T; // NOLINT(readability-identifier-naming): the allocator requirements' name
  PageAllocator() = default;
  template <typename U>
  explicit PageAllocator(const PageAllocator<U>& /*other*/) noexcept
  {
  }
  [[nodiscard]] T* allocate(std::size_t count)
  {
    const std::size_t page = pageBytes();
    const std::size_t bytes = ((count * sizeof(T)) + page - 1U) / page * page;
    void* memory = std::aligned_alloc(page, std::max(bytes, page));
    if (memory == nullptr) {
      throw std::bad_alloc();
    }
    return static_cast<T*>(memory);
  }
  void deallocate(T* memory, std::size_t /*count*/) noexcept { std::free(memory); } // NOLINT(cppcoreguidelines-no-malloc): aligned_alloc's
  [[nodiscard]] static std::size_t pageBytes() { return static_cast<std::size_t>(sysconf(_SC_PAGESIZE)); }
  friend bool operator==(const PageAllocator& /*a*/, const PageAllocator& /*b*/) noexcept { return true; }
};

/// A FlatObject's bytes as they must land on the device: the object, whose
/// internal pointers already point into the device buffer, and its buffer.
/// Kept between ports of the same source: the clone behind the bytes is
/// re-targeted when the device buffer moves, not made again.
struct FlatObjectImage {
  std::vector<std::byte> object;
  std::vector<std::byte, PageAllocator<std::byte>> buffer; ///< whole pages: DeviceFlatObject page-locks it

  const void* source = nullptr; ///< the object the image is of
  char* future = nullptr;       ///< the device buffer its pointers name
  std::shared_ptr<void> clone;  ///< the source's clone over `buffer`, the object's bytes taken from it
};

/// How a FlatObject type is prepared for a device: a specialisation provides
/// `static std::size_t bufferBytes(const T& source)` and
/// `static void image(const T& source, char* futureBuffer, FlatObjectImage& image)`,
/// which makes \a image of \a source for \a futureBuffer or re-targets the
/// one it holds of the same source.
/// O2 compiles a FlatObject's host methods -- cloning, relocation, even its
/// constructor -- only outside GPU translation units, so a specialisation is
/// declared in a header and defined in a host .cxx, usually as
/// flatObjectBufferBytes<T> and flatObjectImage<T> below; DeviceFlatObject
/// then works in any translation unit.
template <typename T>
struct FlatObjectPort;

/// The size of \a source's flat buffer, for host translation units.
template <typename T>
[[nodiscard]] std::size_t flatObjectBufferBytes(const T& source)
{
  return static_cast<std::size_t>(source.getFlatBufferSize());
}

/// FlatObject.h's recipe for a GPU, for host translation units: \a source
/// cloned into \a image's buffer and the clone's pointers relocated to
/// \a futureBuffer with setFutureBufferAddress, the object's bytes taken.
/// An image of the same source is re-targeted instead: its clone's pointers
/// back to the buffer it lives in (setActualBufferAddress), then to the new
/// address; an image already for \a futureBuffer is left as it is. The same
/// source object is taken to hold the same content: one changed in place
/// needs a new image (a fresh FlatObjectImage).
template <typename T>
void flatObjectImage(const T& source, char* futureBuffer, FlatObjectImage& image)
{
  // FlatObject speaks char*; the image holds std::byte.
  auto* const buffer = reinterpret_cast<char*>(image.buffer.data()); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
  auto* clone = static_cast<T*>(image.clone.get());
  if (clone == nullptr || image.source != &source) {
    image.buffer.resize(static_cast<std::size_t>(source.getFlatBufferSize()));
    auto made = std::make_shared<T>();
    made->cloneFromObject(source, reinterpret_cast<char*>(image.buffer.data())); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    clone = made.get();
    image.clone = std::move(made);
    image.source = &source;
  } else if (image.future == futureBuffer) {
    return;
  } else {
    clone->setActualBufferAddress(buffer);
  }
  clone->setFutureBufferAddress(futureBuffer);
  image.future = futureBuffer;
  image.object.resize(sizeof(T));
  std::memcpy(image.object.data(), static_cast<const void*>(clone), sizeof(T));
}

/// An O2 FlatObject -- a MatLayerCylSet, a calibration object -- ported to the
/// device: its image, made by FlatObjectPort<T>
/// for the device buffer's address, copied as bytes, the buffer first. The
/// object and its buffer are one allocation, so the upload is all or nothing,
/// and the object is usable on the device as it lands. The image stays on the
/// host for the next upload of the same source, only re-targeted where the
/// device buffer moved, its buffer page-locked for the copies while it stays.
template <typename T>
class DeviceFlatObject
{
 public:
  DeviceFlatObject() = default;
  DeviceFlatObject(const DeviceFlatObject&) = delete;
  DeviceFlatObject& operator=(const DeviceFlatObject&) = delete;
  DeviceFlatObject(DeviceFlatObject&&) = delete;
  DeviceFlatObject& operator=(DeviceFlatObject&&) = delete;
  ~DeviceFlatObject() { unlock(); }

  /// Port \a source to \a end of \a arena. False, and the arena as it was, if
  /// the arena refuses. Returns once the copies are done.
  [[nodiscard]] bool upload(DeviceArena& arena, const T& source, DeviceStream& stream, ArenaEnd end = ArenaEnd::kBase)
  {
    constexpr std::size_t kAlignment = 256U;
    const std::size_t objectBytes = (sizeof(T) + kAlignment - 1U) / kAlignment * kAlignment;
    const std::size_t bufferBytes = FlatObjectPort<T>::bufferBytes(source);
    DeviceSpan<std::byte> storage;
    if (!storage.allocate(arena, objectBytes + bufferBytes, end)) {
      return false;
    }
    std::byte* const object = storage.data();
    std::byte* const buffer = object + objectBytes;
    auto& runtime = arena.runtime();
    if (mImage.source != &source) {
      unlock(); // the image is made again: its buffer may go
    }
    FlatObjectPort<T>::image(source, reinterpret_cast<char*>(buffer), mImage); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    if (mImage.object.size() != sizeof(T) || mImage.buffer.size() != bufferBytes) {
      throw std::logic_error("a flat object's image does not match its type and buffer size");
    }
    if (mLocked == nullptr && bufferBytes != 0U) {
      runtime.registerHost(mImage.buffer.data(), bufferBytes);
      mLocked = mImage.buffer.data();
      mRuntime = &runtime;
    } else if (mLocked != nullptr && mLocked != mImage.buffer.data()) {
      throw std::logic_error("a flat object's page-locked image buffer moved");
    }
    if (bufferBytes != 0U) {
      runtime.copyToDevice(buffer, mImage.buffer.data(), bufferBytes, stream);
    }
    runtime.copyToDevice(object, mImage.object.data(), sizeof(T), stream);
    stream.synchronize(); // the image may be re-targeted once this returns
    mStorage = storage;
    return true;
  }

  /// The device copy let go (its arena reset), the host image kept for the
  /// next upload.
  void release() noexcept { mStorage = {}; }

  /// The object on the device.
  [[nodiscard]] const T* data() const { return reinterpret_cast<const T*>(mStorage.data()); } // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)

 private:
  /// The image's buffer released from its page lock.
  void unlock() noexcept
  {
    if (mLocked != nullptr) {
      mRuntime->unregisterHost(mLocked);
      mLocked = nullptr;
    }
  }

  DeviceSpan<std::byte> mStorage;
  FlatObjectImage mImage;
  std::byte* mLocked = nullptr; ///< the image's buffer as page-locked
  DeviceRuntime* mRuntime = nullptr;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GPU_DEVICESTORAGE_H_

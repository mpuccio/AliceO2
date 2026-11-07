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

/// \file DeviceContext.h
/// \brief A GPU device on O2's GPU route: the reconstruction that owns it,
///        the blocks of its pool the device trackers work in, and the
///        backend's runtime.
///
/// The reconstruction (`o2::gpu::GPUReconstruction`) initialises the device
/// and allocates the pool; `GPUChainGenericTracking` sets aside one block of
/// it per tracker; the backend library loaded with the reconstruction's
/// backend makes the runtime (RuntimeFactory.h). A tracker made on the
/// context is destroyed before it.
///
/// Or on the GPU reconstruction of a workflow that owns it (the TPC's), not
/// owned: no blocks of its own; each timeframe borrows the pool's free memory
/// after the TPC (`borrow`). The reconstruction's own processing and the
/// borrow are one timeframe at a time on its thread, as its pool requires.

#ifndef GTRK_TRACKING_GPU_DEVICECONTEXT_H_
#define GTRK_TRACKING_GPU_DEVICECONTEXT_H_

#include "GenericTracking/Tracking/GPU/DeviceArena.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string_view>

namespace o2::gpu
{
class GPUReconstruction;
} // namespace o2::gpu

namespace o2::gtrk
{

class DeviceRuntime;
class GPUChainGenericTracking;

class DeviceContext
{
 public:
  /// \a backend ("cuda" or "hip", in any case) on \a device, with one block
  /// of each of \a blockBytes; std::runtime_error if it cannot be had.
  DeviceContext(std::string_view backend, int device, std::span<const std::size_t> blockBytes);
  /// On \a reconstruction, initialised by its owner: its backend, the device
  /// it selected. std::runtime_error if it cannot lend its pool: not a GPU, or
  /// its memory settings do not give back what is borrowed (individual
  /// allocations, keepDisplayMemory, disableMemoryReuse).
  explicit DeviceContext(o2::gpu::GPUReconstruction& reconstruction);
  ~DeviceContext();
  DeviceContext(const DeviceContext&) = delete;
  DeviceContext& operator=(const DeviceContext&) = delete;
  DeviceContext(DeviceContext&&) = delete;
  DeviceContext& operator=(DeviceContext&&) = delete;

  [[nodiscard]] DeviceRuntime& runtime() const noexcept { return *mRuntime; }
  /// Block \a index, in the order of the constructor's sizes; not on a
  /// borrowed reconstruction (std::logic_error).
  [[nodiscard]] DeviceBlock block(std::size_t index) const;

  /// The pool's free memory, borrowed from the reconstruction until the
  /// borrow goes: a stack allocation of all of it between the pool's push and
  /// pop of \a tag. Empty (no data, no bytes) if nothing is free.
  class Borrowed
  {
   public:
    ~Borrowed();
    Borrowed(const Borrowed&) = delete;
    Borrowed& operator=(const Borrowed&) = delete;
    Borrowed(Borrowed&&) = delete;
    Borrowed& operator=(Borrowed&&) = delete;
    [[nodiscard]] DeviceBlock block() const noexcept { return mBlock; }

   private:
    friend class DeviceContext;
    Borrowed(o2::gpu::GPUReconstruction& reconstruction, std::uint64_t tag);
    o2::gpu::GPUReconstruction* mReconstruction;
    std::uint64_t mTag;
    DeviceBlock mBlock;
  };
  /// The timeframe's device memory on a borrowed reconstruction: give it
  /// back (destroy the borrow) only once nothing reads it any more.
  [[nodiscard]] std::unique_ptr<Borrowed> borrow(std::uint64_t tag) const;

 private:
  /// Declared first, destroyed last: the runtime's streams and the blocks
  /// are the reconstruction's device's. Empty on a borrowed one.
  std::unique_ptr<o2::gpu::GPUReconstruction> mOwned;
  o2::gpu::GPUReconstruction* mReconstruction = nullptr;
  GPUChainGenericTracking* mChain = nullptr;
  std::unique_ptr<DeviceRuntime> mRuntime;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GPU_DEVICECONTEXT_H_

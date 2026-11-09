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

/// \file DeviceRuntime.cu
/// \brief CUDA runtime implementation, hipified verbatim for HIP.

#include <cuda_runtime.h>

#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceTracker.h"

#include "Platform.cuh"

#include "GPUCommonHelpers.h"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>

#ifndef GTRK_BACKEND_NAMESPACE
#error "GTRK_BACKEND_NAMESPACE must name this backend's private namespace"
#endif

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE
{

/// In DeviceTracker<Layout>.cu, one translation unit per layout.
std::unique_ptr<DeviceTracker<detectors::ITS, o2::base::MatLayerCylSet>> makeITSTracker(DeviceRuntime& runtime, DeviceBlock block);
std::unique_ptr<DeviceTracker<detectors::MFT, o2::base::MatLayerCylSet>> makeMFTTracker(DeviceRuntime& runtime, DeviceBlock block);
std::unique_ptr<DeviceTracker<detectors::TRKFT3, o2::base::MatLayerCylSet>> makeTRKFT3Tracker(DeviceRuntime& runtime, DeviceBlock block);

namespace
{

#ifdef GTRK_BACKEND_HIP
constexpr std::string_view kBackend = "HIP";
#else
constexpr std::string_view kBackend = "CUDA";
#endif

/// The copies of one copyEachOnDevice launch, each in units of 2^unitShift
/// bytes, the widest of 16, 8, 4 or 1 dividing both addresses, and the
/// `tail` bytes after the last whole unit.
struct CopyBatch {
  static constexpr std::size_t kCopies = 32U;
  struct Copy {
    std::byte* to;
    const std::byte* from;
    std::size_t units;
    unsigned unitShift;
    unsigned tail;
  };
  Copy copies[kCopies];
};

template <std::size_t Bytes>
struct alignas(Bytes) CopyUnit {
  std::byte bytes[Bytes];
};

template <std::size_t Bytes>
GPUdi() void copyUnit(const CopyBatch::Copy& copy, std::size_t unit)
{
  reinterpret_cast<CopyUnit<Bytes>*>(copy.to)[unit] = reinterpret_cast<const CopyUnit<Bytes>*>(copy.from)[unit];
}

/// Copy blockIdx.y of \a batch, one unit a thread, the first `tail` threads
/// also a byte of the tail.
GPUg() void GTRK_KERNEL_BOUNDS copyEachKernel(CopyBatch batch)
{
  const CopyBatch::Copy& copy = batch.copies[blockIdx.y];
  const std::size_t unit = (std::size_t{blockIdx.x} * blockDim.x) + threadIdx.x;
  if (unit < copy.tail) {
    const std::size_t at = (copy.units << copy.unitShift) + unit;
    copy.to[at] = copy.from[at];
  }
  if (unit >= copy.units) {
    return;
  }
  switch (copy.unitShift) {
    case 4U:
      copyUnit<16U>(copy, unit);
      break;
    case 3U:
      copyUnit<8U>(copy, unit);
      break;
    case 2U:
      copyUnit<4U>(copy, unit);
      break;
    default:
      copyUnit<1U>(copy, unit);
      break;
  }
}

class Runtime final : public DeviceRuntime
{
 public:
  explicit Runtime(int device)
  {
    cudaDeviceProp properties{};
    GPUChkErrS(cudaGetDeviceProperties(&properties, device));
    mProperties.name = properties.name;
#ifdef GTRK_BACKEND_HIP
    mProperties.architecture = properties.gcnArchName;
#else
    mProperties.architecture = "sm_" + std::to_string((properties.major * 10) + properties.minor);
#endif
    mProperties.warpSize = properties.warpSize;
  }

  ~Runtime() final = default;

  [[nodiscard]] const char* backend() const noexcept final { return kBackend.data(); }
  [[nodiscard]] const DeviceProperties& properties() const noexcept final { return mProperties; }

  void copyToDevice(void* to, const void* from, std::size_t bytes, DeviceStream& stream) final { copy(to, from, bytes, cudaMemcpyHostToDevice, stream); }
  void copyToHost(void* to, const void* from, std::size_t bytes, DeviceStream& stream) final { copy(to, from, bytes, cudaMemcpyDeviceToHost, stream); }
  void copyOnDevice(void* to, const void* from, std::size_t bytes, DeviceStream& stream) final { copy(to, from, bytes, cudaMemcpyDeviceToDevice, stream); }

  void copyEachOnDevice(std::span<const DeviceCopy> copies, DeviceStream& stream) final
  {
    // Zero-byte copies take no slot: a launch has at least one unit.
    while (!copies.empty()) {
      CopyBatch batch{};
      std::size_t taken = 0U;
      std::size_t largest = 0U;
      for (; !copies.empty() && taken < CopyBatch::kCopies; copies = copies.subspan(1U)) {
        const DeviceCopy& copy = copies.front();
        if (copy.bytes == 0U) {
          continue;
        }
        const auto bits = reinterpret_cast<std::uintptr_t>(copy.to) | reinterpret_cast<std::uintptr_t>(copy.from) | (std::uintptr_t{1} << 4U);
        const int zeros = std::countr_zero(bits);
        const unsigned shift = zeros >= 2 ? static_cast<unsigned>(zeros) : 0U; // no 2-byte unit
        const std::size_t units = copy.bytes >> shift;
        const auto tail = static_cast<unsigned>(copy.bytes - (units << shift));
        batch.copies[taken] = {
          .to = static_cast<std::byte*>(copy.to),
          .from = static_cast<const std::byte*>(copy.from),
          .units = units,
          .unitShift = shift,
          .tail = tail,
        };
        largest = std::max({largest, units, std::size_t{tail}});
        ++taken;
      }
      if (taken == 0U) {
        return;
      }
      const auto threads = static_cast<std::size_t>(kBlockSize);
      if (largest > std::size_t{std::numeric_limits<std::int32_t>::max()} * threads) {
        throw std::length_error("copyEachOnDevice: a copy of more units than a launch counts");
      }
      const dim3 grid{static_cast<unsigned>((largest + threads - 1U) / threads), static_cast<unsigned>(taken)};
      copyEachKernel<<<grid, threads, 0U, nativeStreamOf(stream)>>>(batch);
      checkLaunch();
    }
  }

  void fill(void* to, std::uint8_t value, std::size_t bytes, DeviceStream& stream) final { GPUChkErrS(cudaMemsetAsync(to, static_cast<int>(value), bytes, nativeStreamOf(stream))); }

  void registerHost(void* host, std::size_t bytes) final { GPUChkErrS(cudaHostRegister(host, bytes, cudaHostRegisterDefault)); }
  void unregisterHost(void* host) noexcept final { GPUChkErrSI(cudaHostUnregister(host)); }

 protected:
  [[nodiscard]] void* createStreamHandle() final
  {
    cudaStream_t stream = nullptr;
    GPUChkErrS(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    return stream;
  }

 private:
  [[nodiscard]] static cudaStream_t nativeStreamOf(DeviceStream& stream) { return static_cast<cudaStream_t>(nativeStream(stream)); }

  static void copy(void* to, const void* from, std::size_t bytes, cudaMemcpyKind kind, DeviceStream& stream) { GPUChkErrS(cudaMemcpyAsync(to, from, bytes, kind, nativeStreamOf(stream))); }

  void destroyStreamHandle(void* handle) noexcept final { GPUChkErrSI(cudaStreamDestroy(static_cast<cudaStream_t>(handle))); }

  void synchronizeStreamHandle(void* handle) final { GPUChkErrS(cudaStreamSynchronize(static_cast<cudaStream_t>(handle))); }

  [[nodiscard]] std::unique_ptr<DeviceTracker<detectors::ITS, o2::base::MatLayerCylSet>> makeITSTracker(DeviceBlock block) final { return GTRK_BACKEND_NAMESPACE::makeITSTracker(*this, block); }
  [[nodiscard]] std::unique_ptr<DeviceTracker<detectors::MFT, o2::base::MatLayerCylSet>> makeMFTTracker(DeviceBlock block) final { return GTRK_BACKEND_NAMESPACE::makeMFTTracker(*this, block); }
  [[nodiscard]] std::unique_ptr<DeviceTracker<detectors::TRKFT3, o2::base::MatLayerCylSet>> makeTRKFT3Tracker(DeviceBlock block) final { return GTRK_BACKEND_NAMESPACE::makeTRKFT3Tracker(*this, block); }

  DeviceProperties mProperties;
};

std::unique_ptr<DeviceRuntime> makeRuntime(int device) { return std::make_unique<Runtime>(device); }

/// When the library is loaded: O2's GPU backend of the same name links it.
[[maybe_unused]] const bool registered = (registerRuntimeFactory(kBackend, &makeRuntime), true);

} // namespace

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE

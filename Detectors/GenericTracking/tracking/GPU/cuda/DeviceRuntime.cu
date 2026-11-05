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

#include <cstddef>
#include <cstdint>
#include <memory>
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

  void fill(void* to, std::uint8_t value, std::size_t bytes, DeviceStream& stream) final
  {
    GPUChkErrS(cudaMemsetAsync(to, static_cast<int>(value), bytes, nativeStreamOf(stream)));
  }

 protected:
  [[nodiscard]] void* createStreamHandle() final
  {
    cudaStream_t stream = nullptr;
    GPUChkErrS(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    return stream;
  }

 private:
  [[nodiscard]] static cudaStream_t nativeStreamOf(DeviceStream& stream) { return static_cast<cudaStream_t>(nativeStream(stream)); }

  static void copy(void* to, const void* from, std::size_t bytes, cudaMemcpyKind kind, DeviceStream& stream)
  {
    GPUChkErrS(cudaMemcpyAsync(to, from, bytes, kind, nativeStreamOf(stream)));
  }

  void destroyStreamHandle(void* handle) noexcept final { GPUChkErrSI(cudaStreamDestroy(static_cast<cudaStream_t>(handle))); }

  void synchronizeStreamHandle(void* handle) final { GPUChkErrS(cudaStreamSynchronize(static_cast<cudaStream_t>(handle))); }

  [[nodiscard]] std::unique_ptr<DeviceTracker<detectors::ITS, o2::base::MatLayerCylSet>> makeITSTracker(DeviceBlock block) final
  {
    return GTRK_BACKEND_NAMESPACE::makeITSTracker(*this, block);
  }
  [[nodiscard]] std::unique_ptr<DeviceTracker<detectors::MFT, o2::base::MatLayerCylSet>> makeMFTTracker(DeviceBlock block) final
  {
    return GTRK_BACKEND_NAMESPACE::makeMFTTracker(*this, block);
  }
  [[nodiscard]] std::unique_ptr<DeviceTracker<detectors::TRKFT3, o2::base::MatLayerCylSet>> makeTRKFT3Tracker(DeviceBlock block) final
  {
    return GTRK_BACKEND_NAMESPACE::makeTRKFT3Tracker(*this, block);
  }

  DeviceProperties mProperties;
};

std::unique_ptr<DeviceRuntime> makeRuntime(int device) { return std::make_unique<Runtime>(device); }

/// When the library is loaded: O2's GPU backend of the same name links it.
[[maybe_unused]] const bool registered = (registerRuntimeFactory(kBackend, &makeRuntime), true);

} // namespace

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE

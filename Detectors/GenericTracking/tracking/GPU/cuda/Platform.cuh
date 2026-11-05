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

/// \file Platform.cuh
/// \brief The few spellings that differ between CUDA and HIP, and the wave
///        width and launch shape, for device headers included by the
///        backends' .cu sources.
///
/// hipify-perl translates the .cu files it is given, not the headers they
/// include, so a header with vendor code writes
/// both branches, as ITS's Utils.h does.
/// Triple-chevron launches and the atomic intrinsics are spelled the same on
/// both and need nothing here.

#ifndef GTRK_TRACKING_GPU_PLATFORM_CUH_
#define GTRK_TRACKING_GPU_PLATFORM_CUH_

#ifdef __HIPCC__
#include <hip/hip_runtime.h>
#include <hipcub/hipcub.hpp>
#else
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#endif

#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"

#include "GPUCommonHelpers.h"
#include "GPUDefParametersWrapper.h"

#include <cstddef>
#include <cstdint>

namespace o2::gtrk::gpu
{

inline constexpr int kWarpSize = GPUCA_WARP_SIZE;
static_assert(kWarpSize > 0, "no GPUCA_GPUTYPE_* with a warp size");

/// Every kernel launches blocks of this size (launchFor); on HIP it is also
/// their launch bounds (GTRK_KERNEL_BOUNDS). Two warps of 32, one wavefront of
/// 64: a block of the divergent generate kernels and fits leaves its
/// multiprocessor when its own warps are done, not with the slowest of a larger
/// block.
inline constexpr int kBlockSize = 64;
static_assert(kBlockSize % kWarpSize == 0);

struct LaunchShape {
  std::uint32_t blocks = 0;
  std::uint32_t threads = static_cast<std::uint32_t>(kBlockSize);
};

/// One thread per element: a capped grid would make a thread loop over
/// several elements, and a warp's lanes drift apart on the divergent bodies
/// (the triplet path sweep's grid is smaller on purpose, PathGlue.cuh).
[[nodiscard]] constexpr LaunchShape launchFor(std::uint32_t elements)
{
  constexpr auto threads = static_cast<std::uint32_t>(kBlockSize);
  return {.blocks = (elements / threads) + static_cast<std::uint32_t>(elements % threads != 0U), .threads = threads};
}

#if defined(__CUDA_ARCH__)
static_assert(kWarpSize == 32);
#elif defined(__HIP_DEVICE_COMPILE__)
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-pragma"
#endif
#if defined(__AMDGCN_WAVEFRONT_SIZE)
static_assert(__AMDGCN_WAVEFRONT_SIZE == kWarpSize);
#endif
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#endif

} // namespace o2::gtrk::gpu

// A kernel's launch bounds: on HIP the block size every launch uses
// (launchFor). Without it hipcc compiles a kernel for 1024 threads, 64 VGPRs
// on gfx906, and with ROCm 7.1.1 the refit's spills come back wrong. Not on
// CUDA, where the same bound makes nvcc's register allocation slower.
#ifdef __HIPCC__
#define GTRK_KERNEL_BOUNDS __launch_bounds__(o2::gtrk::gpu::kBlockSize)
#else
#define GTRK_KERNEL_BOUNDS
#endif

#ifndef GTRK_BACKEND_NAMESPACE
#error "GTRK_BACKEND_NAMESPACE must name this backend's private namespace"
#endif

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE
{

#ifdef __HIPCC__
namespace cubapi = ::hipcub;
using NativeStream = hipStream_t;
#else
namespace cubapi = ::cub;
using NativeStream = cudaStream_t;
#endif

/// The native stream behind \a stream, made current on its runtime's device.
inline NativeStream nativeStream(DeviceRuntime& runtime, DeviceStream& stream) { return static_cast<NativeStream>(runtime.nativeStream(stream)); }

/// What the 64-bit atomics take on both vendors.
using Counter = unsigned long long;

/// A failed launch is a runtime error, not a rung of the overflow ladder.
inline void checkLaunch()
{
#ifdef __HIPCC__
  GPUChkErrS(hipGetLastError());
#else
  GPUChkErrS(cudaGetLastError());
#endif
}

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE

#endif // GTRK_TRACKING_GPU_PLATFORM_CUH_

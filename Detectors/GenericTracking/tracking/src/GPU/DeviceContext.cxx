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

#include "GenericTracking/Tracking/GPU/DeviceContext.h"

#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/RuntimeFactory.h"

#include "GPUO2ExternalUser.h" // first of O2's GPU headers outside GPU/
#include "GPUChainGenericTracking.h"
#include "GPUMemoryResource.h"
#include "GPUReconstruction.h"
#include "GPUSettings.h"

#include <algorithm>
#include <cctype>
#include <format>
#include <stdexcept>
#include <string>

namespace o2::gtrk
{

namespace
{
/// The backend's name as GPUReconstruction::CreateInstance spells it.
std::string upperCase(std::string_view backend)
{
  std::string type{backend};
  std::ranges::transform(type, type.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return type;
}
} // namespace

DeviceContext::DeviceContext(std::string_view backend, int device, std::span<const std::size_t> blockBytes)
  : mOwned(o2::gpu::GPUReconstruction::CreateInstance(upperCase(backend).c_str(), true)),
    mReconstruction(mOwned.get())
{
  const std::string type = upperCase(backend);
  if (mReconstruction == nullptr) {
    throw std::runtime_error(std::format("GPU backend {}: no GPU reconstruction", type));
  }
  mChain = mReconstruction->AddChain<GPUChainGenericTracking>();
  for (const auto bytes : blockBytes) {
    static_cast<void>(mChain->AddBlock(bytes));
  }
  const o2::gpu::GPUSettingsGRP grp;
  o2::gpu::GPUSettingsProcessing processing;
  processing.deviceNum = device;
  // The pool is the blocks the chain reports, and the pinned host pool what
  // the reconstruction keeps there itself.
  processing.forceMemoryPoolSize = 0;
  processing.forceHostMemoryPoolSize = GPUChainGenericTracking::HostMemoryBytes();
  mReconstruction->SetSettings(&grp, nullptr, &processing);
  if (mReconstruction->Init() != 0) {
    throw std::runtime_error(std::format("GPU backend {}: device {} does not initialise", type, device));
  }
  // The backend library registered its runtime when O2's GPU backend, which
  // links it, was loaded.
  const auto factory = runtimeFactory(type);
  if (factory == nullptr) {
    throw std::runtime_error(std::format("GPU backend {}: no GenericTracking runtime (O2's GPUTracking{} does not link GTrkTracking{})", type, type, type));
  }
  mRuntime = factory(device);
}

DeviceContext::DeviceContext(o2::gpu::GPUReconstruction& reconstruction) : mReconstruction(&reconstruction)
{
  if (!reconstruction.IsGPU()) {
    throw std::runtime_error("GenericTracking: the GPU reconstruction to borrow from runs on no GPU");
  }
  const auto& processing = reconstruction.GetProcessingSettings();
  if (processing.memoryAllocationStrategy != o2::gpu::GPUMemoryResource::ALLOCATION_GLOBAL || processing.keepDisplayMemory || processing.disableMemoryReuse) {
    throw std::runtime_error("GenericTracking: the GPU reconstruction's pool cannot lend its memory (it needs the global allocation, without keepDisplayMemory or disableMemoryReuse)");
  }
  const std::string type = o2::gpu::gpudatatypes::DEVICE_TYPE_NAMES[reconstruction.GetDeviceType()];
  const auto factory = runtimeFactory(type);
  if (factory == nullptr) {
    throw std::runtime_error(std::format("GPU backend {}: no GenericTracking runtime (O2's GPUTracking{} does not link GTrkTracking{})", type, type, type));
  }
  mRuntime = factory(reconstruction.GetDeviceId());
}

DeviceContext::~DeviceContext() = default;

DeviceContext::Borrowed::Borrowed(o2::gpu::GPUReconstruction& reconstruction, std::uint64_t tag) : mReconstruction(&reconstruction), mTag(tag)
{
  mReconstruction->PushNonPersistentMemory(mTag);
  // All of the free memory, measured first: the allocation cannot fail. None
  // if nothing is free (a zero-size one still moves the pool's end).
  const std::size_t bytes = mReconstruction->DeviceMemoryGap();
  if (bytes == 0U) {
    return;
  }
  try {
    mBlock = {.data = mReconstruction->AllocateDirectMemory(bytes, o2::gpu::GPUMemoryResource::MEMORY_GPU | o2::gpu::GPUMemoryResource::MEMORY_STACK), .bytes = bytes};
  } catch (...) {
    mReconstruction->PopNonPersistentMemory(o2::gpu::gpudatatypes::RecoStep::GenericTracking, mTag);
    throw;
  }
}

DeviceContext::Borrowed::~Borrowed() { mReconstruction->PopNonPersistentMemory(o2::gpu::gpudatatypes::RecoStep::GenericTracking, mTag); }

std::unique_ptr<DeviceContext::Borrowed> DeviceContext::borrow(std::uint64_t tag) const
{
  if (mOwned != nullptr) {
    throw std::logic_error("GenericTracking: borrow on a GPU reconstruction of its own");
  }
  return std::unique_ptr<Borrowed>{new Borrowed(*mReconstruction, tag)};
}

DeviceBlock DeviceContext::block(std::size_t index) const
{
  if (mChain == nullptr) {
    throw std::logic_error("GenericTracking: no blocks on a borrowed GPU reconstruction");
  }
  const auto block = mChain->GetBlock(index);
  return {.data = block.data, .bytes = block.bytes};
}

} // namespace o2::gtrk

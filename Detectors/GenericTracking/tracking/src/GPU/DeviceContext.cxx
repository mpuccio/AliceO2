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
#include "GPUReconstruction.h"
#include "GPUSettings.h"

#include <algorithm>
#include <cctype>
#include <format>
#include <stdexcept>
#include <string>

namespace o2::gtrk
{

DeviceContext::DeviceContext(std::string_view backend, int device, std::span<const std::size_t> blockBytes)
{
  std::string type{backend};
  std::ranges::transform(type, type.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  mReconstruction.reset(o2::gpu::GPUReconstruction::CreateInstance(type.c_str(), true));
  if (mReconstruction == nullptr) {
    throw std::runtime_error(std::format("GPU backend {}: no GPU reconstruction", type));
  }
  mChain = mReconstruction->AddChain<GPUChainGenericTracking>();
  for (const auto bytes : blockBytes) {
    static_cast<void>(mChain->AddBlock(bytes));
  }
  o2::gpu::GPUSettingsGRP grp;
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

DeviceContext::~DeviceContext() = default;

DeviceBlock DeviceContext::block(std::size_t index) const
{
  const auto block = mChain->GetBlock(index);
  return {.data = block.data, .bytes = block.bytes};
}

} // namespace o2::gtrk

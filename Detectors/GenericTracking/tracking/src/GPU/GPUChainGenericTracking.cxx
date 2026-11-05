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

/// \file GPUChainGenericTracking.cxx

#include "GPUO2ExternalUser.h" // first of O2's GPU headers outside GPU/

#include "GPUChainGenericTracking.h"

#include "GPUConstantMem.h"
#include "GPUReconstruction.h"

#include <stdexcept>

namespace o2::gtrk
{

using o2::gpu::GPUMemoryResource;
using o2::gpu::GPUProcessor;

GPUChainGenericTracking::GPUChainGenericTracking(o2::gpu::GPUReconstruction* rec) : GPUChain(rec) {}

GPUChainGenericTracking::~GPUChainGenericTracking() = default;

void* GPUChainGenericTracking::BlockProcessor::SetPointersBlock(void* mem) { return static_cast<char*>(mem) + mBytes; }

size_t GPUChainGenericTracking::AddBlock(size_t bytes)
{
  auto& pair = mBlocks.emplace_back(std::make_unique<BlockPair>());
  pair->host.mBytes = bytes;
  pair->device.mBytes = bytes;
  return mBlocks.size() - 1;
}

void GPUChainGenericTracking::RegisterPermanentMemoryAndProcessors()
{
  if (!mRec->IsGPU()) {
    return;
  }
  for (auto& pair : mBlocks) {
    pair->host.InitGPUProcessor(mRec, GPUProcessor::PROCESSOR_TYPE_SLAVE);
    pair->device.InitGPUProcessor(mRec, GPUProcessor::PROCESSOR_TYPE_DEVICE, &pair->host);
    // Device only. Custom: the per-event allocation pass, which knows a
    // permanent resource as allocated by its host pointer, leaves it alone.
    pair->resource = mRec->RegisterMemoryAllocation(&pair->host, &BlockProcessor::SetPointersBlock,
                                                    GPUMemoryResource::MEMORY_GPU | GPUMemoryResource::MEMORY_PERMANENT | GPUMemoryResource::MEMORY_CUSTOM,
                                                    "GenericTrackingBlock");
  }
}

void GPUChainGenericTracking::MemorySize(size_t& gpuMem, size_t& pageLockedHostMem)
{
  gpuMem = 0;
  for (const auto& pair : mBlocks) {
    gpuMem += GPUProcessor::nextMultipleOf<o2::gpu::constants::GPU_MEMALIGN>(pair->host.mBytes);
  }
  pageLockedHostMem = HostMemoryBytes();
}

size_t GPUChainGenericTracking::HostMemoryBytes() { return GPUProcessor::nextMultipleOf<o2::gpu::constants::GPU_MEMALIGN>(sizeof(o2::gpu::GPUConstantMem)) + o2::gpu::constants::GPU_MEMALIGN; }

GPUChainGenericTracking::Block GPUChainGenericTracking::GetBlock(size_t index) const
{
  if (index >= mBlocks.size() || mBlocks[index]->resource < 0) {
    throw std::out_of_range("GenericTracking: no such device block");
  }
  const auto& resource = mRec->Res(mBlocks[index]->resource);
  return {.data = resource.PtrDevice(), .bytes = mBlocks[index]->host.mBytes};
}

} // namespace o2::gtrk

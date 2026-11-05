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

/// \file GPUChainGenericTracking.h
/// \brief The device memory of GenericTracking's device trackers.
///
/// A holder, as O2's GPUChainITS: the trackers run outside the chain, on
/// their backend's runtime (RuntimeFactory.h). The chain sets aside one
/// permanent block of the reconstruction's device pool per tracker, of the
/// size its creator asks for before Init, and hands it out after Init.
/// DeviceContext's alone: include GPUO2ExternalUser.h before it.

#ifndef GTRK_TRACKING_GPU_GPUCHAINGENERICTRACKING_H_
#define GTRK_TRACKING_GPU_GPUCHAINGENERICTRACKING_H_

#include "GPUChain.h"
#include "GPUProcessor.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace o2::gtrk
{

class GPUChainGenericTracking final : public o2::gpu::GPUChain
{
 public:
  /// As GPUReconstruction::AddChain makes it.
  explicit GPUChainGenericTracking(o2::gpu::GPUReconstruction* rec);

  /// One tracker's block: device memory of the pool, permanent.
  struct Block {
    void* data = nullptr;
    size_t bytes = 0;
  };

  ~GPUChainGenericTracking() final;
  int32_t Init() final { return 0; }
  int32_t PrepareEvent() final { return 0; }
  int32_t Finalize() final { return 0; }
  int32_t RunChain() final { return 0; }

  void RegisterPermanentMemoryAndProcessors() final;
  void RegisterGPUProcessors() final {}
  /// The blocks on the device; on the host the processors' shadow, which the
  /// device reconstruction keeps in its pinned pool.
  void MemorySize(size_t& gpuMem, size_t& pageLockedHostMem) final;
  /// What MemorySize reports on the host: the pinned pool a reconstruction
  /// holding only this chain needs (forceHostMemoryPoolSize).
  static size_t HostMemoryBytes();

  /// A block of \a bytes, before Init; its index for GetBlock.
  size_t AddBlock(size_t bytes);
  /// Block \a index, after Init.
  Block GetBlock(size_t index) const;

 private:
  /// A block's processor: the host-side one registers the resource, its
  /// device mirror places it (a device allocation needs the linked pair).
  struct BlockProcessor : public o2::gpu::GPUProcessor {
    size_t mBytes = 0;
    void* SetPointersBlock(void* mem);
  };
  struct BlockPair {
    BlockProcessor host;
    BlockProcessor device;
    int16_t resource = -1;
  };
  std::vector<std::unique_ptr<BlockPair>> mBlocks;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GPU_GPUCHAINGENERICTRACKING_H_

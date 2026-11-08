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

/// \file DeviceRuntime.h
/// \brief A GPU backend's runtime: streams, copies, fills and the device
///        trackers, vendor-neutral.
///
/// Made for the device O2's GPU reconstruction initialised, by the factory
/// its backend library registered (RuntimeFactory.h); the device memory is
/// the reconstruction's pool, of which each tracker gets a block.

#ifndef GTRK_TRACKING_GPU_DEVICERUNTIME_H_
#define GTRK_TRACKING_GPU_DEVICERUNTIME_H_

#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceStream.h"
#include "GenericTracking/Tracking/GPU/RuntimeFactory.h"

namespace o2::base
{
class MatLayerCylSet;
}

#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace o2::gtrk
{

/// Declared only: the runtime names the trackers it makes, their users
/// include DeviceTracker.h.
template <DetectorLayout Det, typename Material>
class DeviceTracker;

struct DeviceProperties {
  std::string name;
  std::string architecture;
  int warpSize = 0;
};

/// Every call runs on the device current on the calling thread: the one
/// the reconstruction initialised.
class DeviceRuntime
{
 public:
  DeviceRuntime() = default;
  virtual ~DeviceRuntime() = default;
  DeviceRuntime(const DeviceRuntime&) = delete;
  DeviceRuntime& operator=(const DeviceRuntime&) = delete;
  DeviceRuntime(DeviceRuntime&&) = delete;
  DeviceRuntime& operator=(DeviceRuntime&&) = delete;

  [[nodiscard]] virtual const char* backend() const noexcept = 0;
  [[nodiscard]] virtual const DeviceProperties& properties() const noexcept = 0;

  [[nodiscard]] DeviceStream createStream();

  virtual void copyToDevice(void* to, const void* from, std::size_t bytes, DeviceStream& stream) = 0;
  virtual void copyToHost(void* to, const void* from, std::size_t bytes, DeviceStream& stream) = 0;
  virtual void copyOnDevice(void* to, const void* from, std::size_t bytes, DeviceStream& stream) = 0;
  virtual void fill(void* to, std::uint8_t value, std::size_t bytes, DeviceStream& stream) = 0;
  /// Page-lock the \a bytes of host memory at \a host for this device's
  /// copies, which then run without staging: whole pages no other
  /// registration shares (`PageAllocator`), until unregisterHost.
  virtual void registerHost(void* host, std::size_t bytes) = 0;
  virtual void unregisterHost(void* host) noexcept = 0;

  /// Backend launchers use this to recover cudaStream_t/hipStream_t without
  /// exposing either type in a public header.
  [[nodiscard]] static void* nativeStream(DeviceStream& stream);

  /// The device tracker for \a D, its arena on \a block; std::runtime_error
  /// for a layout the backends are not built for.
  template <DetectorLayout D>
  [[nodiscard]] std::unique_ptr<DeviceTracker<D, o2::base::MatLayerCylSet>> makeTracker(DeviceBlock block)
  {
    if constexpr (std::is_same_v<D, detectors::ITS>) {
      return makeITSTracker(block);
    } else if constexpr (std::is_same_v<D, detectors::MFT>) {
      return makeMFTTracker(block);
    } else if constexpr (std::is_same_v<D, detectors::TRKFT3>) {
      return makeTRKFT3Tracker(block);
    } else {
      throw std::runtime_error(std::string{"GPU backend "} + backend() + ": no device tracker is built for " + D::name);
    }
  }

 protected:
  [[nodiscard]] virtual void* createStreamHandle() = 0;

 private:
  friend class DeviceStream;

  virtual void destroyStreamHandle(void* handle) noexcept = 0;
  virtual void synchronizeStreamHandle(void* handle) = 0;

  [[nodiscard]] virtual std::unique_ptr<DeviceTracker<detectors::ITS, o2::base::MatLayerCylSet>> makeITSTracker(DeviceBlock block) = 0;
  [[nodiscard]] virtual std::unique_ptr<DeviceTracker<detectors::MFT, o2::base::MatLayerCylSet>> makeMFTTracker(DeviceBlock block) = 0;
  [[nodiscard]] virtual std::unique_ptr<DeviceTracker<detectors::TRKFT3, o2::base::MatLayerCylSet>> makeTRKFT3Tracker(DeviceBlock block) = 0;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GPU_DEVICERUNTIME_H_

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

/// \file DeviceStream.h
/// \brief A vendor-neutral, move-only stream handle.

#ifndef GTRK_TRACKING_GPU_DEVICESTREAM_H_
#define GTRK_TRACKING_GPU_DEVICESTREAM_H_

namespace o2::gtrk
{

class DeviceRuntime;

class DeviceStream
{
 public:
  DeviceStream() = default;
  ~DeviceStream();
  DeviceStream(const DeviceStream&) = delete;
  DeviceStream& operator=(const DeviceStream&) = delete;
  DeviceStream(DeviceStream&& other) noexcept;
  DeviceStream& operator=(DeviceStream&& other) noexcept;

  void synchronize();

 private:
  friend class DeviceRuntime;
  DeviceStream(DeviceRuntime& runtime, void* handle) noexcept : mRuntime(&runtime), mHandle(handle) {}
  void reset() noexcept;

  DeviceRuntime* mRuntime = nullptr;
  void* mHandle = nullptr;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GPU_DEVICESTREAM_H_

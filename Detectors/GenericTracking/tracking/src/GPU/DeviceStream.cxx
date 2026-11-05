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

#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"

#include <stdexcept>
#include <utility>

namespace o2::gtrk
{

DeviceStream DeviceRuntime::createStream() { return DeviceStream{*this, createStreamHandle()}; }

void* DeviceRuntime::nativeStream(DeviceStream& stream) { return stream.mHandle; }

DeviceStream::~DeviceStream() { reset(); }

DeviceStream::DeviceStream(DeviceStream&& other) noexcept : mRuntime(std::exchange(other.mRuntime, nullptr)), mHandle(std::exchange(other.mHandle, nullptr)) {}

DeviceStream& DeviceStream::operator=(DeviceStream&& other) noexcept
{
  if (this != &other) {
    reset();
    mRuntime = std::exchange(other.mRuntime, nullptr);
    mHandle = std::exchange(other.mHandle, nullptr);
  }
  return *this;
}

void DeviceStream::reset() noexcept
{
  if (mHandle != nullptr) {
    mRuntime->destroyStreamHandle(mHandle);
  }
  mRuntime = nullptr;
  mHandle = nullptr;
}

void DeviceStream::synchronize()
{
  if (mRuntime == nullptr) {
    throw std::logic_error("cannot synchronize an empty device stream");
  }
  mRuntime->synchronizeStreamHandle(mHandle);
}

} // namespace o2::gtrk

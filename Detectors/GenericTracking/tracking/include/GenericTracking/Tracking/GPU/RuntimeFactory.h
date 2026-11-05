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

/// \file RuntimeFactory.h
/// \brief The GPU backends' runtime factories: each backend library
///        registers its own when it is loaded -- O2's GPU backend of the same
///        name links it -- and DeviceContext takes it for the device the
///        reconstruction initialised.

#ifndef GTRK_TRACKING_GPU_RUNTIMEFACTORY_H_
#define GTRK_TRACKING_GPU_RUNTIMEFACTORY_H_

#include <memory>
#include <string_view>

namespace o2::gtrk
{

class DeviceRuntime;

/// A backend's runtime on \a device.
using RuntimeFactory = std::unique_ptr<DeviceRuntime> (*)(int device);

/// At a backend library's load: \a backend ("CUDA", "HIP") is made by \a factory.
void registerRuntimeFactory(std::string_view backend, RuntimeFactory factory);

/// The factory \a backend registered; nullptr if its library is not loaded.
[[nodiscard]] RuntimeFactory runtimeFactory(std::string_view backend);

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GPU_RUNTIMEFACTORY_H_

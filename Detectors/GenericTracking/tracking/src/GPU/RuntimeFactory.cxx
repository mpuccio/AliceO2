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

#include "GenericTracking/Tracking/GPU/RuntimeFactory.h"

#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace o2::gtrk
{

namespace
{

struct Registry {
  std::mutex mutex;
  std::vector<std::pair<std::string, RuntimeFactory>> factories;
};

Registry& registry()
{
  static Registry instance;
  return instance;
}

} // namespace

void registerRuntimeFactory(std::string_view backend, RuntimeFactory factory)
{
  auto& known = registry();
  const std::lock_guard lock{known.mutex};
  known.factories.emplace_back(std::string{backend}, factory);
}

RuntimeFactory runtimeFactory(std::string_view backend)
{
  auto& known = registry();
  const std::lock_guard lock{known.mutex};
  for (const auto& [name, factory] : known.factories) {
    if (name == backend) {
      return factory;
    }
  }
  return nullptr;
}

} // namespace o2::gtrk

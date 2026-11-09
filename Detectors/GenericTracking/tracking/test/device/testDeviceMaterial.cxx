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

/// \file testDeviceMaterial.cxx
/// \brief O2's material LUT ported to the device with DeviceFlatObject, its
///        lookups against the host's.
///
/// Reads the LUT GTRK_TEST_MATLUT names and passes without
/// checking anything if there is none. Built once per backend. Needs the
/// device.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK device material
#include <boost/test/unit_test.hpp>

#include "DeviceTestFixture.h"
#include "DeviceTestKernels.h"

#include "DetectorsBase/MatLayerCylSet.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace o2::gtrk
{

using namespace device_test;

BOOST_AUTO_TEST_CASE(TheLutGivesTheHostsBudgetsOnTheDevice)
{
  const char* path = std::getenv("GTRK_TEST_MATLUT");
  if (path == nullptr || !std::filesystem::exists(path)) {
    BOOST_TEST_MESSAGE("no material LUT: nothing to compare");
    return;
  }
  const std::unique_ptr<o2::base::MatLayerCylSet> lut{o2::base::MatLayerCylSet::loadFromFile(path)};
  BOOST_REQUIRE(lut != nullptr);
  std::mt19937 random{13U};
  std::uniform_real_distribution<float> transverse{-60.F, 60.F};
  std::uniform_real_distribution<float> longitudinal{-150.F, 150.F};
  std::uniform_real_distribution<float> step{-3.F, 3.F};
  std::vector<std::array<float, 6>> segments(20000U);
  for (auto& segment : segments) {
    const float x = transverse(random);
    const float y = transverse(random);
    const float z = longitudinal(random);
    segment = {x, y, z, x + step(random), y + step(random), z + step(random)};
  }
  Device device;
  std::vector<std::array<float, 3>> onDevice;
  BOOST_REQUIRE((kernels().matBudgets(*device.runtime, device.stream, *device.arena, *lut, segments, onDevice) == StageStatus::kComplete));
  std::size_t differing = 0U;
  std::size_t nonZero = 0U;
  double worst = 0.;
  const auto relative = [](float host, float device) { return host == device ? 0. : std::abs(static_cast<double>(device) - static_cast<double>(host)) / std::abs(static_cast<double>(host)); };
  for (std::size_t i = 0; i < segments.size(); ++i) {
    const auto& s = segments[i];
    const auto host = lut->getMatBudget(s[0], s[1], s[2], s[3], s[4], s[5]);
    nonZero += host.meanX2X0 > 0.F ? 1U : 0U;
    const bool same = std::bit_cast<std::uint32_t>(host.meanRho) == std::bit_cast<std::uint32_t>(onDevice[i][0]) && std::bit_cast<std::uint32_t>(host.meanX2X0) == std::bit_cast<std::uint32_t>(onDevice[i][1]) && std::bit_cast<std::uint32_t>(host.length) == std::bit_cast<std::uint32_t>(onDevice[i][2]);
    differing += same ? 0U : 1U;
    worst = std::max({worst, relative(host.meanX2X0, onDevice[i][1]), relative(host.length, onDevice[i][2])});
  }
  BOOST_TEST_MESSAGE(runtime()->backend() << ": " << segments.size() << " segments, " << nonZero << " through material, " << differing << " differing, worst relative difference " << worst);
  BOOST_TEST(nonZero > 0U);
#ifdef GPUCA_DETERMINISTIC_MODE
  BOOST_TEST(differing == 0U);
#else
  // Fast math moves where a segment crosses a cell of the LUT, by a small
  // fraction of a budget; a port that is wrong is off by far more.
  BOOST_TEST(worst < 1.e-2);
#endif
}

} // namespace o2::gtrk

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

/// \file testTrackerMemoryPool.cxx
/// \brief The tracker's memory pool from GTRKRecoITS.memoryLimitGB.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK tracker memory pool
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Workflow/GridConfig.h"
#include "GenericTracking/Workflow/RecoParam.h"

#include "CommonUtils/ConfigurableParam.h"

#include <cstddef>
#include <limits>
#include <new>

namespace o2::gtrk::workflow
{

BOOST_AUTO_TEST_CASE(NoLimitByDefault)
{
  BOOST_TEST(recoParam<detectors::ITS>().memoryLimitGB == 0.);
  BOOST_TEST(trackerMemoryPool(recoParam<detectors::ITS>())->getMaxMemory() == std::numeric_limits<std::size_t>::max());
}

BOOST_AUTO_TEST_CASE(TheLimitBoundsThePool)
{
  o2::conf::ConfigurableParam::updateFromString("GTRKRecoITS.memoryLimitGB=0.5");
  const auto pool = trackerMemoryPool(recoParam<detectors::ITS>());
  BOOST_TEST(pool->getMaxMemory() == std::size_t{512} * 1024U * 1024U);
  BoundedVector<std::byte> past{pool.get()};
  BOOST_CHECK_THROW(past.resize(std::size_t{513} * 1024U * 1024U), std::bad_alloc);
  BOOST_TEST(pool->getUsedMemory() == 0U);
  o2::conf::ConfigurableParam::updateFromString("GTRKRecoITS.memoryLimitGB=0");
}

} // namespace o2::gtrk::workflow

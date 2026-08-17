// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#define BOOST_TEST_MODULE ITSMFT SurfaceMask
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/SurfaceMask.h"

using o2::itsmft::tracking::LayerMask;
using o2::itsmft::tracking::positionalSurfaceMask;
using o2::itsmft::tracking::LayerId;
using o2::itsmft::tracking::SurfaceMask;

/// Gate 4 C5: focused coverage for the shared positionalSurfaceMask() helper
/// that replaced the equivalent private copies formerly duplicated in
/// the common tracking topology builder and the combined workflow application plan. Each set bit
/// in `layerMask` is a *position* in `orderedSurfaces`, never a numeric
/// comparison against the LayerId values it holds.

BOOST_AUTO_TEST_CASE(positional_surface_mask_full_active_count_identity_ids)
{
  const std::vector<LayerId> orderedSurfaces{LayerId{0}, LayerId{1}, LayerId{2}, LayerId{3}};
  const LayerMask layerMask{0, 1, 3}; // positions 0, 1, 3 set

  const auto result = positionalSurfaceMask(layerMask, orderedSurfaces, static_cast<uint32_t>(orderedSurfaces.size()));

  BOOST_CHECK(result.has(LayerId{0}));
  BOOST_CHECK(result.has(LayerId{1}));
  BOOST_CHECK(!result.has(LayerId{2}));
  BOOST_CHECK(result.has(LayerId{3}));
  BOOST_CHECK_EQUAL(result.count(), 3);
}

BOOST_AUTO_TEST_CASE(positional_surface_mask_partial_active_count_ignores_tail_bits)
{
  // Layer positions 4 and 5 are set but fall outside the active range: they
  // must not contribute to the result even though `orderedSurfaces` has
  // entries there.
  const std::vector<LayerId> orderedSurfaces{LayerId{10}, LayerId{11}, LayerId{12}, LayerId{13}, LayerId{14}, LayerId{15}};
  LayerMask layerMask{};
  layerMask.set(1);
  layerMask.set(4);
  layerMask.set(5);

  const auto result = positionalSurfaceMask(layerMask, orderedSurfaces, 3);

  BOOST_CHECK(!result.has(LayerId{10}));
  BOOST_CHECK(result.has(LayerId{11}));
  BOOST_CHECK(!result.has(LayerId{12}));
  BOOST_CHECK(!result.has(LayerId{14}));
  BOOST_CHECK(!result.has(LayerId{15}));
  BOOST_CHECK_EQUAL(result.count(), 1);
}

BOOST_AUTO_TEST_CASE(positional_surface_mask_zero_active_count_is_always_empty)
{
  const std::vector<LayerId> orderedSurfaces{LayerId{0}, LayerId{1}};
  const LayerMask layerMask{0, 1, 1};

  const auto result = positionalSurfaceMask(layerMask, orderedSurfaces, 0);

  BOOST_CHECK(result.empty());
}

BOOST_AUTO_TEST_CASE(positional_surface_mask_sparse_bits_non_identity_global_ids)
{
  // Non-identity, non-contiguous global LayerIds, e.g. an MFT-like span
  // offset behind an ITS block in a combined catalog.
  const std::vector<LayerId> orderedSurfaces{LayerId{7}, LayerId{8}, LayerId{9}, LayerId{10}, LayerId{11}};
  LayerMask layerMask{};
  layerMask.set(0);
  layerMask.set(2);
  layerMask.set(4);

  const auto result = positionalSurfaceMask(layerMask, orderedSurfaces, static_cast<uint32_t>(orderedSurfaces.size()));

  BOOST_CHECK(result.has(LayerId{7}));
  BOOST_CHECK(!result.has(LayerId{8}));
  BOOST_CHECK(result.has(LayerId{9}));
  BOOST_CHECK(!result.has(LayerId{10}));
  BOOST_CHECK(result.has(LayerId{11}));
  BOOST_CHECK_EQUAL(result.count(), 3);
}

BOOST_AUTO_TEST_CASE(positional_surface_mask_empty_layer_mask_is_always_empty)
{
  const std::vector<LayerId> orderedSurfaces{LayerId{2}, LayerId{4}, LayerId{6}};
  const LayerMask layerMask{};

  const auto result = positionalSurfaceMask(layerMask, orderedSurfaces, static_cast<uint32_t>(orderedSurfaces.size()));

  BOOST_CHECK(result.empty());
}

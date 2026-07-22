// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT NominalSurfaceMaterial
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// Deliberately the only project header this translation unit includes: this
// file's own minimal include list is the narrow-header dependency-boundary
// check. NominalSurfaceMaterial.h must not need any STL container, gsl::span,
// owner, provider, material-query API, or TGeo/geometry dependency, so a TU
// that includes nothing beyond it (and the test framework) must compile and
// link cleanly.
#include "ITSMFTTracking/NominalSurfaceMaterial.h"

using namespace o2::itsmft::tracking;

BOOST_AUTO_TEST_CASE(BudgetDefaultsToZero)
{
  NominalSurfaceMaterialBudget budget{};
  BOOST_CHECK_EQUAL(budget.normalXOverX0, 0.f);
  BOOST_CHECK_EQUAL(budget.normalArealDensityGPerCm2, 0.f);
}

BOOST_AUTO_TEST_CASE(EntryHoldsSurfaceAndBudget)
{
  const NominalSurfaceMaterialEntry entry{SurfaceId{3}, NominalSurfaceMaterialBudget{0.01f, 0.5f}};
  BOOST_CHECK_EQUAL(entry.surface.value(), 3);
  BOOST_CHECK_EQUAL(entry.budget.normalXOverX0, 0.01f);
  BOOST_CHECK_EQUAL(entry.budget.normalArealDensityGPerCm2, 0.5f);
}

BOOST_AUTO_TEST_CASE(BothTypesAreCopyable)
{
  constexpr NominalSurfaceMaterialBudget original{0.02f, 1.5f};
  const NominalSurfaceMaterialBudget copy = original; // NOLINT
  BOOST_CHECK_EQUAL(copy.normalXOverX0, original.normalXOverX0);
  BOOST_CHECK_EQUAL(copy.normalArealDensityGPerCm2, original.normalArealDensityGPerCm2);
}

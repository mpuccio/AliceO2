// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT IndexTableConfigBinding
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <cmath>
#include <limits>

#include <boost/test/unit_test.hpp>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IndexTableConfiguration.h"
#include "ITSMFTTracking/IndexTableUtils.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{
constexpr int ITSN = 7;
constexpr int MFTN = 10;

bool isUntouched(const IndexTableUtils<ITSN>& utils)
{
  // Default-constructed sentinel state (IndexTableUtils.h): never mutated.
  return utils.getNrowBins() == 0 && utils.getNcolBins() == 0 && utils.getCoordType() == IndexTableCoordType::PhiZ;
}

bool isUntouched(const IndexTableUtils<MFTN>& utils)
{
  return utils.getNrowBins() == 0 && utils.getNcolBins() == 0 && utils.getCoordType() == IndexTableCoordType::PhiZ;
}
} // namespace

/// Exact CylinderCylinder/PhiZ parity: bindIndexTableConfiguration must
/// reproduce IndexTableUtils::setTrackingParameters() bit-for-bit for the
/// default ITS-shaped TrackingParameters.
BOOST_AUTO_TEST_CASE(CylinderCylinderBindingMatchesSetTrackingParameters)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::ITS);

  IndexTableUtils<ITSN> staged;
  const auto error = bindIndexTableConfiguration<TransitionPolicyTag::CylinderCylinder, ITSN>(staged, params);
  BOOST_REQUIRE(error == IndexTableConfigError::None);

  IndexTableUtils<ITSN> reference;
  reference.setTrackingParameters(params);

  BOOST_CHECK(staged.getCoordType() == reference.getCoordType());
  BOOST_CHECK_EQUAL(staged.getNrowBins(), reference.getNrowBins());
  BOOST_CHECK_EQUAL(staged.getNcolBins(), reference.getNcolBins());
  BOOST_CHECK_EQUAL(staged.getRowOrigin(), reference.getRowOrigin());
  BOOST_CHECK_EQUAL(staged.getRowCoordinateSpan(), reference.getRowCoordinateSpan());
  for (int iLayer = 0; iLayer < ITSN; ++iLayer) {
    BOOST_CHECK_EQUAL(staged.getLayerColHalfExtent(iLayer), reference.getLayerColHalfExtent(iLayer));
  }
}

/// Exact DiskDisk/XY parity with the +/-20 cm default fallback (IndexRowMax
/// left at 0, as a caller bypassing resetDetectorDefaults could do).
BOOST_AUTO_TEST_CASE(DiskDiskBindingMatchesDefaultRowRangeFallback)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::MFT);
  params.IndexRowMin = 0.f;
  params.IndexRowMax = 0.f; // force the "no explicit range" fallback path

  IndexTableUtils<MFTN> staged;
  const auto error = bindIndexTableConfiguration<TransitionPolicyTag::DiskDisk, MFTN>(staged, params);
  BOOST_REQUIRE(error == IndexTableConfigError::None);

  IndexTableUtils<MFTN> reference;
  reference.setTrackingParametersXY(params, -20.f, 20.f);

  BOOST_CHECK(staged.getCoordType() == reference.getCoordType());
  BOOST_CHECK_EQUAL(staged.getNrowBins(), reference.getNrowBins());
  BOOST_CHECK_EQUAL(staged.getNcolBins(), reference.getNcolBins());
  BOOST_CHECK_EQUAL(staged.getRowOrigin(), reference.getRowOrigin());
  BOOST_CHECK_EQUAL(staged.getRowCoordinateSpan(), reference.getRowCoordinateSpan());
  for (int iLayer = 0; iLayer < MFTN; ++iLayer) {
    BOOST_CHECK_EQUAL(staged.getLayerColHalfExtent(iLayer), reference.getLayerColHalfExtent(iLayer));
  }
}

/// Exact DiskDisk/XY parity with an explicit IndexRowMin/Max override.
BOOST_AUTO_TEST_CASE(DiskDiskBindingMatchesExplicitRowRangeOverride)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::MFT);
  params.IndexRowMin = -15.f;
  params.IndexRowMax = 15.f;

  IndexTableUtils<MFTN> staged;
  const auto error = bindIndexTableConfiguration<TransitionPolicyTag::DiskDisk, MFTN>(staged, params);
  BOOST_REQUIRE(error == IndexTableConfigError::None);

  IndexTableUtils<MFTN> reference;
  reference.setTrackingParametersXY(params, -15.f, 15.f);

  BOOST_CHECK_EQUAL(staged.getRowOrigin(), reference.getRowOrigin());
  BOOST_CHECK_EQUAL(staged.getRowCoordinateSpan(), reference.getRowCoordinateSpan());
  BOOST_CHECK_EQUAL(staged.getNrowBins(), reference.getNrowBins());
  BOOST_CHECK_EQUAL(staged.getNcolBins(), reference.getNcolBins());
}

/// params.NLayers may be a smaller active prefix of the template NLayers.
BOOST_AUTO_TEST_CASE(ActiveLayerCountSmallerThanTemplateIsAccepted)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::ITS);
  params.NLayers = 5;

  IndexTableUtils<ITSN> staged;
  const auto error = bindIndexTableConfiguration<TransitionPolicyTag::CylinderCylinder, ITSN>(staged, params);
  BOOST_CHECK(error == IndexTableConfigError::None);
  // Every template slot is still filled, regardless of the active prefix.
  for (int iLayer = 0; iLayer < ITSN; ++iLayer) {
    BOOST_CHECK_GT(staged.getLayerColHalfExtent(iLayer), 0.f);
  }
}

BOOST_AUTO_TEST_CASE(InvalidActiveLayerCountRejectedAndStagedUntouched)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::ITS);

  for (int n : {0, -1, ITSN + 1}) {
    params.NLayers = n;
    IndexTableUtils<ITSN> staged;
    const auto error = bindIndexTableConfiguration<TransitionPolicyTag::CylinderCylinder, ITSN>(staged, params);
    BOOST_CHECK(error == IndexTableConfigError::InvalidActiveLayerCount);
    BOOST_CHECK(isUntouched(staged));
  }
}

BOOST_AUTO_TEST_CASE(NonPositiveBinsRejectedAndStagedUntouched)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::ITS);

  params.RowBins = 0;
  {
    IndexTableUtils<ITSN> staged;
    BOOST_CHECK((bindIndexTableConfiguration<TransitionPolicyTag::CylinderCylinder, ITSN>(staged, params) == IndexTableConfigError::NonPositiveRowBins));
    BOOST_CHECK(isUntouched(staged));
  }
  params.RowBins = 128;
  params.ColBins = -5;
  {
    IndexTableUtils<ITSN> staged;
    BOOST_CHECK((bindIndexTableConfiguration<TransitionPolicyTag::CylinderCylinder, ITSN>(staged, params) == IndexTableConfigError::NonPositiveColBins));
    BOOST_CHECK(isUntouched(staged));
  }
}

/// Exact boundary: 46341*46341 = 2147488281 > INT_MAX; 46340*46340 =
/// 2147395600 <= INT_MAX. No false positive/negative at the boundary.
BOOST_AUTO_TEST_CASE(RowColBinCountExceedsIndexRangeBoundary)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::ITS);

  params.RowBins = 46341;
  params.ColBins = 46341;
  {
    IndexTableUtils<ITSN> staged;
    BOOST_CHECK((bindIndexTableConfiguration<TransitionPolicyTag::CylinderCylinder, ITSN>(staged, params) == IndexTableConfigError::RowColBinCountExceedsIndexRange));
    BOOST_CHECK(isUntouched(staged));
  }
  params.RowBins = 46340;
  params.ColBins = 46340;
  {
    IndexTableUtils<ITSN> staged;
    BOOST_CHECK((bindIndexTableConfiguration<TransitionPolicyTag::CylinderCylinder, ITSN>(staged, params) == IndexTableConfigError::None));
  }
}

/// Conditional source-size rule: when LayerColHalfExtent is non-empty it is
/// authoritative and must itself cover NLayers, even if LayerZ is long
/// enough; when it is empty, LayerZ must cover NLayers instead.
BOOST_AUTO_TEST_CASE(InsufficientExtentSourceRejectedAndStagedUntouched)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::ITS);
  BOOST_REQUIRE_EQUAL(params.LayerZ.size(), static_cast<size_t>(ITSN));

  // LayerColHalfExtent non-empty but short: rejected even though LayerZ is long enough.
  params.LayerColHalfExtent.assign(ITSN - 1, 1.f);
  {
    IndexTableUtils<ITSN> staged;
    BOOST_CHECK((bindIndexTableConfiguration<TransitionPolicyTag::CylinderCylinder, ITSN>(staged, params) == IndexTableConfigError::InsufficientLayerColHalfExtent));
    BOOST_CHECK(isUntouched(staged));
  }

  // LayerColHalfExtent empty, LayerZ short: rejected.
  params.LayerColHalfExtent.clear();
  params.LayerZ.resize(ITSN - 1);
  {
    IndexTableUtils<ITSN> staged;
    BOOST_CHECK((bindIndexTableConfiguration<TransitionPolicyTag::CylinderCylinder, ITSN>(staged, params) == IndexTableConfigError::InsufficientLayerColHalfExtent));
    BOOST_CHECK(isUntouched(staged));
  }
}

/// A naive `!(extent > 0)` predicate lets +Inf through (Inf > 0 is true).
/// NonFiniteColHalfExtent must catch it explicitly, separately from
/// NonPositiveColHalfExtent.
BOOST_AUTO_TEST_CASE(NonFiniteColHalfExtentRejectsNaNAndInfSeparatelyFromNonPositive)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::ITS);

  params.LayerZ[2] = std::numeric_limits<float>::infinity();
  {
    IndexTableUtils<ITSN> staged;
    BOOST_CHECK((bindIndexTableConfiguration<TransitionPolicyTag::CylinderCylinder, ITSN>(staged, params) == IndexTableConfigError::NonFiniteColHalfExtent));
    BOOST_CHECK(isUntouched(staged));
  }

  params.LayerZ[2] = std::numeric_limits<float>::quiet_NaN();
  {
    IndexTableUtils<ITSN> staged;
    BOOST_CHECK((bindIndexTableConfiguration<TransitionPolicyTag::CylinderCylinder, ITSN>(staged, params) == IndexTableConfigError::NonFiniteColHalfExtent));
    BOOST_CHECK(isUntouched(staged));
  }

  params.LayerZ[2] = 0.f;
  {
    IndexTableUtils<ITSN> staged;
    BOOST_CHECK((bindIndexTableConfiguration<TransitionPolicyTag::CylinderCylinder, ITSN>(staged, params) == IndexTableConfigError::NonPositiveColHalfExtent));
    BOOST_CHECK(isUntouched(staged));
  }

  params.LayerZ[2] = -3.f;
  {
    IndexTableUtils<ITSN> staged;
    BOOST_CHECK((bindIndexTableConfiguration<TransitionPolicyTag::CylinderCylinder, ITSN>(staged, params) == IndexTableConfigError::NonPositiveColHalfExtent));
    BOOST_CHECK(isUntouched(staged));
  }
}

BOOST_AUTO_TEST_CASE(RowRangeValidationIsXYOnly)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::MFT);

  params.IndexRowMin = std::numeric_limits<float>::quiet_NaN();
  params.IndexRowMax = 1.f;
  {
    IndexTableUtils<MFTN> staged;
    BOOST_CHECK((bindIndexTableConfiguration<TransitionPolicyTag::DiskDisk, MFTN>(staged, params) == IndexTableConfigError::NonFiniteRowRange));
    BOOST_CHECK(isUntouched(staged));
  }

  params.IndexRowMin = 10.f;
  params.IndexRowMax = 5.f; // rowMax <= rowMin
  {
    IndexTableUtils<MFTN> staged;
    BOOST_CHECK((bindIndexTableConfiguration<TransitionPolicyTag::DiskDisk, MFTN>(staged, params) == IndexTableConfigError::DegenerateRowRange));
    BOOST_CHECK(isUntouched(staged));
  }
}

BOOST_AUTO_TEST_CASE(ConfigurationsMatchIdentifiesEveryStoredField)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::MFT);

  IndexTableUtils<MFTN> a;
  BOOST_REQUIRE((bindIndexTableConfiguration<TransitionPolicyTag::DiskDisk, MFTN>(a, params) == IndexTableConfigError::None));

  IndexTableUtils<MFTN> b;
  BOOST_REQUIRE((bindIndexTableConfiguration<TransitionPolicyTag::DiskDisk, MFTN>(b, params) == IndexTableConfigError::None));
  BOOST_CHECK(indexTableConfigurationsMatch(a, b));

  // Different RowBins/ColBins.
  auto diffBins = params;
  diffBins.RowBins = params.RowBins + 1;
  IndexTableUtils<MFTN> c;
  BOOST_REQUIRE((bindIndexTableConfiguration<TransitionPolicyTag::DiskDisk, MFTN>(c, diffBins) == IndexTableConfigError::None));
  BOOST_CHECK(!indexTableConfigurationsMatch(a, c));

  // Different row range.
  auto diffRange = params;
  diffRange.IndexRowMin = -15.f;
  diffRange.IndexRowMax = 15.f;
  IndexTableUtils<MFTN> d;
  BOOST_REQUIRE((bindIndexTableConfiguration<TransitionPolicyTag::DiskDisk, MFTN>(d, diffRange) == IndexTableConfigError::None));
  BOOST_CHECK(!indexTableConfigurationsMatch(a, d));

  // Different per-layer extent.
  auto diffExtent = params;
  diffExtent.LayerColHalfExtent[0] += 1.f;
  IndexTableUtils<MFTN> e;
  BOOST_REQUIRE((bindIndexTableConfiguration<TransitionPolicyTag::DiskDisk, MFTN>(e, diffExtent) == IndexTableConfigError::None));
  BOOST_CHECK(!indexTableConfigurationsMatch(a, e));
}

BOOST_AUTO_TEST_CASE(CheckedIndexTableSizeProductBoundaries)
{
  std::size_t result = 12345;

  BOOST_CHECK(checkedIndexTableSizeProduct(0, 999999, result));
  BOOST_CHECK_EQUAL(result, 0u);

  BOOST_CHECK(checkedIndexTableSizeProduct(3, 4, result));
  BOOST_CHECK_EQUAL(result, 12u);

  BOOST_CHECK(checkedIndexTableSizeProduct(std::numeric_limits<std::size_t>::max(), 1, result));
  BOOST_CHECK_EQUAL(result, std::numeric_limits<std::size_t>::max());

  BOOST_CHECK(!checkedIndexTableSizeProduct(std::numeric_limits<std::size_t>::max(), 2, result));
  BOOST_CHECK(!checkedIndexTableSizeProduct(std::numeric_limits<std::size_t>::max() / 2 + 2, 2, result));
}

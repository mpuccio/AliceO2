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

bool isUntouched(const IndexTableUtilsCore& utils)
{
  // Default-constructed sentinel state (IndexTableUtils.h): never mutated.
  return utils.getNrowBins() == 0 && utils.getNcolBins() == 0 && utils.getCoordType() == IndexTableCoordType::PhiZ;
}

IndexTableConfigError bindIndexTableConfiguration(IndexTableUtilsCore& staged,
                                                  const TrackingParameters& params,
                                                  int activeSurfaceCount,
                                                  SurfaceKind kind)
{
  std::array<SurfaceChartRange, IndexTableUtilsCore::MaxLayers> ranges{};
  const auto& extents = params.LayerColHalfExtent.empty() ? params.LayerZ : params.LayerColHalfExtent;
  const int rangeCount = std::min(activeSurfaceCount, static_cast<int>(extents.size()));
  for (int iLayer = 0; iLayer < rangeCount; ++iLayer) {
    ranges[iLayer] = kind == SurfaceKind::Disk ? SurfaceChartRange{0.1f, 20.f} : SurfaceChartRange{-extents[iLayer], extents[iLayer]};
  }
  return o2::itsmft::tracking::bindIndexTableConfiguration(staged, params, activeSurfaceCount, kind,
                                                           gsl::span<const SurfaceChartRange>{ranges.data(), static_cast<size_t>(rangeCount)});
}
} // namespace

/// Exact Cylinder/PhiZ parity: bindIndexTableConfiguration must
/// reproduce IndexTableUtils::setTrackingParameters() bit-for-bit for the
/// default ITS-shaped TrackingParameters.
BOOST_AUTO_TEST_CASE(CylinderBindingMatchesSetTrackingParameters)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::ITS);

  IndexTableUtilsCore staged;
  const auto error = bindIndexTableConfiguration(staged, params, ITSN, SurfaceKind::Cylinder);
  BOOST_REQUIRE(error == IndexTableConfigError::None);

  IndexTableUtilsCore reference;
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

/// Disk charts are periodic phi with descriptor-owned radial intervals.
BOOST_AUTO_TEST_CASE(DiskBindingMatchesDefaultRowRangeFallback)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::MFT);
  params.IndexRowMin = 0.f;
  params.IndexRowMax = 0.f; // force the "no explicit range" fallback path

  IndexTableUtilsCore staged;
  const auto error = bindIndexTableConfiguration(staged, params, MFTN, SurfaceKind::Disk);
  BOOST_REQUIRE(error == IndexTableConfigError::None);

  BOOST_CHECK(staged.getCoordType() == IndexTableCoordType::PhiR);
  BOOST_CHECK_EQUAL(staged.getNrowBins(), params.RowBins);
  BOOST_CHECK_EQUAL(staged.getNcolBins(), params.ColBins);
  BOOST_CHECK_EQUAL(staged.getRowOrigin(), 0.f);
  BOOST_CHECK_EQUAL(staged.getRowCoordinateSpan(), o2::constants::math::TwoPI);
  for (int iLayer = 0; iLayer < MFTN; ++iLayer) {
    BOOST_CHECK_EQUAL(staged.getLayerColMin(iLayer), 0.1f);
    BOOST_CHECK_EQUAL(staged.getLayerColMax(iLayer), 20.f);
  }
}

/// Obsolete Cartesian row-range knobs cannot alter a periodic disk chart.
BOOST_AUTO_TEST_CASE(DiskBindingMatchesExplicitRowRangeOverride)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::MFT);
  params.IndexRowMin = -15.f;
  params.IndexRowMax = 15.f;

  IndexTableUtilsCore staged;
  const auto error = bindIndexTableConfiguration(staged, params, MFTN, SurfaceKind::Disk);
  BOOST_REQUIRE(error == IndexTableConfigError::None);

  BOOST_CHECK_EQUAL(staged.getRowOrigin(), 0.f);
  BOOST_CHECK_EQUAL(staged.getRowCoordinateSpan(), o2::constants::math::TwoPI);
  BOOST_CHECK_EQUAL(staged.getNrowBins(), params.RowBins);
  BOOST_CHECK_EQUAL(staged.getNcolBins(), params.ColBins);
}

/// params.NLayers may be a smaller active prefix of the template NLayers.
BOOST_AUTO_TEST_CASE(ActiveLayerCountSmallerThanTemplateIsAccepted)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::ITS);
  params.NLayers = 5;

  IndexTableUtilsCore staged;
  const auto error = bindIndexTableConfiguration(staged, params, 5, SurfaceKind::Cylinder);
  BOOST_CHECK(error == IndexTableConfigError::None);
  // The fixed device arrays retain their capacity, but only the runtime
  // active prefix is populated by the binder.
  for (int iLayer = 0; iLayer < 5; ++iLayer) {
    BOOST_CHECK_GT(staged.getLayerColHalfExtent(iLayer), 0.f);
  }
}

BOOST_AUTO_TEST_CASE(InvalidActiveLayerCountRejectedAndStagedUntouched)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::ITS);

  for (int n : {0, -1, ITSN + 1}) {
    params.NLayers = n;
    IndexTableUtilsCore staged;
    const auto error = bindIndexTableConfiguration(staged, params, ITSN, SurfaceKind::Cylinder);
    BOOST_CHECK(error == IndexTableConfigError::InvalidActiveLayerCount);
    BOOST_CHECK(isUntouched(staged));
  }
}

BOOST_AUTO_TEST_CASE(InvalidSurfaceKindRejectedAndStagedUntouched)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::ITS);
  IndexTableUtilsCore staged;
  const auto invalidKind = static_cast<SurfaceKind>(99);
  BOOST_CHECK(bindIndexTableConfiguration(staged, params, ITSN, invalidKind) == IndexTableConfigError::InvalidSurfaceKind);
  BOOST_CHECK(isUntouched(staged));
}

BOOST_AUTO_TEST_CASE(NonPositiveBinsRejectedAndStagedUntouched)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::ITS);

  params.RowBins = 0;
  {
    IndexTableUtilsCore staged;
    BOOST_CHECK((bindIndexTableConfiguration(staged, params, ITSN, SurfaceKind::Cylinder) == IndexTableConfigError::NonPositiveRowBins));
    BOOST_CHECK(isUntouched(staged));
  }
  params.RowBins = 128;
  params.ColBins = -5;
  {
    IndexTableUtilsCore staged;
    BOOST_CHECK((bindIndexTableConfiguration(staged, params, ITSN, SurfaceKind::Cylinder) == IndexTableConfigError::NonPositiveColBins));
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
    IndexTableUtilsCore staged;
    BOOST_CHECK((bindIndexTableConfiguration(staged, params, ITSN, SurfaceKind::Cylinder) == IndexTableConfigError::RowColBinCountExceedsIndexRange));
    BOOST_CHECK(isUntouched(staged));
  }
  params.RowBins = 46340;
  params.ColBins = 46340;
  {
    IndexTableUtilsCore staged;
    BOOST_CHECK((bindIndexTableConfiguration(staged, params, ITSN, SurfaceKind::Cylinder) == IndexTableConfigError::None));
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
    IndexTableUtilsCore staged;
    BOOST_CHECK((bindIndexTableConfiguration(staged, params, ITSN, SurfaceKind::Cylinder) == IndexTableConfigError::InsufficientChartRanges));
    BOOST_CHECK(isUntouched(staged));
  }

  // LayerColHalfExtent empty, LayerZ short: rejected.
  params.LayerColHalfExtent.clear();
  params.LayerZ.resize(ITSN - 1);
  {
    IndexTableUtilsCore staged;
    BOOST_CHECK((bindIndexTableConfiguration(staged, params, ITSN, SurfaceKind::Cylinder) == IndexTableConfigError::InsufficientChartRanges));
    BOOST_CHECK(isUntouched(staged));
  }
}

/// A naive `!(extent > 0)` predicate lets +Inf through (Inf > 0 is true).
/// NonFiniteChartRange must catch it explicitly, separately from
/// InvalidChartRange.
BOOST_AUTO_TEST_CASE(NonFiniteChartRangeRejectsNaNAndInfSeparatelyFromInvalidRange)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::ITS);

  params.LayerZ[2] = std::numeric_limits<float>::infinity();
  {
    IndexTableUtilsCore staged;
    BOOST_CHECK((bindIndexTableConfiguration(staged, params, ITSN, SurfaceKind::Cylinder) == IndexTableConfigError::NonFiniteChartRange));
    BOOST_CHECK(isUntouched(staged));
  }

  params.LayerZ[2] = std::numeric_limits<float>::quiet_NaN();
  {
    IndexTableUtilsCore staged;
    BOOST_CHECK((bindIndexTableConfiguration(staged, params, ITSN, SurfaceKind::Cylinder) == IndexTableConfigError::NonFiniteChartRange));
    BOOST_CHECK(isUntouched(staged));
  }

  params.LayerZ[2] = 0.f;
  {
    IndexTableUtilsCore staged;
    BOOST_CHECK((bindIndexTableConfiguration(staged, params, ITSN, SurfaceKind::Cylinder) == IndexTableConfigError::InvalidChartRange));
    BOOST_CHECK(isUntouched(staged));
  }

  params.LayerZ[2] = -3.f;
  {
    IndexTableUtilsCore staged;
    BOOST_CHECK((bindIndexTableConfiguration(staged, params, ITSN, SurfaceKind::Cylinder) == IndexTableConfigError::InvalidChartRange));
    BOOST_CHECK(isUntouched(staged));
  }
}

BOOST_AUTO_TEST_CASE(ChartRangeValidationRejectsNonFiniteAndDegenerateIntervals)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::MFT);

  params.IndexRowMin = std::numeric_limits<float>::quiet_NaN();
  params.IndexRowMax = 1.f;
  {
    IndexTableUtilsCore staged;
    std::array<SurfaceChartRange, MFTN> ranges{};
    ranges.fill({0.1f, 20.f});
    ranges[0].min = std::numeric_limits<float>::quiet_NaN();
    BOOST_CHECK((o2::itsmft::tracking::bindIndexTableConfiguration(staged, params, MFTN, SurfaceKind::Disk, ranges) == IndexTableConfigError::NonFiniteChartRange));
    BOOST_CHECK(isUntouched(staged));
  }

  params.IndexRowMin = 10.f;
  params.IndexRowMax = 5.f; // rowMax <= rowMin
  {
    IndexTableUtilsCore staged;
    std::array<SurfaceChartRange, MFTN> ranges{};
    ranges.fill({0.1f, 20.f});
    ranges[0] = {10.f, 5.f};
    BOOST_CHECK((o2::itsmft::tracking::bindIndexTableConfiguration(staged, params, MFTN, SurfaceKind::Disk, ranges) == IndexTableConfigError::InvalidChartRange));
    BOOST_CHECK(isUntouched(staged));
  }
}

BOOST_AUTO_TEST_CASE(ConfigurationsMatchIdentifiesEveryStoredField)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::MFT);

  IndexTableUtilsCore a;
  BOOST_REQUIRE((bindIndexTableConfiguration(a, params, MFTN, SurfaceKind::Disk) == IndexTableConfigError::None));

  IndexTableUtilsCore b;
  BOOST_REQUIRE((bindIndexTableConfiguration(b, params, MFTN, SurfaceKind::Disk) == IndexTableConfigError::None));
  BOOST_CHECK(indexTableConfigurationsMatch(a, b, MFTN));

  // Different RowBins/ColBins.
  auto diffBins = params;
  diffBins.RowBins = params.RowBins + 1;
  IndexTableUtilsCore c;
  BOOST_REQUIRE((bindIndexTableConfiguration(c, diffBins, MFTN, SurfaceKind::Disk) == IndexTableConfigError::None));
  BOOST_CHECK(!indexTableConfigurationsMatch(a, c, MFTN));

  // Legacy row-range knobs do not alter a periodic chart.
  auto diffRange = params;
  diffRange.IndexRowMin = -15.f;
  diffRange.IndexRowMax = 15.f;
  IndexTableUtilsCore d;
  BOOST_REQUIRE((bindIndexTableConfiguration(d, diffRange, MFTN, SurfaceKind::Disk) == IndexTableConfigError::None));
  BOOST_CHECK(indexTableConfigurationsMatch(a, d, MFTN));

  // Obsolete Cartesian extents no longer affect the chart contract.
  auto diffExtent = params;
  diffExtent.LayerColHalfExtent[0] += 1.f;
  IndexTableUtilsCore e;
  BOOST_REQUIRE((bindIndexTableConfiguration(e, diffExtent, MFTN, SurfaceKind::Disk) == IndexTableConfigError::None));
  BOOST_CHECK(indexTableConfigurationsMatch(a, e, MFTN));
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

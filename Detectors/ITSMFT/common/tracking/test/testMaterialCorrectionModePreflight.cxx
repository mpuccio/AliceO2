// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Stage-B Slice C: focused tests for the pure material-correction-mode
// preflight (materialCorrectionModeSupport),
// and the new TraversalFailureReason::UnsupportedMaterialCorrectionMode value
// (TrackerTraits.h). The preflight is wired into
// TrackerTraits::initialiseTimeFrame(); these tests cover its direct result
// and the retained validation contract.

#define BOOST_TEST_MODULE ITSMFT MaterialCorrectionModePreflight
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <type_traits>

#include <gsl/span>

#include "DetectorsBase/Propagator.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/Configuration.h"

using namespace o2::itsmft::tracking;
using o2::itsmft::AttachHitConfigView;
using o2::itsmft::isRecognizedMatCorrType;
using MatCorrType = o2::base::PropagatorF::MatCorrType;

// --- Exact TraversalFailureReason numeric value ----------------------------
//
// Appended without renumbering any existing value (TrackerTraits.h): the
// next value after LegacyMaterialMismatch (10).

static_assert(static_cast<uint8_t>(TraversalFailureReason::MissingLayout) == 0);
static_assert(static_cast<uint8_t>(TraversalFailureReason::StaleLayout) == 1);
static_assert(static_cast<uint8_t>(TraversalFailureReason::IterationOutOfRange) == 2);
static_assert(static_cast<uint8_t>(TraversalFailureReason::SparseTopologyMismatch) == 3);
static_assert(static_cast<uint8_t>(TraversalFailureReason::InvalidTraversalSchedule) == 4);
static_assert(static_cast<uint8_t>(TraversalFailureReason::MixedSurfaceKindLayout) == 5);
static_assert(static_cast<uint8_t>(TraversalFailureReason::SurfaceKindMismatch) == 6);
static_assert(static_cast<uint8_t>(TraversalFailureReason::InvalidSurfaceParameters) == 7);
static_assert(static_cast<uint8_t>(TraversalFailureReason::InvalidIndexTableConfiguration) == 8);
static_assert(static_cast<uint8_t>(TraversalFailureReason::IndexTableConfigurationMismatch) == 9);
static_assert(static_cast<uint8_t>(TraversalFailureReason::LegacyMaterialMismatch) == 10);
static_assert(static_cast<uint8_t>(TraversalFailureReason::UnsupportedMaterialCorrectionMode) == 11);

BOOST_AUTO_TEST_CASE(UnsupportedMaterialCorrectionModeHasExactValue)
{
  BOOST_CHECK_EQUAL(static_cast<int>(TraversalFailureReason::UnsupportedMaterialCorrectionMode), 11);
}

// --- Cylinder: every recognized CorrType --------------------------

BOOST_AUTO_TEST_CASE(CylinderNoneIsSupported)
{
  BOOST_CHECK(materialCorrectionModeSupport(SurfaceKind::Cylinder, MatCorrType::USEMatCorrNONE) ==
              MaterialCorrectionModeSupport::Supported);
}

BOOST_AUTO_TEST_CASE(CylinderLutIsRecognizedButUnsupported)
{
  BOOST_CHECK(materialCorrectionModeSupport(SurfaceKind::Cylinder, MatCorrType::USEMatCorrLUT) ==
              MaterialCorrectionModeSupport::Unsupported);
}

BOOST_AUTO_TEST_CASE(CylinderTGeoIsRecognizedButUnsupported)
{
  BOOST_CHECK(materialCorrectionModeSupport(SurfaceKind::Cylinder, MatCorrType::USEMatCorrTGeo) ==
              MaterialCorrectionModeSupport::Unsupported);
}

BOOST_AUTO_TEST_CASE(CylinderInvalidCorrTypeIsInvalidModeNotUnsupported)
{
  const auto invalid = static_cast<MatCorrType>(99);
  const auto result = materialCorrectionModeSupport(SurfaceKind::Cylinder, invalid);
  BOOST_CHECK(result == MaterialCorrectionModeSupport::InvalidMode);
  BOOST_CHECK(result != MaterialCorrectionModeSupport::Unsupported);
  BOOST_CHECK(result != MaterialCorrectionModeSupport::Supported);
}

// --- Disk: not constrained by this preflight ---------------------------

BOOST_AUTO_TEST_CASE(DiskNoneIsSupported)
{
  BOOST_CHECK(materialCorrectionModeSupport(SurfaceKind::Disk, MatCorrType::USEMatCorrNONE) ==
              MaterialCorrectionModeSupport::Supported);
}

BOOST_AUTO_TEST_CASE(DiskLutIsSupported)
{
  BOOST_CHECK(materialCorrectionModeSupport(SurfaceKind::Disk, MatCorrType::USEMatCorrLUT) ==
              MaterialCorrectionModeSupport::Supported);
}

BOOST_AUTO_TEST_CASE(DiskTGeoIsSupported)
{
  BOOST_CHECK(materialCorrectionModeSupport(SurfaceKind::Disk, MatCorrType::USEMatCorrTGeo) ==
              MaterialCorrectionModeSupport::Supported);
}

BOOST_AUTO_TEST_CASE(DiskInvalidCorrTypeIsInvalidModeNotSilentlySupported)
{
  const auto invalid = static_cast<MatCorrType>(99);
  const auto result = materialCorrectionModeSupport(SurfaceKind::Disk, invalid);
  BOOST_CHECK(result == MaterialCorrectionModeSupport::InvalidMode);
  BOOST_CHECK(result != MaterialCorrectionModeSupport::Supported);
}

BOOST_AUTO_TEST_CASE(InvalidSurfaceKindIsDeferredAndCorrTypeStillHasPrecedence)
{
  const auto invalidKind = static_cast<SurfaceKind>(99);
  BOOST_CHECK(materialCorrectionModeSupport(invalidKind, MatCorrType::USEMatCorrNONE) ==
              MaterialCorrectionModeSupport::InvalidSurfaceKind);
  BOOST_CHECK(materialCorrectionModeSupport(invalidKind, static_cast<MatCorrType>(99)) ==
              MaterialCorrectionModeSupport::InvalidMode);
}

// --- Exact distinction between Unsupported and InvalidMode ------------------

BOOST_AUTO_TEST_CASE(UnsupportedAndInvalidModeAreDistinctValues)
{
  BOOST_CHECK(MaterialCorrectionModeSupport::Supported != MaterialCorrectionModeSupport::Unsupported);
  BOOST_CHECK(MaterialCorrectionModeSupport::Supported != MaterialCorrectionModeSupport::InvalidMode);
  BOOST_CHECK(MaterialCorrectionModeSupport::Unsupported != MaterialCorrectionModeSupport::InvalidMode);

  // The recognized-but-unsupported Cylinder cases must never collapse
  // onto the invalid-mode result, and vice versa.
  BOOST_CHECK(materialCorrectionModeSupport(SurfaceKind::Cylinder, MatCorrType::USEMatCorrLUT) ==
              MaterialCorrectionModeSupport::Unsupported);
  BOOST_CHECK(materialCorrectionModeSupport(SurfaceKind::Cylinder, static_cast<MatCorrType>(-1)) ==
              MaterialCorrectionModeSupport::InvalidMode);
}

// --- Determinism / repeated calls ------------------------------------------

BOOST_AUTO_TEST_CASE(RepeatedCallsAreDeterministic)
{
  for (int i = 0; i < 5; ++i) {
    BOOST_CHECK(materialCorrectionModeSupport(SurfaceKind::Cylinder, MatCorrType::USEMatCorrNONE) ==
                MaterialCorrectionModeSupport::Supported);
    BOOST_CHECK(materialCorrectionModeSupport(SurfaceKind::Cylinder, MatCorrType::USEMatCorrLUT) ==
                MaterialCorrectionModeSupport::Unsupported);
    BOOST_CHECK(materialCorrectionModeSupport(SurfaceKind::Disk, MatCorrType::USEMatCorrTGeo) ==
                MaterialCorrectionModeSupport::Supported);
    BOOST_CHECK(materialCorrectionModeSupport(SurfaceKind::Cylinder, static_cast<MatCorrType>(42)) ==
                MaterialCorrectionModeSupport::InvalidMode);
  }
}

// --- noexcept / signature shape ---------------------------------------------

static_assert(noexcept(materialCorrectionModeSupport(SurfaceKind::Cylinder, MatCorrType::USEMatCorrNONE)));
static_assert(noexcept(materialCorrectionModeSupport(SurfaceKind::Disk, MatCorrType::USEMatCorrNONE)));
static_assert(std::is_trivially_copyable_v<MaterialCorrectionModeSupport>);

// --- isRecognizedMatCorrType is the single shared definition ----------------
//
// AttachHitConfigView::isValid() must still accept all three recognized
// modes and reject an invalid one, using exactly the same predicate this
// preflight uses -- proving the two contracts cannot diverge.

BOOST_AUTO_TEST_CASE(AttachHitConfigViewAcceptsAllRecognizedModes)
{
  const std::array<NominalSurfaceMaterial, 1> material{NominalSurfaceMaterial{0.01f, 0.02f}};
  for (const auto corrType : {MatCorrType::USEMatCorrNONE, MatCorrType::USEMatCorrLUT, MatCorrType::USEMatCorrTGeo}) {
    AttachHitConfigView view{gsl::span<const NominalSurfaceMaterial>(material), corrType};
    BOOST_CHECK(view.isValid(1));
  }
}

BOOST_AUTO_TEST_CASE(AttachHitConfigViewRejectsInvalidCorrType)
{
  const std::array<NominalSurfaceMaterial, 1> material{NominalSurfaceMaterial{0.01f, 0.02f}};
  AttachHitConfigView view{gsl::span<const NominalSurfaceMaterial>(material), static_cast<MatCorrType>(99)};
  BOOST_CHECK(!view.isValid(1));
}

BOOST_AUTO_TEST_CASE(IsRecognizedMatCorrTypeMatchesAttachHitConfigViewContract)
{
  BOOST_CHECK(isRecognizedMatCorrType(MatCorrType::USEMatCorrNONE));
  BOOST_CHECK(isRecognizedMatCorrType(MatCorrType::USEMatCorrLUT));
  BOOST_CHECK(isRecognizedMatCorrType(MatCorrType::USEMatCorrTGeo));
  BOOST_CHECK(!isRecognizedMatCorrType(static_cast<MatCorrType>(99)));
  BOOST_CHECK(!isRecognizedMatCorrType(static_cast<MatCorrType>(-1)));
}

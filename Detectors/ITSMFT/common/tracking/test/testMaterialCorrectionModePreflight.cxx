// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Stage-B Slice C: focused tests for the pure material-correction-mode
// preflight (checkMaterialCorrectionModeSupport<Tag>, TransitionPolicyBinding.h)
// and the new TraversalFailureReason::UnsupportedMaterialCorrectionMode value
// (TrackerTraits.h). This slice is additive and unwired: no production call
// site invokes the preflight, and none of these tests exercise
// TrackerTraits::initialiseTimeFrame() or any TimeFrame/candidate state.

#define BOOST_TEST_MODULE ITSMFT MaterialCorrectionModePreflight
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <type_traits>

#include <gsl/span>

#include "DetectorsBase/Propagator.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/TransitionPolicyBinding.h"

using namespace o2::itsmft::tracking;
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
static_assert(static_cast<uint8_t>(TraversalFailureReason::MixedPolicyLayout) == 5);
static_assert(static_cast<uint8_t>(TraversalFailureReason::StateFamilyMismatch) == 6);
static_assert(static_cast<uint8_t>(TraversalFailureReason::InvalidPolicyParameters) == 7);
static_assert(static_cast<uint8_t>(TraversalFailureReason::InvalidIndexTableConfiguration) == 8);
static_assert(static_cast<uint8_t>(TraversalFailureReason::IndexTableConfigurationMismatch) == 9);
static_assert(static_cast<uint8_t>(TraversalFailureReason::LegacyMaterialMismatch) == 10);
static_assert(static_cast<uint8_t>(TraversalFailureReason::UnsupportedMaterialCorrectionMode) == 11);

BOOST_AUTO_TEST_CASE(UnsupportedMaterialCorrectionModeHasExactValue)
{
  BOOST_CHECK_EQUAL(static_cast<int>(TraversalFailureReason::UnsupportedMaterialCorrectionMode), 11);
}

// --- Compile-time rejection for unsupported policy tags --------------------
//
// checkMaterialCorrectionModeSupport<Tag>'s primary template is `= delete`d
// (TransitionPolicyBinding.h): instantiating it for TransitionPolicyTag::Invalid
// is a hard, SFINAE-observable compile-time error, never a silent fallback to
// Supported/Unsupported. Proven with a dependent-trailing-decltype wrapper
// (the standard idiom std::is_invocable itself relies on to detect deleted
// candidates), not merely asserted in a comment.

namespace
{
template <TransitionPolicyTag Tag>
struct MaterialModeInvoker {
  template <class CorrTypeArg>
  auto operator()(CorrTypeArg corrType) const noexcept -> decltype(checkMaterialCorrectionModeSupport<Tag>(corrType))
  {
    return checkMaterialCorrectionModeSupport<Tag>(corrType);
  }
};
} // namespace

static_assert(std::is_invocable_v<MaterialModeInvoker<TransitionPolicyTag::CylinderCylinder>, MatCorrType>,
              "CylinderCylinder must have a usable specialization");
static_assert(std::is_invocable_v<MaterialModeInvoker<TransitionPolicyTag::DiskDisk>, MatCorrType>,
              "DiskDisk must have a usable specialization");
static_assert(!std::is_invocable_v<MaterialModeInvoker<TransitionPolicyTag::Invalid>, MatCorrType>,
              "Invalid must be rejected at compile time, not silently fall back to a supported family");

// --- CylinderCylinder: every recognized CorrType --------------------------

BOOST_AUTO_TEST_CASE(CylinderCylinderNoneIsSupported)
{
  BOOST_CHECK(checkMaterialCorrectionModeSupport<TransitionPolicyTag::CylinderCylinder>(MatCorrType::USEMatCorrNONE) ==
              MaterialCorrectionModeSupport::Supported);
}

BOOST_AUTO_TEST_CASE(CylinderCylinderLutIsRecognizedButUnsupported)
{
  BOOST_CHECK(checkMaterialCorrectionModeSupport<TransitionPolicyTag::CylinderCylinder>(MatCorrType::USEMatCorrLUT) ==
              MaterialCorrectionModeSupport::Unsupported);
}

BOOST_AUTO_TEST_CASE(CylinderCylinderTGeoIsRecognizedButUnsupported)
{
  BOOST_CHECK(checkMaterialCorrectionModeSupport<TransitionPolicyTag::CylinderCylinder>(MatCorrType::USEMatCorrTGeo) ==
              MaterialCorrectionModeSupport::Unsupported);
}

BOOST_AUTO_TEST_CASE(CylinderCylinderInvalidCorrTypeIsInvalidModeNotUnsupported)
{
  const auto invalid = static_cast<MatCorrType>(99);
  const auto result = checkMaterialCorrectionModeSupport<TransitionPolicyTag::CylinderCylinder>(invalid);
  BOOST_CHECK(result == MaterialCorrectionModeSupport::InvalidMode);
  BOOST_CHECK(result != MaterialCorrectionModeSupport::Unsupported);
  BOOST_CHECK(result != MaterialCorrectionModeSupport::Supported);
}

// --- DiskDisk: not constrained by this preflight ---------------------------

BOOST_AUTO_TEST_CASE(DiskDiskNoneIsSupported)
{
  BOOST_CHECK(checkMaterialCorrectionModeSupport<TransitionPolicyTag::DiskDisk>(MatCorrType::USEMatCorrNONE) ==
              MaterialCorrectionModeSupport::Supported);
}

BOOST_AUTO_TEST_CASE(DiskDiskLutIsSupported)
{
  BOOST_CHECK(checkMaterialCorrectionModeSupport<TransitionPolicyTag::DiskDisk>(MatCorrType::USEMatCorrLUT) ==
              MaterialCorrectionModeSupport::Supported);
}

BOOST_AUTO_TEST_CASE(DiskDiskTGeoIsSupported)
{
  BOOST_CHECK(checkMaterialCorrectionModeSupport<TransitionPolicyTag::DiskDisk>(MatCorrType::USEMatCorrTGeo) ==
              MaterialCorrectionModeSupport::Supported);
}

BOOST_AUTO_TEST_CASE(DiskDiskInvalidCorrTypeIsInvalidModeNotSilentlySupported)
{
  const auto invalid = static_cast<MatCorrType>(99);
  const auto result = checkMaterialCorrectionModeSupport<TransitionPolicyTag::DiskDisk>(invalid);
  BOOST_CHECK(result == MaterialCorrectionModeSupport::InvalidMode);
  BOOST_CHECK(result != MaterialCorrectionModeSupport::Supported);
}

// --- Exact distinction between Unsupported and InvalidMode ------------------

BOOST_AUTO_TEST_CASE(UnsupportedAndInvalidModeAreDistinctValues)
{
  BOOST_CHECK(MaterialCorrectionModeSupport::Supported != MaterialCorrectionModeSupport::Unsupported);
  BOOST_CHECK(MaterialCorrectionModeSupport::Supported != MaterialCorrectionModeSupport::InvalidMode);
  BOOST_CHECK(MaterialCorrectionModeSupport::Unsupported != MaterialCorrectionModeSupport::InvalidMode);

  // The recognized-but-unsupported CylinderCylinder cases must never collapse
  // onto the invalid-mode result, and vice versa.
  BOOST_CHECK(checkMaterialCorrectionModeSupport<TransitionPolicyTag::CylinderCylinder>(MatCorrType::USEMatCorrLUT) ==
              MaterialCorrectionModeSupport::Unsupported);
  BOOST_CHECK(checkMaterialCorrectionModeSupport<TransitionPolicyTag::CylinderCylinder>(static_cast<MatCorrType>(-1)) ==
              MaterialCorrectionModeSupport::InvalidMode);
}

// --- Determinism / repeated calls ------------------------------------------

BOOST_AUTO_TEST_CASE(RepeatedCallsAreDeterministic)
{
  for (int i = 0; i < 5; ++i) {
    BOOST_CHECK(checkMaterialCorrectionModeSupport<TransitionPolicyTag::CylinderCylinder>(MatCorrType::USEMatCorrNONE) ==
                MaterialCorrectionModeSupport::Supported);
    BOOST_CHECK(checkMaterialCorrectionModeSupport<TransitionPolicyTag::CylinderCylinder>(MatCorrType::USEMatCorrLUT) ==
                MaterialCorrectionModeSupport::Unsupported);
    BOOST_CHECK(checkMaterialCorrectionModeSupport<TransitionPolicyTag::DiskDisk>(MatCorrType::USEMatCorrTGeo) ==
                MaterialCorrectionModeSupport::Supported);
    BOOST_CHECK(checkMaterialCorrectionModeSupport<TransitionPolicyTag::CylinderCylinder>(static_cast<MatCorrType>(42)) ==
                MaterialCorrectionModeSupport::InvalidMode);
  }
}

// --- noexcept / signature shape ---------------------------------------------

static_assert(noexcept(checkMaterialCorrectionModeSupport<TransitionPolicyTag::CylinderCylinder>(MatCorrType::USEMatCorrNONE)));
static_assert(noexcept(checkMaterialCorrectionModeSupport<TransitionPolicyTag::DiskDisk>(MatCorrType::USEMatCorrNONE)));
static_assert(std::is_trivially_copyable_v<MaterialCorrectionModeSupport>);

// --- isRecognizedMatCorrType is the single shared definition ----------------
//
// AttachHitPolicyConfigView::isValid() must still accept all three recognized
// modes and reject an invalid one, using exactly the same predicate this
// preflight uses -- proving the two contracts cannot diverge.

BOOST_AUTO_TEST_CASE(AttachHitPolicyConfigViewAcceptsAllRecognizedModes)
{
  const std::array<NominalSurfaceMaterial, 1> material{NominalSurfaceMaterial{0.01f, 0.02f}};
  for (const auto corrType : {MatCorrType::USEMatCorrNONE, MatCorrType::USEMatCorrLUT, MatCorrType::USEMatCorrTGeo}) {
    AttachHitPolicyConfigView view{gsl::span<const NominalSurfaceMaterial>(material), corrType};
    BOOST_CHECK(view.isValid(1));
  }
}

BOOST_AUTO_TEST_CASE(AttachHitPolicyConfigViewRejectsInvalidCorrType)
{
  const std::array<NominalSurfaceMaterial, 1> material{NominalSurfaceMaterial{0.01f, 0.02f}};
  AttachHitPolicyConfigView view{gsl::span<const NominalSurfaceMaterial>(material), static_cast<MatCorrType>(99)};
  BOOST_CHECK(!view.isValid(1));
}

BOOST_AUTO_TEST_CASE(IsRecognizedMatCorrTypeMatchesAttachHitPolicyConfigViewContract)
{
  BOOST_CHECK(isRecognizedMatCorrType(MatCorrType::USEMatCorrNONE));
  BOOST_CHECK(isRecognizedMatCorrType(MatCorrType::USEMatCorrLUT));
  BOOST_CHECK(isRecognizedMatCorrType(MatCorrType::USEMatCorrTGeo));
  BOOST_CHECK(!isRecognizedMatCorrType(static_cast<MatCorrType>(99)));
  BOOST_CHECK(!isRecognizedMatCorrType(static_cast<MatCorrType>(-1)));
}

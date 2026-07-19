// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT TransitionPolicyState
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <cstddef>
#include <limits>
#include <type_traits>

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/TransitionPolicyState.h"

using namespace o2::itsmft::tracking;

/// Proves the cylinder-cylinder and disk-disk policy/state parameter
/// boundaries introduced for Gate 2 are device-compatible PODs with a
/// field-level ABI, and that their tag/family/surface-kind identity agrees
/// with the Gate 1 TransitionPolicyTag foundation rather than duplicating or
/// contradicting it (D007 / Architecture.md §10.1).

BOOST_AUTO_TEST_CASE(PolicyParamsAreDeviceCompatiblePods)
{
  BOOST_CHECK(std::is_standard_layout_v<CylinderCylinderPolicyParams>);
  BOOST_CHECK(std::is_trivially_copyable_v<CylinderCylinderPolicyParams>);
  BOOST_CHECK(std::is_standard_layout_v<DiskDiskPolicyParams>);
  BOOST_CHECK(std::is_trivially_copyable_v<DiskDiskPolicyParams>);
}

BOOST_AUTO_TEST_CASE(PolicyParamDefaultsMatchCurrentProductionValues)
{
  // Mirrors TrackingParameters defaults in Configuration.h; this is a
  // parity check on the new boundary's defaults, not a claim that
  // production tracking reads from these structs yet.
  CylinderCylinderPolicyParams barrel;
  BOOST_CHECK_CLOSE(barrel.trackletMinPt, 0.3f, 1e-6);
  BOOST_CHECK_CLOSE(barrel.cellDeltaTanLambdaSigma, 0.007f, 1e-6);
  BOOST_CHECK_CLOSE(barrel.nSigmaCut, 5.f, 1e-6);
  BOOST_CHECK_CLOSE(barrel.maxChi2ClusterAttachment, 60.f, 1e-6);
  BOOST_CHECK_CLOSE(barrel.maxChi2NDF, 30.f, 1e-6);
  BOOST_CHECK_CLOSE(barrel.pvResolution, 1.e-2f, 1e-6);
  BOOST_CHECK(barrel.isValid());

  // Disk-disk carries its own TrackletMinPt/CellDeltaTanLambdaSigma because
  // the disk path reads them directly (MFT projection and the
  // detector-specific cell-building branch), not shared with the barrel
  // struct. Defaults mirror resetDetectorDefaults(..., MFT) (Configuration.cxx):
  // TrackletMinPt/CellDeltaTanLambdaSigma are left at the TrackingParameters
  // struct defaults (unset by the MFT branch), TrackletMinAbsX is explicitly
  // set to 0.05f there.
  DiskDiskPolicyParams forward;
  BOOST_CHECK_CLOSE(forward.trackletMinPt, 0.3f, 1e-6);
  BOOST_CHECK_CLOSE(forward.cellDeltaTanLambdaSigma, 0.007f, 1e-6);
  BOOST_CHECK_CLOSE(forward.cellRoadRCut, 0.05f, 1e-6);
  BOOST_CHECK_CLOSE(forward.trackletMinAbsX, 0.05f, 1e-6);
  BOOST_CHECK_CLOSE(forward.nSigmaCut, 5.f, 1e-6);
  BOOST_CHECK_CLOSE(forward.maxChi2ClusterAttachment, 60.f, 1e-6);
  BOOST_CHECK_CLOSE(forward.maxChi2NDF, 30.f, 1e-6);
  BOOST_CHECK(forward.isValid());
}

BOOST_AUTO_TEST_CASE(PolicyParamsAbiIsLocked)
{
  BOOST_CHECK_EQUAL(sizeof(CylinderCylinderPolicyParams), 24u);
  BOOST_CHECK_EQUAL(alignof(CylinderCylinderPolicyParams), alignof(float));
  BOOST_CHECK_EQUAL(offsetof(CylinderCylinderPolicyParams, trackletMinPt), 0u);
  BOOST_CHECK_EQUAL(offsetof(CylinderCylinderPolicyParams, cellDeltaTanLambdaSigma), 4u);
  BOOST_CHECK_EQUAL(offsetof(CylinderCylinderPolicyParams, nSigmaCut), 8u);
  BOOST_CHECK_EQUAL(offsetof(CylinderCylinderPolicyParams, maxChi2ClusterAttachment), 12u);
  BOOST_CHECK_EQUAL(offsetof(CylinderCylinderPolicyParams, maxChi2NDF), 16u);
  BOOST_CHECK_EQUAL(offsetof(CylinderCylinderPolicyParams, pvResolution), 20u);

  BOOST_CHECK_EQUAL(sizeof(DiskDiskPolicyParams), 28u);
  BOOST_CHECK_EQUAL(alignof(DiskDiskPolicyParams), alignof(float));
  BOOST_CHECK_EQUAL(offsetof(DiskDiskPolicyParams, trackletMinPt), 0u);
  BOOST_CHECK_EQUAL(offsetof(DiskDiskPolicyParams, cellDeltaTanLambdaSigma), 4u);
  BOOST_CHECK_EQUAL(offsetof(DiskDiskPolicyParams, cellRoadRCut), 8u);
  BOOST_CHECK_EQUAL(offsetof(DiskDiskPolicyParams, trackletMinAbsX), 12u);
  BOOST_CHECK_EQUAL(offsetof(DiskDiskPolicyParams, nSigmaCut), 16u);
  BOOST_CHECK_EQUAL(offsetof(DiskDiskPolicyParams, maxChi2ClusterAttachment), 20u);
  BOOST_CHECK_EQUAL(offsetof(DiskDiskPolicyParams, maxChi2NDF), 24u);
}

BOOST_AUTO_TEST_CASE(PolicyParamBoundsAreValidated)
{
  CylinderCylinderPolicyParams barrel;
  barrel.trackletMinPt = 0.f;
  BOOST_CHECK(!barrel.isValid());

  barrel = CylinderCylinderPolicyParams{};
  barrel.cellDeltaTanLambdaSigma = 0.f;
  BOOST_CHECK(!barrel.isValid());

  barrel = CylinderCylinderPolicyParams{};
  barrel.nSigmaCut = -1.f;
  BOOST_CHECK(!barrel.isValid());

  barrel = CylinderCylinderPolicyParams{};
  barrel.maxChi2ClusterAttachment = 0.f;
  BOOST_CHECK(!barrel.isValid());

  barrel = CylinderCylinderPolicyParams{};
  barrel.maxChi2NDF = -1.f;
  BOOST_CHECK(!barrel.isValid());

  barrel = CylinderCylinderPolicyParams{};
  barrel.pvResolution = 0.f;
  BOOST_CHECK(barrel.isValid());

  barrel.pvResolution = -1.f;
  BOOST_CHECK(!barrel.isValid());

  DiskDiskPolicyParams forward;
  forward.trackletMinPt = 0.f;
  BOOST_CHECK(!forward.isValid());

  forward = DiskDiskPolicyParams{};
  forward.cellDeltaTanLambdaSigma = 0.f;
  BOOST_CHECK(!forward.isValid());

  forward = DiskDiskPolicyParams{};
  forward.cellRoadRCut = 0.f;
  BOOST_CHECK(!forward.isValid());

  forward = DiskDiskPolicyParams{};
  forward.trackletMinAbsX = -1.f;
  BOOST_CHECK(!forward.isValid());

  forward = DiskDiskPolicyParams{};
  forward.nSigmaCut = -1.f;
  BOOST_CHECK(!forward.isValid());

  forward = DiskDiskPolicyParams{};
  forward.maxChi2ClusterAttachment = 0.f;
  BOOST_CHECK(!forward.isValid());

  forward = DiskDiskPolicyParams{};
  forward.maxChi2NDF = -1.f;
  BOOST_CHECK(!forward.isValid());
}

BOOST_AUTO_TEST_CASE(PolicyParamsRejectNonFiniteValues)
{
  const float nan = std::numeric_limits<float>::quiet_NaN();
  const float inf = std::numeric_limits<float>::infinity();

  CylinderCylinderPolicyParams barrel;
  barrel.trackletMinPt = nan;
  BOOST_CHECK(!barrel.isValid());

  barrel = CylinderCylinderPolicyParams{};
  barrel.cellDeltaTanLambdaSigma = inf;
  BOOST_CHECK(!barrel.isValid());

  barrel = CylinderCylinderPolicyParams{};
  barrel.nSigmaCut = -inf;
  BOOST_CHECK(!barrel.isValid());

  barrel = CylinderCylinderPolicyParams{};
  barrel.maxChi2ClusterAttachment = nan;
  BOOST_CHECK(!barrel.isValid());

  barrel = CylinderCylinderPolicyParams{};
  barrel.maxChi2NDF = inf;
  BOOST_CHECK(!barrel.isValid());

  barrel = CylinderCylinderPolicyParams{};
  barrel.pvResolution = nan;
  BOOST_CHECK(!barrel.isValid());

  barrel.pvResolution = inf;
  BOOST_CHECK(!barrel.isValid());

  DiskDiskPolicyParams forward;
  forward.trackletMinPt = nan;
  BOOST_CHECK(!forward.isValid());

  forward = DiskDiskPolicyParams{};
  forward.cellDeltaTanLambdaSigma = inf;
  BOOST_CHECK(!forward.isValid());

  forward = DiskDiskPolicyParams{};
  forward.cellRoadRCut = nan;
  BOOST_CHECK(!forward.isValid());

  forward = DiskDiskPolicyParams{};
  forward.trackletMinAbsX = inf;
  BOOST_CHECK(!forward.isValid());

  forward = DiskDiskPolicyParams{};
  forward.nSigmaCut = nan;
  BOOST_CHECK(!forward.isValid());

  forward = DiskDiskPolicyParams{};
  forward.maxChi2ClusterAttachment = inf;
  BOOST_CHECK(!forward.isValid());

  forward = DiskDiskPolicyParams{};
  forward.maxChi2NDF = nan;
  BOOST_CHECK(!forward.isValid());
}

BOOST_AUTO_TEST_CASE(TraitsAgreeWithStageAStateFamilyFoundation)
{
  using Barrel = TransitionPolicyTraits<TransitionPolicyTag::CylinderCylinder>;
  using Forward = TransitionPolicyTraits<TransitionPolicyTag::DiskDisk>;

  BOOST_CHECK(Barrel::Family == stateFamilyOf(TransitionPolicyTag::CylinderCylinder));
  BOOST_CHECK(Forward::Family == stateFamilyOf(TransitionPolicyTag::DiskDisk));
  BOOST_CHECK(Barrel::Family == StateFamily::Barrel);
  BOOST_CHECK(Forward::Family == StateFamily::Forward);
}

BOOST_AUTO_TEST_CASE(TraitsSelectDistinctStageAStateTypes)
{
  using Barrel = TransitionPolicyTraits<TransitionPolicyTag::CylinderCylinder>;
  using Forward = TransitionPolicyTraits<TransitionPolicyTag::DiskDisk>;

  // Barrel paths use the barrel track state, forward (disk) paths use the
  // forward track state (Architecture.md Stage A); the two must not collapse
  // onto one shared representation.
  BOOST_CHECK(!(std::is_same_v<Barrel::SeedState, Forward::SeedState>));
  BOOST_CHECK(!(std::is_same_v<Barrel::Params, Forward::Params>));
}

BOOST_AUTO_TEST_CASE(PolicySurfaceKindCompatibilityIsExplicit)
{
  BOOST_CHECK(isSurfaceKindCompatible(TransitionPolicyTag::CylinderCylinder, SurfaceKind::Cylinder));
  BOOST_CHECK(!isSurfaceKindCompatible(TransitionPolicyTag::CylinderCylinder, SurfaceKind::Disk));
  BOOST_CHECK(isSurfaceKindCompatible(TransitionPolicyTag::DiskDisk, SurfaceKind::Disk));
  BOOST_CHECK(!isSurfaceKindCompatible(TransitionPolicyTag::DiskDisk, SurfaceKind::Cylinder));
  BOOST_CHECK(!isSurfaceKindCompatible(TransitionPolicyTag::Invalid, SurfaceKind::Cylinder));
  BOOST_CHECK(!isSurfaceKindCompatible(TransitionPolicyTag::Invalid, SurfaceKind::Disk));

  using Barrel = TransitionPolicyTraits<TransitionPolicyTag::CylinderCylinder>;
  using Forward = TransitionPolicyTraits<TransitionPolicyTag::DiskDisk>;
  BOOST_CHECK(isSurfaceKindCompatible(Barrel::Tag, Barrel::ExpectedSurfaceKind));
  BOOST_CHECK(isSurfaceKindCompatible(Forward::Tag, Forward::ExpectedSurfaceKind));
}

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT TrackingKernelParameters
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <cstddef>
#include <limits>
#include <type_traits>

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/detail/TrackingKernelParameters.h"

using namespace o2::itsmft::tracking;

/// Proves the single SurfaceKind-keyed parameter record is a device-compatible
/// POD with a field-level ABI and preserves the old common and family-specific
/// values and finite-value validation contracts.

BOOST_AUTO_TEST_CASE(TrackingKernelParametersAreDeviceCompatiblePods)
{
  BOOST_CHECK(std::is_standard_layout_v<TrackingKernelParameters>);
  BOOST_CHECK(std::is_trivially_copyable_v<TrackingKernelParameters>);
}

BOOST_AUTO_TEST_CASE(TrackingKernelParametersDefaultsMatchCurrentProductionValues)
{
  // Mirrors TrackingParameters defaults in Configuration.h; this is a
  // parity check on the new boundary's defaults, not a claim that
  // production tracking reads from these structs yet.
  TrackingKernelParameters barrel;
  barrel.kind = SurfaceKind::Cylinder;
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
  TrackingKernelParameters forward;
  forward.kind = SurfaceKind::Disk;
  BOOST_CHECK_CLOSE(forward.trackletMinPt, 0.3f, 1e-6);
  BOOST_CHECK_CLOSE(forward.cellDeltaTanLambdaSigma, 0.007f, 1e-6);
  BOOST_CHECK_CLOSE(forward.cellRoadRCut, 0.05f, 1e-6);
  BOOST_CHECK_CLOSE(forward.trackletMinAbsX, 0.05f, 1e-6);
  BOOST_CHECK_CLOSE(forward.nSigmaCut, 5.f, 1e-6);
  BOOST_CHECK_CLOSE(forward.maxChi2ClusterAttachment, 60.f, 1e-6);
  BOOST_CHECK_CLOSE(forward.maxChi2NDF, 30.f, 1e-6);
  BOOST_CHECK(forward.isValid());
}

BOOST_AUTO_TEST_CASE(TrackingKernelParametersAbiIsLocked)
{
  BOOST_CHECK_EQUAL(sizeof(TrackingKernelParameters), 36u);
  BOOST_CHECK_EQUAL(alignof(TrackingKernelParameters), alignof(float));
  BOOST_CHECK_EQUAL(offsetof(TrackingKernelParameters, kind), 0u);
  BOOST_CHECK_EQUAL(offsetof(TrackingKernelParameters, trackletMinPt), 4u);
  BOOST_CHECK_EQUAL(offsetof(TrackingKernelParameters, cellDeltaTanLambdaSigma), 8u);
  BOOST_CHECK_EQUAL(offsetof(TrackingKernelParameters, nSigmaCut), 12u);
  BOOST_CHECK_EQUAL(offsetof(TrackingKernelParameters, maxChi2ClusterAttachment), 16u);
  BOOST_CHECK_EQUAL(offsetof(TrackingKernelParameters, maxChi2NDF), 20u);
  BOOST_CHECK_EQUAL(offsetof(TrackingKernelParameters, pvResolution), 24u);
  BOOST_CHECK_EQUAL(offsetof(TrackingKernelParameters, cellRoadRCut), 28u);
  BOOST_CHECK_EQUAL(offsetof(TrackingKernelParameters, trackletMinAbsX), 32u);
}

BOOST_AUTO_TEST_CASE(TrackingKernelParametersBoundsAreValidated)
{
  TrackingKernelParameters barrel;
  barrel.trackletMinPt = 0.f;
  BOOST_CHECK(!barrel.isValid());

  barrel = TrackingKernelParameters{};
  barrel.cellDeltaTanLambdaSigma = 0.f;
  BOOST_CHECK(!barrel.isValid());

  barrel = TrackingKernelParameters{};
  barrel.nSigmaCut = -1.f;
  BOOST_CHECK(!barrel.isValid());

  barrel = TrackingKernelParameters{};
  barrel.maxChi2ClusterAttachment = 0.f;
  BOOST_CHECK(!barrel.isValid());

  barrel = TrackingKernelParameters{};
  barrel.maxChi2NDF = -1.f;
  BOOST_CHECK(!barrel.isValid());

  barrel = TrackingKernelParameters{};
  barrel.pvResolution = 0.f;
  BOOST_CHECK(barrel.isValid());

  barrel.pvResolution = -1.f;
  BOOST_CHECK(!barrel.isValid());

  TrackingKernelParameters forward;
  forward.kind = SurfaceKind::Disk;
  forward.trackletMinPt = 0.f;
  BOOST_CHECK(!forward.isValid());

  forward = TrackingKernelParameters{};
  forward.kind = SurfaceKind::Disk;
  forward.cellDeltaTanLambdaSigma = 0.f;
  BOOST_CHECK(!forward.isValid());

  forward = TrackingKernelParameters{};
  forward.kind = SurfaceKind::Disk;
  forward.cellRoadRCut = 0.f;
  BOOST_CHECK(!forward.isValid());

  forward = TrackingKernelParameters{};
  forward.kind = SurfaceKind::Disk;
  forward.trackletMinAbsX = -1.f;
  BOOST_CHECK(!forward.isValid());

  forward = TrackingKernelParameters{};
  forward.kind = SurfaceKind::Disk;
  forward.nSigmaCut = -1.f;
  BOOST_CHECK(!forward.isValid());

  forward = TrackingKernelParameters{};
  forward.kind = SurfaceKind::Disk;
  forward.maxChi2ClusterAttachment = 0.f;
  BOOST_CHECK(!forward.isValid());

  forward = TrackingKernelParameters{};
  forward.kind = SurfaceKind::Disk;
  forward.maxChi2NDF = -1.f;
  BOOST_CHECK(!forward.isValid());
}

BOOST_AUTO_TEST_CASE(TrackingKernelParametersRejectNonFiniteValues)
{
  const float nan = std::numeric_limits<float>::quiet_NaN();
  const float inf = std::numeric_limits<float>::infinity();

  TrackingKernelParameters barrel;
  barrel.trackletMinPt = nan;
  BOOST_CHECK(!barrel.isValid());

  barrel = TrackingKernelParameters{};
  barrel.cellDeltaTanLambdaSigma = inf;
  BOOST_CHECK(!barrel.isValid());

  barrel = TrackingKernelParameters{};
  barrel.nSigmaCut = -inf;
  BOOST_CHECK(!barrel.isValid());

  barrel = TrackingKernelParameters{};
  barrel.maxChi2ClusterAttachment = nan;
  BOOST_CHECK(!barrel.isValid());

  barrel = TrackingKernelParameters{};
  barrel.maxChi2NDF = inf;
  BOOST_CHECK(!barrel.isValid());

  barrel = TrackingKernelParameters{};
  barrel.pvResolution = nan;
  BOOST_CHECK(!barrel.isValid());

  barrel.pvResolution = inf;
  BOOST_CHECK(!barrel.isValid());

  TrackingKernelParameters forward;
  forward.kind = SurfaceKind::Disk;
  forward.trackletMinPt = nan;
  BOOST_CHECK(!forward.isValid());

  forward = TrackingKernelParameters{};
  forward.kind = SurfaceKind::Disk;
  forward.cellDeltaTanLambdaSigma = inf;
  BOOST_CHECK(!forward.isValid());

  forward = TrackingKernelParameters{};
  forward.kind = SurfaceKind::Disk;
  forward.cellRoadRCut = nan;
  BOOST_CHECK(!forward.isValid());

  forward = TrackingKernelParameters{};
  forward.kind = SurfaceKind::Disk;
  forward.trackletMinAbsX = inf;
  BOOST_CHECK(!forward.isValid());

  forward = TrackingKernelParameters{};
  forward.kind = SurfaceKind::Disk;
  forward.nSigmaCut = nan;
  BOOST_CHECK(!forward.isValid());

  forward = TrackingKernelParameters{};
  forward.kind = SurfaceKind::Disk;
  forward.maxChi2ClusterAttachment = inf;
  BOOST_CHECK(!forward.isValid());

  forward = TrackingKernelParameters{};
  forward.kind = SurfaceKind::Disk;
  forward.maxChi2NDF = nan;
  BOOST_CHECK(!forward.isValid());
}

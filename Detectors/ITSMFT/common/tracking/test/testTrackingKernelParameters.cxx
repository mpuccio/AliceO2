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

/// Proves the single parameter record is a device-compatible POD with a
/// field-level ABI and common finite-value validation contracts.

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
  TrackingKernelParameters parameters;
  BOOST_CHECK_CLOSE(parameters.trackletMinPt, 0.3f, 1e-6);
  BOOST_CHECK_CLOSE(parameters.nSigmaCut, 5.f, 1e-6);
  BOOST_CHECK_CLOSE(parameters.maxChi2ClusterAttachment, 60.f, 1e-6);
  BOOST_CHECK_CLOSE(parameters.maxChi2NDF, 30.f, 1e-6);
  BOOST_CHECK_CLOSE(parameters.pvResolution, 1.e-2f, 1e-6);
  BOOST_CHECK(parameters.isValid());
}

BOOST_AUTO_TEST_CASE(TrackingKernelParametersAbiIsLocked)
{
  BOOST_CHECK_EQUAL(sizeof(TrackingKernelParameters), 20u);
  BOOST_CHECK_EQUAL(alignof(TrackingKernelParameters), alignof(float));
  BOOST_CHECK_EQUAL(offsetof(TrackingKernelParameters, trackletMinPt), 0u);
  BOOST_CHECK_EQUAL(offsetof(TrackingKernelParameters, nSigmaCut), 4u);
  BOOST_CHECK_EQUAL(offsetof(TrackingKernelParameters, maxChi2ClusterAttachment), 8u);
  BOOST_CHECK_EQUAL(offsetof(TrackingKernelParameters, maxChi2NDF), 12u);
  BOOST_CHECK_EQUAL(offsetof(TrackingKernelParameters, pvResolution), 16u);
}

BOOST_AUTO_TEST_CASE(TrackingKernelParametersBoundsAreValidated)
{
  TrackingKernelParameters parameters;
  parameters.trackletMinPt = 0.f;
  BOOST_CHECK(!parameters.isValid());

  parameters = TrackingKernelParameters{};
  parameters.nSigmaCut = -1.f;
  BOOST_CHECK(!parameters.isValid());

  parameters = TrackingKernelParameters{};
  parameters.maxChi2ClusterAttachment = 0.f;
  BOOST_CHECK(!parameters.isValid());

  parameters = TrackingKernelParameters{};
  parameters.maxChi2NDF = -1.f;
  BOOST_CHECK(!parameters.isValid());

  parameters = TrackingKernelParameters{};
  parameters.pvResolution = 0.f;
  BOOST_CHECK(parameters.isValid());

  parameters.pvResolution = -1.f;
  BOOST_CHECK(!parameters.isValid());
}

BOOST_AUTO_TEST_CASE(TrackingKernelParametersRejectNonFiniteValues)
{
  const float nan = std::numeric_limits<float>::quiet_NaN();
  const float inf = std::numeric_limits<float>::infinity();

  TrackingKernelParameters parameters;
  parameters.trackletMinPt = nan;
  BOOST_CHECK(!parameters.isValid());

  parameters = TrackingKernelParameters{};
  parameters.nSigmaCut = -inf;
  BOOST_CHECK(!parameters.isValid());

  parameters = TrackingKernelParameters{};
  parameters.maxChi2ClusterAttachment = nan;
  BOOST_CHECK(!parameters.isValid());

  parameters = TrackingKernelParameters{};
  parameters.maxChi2NDF = inf;
  BOOST_CHECK(!parameters.isValid());

  parameters = TrackingKernelParameters{};
  parameters.pvResolution = nan;
  BOOST_CHECK(!parameters.isValid());

  parameters.pvResolution = inf;
  BOOST_CHECK(!parameters.isValid());
}

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file TimeFrameLoadFailure.h
/// \brief Typed failures at the ITSMFTTrackingInterface loading boundary
///

#ifndef ALICEO2_ITSMFT_TRACKING_TIMEFRAMELOADFAILURE_H_
#define ALICEO2_ITSMFT_TRACKING_TIMEFRAMELOADFAILURE_H_

// Host loading-boundary only (Architecture.md 7.1): these types, and the
// classification function below, must never enter device views or CA loops.
#ifndef GPUCA_GPUCODE

#include <cstdint>
#include <stdexcept>
#include <string>

#include "ITSMFTTracking/MultiSourceLoading.h"
#include "ITSMFTTracking/SurfaceTiming.h"

namespace o2::itsmft::tracking
{

// Typed, per-TF malformed input at the ITSMFTTrackingInterface loading
// boundary (e.g. a malformed cluster pattern, an invalid ROF range, a
// genuine per-TF BC-arithmetic timing overflow). DropTFUponFailure applies
// to this type and to the two concrete resource-exhaustion exceptions
// (BoundedMemoryResource::MemoryLimitExceeded, std::bad_alloc) -- never to
// TimeFrameLoadException below.
class RecoverableLoadFailure final : public std::runtime_error
{
 public:
  explicit RecoverableLoadFailure(const LoadSourcesResult& result);

  MultiSourceLoadError error() const noexcept { return mResult.error; }
  const LoadSourcesResult& result() const noexcept { return mResult; }

 private:
  LoadSourcesResult mResult;
};

// Why a structural loading-boundary failure occurred, exposed as a typed
// discriminant so callers do not need to string-match TimeFrameLoadException
// ::what(). This is not the discarded, unused v1 FailureKind: every value
// here is actually produced and actually distinguished by a construction
// path below.
enum class TimeFrameLoadFailureReason : uint8_t {
  DictionaryNotConfigured, // ITSMFTTrackingInterface::mDict == nullptr
  NonUniformROFTiming,     // deriveUniformROFTimingConfig() found divergent per-layer length/delay/bias/addTimeErr
  ZeroROFCount,            // configured per-layer timing yields mNROFsTF == 0 (e.g. rofLengthInBC > LHCMaxBunches, or 0 orbits/TF)
  LoadSourcesFailure       // a structural LoadSourcesResult; see loadResult()
};

// Structural loading-boundary failure: catalog/layout not configured, an
// interface-level configuration gap, or a LoadSourcesResult classified
// structural by isRecoverableLoadError(). Never gated by DropTFUponFailure --
// always propagates.
class TimeFrameLoadException final : public std::runtime_error
{
 public:
  // DictionaryNotConfigured / NonUniformROFTiming: no LoadSourcesResult
  // exists yet at the point these are detected, so loadResult() returns a
  // default-constructed (error == None) result for them.
  TimeFrameLoadException(TimeFrameLoadFailureReason reason, std::string message);
  // LoadSourcesFailure: retains the complete LoadSourcesResult that caused
  // it, not just a formatted message.
  explicit TimeFrameLoadException(const LoadSourcesResult& result);

  TimeFrameLoadFailureReason reason() const noexcept { return mReason; }
  // Only meaningful when reason() == LoadSourcesFailure; a default-
  // constructed LoadSourcesResult (error == None) for every other reason.
  const LoadSourcesResult& loadResult() const noexcept { return mLoadResult; }

 private:
  TimeFrameLoadFailureReason mReason;
  LoadSourcesResult mLoadResult{};
};

// Exhaustive, no default case. Omitting `default:` only produces a
// -Wswitch warning for an unhandled enumerator under common build
// configurations, not a hard compile error, so this switch alone does not
// guarantee coverage of a future MultiSourceLoadError/TimingBuildError
// enumerator -- the dedicated classification test
// (ClassifyEveryMultiSourceLoadError) is the actual, checked guarantee.
bool isRecoverableLoadError(MultiSourceLoadError error, TimingBuildError timingDetail) noexcept;

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE
#endif /* ALICEO2_ITSMFT_TRACKING_TIMEFRAMELOADFAILURE_H_ */

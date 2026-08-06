// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/TimeFrameLoadFailure.h"

#include <format>

namespace o2::itsmft::tracking
{

namespace
{
std::string formatLoadSourcesResult(const char* label, const LoadSourcesResult& result)
{
  return std::format("{}: error={} source={} rof={} clusterIndex={} timingDetail={}",
                      label, static_cast<int>(result.error), result.source.value(),
                      result.rof, result.clusterIndex, static_cast<int>(result.timingDetail));
}
} // namespace

RecoverableLoadFailure::RecoverableLoadFailure(const LoadSourcesResult& result)
  : std::runtime_error(formatLoadSourcesResult("TimeFrame loading boundary: recoverable data failure", result)),
    mResult(result)
{
}

TimeFrameLoadException::TimeFrameLoadException(TimeFrameLoadFailureReason reason, std::string message)
  : std::runtime_error(std::move(message)), mReason(reason)
{
}

TimeFrameLoadException::TimeFrameLoadException(const LoadSourcesResult& result)
  : std::runtime_error(formatLoadSourcesResult("TimeFrame loading boundary: structural failure", result)),
    mReason(TimeFrameLoadFailureReason::LoadSourcesFailure),
    mLoadResult(result)
{
}

bool isRecoverableLoadError(MultiSourceLoadError error, TimingBuildError timingDetail) noexcept
{
  switch (error) {
    case MultiSourceLoadError::InvalidROFRange:
    case MultiSourceLoadError::TruncatedExplicitPattern:
    case MultiSourceLoadError::MalformedExplicitPattern:
    case MultiSourceLoadError::InvalidPatternId:
    case MultiSourceLoadError::InvalidSensor:
    case MultiSourceLoadError::InvalidDecodedLayer:
    case MultiSourceLoadError::OtherMalformedInput:
    case MultiSourceLoadError::TrailingPatternData:
      return true;
    case MultiSourceLoadError::TimingError:
      // Overflow is a genuine per-TF BC-arithmetic overflow caused by the
      // incoming ROF data; InvalidROFLength/InvalidSourceROF are
      // configuration problems (an invalid ROFTimingConfig or an
      // unreachable-in-practice source-ROF ordinal), not per-TF data.
      return timingDetail == TimingBuildError::Overflow;
    case MultiSourceLoadError::None:
    case MultiSourceLoadError::NonDenseSourceIds:
    case MultiSourceLoadError::DuplicateSourceId:
    case MultiSourceLoadError::UnsupportedDetector:
    case MultiSourceLoadError::MissingDecoder:
    case MultiSourceLoadError::InvalidLayerMapping:
    case MultiSourceLoadError::DetectorSurfaceMismatch:
    case MultiSourceLoadError::InconsistentDecoderMetadata:
    case MultiSourceLoadError::SurfaceKindMismatch:
    case MultiSourceLoadError::SurfaceCatalogNotConfigured:
    case MultiSourceLoadError::SurfaceCatalogStale:
    case MultiSourceLoadError::MissingDictionary:
    case MultiSourceLoadError::GeometryUnavailable:
    case MultiSourceLoadError::FrameNotConfigured:
      return false;
  }
  return false; // unreachable for any enumerator listed above
}

} // namespace o2::itsmft::tracking

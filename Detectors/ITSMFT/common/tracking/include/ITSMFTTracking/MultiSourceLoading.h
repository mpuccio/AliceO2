// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_MULTISOURCELOADING_H_
#define ALICEO2_ITSMFT_TRACKING_MULTISOURCELOADING_H_

#include <cstdint>
#include <limits>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/SurfaceCatalogView.h"
#include "ITSMFTTracking/SurfaceTiming.h"

namespace o2::itsmft::tracking
{

enum class MultiSourceLoadError : uint8_t {
  None,
  NonDenseSourceIds,
  DuplicateSourceId,
  UnsupportedDetector,
  MissingDecoder,
  InvalidROFRange,
  InvalidLayerMapping,
  DetectorSurfaceMismatch,
  // The decoder-produced SurfaceMeasurement disagrees with the request that
  // produced it (surface not equal to layerToSurface[layer], or cluster
  // source/index/sourceROF/sensor-detector not equal to what was requested).
  // This guards cluster identity against a buggy host adapter; it is
  // distinct from DetectorSurfaceMismatch, which is a layout/source
  // detector-qualification mismatch rather than a decoder self-consistency
  // failure.
  InconsistentDecoderMetadata,
  // The decoder's declared SurfaceKind (cylinder/disk) does not match the
  // target surface's descriptor. Never inferred from surface count.
  SurfaceKindMismatch,
  TimingError,
  // Normalized loading preflight (loadNormalizedSource()): the caller-
  // supplied SurfaceCatalogView is empty (no surfaces) -- TimeFrame owns no
  // catalog of its own (Gate 4 B2 Slice 2), so this means the caller (the
  // plan's owner) had no plan to pass.
  SurfaceCatalogNotConfigured,
  // Retained for enum-value compatibility with the pre-Slice-2 contract.
  // Unreachable via loadNormalizedSource() itself: a caller-supplied
  // SurfaceCatalogView has no currency concept to go stale against (there is
  // no TimeFrame-owned epoch/rebuild-on-demand any more).
  SurfaceCatalogStale,
  // Appended to preserve the numeric values of the pre-existing contract.
  MissingDictionary,
  TruncatedExplicitPattern,
  MalformedExplicitPattern,
  InvalidPatternId,
  InvalidSensor,
  InvalidDecodedLayer,
  GeometryUnavailable,
  OtherMalformedInput,
  TrailingPatternData,
  FrameNotConfigured
};

struct LoadSourcesResult {
  MultiSourceLoadError error{MultiSourceLoadError::None};
  ClusterSourceId source{};
  uint32_t rof{std::numeric_limits<uint32_t>::max()};
  uint32_t clusterIndex{std::numeric_limits<uint32_t>::max()};
  // Meaningful only when error == MultiSourceLoadError::TimingError: the
  // underlying computeROFIntervalBC() failure. None for every other error
  // value, including MultiSourceLoadError::None on success. Callers must not
  // classify TimingError without inspecting this field -- InvalidROFLength
  // and InvalidSourceROF are configuration problems, while Overflow is a
  // genuine per-TF BC-arithmetic overflow caused by the incoming ROF data.
  TimingBuildError timingDetail{TimingBuildError::None};

  bool ok() const noexcept { return error == MultiSourceLoadError::None; }
};

// Loads one or more detector-qualified sources into `frame`. Geometry
// lookup, pattern consumption, covariance calculation and systematic-error
// application occur exactly once per cluster. Transactional: on any
// validation failure `frame` is left completely unmodified (whatever it
// held before the call) and the returned result identifies the first
// offending source/ROF/cluster. Pattern buffers must be consumed exactly;
// trailing data is reported at the after-last ROF and cluster ordinals.
// `origin` is the single explicit
// InteractionRecord chosen for the loaded frame; every source ROF is
// converted to a signed TF-relative BC interval against it.
LoadSourcesResult loadSources(MultiSourceFrame& frame,
                              const SurfaceCatalogView& catalog,
                              gsl::span<const ClusterSourceInput> sources,
                              const o2::InteractionRecord& origin);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_MULTISOURCELOADING_H_ */

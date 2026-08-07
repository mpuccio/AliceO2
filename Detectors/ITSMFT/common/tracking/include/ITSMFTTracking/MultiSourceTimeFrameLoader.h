// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_MULTISOURCETIMEFRAMELOADER_H_
#define ALICEO2_ITSMFT_TRACKING_MULTISOURCETIMEFRAMELOADER_H_

#include <cstdint>
#include <limits>

#ifndef GPUCA_GPUCODE
#include <stdexcept>
#include <string>
#endif

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/ROFViews.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceId.h"
#include "ITSMFTTracking/SurfaceTiming.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

namespace o2::itsmft::tracking
{

class TimeFrame;

// Non-owning, host-side description of one detector-qualified input stream.
// Source IDs are unique, dense, TimeFrame-local, and stable for the loaded
// frame; the all-ones value is reserved as invalid. The source retains
// ownership of detector-specific inputs; loading copies only normalized
// measurements and timing into the staged frame.
struct ClusterSourceInput {
  ClusterSourceId id{};
  o2::detectors::DetID::ID detector{o2::detectors::DetID::ITS};
  gsl::span<const o2::itsmft::CompClusterExt> clusters{};
  // Exact source-local serialization of explicit and grouped patterns.
  gsl::span<const unsigned char> patterns{};
  // `rofs` must form an exact, ordered partition of `clusters`. An empty
  // `rofs` is valid only when `clusters` is also empty.
  gsl::span<const o2::itsmft::ROFRecord> rofs{};
  const o2::itsmft::TopologyDictionary* dictionary{nullptr};
  const o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels{nullptr};
  // Detector-local layer (as discovered by the decoder) -> global SurfaceId.
  gsl::span<const SurfaceId> layerToSurface{};
  ROFTimingConfig timing{};
  const ClusterDecoder* decoder{nullptr};
  bool applySysErrors{true};
  // Non-owning timing/mask context assembled by the adapter for this event.
  RuntimeROFViews rofViews{};
};

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
  // produced it (surface, cluster identity, source ROF, or sensor detector).
  InconsistentDecoderMetadata,
  // The decoder's declared SurfaceKind does not match the target descriptor.
  SurfaceKindMismatch,
  TimingError,
  // The caller supplied no SurfaceCatalogView (there is no frame-owned
  // catalog to rebuild or compare against).
  SurfaceCatalogNotConfigured,
  // Retained for enum-value compatibility; unreachable for a caller-supplied
  // SurfaceCatalogView, which has no currency concept.
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
  // Meaningful only when error == TimingError. None for all other values.
  TimingBuildError timingDetail{TimingBuildError::None};

  bool ok() const noexcept { return error == MultiSourceLoadError::None; }
};

// Decode and stage all detector-qualified sources into a normalized frame.
// Geometry lookup, pattern consumption, covariance calculation, and
// systematic-error application occur once per cluster. On failure, `frame`
// remains unmodified and the result identifies the first offending input.
LoadSourcesResult loadSources(MultiSourceFrame& frame,
                              const SurfaceCatalogView& catalog,
                              gsl::span<const ClusterSourceInput> sources,
                              const o2::InteractionRecord& origin);

// Non-owning direct loading component. It stages every source and installs
// one complete event into a configured TimeFrame; it does not track, publish,
// own raw ROFs, or make workflow decisions. A failure leaves the frame's
// previous event and configuration untouched.
class MultiSourceTimeFrameLoader
{
 public:
  static LoadSourcesResult load(TimeFrame& frame, gsl::span<const ClusterSourceInput> sources,
                                SurfaceCatalogView catalog, const o2::InteractionRecord& origin);
};

#ifndef GPUCA_GPUCODE

// Typed, per-TimeFrame malformed input at the workflow loading boundary.
// DropTFUponFailure applies to this type and resource-exhaustion exceptions,
// never to TimeFrameLoadException below.
class RecoverableLoadFailure final : public std::runtime_error
{
 public:
  explicit RecoverableLoadFailure(const LoadSourcesResult& result);

  MultiSourceLoadError error() const noexcept { return mResult.error; }
  const LoadSourcesResult& result() const noexcept { return mResult; }

 private:
  LoadSourcesResult mResult;
};

enum class TimeFrameLoadFailureReason : uint8_t {
  DictionaryNotConfigured,
  NonUniformROFTiming,
  ZeroROFCount,
  LoadSourcesFailure
};

// Structural loading-boundary failure: configuration gaps or a structural
// LoadSourcesResult. It always propagates and is never gated by drop policy.
class TimeFrameLoadException final : public std::runtime_error
{
 public:
  TimeFrameLoadException(TimeFrameLoadFailureReason reason, std::string message);
  explicit TimeFrameLoadException(const LoadSourcesResult& result);

  TimeFrameLoadFailureReason reason() const noexcept { return mReason; }
  const LoadSourcesResult& loadResult() const noexcept { return mLoadResult; }

 private:
  TimeFrameLoadFailureReason mReason;
  LoadSourcesResult mLoadResult{};
};

// Classify a decoder/loading result for workflow drop-versus-propagate policy.
bool isRecoverableLoadError(MultiSourceLoadError error, TimingBuildError timingDetail) noexcept;

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking

#endif

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.

#ifndef ALICEO2_ITSMFT_TRACKING_COMMONTRACKOUTPUTADAPTER_H_
#define ALICEO2_ITSMFT_TRACKING_COMMONTRACKOUTPUTADAPTER_H_

// A deliberately pure, host-side boundary for the later opt-in DPL adapters.
// It consumes immutable owner data plus workflow-owned ROF context and returns
// only fully staged vectors.  In particular this file has no Framework/DPL API.

#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <vector>

#include <gsl/span>

#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/ClockTimingPublicationView.h"
#include "ITSMFTTracking/MCLabelAccumulator.h"
#include "ITSMFTTracking/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/SurfaceKinematicStateLegacyAdapters.h"
#include "ITSMFTTracking/TimeFrame.h"

namespace o2::itsmft::tracking
{

enum class CommonTrackOutputAdapterError : uint8_t {
  None,
  TooManyCommonTracks,
  InvalidTrackRange,
  UnresolvedReference,
  MixedDetector,
  MixedSources,
  InvalidExternalClusterIndex,
  InvalidLayerLayout,
  InvalidTimestamp,
  InvalidROF,
  InvalidState,
  MissingCompatibility
};

struct CommonTrackOutputAdapterSelection {
  std::vector<uint32_t> globalIndices;
};

struct CommonTrackOutputOrderEntry {
  uint32_t globalIndex{};
  o2::its::TimeStamp timestamp{};
};

// This context is intentionally source-local.  ROFRecord payload is copied
// only into the returned publication product, never into TimeFrame.
struct CommonTrackOutputTimingContext {
  gsl::span<const o2::itsmft::ROFRecord> inputROFs;
  ClockTimingPublicationView clock;
};

struct CommonTrackPublicationContext {
  o2::detectors::DetID::ID detector{};
  ClusterSourceId source{};
  gsl::span<const o2::itsmft::ROFRecord> inputROFs;
  ClockTimingPublicationView clock;
  gsl::span<const SurfaceId> orderedSurfaces;
};

struct ITSCommonTrackOutput {
  std::vector<o2::its::TrackITS> tracks;
  std::vector<int> clusterIndices;
  std::vector<o2::itsmft::ROFRecord> trackROFs;
  std::vector<o2::MCCompLabel> labels;
};

struct MFTCommonTrackOutput {
  std::vector<o2::mft::TrackMFT> tracks;
  std::vector<int> clusterIndices;
  std::vector<o2::itsmft::ROFRecord> trackROFs;
  std::vector<uint16_t> seedPatterns;
  std::vector<o2::MCCompLabel> labels;
};

inline std::optional<CommonTrackOutputAdapterSelection> selectCommonTracksForSource(
  const TimeFrame& frame, o2::detectors::DetID::ID detector, ClusterSourceId source,
  CommonTrackOutputAdapterError& error)
{
  error = CommonTrackOutputAdapterError::None;
  const auto& tracks = frame.getCommonTracks();
  if (tracks.size() > std::numeric_limits<uint32_t>::max()) {
    error = CommonTrackOutputAdapterError::TooManyCommonTracks;
    return std::nullopt;
  }
  CommonTrackOutputAdapterSelection selection;
  const auto& references = frame.getTrackClusterIndices();
  const auto& normalized = frame.getNormalizedFrame();
  selection.globalIndices.reserve(tracks.size());
  for (uint32_t globalIndex = 0; globalIndex < tracks.size(); ++globalIndex) {
    const auto& track = tracks[globalIndex];
    if (!isValidTrackRange(track, static_cast<uint32_t>(references.size()))) {
      error = CommonTrackOutputAdapterError::InvalidTrackRange;
      return std::nullopt;
    }
    bool requested = false;
    bool foreign = false;
    for (uint32_t i = track.firstClusterRef; i < track.clusterRefEnd; ++i) {
      const auto& reference = references[i];
      const auto* measurement = normalized.getMeasurement(reference.surface, reference.index);
      if (measurement == nullptr || measurement->surface != reference.surface) {
        error = CommonTrackOutputAdapterError::UnresolvedReference;
        return std::nullopt;
      }
      const bool match = measurement->sensor.detector == static_cast<uint32_t>(detector) && measurement->cluster.source == source;
      requested |= match;
      foreign |= !match;
      if (measurement->cluster.index > static_cast<uint32_t>(std::numeric_limits<int>::max())) {
        error = CommonTrackOutputAdapterError::InvalidExternalClusterIndex;
        return std::nullopt;
      }
    }
    if (requested && foreign) {
      error = CommonTrackOutputAdapterError::MixedDetector;
      return std::nullopt;
    }
    if (requested) {
      selection.globalIndices.push_back(globalIndex);
    }
  }
  return selection;
}

inline std::optional<o2::its::TimeStamp> makeOutputTimestamp(const CommonTrackTimestamp& timestamp,
                                                             const ClockTimingPublicationView& clock,
                                                             CommonTrackOutputAdapterError& error)
{
  const auto result = clock.makeOutputTimestamp(timestamp);
  if (!result) {
    error = CommonTrackOutputAdapterError::InvalidTimestamp;
    return std::nullopt;
  }
  return result;
}

inline std::optional<std::vector<CommonTrackOutputOrderEntry>> makeLegacyOutputOrder(
  const TimeFrame& frame, const CommonTrackOutputAdapterSelection& selection,
  const ClockTimingPublicationView& clock, CommonTrackOutputAdapterError& error)
{
  std::vector<CommonTrackOutputOrderEntry> ordered;
  ordered.reserve(selection.globalIndices.size());
  for (const auto index : selection.globalIndices) {
    const auto timestamp = makeOutputTimestamp(frame.getCommonTracks()[index].timestamp, clock, error);
    if (!timestamp) {
      return std::nullopt;
    }
    ordered.push_back({index, *timestamp});
  }
  // Mirror Tracker<NLayers>::sortTracks(): it reorders the scratch results
  // after all accepted CommonTrack shadows were appended, by the lower edge
  // of the legacy symmetric/clamped timestamp and then chi2.
  std::sort(ordered.begin(), ordered.end(), [&frame](const auto& left, const auto& right) {
    const auto& leftTrack = frame.getCommonTracks()[left.globalIndex];
    const auto& rightTrack = frame.getCommonTracks()[right.globalIndex];
    const auto leftLower = left.timestamp.getTimeStamp() - left.timestamp.getTimeStampError();
    const auto rightLower = right.timestamp.getTimeStamp() - right.timestamp.getTimeStampError();
    if (leftLower != rightLower) {
      return leftLower < rightLower;
    }
    return leftTrack.chi2 < rightTrack.chi2;
  });
  return ordered;
}

inline bool finalizeROFs(std::vector<o2::itsmft::ROFRecord>& rofs, const std::vector<o2::its::TimeStamp>& times,
                         const CommonTrackOutputTimingContext& context, CommonTrackOutputAdapterError& error)
{
  for (auto& rof : rofs) {
    rof.setFirstEntry(0);
    rof.setNEntries(0);
  }
  for (const auto& time : times) {
    const int rof = context.clock.getROF(time);
    if (rof < 0 || static_cast<size_t>(rof) >= rofs.size()) {
      // Preserve fillITSOutputs()/fillMFTOutputs(): a legacy accepted track
      // is still published when its reconstructed timestamp lies outside
      // the workflow ROF span; only its TrackROF entry is omitted.
      continue;
    }
    rofs[rof].setNEntries(rofs[rof].getNEntries() + 1);
  }
  std::vector<int> counts(rofs.size());
  for (size_t i = 0; i < rofs.size(); ++i) {
    counts[i] = rofs[i].getNEntries();
  }
  std::exclusive_scan(counts.begin(), counts.end(), counts.begin(), 0);
  for (size_t i = 0; i < rofs.size(); ++i) {
    rofs[i].setFirstEntry(counts[i]);
  }
  return true;
}

inline void setOutputClusterRange(o2::its::TrackITS& track, int first, int count)
{
  track.setClusterRefs(first, count);
}

inline void setOutputClusterRange(o2::mft::TrackMFT& track, int first, int count)
{
  track.setExternalClusterIndexOffset(first);
  track.setNumberOfPoints(count);
}

template <typename OutputTrack>
inline bool collectReferences(const TimeFrame& frame, const CommonTrack& common, gsl::span<const SurfaceId> orderedSurfaces,
                              uint32_t maxLayers, std::vector<int>& outputIndices, OutputTrack& output,
                              MCLabelAccumulator* labels, uint32_t& pattern, CommonTrackOutputAdapterError& error)
{
  const auto& references = frame.getTrackClusterIndices();
  const auto& normalized = frame.getNormalizedFrame();
  std::vector<const SurfaceMeasurement*> byLayer(maxLayers, nullptr);
  for (uint32_t ref = common.firstClusterRef; ref < common.clusterRefEnd; ++ref) {
    const auto& key = references[ref];
    const auto* measurement = normalized.getMeasurement(key.surface, key.index);
    if (measurement == nullptr) {
      error = CommonTrackOutputAdapterError::UnresolvedReference;
      return false;
    }
    const auto where = std::find(orderedSurfaces.begin(), orderedSurfaces.end(), key.surface);
    if (where == orderedSurfaces.end() || static_cast<uint32_t>(where - orderedSurfaces.begin()) >= maxLayers) {
      error = CommonTrackOutputAdapterError::InvalidLayerLayout;
      return false;
    }
    const auto layer = static_cast<uint32_t>(where - orderedSurfaces.begin());
    if (byLayer[layer] != nullptr) {
      error = CommonTrackOutputAdapterError::InvalidLayerLayout;
      return false;
    }
    byLayer[layer] = measurement;
  }
  const int first = static_cast<int>(outputIndices.size());
  uint32_t count = 0;
  for (uint32_t layer = maxLayers; layer-- > 0;) {
    const auto* measurement = byLayer[layer];
    if (measurement == nullptr) {
      continue;
    }
    outputIndices.push_back(static_cast<int>(measurement->cluster.index));
    output.setClusterSize(layer, measurement->shape.nPixels);
    pattern |= 1u << layer;
    ++count;
  }
  setOutputClusterRange(output, first, static_cast<int>(count));
  if (labels != nullptr) {
    for (const auto* measurement : byLayer) {
      if (measurement != nullptr) {
        labels->addCluster(normalized.getLabels(measurement->cluster));
      }
    }
  }
  return true;
}

inline std::optional<ITSCommonTrackOutput> stageITSCommonTrackOutput(const TimeFrame& frame, ClusterSourceId source,
                                                                     gsl::span<const SurfaceId> surfaces,
                                                                     const CommonTrackOutputTimingContext& context,
                                                                     const ITSSharedClusterCompatibility& compatibility,
                                                                     bool withMC, CommonTrackOutputAdapterError& error)
{
  const auto selection = selectCommonTracksForSource(frame, o2::detectors::DetID::ITS, source, error);
  if (!selection || (!selection->globalIndices.empty() && !compatibility.isSealed())) {
    if (error == CommonTrackOutputAdapterError::None)
      error = CommonTrackOutputAdapterError::MissingCompatibility;
    return std::nullopt;
  }
  const auto ordered = makeLegacyOutputOrder(frame, *selection, context.clock, error);
  if (!ordered) {
    return std::nullopt;
  }
  ITSCommonTrackOutput staged;
  staged.trackROFs.assign(context.inputROFs.begin(), context.inputROFs.end());
  staged.tracks.reserve(ordered->size());
  staged.labels.reserve(withMC ? ordered->size() : 0);
  std::vector<o2::its::TimeStamp> times;
  times.reserve(ordered->size());
  for (const auto& orderedTrack : *ordered) {
    const auto index = orderedTrack.globalIndex;
    o2::track::TrackParCovF inner, outer;
    const auto& common = frame.getCommonTracks()[index];
    if (!legacy::exportBarrelTrackParCov(common.innerState, inner) || !legacy::exportBarrelTrackParCov(common.outerState, outer)) {
      error = CommonTrackOutputAdapterError::InvalidState;
      return std::nullopt;
    }
    const auto it = std::lower_bound(compatibility.entries().begin(), compatibility.entries().end(), index,
                                     [](const auto& entry, uint32_t value) { return entry.commonTrackIndex < value; });
    if (it == compatibility.entries().end() || it->commonTrackIndex != index) {
      error = CommonTrackOutputAdapterError::MissingCompatibility;
      return std::nullopt;
    }
    o2::its::TrackITS output{inner, common.chi2, outer};
    MCLabelAccumulator accumulator;
    uint32_t pattern = 0;
    if (!collectReferences(frame, common, surfaces, 7, staged.clusterIndices, output, withMC ? &accumulator : nullptr, pattern, error))
      return std::nullopt;
    output.setPattern(pattern);
    output.setSharedClusters(it->hasSharedClusters);
    output.getTimeStamp() = orderedTrack.timestamp;
    staged.tracks.push_back(std::move(output));
    times.push_back(orderedTrack.timestamp);
    if (withMC)
      staged.labels.push_back(accumulator.finalize());
  }
  if (!finalizeROFs(staged.trackROFs, times, context, error))
    return std::nullopt;
  return staged;
}

inline std::optional<MFTCommonTrackOutput> stageMFTCommonTrackOutput(const TimeFrame& frame, ClusterSourceId source,
                                                                     gsl::span<const SurfaceId> surfaces,
                                                                     const CommonTrackOutputTimingContext& context,
                                                                     const MFTPublicationCompatibility& compatibility,
                                                                     bool withMC, CommonTrackOutputAdapterError& error)
{
  const auto selection = selectCommonTracksForSource(frame, o2::detectors::DetID::MFT, source, error);
  if (!selection)
    return std::nullopt;
  const auto ordered = makeLegacyOutputOrder(frame, *selection, context.clock, error);
  if (!ordered) {
    return std::nullopt;
  }
  MFTCommonTrackOutput staged;
  staged.trackROFs.assign(context.inputROFs.begin(), context.inputROFs.end());
  staged.tracks.reserve(ordered->size());
  staged.seedPatterns.reserve(ordered->size());
  std::vector<o2::its::TimeStamp> times;
  times.reserve(ordered->size());
  for (const auto& orderedTrack : *ordered) {
    const auto index = orderedTrack.globalIndex;
    const auto& common = frame.getCommonTracks()[index];
    const auto* sidecar = compatibility.find(index, frame.getCommonTracks().size());
    if (sidecar == nullptr) {
      error = CommonTrackOutputAdapterError::MissingCompatibility;
      return std::nullopt;
    }
    o2::track::TrackParCovFwd inner, outer;
    if (!legacy::exportLegacyForwardTrackParCov(common.innerState, inner) || !legacy::exportLegacyForwardTrackParCov(common.outerState, outer)) {
      error = CommonTrackOutputAdapterError::InvalidState;
      return std::nullopt;
    }
    // TrackMFT persists the fitted chi2 in both its base and out-parameter
    // TrackParCovFwd objects. CommonTrack owns that one fitted value, so the
    // output reconstruction must restore it on the exported outer state too.
    outer.setTrackChi2(sidecar->outParamChi2);
    o2::mft::TrackMFT output;
    static_cast<o2::track::TrackParCovFwd&>(output) = inner;
    output.setOutParam(outer);
    output.setTrackChi2(common.chi2);
    output.setCA(true);
    output.setInvQPtSeed(sidecar->invQPtSeed);
    output.setChi2QPtSeed(sidecar->chi2QPtSeed);
    MCLabelAccumulator accumulator;
    uint32_t ignoredPattern = 0;
    if (!collectReferences(frame, common, surfaces, 10, staged.clusterIndices, output, withMC ? &accumulator : nullptr, ignoredPattern, error))
      return std::nullopt;
    staged.tracks.push_back(std::move(output));
    staged.seedPatterns.push_back(sidecar->seedPattern);
    times.push_back(orderedTrack.timestamp);
    if (withMC)
      staged.labels.push_back(accumulator.finalize());
  }
  if (!finalizeROFs(staged.trackROFs, times, context, error))
    return std::nullopt;
  return staged;
}

inline std::optional<ITSCommonTrackOutput> stageITSCommonTrackOutput(const TimeFrame& frame, const CommonTrackPublicationContext& context,
                                                                     const ITSSharedClusterCompatibility& compatibility, bool withMC,
                                                                     CommonTrackOutputAdapterError& error)
{
  if (context.detector != o2::detectors::DetID::ITS) {
    error = CommonTrackOutputAdapterError::MixedDetector;
    return std::nullopt;
  }
  return stageITSCommonTrackOutput(frame, context.source, context.orderedSurfaces, {context.inputROFs, context.clock}, compatibility, withMC, error);
}

inline std::optional<MFTCommonTrackOutput> stageMFTCommonTrackOutput(const TimeFrame& frame, const CommonTrackPublicationContext& context,
                                                                     const MFTPublicationCompatibility& compatibility, bool withMC,
                                                                     CommonTrackOutputAdapterError& error)
{
  if (context.detector != o2::detectors::DetID::MFT) {
    error = CommonTrackOutputAdapterError::MixedDetector;
    return std::nullopt;
  }
  return stageMFTCommonTrackOutput(frame, context.source, context.orderedSurfaces, {context.inputROFs, context.clock}, compatibility, withMC, error);
}

} // namespace o2::itsmft::tracking

#endif // ALICEO2_ITSMFT_TRACKING_COMMONTRACKOUTPUTADAPTER_H_

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/MultiSourceLoading.h"

#include <vector>

#include "DataFormatsITSMFT/ROFRecord.h"

namespace o2::itsmft::tracking
{

namespace
{
bool isSupportedDetector(o2::detectors::DetID::ID det) noexcept
{
  return det == o2::detectors::DetID::ITS || det == o2::detectors::DetID::MFT;
}

MultiSourceLoadError mapDecodeError(ClusterDecodeError error) noexcept
{
  switch (error) {
    case ClusterDecodeError::None:
      return MultiSourceLoadError::None;
    case ClusterDecodeError::MissingDictionary:
      return MultiSourceLoadError::MissingDictionary;
    case ClusterDecodeError::TruncatedExplicitPattern:
      return MultiSourceLoadError::TruncatedExplicitPattern;
    case ClusterDecodeError::MalformedExplicitPattern:
      return MultiSourceLoadError::MalformedExplicitPattern;
    case ClusterDecodeError::InvalidPatternId:
      return MultiSourceLoadError::InvalidPatternId;
    case ClusterDecodeError::InvalidSensor:
      return MultiSourceLoadError::InvalidSensor;
    case ClusterDecodeError::InvalidLayer:
      return MultiSourceLoadError::InvalidDecodedLayer;
    case ClusterDecodeError::InvalidLayerMapping:
      return MultiSourceLoadError::InvalidLayerMapping;
    case ClusterDecodeError::GeometryUnavailable:
      return MultiSourceLoadError::GeometryUnavailable;
    case ClusterDecodeError::OtherMalformedInput:
      return MultiSourceLoadError::OtherMalformedInput;
  }
  return MultiSourceLoadError::OtherMalformedInput;
}
} // namespace

LoadSourcesResult loadSources(MultiSourceFrame& frame,
                              const SurfaceCatalogView& catalog,
                              gsl::span<const ClusterSourceInput> sources,
                              const o2::InteractionRecord& origin)
{
  const auto nSources = static_cast<uint32_t>(sources.size());

  // 1. Validate source-level contract: dense/unique IDs, supported detector,
  //    a decoder is present. Nothing below this point touches `frame`.
  std::vector<bool> seen(nSources, false);
  for (const auto& src : sources) {
    if (!src.id.isValid() || src.id.value() >= nSources) {
      return {MultiSourceLoadError::NonDenseSourceIds, src.id};
    }
    if (seen[src.id.value()]) {
      return {MultiSourceLoadError::DuplicateSourceId, src.id};
    }
    seen[src.id.value()] = true;
    if (!isSupportedDetector(src.detector)) {
      return {MultiSourceLoadError::UnsupportedDetector, src.id};
    }
    if (src.decoder == nullptr) {
      return {MultiSourceLoadError::MissingDecoder, src.id};
    }
    if (!src.clusters.empty() && src.dictionary == nullptr) {
      return {MultiSourceLoadError::MissingDictionary, src.id, 0, 0};
    }
  }

  // Staging storage: nothing here is visible to `frame` until every source
  // has been fully validated and decoded.
  std::vector<std::vector<SurfaceMeasurement>> perSurface(catalog.nSurfaces);
  std::vector<SourceMetadata> sourcesMeta(nSources);
  std::vector<std::vector<ROFIntervalBC>> perSourceIntervals(nSources);
  std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*> labelSources(nSources, nullptr);

  for (const auto& src : sources) {
    const auto srcIdx = src.id.value();
    sourcesMeta[srcIdx] = SourceMetadata{src.id, src.detector, static_cast<uint32_t>(src.rofs.size())};
    labelSources[srcIdx] = src.labels;

    // 2. Validate ROF cluster ranges against this source's cluster span
    //    before any decoding, so an out-of-range access can never happen.
    //    The contract (documented on ClusterSourceInput::rofs) is an exact,
    //    ordered partition of the cluster span: the first range begins at
    //    zero, each range begins exactly where the previous one ended, and
    //    the final range ends exactly at clusters.size(). Leading gaps,
    //    internal gaps and trailing unreferenced clusters are all rejected,
    //    because explicit-pattern consumption follows ROF traversal order
    //    and a gap would silently associate the wrong pattern bytes with
    //    the clusters that follow it.
    int64_t expectedNext = 0;
    for (uint32_t r = 0; r < src.rofs.size(); ++r) {
      const auto& rof = src.rofs[r];
      const int64_t first = rof.getFirstEntry();
      const int64_t n = rof.getNEntries();
      if (n < 0 || first != expectedNext) {
        return {MultiSourceLoadError::InvalidROFRange, src.id, r};
      }
      expectedNext = first + n;
      if (expectedNext > static_cast<int64_t>(src.clusters.size())) {
        return {MultiSourceLoadError::InvalidROFRange, src.id, r};
      }
    }
    if (expectedNext != static_cast<int64_t>(src.clusters.size())) {
      return {MultiSourceLoadError::InvalidROFRange, src.id, static_cast<uint32_t>(src.rofs.size())};
    }

    // 3. Timing: one independent config per source; sourceROF is the
    //    source-local ordinal position of the ROF within `src.rofs`.
    auto& intervals = perSourceIntervals[srcIdx];
    intervals.reserve(src.rofs.size());
    for (uint32_t r = 0; r < src.rofs.size(); ++r) {
      const auto built = computeROFIntervalBC(src.rofs[r].getBCData(), origin, src.timing, r);
      if (!built.ok()) {
        // Designated initialization, not {error, id, r, {}}: the latter
        // would explicitly zero clusterIndex instead of leaving its
        // invalid-sentinel default member initializer in place, which is
        // not meaningful for a timing failure (timing is validated before
        // any cluster is decoded).
        return LoadSourcesResult{.error = MultiSourceLoadError::TimingError, .source = src.id, .rof = r, .timingDetail = built.error};
      }
      intervals.push_back(built.interval);
    }

    // 4. Decode: one independent pattern cursor per source. Each cluster is
    //    decoded exactly once by the source's decoder.
    src.decoder->prepare();
    BoundedPatternCursor patterns{src.patterns};
    for (uint32_t r = 0; r < src.rofs.size(); ++r) {
      const auto& rof = src.rofs[r];
      const auto firstEntry = rof.getFirstEntry();
      const auto nEntries = rof.getNEntries();
      for (int32_t clusterId = firstEntry; clusterId < firstEntry + nEntries; ++clusterId) {
        const auto& cluster = src.clusters[clusterId];
        const auto externalIndex = static_cast<uint32_t>(clusterId);
        const auto decoded = src.decoder->decode(cluster, patterns, src.dictionary, src.layerToSurface,
                                                 src.id, externalIndex, r, src.applySysErrors);
        if (!decoded.ok()) {
          return {mapDecodeError(decoded.error), src.id, r, externalIndex};
        }
        // A buggy decoder could report layerMapped=true while its own
        // `layer` is negative or out of range for `src.layerToSurface`;
        // layerMapped is never trusted on its own before indexing.
        if (!decoded.layerMapped || decoded.layer < 0 ||
            static_cast<size_t>(decoded.layer) >= src.layerToSurface.size()) {
          return {MultiSourceLoadError::InvalidLayerMapping, src.id, r, externalIndex};
        }
        const auto expectedSurface = src.layerToSurface[decoded.layer];
        if (!expectedSurface.isValid() || expectedSurface.value() >= catalog.nSurfaces) {
          return {MultiSourceLoadError::InvalidLayerMapping, src.id, r, externalIndex};
        }
        // Authoritative cluster identity: a buggy host adapter must not be
        // able to substitute a different surface, source, external index,
        // source-local ROF or sensor detector than the ones it was asked to
        // decode against.
        const auto& measurement = decoded.measurement;
        if (measurement.surface != expectedSurface ||
            measurement.cluster.source != src.id ||
            measurement.cluster.index != externalIndex ||
            measurement.sourceROF != r ||
            measurement.sensor.detector != static_cast<uint32_t>(src.detector)) {
          return {MultiSourceLoadError::InconsistentDecoderMetadata, src.id, r, externalIndex};
        }
        const auto& surfaceDescriptor = catalog.getSurface(expectedSurface);
        if (surfaceDescriptor.detectorId != static_cast<uint8_t>(src.detector)) {
          return {MultiSourceLoadError::DetectorSurfaceMismatch, src.id, r, externalIndex};
        }
        // Explicit surface-kind check: geometry kind is never inferred from
        // surface count, only compared between the decoder's own declared
        // kind and the layout's explicit descriptor.
        if (surfaceDescriptor.kind != decoded.kind) {
          return {MultiSourceLoadError::SurfaceKindMismatch, src.id, r, externalIndex};
        }
        perSurface[expectedSurface.value()].push_back(measurement);
      }
    }
    // The source pattern buffer is an exact serialization of the explicit
    // and grouped patterns encountered above. Any bytes left after the last
    // cluster are inconsistent input and are rejected, with the after-last
    // ROF/cluster ordinals identifying the boundary where consumption ended.
    if (!patterns.empty()) {
      return {MultiSourceLoadError::TrailingPatternData, src.id,
              static_cast<uint32_t>(src.rofs.size()), static_cast<uint32_t>(src.clusters.size())};
    }
  }

  // 5. Commit: flatten per-surface buckets and per-source timing into the
  //    frame's contiguous storage in a single atomic assignment.
  std::vector<SurfaceMeasurement> flatMeasurements;
  std::vector<SurfaceMeasurementRange> ranges(catalog.nSurfaces);
  size_t total = 0;
  for (uint32_t s = 0; s < catalog.nSurfaces; ++s) {
    total += perSurface[s].size();
  }
  flatMeasurements.reserve(total);
  for (uint32_t s = 0; s < catalog.nSurfaces; ++s) {
    ranges[s] = SurfaceMeasurementRange{static_cast<uint32_t>(flatMeasurements.size()), static_cast<uint32_t>(perSurface[s].size())};
    flatMeasurements.insert(flatMeasurements.end(), perSurface[s].begin(), perSurface[s].end());
  }

  std::vector<uint32_t> sourceOffsets(nSources + 1, 0);
  for (uint32_t s = 0; s < nSources; ++s) {
    sourceOffsets[s + 1] = sourceOffsets[s] + static_cast<uint32_t>(perSourceIntervals[s].size());
  }
  std::vector<ROFIntervalBC> flatIntervals;
  flatIntervals.reserve(sourceOffsets.back());
  for (uint32_t s = 0; s < nSources; ++s) {
    flatIntervals.insert(flatIntervals.end(), perSourceIntervals[s].begin(), perSourceIntervals[s].end());
  }

  frame.assignLoadedData(std::move(flatMeasurements), std::move(ranges), std::move(sourcesMeta),
                         std::move(flatIntervals), std::move(sourceOffsets), std::move(labelSources));
  return {};
}

} // namespace o2::itsmft::tracking

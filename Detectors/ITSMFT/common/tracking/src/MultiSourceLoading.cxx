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
} // namespace

LoadSourcesResult loadSources(MultiSourceFrame& frame,
                              const DetectorLayoutView& layout,
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
  }

  // Staging storage: nothing here is visible to `frame` until every source
  // has been fully validated and decoded.
  std::vector<std::vector<SurfaceMeasurement>> perSurface(layout.nSurfaces);
  std::vector<SourceMetadata> sourcesMeta(nSources);
  std::vector<std::vector<ROFIntervalBC>> perSourceIntervals(nSources);
  std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*> labelSources(nSources, nullptr);

  for (const auto& src : sources) {
    const auto srcIdx = src.id.value();
    sourcesMeta[srcIdx] = SourceMetadata{src.id, src.detector, static_cast<uint32_t>(src.rofs.size())};
    labelSources[srcIdx] = src.labels;

    // 2. Validate ROF cluster ranges against this source's cluster span
    //    before any decoding, so an out-of-range access can never happen.
    int64_t previousEnd = 0;
    for (uint32_t r = 0; r < src.rofs.size(); ++r) {
      const auto& rof = src.rofs[r];
      const int64_t first = rof.getFirstEntry();
      const int64_t n = rof.getNEntries();
      if (first < 0 || n < 0 || first + n > static_cast<int64_t>(src.clusters.size()) || first < previousEnd) {
        return {MultiSourceLoadError::InvalidROFRange, src.id, r};
      }
      previousEnd = first + n;
    }

    // 3. Timing: one independent config per source; sourceROF is the
    //    source-local ordinal position of the ROF within `src.rofs`.
    auto& intervals = perSourceIntervals[srcIdx];
    intervals.reserve(src.rofs.size());
    for (uint32_t r = 0; r < src.rofs.size(); ++r) {
      const auto built = computeROFIntervalBC(src.rofs[r].getBCData(), origin, src.timing, r);
      if (!built.ok()) {
        return {MultiSourceLoadError::TimingError, src.id, r};
      }
      intervals.push_back(built.interval);
    }

    // 4. Decode: one independent pattern cursor per source. Each cluster is
    //    decoded exactly once by the source's decoder.
    src.decoder->prepare();
    auto pattIt = src.patterns.begin();
    for (uint32_t r = 0; r < src.rofs.size(); ++r) {
      const auto& rof = src.rofs[r];
      const auto firstEntry = rof.getFirstEntry();
      const auto nEntries = rof.getNEntries();
      for (int32_t clusterId = firstEntry; clusterId < firstEntry + nEntries; ++clusterId) {
        const auto& cluster = src.clusters[clusterId];
        const auto decoded = src.decoder->decode(cluster, pattIt, src.dictionary, src.layerToSurface,
                                                  src.id, static_cast<uint32_t>(clusterId), r, src.applySysErrors);
        if (!decoded.layerMapped) {
          return {MultiSourceLoadError::InvalidLayerMapping, src.id, r, static_cast<uint32_t>(clusterId)};
        }
        const auto surface = decoded.measurement.surface;
        if (!surface.isValid() || surface.value() >= layout.nSurfaces) {
          return {MultiSourceLoadError::InvalidLayerMapping, src.id, r, static_cast<uint32_t>(clusterId)};
        }
        if (layout.getSurface(surface).detectorId != static_cast<uint8_t>(src.detector)) {
          return {MultiSourceLoadError::DetectorSurfaceMismatch, src.id, r, static_cast<uint32_t>(clusterId)};
        }
        perSurface[surface.value()].push_back(decoded.measurement);
      }
    }
  }

  // 5. Commit: flatten per-surface buckets and per-source timing into the
  //    frame's contiguous storage in a single atomic assignment.
  std::vector<SurfaceMeasurement> flatMeasurements;
  std::vector<SurfaceMeasurementRange> ranges(layout.nSurfaces);
  size_t total = 0;
  for (uint32_t s = 0; s < layout.nSurfaces; ++s) {
    total += perSurface[s].size();
  }
  flatMeasurements.reserve(total);
  for (uint32_t s = 0; s < layout.nSurfaces; ++s) {
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

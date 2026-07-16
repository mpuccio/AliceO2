// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_MULTISOURCEFRAME_H_
#define ALICEO2_ITSMFT_TRACKING_MULTISOURCEFRAME_H_

#include <cstdint>
#include <type_traits>

#include "GPUCommonDef.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceTiming.h"

#ifndef GPUCA_GPUCODE
#include <vector>

#include <gsl/gsl>

#include "DetectorsCommonDataFormats/DetID.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#endif

namespace o2::itsmft::tracking
{

struct SurfaceMeasurementRange {
  uint32_t firstEntry{0};
  uint32_t entries{0};

  GPUhdi() constexpr uint32_t getFirstEntry() const noexcept { return firstEntry; }
  GPUhdi() constexpr uint32_t getEntries() const noexcept { return entries; }
};

static_assert(std::is_standard_layout_v<SurfaceMeasurementRange>);
static_assert(std::is_trivially_copyable_v<SurfaceMeasurementRange>);

// Read-only, device-facing view: pointer/count pairs only, no STL containers
// or gsl::span (Architecture.md section 6/14).
struct MultiSourceFrameView {
  const SurfaceMeasurement* measurements{nullptr};
  uint32_t nMeasurements{0};
  const SurfaceMeasurementRange* surfaceRanges{nullptr}; // indexed by SurfaceId, size == nSurfaces
  uint32_t nSurfaces{0};
  const ROFIntervalBC* rofIntervals{nullptr}; // flattened per source, in sourceROF order
  const uint32_t* sourceROFOffsets{nullptr};  // size == nSources + 1
  uint32_t nSources{0};

  GPUhdi() const SurfaceMeasurement* getSurfaceMeasurements(SurfaceId surface) const noexcept
  {
    return measurements + surfaceRanges[surface.value()].firstEntry;
  }
  GPUhdi() uint32_t getSurfaceMeasurementCount(SurfaceId surface) const noexcept
  {
    return surfaceRanges[surface.value()].entries;
  }
  GPUhdi() const ROFIntervalBC* getSourceIntervals(ClusterSourceId source) const noexcept
  {
    return rofIntervals + sourceROFOffsets[source.value()];
  }
  GPUhdi() uint32_t getSourceIntervalCount(ClusterSourceId source) const noexcept
  {
    return sourceROFOffsets[source.value() + 1] - sourceROFOffsets[source.value()];
  }
};

static_assert(std::is_standard_layout_v<MultiSourceFrameView>);
static_assert(std::is_trivially_copyable_v<MultiSourceFrameView>);

#ifndef GPUCA_GPUCODE

struct SourceMetadata {
  ClusterSourceId id{};
  o2::detectors::DetID::ID detector{o2::detectors::DetID::ITS};
  uint32_t nROFs{0};
};

// Standalone host owner (Architecture.md section 7.1): normalized
// measurements partitioned by global SurfaceId, source metadata, source
// timing intervals, and label-lookup metadata. This is not the tracking
// TimeFrame: no CA artefacts, sorting, index tables, cells, roads or tracks
// are stored here. It retains no geometry singletons, dictionaries, compact
// clusters or detector output types after loading.
class MultiSourceFrame
{
 public:
  MultiSourceFrame() = default;

  void clear() noexcept;

  MultiSourceFrameView getView() const noexcept;

  gsl::span<const SurfaceMeasurement> getSurfaceMeasurements(SurfaceId surface) const;
  gsl::span<const ROFIntervalBC> getSourceIntervals(ClusterSourceId source) const;
  // MC labels stay external; resolved through ClusterRef. Absent labels
  // (no container for that source, or an out-of-range index) are legal and
  // yield an empty span rather than an error.
  gsl::span<const o2::MCCompLabel> getLabels(ClusterRef cluster) const;

  const std::vector<SourceMetadata>& getSources() const noexcept { return mSources; }
  uint32_t getNSurfaces() const noexcept { return static_cast<uint32_t>(mSurfaceRanges.size()); }
  size_t getTotalMeasurements() const noexcept { return mMeasurements.size(); }

  // Assembled exclusively by loadSources() once every source has been fully
  // validated and decoded; replaces all prior content atomically so that a
  // failed load (which never calls this) leaves the frame untouched.
  void assignLoadedData(std::vector<SurfaceMeasurement>&& measurements,
                        std::vector<SurfaceMeasurementRange>&& surfaceRanges,
                        std::vector<SourceMetadata>&& sources,
                        std::vector<ROFIntervalBC>&& rofIntervals,
                        std::vector<uint32_t>&& sourceROFOffsets,
                        std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>&& labelSources);

 private:
  std::vector<SurfaceMeasurement> mMeasurements;
  std::vector<SurfaceMeasurementRange> mSurfaceRanges;
  std::vector<SourceMetadata> mSources;
  std::vector<ROFIntervalBC> mROFIntervals;
  std::vector<uint32_t> mSourceROFOffsets;
  std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*> mLabelSources;
};

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_MULTISOURCEFRAME_H_ */

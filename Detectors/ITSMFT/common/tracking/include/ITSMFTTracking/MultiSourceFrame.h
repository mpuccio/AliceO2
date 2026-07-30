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
#include "ITSMFTTracking/SurfaceMeasurementIndex.h"
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
// or gsl::span (Architecture.md section 6/14). Every pointer here is
// non-owning and remains valid only while the MultiSourceFrame that produced
// it is alive and has not been cleared, reloaded, moved from, or destroyed;
// this view does not extend that owner's lifetime.
struct MultiSourceFrameView {
  const SurfaceMeasurement* measurements{nullptr};
  uint32_t nMeasurements{0};
  const SurfaceMeasurementRange* surfaceRanges{nullptr}; // indexed by SurfaceId, size == nSurfaces
  uint32_t nSurfaces{0};
  const ROFIntervalBC* rofIntervals{nullptr}; // flattened per source, in sourceROF order
  const uint32_t* sourceROFOffsets{nullptr};  // size == nSources + 1
  uint32_t nSources{0};

  // Returns nullptr (never a computed nullptr+0) when the surface has no
  // measurements.
  GPUhdi() const SurfaceMeasurement* getSurfaceMeasurements(SurfaceId surface) const noexcept
  {
    const auto& range = surfaceRanges[surface.value()];
    return range.entries == 0 ? nullptr : measurements + range.firstEntry;
  }
  GPUhdi() uint32_t getSurfaceMeasurementCount(SurfaceId surface) const noexcept
  {
    return surfaceRanges[surface.value()].entries;
  }
  // Bounds-checked lookup of one measurement by its canonical, flattened
  // position (SurfaceMeasurementIndex is a position in `measurements`, the
  // single flat array every surface's span is carved from -- not a
  // per-surface-local position). Returns nullptr for an invalid index or one
  // at/past `nMeasurements`. This is the minimal support CommonTrack's
  // trackClusterIndices (ITSMFTTracking/CommonTrack.h) needs to validate and
  // dereference a stored SurfaceMeasurementIndex.
  GPUhdi() const SurfaceMeasurement* getMeasurement(SurfaceMeasurementIndex index) const noexcept
  {
    return (index.isValid() && index.value() < nMeasurements) ? measurements + index.value() : nullptr;
  }
  // Returns nullptr (never a computed nullptr+0) when the source has no
  // ROF intervals.
  GPUhdi() const ROFIntervalBC* getSourceIntervals(ClusterSourceId source) const noexcept
  {
    const auto first = sourceROFOffsets[source.value()];
    const auto last = sourceROFOffsets[source.value() + 1];
    return first == last ? nullptr : rofIntervals + first;
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
//
// Non-owning label lifetime: MC labels are never copied into this owner.
// `assignLoadedData()` (called only by loadSources()) stores raw, non-owning
// MCTruthContainer pointers taken directly from each source's
// ClusterSourceInput::labels. Every such pointer -- and therefore every
// span returned by getLabels() -- is valid only as long as the caller's
// label container outlives it. clear(), a subsequent successful
// loadSources() call, moving from or destroying this owner, or destroying a
// source's label container all invalidate the corresponding lookups.
// Likewise, a MultiSourceFrameView obtained from getView() only remains
// valid while this owner is alive and has not been cleared or reloaded.
class MultiSourceFrame
{
 public:
  MultiSourceFrame() = default;

  void clear() noexcept;

  MultiSourceFrameView getView() const noexcept;

  gsl::span<const SurfaceMeasurement> getSurfaceMeasurements(SurfaceId surface) const;
  // Bounds-checked lookup by canonical flattened position -- see
  // MultiSourceFrameView::getMeasurement()'s doc above; same contract,
  // host-owner-side accessor.
  const SurfaceMeasurement* getMeasurement(SurfaceMeasurementIndex index) const noexcept;
  gsl::span<const ROFIntervalBC> getSourceIntervals(ClusterSourceId source) const;
  // MC labels stay external -- never copied here -- and are resolved
  // through ClusterRef via the non-owning per-source container pointer
  // recorded by assignLoadedData(). Absent labels (no container for that
  // source, or an out-of-range index) are legal and yield an empty span
  // rather than an error. The returned span is only valid while the
  // referenced label container is still alive.
  gsl::span<const o2::MCCompLabel> getLabels(ClusterRef cluster) const;

  const std::vector<SourceMetadata>& getSources() const noexcept { return mSources; }
  uint32_t getNSurfaces() const noexcept { return static_cast<uint32_t>(mSurfaceRanges.size()); }
  size_t getTotalMeasurements() const noexcept { return mMeasurements.size(); }

  // Assembled exclusively by loadSources() once every source has been fully
  // validated and decoded; replaces all prior content atomically so that a
  // failed load (which never calls this) leaves the frame untouched.
  // `labelSources` are stored as-is (non-owning): see the class-level note
  // on label lifetime above.
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
  // Non-owning: see the class-level note on label lifetime above.
  std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*> mLabelSources;
};

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_MULTISOURCEFRAME_H_ */

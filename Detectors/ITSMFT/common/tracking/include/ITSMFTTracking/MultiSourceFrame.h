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
#include "ITSMFTTracking/SurfaceId.h"
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

// Non-owning pointer/count pair into one surface's own measurement array.
// Measurements are owned per surface (one array per SurfaceId), never as a
// single flattened array with per-surface offsets, so this is a direct
// pointer into that surface's own storage -- never an offset into a shared
// array. Indexed by SurfaceId in MultiSourceFrameView::surfaces.
struct SurfaceMeasurementSpan {
  const SurfaceMeasurement* data{nullptr};
  uint32_t count{0};
};

static_assert(std::is_standard_layout_v<SurfaceMeasurementSpan>);
static_assert(std::is_trivially_copyable_v<SurfaceMeasurementSpan>);

// Read-only, device-facing view: pointer/count pairs only, no STL containers
// or gsl::span (Architecture.md section 6/14). Every pointer here is
// non-owning and remains valid only while the MultiSourceFrame that produced
// it is alive and has not been cleared, reloaded, moved from, or destroyed;
// this view does not extend that owner's lifetime.
//
// `surfaces` is a POD array of per-surface pointer/count spans, indexed
// directly by SurfaceId -- there is no flattened/global measurement index
// domain anywhere in this view or its owner: a SurfaceMeasurementIndex is
// always a position within one specific surface's own array, never a
// position in some larger shared array.
struct MultiSourceFrameView {
  const SurfaceMeasurementSpan* surfaces{nullptr}; // indexed by SurfaceId, size == nSurfaces
  uint32_t nSurfaces{0};
  const ROFIntervalBC* rofIntervals{nullptr}; // flattened per source, in sourceROF order
  const uint32_t* sourceROFOffsets{nullptr};  // size == nSources + 1
  uint32_t nSources{0};

  // Bounds-unchecked (matching the per-surface accessor convention every
  // other device view in this library uses: callers must already know
  // `surface` is valid for this
  // view. getMeasurement() below is the bounds-checked counterpart. Returns
  // nullptr (never a computed nullptr+0) when the surface has no
  // measurements.
  GPUhdi() const SurfaceMeasurement* getSurfaceMeasurements(SurfaceId surface) const noexcept
  {
    const auto& span = surfaces[surface.value()];
    return span.count == 0 ? nullptr : span.data;
  }
  GPUhdi() uint32_t getSurfaceMeasurementCount(SurfaceId surface) const noexcept
  {
    return surfaces[surface.value()].count;
  }
  // Bounds-checked lookup of one measurement by the surface it belongs to
  // and its position within that surface's own measurement array (a
  // SurfaceMeasurementIndex is always surface-local; it is never valid to
  // resolve one against a different surface than the one it was paired
  // with, e.g. via TrackClusterReference, ITSMFTTracking/CommonTrack.h).
  // Returns nullptr for an invalid/out-of-range surface or an
  // invalid/out-of-range index.
  GPUhdi() const SurfaceMeasurement* getMeasurement(SurfaceId surface, SurfaceMeasurementIndex index) const noexcept
  {
    if (!surface.isValid() || surface.value() >= nSurfaces || !index.isValid()) {
      return nullptr;
    }
    const auto& span = surfaces[surface.value()];
    return index.value() < span.count ? span.data + index.value() : nullptr;
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
// measurements owned per surface (one array per SurfaceId, never a single
// flattened array with per-surface offsets), source metadata, source timing
// intervals, and label-lookup metadata. This is not the tracking TimeFrame:
// no CA artefacts, sorting, index tables, cells, roads or tracks are stored
// here. It retains no geometry singletons, dictionaries, compact clusters or
// detector output types after loading.
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

  // mSurfaceSpans caches raw pointers into mPerSurfaceMeasurements. A
  // member-wise (default) copy would copy those pointers as-is, leaving the
  // copy's mSurfaceSpans pointing into the *original* object's storage
  // instead of its own already-duplicated one -- silently wrong, not merely
  // unsafe, immediately after the copy completes. Correctly supporting copy
  // would require rebuilding mSurfaceSpans from the copy's own
  // mPerSurfaceMeasurements after duplicating it; no caller of this type
  // needs that today (TimeFrame owns exactly one MultiSourceFrame per
  // normalized load and only ever moves it, never copies it), so copy is
  // deleted rather than implemented incorrectly by default.
  MultiSourceFrame(const MultiSourceFrame&) = delete;
  MultiSourceFrame& operator=(const MultiSourceFrame&) = delete;

  // Move is implemented via swap() rather than relying on the member-wise
  // move std::vector would otherwise generate implicitly. Both are
  // semantically equivalent here (moving a vector-of-vectors never touches
  // the innermost heap buffers mSurfaceSpans' cached pointers reference, so
  // plain member-wise move would also leave every span correctly pointing
  // into this object's own, now-moved-in storage) -- but swap() makes the
  // noexcept guarantee provable from this operator's own definition
  // (std::vector::swap is unconditionally noexcept: a pointer/size/capacity
  // exchange only, never a reallocation or per-element operation) rather
  // than from a conditional trait computed over every member's move
  // operation, and it leaves `other` holding a fully self-consistent state
  // (`*this`'s old content, spans included) rather than an unspecified-but-
  // valid moved-from state. This is what guarantees a successful
  // TimeFrame::loadNormalizedSource() commit (`mNormalizedFrame =
  // std::move(staged)`) can never leave mNormalizedFrame's cached spans
  // pointing into `staged` (or vice versa).
  MultiSourceFrame(MultiSourceFrame&& other) noexcept { swap(other); }
  MultiSourceFrame& operator=(MultiSourceFrame&& other) noexcept
  {
    swap(other);
    return *this;
  }

  // Exchanges all owned storage, including the mSurfaceSpans cache, in
  // lockstep -- every member swap below is a std::vector::swap, so this
  // can never throw and can never leave either side's mSurfaceSpans
  // pointing into the other side's mPerSurfaceMeasurements.
  void swap(MultiSourceFrame& other) noexcept
  {
    mPerSurfaceMeasurements.swap(other.mPerSurfaceMeasurements);
    mSurfaceSpans.swap(other.mSurfaceSpans);
    mSources.swap(other.mSources);
    mROFIntervals.swap(other.mROFIntervals);
    mSourceROFOffsets.swap(other.mSourceROFOffsets);
    mLabelSources.swap(other.mLabelSources);
  }

  void clear() noexcept;

  MultiSourceFrameView getView() const noexcept;

  gsl::span<const SurfaceMeasurement> getSurfaceMeasurements(SurfaceId surface) const;
  // Bounds-checked lookup by the surface it belongs to and its position
  // within that surface's own measurement array -- see
  // MultiSourceFrameView::getMeasurement()'s doc above; same contract,
  // host-owner-side accessor.
  const SurfaceMeasurement* getMeasurement(SurfaceId surface, SurfaceMeasurementIndex index) const noexcept;
  gsl::span<const ROFIntervalBC> getSourceIntervals(ClusterSourceId source) const;
  // MC labels stay external -- never copied here -- and are resolved
  // through ClusterRef via the non-owning per-source container pointer
  // recorded by assignLoadedData(). Absent labels (no container for that
  // source, or an out-of-range index) are legal and yield an empty span
  // rather than an error. The returned span is only valid while the
  // referenced label container is still alive.
  gsl::span<const o2::MCCompLabel> getLabels(ClusterRef cluster) const;

  const std::vector<SourceMetadata>& getSources() const noexcept { return mSources; }
  uint32_t getNSurfaces() const noexcept { return static_cast<uint32_t>(mPerSurfaceMeasurements.size()); }
  size_t getTotalMeasurements() const noexcept
  {
    size_t total = 0;
    for (const auto& measurements : mPerSurfaceMeasurements) {
      total += measurements.size();
    }
    return total;
  }

  // Assembled exclusively by loadSources() once every source has been fully
  // validated and decoded; replaces all prior content atomically so that a
  // failed load (which never calls this) leaves the frame untouched.
  // `perSurfaceMeasurements` is indexed by SurfaceId (size == the resolved
  // catalog's surface count); each entry is that surface's own, already
  // complete measurement array -- never flattened into, or reconstructed
  // from, a shared array. `labelSources` are stored as-is (non-owning): see
  // the class-level note on label lifetime above.
  void assignLoadedData(std::vector<std::vector<SurfaceMeasurement>>&& perSurfaceMeasurements,
                        std::vector<SourceMetadata>&& sources,
                        std::vector<ROFIntervalBC>&& rofIntervals,
                        std::vector<uint32_t>&& sourceROFOffsets,
                        std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>&& labelSources);

 private:
  void rebuildSurfaceSpans();

  std::vector<std::vector<SurfaceMeasurement>> mPerSurfaceMeasurements;
  // Cached pointer/count view over mPerSurfaceMeasurements, rebuilt whenever
  // that storage changes (assignLoadedData()/clear()); this is what getView()
  // hands out and is never itself the owner of any measurement.
  std::vector<SurfaceMeasurementSpan> mSurfaceSpans;
  std::vector<SourceMetadata> mSources;
  std::vector<ROFIntervalBC> mROFIntervals;
  std::vector<uint32_t> mSourceROFOffsets;
  // Non-owning: see the class-level note on label lifetime above.
  std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*> mLabelSources;
};

static_assert(std::is_nothrow_move_constructible_v<MultiSourceFrame>);
static_assert(std::is_nothrow_move_assignable_v<MultiSourceFrame>);
static_assert(!std::is_copy_constructible_v<MultiSourceFrame>);
static_assert(!std::is_copy_assignable_v<MultiSourceFrame>);

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_MULTISOURCEFRAME_H_ */

// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TrackRefit.h
/// \brief The refit of the selected tracks into the output record.
///
/// The candidate fit kept a road's chi2 and its outward pass's end state; a
/// selected track gets the fit's inward pass (`fitInward`) on its rows, and with it
/// its state and covariance on its innermost surface, the outward state taken
/// from the record. One element per
/// selected track, the output one record per element: nothing to scan.

#ifndef GTRK_TRACKING_TRACKREFIT_H_
#define GTRK_TRACKING_TRACKREFIT_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/Columns.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/Selection.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackFitter.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include <memlayout/wrapper.h>

#include <cstddef>
#include <cstdint>
#include <memory>

namespace o2::gtrk
{

/// A track as the iteration outputs it: the inward pass's state on its
/// innermost surface, the outward pass's state on its outermost surface, its
/// time box `[timePS, timePS + widthPS)`, its
/// fit quality, and its measurement rows, innermost first. \a K is the
/// detector's longest chain.
template <std::size_t K>
struct TrackOf {
  // As for TrackletT: upstream's `apply` macro also uses the member names as
  // enum constants.
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
  template <template <class> class C>
  struct RecordT {
    MEMLAYOUT_APPLY(RecordT, state, outwardState, timePS, widthPS, chi2, outwardChi2, ndf, status, clusters, count, iteration, addedOutward, addedInward, kinks, kinkCluster)
    C<BoundTrackState> state;
    C<BoundTrackState> outwardState;
    C<std::uint32_t> timePS;
    C<std::uint32_t> widthPS;
    C<Scalar> chi2;        ///< the inward pass's
    C<Scalar> outwardChi2; ///< the outward pass's, which `outwardState` ends
    C<std::int16_t> ndf;
    C<PropagationStatus> status; ///< the refit's: kOk, or why it failed
    C<ClusterList<K>> clusters;
    C<std::uint8_t> count;
    C<std::uint8_t> iteration; ///< the iteration that found it
    C<std::uint8_t> addedOutward;
    C<std::uint8_t> addedInward;
    C<std::uint8_t> kinks;       ///< the inward pass's
    C<std::uint8_t> kinkCluster; ///< its last kink's position in `clusters`; kNoKink without
  };
  static constexpr std::uint8_t kNoKink = 0xFFU;
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
};

/// A track's time box: the intersection of its \a count rows' boxes, not
/// empty -- the roads' time test and the extension keep it so.
template <typename Columns, std::size_t K>
[[nodiscard]] GPUhdi() GPUAlwaysInline() TimeEstPS trackTimeBox(const Columns& columns, const ClusterList<K>& clusters, std::uint8_t count)
{
  auto box = timeBoxAt(columns, clusters.index[0]);
  for (std::uint8_t i = 1; i < count; ++i) {
    box += timeBoxAt(columns, clusters.index[i]);
  }
  return box;
}

/// Write \a fit of the rows \a clusters, \a count of them, as track \a track
/// of \a tracks with the time box \a box (`trackTimeBox`): every column but
/// the iteration and the clusters an extension added, which the caller knows
/// (and the outward chi2 after `fitInward`, which has none).
template <typename Tracks, std::size_t K>
GPUhdi() GPUAlwaysInline() void writeFit(const Tracks& tracks, std::uint32_t track, const TrackFitResult& fit, const ClusterList<K>& clusters, std::uint8_t count, const TimeEstPS& box)
{
  tracks.state[track] = fit.state;
  tracks.outwardState[track] = fit.outwardState;
  tracks.timePS[track] = box.start;
  tracks.widthPS[track] = box.width;
  tracks.chi2[track] = fit.chi2;
  tracks.outwardChi2[track] = fit.outwardChi2;
  tracks.ndf[track] = fit.ndf;
  tracks.status[track] = fit.status;
  tracks.clusters[track] = clusters;
  tracks.count[track] = count;
  tracks.kinks[track] = fit.kinks;
  tracks.kinkCluster[track] = fit.kinkMeasurement < 0 ? TrackOf<K>::kNoKink : static_cast<std::uint8_t>(fit.kinkMeasurement);
}

/// One element per selected record: the inward pass on its rows from the road
/// fit's outward state (`fitInward`), and the track's time from their boxes;
/// the pass may take the iteration's kinks. Bidirectional reference mode starts
/// its inward reference from that outward state. Its frame, geometry and material
/// are pointers, which may be a device's.
template <DetectorLayout Det, typename Material>
class TrackRefitElement
{
 public:
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using View = typename BoundedColumns<TrackOf<kMaxClusters>::template RecordT>::View;

  TrackRefitElement(const TimeFrame<Det>* frame, const TrackingGeometry<Det>* geometry, const Material* material, const Propagator& propagator, const ClusterList<kMaxClusters>* clusters, const std::uint8_t* count, const BoundTrackState* outward, const Scalar* outwardChi2, const std::uint32_t* selected, View tracks, std::uint32_t offset, std::uint8_t iteration, KinkBudget kinks, FitReference reference = FitReference::kFilteredState)
    : mFrame(frame),
      mGeometry(geometry),
      mMaterial(material),
      mPropagator(propagator),
      mClusters(clusters),
      mCount(count),
      mOutward(outward),
      mOutwardChi2(outwardChi2),
      mSelected(selected),
      mTracks(tracks),
      mOffset(offset),
      mIteration(iteration),
      mKinks(kinks),
      mReference(reference)
  {
  }

  GPUhdi() void operator()(std::uint32_t element) const
  {
    const std::uint32_t record = mSelected[element];
    const std::uint32_t track = mOffset + element;
    const auto& rows = mClusters[record];
    const std::uint8_t count = mCount[record];
    // The road fit ran the outward pass already and kept its end state.
    const auto chain = RowChain<Det>{.frame = mFrame, .rows = rows.index};
    const auto fit = fitInwardReference(mPropagator, *mGeometry, *mMaterial, chain, static_cast<int>(count), mOutward[record], mReference, mKinks);
    writeFit(mTracks, track, fit, rows, count, trackTimeBox(mFrame->measurements(), rows, count));
    mTracks.outwardChi2[track] = mOutwardChi2[record];
    mTracks.iteration[track] = mIteration;
    mTracks.addedOutward[track] = 0U;
    mTracks.addedInward[track] = 0U;
  }

 private:
  const TimeFrame<Det>* mFrame = nullptr;
  const TrackingGeometry<Det>* mGeometry = nullptr;
  const Material* mMaterial = nullptr;
  Propagator mPropagator;
  const ClusterList<kMaxClusters>* mClusters = nullptr;
  const std::uint8_t* mCount = nullptr;
  const BoundTrackState* mOutward = nullptr;
  const Scalar* mOutwardChi2 = nullptr; ///< the road fit's, whose end state `mOutward` is
  const std::uint32_t* mSelected = nullptr;
  View mTracks;
  std::uint32_t mOffset = 0U;
  std::uint8_t mIteration = 0U;
  KinkBudget mKinks;
  FitReference mReference;
};

/// The timeframe's tracks: each iteration's appended after the ones before,
/// one per selected record, in the selection's order, best first.
template <DetectorLayout Det>
class TrackRefit
{
 public:
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using Tracks = BoundedColumns<TrackOf<kMaxClusters>::template RecordT>;

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool.swap(pool);
    clearStorage();
  }

  /// Start a timeframe: no tracks.
  void clear() { mTracks.resize(0U); }

  /// The tracks' storage back to the pool.
  void clearStorage() noexcept { mTracks.setMemoryResource(mMemoryPool.get()); }

  /// Refit \a selection's selected records, iteration \a iteration's, after
  /// the tracks already held. \a forEach(elements, body) runs body(element)
  /// for every element, serially or in parallel.
  template <typename Material, typename ForEach>
  [[nodiscard]] StageStatus run(ForEach forEach, const GlobalSelection<kMaxClusters>& selection, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator, std::uint8_t iteration)
  {
    return guardOutOfMemory([&] {
      const auto selected = selection.selected();
      const auto offset = static_cast<std::uint32_t>(mTracks.size());
      mTracks.resize(offset + selected.size());
      const auto records = selection.records();
      const TrackRefitElement<Det, Material> element{&frame, &geometry, &material, propagator, records.clusters, records.count, records.outwardState, records.chi2, selected.data(), mTracks.view(), offset, iteration, kinkBudget(parameters), parameters.fitReference};
      forEach(static_cast<std::uint32_t>(selected.size()), element);
      return StageStatus::kComplete;
    });
  }

  [[nodiscard]] const Tracks& tracks() const { return mTracks; }
  [[nodiscard]] Tracks& tracks() { return mTracks; }
  [[nodiscard]] std::size_t size() const { return mTracks.size(); }

 private:
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  Tracks mTracks;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRACKREFIT_H_

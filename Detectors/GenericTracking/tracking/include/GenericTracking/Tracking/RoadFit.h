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

/// \file RoadFit.h
/// \brief The road fit's per-road function and the selection record it
///        emits.

#ifndef GTRK_TRACKING_ROADFIT_H_
#define GTRK_TRACKING_ROADFIT_H_

#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/BeamConstraint.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/RoadCount.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackFitter.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/Triplet.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include "GPUCommonMath.h"

#include <memlayout/wrapper.h>

#include <cstddef>
#include <cstdint>
#include <limits>

namespace o2::gtrk
{

/// A road's measurement rows, innermost first; the slots past its count hold
/// `kNoCluster`. \a K is the detector's longest chain. Rows are
/// timeframe-global and unique like cluster indices, and what the used mask
/// indexes; cluster indices are read only for the accepted tracks' output.
template <std::size_t K>
struct ClusterList {
  static constexpr std::uint32_t kNoCluster = std::numeric_limits<std::uint32_t>::max();
  std::uint32_t index[K];
};

/// What the global selection needs of a road: its clusters' rows, in
/// road order -- a function of the set, the link table being acyclic -- and
/// its chi2; and the outward pass's end state, which the refit of a selected
/// road takes instead of running that pass again, and runs only the inward one.
template <std::size_t K>
struct SelectionRecordOf {
  // As for TrackletT: upstream's `apply` macro also uses the member names as
  // enum constants.
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
  template <template <class> class C>
  struct RecordT {
    MEMLAYOUT_APPLY(RecordT, clusters, count, chi2, kinks, outwardState)
    C<ClusterList<K>> clusters;
    C<std::uint8_t> count;           ///< clusters on the road
    C<Scalar> chi2;                  ///< the outward pass's, over all of them
    C<std::uint8_t> kinks;           ///< the outward pass's
    C<BoundTrackState> outwardState; ///< the outward pass's, on the outermost surface
  };
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
};

/// Row \a row of \a frame as a fit reads it: only the nine columns it uses
/// -- a whole row through the view's `operator[]` is a call that reads all
/// nineteen.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() Measurement measurementAt(const TimeFrame<Det>& frame, std::uint32_t row)
{
  const auto columns = frame.measurements();
  return Measurement{.loc0 = columns.loc0[row],
                     .loc1 = columns.loc1[row],
                     .timePS = columns.timePS[row],
                     .widthPS = columns.widthPS[row],
                     .varLoc0 = columns.varLoc0[row],
                     .covLoc01 = columns.covLoc01[row],
                     .varLoc1 = columns.varLoc1[row],
                     .clusterIndex = columns.clusterIndex[row],
                     .surface = columns.surface[row]};
}

/// A road's measurements by their rows in the frame's columns: the fit reads
/// one when it needs it instead of copying the road's worth into the element.
template <DetectorLayout Det>
struct RowChain {
  const TimeFrame<Det>* frame = nullptr;
  const std::uint32_t* rows = nullptr;

  [[nodiscard]] GPUhdi() GPUAlwaysInline() Measurement operator[](int i) const { return measurementAt(*frame, rows[i]); }
};

/// Whether \a state, a fit of \a count clusters, is at or above the
/// iteration's pT floor for that count.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool clearsPtFloor(const IterationParameters<Det>& parameters, std::uint32_t count, const BoundTrackState& state)
{
  const Scalar floor = parameters.minPtByClusters[count];
  return !(floor > 0.F) || o2::gpu::CAMath::Abs(o2::gpu::CAMath::Sin(state.theta())) >= floor * o2::gpu::CAMath::Abs(state.qOverP());
}

/// Whether a road of \a count clusters the fit accepted clears the
/// iteration's pT floor, at its outermost measurement. The floors are the
/// final track's: a road is one only where the
/// iteration does not follow; where it does, they act after the follow
/// (`keptAfterFollow`, Follow.h) and every road clears them here. O2 applies its
/// floor at the innermost, after the inward pass; between the two only the
/// energy lost in the road's own layers differs.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool clearsPtFloor(const IterationParameters<Det>& parameters, std::uint32_t count, const RoadFitResult& fit)
{
  return parameters.follow || clearsPtFloor(parameters, count, fit.state);
}

/// The beam-line test of a road the fit accepted: its inward pass, as the refit runs it from \a outwardState, and
/// the update of the inward state with the beam line (`updateWithBeamLine`);
/// with a kink allowance, again with the scattering inside scaled by it if
/// the nominal update is above the gate. False if either fails or the update
/// is above the iteration's gate; \a chi2 is the update's.
template <DetectorLayout Det, typename Material>
[[nodiscard]] GPUdi() bool clearsBeamLine(const Propagator& propagator, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame,
                                          const TrackingGeometry<Det>& geometry, const Material& material, const std::uint32_t* rows, std::uint32_t count,
                                          const BoundTrackState& outwardState, Scalar& chi2)
{
  const auto inward =
    fitInward(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = rows}, static_cast<int>(count), outwardState, kinkBudget(parameters));
  if (inward.status != PropagationStatus::kOk) {
    return false;
  }
  BoundTrackState state = inward.state;
  const auto& plane = geometry.planeOf(state);
  const LayerID layer = geometry.surface(state.surface).layer;
  const auto point = pairPointAt(frame.measurements(), rows[0]);
  auto update = updateWithBeamLine(propagator, parameters, plane, layer, point, state);
  if (update.ok && !(update.chi2 <= parameters.beamChi2Max) && parameters.kinkScale2 > 1.F) {
    // The kink allowance: the scattering inside scaled by it.
    state = inward.state;
    update = updateWithBeamLine(propagator, parameters, plane, layer, point, state, parameters.kinkScale2);
  }
  chi2 = update.chi2;
  return update.ok && update.chi2 <= parameters.beamChi2Max;
}

/// One element per road, by its index: walk to it, drop it if its boxes share
/// no time (with widened boxes a road can drift a frame per neighbour), fit its measurements with
/// the outward pass, and emit its selection record if the chi2 is within the
/// iteration's bound for its cluster count and its pT at or above its floor. Everything in the element has a
/// size fixed by `maxChainLayers<Det>`. With \a BeamLine, for an iteration
/// with `beamConstraint`, a road must also clear the beam line
/// (`clearsBeamLine`), and its record's chi2 takes the update's. With
/// \a NominalOnly the stage takes the roads of nominal
/// triplets only: of a graph built with the next iteration's kink allowance,
/// the roads its own iteration would have made (one graph for both rounds;
/// the next iteration fits every road). With
/// \a LastStepAllowance, for an iteration with a kink allowance, the fit's
/// last step takes it (`fitOutward`, against the road-end fakes): the host and the device trackers instantiate it for such an
/// iteration. With \a InnerTrim, for
/// an iteration with `innerTrim`, a road the stage refuses -- no common
/// time, the fit, the floor, the beam line -- is tried again without its
/// innermost cluster, then without two, down to `innerTrim` of them and the
/// iteration's `minClusters`, and the first the stage takes is the road's
/// record (a decay daughter's chain walked on below it into a triplet of
/// another's cluster).
template <DetectorLayout Det, typename Material, bool BeamLine, bool NominalOnly, bool LastStepAllowance, bool InnerTrim>
class RoadFitStage
{
 public:
  using Index = StageOutput<TripletT>::Index;
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using Record = memlayout::Wrapper<SelectionRecordOf<kMaxClusters>::template RecordT, memlayout::value>;

  RoadFitStage(const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material,
               const Propagator& propagator, const StageOutput<TrackletT>& tracklets, const StageOutput<TripletT>& triplets, const RoadWalk& walk)
    : RoadFitStage(&parameters, &frame, &geometry, &material, propagator, tracklets.view(), triplets.view(), walk)
  {
  }

  /// The same stage over inputs wherever they are -- a device's copies:
  /// every argument is only stored.
  RoadFitStage(const IterationParameters<Det>* parameters, const TimeFrame<Det>* frame, const TrackingGeometry<Det>* geometry, const Material* material,
               const Propagator& propagator, TrackletConstView tracklets, TripletConstView triplets, const RoadWalk& walk)
    : mParameters(parameters),
      mFrame(frame),
      mGeometry(geometry),
      mMaterial(material),
      mPropagator(propagator),
      mTracklets(tracklets),
      mTriplets(triplets),
      mWalk(walk)
  {
  }

  template <typename Emitter>
  GPUhdi() void operator()(Index road, Emitter& emitter) const
  {
    Index triplets[kMaxClusters - 2U];
    const std::uint32_t length = mWalk(road, triplets);
    if constexpr (NominalOnly) {
      for (std::uint32_t i = 0; i < length; ++i) {
        if (mTriplets.kinked[triplets[i]] != 0U) {
          return;
        }
      }
    }
    std::uint32_t rows[kMaxClusters];
    const Index first = mTriplets.first[triplets[0]];
    rows[0] = mTracklets.source[first];
    rows[1] = mTracklets.target[first];
    for (std::uint32_t i = 0; i < length; ++i) {
      rows[i + 2U] = mTracklets.target[mTriplets.second[triplets[i]]];
    }
    const std::uint32_t count = length + 2U;
    if constexpr (InnerTrim) {
      const std::uint32_t shortest = mParameters->minClusters;
      const std::uint32_t trims = count > shortest ? o2::gpu::CAMath::Min(std::uint32_t{mParameters->innerTrim}, count - shortest) : 0U;
      for (std::uint32_t dropped = 0; dropped <= trims; ++dropped) {
        if (dropped != 0U) {
          countWork(WorkCounter::kRoadTrimFits);
        }
        if (fitRows(rows + dropped, count - dropped, emitter)) {
          if (dropped != 0U) {
            countWork(WorkCounter::kRoadTrimAccepted);
          }
          return;
        }
      }
    } else {
      fitRows(rows, count, emitter);
    }
  }

 private:
  /// The road of \a count rows \a rows, innermost first: dropped if its
  /// boxes share no time, else fitted with the outward pass, tested against
  /// its floor and, with \a BeamLine, the beam line, and its record emitted.
  /// Whether it was.
  template <typename Emitter>
  GPUhdi() GPUAlwaysInline() bool fitRows(const std::uint32_t* rows, std::uint32_t count, Emitter& emitter) const
  {
    const auto columns = mFrame->measurements();
    auto common = timeBoxAt(columns, rows[0]);
    for (std::uint32_t i = 1; i < count; ++i) {
      const auto box = timeBoxAt(columns, rows[i]);
      if (!common.isCompatible(box)) {
        countWork(WorkCounter::kRoadTime);
        return false;
      }
      common += box;
    }
    const auto fit =
      fitRoad<LastStepAllowance>(*mParameters, mPropagator, *mGeometry, *mMaterial, RowChain<Det>{.frame = mFrame, .rows = rows}, static_cast<int>(count));
    countWork(WorkCounter::kRoadFits);
    if (fit.status != RoadFitStatus::kAccepted) {
      return false;
    }
    if (!clearsPtFloor(*mParameters, count, fit)) {
      countWork(WorkCounter::kRoadFloor);
      return false;
    }
    if constexpr (BeamLine) {
      Scalar beamChi2 = 0.F;
      if (!clearsBeamLine(mPropagator, *mParameters, *mFrame, *mGeometry, *mMaterial, rows, count, fit.state, beamChi2)) {
        countWork(WorkCounter::kRoadBeam);
        return false;
      }
      emitRecord(emitter, rows, count, fit, fit.chi2 + beamChi2);
    } else {
      emitRecord(emitter, rows, count, fit, fit.chi2);
    }
    return true;
  }

  /// The accepted road's record, of chi2 \a chi2.
  template <typename Emitter>
  GPUhdi() GPUAlwaysInline() static void emitRecord(Emitter& emitter, const std::uint32_t* rows, std::uint32_t count, const RoadFitResult& fit, Scalar chi2)
  {
    countWork(WorkCounter::kRoadAccepted);
    Record record;
    for (std::uint32_t i = 0; i < kMaxClusters; ++i) {
      record.clusters.index[i] = i < count ? rows[i] : ClusterList<kMaxClusters>::kNoCluster;
    }
    record.count = static_cast<std::uint8_t>(count);
    record.chi2 = chi2;
    record.kinks = fit.kinks;
    record.outwardState = fit.state;
    emitter.emit(record);
  }

  const IterationParameters<Det>* mParameters = nullptr;
  const TimeFrame<Det>* mFrame = nullptr;
  const TrackingGeometry<Det>* mGeometry = nullptr;
  const Material* mMaterial = nullptr;
  Propagator mPropagator;
  TrackletConstView mTracklets;
  TripletConstView mTriplets;
  RoadWalk mWalk;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_ROADFIT_H_

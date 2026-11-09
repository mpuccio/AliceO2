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

/// \file VertexingPass.h
/// \brief The own vertexer before the tracking, on the host.
///
/// The vertexing pass is a tracker of one iteration at the beam line whose
/// links join the innermost layers: it runs over the timeframe in the
/// tracking's slices and stops at their triplets; the vertexer turns each
/// slice's into the vertices of the frames its core owns (`Vertexer::find`,
/// `OwnedFrames`), so that the slices' lists make the timeframe's, each
/// vertex once, and the pass's memory is a slice's. The list
/// takes the place of the MC's for the iterations that seed from it. The
/// pass keeps its own used mask, all free: it takes no clusters from the
/// tracking.

#ifndef GTRK_TRACKING_VERTEXINGPASS_H_
#define GTRK_TRACKING_VERTEXINGPASS_H_

#include "GenericTracking/Tracking/Tracker.h"
#include "GenericTracking/Tracking/VertexStage.h"

#include <memory>
#include <utility>

namespace o2::gtrk
{

template <DetectorLayout Det, typename Material>
class VertexingPass
{
 public:
  using Pass = Tracker<Det, Material>;
  using Frame = typename Pass::Frame;

  /// What `find` gives: the vertexer's result, or the status the pass
  /// failed with -- out of memory, as any stage can. The pass's tracklets and
  /// triplets stay valid until the next `find`, for the QA.
  struct Result {
    StageStatus status = StageStatus::kComplete;
    VertexerResult found;
    /// With claims asked for, every triplet of the frames owned, its rows and
    /// the vertex in `found` it went to (`VertexerResult::kFree`, `kNone`),
    /// for the QA; the not-owned triplets are left out.
    std::vector<TripletClaim> claims;
  };

  /// \a pass the pass's iteration (GridConfig.h, `vertexingPassSettings`),
  /// \a vertexer the vertexer's settings; the rest as for a Tracker.
  VertexingPass(const IterationSettings& pass, const VertexerSettings& vertexer, const VertexConstraint& beamSpot, const ScatteringTable<Det>& scattering, std::shared_ptr<BoundedMemoryResource> memoryPool, typename Pass::Index chunkSize, int maxConcurrency = 0)
    : mPass(pass, beamSpot, scattering, memoryPool, chunkSize, maxConcurrency),
      mVertexer(vertexer, beamSpot),
      mMemoryPool(std::move(memoryPool))
  {
  }

  /// The pass's memory back to the pool, after the caller has read `find`'s
  /// result.
  void endTimeframe() noexcept { mPass.endTimeframe(); }

  /// The pass's iteration parameters, as a device run configures the pass
  /// after the tracking's iterations.
  [[nodiscard]] const IterationParameters<Det>& passParameters() const { return mPass.iteration(0U); }
  [[nodiscard]] const VertexerSettings& settings() const { return mVertexer.settings(); }

  /// The vertices of \a frame, sorted by time, the pass in slices of
  /// \a sliceWidthBlocks blocks with at least \a haloPS of halo, as
  /// `planTimeframe` plans the tracking's; with \a keepClaims the claims too.
  [[nodiscard]] Result find(const Frame& frame, std::int32_t sliceWidthBlocks, std::int64_t haloPS, bool keepClaims = false)
  {
    Result result;
    result.status = mPass.beginTimeframe(frame);
    if (result.status != StageStatus::kComplete) {
      return result;
    }
    const auto plan = mPass.planTimeframe(frame, {}, sliceWidthBlocks, haloPS);
    // Each candidate's z0 error from the pass's material at its own momentum
    // and the beam's transverse spread.
    const ScatteringErrors<Det> errors{
      .table = &mPass.iteration(0U).scattering,
      .cut = mPass.iteration(0U).pairCut,
      .varTransverse = mPass.beamLine().varTransverse,
    };
    for (const auto& planned : plan.slices) {
      const auto slice = mPass.runSlice(0U, frame, plan, planned, {}, plan.vertexWindows[0]);
      result.status = slice.status;
      if (result.status != StageStatus::kComplete) {
        return result;
      }
      const auto tracklets = slice.tracklets->view();
      const auto triplets = slice.triplets->view();
      const auto found = mVertexer.find(frame.measurements(), tracklets, triplets, slice.triplets->size(), &mPass.arena(), keepClaims, errors, mMemoryPool.get(), ownedFrames(frame.grid().geometry(), planned.slice.coreBlocks));
      const auto offset = static_cast<std::int32_t>(result.found.vertices.size());
      for (std::size_t t = 0; t < found.tripletVertex.size(); ++t) {
        if (const auto vertex = found.tripletVertex[t]; vertex != VertexerResult::kNotOwned) {
          result.claims.push_back({.rows = tripletRows(tracklets, triplets, static_cast<std::uint32_t>(t)), .vertex = vertex < 0 ? vertex : offset + vertex});
        }
      }
      appendVertices(result.found, found);
    }
    return result;
  }

 private:
  Pass mPass;
  Vertexer mVertexer;
  std::shared_ptr<BoundedMemoryResource> mMemoryPool; ///< the pass's, for the vertexer's work arrays
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_VERTEXINGPASS_H_

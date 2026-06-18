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
///
/// \file MFTFwdTrackHelpers.h
/// \brief Forward-track helpers for MFT CA cell fitting (TrackParCovFwd)
///

#ifndef ALICEO2_ITSMFT_TRACKING_MFTFWDTRACKHELPERS_H_
#define ALICEO2_ITSMFT_TRACKING_MFTFWDTRACKHELPERS_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "CommonConstants/MathConstants.h"
#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/CATrackTypes.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/LayerMask.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Constants.h"
#include "MFTTracking/Constants.h"
#include "ReconstructionDataFormats/TrackFwd.h"

namespace o2::itsmft::tracking::detail
{

inline int mftLayerAtSlot(LayerMask mask, int requestedSlot)
{
  int slot = 0;
  for (int layer = mask.first(); layer <= mask.last(); ++layer) {
    if (!mask.has(layer)) {
      continue;
    }
    if (slot++ == requestedSlot) {
      return layer;
    }
  }
  return o2::its::constants::UnusedIndex;
}

namespace mftFwdSeedDefaults
{
constexpr float kMinDr = 1.e-6f;
constexpr float kMinDz = 1.e-6f;
constexpr float kMinBz = 0.01f;
constexpr float kDefaultSigma2 = 1.f;
} // namespace mftFwdSeedDefaults

/// Returns true when three clusters form a usable inward MFT seed (inner z > outer z, non-degenerate dr).
inline bool mftFwdCellSeedGeometryValid(const o2::its::Cluster& cInner,
                                        const o2::its::Cluster& cMid,
                                        const o2::its::Cluster& cOuter)
{
  using namespace mftFwdSeedDefaults;
  if (cInner.zCoordinate <= cOuter.zCoordinate + kMinDz) {
    return false;
  }
  const float dxTan = cMid.xCoordinate - cInner.xCoordinate;
  const float dyTan = cMid.yCoordinate - cInner.yCoordinate;
  const float dzTan = cMid.zCoordinate - cInner.zCoordinate;
  const float drTan = std::sqrt(dxTan * dxTan + dyTan * dyTan);
  if (drTan < kMinDr || std::abs(dzTan) < kMinDz) {
    return false;
  }
  const float dxPhi = cOuter.xCoordinate - cInner.xCoordinate;
  const float dyPhi = cOuter.yCoordinate - cInner.yCoordinate;
  const float dzPhi = cOuter.zCoordinate - cInner.zCoordinate;
  const float drPhi = std::sqrt(dxPhi * dxPhi + dyPhi * dyPhi);
  return drPhi >= kMinDr && std::abs(dzPhi) >= kMinDz;
}

/// Build a TrackParCovFwd seed at the outer cluster, inward-fitting convention (TrackFitter::initTrack).
/// Clusters must be ordered inner → mid → outer (MFT: inner has less negative z).
/// Reference point (x, y, z) is at cOuter; invQPt defaults to 1/minPt when B is on.
inline bool mftBuildFwdTrackSeedInward(const o2::its::Cluster& cInner,
                                       const o2::its::Cluster& cMid,
                                       const o2::its::Cluster& cOuter,
                                       float bz,
                                       float minPt,
                                       o2::track::TrackParCovFwd& track,
                                       float sigma2XOuter = mftFwdSeedDefaults::kDefaultSigma2,
                                       float sigma2YOuter = mftFwdSeedDefaults::kDefaultSigma2)
{
  using namespace mftFwdSeedDefaults;
  if (!mftFwdCellSeedGeometryValid(cInner, cMid, cOuter)) {
    return false;
  }

  const float x0 = cOuter.xCoordinate;
  const float y0 = cOuter.yCoordinate;
  const float z0 = cOuter.zCoordinate;

  const float dxTan = cMid.xCoordinate - cInner.xCoordinate;
  const float dyTan = cMid.yCoordinate - cInner.yCoordinate;
  const float dzTan = cMid.zCoordinate - cInner.zCoordinate;
  const float drTan = std::sqrt(dxTan * dxTan + dyTan * dyTan);

  const float dxPhi = cOuter.xCoordinate - cInner.xCoordinate;
  const float dyPhi = cOuter.yCoordinate - cInner.yCoordinate;
  const float dzPhi = cOuter.zCoordinate - cInner.zCoordinate;
  const float drPhi = std::sqrt(dxPhi * dxPhi + dyPhi * dyPhi);

  const float invQPt = (minPt > 0.f) ? 1.f / minPt : 0.f;

  float tanl{0.f};
  float phi{0.f};
  if (std::abs(bz) > kMinBz) {
    // Magnet on: tan(lambda) from inner→mid, phi from inner→outer with bending correction (TrackFitter inward).
    tanl = -std::abs(dzTan) / drTan;
    phi = std::atan2(dyPhi, dxPhi);
    if (std::abs(tanl) > kMinDz) {
      const float k = std::abs(o2::constants::math::B2C * bz);
      const float hz = (bz > 0.f) ? 1.f : -1.f;
      phi -= 0.5f * hz * invQPt * dzPhi * k / tanl;
    }
  } else {
    // Magnet off: both slopes from inner→outer (TrackFitter linear branch).
    tanl = -std::abs(dzPhi) / drPhi;
    phi = std::atan2(dyPhi, dxPhi);
  }

  ROOT::Math::SVector<double, 5> params{x0, y0, phi, tanl, invQPt};
  ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>> cov{};
  cov(0, 0) = sigma2XOuter > 0.f ? sigma2XOuter : kDefaultSigma2;
  cov(1, 1) = sigma2YOuter > 0.f ? sigma2YOuter : kDefaultSigma2;
  cov(2, 2) = 1.;
  cov(3, 3) = 1.;
  const double qptSigma = std::clamp(static_cast<double>(std::abs(invQPt)), 1., 10.);
  cov(4, 4) = qptSigma * qptSigma;

  track = {z0, params, cov, 0.};
  return true;
}

/// Convenience wrapper for CA cell building: clusters indexed inner → mid → outer in the time frame.
template <int NLayers>
inline bool mftBuildFwdTrackSeedFromCellClusters(const int hitLayers[3],
                                                 const int clusIds[3],
                                                 const TimeFrame<NLayers>& tf,
                                                 float bz,
                                                 float minPt,
                                                 o2::track::TrackParCovFwd& track)
{
  const auto& cInner = tf.getUnsortedClusters()[hitLayers[0]][clusIds[0]];
  const auto& cMid = tf.getUnsortedClusters()[hitLayers[1]][clusIds[1]];
  const auto& cOuter = tf.getUnsortedClusters()[hitLayers[2]][clusIds[2]];
  const auto& tfOuter = tf.getTrackingFrameInfoOnLayer(hitLayers[2])[clusIds[2]];
  return mftBuildFwdTrackSeedInward(cInner,
                                    cMid,
                                    cOuter,
                                    bz,
                                    minPt,
                                    track,
                                    tfOuter.covarianceTrackingFrame[0],
                                    tfOuter.covarianceTrackingFrame[2]);
}

/// Legacy helper: always returns a track object; invalid geometry yields a default-constructed state.
inline o2::track::TrackParCovFwd mftBuildFwdTrackSeed(const o2::its::Cluster& cInner,
                                                      const o2::its::Cluster& cMid,
                                                      const o2::its::Cluster& cOuter,
                                                      float bz,
                                                      float minPt)
{
  o2::track::TrackParCovFwd track;
  (void)mftBuildFwdTrackSeedInward(cInner, cMid, cOuter, bz, minPt, track);
  return track;
}

inline void mftFwdPropagateToZ(o2::track::TrackParCovFwd& track, float z, float bz)
{
  if (std::abs(bz) > 0.01f) {
    track.propagateToZhelix(z, bz);
  } else {
    track.propagateToZlinear(z);
  }
}

inline float mftFwdPredictedChi2Quiet(const o2::track::TrackParCovFwd& track, float x, float y, float sigma2X, float sigma2Y)
{
  const float dx = x - static_cast<float>(track.getX());
  const float dy = y - static_cast<float>(track.getY());
  const float vx = static_cast<float>(track.getSigma2X()) + sigma2X;
  const float vy = static_cast<float>(track.getSigma2Y()) + sigma2Y;
  if (vx <= 0.f || vy <= 0.f) {
    return o2::constants::math::VeryBig;
  }
  return dx * dx / vx + dy * dy / vy;
}

inline float mftFwdGetPredictedChi2Quiet(const o2::track::TrackParCovFwd& current, const o2::track::TrackParCovFwd& rhs)
{
  ROOT::Math::SVector<double, 5> diff{
    rhs.getX() - current.getX(),
    rhs.getY() - current.getY(),
    rhs.getPhi() - current.getPhi(),
    rhs.getTanl() - current.getTanl(),
    rhs.getInvQPt() - current.getInvQPt()};
  auto cov = current.getCovariances();
  cov += rhs.getCovariances();
  if (!cov.Invert()) {
    return o2::constants::math::VeryBig;
  }
  return static_cast<float>(ROOT::Math::Similarity(cov, diff));
}

/// ITS getPredictedChi2 analogue: propagate next cell track to current z, compare 5-param states.
inline float mftFwdCellGetPredictedChi2(const o2::track::TrackParCovFwd& current,
                                        o2::track::TrackParCovFwd& next,
                                        float bz)
{
  mftFwdPropagateToZ(next, static_cast<float>(current.getZ()), bz);
  return mftFwdGetPredictedChi2Quiet(current, next);
}

inline bool mftFwdAttachCluster(o2::track::TrackParCovFwd& track,
                                float z,
                                float x,
                                float y,
                                float sigma2X,
                                float sigma2Y,
                                float xOverX0,
                                float bz,
                                float maxChi2,
                                float& chi2,
                                bool checkChi2OnLast = false)
{
  mftFwdPropagateToZ(track, z, bz);
  if (xOverX0 > 0.f) {
    track.addMCSEffect(xOverX0);
  }
  const float predChi2 = mftFwdPredictedChi2Quiet(track, x, y, sigma2X, sigma2Y);
  if (checkChi2OnLast && predChi2 > maxChi2) {
    return false;
  }
  const std::array<float, 2> p{x, y};
  const std::array<float, 2> cov{sigma2X, sigma2Y};
  if (!track.update(p, cov)) {
    return false;
  }
  chi2 += predChi2;
  return true;
}

template <int NLayers>
inline bool mftFwdFitThreeClusters(o2::track::TrackParCovFwd& track,
                                   const std::array<int, 3>& hitLayers,
                                   const std::array<int, 3>& clusIds,
                                   const TimeFrame<NLayers>& tf,
                                   const TrackingParameters& params,
                                   float bz,
                                   float& chi2)
{
  chi2 = 0.f;
  for (int iC{2}; iC >= 0; --iC) {
    const int layer = hitLayers[iC];
    const int clIdx = clusIds[iC];
    const auto& tfInfo = tf.getTrackingFrameInfoOnLayer(layer)[clIdx];
    if (!mftFwdAttachCluster(track,
                             tfInfo.zCoordinate,
                             tfInfo.xCoordinate,
                             tfInfo.yCoordinate,
                             tfInfo.covarianceTrackingFrame[0],
                             tfInfo.covarianceTrackingFrame[2],
                             params.LayerxX0[layer],
                             bz,
                             params.MaxChi2ClusterAttachment,
                             chi2,
                             iC == 0)) {
      return false;
    }
  }
  return true;
}

/// Build a forward seed from three cell clusters and fit with propagate / MCS / χ² / update.
template <int NLayers>
inline bool mftFwdFitCellClusters(const int hitLayers[3],
                                  const int clusIds[3],
                                  const TimeFrame<NLayers>& tf,
                                  const TrackingParameters& params,
                                  float bz,
                                  o2::track::TrackParCovFwd& track,
                                  float& chi2)
{
  if (!mftBuildFwdTrackSeedFromCellClusters(hitLayers, clusIds, tf, bz, params.TrackletMinPt, track)) {
    return false;
  }
  const std::array<int, 3> layersArr{hitLayers[0], hitLayers[1], hitLayers[2]};
  const std::array<int, 3> clusArr{clusIds[0], clusIds[1], clusIds[2]};
  return mftFwdFitThreeClusters(track, layersArr, clusArr, tf, params, bz, chi2);
}

template <int NLayers>
inline bool mftFwdFitCellSeedWithChi2(const CellSeedTpl<o2::track::TrackParCovFwd>& cell,
                                      const TimeFrame<NLayers>& tf,
                                      const TrackingParameters& params,
                                      float bz,
                                      o2::track::TrackParCovFwd& fwdTrack,
                                      float& chi2)
{
  const int hitLayers[3]{
    mftLayerAtSlot(cell.getHitLayerMask(), 0),
    mftLayerAtSlot(cell.getHitLayerMask(), 1),
    mftLayerAtSlot(cell.getHitLayerMask(), 2)};
  const int clusId[3]{cell.getFirstClusterIndex(), cell.getSecondClusterIndex(), cell.getThirdClusterIndex()};
  return mftFwdFitCellClusters(hitLayers, clusId, tf, params, bz, fwdTrack, chi2);
}

template <int NLayers>
inline bool mftFwdFitCellSeed(const CellSeedTpl<o2::track::TrackParCovFwd>& cell,
                              const TimeFrame<NLayers>& tf,
                              const TrackingParameters& params,
                              float bz,
                              o2::track::TrackParCovFwd& fwdTrack)
{
  float chi2{0.f};
  return mftFwdFitCellSeedWithChi2(cell, tf, params, bz, fwdTrack, chi2);
}

inline o2::track::TrackParCovF mftFwdToBarrelTrack(const o2::track::TrackParCovFwd& fwd)
{
  o2::track::TrackParCovF barrel;
  fwd.toBarrelTrackParCov(barrel);
  return barrel;
}

template <int NLayers, typename SeedT>
inline std::vector<std::pair<int, int>> mftSortedHitsFromSeed(const SeedT& seed)
{
  std::vector<std::pair<int, int>> hits;
  const auto mask = seed.getHitLayerMask();
  for (int layer = mask.first(); layer <= mask.last(); ++layer) {
    if (!mask.has(layer)) {
      continue;
    }
    const int clIdx = seed.getCluster(layer);
    if (clIdx != o2::its::constants::UnusedIndex) {
      hits.emplace_back(layer, clIdx);
    }
  }
  std::sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
  return hits;
}

template <int NLayers, typename SeedT>
inline bool mftFwdFitSeedClusters(const SeedT& seed,
                                const TimeFrame<NLayers>& tf,
                                const TrackingParameters& params,
                                float bz,
                                o2::track::TrackParCovFwd& track,
                                float& chi2)
{
  const auto hits = mftSortedHitsFromSeed<NLayers>(seed);
  if (hits.size() < 3) {
    return false;
  }
  const auto& cOuter = tf.getUnsortedClusters()[hits[0].first][hits[0].second];
  const auto& cInner = tf.getUnsortedClusters()[hits.back().first][hits.back().second];
  const size_t midIdx = (hits.size() <= 3) ? 1u : (hits.size() / 2u);
  const auto& cMid = tf.getUnsortedClusters()[hits[midIdx].first][hits[midIdx].second];
  std::array<int, 3> hitLayers{hits.back().first, hits[midIdx].first, hits[0].first};
  std::array<int, 3> clusIds{hits.back().second, hits[midIdx].second, hits[0].second};
  if (!mftBuildFwdTrackSeedInward(cInner, cMid, cOuter, bz, params.TrackletMinPt, track)) {
    return false;
  }
  chi2 = 0.f;
  if (!mftFwdFitThreeClusters(track, hitLayers, clusIds, tf, params, bz, chi2)) {
    return false;
  }
  for (size_t i{3}; i < hits.size(); ++i) {
    const auto [layer, clIdx] = hits[i];
    const auto& tfInfo = tf.getTrackingFrameInfoOnLayer(layer)[clIdx];
    if (!mftFwdAttachCluster(track,
                             tfInfo.zCoordinate,
                             tfInfo.xCoordinate,
                             tfInfo.yCoordinate,
                             tfInfo.covarianceTrackingFrame[0],
                             tfInfo.covarianceTrackingFrame[2],
                             params.LayerxX0[layer],
                             bz,
                             params.MaxChi2ClusterAttachment,
                             chi2,
                             true)) {
      return false;
    }
  }
  return true;
}

/// Forward-fit neighbour attachment: propagate stored seed track, χ²-update at neighbour (ITS processNeighbours analogue).
template <int NLayers, typename SeedT>
inline bool mftFwdAttachNeighbourToSeed(const SeedT& seed,
                                        int neighbourLayer,
                                        int neighbourClIdx,
                                        const TimeFrame<NLayers>& tf,
                                        const TrackingParameters& params,
                                        float bz,
                                        o2::track::TrackParCovFwd& fwdTrack,
                                        float& chi2)
{
  fwdTrack = static_cast<const o2::track::TrackParCovFwd&>(seed);
  chi2 = seed.getChi2();
  const auto& trHit = tf.getTrackingFrameInfoOnLayer(neighbourLayer)[neighbourClIdx];
  return mftFwdAttachCluster(fwdTrack,
                             trHit.zCoordinate,
                             trHit.xCoordinate,
                             trHit.yCoordinate,
                             trHit.covarianceTrackingFrame[0],
                             trHit.covarianceTrackingFrame[2],
                             params.LayerxX0[neighbourLayer],
                             bz,
                             params.MaxChi2ClusterAttachment,
                             chi2,
                             true);
}

/// χ² compatibility between connected cells using stored forward track states (ITS cellsConnect analogue).
inline bool mftFwdCellsAreCompatible(const CellSeedTpl<o2::track::TrackParCovFwd>& currentCell,
                                     const CellSeedTpl<o2::track::TrackParCovFwd>& nextCell,
                                     float bz,
                                     float maxChi2)
{
  if (currentCell.getSecondClusterIndex() != nextCell.getFirstClusterIndex()) {
    return false;
  }
  o2::track::TrackParCovFwd nextFwd = static_cast<const o2::track::TrackParCovFwd&>(nextCell);
  const auto& currentFwd = static_cast<const o2::track::TrackParCovFwd&>(currentCell);
  const float chi2 = mftFwdCellGetPredictedChi2(currentFwd, nextFwd, bz);
  return chi2 <= maxChi2;
}

/// Legacy helper: full refit of both cells then compare (expensive; prefer mftFwdCellsAreCompatible above).
template <int NLayers>
inline bool mftFwdCellsAreCompatibleRefit(const CellSeedTpl<o2::track::TrackParCovFwd>& currentCell,
                                          const CellSeedTpl<o2::track::TrackParCovFwd>& nextCell,
                                          const TimeFrame<NLayers>& tf,
                                          const TrackingParameters& params,
                                          float bz)
{
  o2::track::TrackParCovFwd currentFwd;
  o2::track::TrackParCovFwd nextFwd;
  float chi2Current{0.f};
  float chi2Next{0.f};
  if (!mftFwdFitCellSeedWithChi2(currentCell, tf, params, bz, currentFwd, chi2Current) ||
      !mftFwdFitCellSeedWithChi2(nextCell, tf, params, bz, nextFwd, chi2Next)) {
    return false;
  }
  mftFwdPropagateToZ(nextFwd, static_cast<float>(currentFwd.getZ()), bz);
  const float chi2 = mftFwdGetPredictedChi2Quiet(currentFwd, nextFwd);
  return chi2 <= params.MaxChi2ClusterAttachment;
}

/// Forward-fit road rescue: extrapolate from seed hits, attach neighbour only.
template <int NLayers, typename SeedT>
inline bool mftFwdRoadAttachNeighbour(const SeedT& seed,
                                      int neighbourLayer,
                                      int neighbourClIdx,
                                      const TimeFrame<NLayers>& tf,
                                      const TrackingParameters& params,
                                      float bz,
                                      o2::track::TrackParCovFwd& fwdTrack,
                                      float& attachChi2)
{
  return mftFwdAttachNeighbourToSeed(seed, neighbourLayer, neighbourClIdx, tf, params, bz, fwdTrack, attachChi2);
}

} // namespace o2::itsmft::tracking::detail

#endif /* ALICEO2_ITSMFT_TRACKING_MFTFWDTRACKHELPERS_H_ */

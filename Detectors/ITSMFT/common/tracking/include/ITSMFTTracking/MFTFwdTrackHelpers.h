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
/// \brief Forward-track helpers for MFT CA cell fitting and final track refit
///

#ifndef ALICEO2_ITSMFT_TRACKING_MFTFWDTRACKHELPERS_H_
#define ALICEO2_ITSMFT_TRACKING_MFTFWDTRACKHELPERS_H_

#include <array>
#include <cmath>

#include "CommonConstants/MathConstants.h"
#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/MFTCATrack.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Constants.h"
#include "MFTTracking/Constants.h"
#include "ReconstructionDataFormats/TrackFwd.h"

namespace o2::itsmft::tracking::detail
{

/// MFT CA uses o2::mft::constants::mft::LayersNumber half-disk layers (same index as GeometryTGeo::getLayer).
/// Physical disk index is halfLayer / 2; ROFOverlapTable stores one LayerTiming per half-layer.

inline float mftLayerZ(int layer)
{
  return o2::mft::constants::mft::LayerZCoordinate()[layer];
}

inline float mftLayerMSAngle(int layer, const TrackingParameters& params)
{
  const float invP = 1.f / params.TrackletMinPt;
  const float zLayer = mftLayerZ(layer);
  const float rRef = params.LayerRadii[layer];
  const float tanlRef = (std::abs(rRef) > 1e-6f) ? zLayer / rRef : 0.f;
  const float absTanl = std::abs(tanlRef);
  const float cscLambda = (absTanl > 1e-6f) ? std::sqrt(1.f + tanlRef * tanlRef) / absTanl : 1e6f;
  return 0.0136f * invP * std::sqrt(params.LayerxX0[layer] * cscLambda);
}

inline void mftTrackletProject(float xCl, float yCl, float zCl, float pvX, float pvY, float pvZ,
                               int fromLayer, int toLayer, float bz, float minPt,
                               float& xProj, float& yProj)
{
  if (std::abs(bz) > 0.01f && minPt > 0.f) {
    const float zTo = mftLayerZ(toLayer);
    const float dxTan = xCl - pvX;
    const float dyTan = yCl - pvY;
    const float dzTan = zCl - pvZ;
    const float drTan = std::sqrt(dxTan * dxTan + dyTan * dyTan);
    float invQPt = 1.f / minPt;
    float tanl = (drTan > 1e-6f) ? -std::abs(dzTan) / drTan : -1.f;
    float phi = (drTan > 1e-6f) ? std::atan2(dyTan, dxTan) : 0.f;
    if (std::abs(tanl) > 1e-6f) {
      const float k = std::abs(o2::constants::math::B2C * bz);
      const float hz = (bz > 0.f) ? 1.f : -1.f;
      phi -= 0.5f * hz * invQPt * dzTan * k / tanl;
    }
    ROOT::Math::SVector<double, 5> params{xCl, yCl, phi, tanl, invQPt};
    ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>> cov{};
    cov(0, 0) = cov(1, 1) = cov(2, 2) = cov(3, 3) = 1.;
    const double qptSigma = std::clamp(static_cast<double>(std::abs(invQPt)), 1., 10.);
    cov(4, 4) = qptSigma * qptSigma;
    o2::track::TrackParCovFwd track{zCl, params, cov, 0.};
    track.propagateToZhelix(zTo, bz);
    xProj = static_cast<float>(track.getX());
    yProj = static_cast<float>(track.getY());
  } else {
    const float dz0 = mftLayerZ(fromLayer) - pvZ;
    if (std::abs(dz0) < 1e-6f) {
      xProj = xCl;
      yProj = yCl;
      return;
    }
    const float w = (mftLayerZ(toLayer) - pvZ) / dz0;
    xProj = pvX + w * (xCl - pvX);
    yProj = pvY + w * (yCl - pvY);
  }
}

inline void mftTrackletSigmaXY(float x0, float y0, float pvX, float pvY, float pvZ,
                               float sigma2X0, float sigma2Y0, float sigma2PvX, float sigma2PvY, float sigma2PvZ,
                               int fromLayer, int toLayer, float rLayerFrom, float meanDeltaZ, float msAngle,
                               float bendingAngle, float xProj, float yProj, float& sigmaX, float& sigmaY)
{
  const float zFrom = mftLayerZ(fromLayer);
  const float zTo = mftLayerZ(toLayer);
  const float dz0 = zFrom - pvZ;
  const float tanlRef = (std::abs(rLayerFrom) > 1e-6f) ? zFrom / rLayerFrom : 0.f;
  const float sigma2MS = meanDeltaZ * meanDeltaZ * msAngle * msAngle * (tanlRef * tanlRef + 1.f);
  if (std::abs(dz0) < o2::its::constants::Tolerance) {
    sigmaX = std::sqrt(sigma2X0 + sigma2PvX + sigma2MS);
    sigmaY = std::sqrt(sigma2Y0 + sigma2PvY + sigma2MS);
  } else {
    const float w = (zTo - pvZ) / dz0;
    const float invDz0 = w / dz0;
    const float sigma2W = invDz0 * invDz0 * sigma2PvZ;
    const float dx0 = x0 - pvX;
    const float dy0 = y0 - pvY;
    const float oneMinusW = 1.f - w;
    sigmaX = std::sqrt(oneMinusW * oneMinusW * sigma2PvX + w * w * sigma2X0 + dx0 * dx0 * sigma2W + sigma2MS);
    sigmaY = std::sqrt(oneMinusW * oneMinusW * sigma2PvY + w * w * sigma2Y0 + dy0 * dy0 * sigma2W + sigma2MS);
  }
  const float rProj = std::hypot(xProj, yProj);
  if (rProj > 1e-6f && bendingAngle > 0.f) {
    const float dr = rProj * bendingAngle;
    const float invR = 1.f / rProj;
    const float sinPhi = yProj * invR;
    const float cosPhi = xProj * invR;
    sigmaX = std::sqrt(sigmaX * sigmaX + dr * dr * sinPhi * sinPhi);
    sigmaY = std::sqrt(sigmaY * sigmaY + dr * dr * cosPhi * cosPhi);
  }
}

inline void mftFwdPropagateToZ(o2::track::TrackParCovFwd& track, float z, float bz)
{
  if (std::abs(bz) > 0.01f) {
    track.propagateToZhelix(z, bz);
  } else {
    track.propagateToZlinear(z);
  }
}

inline float mftFwdPredictedChi2(const o2::track::TrackParCovFwd& track, float x, float y, float sigma2X, float sigma2Y)
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

inline float mftFwdStateChi2(const o2::track::TrackParCovFwd& current, const o2::track::TrackParCovFwd& rhs)
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

inline bool mftFwdAttachCluster(o2::track::TrackParCovFwd& track, float z, float x, float y,
                                float sigma2X, float sigma2Y, float xOverX0, float bz, float maxChi2,
                                float& chi2, bool checkChi2OnLast = false)
{
  mftFwdPropagateToZ(track, z, bz);
  if (xOverX0 > 0.f) {
    track.addMCSEffect(xOverX0);
  }
  const float predChi2 = mftFwdPredictedChi2(track, x, y, sigma2X, sigma2Y);
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

/// Build inward forward seed at the outer cluster and Kalman-fit the three cell clusters.
template <int NLayers>
inline bool mftFwdFitCellClusters(const int hitLayers[3], const int clusIds[3],
                                  const TimeFrame<NLayers>& tf, const TrackingParameters& params,
                                  float bz, o2::track::TrackParCovFwd& track, float& chi2)
{
  const auto& cInner = tf.getUnsortedClusters()[hitLayers[0]][clusIds[0]];
  const auto& cMid = tf.getUnsortedClusters()[hitLayers[1]][clusIds[1]];
  const auto& cOuter = tf.getUnsortedClusters()[hitLayers[2]][clusIds[2]];
  if (cInner.zCoordinate <= cOuter.zCoordinate + 1.e-6f) {
    return false;
  }

  const float dxTan = cMid.xCoordinate - cInner.xCoordinate;
  const float dyTan = cMid.yCoordinate - cInner.yCoordinate;
  const float dzTan = cMid.zCoordinate - cInner.zCoordinate;
  const float drTan = std::sqrt(dxTan * dxTan + dyTan * dyTan);
  const float dxPhi = cOuter.xCoordinate - cInner.xCoordinate;
  const float dyPhi = cOuter.yCoordinate - cInner.yCoordinate;
  const float dzPhi = cOuter.zCoordinate - cInner.zCoordinate;
  const float drPhi = std::sqrt(dxPhi * dxPhi + dyPhi * dyPhi);
  if (drTan < 1.e-6f || std::abs(dzTan) < 1.e-6f || drPhi < 1.e-6f || std::abs(dzPhi) < 1.e-6f) {
    return false;
  }

  const float minPt = params.TrackletMinPt;
  const float invQPt = (minPt > 0.f) ? 1.f / minPt : 0.f;
  float tanl{0.f};
  float phi{0.f};
  if (std::abs(bz) > 0.01f) {
    tanl = -std::abs(dzTan) / drTan;
    phi = std::atan2(dyPhi, dxPhi);
    if (std::abs(tanl) > 1.e-6f) {
      const float k = std::abs(o2::constants::math::B2C * bz);
      const float hz = (bz > 0.f) ? 1.f : -1.f;
      phi -= 0.5f * hz * invQPt * dzPhi * k / tanl;
    }
  } else {
    tanl = -std::abs(dzPhi) / drPhi;
    phi = std::atan2(dyPhi, dxPhi);
  }

  const auto& tfOuter = tf.getTrackingFrameInfoOnLayer(hitLayers[2])[clusIds[2]];
  ROOT::Math::SVector<double, 5> seedParams{cOuter.xCoordinate, cOuter.yCoordinate, phi, tanl, invQPt};
  ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>> seedCov{};
  seedCov(0, 0) = tfOuter.covarianceTrackingFrame[0] > 0.f ? tfOuter.covarianceTrackingFrame[0] : 1.f;
  seedCov(1, 1) = tfOuter.covarianceTrackingFrame[2] > 0.f ? tfOuter.covarianceTrackingFrame[2] : 1.f;
  seedCov(2, 2) = seedCov(3, 3) = 1.;
  const double qptSigma = std::clamp(static_cast<double>(std::abs(invQPt)), 1., 10.);
  seedCov(4, 4) = qptSigma * qptSigma;
  track = {cOuter.zCoordinate, seedParams, seedCov, 0.};

  chi2 = 0.f;
  for (int iC{2}; iC >= 0; --iC) {
    const int layer = hitLayers[iC];
    const int clIdx = clusIds[iC];
    const auto& tfInfo = tf.getTrackingFrameInfoOnLayer(layer)[clIdx];
    if (!mftFwdAttachCluster(track, tfInfo.zCoordinate, tfInfo.xCoordinate, tfInfo.yCoordinate,
                             tfInfo.covarianceTrackingFrame[0], tfInfo.covarianceTrackingFrame[2],
                             params.LayerxX0[layer], bz, params.MaxChi2ClusterAttachment, chi2, iC == 0)) {
      return false;
    }
  }
  return true;
}

template <int NLayers, typename SeedT>
inline bool mftFwdAttachNeighbourToSeed(const SeedT& seed, int neighbourLayer, int neighbourClIdx,
                                        const TimeFrame<NLayers>& tf, const TrackingParameters& params,
                                        float bz, o2::track::TrackParCovFwd& fwdTrack, float& chi2)
{
  fwdTrack = static_cast<const o2::track::TrackParCovFwd&>(seed);
  chi2 = seed.getChi2();
  const auto& trHit = tf.getTrackingFrameInfoOnLayer(neighbourLayer)[neighbourClIdx];
  return mftFwdAttachCluster(fwdTrack, trHit.zCoordinate, trHit.xCoordinate, trHit.yCoordinate,
                             trHit.covarianceTrackingFrame[0], trHit.covarianceTrackingFrame[2],
                             params.LayerxX0[neighbourLayer], bz, params.MaxChi2ClusterAttachment, chi2, true);
}

inline bool mftFwdCellsAreCompatible(const CellSeedTpl<o2::track::TrackParCovFwd>& currentCell,
                                   const CellSeedTpl<o2::track::TrackParCovFwd>& nextCell,
                                   float bz, float maxChi2)
{
  if (currentCell.getSecondClusterIndex() != nextCell.getFirstClusterIndex()) {
    return false;
  }
  o2::track::TrackParCovFwd nextFwd = static_cast<const o2::track::TrackParCovFwd&>(nextCell);
  const auto& currentFwd = static_cast<const o2::track::TrackParCovFwd&>(currentCell);
  mftFwdPropagateToZ(nextFwd, static_cast<float>(currentFwd.getZ()), bz);
  return mftFwdStateChi2(currentFwd, nextFwd) <= maxChi2;
}

} // namespace o2::itsmft::tracking::detail

namespace o2::itsmft::tracking
{

bool refitTrackFwd(const TrackSeedN<o2::mft::constants::mft::LayersNumber>& seed,
                   MFTCATrack& track,
                   const TimeFrame<o2::mft::constants::mft::LayersNumber>& tf,
                   const TrackingParameters& params,
                   float bz);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_MFTFWDTRACKHELPERS_H_ */

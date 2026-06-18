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
/// \file DetectorTraits.h
/// \brief Minimal detector-specific hooks for shared CA tracking (ITS vs MFT)
///

#ifndef ALICEO2_ITSMFT_TRACKING_DETECTOR_TRAITS_H_
#define ALICEO2_ITSMFT_TRACKING_DETECTOR_TRAITS_H_

#include "DetectorsBase/Propagator.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "ITSMFTTracking/CATrackTypes.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/Constants.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/Cell.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Tracklet.h"
#include "ReconstructionDataFormats/Track.h"
#include "MFTTracking/Constants.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include <cmath>
#include <limits>

namespace o2::itsmft::tracking
{

namespace detail
{
inline float mftDistanceToSeedSquared(const o2::its::Cluster& c1, const o2::its::Cluster& c2, const o2::its::Cluster& c)
{
  const float dxSeed = c2.xCoordinate - c1.xCoordinate;
  const float dySeed = c2.yCoordinate - c1.yCoordinate;
  const float dzSeed = c2.zCoordinate - c1.zCoordinate;
  if (std::abs(dzSeed) < 1e-9f) {
    return std::numeric_limits<float>::max();
  }
  const float invdzSeed = (c.zCoordinate - c1.zCoordinate) / dzSeed;
  const float xSeed = c1.xCoordinate + dxSeed * invdzSeed;
  const float ySeed = c1.yCoordinate + dySeed * invdzSeed;
  const float dx = c.xCoordinate - xSeed;
  const float dy = c.yCoordinate - ySeed;
  return dx * dx + dy * dy;
}

inline float mftLayerZ(int layer)
{
  return o2::mft::constants::mft::LayerZCoordinate()[layer];
}

/// Expected cluster z on toLayer for a straight line from the vertex through cluster on fromLayer.
inline float mftExpectedZAtLayer(float clusterZ, int fromLayer, int toLayer, float pvZ)
{
  const float zFrom = mftLayerZ(fromLayer);
  const float dzLayers = mftLayerZ(toLayer) - zFrom;
  const float denom = zFrom - pvZ;
  if (std::abs(denom) < 1e-6f) {
    return clusterZ + dzLayers;
  }
  return clusterZ + dzLayers * ((clusterZ - pvZ) / denom);
}

/// Per-layer forward MCS angle (TrackParCovFwd::addMCSEffect-style), for quadrature sum over layers.
inline float mftLayerMSAngle(int layer, const TrackingParameters& params)
{
  const float invP = 1.f / params.TrackletMinPt;
  const float zLayer = mftLayerZ(layer);
  const float rRef = params.LayerRadii[layer];
  const float tanlRef = (std::abs(rRef) > 1e-6f) ? zLayer / rRef : 0.f;
  const float absTanl = std::abs(tanlRef);
  const float cscLambda = (absTanl > 1e-6f) ? std::sqrt(1.f + tanlRef * tanlRef) / absTanl : 1e6f;
  const float pathLengthOverX0 = params.LayerxX0[layer] * cscLambda;
  return 0.0136f * invP * std::sqrt(pathLengthOverX0);
}

/// Straight-line extrapolation of cluster to toLayer through the primary vertex.
inline void mftLineProject(float xCl, float yCl, float zCl, float pvX, float pvY, float pvZ,
                           int fromLayer, int toLayer, float& xProj, float& yProj)
{
  const float zFrom = mftLayerZ(fromLayer);
  const float zTo = mftLayerZ(toLayer);
  const float dz0 = zFrom - pvZ;
  if (std::abs(dz0) < 1e-6f) {
    xProj = xCl;
    yProj = yCl;
    return;
  }
  const float w = (zTo - pvZ) / dz0;
  xProj = pvX + w * (xCl - pvX);
  yProj = pvY + w * (yCl - pvY);
}

/// Helix extrapolation at TrackletMinPt from cluster on fromLayer to toLayer (forward parametrization).
inline void mftHelixProject(float xCl, float yCl, float zCl, float pvX, float pvY, float pvZ,
                            int toLayer, float bz, float minPt, float& xProj, float& yProj)
{
  const float zTo = mftLayerZ(toLayer);
  const float dxTan = xCl - pvX;
  const float dyTan = yCl - pvY;
  const float dzTan = zCl - pvZ;
  const float drTan = std::sqrt(dxTan * dxTan + dyTan * dyTan);
  float invQPt = (minPt > 0.f) ? 1.f / minPt : 0.f;
  float tanl = (drTan > 1e-6f) ? -std::abs(dzTan) / drTan : -1.f;
  float phi = (drTan > 1e-6f) ? std::atan2(dyTan, dxTan) : 0.f;
  if (std::abs(bz) > 0.01f && std::abs(tanl) > 1e-6f) {
    const float k = std::abs(o2::constants::math::B2C * bz);
    const float hz = (bz > 0.f) ? 1.f : -1.f;
    phi -= 0.5f * hz * invQPt * dzTan * k / tanl;
  }
  ROOT::Math::SVector<double, 5> params{xCl, yCl, phi, tanl, invQPt};
  ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>> cov{};
  cov(0, 0) = 1.;
  cov(1, 1) = 1.;
  cov(2, 2) = 1.;
  cov(3, 3) = 1.;
  const double qptSigma = std::clamp(static_cast<double>(std::abs(invQPt)), 1., 10.);
  cov(4, 4) = qptSigma * qptSigma;
  o2::track::TrackParCovFwd track{zCl, params, cov, 0.};
  if (std::abs(bz) > 0.01f) {
    track.propagateToZhelix(zTo, bz);
  } else {
    track.propagateToZlinear(zTo);
  }
  xProj = static_cast<float>(track.getX());
  yProj = static_cast<float>(track.getY());
}

/// Project cluster to toLayer: helix at min pT when B is on, otherwise vertex line.
inline void mftTrackletProject(float xCl, float yCl, float zCl, float pvX, float pvY, float pvZ,
                               int fromLayer, int toLayer, float bz, float minPt,
                               float& xProj, float& yProj)
{
  if (std::abs(bz) > 0.01f && minPt > 0.f) {
    mftHelixProject(xCl, yCl, zCl, pvX, pvY, pvZ, toLayer, bz, minPt, xProj, yProj);
  } else {
    mftLineProject(xCl, yCl, zCl, pvX, pvY, pvZ, fromLayer, toLayer, xProj, yProj);
  }
}

/// Tracklet x/y resolution at toLayer: vertex line errors + MS + helix bending at min pT.
inline void mftTrackletSigmaXY(float x0,
                               float y0,
                               float pvX,
                               float pvY,
                               float pvZ,
                               float sigma2X0,
                               float sigma2Y0,
                               float sigma2PvX,
                               float sigma2PvY,
                               float sigma2PvZ,
                               int fromLayer,
                               int toLayer,
                               float rLayerFrom,
                               float meanDeltaZ,
                               float msAngle,
                               float bendingAngle,
                               float xProj,
                               float yProj,
                               float& sigmaX,
                               float& sigmaY)
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
    const float sigma2BendX = dr * dr * sinPhi * sinPhi;
    const float sigma2BendY = dr * dr * cosPhi * cosPhi;
    sigmaX = std::sqrt(sigmaX * sigmaX + sigma2BendX);
    sigmaY = std::sqrt(sigmaY * sigmaY + sigma2BendY);
  }
}

/// Diamond vertex z extent used for conical R spread (legacy MFT CA indexing).
inline void mftDiamondZExtents(float pvZ, float pvSigmaZ, float nSigmaCut, float& zVtxMin, float& zVtxMax)
{
  const float zSpread = nSigmaCut * pvSigmaZ;
  zVtxMin = pvZ - zSpread;
  zVtxMax = pvZ + zSpread;
}

/// Radial spread at toLayer from seed radius at fromLayer and diamond z extremes.
inline void mftRSpreadAtLayer(float seedRadius,
                              float zLayerFrom,
                              float zLayerTo,
                              float zVtxMin,
                              float zVtxMax,
                              float& rMin,
                              float& rMax)
{
  const float absZFrom = std::abs(zLayerFrom);
  const float absZTo = std::abs(zLayerTo);
  const float denomMin = zVtxMax + absZFrom;
  const float denomMax = absZFrom + zVtxMin;
  rMin = (std::abs(denomMin) > 1.e-6f) ? seedRadius * (zVtxMax + absZTo) / denomMin : seedRadius;
  rMax = (std::abs(denomMax) > 1.e-6f) ? seedRadius * (absZTo + zVtxMin) / denomMax : seedRadius;
  if (rMin > rMax) {
    const float tmp = rMin;
    rMin = rMax;
    rMax = tmp;
  }
}

/// Normalised transverse chi2 in cone-projected x-y using per-axis sigmas.
inline float mftTrackletTransverseChi2(float dx, float dy, float sigmaX, float sigmaY)
{
  const float invSigmaX2 = (sigmaX > 0.f) ? 1.f / (sigmaX * sigmaX) : 0.f;
  const float invSigmaY2 = (sigmaY > 0.f) ? 1.f / (sigmaY * sigmaY) : 0.f;
  return dx * dx * invSigmaX2 + dy * dy * invSigmaY2;
}

/// Legacy MFT conical road scale: (1 + dz/z_ref)^2 for layerFrom -> layerTo.
inline float mftConicalRoadR2Scale(int layerFrom, int layerTo)
{
  const float zFrom = mftLayerZ(layerFrom);
  if (std::abs(zFrom) < 1e-6f) {
    return 1.f;
  }
  const float dCone = 1.f + (mftLayerZ(layerTo) - zFrom) / zFrom;
  return dCone * dCone;
}

/// Geometric road acceptance with conical scaling of CellRoadRCut (legacy ROADclsRCut behaviour).
inline bool mftGeometricRoadAccept(const o2::its::Cluster& cA,
                                   const o2::its::Cluster& cB,
                                   const o2::its::Cluster& cN,
                                   int refLayerInner,
                                   int neighbourLayer,
                                   float roadRCut)
{
  const float r2Cut = roadRCut * roadRCut * mftConicalRoadR2Scale(refLayerInner, neighbourLayer);
  return mftDistanceToSeedSquared(cA, cB, cN) < r2Cut;
}
} // namespace detail

/// Per-detector differences in refit, track acceptance, and index-table setup.
/// Everything else stays in TrackerTraits and matches ITS line-for-line.
template <int NLayers>
struct DetectorTraits {
  using TrackType = CATrackType<NLayers>;
  using TrackSeedN = o2::itsmft::tracking::TrackSeedN<NLayers>;
  using CellSeedN = o2::itsmft::tracking::CellSeedN<NLayers>;
  using TimeFrameN = TimeFrame<NLayers>;
  static constexpr o2::detectors::DetID::ID DetId = constants::detIdFromNLayers<NLayers>();

  static bool refitSeed(const TrackSeedN& seed,
                        TrackType& track,
                        const TrackingParameters& params,
                        float bz,
                        TimeFrameN& tf,
                        const o2::its::TrackingFrameInfo* const tfInfos[NLayers],
                        const o2::its::Cluster* const unsortedClusters[NLayers],
                        const o2::base::PropagatorImpl<float>* propagator);

  static void sortRefittedTracks(bounded_vector<TrackType>& tracks);
  static void finalizeAcceptedTrack(TrackType& track);
  static bool sameTrackSign(const TrackType& t1, const TrackType& t2);

  static void configureIndexTableUtils(IndexTableUtils<NLayers>& utils, const TrackingParameters& params);

  /// ITS: helix fit on three clusters; MFT: forward fit (TrackParCovFwd) on three clusters.
  static bool fitCellClusters(const int clusId[3],
                              const int hitLayers[3],
                              const TimeFrameN& tf,
                              const TrackingParameters& params,
                              float bz,
                              const o2::base::PropagatorImpl<float>* propagator,
                              o2::track::TrackParCovF& track,
                              float& chi2);

  /// ITS: propagate barrel seed and χ²-update; MFT: forward-fit seed and attach neighbour cluster.
  static bool attachNeighbourToSeed(TrackSeedN& seed,
                                    int neighbourLayer,
                                    int neighbourCluster,
                                    const TimeFrameN& tf,
                                    const TrackingParameters& params,
                                    float bz,
                                    const o2::base::PropagatorImpl<float>* propagator);

  /// ITS: barrel χ² between connected cells; MFT: forward χ² between connected cells.
  static bool cellsAreCompatible(const CellSeedN& currentCell,
                                 const CellSeedN& nextCell,
                                 const TimeFrameN& tf,
                                 const TrackingParameters& params,
                                 float bz);

  static bool validateMFTCellClusters(const o2::its::Cluster& c0, int layer0,
                                      const o2::its::Cluster& c1, int layer1,
                                      const o2::its::Cluster& c2, int layer2,
                                      float r2Cut);
  static bool mftCellsConnect(const o2::its::Cluster& cEnd, const o2::its::Cluster& cStart, float r2Cut);
};

template <o2::detectors::DetID::ID DetId, int NLayers>
struct TrackingLoadPolicy {
  static void configureBeamPosition(TimeFrame<NLayers>& tf,
                                    const TrackingParameters& p,
                                    const o2::dataformats::MeanVertexObject* meanVertex,
                                    bool overrideBeamEstimation);
};

template <int NLayers>
using TrackingLoadPolicyN = TrackingLoadPolicy<constants::detIdFromNLayers<NLayers>(), NLayers>;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_DETECTOR_TRAITS_H_ */

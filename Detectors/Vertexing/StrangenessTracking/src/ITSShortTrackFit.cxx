// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "StrangenessTracking/ITSShortTrackFit.h"

#include "ITSBase/GeometryTGeo.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/TrackHelpers.h"
#include <cmath>

namespace o2::strangeness_tracking
{
namespace
{
template <int N, typename T>
bool positiveDefinite(const std::array<T, N*(N + 1) / 2>& cov)
{
  // Check the packed lower triangle without changing/clamping the covariance.
  double l[N][N]{};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j <= i; ++j) {
      double v = cov[i * (i + 1) / 2 + j];
      for (int k = 0; k < j; ++k) {
        v -= l[i][k] * l[j][k];
      }
      if (i == j) {
        if (!(v > 0.)) {
          return false;
        }
        l[i][j] = std::sqrt(v);
      } else {
        l[i][j] = v / l[j][j];
      }
    }
  }
  return true;
}

bool validTrack(const o2::track::TrackParCovF& track)
{
  return std::abs(track.getSnp()) < 1.f && positiveDefinite<5>(track.getCov());
}
struct PVPointingMeasurement {
  std::array<float, 2> position{};
  std::array<float, 3> covariance{};
};

// PV uncertainty enters the compatibility residual only. Fixed preselection
// slopes project its full 3D covariance; it is never assimilated into the track.
bool projectPVForPointing(const std::array<float, 3>& position,
                          const std::array<float, 6>& covariance,
                          float alpha, float referenceX, float snp, float tgl,
                          PVPointingMeasurement& result)
{
  if (!(std::abs(snp) < 1.f) || !positiveDefinite<3>(covariance)) {
    return false;
  }
  const double sn = std::sin(alpha), cs = std::cos(alpha);
  const double csp = std::sqrt((1. - snp) * (1. + snp));
  const double a = snp / csp, b = tgl / csp;
  const double x = cs * position[0] + sn * position[1];
  const double y = -sn * position[0] + cs * position[1];
  PVPointingMeasurement projected;
  projected.position = {float(y + a * (referenceX - x)), float(position[2] + b * (referenceX - x))};
  // Rotation to the track frame followed by J=[[-a,1,0],[-b,0,1]].
  const double jacobian[2][3]{{-sn - a * cs, cs - a * sn, 0.}, {-b * cs, -b * sn, 1.}};
  constexpr int index[3][3]{{0, 1, 3}, {1, 2, 4}, {3, 4, 5}};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j <= i; ++j) {
      double value = 0.;
      for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
          value += jacobian[i][k] * covariance[index[k][l]] * jacobian[j][l];
        }
      }
      projected.covariance[i * (i + 1) / 2 + j] = value;
    }
  }
  if (!positiveDefinite<2>(projected.covariance)) {
    return false;
  }
  result = projected;
  return true;
}

} // namespace

bool fitITSShortTrack(
  const std::array<o2::BaseCluster<float>, 3>& trackingClusters,
  const o2::dataformats::PrimaryVertex& pv, float bz,
  o2::base::Propagator::MatCorrType corr,
  const std::array<float, 3>& layerxX0, float minPt, float maxPointingChi2,
  float maxChi2NDF, o2::track::TrackParCovF& out, float& hitChi2,
  float& pointingChi2, int& ndf, std::array<float, 2>* pointingDCA,
  o2::track::TrackParCovF* outerOut)
{
  if (std::abs(bz) < 0.01f || minPt < 0.f || !(maxPointingChi2 > 0.f) ||
      !(maxChi2NDF > 0.f) ||
      !positiveDefinite<3>(std::array<float, 6>{
        pv.getSigmaX2(), pv.getSigmaXY(), pv.getSigmaY2(), pv.getSigmaXZ(),
        pv.getSigmaYZ(), pv.getSigmaZ2()})) {
    return false;
  }

  auto* geometry = o2::its::GeometryTGeo::Instance();
  std::array<o2::its::Cluster, 3> global;
  std::array<o2::its::TrackingFrameInfo, 3> frame;
  for (int layer = 0; layer < 3; ++layer) {
    const auto& cluster = trackingClusters[layer];
    if (cluster.getSensorID() < 0 ||
        cluster.getSensorID() >= geometry->getNumberOfChips() ||
        geometry->getLayer(cluster.getSensorID()) != layer ||
        !positiveDefinite<2>(std::array<float, 3>{cluster.getSigmaY2(),
                                                  cluster.getSigmaYZ(),
                                                  cluster.getSigmaZ2()}) ||
        layerxX0[layer] < 0.f) {
      return false;
    }
    const auto xyz = cluster.getXYZGloRot(*geometry);
    global[layer] = o2::its::Cluster{xyz.X(), xyz.Y(), xyz.Z(), layer};
    frame[layer] = o2::its::TrackingFrameInfo{
      xyz.X(),
      xyz.Y(),
      xyz.Z(),
      cluster.getX(),
      geometry->getSensorRefAlpha(cluster.getSensorID()),
      {cluster.getY(), cluster.getZ()},
      {cluster.getSigmaY2(), cluster.getSigmaYZ(), cluster.getSigmaZ2()}};
  }

  // Match the ITS seed's outward momentum convention, starting at L2 and
  // transporting inward. The seed covariance contains hit information and is
  // deliberately discarded before the actual measurement fit.
  auto fitted =
    o2::its::track::buildTrackSeed(global[0], global[1], frame[2], bz);
  const double dx = global[2].xCoordinate - global[0].xCoordinate;
  const double dy = global[2].yCoordinate - global[0].yCoordinate;
  const double distance = std::hypot(dx, dy);
  const double area = (global[1].xCoordinate - global[0].xCoordinate) * dy -
                      (global[1].yCoordinate - global[0].yCoordinate) * dx;
  if (!(distance > 1.e-5)) {
    return false;
  }
  if (std::abs(area) < o2::its::constants::Tolerance) {
    fitted.setQ2Pt(0.f);
  }
  // Obtain the outer tangent from the chord and curvature, avoiding the ITS
  // circle-centre infinity sentinel and cancellation for almost straight hits.
  const float alpha = frame[2].alphaTrackingFrame;
  const double localDx = dx * std::cos(alpha) + dy * std::sin(alpha);
  const double localDy = -dx * std::sin(alpha) + dy * std::cos(alpha);
  const double sinBend =
    0.5 * fitted.getQ2Pt() * bz * o2::constants::math::B2C * distance;
  if (!(std::abs(sinBend) < 1.)) {
    return false;
  }
  const double cosBend = std::sqrt((1. - sinBend) * (1. + sinBend));
  if (!(localDx * cosBend - localDy * sinBend > 0.)) {
    return false;
  }
  fitted.setSnp((localDy * cosBend + localDx * sinBend) / distance);
  fitted.resetCovariance();
  if (!validTrack(fitted)) {
    return false;
  }

  const auto* propagator = o2::base::Propagator::Instance();
  float hitSum = 0.f;
  for (int layer = 2; layer >= 0; --layer) {
    const auto& measurement = frame[layer];
    if (!fitted.rotate(measurement.alphaTrackingFrame) ||
        !propagator->propagateToX(fitted, measurement.xTrackingFrame, bz,
                                  o2::base::Propagator::MAX_SIN_PHI,
                                  o2::base::Propagator::MAX_STEP, corr) ||
        (corr == o2::base::Propagator::MatCorrType::USEMatCorrNONE &&
         !fitted.correctForMaterial(layerxX0[layer],
                                    layerxX0[layer] * o2::its::constants::Radl *
                                      o2::its::constants::Rho,
                                    true))) {
      return false;
    }
    const float chi2 = fitted.getPredictedChi2Quiet(
      measurement.positionTrackingFrame, measurement.covarianceTrackingFrame);
    if (chi2 < 0.f ||
        !fitted.update(measurement.positionTrackingFrame,
                       measurement.covarianceTrackingFrame) ||
        !validTrack(fitted)) {
      return false;
    }
    hitSum += chi2;
  }

  constexpr int fitNDF = 1; // six ITS coordinates, five fitted parameters
  if (!(hitSum < fitNDF * maxChi2NDF) ||
      (minPt > 0.f && !(std::abs(fitted.getQ2Pt()) < 1.f / minPt))) {
    return false;
  }

  // Persist the L0 hit fit, entirely independent of the PV. Pointing changes
  // only the reference surface of a temporary copy and never updates it.
  auto pointing = fitted;
  if (!propagator->propagateToDCA(pv, pointing, bz,
                                  o2::base::Propagator::MAX_STEP, corr) ||
      !validTrack(pointing)) {
    return false;
  }
  PVPointingMeasurement pvMeasurement;
  if (!projectPVForPointing({pv.getX(), pv.getY(), pv.getZ()},
                            {pv.getSigmaX2(), pv.getSigmaXY(), pv.getSigmaY2(), pv.getSigmaXZ(), pv.getSigmaYZ(), pv.getSigmaZ2()},
                            pointing.getAlpha(), pointing.getX(), pointing.getSnp(), pointing.getTgl(), pvMeasurement)) {
    return false;
  }
  const auto& pvYZ = pvMeasurement.position;
  const auto& pvCov = pvMeasurement.covariance;
  const std::array<float, 2> dca{pointing.getY() - pvYZ[0], pointing.getZ() - pvYZ[1]};
  const float pvScore = pointing.getPredictedChi2Quiet(pvYZ, pvCov);
  // Pointing is a separate selection and contributes no fit chi2/covariance.
  if (pvScore < 0.f || !(pvScore < maxPointingChi2)) {
    return false;
  }
  auto outer = fitted;
  if (outerOut) {
    for (int layer = 1; layer <= 2; ++layer) {
      const auto& surface = frame[layer];
      if (!outer.rotate(surface.alphaTrackingFrame) ||
          !propagator->propagateToX(outer, surface.xTrackingFrame, bz,
                                    o2::base::Propagator::MAX_SIN_PHI,
                                    o2::base::Propagator::MAX_STEP, corr) ||
          (corr == o2::base::Propagator::MatCorrType::USEMatCorrNONE &&
           !outer.correctForMaterial(layerxX0[layer],
                                     -layerxX0[layer] * o2::its::constants::Radl * o2::its::constants::Rho, true)) ||
          !validTrack(outer)) {
        return false;
      }
    }
  }
  out = fitted;
  if (outerOut) {
    *outerOut = outer;
  }
  hitChi2 = hitSum;
  pointingChi2 = pvScore;
  ndf = fitNDF;
  if (pointingDCA) {
    *pointingDCA = dca;
  }
  return true;
}
} // namespace o2::strangeness_tracking

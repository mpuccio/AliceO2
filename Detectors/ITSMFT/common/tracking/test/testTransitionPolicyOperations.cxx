// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT TransitionPolicyOperations
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <utility>

#include <boost/test/unit_test.hpp>

#include <TGeoGlobalMagField.h>

#include "DetectorsCommonDataFormats/DetID.h"
#include "Field/MagneticField.h"
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/TransitionPolicyBinding.h"
#include "ITSMFTTracking/TransitionPolicyOperations.h"
#include "ITStracking/TrackHelpers.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

struct PropagatorFieldFixture {
  PropagatorFieldFixture()
  {
    if (!TGeoGlobalMagField::Instance()->GetField()) {
      TGeoGlobalMagField::Instance()->SetField(o2::field::MagneticField::createNominalField(5, true));
      TGeoGlobalMagField::Instance()->Lock();
    }
  }
};

BOOST_GLOBAL_FIXTURE(PropagatorFieldFixture);

/// Focused numerical-parity coverage for the first D007 policy-boundary
/// operation migrated off the legacy per-detector branch (Architecture.md
/// §10, cellsAreCompatible). These tests do not exercise TrackerTraits'
/// production traversal -- see the handoff note on scope.

namespace
{

o2::track::TrackParCovF makeBarrelState(float x, float alpha, float y, float z, float snp, float tgl, float q2pt)
{
  const o2::track::TrackParCovF::params_t par{y, z, snp, tgl, q2pt};
  o2::track::TrackParCovF::covMat_t cov{};
  cov[0] = 1.e-4f;  // sigma_Y^2
  cov[2] = 1.e-4f;  // sigma_Z^2
  cov[5] = 1.e-4f;  // sigma_Snp^2
  cov[9] = 1.e-4f;  // sigma_Tgl^2
  cov[14] = 1.e-2f; // sigma_Q2Pt^2
  return o2::track::TrackParCovF{x, alpha, par, cov};
}

CellSeedTpl<o2::track::TrackParCovF> makeBarrelCell(const o2::track::TrackParCovF& state, int cl0, int cl1, int cl2)
{
  return CellSeedTpl<o2::track::TrackParCovF>{0, cl0, cl1, cl2, 0, 1, state, 0.f, o2::its::TimeEstBC{}};
}

o2::track::TrackParCovFwd makeForwardState(float z, float x, float y, float phi, float tanl, float invQPt)
{
  const o2::track::SMatrix5 par{x, y, phi, tanl, invQPt};
  o2::track::SMatrix55Sym cov{};
  cov(0, 0) = 1.e-2;
  cov(1, 1) = 1.e-2;
  cov(2, 2) = 1.e-3;
  cov(3, 3) = 1.e-3;
  cov(4, 4) = 1.e-2;
  return o2::track::TrackParCovFwd{z, par, cov, 0.};
}

CellSeedTpl<o2::track::TrackParCovFwd> makeForwardCell(const o2::track::TrackParCovFwd& state, int cl0, int cl1, int cl2)
{
  return CellSeedTpl<o2::track::TrackParCovFwd>{0, cl0, cl1, cl2, 0, 1, state, 0.f, o2::its::TimeEstBC{}};
}

constexpr float Bz = 0.5f;
constexpr auto NoMaterialCorrection = o2::base::PropagatorF::MatCorrType::USEMatCorrNONE;

o2::its::TrackingFrameInfo makeBarrelHit(float xTF, float alpha, float y, float z, float sigma2Y = 1.e-4f, float sigma2Z = 1.e-4f)
{
  return o2::its::TrackingFrameInfo{xTF, y, z, xTF, alpha, {y, z}, {sigma2Y, 0.f, sigma2Z}};
}

o2::its::TrackingFrameInfo makeDiskHit(float z, float x, float y, float sigma2X = 1.e-2f, float sigma2Y = 1.e-2f)
{
  return o2::its::TrackingFrameInfo{x, y, z, 0.f, 0.f, {x, y}, {sigma2X, 0.f, sigma2Y}};
}

o2::its::Vertex makeVertex(float x, float y, float z,
                           float sigma2X, float sigma2Y, float sigma2Z,
                           unsigned short contributors = 1)
{
  const float position[3]{x, y, z};
  const float covariance[6]{sigma2X, 0.f, sigma2Y, 0.f, 0.f, sigma2Z};
  return o2::its::Vertex{position, covariance, contributors, 1.f};
}

bool legacyBarrelAttach(o2::track::TrackParCovF& state, const o2::its::TrackingFrameInfo& hit,
                        float xOverX0, o2::base::PropagatorF::MatCorrType corrType,
                        float bz, float maxChi2, float& chi2)
{
  if (!state.rotate(hit.alphaTrackingFrame)) {
    return false;
  }
  if (!o2::base::Propagator::Instance()->propagateToX(state, hit.xTrackingFrame, bz,
                                                      o2::base::PropagatorImpl<float>::MAX_SIN_PHI,
                                                      o2::base::PropagatorImpl<float>::MAX_STEP, corrType)) {
    return false;
  }
  if (corrType == NoMaterialCorrection &&
      !state.correctForMaterial(xOverX0, xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho, true)) {
    return false;
  }
  const float predictedChi2 = state.getPredictedChi2Quiet(hit.positionTrackingFrame, hit.covarianceTrackingFrame);
  if (predictedChi2 > maxChi2 || predictedChi2 < 0.f) {
    return false;
  }
  chi2 += predictedChi2;
  return state.o2::track::TrackParCov::update(hit.positionTrackingFrame, hit.covarianceTrackingFrame);
}

bool legacyDiskAttach(o2::track::TrackParCovFwd& state, const o2::its::TrackingFrameInfo& hit,
                      float xOverX0, float bz, float maxChi2, float& chi2)
{
  auto updated = state;
  float updatedChi2 = chi2;
  if (!detail::mftFwdAttachCluster(updated, hit.zCoordinate, hit.xCoordinate, hit.yCoordinate,
                                   hit.covarianceTrackingFrame[0], hit.covarianceTrackingFrame[2],
                                   xOverX0, bz, maxChi2, updatedChi2, true)) {
    return false;
  }
  state = updated;
  chi2 = updatedChi2;
  return true;
}

/// buildCellSeed<Tag> equivalence coverage (D007, Architecture.md Sec 10).
/// These reference helpers independently re-transcribe the pre-existing
/// legacy formulas (TrackerTraits::computeLayerCells' barrel branch and
/// detail::mftFwdFitCellClusters) directly from clusters/hits, so a
/// transcription mistake in either the operation or the test would show up
/// as a mismatch rather than being masked by sharing code.

o2::its::Cluster makeGlobalCluster(float x, float y, float z, int id = 0)
{
  return o2::its::Cluster{x, y, z, id};
}

template <TransitionPolicyTag Tag>
void checkSearchWindowEqual(const TrackletSearchWindow<Tag>& lhs, const TrackletSearchWindow<Tag>& rhs)
{
  BOOST_CHECK_EQUAL(lhs.bins.x, rhs.bins.x);
  BOOST_CHECK_EQUAL(lhs.bins.y, rhs.bins.y);
  BOOST_CHECK_EQUAL(lhs.bins.z, rhs.bins.z);
  BOOST_CHECK_EQUAL(lhs.bins.w, rhs.bins.w);
  if constexpr (Tag == TransitionPolicyTag::CylinderCylinder) {
    BOOST_CHECK_EQUAL(lhs.tanLambda, rhs.tanLambda);
    BOOST_CHECK_EQUAL(lhs.sigmaZ, rhs.sigmaZ);
    BOOST_CHECK_EQUAL(lhs.phiCut, rhs.phiCut);
    BOOST_CHECK_EQUAL(lhs.nSigmaCut, rhs.nSigmaCut);
  } else {
    BOOST_CHECK_EQUAL(lhs.xProj, rhs.xProj);
    BOOST_CHECK_EQUAL(lhs.yProj, rhs.yProj);
    BOOST_CHECK_EQUAL(lhs.sigmaX, rhs.sigmaX);
    BOOST_CHECK_EQUAL(lhs.sigmaY, rhs.sigmaY);
    BOOST_CHECK_EQUAL(lhs.meanDeltaZ, rhs.meanDeltaZ);
    BOOST_CHECK_EQUAL(lhs.nSigmaCut, rhs.nSigmaCut);
  }
}

/// Independent re-transcription of TrackerTraits::computeLayerCells' barrel
/// branch: buildTrackSeed(clusterInner, clusterMiddle, hitOuter, bz) then
/// middle-then-inner rotate/propagateTo/correctForMaterial/update, chi2 cut
/// only on the last (inner) step. `nSteps==1` runs only the middle step, for
/// deriving the inner step's exact marginal chi2 in boundary tests.
bool legacyBarrelCellFit(const o2::its::Cluster& clusterInner, const o2::its::Cluster& clusterMiddle,
                         const o2::its::TrackingFrameInfo& hitMiddle, const o2::its::TrackingFrameInfo& hitInner,
                         const o2::its::TrackingFrameInfo& hitOuter,
                         float xOverX0Middle, float xOverX0Inner, float bz, float maxChi2,
                         o2::track::TrackParCovF& outTrack, float& outChi2, int nSteps = 2)
{
  auto track = o2::its::track::buildTrackSeed(clusterInner, clusterMiddle, hitOuter, bz);
  const std::array<const o2::its::TrackingFrameInfo*, 2> hits{&hitMiddle, &hitInner};
  const std::array<float, 2> x0{xOverX0Middle, xOverX0Inner};
  float chi2 = 0.f;
  for (int step = 0; step < nSteps; ++step) {
    const bool isLast = (nSteps == 2) && (step == 1);
    const auto& h = *hits[step];
    if (!track.rotate(h.alphaTrackingFrame)) {
      return false;
    }
    if (!track.propagateTo(h.xTrackingFrame, bz)) {
      return false;
    }
    if (!track.correctForMaterial(x0[step], x0[step] * o2::its::constants::Radl * o2::its::constants::Rho, true)) {
      return false;
    }
    const float predChi2 = track.getPredictedChi2Quiet(h.positionTrackingFrame, h.covarianceTrackingFrame);
    if (isLast && predChi2 > maxChi2) {
      return false;
    }
    if (!track.o2::track::TrackParCov::update(h.positionTrackingFrame, h.covarianceTrackingFrame)) {
      return false;
    }
    chi2 += predChi2;
  }
  outTrack = track;
  outChi2 = chi2;
  return true;
}

/// Independent re-transcription of detail::mftFwdFitCellClusters, reading
/// its three clusters/hits directly instead of through a TimeFrame, and
/// deliberately NOT applying the geometric road pre-cut (that guard is
/// TrackerTraits-owned, not part of the operation under test).
/// `nSteps` limits how many of the outer/middle/inner attach steps run (used
/// to derive the inner step's exact marginal chi2 in boundary tests); the
/// chi2 cut is only active when all three steps run.
bool legacyDiskCellFit(const o2::its::Cluster& clusterInner, const o2::its::Cluster& clusterMiddle, const o2::its::Cluster& clusterOuter,
                       const o2::its::TrackingFrameInfo& hitOuter, const o2::its::TrackingFrameInfo& hitMiddle, const o2::its::TrackingFrameInfo& hitInner,
                       float xOverX0Outer, float xOverX0Middle, float xOverX0Inner,
                       float bz, float trackletMinPt, float maxChi2,
                       o2::track::TrackParCovFwd& outTrack, float& outChi2, int nSteps = 3)
{
  if (clusterInner.zCoordinate <= clusterOuter.zCoordinate + 1.e-6f) {
    return false;
  }

  const float dxTan = clusterMiddle.xCoordinate - clusterInner.xCoordinate;
  const float dyTan = clusterMiddle.yCoordinate - clusterInner.yCoordinate;
  const float dzTan = clusterMiddle.zCoordinate - clusterInner.zCoordinate;
  const float drTan = std::sqrt(dxTan * dxTan + dyTan * dyTan);
  const float dxPhi = clusterOuter.xCoordinate - clusterInner.xCoordinate;
  const float dyPhi = clusterOuter.yCoordinate - clusterInner.yCoordinate;
  const float dzPhi = clusterOuter.zCoordinate - clusterInner.zCoordinate;
  const float drPhi = std::sqrt(dxPhi * dxPhi + dyPhi * dyPhi);
  if (drTan < 1.e-6f || std::abs(dzTan) < 1.e-6f || drPhi < 1.e-6f || std::abs(dzPhi) < 1.e-6f) {
    return false;
  }

  const float invQPt = (trackletMinPt > 0.f) ? 1.f / trackletMinPt : 0.f;
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

  ROOT::Math::SVector<double, 5> seedParams{clusterOuter.xCoordinate, clusterOuter.yCoordinate, phi, tanl, invQPt};
  ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>> seedCov{};
  seedCov(0, 0) = hitOuter.covarianceTrackingFrame[0] > 0.f ? hitOuter.covarianceTrackingFrame[0] : 1.f;
  seedCov(1, 1) = hitOuter.covarianceTrackingFrame[2] > 0.f ? hitOuter.covarianceTrackingFrame[2] : 1.f;
  seedCov(2, 2) = seedCov(3, 3) = 1.;
  const double qptSigma = std::clamp(static_cast<double>(std::abs(invQPt)), 1., 10.);
  seedCov(4, 4) = qptSigma * qptSigma;

  o2::track::TrackParCovFwd track{clusterOuter.zCoordinate, seedParams, seedCov, 0.};
  float chi2 = 0.f;
  const std::array<const o2::its::TrackingFrameInfo*, 3> hits{&hitOuter, &hitMiddle, &hitInner};
  const std::array<float, 3> x0{xOverX0Outer, xOverX0Middle, xOverX0Inner};
  for (int step = 0; step < nSteps; ++step) {
    const bool checkLast = (nSteps == 3) && (step == 2);
    const auto& h = *hits[step];
    if (!detail::mftFwdAttachCluster(track, h.zCoordinate, h.xCoordinate, h.yCoordinate,
                                     h.covarianceTrackingFrame[0], h.covarianceTrackingFrame[2],
                                     x0[step], bz, maxChi2, chi2, checkLast)) {
      return false;
    }
  }
  outTrack = track;
  outChi2 = chi2;
  return true;
}

NominalSurfaceMaterial toMaterial(float xOverX0)
{
  return NominalSurfaceMaterial{xOverX0, xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho};
}

std::array<NominalSurfaceMaterial, 3> toMaterial(const std::array<float, 3>& xOverX0)
{
  return {toMaterial(xOverX0[0]), toMaterial(xOverX0[1]), toMaterial(xOverX0[2])};
}

std::vector<NominalSurfaceMaterial> toMaterial(const std::vector<float>& xOverX0)
{
  std::vector<NominalSurfaceMaterial> material;
  material.reserve(xOverX0.size());
  for (const float x0 : xOverX0) {
    material.push_back(toMaterial(x0));
  }
  return material;
}

void checkBarrelStateEqual(const o2::track::TrackParCovF& lhs, const o2::track::TrackParCovF& rhs)
{
  BOOST_CHECK_EQUAL(lhs.getX(), rhs.getX());
  BOOST_CHECK_EQUAL(lhs.getAlpha(), rhs.getAlpha());
  BOOST_CHECK_EQUAL(lhs.getY(), rhs.getY());
  BOOST_CHECK_EQUAL(lhs.getZ(), rhs.getZ());
  BOOST_CHECK_EQUAL(lhs.getSnp(), rhs.getSnp());
  BOOST_CHECK_EQUAL(lhs.getTgl(), rhs.getTgl());
  BOOST_CHECK_EQUAL(lhs.getQ2Pt(), rhs.getQ2Pt());
  const auto& lhsCov = lhs.getCov();
  const auto& rhsCov = rhs.getCov();
  for (size_t element = 0; element < lhsCov.size(); ++element) {
    BOOST_CHECK_EQUAL(lhsCov[element], rhsCov[element]);
  }
}

void checkDiskStateEqual(const o2::track::TrackParCovFwd& lhs, const o2::track::TrackParCovFwd& rhs)
{
  BOOST_CHECK_EQUAL(lhs.getZ(), rhs.getZ());
  BOOST_CHECK_EQUAL(lhs.getX(), rhs.getX());
  BOOST_CHECK_EQUAL(lhs.getY(), rhs.getY());
  BOOST_CHECK_EQUAL(lhs.getPhi(), rhs.getPhi());
  BOOST_CHECK_EQUAL(lhs.getTanl(), rhs.getTanl());
  BOOST_CHECK_EQUAL(lhs.getInvQPt(), rhs.getInvQPt());
  BOOST_CHECK_EQUAL(lhs.getTrackChi2(), rhs.getTrackChi2());
  const auto& lhsCov = lhs.getCovariances();
  const auto& rhsCov = rhs.getCovariances();
  for (size_t row = 0; row < 5; ++row) {
    for (size_t column = 0; column < 5; ++column) {
      BOOST_CHECK_EQUAL(lhsCov(row, column), rhsCov(row, column));
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(BarrelAcceptsAtOrBelowTheExactChi2Threshold)
{
  const auto currentState = makeBarrelState(5.f, 0.3f, 0.1f, 1.f, 0.05f, 0.4f, 0.2f);
  const auto nextState = makeBarrelState(4.f, 0.3f, 0.2f, 1.1f, 0.06f, 0.42f, 0.21f);
  const auto currentCell = makeBarrelCell(currentState, 10, 20, 30);
  const auto nextCell = makeBarrelCell(nextState, 20, 30, 40);

  // Reference chi2, computed with the same unmodified rotate/propagateTo/
  // getPredictedChi2 primitives the operation itself uses on a local copy.
  auto reference = nextCell;
  BOOST_REQUIRE(reference.rotate(currentCell.getAlpha()));
  BOOST_REQUIRE(reference.propagateTo(currentCell.getX(), Bz));
  const float refChi2 = currentCell.getPredictedChi2(reference);
  BOOST_REQUIRE_GT(refChi2, 0.f);

  CylinderCylinderPolicyParams accept;
  accept.maxChi2ClusterAttachment = refChi2; // exact threshold: "<=" must accept
  BOOST_CHECK(cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(currentCell, nextCell, Bz, accept));

  CylinderCylinderPolicyParams reject;
  reject.maxChi2ClusterAttachment = refChi2 * 0.999f; // just under the exact chi2: must reject
  BOOST_CHECK(!cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(currentCell, nextCell, Bz, reject));
}

BOOST_AUTO_TEST_CASE(BarrelRotationFailureIsRejectedNotThrown)
{
  // nextCell's frame is far enough from currentCell's that rotating it into
  // currentCell's alpha fails the local cos(phi) >= 0 precondition
  // (TrackParametrization::rotateParam); the operation must report this as
  // "not compatible", not crash or silently ignore the failure.
  const auto currentState = makeBarrelState(5.f, 0.f, 0.1f, 1.f, 0.f, 0.4f, 0.2f);
  const auto nextState = makeBarrelState(5.f, 3.f, 0.1f, 1.f, 0.f, 0.4f, 0.2f);
  const auto currentCell = makeBarrelCell(currentState, 10, 20, 30);
  const auto nextCell = makeBarrelCell(nextState, 20, 30, 40);

  auto reference = nextCell;
  BOOST_REQUIRE(!reference.rotate(currentCell.getAlpha()));

  CylinderCylinderPolicyParams params; // permissive threshold: only rotation should matter
  params.maxChi2ClusterAttachment = 1.e6f;
  BOOST_CHECK(!cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(currentCell, nextCell, Bz, params));
}

BOOST_AUTO_TEST_CASE(BarrelPropagationFailureIsRejectedNotThrown)
{
  // Distinct from the rotation-failure case above: here rotate() succeeds
  // (both cells share alpha, so the local cos(phi)>=0 precondition trivially
  // holds) and propagateTo() itself is the one that fails.
  //
  // TrackParametrization::propagateParamTo(xk, b) (the scalar-field overload
  // CellSeedTpl::propagateTo resolves to) fails exactly when the requested
  // step would push the local sin(phi) parameter (snp) outside (-1, 1):
  // with curvature crv = q2pt * b * B2C and step x2r = crv * dx (dx = xk -
  // x_current), it requires |snp + x2r| <= Almost1 (~0.999999). Choosing an
  // extreme q2pt (very low pt) for `nextCell` and a several-cm dx makes x2r
  // exceed that bound by a comfortable margin regardless of B2C's exact
  // value, without needing any other precondition to be near a boundary.
  const auto currentState = makeBarrelState(5.f, 0.3f, 0.1f, 1.f, 0.f, 0.4f, 0.2f);
  const auto nextState = makeBarrelState(0.f, 0.3f, 0.1f, 1.f, 0.f, 0.4f, 2000.f); // same alpha as currentState; extreme q2pt
  const auto currentCell = makeBarrelCell(currentState, 10, 20, 30);
  const auto nextCell = makeBarrelCell(nextState, 20, 30, 40);

  auto reference = nextCell;
  BOOST_REQUIRE(reference.rotate(currentCell.getAlpha()));
  BOOST_REQUIRE(!reference.propagateTo(currentCell.getX(), Bz));

  CylinderCylinderPolicyParams params; // permissive threshold: only propagation should matter
  params.maxChi2ClusterAttachment = 1.e6f;
  BOOST_CHECK(!cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(currentCell, nextCell, Bz, params));
}

BOOST_AUTO_TEST_CASE(BarrelInputCellsAreNotMutated)
{
  const auto currentState = makeBarrelState(5.f, 0.3f, 0.1f, 1.f, 0.05f, 0.4f, 0.2f);
  const auto nextState = makeBarrelState(4.f, 0.3f, 0.2f, 1.1f, 0.06f, 0.42f, 0.21f);
  const auto currentCell = makeBarrelCell(currentState, 10, 20, 30);
  const auto nextCell = makeBarrelCell(nextState, 20, 30, 40);
  const auto currentBefore = currentCell;
  const auto nextBefore = nextCell;

  CylinderCylinderPolicyParams params;
  params.maxChi2ClusterAttachment = 1.e6f;
  cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(currentCell, nextCell, Bz, params);

  BOOST_CHECK_EQUAL(currentCell.getX(), currentBefore.getX());
  BOOST_CHECK_EQUAL(currentCell.getAlpha(), currentBefore.getAlpha());
  BOOST_CHECK_EQUAL(currentCell.getY(), currentBefore.getY());
  BOOST_CHECK_EQUAL(currentCell.getZ(), currentBefore.getZ());
  BOOST_CHECK_EQUAL(currentCell.getFirstClusterIndex(), currentBefore.getFirstClusterIndex());

  BOOST_CHECK_EQUAL(nextCell.getX(), nextBefore.getX());
  BOOST_CHECK_EQUAL(nextCell.getAlpha(), nextBefore.getAlpha());
  BOOST_CHECK_EQUAL(nextCell.getY(), nextBefore.getY());
  BOOST_CHECK_EQUAL(nextCell.getZ(), nextBefore.getZ());
  BOOST_CHECK_EQUAL(nextCell.getFirstClusterIndex(), nextBefore.getFirstClusterIndex());
}

BOOST_AUTO_TEST_CASE(DiskAcceptsAtOrBelowTheExactChi2Threshold)
{
  const auto currentState = makeForwardState(0.f, 0.f, 0.f, 0.f, -1.2f, 0.1f);
  const auto nextState = makeForwardState(-1.f, 0.05f, -0.03f, 0.02f, -1.19f, 0.11f);
  const auto currentCell = makeForwardCell(currentState, 10, 20, 30);
  const auto nextCell = makeForwardCell(nextState, 20, 30, 40); // firstClusterIndex==20==currentCell.getSecondClusterIndex()

  auto reference = nextState;
  detail::mftFwdPropagateToZ(reference, static_cast<float>(currentState.getZ()), Bz);
  const float refChi2 = detail::mftFwdStateChi2(currentState, reference);
  BOOST_REQUIRE_GT(refChi2, 0.f);

  DiskDiskPolicyParams accept;
  accept.maxChi2ClusterAttachment = refChi2;
  BOOST_CHECK(cellsAreCompatible<TransitionPolicyTag::DiskDisk>(currentCell, nextCell, Bz, accept));
  BOOST_CHECK(detail::mftFwdCellsAreCompatible(currentCell, nextCell, Bz, accept.maxChi2ClusterAttachment));

  DiskDiskPolicyParams reject;
  reject.maxChi2ClusterAttachment = refChi2 * 0.999f;
  BOOST_CHECK(!cellsAreCompatible<TransitionPolicyTag::DiskDisk>(currentCell, nextCell, Bz, reject));
  BOOST_CHECK(!detail::mftFwdCellsAreCompatible(currentCell, nextCell, Bz, reject.maxChi2ClusterAttachment));
}

BOOST_AUTO_TEST_CASE(DiskRejectsDiscontinuousClusterIndices)
{
  // mftFwdCellsAreCompatible's own precondition: nextCell must continue
  // currentCell's second cluster. A generous chi2 bound must not override it.
  const auto currentState = makeForwardState(0.f, 0.f, 0.f, 0.f, -1.2f, 0.1f);
  const auto nextState = makeForwardState(-1.f, 0.f, 0.f, 0.f, -1.2f, 0.1f); // identical state, chi2 ~ 0
  const auto currentCell = makeForwardCell(currentState, 10, 20, 30);
  const auto nextCell = makeForwardCell(nextState, 999, 30, 40); // firstClusterIndex != 20

  DiskDiskPolicyParams params;
  params.maxChi2ClusterAttachment = 1.e6f;
  BOOST_CHECK(!cellsAreCompatible<TransitionPolicyTag::DiskDisk>(currentCell, nextCell, Bz, params));
}

BOOST_AUTO_TEST_CASE(DiskInputCellsAreNotMutated)
{
  const auto currentState = makeForwardState(0.f, 0.f, 0.f, 0.f, -1.2f, 0.1f);
  const auto nextState = makeForwardState(-1.f, 0.05f, -0.03f, 0.02f, -1.19f, 0.11f);
  const auto currentCell = makeForwardCell(currentState, 10, 20, 30);
  const auto nextCell = makeForwardCell(nextState, 20, 30, 40);
  const auto currentBefore = currentCell;
  const auto nextBefore = nextCell;

  DiskDiskPolicyParams params;
  params.maxChi2ClusterAttachment = 1.e6f;
  cellsAreCompatible<TransitionPolicyTag::DiskDisk>(currentCell, nextCell, Bz, params);

  BOOST_CHECK_EQUAL(currentCell.getX(), currentBefore.getX());
  BOOST_CHECK_EQUAL(currentCell.getZ(), currentBefore.getZ());
  BOOST_CHECK_EQUAL(currentCell.getFirstClusterIndex(), currentBefore.getFirstClusterIndex());

  BOOST_CHECK_EQUAL(nextCell.getX(), nextBefore.getX());
  BOOST_CHECK_EQUAL(nextCell.getZ(), nextBefore.getZ());
  BOOST_CHECK_EQUAL(nextCell.getFirstClusterIndex(), nextBefore.getFirstClusterIndex());
}

BOOST_AUTO_TEST_CASE(BindingCopiesEveryFieldToTheCorrectSlot)
{
  // Distinct sentinel per field so a field-swap bug in the binding is caught.
  TrackingParameters legacy;
  legacy.TrackletMinPt = 1.11f;
  legacy.CellDeltaTanLambdaSigma = 2.22f;
  legacy.NSigmaCut = 3.33f;
  legacy.MaxChi2ClusterAttachment = 4.44f;
  legacy.MaxChi2NDF = 5.55f;
  legacy.PVres = 8.88f;
  legacy.CellRoadRCut = 6.66f;
  legacy.TrackletMinAbsX = 7.77f;
  legacy.LayerxX0 = {0.011f, 0.022f, 0.033f};
  legacy.CorrType = o2::base::PropagatorF::MatCorrType::USEMatCorrLUT;

  const auto barrel = bindTransitionPolicyParams<TransitionPolicyTag::CylinderCylinder>(legacy);
  BOOST_CHECK_CLOSE(barrel.trackletMinPt, 1.11f, 1e-6);
  BOOST_CHECK_CLOSE(barrel.cellDeltaTanLambdaSigma, 2.22f, 1e-6);
  BOOST_CHECK_CLOSE(barrel.nSigmaCut, 3.33f, 1e-6);
  BOOST_CHECK_CLOSE(barrel.maxChi2ClusterAttachment, 4.44f, 1e-6);
  BOOST_CHECK_CLOSE(barrel.maxChi2NDF, 5.55f, 1e-6);
  BOOST_CHECK_CLOSE(barrel.pvResolution, 8.88f, 1e-6);
  BOOST_CHECK(barrel.isValid());

  const auto disk = bindTransitionPolicyParams<TransitionPolicyTag::DiskDisk>(legacy);
  BOOST_CHECK_CLOSE(disk.trackletMinPt, 1.11f, 1e-6);
  BOOST_CHECK_CLOSE(disk.cellDeltaTanLambdaSigma, 2.22f, 1e-6);
  BOOST_CHECK_CLOSE(disk.cellRoadRCut, 6.66f, 1e-6);
  BOOST_CHECK_CLOSE(disk.trackletMinAbsX, 7.77f, 1e-6);
  BOOST_CHECK_CLOSE(disk.nSigmaCut, 3.33f, 1e-6);
  BOOST_CHECK_CLOSE(disk.maxChi2ClusterAttachment, 4.44f, 1e-6);
  BOOST_CHECK_CLOSE(disk.maxChi2NDF, 5.55f, 1e-6);
  BOOST_CHECK(disk.isValid());

  const auto legacyMaterial = toMaterial(legacy.LayerxX0);
  const auto attach = bindAttachHitPolicyConfig(gsl::span<const NominalSurfaceMaterial>(legacyMaterial), legacy);
  BOOST_REQUIRE_EQUAL(attach.layerMaterial.size(), 3u);
  BOOST_CHECK_CLOSE(attach.layerMaterial[0].xOverX0, 0.011f, 1e-6);
  BOOST_CHECK_CLOSE(attach.layerMaterial[1].xOverX0, 0.022f, 1e-6);
  BOOST_CHECK_CLOSE(attach.layerMaterial[2].xOverX0, 0.033f, 1e-6);
  BOOST_CHECK(attach.corrType == o2::base::PropagatorF::MatCorrType::USEMatCorrLUT);
  BOOST_CHECK(attach.isValid(3));
  BOOST_CHECK(!attach.isValid(4));
}

BOOST_AUTO_TEST_CASE(BoundNonFiniteParametersAreDetectableThroughIsValid)
{
  const float nan = std::numeric_limits<float>::quiet_NaN();
  const float inf = std::numeric_limits<float>::infinity();

  TrackingParameters legacy;
  legacy.TrackletMinPt = 1.11f;
  legacy.CellDeltaTanLambdaSigma = 2.22f;
  legacy.NSigmaCut = 3.33f;
  legacy.MaxChi2ClusterAttachment = 4.44f;
  legacy.MaxChi2NDF = 5.55f;
  legacy.CellRoadRCut = 6.66f;
  legacy.TrackletMinAbsX = 7.77f;

  auto barrelHealthy = legacy;
  BOOST_CHECK(bindTransitionPolicyParams<TransitionPolicyTag::CylinderCylinder>(barrelHealthy).isValid());
  auto barrelNaN = legacy;
  barrelNaN.MaxChi2ClusterAttachment = nan;
  BOOST_CHECK(!bindTransitionPolicyParams<TransitionPolicyTag::CylinderCylinder>(barrelNaN).isValid());
  auto barrelInf = legacy;
  barrelInf.NSigmaCut = inf;
  BOOST_CHECK(!bindTransitionPolicyParams<TransitionPolicyTag::CylinderCylinder>(barrelInf).isValid());

  auto diskHealthy = legacy;
  BOOST_CHECK(bindTransitionPolicyParams<TransitionPolicyTag::DiskDisk>(diskHealthy).isValid());
  auto diskNaN = legacy;
  diskNaN.CellRoadRCut = nan;
  BOOST_CHECK(!bindTransitionPolicyParams<TransitionPolicyTag::DiskDisk>(diskNaN).isValid());
  auto diskInf = legacy;
  diskInf.TrackletMinAbsX = inf;
  BOOST_CHECK(!bindTransitionPolicyParams<TransitionPolicyTag::DiskDisk>(diskInf).isValid());

  auto materialNaN = legacy;
  materialNaN.LayerxX0[2] = nan;
  const auto materialNaNSpan = toMaterial(materialNaN.LayerxX0);
  BOOST_CHECK(!bindAttachHitPolicyConfig(gsl::span<const NominalSurfaceMaterial>(materialNaNSpan), materialNaN)
                 .isValid(materialNaN.LayerxX0.size()));
  auto materialInf = legacy;
  materialInf.LayerxX0[4] = inf;
  const auto materialInfSpan = toMaterial(materialInf.LayerxX0);
  BOOST_CHECK(!bindAttachHitPolicyConfig(gsl::span<const NominalSurfaceMaterial>(materialInfSpan), materialInf)
                 .isValid(materialInf.LayerxX0.size()));
  auto invalidCorrection = legacy;
  invalidCorrection.CorrType = static_cast<o2::base::PropagatorF::MatCorrType>(99);
  const auto invalidCorrectionMaterial = toMaterial(invalidCorrection.LayerxX0);
  BOOST_CHECK(!bindAttachHitPolicyConfig(gsl::span<const NominalSurfaceMaterial>(invalidCorrectionMaterial), invalidCorrection)
                 .isValid(invalidCorrection.LayerxX0.size()));
}

BOOST_AUTO_TEST_CASE(CylinderProjectSearchWindowMatchesInlineFormulaAndDirectPhiZBins)
{
  TrackingParameters legacy;
  legacy.PVres = 0.f; // valid: disables only the primary-vertex resolution term
  const auto params = bindTransitionPolicyParams<TransitionPolicyTag::CylinderCylinder>(legacy);
  BOOST_REQUIRE(params.isValid());

  IndexTableUtils<7> indexUtils;
  indexUtils.setTrackingParameters(legacy);

  const auto source = makeGlobalCluster(2.f, 0.f, 0.5f);
  const auto sourceHit = makeBarrelHit(2.f, 0.f, 0.f, 0.5f);
  const auto vertex = makeVertex(0.f, 0.f, 0.f, 1.e-4f, 1.e-4f, 4.e-4f, 4);
  const TrackletProjectionState<TransitionPolicyTag::CylinderCylinder> state{
    0, 3, 2.f, 3.8f, 4.2f, 5.e-4f, 2.e-3f, 0.08f};

  TrackletSearchWindow<TransitionPolicyTag::CylinderCylinder> window{};
  BOOST_REQUIRE((projectSearchWindow<TransitionPolicyTag::CylinderCylinder, 7>(
    source, sourceHit, vertex, state, Bz, indexUtils, params, window)));

  const float inverseR0 = 1.f / source.radius;
  const float resolution = o2::gpu::CAMath::Sqrt(o2::its::math_utils::Sq(state.sourcePositionResolution) +
                                                 o2::its::math_utils::Sq(params.pvResolution) / float(vertex.getNContributors()));
  const float tanLambda = (source.zCoordinate - vertex.getZ()) * inverseR0;
  const float zAtTargetMinR = tanLambda * (state.targetMinR - source.radius) + source.zCoordinate;
  const float zAtTargetMaxR = tanLambda * (state.targetMaxR - source.radius) + source.zCoordinate;
  const float sqInvDeltaZ0 = 1.f / (o2::its::math_utils::Sq(source.zCoordinate - vertex.getZ()) + o2::its::constants::Tolerance);
  const float sigmaZ = o2::gpu::CAMath::Sqrt((o2::its::math_utils::Sq(resolution) * o2::its::math_utils::Sq(tanLambda) *
                                              ((o2::its::math_utils::Sq(inverseR0) + sqInvDeltaZ0) * o2::its::math_utils::Sq(state.meanDeltaR) + 1.f)) +
                                             o2::its::math_utils::Sq(state.meanDeltaR * state.transitionMSAngle));
  const auto directBins = getBinsPhiZ(source.phi, state.toLayer, zAtTargetMinR, zAtTargetMaxR,
                                      sigmaZ * params.nSigmaCut, state.transitionPhiCut, indexUtils);

  BOOST_CHECK_EQUAL(window.bins.x, directBins.x);
  BOOST_CHECK_EQUAL(window.bins.y, directBins.y);
  BOOST_CHECK_EQUAL(window.bins.z, directBins.z);
  BOOST_CHECK_EQUAL(window.bins.w, directBins.w);
  BOOST_CHECK_EQUAL(window.tanLambda, tanLambda);
  BOOST_CHECK_EQUAL(window.sigmaZ, sigmaZ);

  legacy.PVres = 0.025f;
  const auto positivePVParams = bindTransitionPolicyParams<TransitionPolicyTag::CylinderCylinder>(legacy);
  BOOST_REQUIRE(positivePVParams.isValid());
  TrackletSearchWindow<TransitionPolicyTag::CylinderCylinder> positivePVWindow{};
  BOOST_REQUIRE((projectSearchWindow<TransitionPolicyTag::CylinderCylinder, 7>(
    source, sourceHit, vertex, state, Bz, indexUtils, positivePVParams, positivePVWindow)));
  const float positivePVResolution = o2::gpu::CAMath::Sqrt(o2::its::math_utils::Sq(state.sourcePositionResolution) +
                                                           o2::its::math_utils::Sq(positivePVParams.pvResolution) / float(vertex.getNContributors()));
  const float positivePVSigmaZ = o2::gpu::CAMath::Sqrt((o2::its::math_utils::Sq(positivePVResolution) * o2::its::math_utils::Sq(tanLambda) *
                                                        ((o2::its::math_utils::Sq(inverseR0) + sqInvDeltaZ0) * o2::its::math_utils::Sq(state.meanDeltaR) + 1.f)) +
                                                       o2::its::math_utils::Sq(state.meanDeltaR * state.transitionMSAngle));
  BOOST_CHECK_EQUAL(positivePVWindow.sigmaZ, positivePVSigmaZ);
  BOOST_CHECK_GT(positivePVWindow.sigmaZ, window.sigmaZ);

  const float targetRadius = 4.f;
  const float targetZ = tanLambda * (targetRadius - source.radius) + source.zCoordinate;
  const auto acceptedTarget = makeGlobalCluster(targetRadius, 0.f, targetZ);
  float acceptedTanLambda = -999.f;
  BOOST_CHECK(window.acceptCandidate(source, acceptedTarget, acceptedTanLambda));
  BOOST_CHECK_EQUAL(acceptedTanLambda, (source.zCoordinate - acceptedTarget.zCoordinate) / (source.radius - acceptedTarget.radius));

  const auto rejectedTarget = makeGlobalCluster(-targetRadius, 0.f, targetZ);
  float rejectedTanLambda = 123.f;
  BOOST_CHECK(!window.acceptCandidate(source, rejectedTarget, rejectedTanLambda));
  BOOST_CHECK_EQUAL(rejectedTanLambda, 123.f);
}

BOOST_AUTO_TEST_CASE(DiskProjectSearchWindowReusesHelpersAndDirectProjectedXYBins)
{
  TrackingParameters legacy;
  const auto params = bindTransitionPolicyParams<TransitionPolicyTag::DiskDisk>(legacy);
  BOOST_REQUIRE(params.isValid());

  IndexTableUtils<10> indexUtils;
  std::array<float, 10> halfExtents{};
  halfExtents.fill(20.f);
  indexUtils.setIndexTableParams(IndexTableCoordType::XY, legacy.RowBins, legacy.ColBins, -20.f, 20.f, halfExtents);

  constexpr int fromLayer = 1;
  constexpr int toLayer = 4; // deliberately skipped/nonadjacent transition
  const float fromZ = detail::mftLayerZ(fromLayer);
  const float toZ = detail::mftLayerZ(toLayer);
  const auto source = makeGlobalCluster(1.2f, 0.7f, fromZ);
  const auto sourceHit = makeDiskHit(fromZ, source.xCoordinate, source.yCoordinate, 2.e-4f, 3.e-4f);
  const auto vertex = makeVertex(0.01f, -0.02f, 0.1f, 4.e-4f, 5.e-4f, 0.04f, 3);
  const TrackletProjectionState<TransitionPolicyTag::DiskDisk> state{
    fromLayer, toLayer, fromZ, toZ, toZ - fromZ, 2.f, 3.e-3f, 0.04f};

  TrackletSearchWindow<TransitionPolicyTag::DiskDisk> window{};
  BOOST_REQUIRE((projectSearchWindow<TransitionPolicyTag::DiskDisk, 10>(
    source, sourceHit, vertex, state, Bz, indexUtils, params, window)));

  float expectedX = 0.f;
  float expectedY = 0.f;
  detail::mftTrackletProject(source.xCoordinate, source.yCoordinate, source.zCoordinate,
                             vertex.getX(), vertex.getY(), vertex.getZ(),
                             fromLayer, toLayer, Bz, params.trackletMinPt, expectedX, expectedY);
  float expectedSigmaX = 0.f;
  float expectedSigmaY = 0.f;
  detail::mftTrackletSigmaXY(source.xCoordinate, source.yCoordinate,
                             vertex.getX(), vertex.getY(), vertex.getZ(),
                             sourceHit.covarianceTrackingFrame[0], sourceHit.covarianceTrackingFrame[2],
                             vertex.getSigmaX2(), vertex.getSigmaY2(), vertex.getSigmaZ2(),
                             fromLayer, toLayer, state.sourceReferenceRadius, state.meanDeltaZ,
                             state.transitionMSAngle, state.transitionBendingAngle,
                             expectedX, expectedY, expectedSigmaX, expectedSigmaY);

  const float zSpread = params.nSigmaCut * vertex.getSigmaZ();
  const float zVtxMin = vertex.getZ() - zSpread;
  const float zVtxMax = vertex.getZ() + zSpread;
  const float absZFrom = std::abs(fromZ);
  const float absZTo = std::abs(toZ);
  const float denomMin = zVtxMax + absZFrom;
  const float denomMax = absZFrom + zVtxMin;
  float radialRangeMin = (std::abs(denomMin) > 1.e-6f) ? source.radius * (zVtxMax + absZTo) / denomMin : source.radius;
  float radialRangeMax = (std::abs(denomMax) > 1.e-6f) ? source.radius * (absZTo + zVtxMin) / denomMax : source.radius;
  if (radialRangeMin > radialRangeMax) {
    std::swap(radialRangeMin, radialRangeMax);
  }
  const auto directBins = getBinsRectClusterAtProj<10>(expectedX, expectedY, toLayer,
                                                       radialRangeMin, radialRangeMax,
                                                       expectedSigmaX * params.nSigmaCut,
                                                       expectedSigmaY * params.nSigmaCut,
                                                       indexUtils);

  BOOST_CHECK_EQUAL(window.bins.x, directBins.x);
  BOOST_CHECK_EQUAL(window.bins.y, directBins.y);
  BOOST_CHECK_EQUAL(window.bins.z, directBins.z);
  BOOST_CHECK_EQUAL(window.bins.w, directBins.w);
  BOOST_CHECK_EQUAL(window.xProj, expectedX);
  BOOST_CHECK_EQUAL(window.yProj, expectedY);
  BOOST_CHECK_EQUAL(window.sigmaX, expectedSigmaX);
  BOOST_CHECK_EQUAL(window.sigmaY, expectedSigmaY);

  const auto acceptedTarget = makeGlobalCluster(expectedX, expectedY, toZ);
  float acceptedTanLambda = -999.f;
  BOOST_CHECK(window.acceptCandidate(source, acceptedTarget, acceptedTanLambda));
  BOOST_CHECK_EQUAL(acceptedTanLambda, (source.zCoordinate - acceptedTarget.zCoordinate) / state.meanDeltaZ);

  auto rejectedWindow = window;
  rejectedWindow.meanDeltaZ = 0.f;
  float rejectedTanLambda = 123.f;
  BOOST_CHECK(!rejectedWindow.acceptCandidate(source, acceptedTarget, rejectedTanLambda));
  BOOST_CHECK_EQUAL(rejectedTanLambda, 123.f);
}

BOOST_AUTO_TEST_CASE(ProjectSearchWindowInvalidBinsLeaveEveryOutputFieldUnchanged)
{
  TrackingParameters legacy;

  IndexTableUtils<7> cylinderIndexUtils;
  cylinderIndexUtils.setTrackingParameters(legacy);
  const auto cylinderParams = bindTransitionPolicyParams<TransitionPolicyTag::CylinderCylinder>(legacy);
  const auto cylinderSource = makeGlobalCluster(2.f, 0.f, 100.f);
  const auto cylinderHit = makeBarrelHit(2.f, 0.f, 0.f, 100.f);
  const auto cylinderVertex = makeVertex(0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
  const TrackletProjectionState<TransitionPolicyTag::CylinderCylinder> cylinderState{
    0, 3, 2.f, 3.8f, 4.2f, 5.e-4f, 2.e-3f, 0.08f};
  const TrackletSearchWindow<TransitionPolicyTag::CylinderCylinder> cylinderSentinel{
    {101, 102, 103, 104}, 105.f, 106.f, 107.f, 108.f};
  auto cylinderOut = cylinderSentinel;
  BOOST_CHECK(!(projectSearchWindow<TransitionPolicyTag::CylinderCylinder, 7>(
    cylinderSource, cylinderHit, cylinderVertex, cylinderState, Bz, cylinderIndexUtils, cylinderParams, cylinderOut)));
  checkSearchWindowEqual(cylinderOut, cylinderSentinel);

  IndexTableUtils<10> diskIndexUtils;
  std::array<float, 10> tinyHalfExtents{};
  tinyHalfExtents.fill(0.01f);
  diskIndexUtils.setIndexTableParams(IndexTableCoordType::XY, legacy.RowBins, legacy.ColBins, -0.01f, 0.01f, tinyHalfExtents);
  const auto diskParams = bindTransitionPolicyParams<TransitionPolicyTag::DiskDisk>(legacy);
  constexpr int fromLayer = 0;
  constexpr int toLayer = 1;
  const float fromZ = detail::mftLayerZ(fromLayer);
  const float toZ = detail::mftLayerZ(toLayer);
  const auto diskSource = makeGlobalCluster(1.f, 0.5f, fromZ);
  const auto diskHit = makeDiskHit(fromZ, diskSource.xCoordinate, diskSource.yCoordinate);
  const auto diskVertex = makeVertex(0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
  const TrackletProjectionState<TransitionPolicyTag::DiskDisk> diskState{
    fromLayer, toLayer, fromZ, toZ, toZ - fromZ, 2.f, 3.e-3f, 0.04f};
  const TrackletSearchWindow<TransitionPolicyTag::DiskDisk> diskSentinel{
    {201, 202, 203, 204}, 205.f, 206.f, 207.f, 208.f, 209.f, 210.f};
  auto diskOut = diskSentinel;
  BOOST_CHECK(!(projectSearchWindow<TransitionPolicyTag::DiskDisk, 10>(
    diskSource, diskHit, diskVertex, diskState, Bz, diskIndexUtils, diskParams, diskOut)));
  checkSearchWindowEqual(diskOut, diskSentinel);
}

BOOST_AUTO_TEST_CASE(CylinderCandidateUsesPeriodicPhiAndStrictSigmaAndPhiBoundaries)
{
  TrackletSearchWindow<TransitionPolicyTag::CylinderCylinder> window{
    {}, 0.f, 1.f, 0.02f, 5.f};

  const float wrapEpsilon = 0.005f;
  const auto wrappedSource = makeGlobalCluster(std::cos(wrapEpsilon), -std::sin(wrapEpsilon), 0.f);
  const auto wrappedTarget = makeGlobalCluster(2.f * std::cos(wrapEpsilon), 2.f * std::sin(wrapEpsilon), 0.f);
  float tanLambda = -9.f;
  BOOST_REQUIRE(o2::its::math_utils::isPhiDifferenceBelow(wrappedSource.phi, wrappedTarget.phi, window.phiCut));
  BOOST_CHECK(window.acceptCandidate(wrappedSource, wrappedTarget, tanLambda));

  const auto source = makeGlobalCluster(1.f, 0.f, 0.f);
  const auto exactSigmaTarget = makeGlobalCluster(2.f, 0.f, 5.f);
  tanLambda = 71.f;
  BOOST_CHECK(!window.acceptCandidate(source, exactSigmaTarget, tanLambda));
  BOOST_CHECK_EQUAL(tanLambda, 71.f);
  const auto insideSigmaTarget = makeGlobalCluster(2.f, 0.f, std::nextafter(5.f, 0.f));
  BOOST_CHECK(window.acceptCandidate(source, insideSigmaTarget, tanLambda));

  const auto phiTarget = makeGlobalCluster(2.f * std::cos(0.125f), 2.f * std::sin(0.125f), 0.f);
  window.phiCut = o2::gpu::CAMath::Abs(source.phi - phiTarget.phi);
  const bool directPhiDecision = o2::its::math_utils::isPhiDifferenceBelow(source.phi, phiTarget.phi, window.phiCut);
  tanLambda = 72.f;
  BOOST_CHECK_EQUAL(window.acceptCandidate(source, phiTarget, tanLambda), directPhiDecision);
  BOOST_CHECK(!directPhiDecision);
  BOOST_CHECK_EQUAL(tanLambda, 72.f);
}

BOOST_AUTO_TEST_CASE(DiskProjectionCoversStraightLineNearZeroDenominatorAndRadialSwap)
{
  TrackingParameters legacy;
  const auto params = bindTransitionPolicyParams<TransitionPolicyTag::DiskDisk>(legacy);
  constexpr int fromLayer = 0;
  constexpr int toLayer = 1;
  const float fromZ = detail::mftLayerZ(fromLayer);
  const float toZ = detail::mftLayerZ(toLayer);
  const auto source = makeGlobalCluster(1.f, 0.5f, fromZ);
  const auto sourceHit = makeDiskHit(fromZ, source.xCoordinate, source.yCoordinate);
  const TrackletProjectionState<TransitionPolicyTag::DiskDisk> state{
    fromLayer, toLayer, fromZ, toZ, toZ - fromZ, 2.f, 3.e-3f, 0.04f};

  IndexTableUtils<10> indexUtils;
  std::array<float, 10> halfExtents{};
  halfExtents.fill(200.f);
  indexUtils.setIndexTableParams(IndexTableCoordType::XY, legacy.RowBins, legacy.ColBins, -200.f, 200.f, halfExtents);

  const auto straightVertex = makeVertex(0.1f, -0.2f, 0.3f, 4.e-4f, 5.e-4f, 0.04f);
  TrackletSearchWindow<TransitionPolicyTag::DiskDisk> straightWindow{};
  BOOST_REQUIRE((projectSearchWindow<TransitionPolicyTag::DiskDisk, 10>(
    source, sourceHit, straightVertex, state, 0.f, indexUtils, params, straightWindow)));
  float expectedX = 0.f;
  float expectedY = 0.f;
  detail::mftTrackletProject(source.xCoordinate, source.yCoordinate, source.zCoordinate,
                             straightVertex.getX(), straightVertex.getY(), straightVertex.getZ(),
                             fromLayer, toLayer, 0.f, params.trackletMinPt, expectedX, expectedY);
  BOOST_CHECK_EQUAL(straightWindow.xProj, expectedX);
  BOOST_CHECK_EQUAL(straightWindow.yProj, expectedY);

  const auto fallbackVertex = makeVertex(0.1f, -0.2f, fromZ, 4.e-4f, 5.e-4f, 0.f);
  TrackletSearchWindow<TransitionPolicyTag::DiskDisk> fallbackWindow{};
  BOOST_REQUIRE((projectSearchWindow<TransitionPolicyTag::DiskDisk, 10>(
    source, sourceHit, fallbackVertex, state, 0.f, indexUtils, params, fallbackWindow)));
  expectedX = expectedY = 0.f;
  detail::mftTrackletProject(source.xCoordinate, source.yCoordinate, source.zCoordinate,
                             fallbackVertex.getX(), fallbackVertex.getY(), fallbackVertex.getZ(),
                             fromLayer, toLayer, 0.f, params.trackletMinPt, expectedX, expectedY);
  BOOST_CHECK_EQUAL(expectedX, source.xCoordinate);
  BOOST_CHECK_EQUAL(expectedY, source.yCoordinate);
  BOOST_CHECK_EQUAL(fallbackWindow.xProj, expectedX);
  BOOST_CHECK_EQUAL(fallbackWindow.yProj, expectedY);

  const auto swapVertex = makeVertex(0.f, 0.f, fromZ + 0.1f, 0.f, 0.f, 0.01f);
  const float zSpread = params.nSigmaCut * swapVertex.getSigmaZ();
  const float zVtxMin = swapVertex.getZ() - zSpread;
  const float zVtxMax = swapVertex.getZ() + zSpread;
  const float absZFrom = std::abs(fromZ);
  const float absZTo = std::abs(toZ);
  const float rawRadialMin = source.radius * (zVtxMax + absZTo) / (zVtxMax + absZFrom);
  const float rawRadialMax = source.radius * (absZTo + zVtxMin) / (absZFrom + zVtxMin);
  BOOST_REQUIRE_GT(rawRadialMin, rawRadialMax);
  TrackletSearchWindow<TransitionPolicyTag::DiskDisk> swapWindow{};
  BOOST_REQUIRE((projectSearchWindow<TransitionPolicyTag::DiskDisk, 10>(
    source, sourceHit, swapVertex, state, 0.f, indexUtils, params, swapWindow)));
  expectedX = expectedY = 0.f;
  detail::mftTrackletProject(source.xCoordinate, source.yCoordinate, source.zCoordinate,
                             swapVertex.getX(), swapVertex.getY(), swapVertex.getZ(),
                             fromLayer, toLayer, 0.f, params.trackletMinPt, expectedX, expectedY);
  float expectedSigmaX = 0.f;
  float expectedSigmaY = 0.f;
  detail::mftTrackletSigmaXY(source.xCoordinate, source.yCoordinate,
                             swapVertex.getX(), swapVertex.getY(), swapVertex.getZ(),
                             sourceHit.covarianceTrackingFrame[0], sourceHit.covarianceTrackingFrame[2],
                             swapVertex.getSigmaX2(), swapVertex.getSigmaY2(), swapVertex.getSigmaZ2(),
                             fromLayer, toLayer, state.sourceReferenceRadius, state.meanDeltaZ,
                             state.transitionMSAngle, state.transitionBendingAngle,
                             expectedX, expectedY, expectedSigmaX, expectedSigmaY);
  const auto directBins = getBinsRectClusterAtProj<10>(expectedX, expectedY, toLayer,
                                                       rawRadialMax, rawRadialMin,
                                                       expectedSigmaX * params.nSigmaCut,
                                                       expectedSigmaY * params.nSigmaCut,
                                                       indexUtils);
  BOOST_CHECK_EQUAL(swapWindow.bins.x, directBins.x);
  BOOST_CHECK_EQUAL(swapWindow.bins.y, directBins.y);
  BOOST_CHECK_EQUAL(swapWindow.bins.z, directBins.z);
  BOOST_CHECK_EQUAL(swapWindow.bins.w, directBins.w);
}

BOOST_AUTO_TEST_CASE(DiskCandidatePreservesInverseVarianceAndStrictBoundarySemantics)
{
  const auto source = makeGlobalCluster(1.f, 0.f, -45.f);
  const auto distantTarget = makeGlobalCluster(100.f, -80.f, -47.f);
  TrackletSearchWindow<TransitionPolicyTag::DiskDisk> zeroSigmaWindow{
    {}, 0.f, 0.f, 0.f, -1.f, 2.f, 5.f};
  float tanLambda = -8.f;
  BOOST_CHECK(zeroSigmaWindow.acceptCandidate(source, distantTarget, tanLambda));
  BOOST_CHECK_EQUAL(tanLambda, 1.f);

  TrackletSearchWindow<TransitionPolicyTag::DiskDisk> chi2Window{
    {}, 0.f, 0.f, 1.f, 1.f, 2.f, 5.f};
  const auto exactChi2Target = makeGlobalCluster(5.f, 0.f, -47.f);
  tanLambda = 81.f;
  BOOST_CHECK(!chi2Window.acceptCandidate(source, exactChi2Target, tanLambda));
  BOOST_CHECK_EQUAL(tanLambda, 81.f);
  const auto insideChi2Target = makeGlobalCluster(std::nextafter(5.f, 0.f), 0.f, -47.f);
  BOOST_CHECK(chi2Window.acceptCandidate(source, insideChi2Target, tanLambda));

  const auto centeredTarget = makeGlobalCluster(0.f, 0.f, -47.f);
  chi2Window.meanDeltaZ = 1.e-6f;
  tanLambda = 82.f;
  BOOST_CHECK(!chi2Window.acceptCandidate(source, centeredTarget, tanLambda));
  BOOST_CHECK_EQUAL(tanLambda, 82.f);
  chi2Window.meanDeltaZ = std::nextafter(1.e-6f, std::numeric_limits<float>::infinity());
  BOOST_CHECK(chi2Window.acceptCandidate(source, centeredTarget, tanLambda));
}

BOOST_AUTO_TEST_CASE(BarrelAttachHitExactBoundaryAndInlineEquivalence)
{
  const auto initial = makeBarrelState(4.f, 0.3f, 0.2f, 1.1f, 0.06f, 0.42f, 0.21f);
  const auto hit = makeBarrelHit(5.f, 0.3f, 0.35f, 1.25f);
  auto reference = initial;
  BOOST_REQUIRE(reference.rotate(hit.alphaTrackingFrame));
  BOOST_REQUIRE(o2::base::Propagator::Instance()->propagateToX(reference, hit.xTrackingFrame, Bz,
                                                               o2::base::PropagatorImpl<float>::MAX_SIN_PHI,
                                                               o2::base::PropagatorImpl<float>::MAX_STEP,
                                                               NoMaterialCorrection));
  BOOST_REQUIRE(reference.correctForMaterial(0.f, 0.f, true));
  const float exactChi2 = reference.getPredictedChi2Quiet(hit.positionTrackingFrame, hit.covarianceTrackingFrame);
  BOOST_REQUIRE_GT(exactChi2, 0.f);

  CylinderCylinderPolicyParams exact;
  exact.maxChi2ClusterAttachment = exactChi2;
  auto policyState = initial;
  auto inlineState = initial;
  float policyChi2 = 7.25f;
  float inlineChi2 = policyChi2;
  BOOST_CHECK(attachHit<TransitionPolicyTag::CylinderCylinder>(policyState, hit, toMaterial(0.f), NoMaterialCorrection, Bz, policyChi2, exact));
  BOOST_CHECK(legacyBarrelAttach(inlineState, hit, 0.f, NoMaterialCorrection, Bz, exact.maxChi2ClusterAttachment, inlineChi2));
  checkBarrelStateEqual(policyState, inlineState);
  BOOST_CHECK_EQUAL(policyChi2, inlineChi2);

  CylinderCylinderPolicyParams below = exact;
  below.maxChi2ClusterAttachment = std::nextafter(exactChi2, 0.f);
  policyState = initial;
  inlineState = initial;
  policyChi2 = inlineChi2 = 7.25f;
  BOOST_CHECK(!attachHit<TransitionPolicyTag::CylinderCylinder>(policyState, hit, toMaterial(0.f), NoMaterialCorrection, Bz, policyChi2, below));
  BOOST_CHECK(!legacyBarrelAttach(inlineState, hit, 0.f, NoMaterialCorrection, Bz, below.maxChi2ClusterAttachment, inlineChi2));
  checkBarrelStateEqual(policyState, inlineState);
  BOOST_CHECK_EQUAL(policyChi2, inlineChi2);
  BOOST_CHECK_NE(policyState.getX(), initial.getX()); // legacy barrel failure retains successful propagation
}

BOOST_AUTO_TEST_CASE(BarrelAttachHitPreservesRotationPropagationAndUpdateFailures)
{
  CylinderCylinderPolicyParams params;
  params.maxChi2ClusterAttachment = std::numeric_limits<float>::max();

  auto rotationState = makeBarrelState(5.f, 3.f, 0.1f, 1.f, 0.f, 0.4f, 0.2f);
  const auto rotationBefore = rotationState;
  const auto rotationHit = makeBarrelHit(5.f, 0.f, 0.1f, 1.f);
  float chi2 = 3.f;
  BOOST_CHECK(!attachHit<TransitionPolicyTag::CylinderCylinder>(rotationState, rotationHit, toMaterial(0.f), NoMaterialCorrection, Bz, chi2, params));
  checkBarrelStateEqual(rotationState, rotationBefore);
  BOOST_CHECK_EQUAL(chi2, 3.f);

  auto propagationState = makeBarrelState(0.f, 0.3f, 0.1f, 1.f, 0.f, 0.4f, 2000.f);
  const auto propagationHit = makeBarrelHit(5.f, 0.3f, 0.1f, 1.f);
  chi2 = 4.f;
  BOOST_CHECK(!attachHit<TransitionPolicyTag::CylinderCylinder>(propagationState, propagationHit, toMaterial(0.f), NoMaterialCorrection, Bz, chi2, params));
  BOOST_CHECK_EQUAL(chi2, 4.f);

  const o2::track::TrackParCovF::params_t zeroPar{0.f, 0.f, 0.f, 0.4f, 0.2f};
  o2::track::TrackParCovF::covMat_t singularCov{};
  singularCov[5] = 1.e-4f;
  singularCov[9] = 1.e-4f;
  singularCov[14] = 1.e-2f;
  auto updateState = o2::track::TrackParCovF{5.f, 0.3f, zeroPar, singularCov};
  auto inlineUpdateState = updateState;
  const auto singularHit = makeBarrelHit(5.f, 0.3f, 0.f, 0.f, 0.f, 0.f);
  chi2 = 5.f;
  float inlineChi2 = chi2;
  BOOST_CHECK(!attachHit<TransitionPolicyTag::CylinderCylinder>(updateState, singularHit, toMaterial(0.f), NoMaterialCorrection, Bz, chi2, params));
  BOOST_CHECK(!legacyBarrelAttach(inlineUpdateState, singularHit, 0.f, NoMaterialCorrection, Bz,
                                  params.maxChi2ClusterAttachment, inlineChi2));
  checkBarrelStateEqual(updateState, inlineUpdateState);
  BOOST_CHECK_EQUAL(chi2, inlineChi2);
  BOOST_CHECK_GT(chi2, 5.f); // legacy increments chi2 before the failed update
}

BOOST_AUTO_TEST_CASE(DiskAttachHitExactBoundaryFailureImmutabilityAndInlineEquivalence)
{
  const auto initial = makeForwardState(-1.f, 0.05f, -0.03f, 0.02f, -1.19f, 0.11f);
  const auto hit = makeDiskHit(0.f, 0.2f, -0.1f);
  auto propagated = initial;
  detail::mftFwdPropagateToZ(propagated, hit.zCoordinate, Bz);
  const float exactChi2 = detail::mftFwdPredictedChi2(propagated, hit.xCoordinate, hit.yCoordinate,
                                                      hit.covarianceTrackingFrame[0], hit.covarianceTrackingFrame[2]);
  BOOST_REQUIRE_GT(exactChi2, 0.f);

  DiskDiskPolicyParams exact;
  exact.maxChi2ClusterAttachment = exactChi2;
  auto policyState = initial;
  auto inlineState = initial;
  float policyChi2 = 8.5f;
  float inlineChi2 = policyChi2;
  BOOST_CHECK(attachHit<TransitionPolicyTag::DiskDisk>(policyState, hit, toMaterial(0.017f), NoMaterialCorrection, Bz, policyChi2, exact));
  BOOST_CHECK(legacyDiskAttach(inlineState, hit, 0.017f, Bz, exact.maxChi2ClusterAttachment, inlineChi2));
  checkDiskStateEqual(policyState, inlineState);
  BOOST_CHECK_EQUAL(policyChi2, inlineChi2);
  BOOST_CHECK_NE(policyState.getZ(), initial.getZ());

  DiskDiskPolicyParams below = exact;
  below.maxChi2ClusterAttachment = std::nextafter(exactChi2, 0.f);
  policyState = initial;
  inlineState = initial;
  const auto before = policyState;
  policyChi2 = inlineChi2 = 8.5f;
  BOOST_CHECK(!attachHit<TransitionPolicyTag::DiskDisk>(policyState, hit, toMaterial(0.017f), NoMaterialCorrection, Bz, policyChi2, below));
  BOOST_CHECK(!legacyDiskAttach(inlineState, hit, 0.017f, Bz, below.maxChi2ClusterAttachment, inlineChi2));
  checkDiskStateEqual(policyState, before);
  checkDiskStateEqual(policyState, inlineState);
  BOOST_CHECK_EQUAL(policyChi2, 8.5f);
  BOOST_CHECK_EQUAL(policyChi2, inlineChi2);
}

BOOST_AUTO_TEST_CASE(DiskAttachHitIgnoresArealDensityMCSOnly)
{
  // attachHit<DiskDisk> passes only material.xOverX0 to the unchanged
  // legacy MCS-only MFT path; arealDensityGPerCm2 must never reach it, so two
  // material values sharing xOverX0 but wildly differing in
  // arealDensityGPerCm2 must produce byte-identical state and chi2.
  const auto initial = makeForwardState(-1.f, 0.05f, -0.03f, 0.02f, -1.19f, 0.11f);
  const auto hit = makeDiskHit(0.f, 0.2f, -0.1f);

  DiskDiskPolicyParams params;
  params.maxChi2ClusterAttachment = 1.e6f;

  auto lowDensityState = initial;
  float lowDensityChi2 = 0.f;
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::DiskDisk>(lowDensityState, hit, NominalSurfaceMaterial{0.017f, 0.001f},
                                                         NoMaterialCorrection, Bz, lowDensityChi2, params));

  auto highDensityState = initial;
  float highDensityChi2 = 0.f;
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::DiskDisk>(highDensityState, hit, NominalSurfaceMaterial{0.017f, 999.f},
                                                         NoMaterialCorrection, Bz, highDensityChi2, params));

  checkDiskStateEqual(lowDensityState, highDensityState);
  BOOST_CHECK_EQUAL(lowDensityChi2, highDensityChi2);
}

// ---------------------------------------------------------------------------
// buildCellSeed<Tag>
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(CylinderBuildCellSeedAcceptsAtExactChi2ThresholdAndMatchesInlineEquivalence)
{
  // Nearly-collinear points: a realistic high-pt (gentle curvature) barrel
  // trajectory, so rotate/propagateTo/correctForMaterial/update all succeed
  // on both steps.
  const auto clusterInner = makeGlobalCluster(3.0f, 0.100f, 0.9f, 10);
  const auto clusterMiddle = makeGlobalCluster(4.0f, 0.150f, 1.05f, 20);
  const auto clusterOuter = makeGlobalCluster(5.0f, 0.201f, 1.25f, 30);
  const auto hitOuter = makeBarrelHit(5.f, 0.f, 0.201f, 1.25f);
  const auto hitMiddle = makeBarrelHit(4.f, 0.f, 0.150f, 1.05f);
  const auto hitInner = makeBarrelHit(3.f, 0.f, 0.100f, 0.9f);
  const std::array<float, 3> xOverX0{0.005f, 0.005f, 0.005f}; // inner, middle, outer
  const auto material = toMaterial(xOverX0);

  // Independently derive the exact marginal chi2 the inner (last) step
  // contributes: run only the middle step generously, then perform the
  // inner step's rotate/propagateTo/correctForMaterial without updating.
  o2::track::TrackParCovF afterMiddle;
  float chi2AfterMiddle = 0.f;
  BOOST_REQUIRE(legacyBarrelCellFit(clusterInner, clusterMiddle, hitMiddle, hitInner, hitOuter,
                                    xOverX0[1], xOverX0[0], Bz, 1.e6f, afterMiddle, chi2AfterMiddle, 1));
  auto atInner = afterMiddle;
  BOOST_REQUIRE(atInner.rotate(hitInner.alphaTrackingFrame));
  BOOST_REQUIRE(atInner.propagateTo(hitInner.xTrackingFrame, Bz));
  BOOST_REQUIRE(atInner.correctForMaterial(xOverX0[0], xOverX0[0] * o2::its::constants::Radl * o2::its::constants::Rho, true));
  const float exactChi2 = atInner.getPredictedChi2Quiet(hitInner.positionTrackingFrame, hitInner.covarianceTrackingFrame);
  BOOST_REQUIRE_GT(exactChi2, 0.f);

  CylinderCylinderPolicyParams accept;
  accept.maxChi2ClusterAttachment = exactChi2;

  o2::track::TrackParCovF policyState{};
  float policyChi2 = 6.25f;
  BOOST_CHECK(buildCellSeed<TransitionPolicyTag::CylinderCylinder>(clusterInner, clusterMiddle, clusterOuter,
                                                                   hitInner, hitMiddle, hitOuter, material, Bz,
                                                                   policyState, policyChi2, accept));

  o2::track::TrackParCovF inlineState{};
  float inlineChi2 = 0.f;
  BOOST_CHECK(legacyBarrelCellFit(clusterInner, clusterMiddle, hitMiddle, hitInner, hitOuter,
                                  xOverX0[1], xOverX0[0], Bz, accept.maxChi2ClusterAttachment, inlineState, inlineChi2));
  checkBarrelStateEqual(policyState, inlineState);
  BOOST_CHECK_EQUAL(policyChi2, inlineChi2);

  CylinderCylinderPolicyParams reject = accept;
  reject.maxChi2ClusterAttachment = std::nextafter(exactChi2, 0.f);
  policyState = o2::track::TrackParCovF{};
  policyChi2 = 6.25f;
  BOOST_CHECK(!buildCellSeed<TransitionPolicyTag::CylinderCylinder>(clusterInner, clusterMiddle, clusterOuter,
                                                                    hitInner, hitMiddle, hitOuter, material, Bz,
                                                                    policyState, policyChi2, reject));
  checkBarrelStateEqual(policyState, o2::track::TrackParCovF{});
  BOOST_CHECK_EQUAL(policyChi2, 6.25f);
}

BOOST_AUTO_TEST_CASE(CylinderBuildCellSeedRotationFailureLeavesOutputUnchanged)
{
  const auto clusterInner = makeGlobalCluster(3.0f, 0.05f, 0.9f, 10);
  const auto clusterMiddle = makeGlobalCluster(4.0f, 0.15f, 1.05f, 20);
  const auto clusterOuter = makeGlobalCluster(5.0f, 0.30f, 1.25f, 30);
  const auto hitOuter = makeBarrelHit(5.f, 0.f, 0.30f, 1.25f);
  const auto hitMiddle = makeBarrelHit(4.f, 3.f, 0.15f, 1.05f); // far alpha: rotate must fail
  const auto hitInner = makeBarrelHit(3.f, 0.f, 0.05f, 0.9f);
  const std::array<float, 3> xOverX0{0.005f, 0.005f, 0.005f};
  const auto material = toMaterial(xOverX0);

  auto seed = o2::its::track::buildTrackSeed(clusterInner, clusterMiddle, hitOuter, Bz);
  BOOST_REQUIRE(!seed.rotate(hitMiddle.alphaTrackingFrame));

  CylinderCylinderPolicyParams params;
  params.maxChi2ClusterAttachment = 1.e6f;
  auto outState = makeBarrelState(9.f, 9.f, 9.f, 9.f, 9.f, 9.f, 9.f);
  const auto before = outState;
  float chi2 = 12.5f;
  BOOST_CHECK(!buildCellSeed<TransitionPolicyTag::CylinderCylinder>(clusterInner, clusterMiddle, clusterOuter,
                                                                    hitInner, hitMiddle, hitOuter, material, Bz,
                                                                    outState, chi2, params));
  checkBarrelStateEqual(outState, before);
  BOOST_CHECK_EQUAL(chi2, 12.5f);
}

BOOST_AUTO_TEST_CASE(CylinderBuildCellSeedPropagationFailureLeavesOutputUnchanged)
{
  const auto clusterInner = makeGlobalCluster(3.0f, 0.05f, 0.9f, 10);
  const auto clusterMiddle = makeGlobalCluster(4.0f, 0.15f, 1.05f, 20);
  const auto clusterOuter = makeGlobalCluster(5.0f, 0.30f, 1.25f, 30);
  const auto hitOuter = makeBarrelHit(5.f, 0.f, 0.30f, 1.25f);
  const auto hitMiddle = makeBarrelHit(500.f, 0.f, 0.15f, 1.05f); // huge step: propagateTo must fail
  const auto hitInner = makeBarrelHit(3.f, 0.f, 0.05f, 0.9f);
  const std::array<float, 3> xOverX0{0.005f, 0.005f, 0.005f};
  const auto material = toMaterial(xOverX0);

  auto seed = o2::its::track::buildTrackSeed(clusterInner, clusterMiddle, hitOuter, Bz);
  BOOST_REQUIRE(seed.rotate(hitMiddle.alphaTrackingFrame));
  BOOST_REQUIRE(!seed.propagateTo(hitMiddle.xTrackingFrame, Bz));

  CylinderCylinderPolicyParams params;
  params.maxChi2ClusterAttachment = 1.e6f;
  auto outState = makeBarrelState(9.f, 9.f, 9.f, 9.f, 9.f, 9.f, 9.f);
  const auto before = outState;
  float chi2 = 7.5f;
  BOOST_CHECK(!buildCellSeed<TransitionPolicyTag::CylinderCylinder>(clusterInner, clusterMiddle, clusterOuter,
                                                                    hitInner, hitMiddle, hitOuter, material, Bz,
                                                                    outState, chi2, params));
  checkBarrelStateEqual(outState, before);
  BOOST_CHECK_EQUAL(chi2, 7.5f);
}

BOOST_AUTO_TEST_CASE(CylinderBuildCellSeedUpdateFailureLeavesOutputUnchanged)
{
  // Force a singular innovation covariance on the first (middle) attach
  // step: the seed's Y/Z covariance comes directly from hitOuter's
  // covariance (o2::its::track::buildTrackSeed), and a zero-distance
  // propagation (hitMiddle at the seed's own xTrackingFrame) leaves it
  // untouched; pairing a zero-covariance outer seed with a zero-covariance
  // middle hit makes the combined innovation covariance exactly singular, so
  // update() must fail -- same mechanism as the existing singular-covariance
  // attachHit coverage.
  const auto clusterInner = makeGlobalCluster(3.0f, 0.05f, 0.9f, 10);
  const auto clusterMiddle = makeGlobalCluster(4.0f, 0.15f, 1.05f, 20);
  const auto clusterOuter = makeGlobalCluster(5.0f, 0.30f, 1.25f, 30);
  const auto hitOuter = makeBarrelHit(5.f, 0.f, 0.30f, 1.25f, 0.f, 0.f);
  const auto hitMiddle = makeBarrelHit(5.f, 0.f, 0.15f, 1.05f, 0.f, 0.f);
  const auto hitInner = makeBarrelHit(3.f, 0.f, 0.05f, 0.9f);
  const std::array<float, 3> xOverX0{0.005f, 0.f, 0.005f};
  const auto material = toMaterial(xOverX0);

  auto seed = o2::its::track::buildTrackSeed(clusterInner, clusterMiddle, hitOuter, Bz);
  BOOST_REQUIRE(seed.rotate(hitMiddle.alphaTrackingFrame));
  BOOST_REQUIRE(seed.propagateTo(hitMiddle.xTrackingFrame, Bz));
  BOOST_REQUIRE(seed.correctForMaterial(0.f, 0.f, true));
  BOOST_REQUIRE(!seed.o2::track::TrackParCov::update(hitMiddle.positionTrackingFrame, hitMiddle.covarianceTrackingFrame));

  CylinderCylinderPolicyParams params;
  params.maxChi2ClusterAttachment = 1.e6f;
  auto outState = makeBarrelState(9.f, 9.f, 9.f, 9.f, 9.f, 9.f, 9.f);
  const auto before = outState;
  float chi2 = 3.5f;
  BOOST_CHECK(!buildCellSeed<TransitionPolicyTag::CylinderCylinder>(clusterInner, clusterMiddle, clusterOuter,
                                                                    hitInner, hitMiddle, hitOuter, material, Bz,
                                                                    outState, chi2, params));
  checkBarrelStateEqual(outState, before);
  BOOST_CHECK_EQUAL(chi2, 3.5f);
}

BOOST_AUTO_TEST_CASE(CylinderBuildCellSeedInputsAreNotMutated)
{
  const auto clusterInner = makeGlobalCluster(3.0f, 0.100f, 0.9f, 10);
  const auto clusterMiddle = makeGlobalCluster(4.0f, 0.150f, 1.05f, 20);
  const auto clusterOuter = makeGlobalCluster(5.0f, 0.201f, 1.25f, 30);
  const auto hitOuter = makeBarrelHit(5.f, 0.f, 0.201f, 1.25f);
  const auto hitMiddle = makeBarrelHit(4.f, 0.f, 0.150f, 1.05f);
  const auto hitInner = makeBarrelHit(3.f, 0.f, 0.100f, 0.9f);
  const std::array<float, 3> xOverX0{0.005f, 0.005f, 0.005f};
  const auto material = toMaterial(xOverX0);
  const auto clusterInnerBefore = clusterInner;
  const auto clusterMiddleBefore = clusterMiddle;
  const auto clusterOuterBefore = clusterOuter;
  const auto hitOuterBefore = hitOuter;
  const auto hitMiddleBefore = hitMiddle;
  const auto hitInnerBefore = hitInner;

  CylinderCylinderPolicyParams params;
  params.maxChi2ClusterAttachment = 1.e6f;
  o2::track::TrackParCovF outState{};
  float chi2 = 0.f;
  BOOST_CHECK(buildCellSeed<TransitionPolicyTag::CylinderCylinder>(clusterInner, clusterMiddle, clusterOuter,
                                                                   hitInner, hitMiddle, hitOuter, material, Bz,
                                                                   outState, chi2, params));

  BOOST_CHECK_EQUAL(clusterInner.xCoordinate, clusterInnerBefore.xCoordinate);
  BOOST_CHECK_EQUAL(clusterInner.yCoordinate, clusterInnerBefore.yCoordinate);
  BOOST_CHECK_EQUAL(clusterInner.zCoordinate, clusterInnerBefore.zCoordinate);
  BOOST_CHECK_EQUAL(clusterMiddle.xCoordinate, clusterMiddleBefore.xCoordinate);
  BOOST_CHECK_EQUAL(clusterMiddle.yCoordinate, clusterMiddleBefore.yCoordinate);
  BOOST_CHECK_EQUAL(clusterMiddle.zCoordinate, clusterMiddleBefore.zCoordinate);
  BOOST_CHECK_EQUAL(clusterOuter.xCoordinate, clusterOuterBefore.xCoordinate);
  BOOST_CHECK_EQUAL(clusterOuter.yCoordinate, clusterOuterBefore.yCoordinate);
  BOOST_CHECK_EQUAL(clusterOuter.zCoordinate, clusterOuterBefore.zCoordinate);
  BOOST_CHECK_EQUAL(hitOuter.xTrackingFrame, hitOuterBefore.xTrackingFrame);
  BOOST_CHECK_EQUAL(hitOuter.alphaTrackingFrame, hitOuterBefore.alphaTrackingFrame);
  BOOST_CHECK_EQUAL(hitMiddle.xTrackingFrame, hitMiddleBefore.xTrackingFrame);
  BOOST_CHECK_EQUAL(hitMiddle.alphaTrackingFrame, hitMiddleBefore.alphaTrackingFrame);
  BOOST_CHECK_EQUAL(hitInner.xTrackingFrame, hitInnerBefore.xTrackingFrame);
  BOOST_CHECK_EQUAL(hitInner.alphaTrackingFrame, hitInnerBefore.alphaTrackingFrame);
}

BOOST_AUTO_TEST_CASE(CylinderBuildCellSeedUsesMaterialSlotsOneThenZero)
{
  const auto clusterInner = makeGlobalCluster(3.0f, 0.100f, 0.9f, 10);
  const auto clusterMiddle = makeGlobalCluster(4.0f, 0.150f, 1.05f, 20);
  const auto clusterOuter = makeGlobalCluster(5.0f, 0.201f, 1.25f, 30);
  const auto hitOuter = makeBarrelHit(5.f, 0.f, 0.201f, 1.25f);
  const auto hitMiddle = makeBarrelHit(4.f, 0.f, 0.150f, 1.05f);
  const auto hitInner = makeBarrelHit(3.f, 0.f, 0.100f, 0.9f);
  // Distinct values in the two real slots; NaN in the unused outer slot
  // proves it is never read (a real read would poison state/chi2 with NaN
  // and desynchronize them from the reference, which never touches index 2).
  const std::array<float, 3> xOverX0{0.05f, 0.005f, std::numeric_limits<float>::quiet_NaN()};
  const auto material = toMaterial(xOverX0);

  CylinderCylinderPolicyParams params;
  params.maxChi2ClusterAttachment = 1.e6f;

  o2::track::TrackParCovF policyState{};
  float policyChi2 = 0.f;
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::CylinderCylinder>(clusterInner, clusterMiddle, clusterOuter,
                                                                     hitInner, hitMiddle, hitOuter, material, Bz,
                                                                     policyState, policyChi2, params));
  BOOST_CHECK(std::isfinite(policyChi2));

  o2::track::TrackParCovF referenceState{};
  float referenceChi2 = 0.f;
  BOOST_REQUIRE(legacyBarrelCellFit(clusterInner, clusterMiddle, hitMiddle, hitInner, hitOuter,
                                    xOverX0[1], xOverX0[0], Bz, params.maxChi2ClusterAttachment,
                                    referenceState, referenceChi2));
  checkBarrelStateEqual(policyState, referenceState);
  BOOST_CHECK_EQUAL(policyChi2, referenceChi2);
}

BOOST_AUTO_TEST_CASE(CylinderBuildCellSeedRepeatedCallsAreDeterministic)
{
  const auto clusterInner = makeGlobalCluster(3.0f, 0.100f, 0.9f, 10);
  const auto clusterMiddle = makeGlobalCluster(4.0f, 0.150f, 1.05f, 20);
  const auto clusterOuter = makeGlobalCluster(5.0f, 0.201f, 1.25f, 30);
  const auto hitOuter = makeBarrelHit(5.f, 0.f, 0.201f, 1.25f);
  const auto hitMiddle = makeBarrelHit(4.f, 0.f, 0.150f, 1.05f);
  const auto hitInner = makeBarrelHit(3.f, 0.f, 0.100f, 0.9f);
  const std::array<float, 3> xOverX0{0.005f, 0.005f, 0.005f};
  const auto material = toMaterial(xOverX0);

  CylinderCylinderPolicyParams params;
  params.maxChi2ClusterAttachment = 1.e6f;

  o2::track::TrackParCovF firstState{};
  float firstChi2 = 0.f;
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::CylinderCylinder>(clusterInner, clusterMiddle, clusterOuter,
                                                                     hitInner, hitMiddle, hitOuter, material, Bz,
                                                                     firstState, firstChi2, params));

  o2::track::TrackParCovF secondState{};
  float secondChi2 = 0.f;
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::CylinderCylinder>(clusterInner, clusterMiddle, clusterOuter,
                                                                     hitInner, hitMiddle, hitOuter, material, Bz,
                                                                     secondState, secondChi2, params));

  checkBarrelStateEqual(firstState, secondState);
  BOOST_CHECK_EQUAL(firstChi2, secondChi2);
}

BOOST_AUTO_TEST_CASE(DiskBuildCellSeedAcceptsAtExactChi2ThresholdAndMatchesInlineEquivalence)
{
  const auto clusterInner = makeGlobalCluster(1.0f, 0.5f, -0.4f, 10);
  const auto clusterMiddle = makeGlobalCluster(1.3f, 0.62f, -0.6f, 20);
  const auto clusterOuter = makeGlobalCluster(1.7f, 0.78f, -0.9f, 30);
  const auto hitInner = makeDiskHit(-0.4f, 1.0f, 0.5f);
  const auto hitMiddle = makeDiskHit(-0.6f, 1.3f, 0.62f);
  const auto hitOuter = makeDiskHit(-0.9f, 1.7f, 0.78f);
  const std::array<float, 3> xOverX0{0.015f, 0.017f, 0.02f}; // inner, middle, outer
  const auto material = toMaterial(xOverX0);

  const float trackletMinPt = 0.3f;

  // Independently derive the exact marginal chi2 the inner (last) step
  // contributes: run only the outer+middle steps generously, then propagate
  // to the inner hit and read its predicted chi2 without updating.
  o2::track::TrackParCovFwd afterMiddle;
  float chi2AfterMiddle = 0.f;
  BOOST_REQUIRE(legacyDiskCellFit(clusterInner, clusterMiddle, clusterOuter, hitOuter, hitMiddle, hitInner,
                                  xOverX0[2], xOverX0[1], xOverX0[0], Bz, trackletMinPt, 1.e6f,
                                  afterMiddle, chi2AfterMiddle, 2));
  auto atInner = afterMiddle;
  detail::mftFwdPropagateToZ(atInner, hitInner.zCoordinate, Bz);
  const float exactChi2 = detail::mftFwdPredictedChi2(atInner, hitInner.xCoordinate, hitInner.yCoordinate,
                                                      hitInner.covarianceTrackingFrame[0], hitInner.covarianceTrackingFrame[2]);
  BOOST_REQUIRE_GT(exactChi2, 0.f);

  DiskDiskPolicyParams accept;
  accept.trackletMinPt = trackletMinPt;
  accept.maxChi2ClusterAttachment = exactChi2;

  o2::track::TrackParCovFwd policyState{};
  float policyChi2 = 9.25f;
  BOOST_CHECK(buildCellSeed<TransitionPolicyTag::DiskDisk>(clusterInner, clusterMiddle, clusterOuter,
                                                           hitInner, hitMiddle, hitOuter, material, Bz,
                                                           policyState, policyChi2, accept));

  o2::track::TrackParCovFwd inlineState{};
  float inlineChi2 = 0.f;
  BOOST_CHECK(legacyDiskCellFit(clusterInner, clusterMiddle, clusterOuter, hitOuter, hitMiddle, hitInner,
                                xOverX0[2], xOverX0[1], xOverX0[0], Bz, accept.trackletMinPt,
                                accept.maxChi2ClusterAttachment, inlineState, inlineChi2));
  checkDiskStateEqual(policyState, inlineState);
  BOOST_CHECK_EQUAL(policyChi2, inlineChi2);

  DiskDiskPolicyParams reject = accept;
  reject.maxChi2ClusterAttachment = std::nextafter(exactChi2, 0.f);
  policyState = o2::track::TrackParCovFwd{};
  policyChi2 = 9.25f;
  BOOST_CHECK(!buildCellSeed<TransitionPolicyTag::DiskDisk>(clusterInner, clusterMiddle, clusterOuter,
                                                            hitInner, hitMiddle, hitOuter, material, Bz,
                                                            policyState, policyChi2, reject));
  checkDiskStateEqual(policyState, o2::track::TrackParCovFwd{});
  BOOST_CHECK_EQUAL(policyChi2, 9.25f);
}

BOOST_AUTO_TEST_CASE(DiskBuildCellSeedZOrderingRejectionLeavesOutputUnchanged)
{
  // clusterInner.z is not strictly greater than clusterOuter.z: the legacy
  // z-ordering precondition (detail::mftFwdFitCellClusters) must reject.
  const auto clusterInner = makeGlobalCluster(1.0f, 0.5f, -0.9f, 10);
  const auto clusterMiddle = makeGlobalCluster(1.3f, 0.62f, -0.6f, 20);
  const auto clusterOuter = makeGlobalCluster(1.7f, 0.78f, -0.9f, 30);
  const auto hitInner = makeDiskHit(-0.9f, 1.0f, 0.5f);
  const auto hitMiddle = makeDiskHit(-0.6f, 1.3f, 0.62f);
  const auto hitOuter = makeDiskHit(-0.9f, 1.7f, 0.78f);
  const std::array<float, 3> xOverX0{0.015f, 0.017f, 0.02f};
  const auto material = toMaterial(xOverX0);

  DiskDiskPolicyParams params;
  params.trackletMinPt = 0.3f;
  params.maxChi2ClusterAttachment = 1.e6f;
  auto outState = makeForwardState(1.f, 2.f, 3.f, 4.f, 5.f, 6.f);
  const auto before = outState;
  float chi2 = 4.25f;
  BOOST_CHECK(!buildCellSeed<TransitionPolicyTag::DiskDisk>(clusterInner, clusterMiddle, clusterOuter,
                                                            hitInner, hitMiddle, hitOuter, material, Bz,
                                                            outState, chi2, params));
  checkDiskStateEqual(outState, before);
  BOOST_CHECK_EQUAL(chi2, 4.25f);
}

BOOST_AUTO_TEST_CASE(DiskBuildCellSeedDegenerateGeometryRejectionLeavesOutputUnchanged)
{
  // clusterMiddle coincides with clusterInner in x/y: drTan collapses to
  // zero, the legacy degenerate-geometry precondition must reject.
  const auto clusterInner = makeGlobalCluster(1.0f, 0.5f, -0.4f, 10);
  const auto clusterMiddle = makeGlobalCluster(1.0f, 0.5f, -0.6f, 20);
  const auto clusterOuter = makeGlobalCluster(1.7f, 0.78f, -0.9f, 30);
  const auto hitInner = makeDiskHit(-0.4f, 1.0f, 0.5f);
  const auto hitMiddle = makeDiskHit(-0.6f, 1.0f, 0.5f);
  const auto hitOuter = makeDiskHit(-0.9f, 1.7f, 0.78f);
  const std::array<float, 3> xOverX0{0.015f, 0.017f, 0.02f};
  const auto material = toMaterial(xOverX0);

  DiskDiskPolicyParams params;
  params.trackletMinPt = 0.3f;
  params.maxChi2ClusterAttachment = 1.e6f;
  auto outState = makeForwardState(1.f, 2.f, 3.f, 4.f, 5.f, 6.f);
  const auto before = outState;
  float chi2 = 2.75f;
  BOOST_CHECK(!buildCellSeed<TransitionPolicyTag::DiskDisk>(clusterInner, clusterMiddle, clusterOuter,
                                                            hitInner, hitMiddle, hitOuter, material, Bz,
                                                            outState, chi2, params));
  checkDiskStateEqual(outState, before);
  BOOST_CHECK_EQUAL(chi2, 2.75f);
}

BOOST_AUTO_TEST_CASE(DiskBuildCellSeedDoesNotApplyRoadPreCut)
{
  // clusterMiddle is displaced far off the inner-outer seed line --
  // passesCellRoadPrecut<DiskDisk> (CellRoadRCut) would reject this, but
  // that guard is TrackerTraits-owned, not part of buildCellSeed. With a
  // generous chi2 bound and otherwise valid geometry, the fit must still
  // succeed, proving no road pre-cut is applied inside the operation.
  const auto clusterInner = makeGlobalCluster(1.0f, 0.5f, -0.4f, 10);
  const auto clusterMiddle = makeGlobalCluster(5.0f, 4.5f, -0.6f, 20); // far off the seed line
  const auto clusterOuter = makeGlobalCluster(1.7f, 0.78f, -0.9f, 30);
  const auto hitInner = makeDiskHit(-0.4f, 1.0f, 0.5f);
  const auto hitMiddle = makeDiskHit(-0.6f, 5.0f, 4.5f);
  const auto hitOuter = makeDiskHit(-0.9f, 1.7f, 0.78f);
  const std::array<float, 3> xOverX0{0.015f, 0.017f, 0.02f};
  const auto material = toMaterial(xOverX0);

  DiskDiskPolicyParams params;
  params.trackletMinPt = 0.3f;
  params.maxChi2ClusterAttachment = 1.e9f;

  o2::track::TrackParCovFwd outState{};
  float chi2 = 0.f;
  BOOST_CHECK(buildCellSeed<TransitionPolicyTag::DiskDisk>(clusterInner, clusterMiddle, clusterOuter,
                                                           hitInner, hitMiddle, hitOuter, material, Bz,
                                                           outState, chi2, params));
}

BOOST_AUTO_TEST_CASE(DiskBuildCellSeedInputsAreNotMutated)
{
  const auto clusterInner = makeGlobalCluster(1.0f, 0.5f, -0.4f, 10);
  const auto clusterMiddle = makeGlobalCluster(1.3f, 0.62f, -0.6f, 20);
  const auto clusterOuter = makeGlobalCluster(1.7f, 0.78f, -0.9f, 30);
  const auto hitInner = makeDiskHit(-0.4f, 1.0f, 0.5f);
  const auto hitMiddle = makeDiskHit(-0.6f, 1.3f, 0.62f);
  const auto hitOuter = makeDiskHit(-0.9f, 1.7f, 0.78f);
  const std::array<float, 3> xOverX0{0.015f, 0.017f, 0.02f};
  const auto material = toMaterial(xOverX0);
  const auto clusterInnerBefore = clusterInner;
  const auto clusterMiddleBefore = clusterMiddle;
  const auto clusterOuterBefore = clusterOuter;
  const auto hitInnerBefore = hitInner;
  const auto hitMiddleBefore = hitMiddle;
  const auto hitOuterBefore = hitOuter;

  DiskDiskPolicyParams params;
  params.trackletMinPt = 0.3f;
  params.maxChi2ClusterAttachment = 1.e6f;
  o2::track::TrackParCovFwd outState{};
  float chi2 = 0.f;
  BOOST_CHECK(buildCellSeed<TransitionPolicyTag::DiskDisk>(clusterInner, clusterMiddle, clusterOuter,
                                                           hitInner, hitMiddle, hitOuter, material, Bz,
                                                           outState, chi2, params));

  BOOST_CHECK_EQUAL(clusterInner.xCoordinate, clusterInnerBefore.xCoordinate);
  BOOST_CHECK_EQUAL(clusterInner.zCoordinate, clusterInnerBefore.zCoordinate);
  BOOST_CHECK_EQUAL(clusterMiddle.xCoordinate, clusterMiddleBefore.xCoordinate);
  BOOST_CHECK_EQUAL(clusterMiddle.zCoordinate, clusterMiddleBefore.zCoordinate);
  BOOST_CHECK_EQUAL(clusterOuter.xCoordinate, clusterOuterBefore.xCoordinate);
  BOOST_CHECK_EQUAL(clusterOuter.zCoordinate, clusterOuterBefore.zCoordinate);
  BOOST_CHECK_EQUAL(hitInner.zCoordinate, hitInnerBefore.zCoordinate);
  BOOST_CHECK_EQUAL(hitMiddle.zCoordinate, hitMiddleBefore.zCoordinate);
  BOOST_CHECK_EQUAL(hitOuter.zCoordinate, hitOuterBefore.zCoordinate);
}

BOOST_AUTO_TEST_CASE(DiskBuildCellSeedUsesMaterialSlotsTwoOneZero)
{
  const auto clusterInner = makeGlobalCluster(1.0f, 0.5f, -0.4f, 10);
  const auto clusterMiddle = makeGlobalCluster(1.3f, 0.62f, -0.6f, 20);
  const auto clusterOuter = makeGlobalCluster(1.7f, 0.78f, -0.9f, 30);
  const auto hitInner = makeDiskHit(-0.4f, 1.0f, 0.5f);
  const auto hitMiddle = makeDiskHit(-0.6f, 1.3f, 0.62f);
  const auto hitOuter = makeDiskHit(-0.9f, 1.7f, 0.78f);
  // All three slots distinct: an index-order bug (e.g. inner/outer swapped)
  // would desynchronize the policy result from the correctly-ordered
  // reference below.
  const std::array<float, 3> xOverX0{0.005f, 0.05f, 0.1f}; // inner, middle, outer
  const auto material = toMaterial(xOverX0);

  DiskDiskPolicyParams params;
  params.trackletMinPt = 0.3f;
  params.maxChi2ClusterAttachment = 1.e6f;

  o2::track::TrackParCovFwd policyState{};
  float policyChi2 = 0.f;
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::DiskDisk>(clusterInner, clusterMiddle, clusterOuter,
                                                             hitInner, hitMiddle, hitOuter, material, Bz,
                                                             policyState, policyChi2, params));

  o2::track::TrackParCovFwd referenceState{};
  float referenceChi2 = 0.f;
  BOOST_REQUIRE(legacyDiskCellFit(clusterInner, clusterMiddle, clusterOuter, hitOuter, hitMiddle, hitInner,
                                  xOverX0[2], xOverX0[1], xOverX0[0], Bz, params.trackletMinPt,
                                  params.maxChi2ClusterAttachment, referenceState, referenceChi2));
  checkDiskStateEqual(policyState, referenceState);
  BOOST_CHECK_EQUAL(policyChi2, referenceChi2);
}

BOOST_AUTO_TEST_CASE(DiskBuildCellSeedIgnoresArealDensityMCSOnly)
{
  // DiskDisk passes only material.xOverX0 to the unchanged legacy MCS-only
  // MFT path (detail::mftFwdAttachCluster); arealDensityGPerCm2 must never
  // reach it. Two material sets sharing xOverX0 but wildly differing in
  // arealDensityGPerCm2 must therefore produce byte-identical state and chi2.
  const auto clusterInner = makeGlobalCluster(1.0f, 0.5f, -0.4f, 10);
  const auto clusterMiddle = makeGlobalCluster(1.3f, 0.62f, -0.6f, 20);
  const auto clusterOuter = makeGlobalCluster(1.7f, 0.78f, -0.9f, 30);
  const auto hitInner = makeDiskHit(-0.4f, 1.0f, 0.5f);
  const auto hitMiddle = makeDiskHit(-0.6f, 1.3f, 0.62f);
  const auto hitOuter = makeDiskHit(-0.9f, 1.7f, 0.78f);

  DiskDiskPolicyParams params;
  params.trackletMinPt = 0.3f;
  params.maxChi2ClusterAttachment = 1.e6f;

  const std::array<NominalSurfaceMaterial, 3> lowDensity{
    NominalSurfaceMaterial{0.005f, 0.001f}, NominalSurfaceMaterial{0.05f, 0.002f}, NominalSurfaceMaterial{0.1f, 0.003f}};
  const std::array<NominalSurfaceMaterial, 3> highDensity{
    NominalSurfaceMaterial{0.005f, 500.f}, NominalSurfaceMaterial{0.05f, 800.f}, NominalSurfaceMaterial{0.1f, 999.f}};

  o2::track::TrackParCovFwd lowDensityState{};
  float lowDensityChi2 = 0.f;
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::DiskDisk>(clusterInner, clusterMiddle, clusterOuter,
                                                             hitInner, hitMiddle, hitOuter, lowDensity, Bz,
                                                             lowDensityState, lowDensityChi2, params));

  o2::track::TrackParCovFwd highDensityState{};
  float highDensityChi2 = 0.f;
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::DiskDisk>(clusterInner, clusterMiddle, clusterOuter,
                                                             hitInner, hitMiddle, hitOuter, highDensity, Bz,
                                                             highDensityState, highDensityChi2, params));

  checkDiskStateEqual(lowDensityState, highDensityState);
  BOOST_CHECK_EQUAL(lowDensityChi2, highDensityChi2);
}

BOOST_AUTO_TEST_CASE(DiskBuildCellSeedRepeatedCallsAreDeterministic)
{
  const auto clusterInner = makeGlobalCluster(1.0f, 0.5f, -0.4f, 10);
  const auto clusterMiddle = makeGlobalCluster(1.3f, 0.62f, -0.6f, 20);
  const auto clusterOuter = makeGlobalCluster(1.7f, 0.78f, -0.9f, 30);
  const auto hitInner = makeDiskHit(-0.4f, 1.0f, 0.5f);
  const auto hitMiddle = makeDiskHit(-0.6f, 1.3f, 0.62f);
  const auto hitOuter = makeDiskHit(-0.9f, 1.7f, 0.78f);
  const std::array<float, 3> xOverX0{0.015f, 0.017f, 0.02f};
  const auto material = toMaterial(xOverX0);

  DiskDiskPolicyParams params;
  params.trackletMinPt = 0.3f;
  params.maxChi2ClusterAttachment = 1.e6f;

  o2::track::TrackParCovFwd firstState{};
  float firstChi2 = 0.f;
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::DiskDisk>(clusterInner, clusterMiddle, clusterOuter,
                                                             hitInner, hitMiddle, hitOuter, material, Bz,
                                                             firstState, firstChi2, params));

  o2::track::TrackParCovFwd secondState{};
  float secondChi2 = 0.f;
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::DiskDisk>(clusterInner, clusterMiddle, clusterOuter,
                                                             hitInner, hitMiddle, hitOuter, material, Bz,
                                                             secondState, secondChi2, params));

  checkDiskStateEqual(firstState, secondState);
  BOOST_CHECK_EQUAL(firstChi2, secondChi2);
}

/// Gate 3 cell-road pre-cut slice coverage: passesCellRoadPrecut<Tag>
/// (TransitionPolicyOperations.h), the last detector branch removed from
/// TrackerTraits::computeLayerCellsForPolicy's candidate loop. These tests
/// use an independent re-derivation of the legacy formula (formerly
/// detail::validateMFTCellClusters/mftDistanceToSeedSquared/
/// mftConicalRoadR2Scale, all three now deleted from MFTFwdTrackHelpers.h)
/// as the oracle, never the operation itself, so a transcription mistake in
/// either would show up as a mismatch.

namespace
{
/// Independent re-derivation of the legacy squared transverse distance from
/// `point` to the seed line `from`->`to` (formerly
/// detail::mftDistanceToSeedSquared). |dzSeed| < 1e-9f returns FLT_MAX.
float referenceDistanceToSeedLineSquared(const GlobalPoint3F& from, const GlobalPoint3F& to, const GlobalPoint3F& point)
{
  const float dz = to.z - from.z;
  if (std::abs(dz) < 1e-9f) {
    return std::numeric_limits<float>::max();
  }
  const float t = (point.z - from.z) / dz;
  const float xOnLine = from.x + t * (to.x - from.x);
  const float yOnLine = from.y + t * (to.y - from.y);
  const float dx = point.x - xOnLine;
  const float dy = point.y - yOnLine;
  return dx * dx + dy * dy;
}

/// Independent re-derivation of the legacy conical road scale (formerly
/// detail::mftConicalRoadR2Scale), taking zFrom/zTo explicitly. |zFrom| <
/// 1e-6f falls back to 1.f. Preserved as `1.f + (zTo - zFrom) / zFrom`,
/// matching the operation -- not simplified to zTo/zFrom.
float referenceConicalRoadR2Scale(float zFrom, float zTo)
{
  if (std::abs(zFrom) < 1e-6f) {
    return 1.f;
  }
  const float scale = 1.f + (zTo - zFrom) / zFrom;
  return scale * scale;
}

/// Independent re-derivation of the combined three-check road pre-cut
/// (formerly detail::validateMFTCellClusters), used only as the oracle for
/// passesCellRoadPrecut<DiskDisk> -- never called by the operation under test.
bool referenceCellRoadPrecut(const GlobalPoint3F& inner, const GlobalPoint3F& middle, const GlobalPoint3F& outer,
                             float zInner, float zMiddle, float zOuter, float cellRoadRCut)
{
  const float r2Cut = cellRoadRCut * cellRoadRCut;
  const bool checkMiddle = referenceDistanceToSeedLineSquared(inner, outer, middle) < r2Cut * referenceConicalRoadR2Scale(zInner, zMiddle);
  const bool checkOuter = referenceDistanceToSeedLineSquared(inner, middle, outer) < r2Cut * referenceConicalRoadR2Scale(zInner, zOuter);
  const bool checkInner = referenceDistanceToSeedLineSquared(middle, outer, inner) < r2Cut * referenceConicalRoadR2Scale(zMiddle, zInner);
  return checkMiddle && checkOuter && checkInner;
}
} // namespace

BOOST_AUTO_TEST_CASE(CylinderCellRoadPrecutAlwaysAcceptsAndIgnoresEmptyReferenceSpan)
{
  // Garbage-ish points/layers and an empty reference span (exactly the state
  // TrackerTraits::mDiskLayerReferenceZ is left in for a CylinderCylinder
  // iteration): the CylinderCylinder specialization must not read any of
  // them and must still return true.
  const GlobalPoint3F garbage{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::infinity(), -1.f};
  BOOST_CHECK(passesCellRoadPrecut<TransitionPolicyTag::CylinderCylinder>(
    garbage, garbage, garbage, 0, 1, 2, gsl::span<const float>{}, CylinderCylinderPolicyParams{}));
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutMatchesIndependentOracleAcceptAndReject)
{
  const GlobalPoint3F pointInner{1.0f, 0.5f, -0.4f};
  const GlobalPoint3F pointMiddle{1.3f, 0.62f, -0.6f};
  const GlobalPoint3F pointOuter{1.7f, 0.78f, -0.9f};
  const std::array<float, 3> referenceZ{-45.3f, -46.7f, -48.6f};

  DiskDiskPolicyParams generous;
  generous.cellRoadRCut = 1000.f;
  const bool oracleAccepts = referenceCellRoadPrecut(pointInner, pointMiddle, pointOuter,
                                                     referenceZ[0], referenceZ[1], referenceZ[2], generous.cellRoadRCut);
  BOOST_REQUIRE(oracleAccepts);
  BOOST_CHECK_EQUAL(oracleAccepts,
                    passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
                      pointInner, pointMiddle, pointOuter, 0, 1, 2,
                      gsl::span<const float>(referenceZ), generous));

  DiskDiskPolicyParams tight;
  tight.cellRoadRCut = 1.e-6f;
  const bool oracleRejects = referenceCellRoadPrecut(pointInner, pointMiddle, pointOuter,
                                                     referenceZ[0], referenceZ[1], referenceZ[2], tight.cellRoadRCut);
  BOOST_REQUIRE(!oracleRejects);
  BOOST_CHECK_EQUAL(oracleRejects,
                    passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
                      pointInner, pointMiddle, pointOuter, 0, 1, 2,
                      gsl::span<const float>(referenceZ), tight));
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutMiddleCheckRejectsIndependently)
{
  // Non-collinear geometry (middle bulged in y) gives all three checks a
  // genuine positive distance. zMiddle tiny relative to zInner/zOuter makes
  // scale(zInner, zMiddle) ~ (zMiddle/zInner)^2 tiny -- checkMiddle's
  // threshold collapses to ~0 -- while scale(zInner, zOuter) and
  // scale(zMiddle, zInner) are correspondingly huge, so checkOuter/checkInner
  // pass regardless of their own (unexamined) geometric distance.
  const GlobalPoint3F pointInner{1.0f, 0.5f, -0.4f};
  const GlobalPoint3F pointMiddle{1.3f, 0.9f, -0.6f};
  const GlobalPoint3F pointOuter{1.7f, 0.78f, -0.9f};
  const float zInner = -40.f;
  const float zMiddle = 1.e-4f;
  const float zOuter = -42.f;

  const bool checkMiddle = referenceDistanceToSeedLineSquared(pointInner, pointOuter, pointMiddle) <
                           1.f * referenceConicalRoadR2Scale(zInner, zMiddle);
  const bool checkOuter = referenceDistanceToSeedLineSquared(pointInner, pointMiddle, pointOuter) <
                          1.f * referenceConicalRoadR2Scale(zInner, zOuter);
  const bool checkInner = referenceDistanceToSeedLineSquared(pointMiddle, pointOuter, pointInner) <
                          1.f * referenceConicalRoadR2Scale(zMiddle, zInner);
  BOOST_REQUIRE(!checkMiddle);
  BOOST_REQUIRE(checkOuter);
  BOOST_REQUIRE(checkInner);

  const std::array<float, 3> referenceZ{zInner, zMiddle, zOuter};
  DiskDiskPolicyParams params;
  params.cellRoadRCut = 1.f;
  BOOST_CHECK(!passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
    pointInner, pointMiddle, pointOuter, 0, 1, 2,
    gsl::span<const float>(referenceZ), params));
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutOuterCheckRejectsIndependently)
{
  const GlobalPoint3F pointInner{1.0f, 0.5f, -0.4f};
  const GlobalPoint3F pointMiddle{1.3f, 0.9f, -0.6f};
  const GlobalPoint3F pointOuter{1.7f, 0.78f, -0.9f};
  const float zInner = -40.f;
  const float zMiddle = -40.f; // scale(zInner, zMiddle) == scale(zMiddle, zInner) == 1 exactly
  const float zOuter = 1.e-4f; // scale(zInner, zOuter) ~ (zOuter/zInner)^2, tiny -> checkOuter collapses

  const bool checkMiddle = referenceDistanceToSeedLineSquared(pointInner, pointOuter, pointMiddle) <
                           100.f * referenceConicalRoadR2Scale(zInner, zMiddle);
  const bool checkOuter = referenceDistanceToSeedLineSquared(pointInner, pointMiddle, pointOuter) <
                          100.f * referenceConicalRoadR2Scale(zInner, zOuter);
  const bool checkInner = referenceDistanceToSeedLineSquared(pointMiddle, pointOuter, pointInner) <
                          100.f * referenceConicalRoadR2Scale(zMiddle, zInner);
  BOOST_REQUIRE(checkMiddle);
  BOOST_REQUIRE(!checkOuter);
  BOOST_REQUIRE(checkInner);

  const std::array<float, 3> referenceZ{zInner, zMiddle, zOuter};
  DiskDiskPolicyParams params;
  params.cellRoadRCut = 10.f; // r2Cut == 100.f, matching the checks above
  BOOST_CHECK(!passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
    pointInner, pointMiddle, pointOuter, 0, 1, 2,
    gsl::span<const float>(referenceZ), params));
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutInnerCheckRejectsIndependently)
{
  const GlobalPoint3F pointInner{1.0f, 0.5f, -0.4f};
  const GlobalPoint3F pointMiddle{1.3f, 0.9f, -0.6f};
  const GlobalPoint3F pointOuter{1.7f, 0.78f, -0.9f};
  const float zInner = 1.e-4f; // scale(zMiddle, zInner) ~ (zInner/zMiddle)^2, tiny -> checkInner collapses
  const float zMiddle = -40.f;
  const float zOuter = -42.f;

  const bool checkMiddle = referenceDistanceToSeedLineSquared(pointInner, pointOuter, pointMiddle) <
                           1.f * referenceConicalRoadR2Scale(zInner, zMiddle);
  const bool checkOuter = referenceDistanceToSeedLineSquared(pointInner, pointMiddle, pointOuter) <
                          1.f * referenceConicalRoadR2Scale(zInner, zOuter);
  const bool checkInner = referenceDistanceToSeedLineSquared(pointMiddle, pointOuter, pointInner) <
                          1.f * referenceConicalRoadR2Scale(zMiddle, zInner);
  BOOST_REQUIRE(checkMiddle);
  BOOST_REQUIRE(checkOuter);
  BOOST_REQUIRE(!checkInner);

  const std::array<float, 3> referenceZ{zInner, zMiddle, zOuter};
  DiskDiskPolicyParams params;
  params.cellRoadRCut = 1.f;
  BOOST_CHECK(!passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
    pointInner, pointMiddle, pointOuter, 0, 1, 2,
    gsl::span<const float>(referenceZ), params));
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutExactBoundaryAndNextRepresentableInsideOutside)
{
  // Degenerate xy-collinear inner/middle (same x,y, different z) makes the
  // checkOuter (line inner->middle) and checkInner (line middle->outer)
  // distances both exactly D*D by construction; checkMiddle (line
  // inner->outer) is comfortably smaller (0.25*D*D). Uniform reference z
  // gives scale == 1.f exactly for all three checks, so r2Cut == cut*cut
  // directly gates the outcome.
  constexpr float D = 0.25f;
  const GlobalPoint3F pointInner{0.f, 0.f, -10.f};
  const GlobalPoint3F pointMiddle{0.f, 0.f, -20.f};
  const GlobalPoint3F pointOuter{0.f, D, -30.f};
  const std::array<float, 3> referenceZ{-1.f, -1.f, -1.f};

  const float distMiddle = referenceDistanceToSeedLineSquared(pointInner, pointOuter, pointMiddle);
  const float distOuter = referenceDistanceToSeedLineSquared(pointInner, pointMiddle, pointOuter);
  const float distInner = referenceDistanceToSeedLineSquared(pointMiddle, pointOuter, pointInner);
  BOOST_REQUIRE_EQUAL(distOuter, D * D);
  BOOST_REQUIRE_EQUAL(distInner, D * D);
  BOOST_REQUIRE_LT(distMiddle, D * D);

  const float atBoundary = D;
  BOOST_REQUIRE_EQUAL(atBoundary * atBoundary, D * D); // same multiplication, bit-identical to distOuter/distInner

  DiskDiskPolicyParams params;
  params.cellRoadRCut = atBoundary;
  BOOST_CHECK(!passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
    pointInner, pointMiddle, pointOuter, 0, 1, 2, gsl::span<const float>(referenceZ), params));

  const float justAbove = std::nextafter(atBoundary, std::numeric_limits<float>::max());
  const float justBelow = std::nextafter(atBoundary, 0.f);
  BOOST_REQUIRE_GT(justAbove * justAbove, D * D);
  BOOST_REQUIRE_LT(justBelow * justBelow, D * D);

  DiskDiskPolicyParams paramsAbove;
  paramsAbove.cellRoadRCut = justAbove;
  BOOST_CHECK(passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
    pointInner, pointMiddle, pointOuter, 0, 1, 2, gsl::span<const float>(referenceZ), paramsAbove));

  DiskDiskPolicyParams paramsBelow;
  paramsBelow.cellRoadRCut = justBelow;
  BOOST_CHECK(!passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
    pointInner, pointMiddle, pointOuter, 0, 1, 2, gsl::span<const float>(referenceZ), paramsBelow));
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutNearZeroSeedDzRejects)
{
  // inner and outer share the same z: the checkMiddle line (inner->outer) has
  // |dzSeed| < 1e-9f, so referenceDistanceToSeedLineSquared/the operation's
  // internal equivalent both return FLT_MAX for that check, which cannot be
  // less than any finite threshold.
  const GlobalPoint3F pointInner{1.0f, 0.5f, -40.f};
  const GlobalPoint3F pointMiddle{1.3f, 0.62f, -41.f};
  const GlobalPoint3F pointOuter{1.7f, 0.78f, -40.f}; // same z as inner
  const std::array<float, 3> referenceZ{-45.f, -46.f, -47.f};

  BOOST_REQUIRE_EQUAL(referenceDistanceToSeedLineSquared(pointInner, pointOuter, pointMiddle),
                      std::numeric_limits<float>::max());

  DiskDiskPolicyParams params;
  params.cellRoadRCut = 1.e6f; // even a very generous cut cannot beat FLT_MAX
  BOOST_CHECK(!passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
    pointInner, pointMiddle, pointOuter, 0, 1, 2, gsl::span<const float>(referenceZ), params));
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutZeroReferenceZFallbackMatchesUnitScale)
{
  // zInner == 0.f triggers the |zFrom| < 1e-6f fallback (scale == 1.f) for
  // the two pairs where zInner is zFrom: checkMiddle's (zInner, zMiddle) and
  // checkOuter's (zInner, zOuter). The third pair, checkInner's (zMiddle,
  // zInner), uses zInner as zTo, not zFrom, so it is NOT protected by the
  // fallback: with zTo == 0 exactly, dCone == 1 + (0 - zMiddle) / zMiddle ==
  // 0 exactly for any nonzero zMiddle, giving scale == 0.f, not 1.f. This
  // documents that the fallback applies only where the formula defines it
  // (zFrom itself near zero), not to every pair that merely involves zInner,
  // and cross-checks the full predicate against the independent oracle for
  // this exact mixed-scale configuration.
  const GlobalPoint3F pointInner{1.0f, 0.5f, -0.4f};
  const GlobalPoint3F pointMiddle{1.3f, 0.62f, -0.6f};
  const GlobalPoint3F pointOuter{1.7f, 0.78f, -0.9f};
  const float zInner = 0.f;
  const float zMiddle = -40.f;
  const float zOuter = -46.f;

  BOOST_CHECK_EQUAL(referenceConicalRoadR2Scale(zInner, zMiddle), 1.f); // checkMiddle: zFrom == zInner == 0 -> fallback
  BOOST_CHECK_EQUAL(referenceConicalRoadR2Scale(zInner, zOuter), 1.f);  // checkOuter: zFrom == zInner == 0 -> fallback
  BOOST_CHECK_EQUAL(referenceConicalRoadR2Scale(zMiddle, zInner), 0.f); // checkInner: zFrom == zMiddle != 0 -> real formula, not fallback

  const std::array<float, 3> referenceZ{zInner, zMiddle, zOuter};
  DiskDiskPolicyParams params;
  params.cellRoadRCut = 0.2f;

  const bool oracle = referenceCellRoadPrecut(pointInner, pointMiddle, pointOuter, zInner, zMiddle, zOuter, params.cellRoadRCut);
  BOOST_CHECK_EQUAL(oracle,
                    passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
                      pointInner, pointMiddle, pointOuter, 0, 1, 2,
                      gsl::span<const float>(referenceZ), params));
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutNonAdjacentLayerIndices)
{
  // perLayerReferenceZ has entries for 6 legacy layout-local layers with
  // distinct values; layerInner/layerMiddle/layerOuter (0, 2, 5) skip several
  // indices, exactly as a road spanning non-adjacent transitions would.
  const std::array<float, 6> referenceZ{-10.f, -999.f, -30.f, -999.f, -999.f, -50.f};
  const GlobalPoint3F pointInner{1.0f, 0.5f, -0.4f};
  const GlobalPoint3F pointMiddle{1.3f, 0.62f, -0.6f};
  const GlobalPoint3F pointOuter{1.7f, 0.78f, -0.9f};
  DiskDiskPolicyParams params;
  params.cellRoadRCut = 1000.f;

  const bool oracle = referenceCellRoadPrecut(pointInner, pointMiddle, pointOuter,
                                              referenceZ[0], referenceZ[2], referenceZ[5], params.cellRoadRCut);
  BOOST_CHECK_EQUAL(oracle,
                    passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
                      pointInner, pointMiddle, pointOuter, 0, 2, 5,
                      gsl::span<const float>(referenceZ), params));
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutNegativeDiskZRepresentative)
{
  // Representative realistic MFT-scale negative half-disk z values (the
  // production convention); no sign-based branch exists in the formula, so
  // this is exercised identically to every other test above -- kept as its
  // own case so the negative-z convention is explicitly documented/covered.
  const GlobalPoint3F pointInner{2.5f, -1.2f, -45.3f};
  const GlobalPoint3F pointMiddle{2.6f, -1.25f, -46.7f};
  const GlobalPoint3F pointOuter{2.75f, -1.32f, -48.6f};
  const std::array<float, 3> referenceZ{-45.3f, -46.7f, -48.6f};
  DiskDiskPolicyParams params;
  params.cellRoadRCut = 0.5f;

  const bool oracle = referenceCellRoadPrecut(pointInner, pointMiddle, pointOuter,
                                              referenceZ[0], referenceZ[1], referenceZ[2], params.cellRoadRCut);
  BOOST_CHECK_EQUAL(oracle,
                    passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
                      pointInner, pointMiddle, pointOuter, 0, 1, 2,
                      gsl::span<const float>(referenceZ), params));
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutInputsAreNotMutated)
{
  GlobalPoint3F pointInner{1.0f, 0.5f, -0.4f};
  GlobalPoint3F pointMiddle{1.3f, 0.62f, -0.6f};
  GlobalPoint3F pointOuter{1.7f, 0.78f, -0.9f};
  const auto pointInnerBefore = pointInner;
  const auto pointMiddleBefore = pointMiddle;
  const auto pointOuterBefore = pointOuter;
  std::array<float, 3> referenceZ{-45.3f, -46.7f, -48.6f};
  const auto referenceZBefore = referenceZ;
  DiskDiskPolicyParams params;
  params.cellRoadRCut = 1000.f;
  const auto paramsBefore = params;

  passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
    pointInner, pointMiddle, pointOuter, 0, 1, 2, gsl::span<const float>(referenceZ), params);

  BOOST_CHECK_EQUAL(pointInner.x, pointInnerBefore.x);
  BOOST_CHECK_EQUAL(pointInner.y, pointInnerBefore.y);
  BOOST_CHECK_EQUAL(pointInner.z, pointInnerBefore.z);
  BOOST_CHECK_EQUAL(pointMiddle.x, pointMiddleBefore.x);
  BOOST_CHECK_EQUAL(pointOuter.x, pointOuterBefore.x);
  for (size_t i = 0; i < referenceZ.size(); ++i) {
    BOOST_CHECK_EQUAL(referenceZ[i], referenceZBefore[i]);
  }
  BOOST_CHECK_EQUAL(params.cellRoadRCut, paramsBefore.cellRoadRCut);
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutNaNClusterCoordinateRejects)
{
  const GlobalPoint3F pointInner{1.0f, 0.5f, -0.4f};
  const GlobalPoint3F pointMiddle{std::numeric_limits<float>::quiet_NaN(), 0.62f, -0.6f};
  const GlobalPoint3F pointOuter{1.7f, 0.78f, -0.9f};
  const std::array<float, 3> referenceZ{-45.f, -46.f, -47.f};
  DiskDiskPolicyParams params;
  params.cellRoadRCut = 1000.f; // generous: only the NaN can cause rejection
  BOOST_CHECK(!passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
    pointInner, pointMiddle, pointOuter, 0, 1, 2,
    gsl::span<const float>(referenceZ), params));
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutNaNReferenceZRejects)
{
  const GlobalPoint3F pointInner{1.0f, 0.5f, -0.4f};
  const GlobalPoint3F pointMiddle{1.3f, 0.62f, -0.6f};
  const GlobalPoint3F pointOuter{1.7f, 0.78f, -0.9f};
  const std::array<float, 3> referenceZ{-45.f, std::numeric_limits<float>::quiet_NaN(), -47.f};
  DiskDiskPolicyParams params;
  params.cellRoadRCut = 1000.f;
  BOOST_CHECK(!passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
    pointInner, pointMiddle, pointOuter, 0, 1, 2,
    gsl::span<const float>(referenceZ), params));
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutNaNCutRejectsInIsolatedPredicate)
{
  // Exercises the isolated predicate directly with a non-finite cut --
  // DiskDiskPolicyParams::isValid() already prevents this from reaching
  // traversal in production (see the dedicated isValid() test below); this
  // test documents the predicate's own, otherwise-untested, arithmetic
  // behavior for this input.
  const GlobalPoint3F pointInner{1.0f, 0.5f, -0.4f};
  const GlobalPoint3F pointMiddle{1.3f, 0.62f, -0.6f};
  const GlobalPoint3F pointOuter{1.7f, 0.78f, -0.9f};
  const std::array<float, 3> referenceZ{-45.f, -46.f, -47.f};
  DiskDiskPolicyParams params;
  params.cellRoadRCut = std::numeric_limits<float>::quiet_NaN();
  BOOST_CHECK(!passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
    pointInner, pointMiddle, pointOuter, 0, 1, 2,
    gsl::span<const float>(referenceZ), params));
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutInfiniteCutAcceptsInIsolatedPredicate)
{
  // An infinite cut makes every finite-scale threshold infinite, so every
  // finite distance passes -- this is a real, different-from-NaN behavior of
  // the isolated formula (not "all non-finite input rejects"). Reference z
  // values stay finite and comparable in magnitude so every scale() call
  // stays finite. DiskDiskPolicyParams::isValid() already prevents this input
  // from reaching production traversal (see the dedicated isValid() test).
  const GlobalPoint3F pointInner{1.0f, 0.5f, -0.4f};
  const GlobalPoint3F pointMiddle{1.3f, 0.62f, -0.6f};
  const GlobalPoint3F pointOuter{1.7f, 0.78f, -0.9f};
  const std::array<float, 3> referenceZ{-45.f, -46.f, -47.f};
  DiskDiskPolicyParams params;
  params.cellRoadRCut = std::numeric_limits<float>::infinity();
  BOOST_CHECK(passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
    pointInner, pointMiddle, pointOuter, 0, 1, 2,
    gsl::span<const float>(referenceZ), params));
}

BOOST_AUTO_TEST_CASE(DiskCellRoadPrecutInfiniteReferenceZRejects)
{
  // zInner == Inf serves as zFrom in both the checkMiddle and checkOuter
  // scale() calls: dCone = 1 + (zTo - Inf)/Inf is the Inf/Inf indeterminate
  // form -> NaN -> that scale is NaN -> those checks' `<` comparisons are
  // false regardless of geometry, so the combined AND rejects even though
  // zInner also serves as zTo in checkInner's pair, where the same Inf
  // produces an accepting (loose) scale -- the rejection is not universal
  // across all three checks, only sufficient via short-circuit AND.
  const GlobalPoint3F pointInner{1.0f, 0.5f, -0.4f};
  const GlobalPoint3F pointMiddle{1.3f, 0.62f, -0.6f};
  const GlobalPoint3F pointOuter{1.7f, 0.78f, -0.9f};
  const float zInner = std::numeric_limits<float>::infinity();
  const float zMiddle = -46.f;
  const float zOuter = -47.f;

  BOOST_CHECK(std::isnan(referenceConicalRoadR2Scale(zInner, zMiddle)));
  BOOST_CHECK(std::isnan(referenceConicalRoadR2Scale(zInner, zOuter)));

  const std::array<float, 3> referenceZ{zInner, zMiddle, zOuter};
  DiskDiskPolicyParams params;
  params.cellRoadRCut = 1000.f;
  BOOST_CHECK(!passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
    pointInner, pointMiddle, pointOuter, 0, 1, 2,
    gsl::span<const float>(referenceZ), params));
}

BOOST_AUTO_TEST_CASE(DiskDiskPolicyParamsIsValidRejectsNonFiniteCellRoadRCut)
{
  // Production guard (TransitionPolicyState.h): initialiseTimeFrame() calls
  // this before any candidate is evaluated, so the infinite/NaN-cut behavior
  // exercised in isolation above never reaches traversal in practice.
  DiskDiskPolicyParams infCut;
  infCut.cellRoadRCut = std::numeric_limits<float>::infinity();
  BOOST_CHECK(!infCut.isValid());

  DiskDiskPolicyParams nanCut;
  nanCut.cellRoadRCut = std::numeric_limits<float>::quiet_NaN();
  BOOST_CHECK(!nanCut.isValid());

  DiskDiskPolicyParams negativeCut;
  negativeCut.cellRoadRCut = -0.05f;
  BOOST_CHECK(!negativeCut.isValid());

  DiskDiskPolicyParams valid;
  valid.cellRoadRCut = 0.05f;
  BOOST_CHECK(valid.isValid());
}

/// Gate 3 transition-preparation slice coverage (relocated from
/// TimeFrame::initialise() into TrackerTraits::initialiseTimeFrame(); see
/// TransitionPolicyOperations.h layerMultipleScatteringAngle<Tag>,
/// clampTransitionCurvature<Tag>, prepareTransitionScatteringAndBending, and
/// TransitionPolicyBinding.h LayerGeometryConfigView). These tests verify
/// exact legacy-formula parity, the family-specific arithmetic literal that
/// integration review required preserved (not canonicalized), and the
/// order-sensitive oneOverR ratchet -- independently of TrackerTraits'
/// production traversal (covered separately in
/// testComputeLayerTrackletsOrchestration.cxx).

namespace
{
/// Independent re-transcription of the shared arithmetic in the frozen
/// ITS-only TimeFrame::initialise() (ITS/tracking/src/TimeFrame.cxx:352-370),
/// which the (now-removed) common-CA non-MFT branch reproduced verbatim.
/// Deliberately re-derived here rather than calling
/// prepareTransitionScatteringAndBending, so a transcription mistake in
/// either the operation or this reference would show up as a mismatch.
TransitionScatteringBendingPrep referenceTransitionScatteringAndBending(
  gsl::span<const float> perLayerMSAngle, int fromLayer, int toLayer,
  float r1, float r2, float clampedOneOverR, float res1, float res2)
{
  float ms2 = 0.f;
  for (int layer = fromLayer; layer < toLayer; ++layer) {
    ms2 += o2::its::math_utils::Sq(perLayerMSAngle[layer]);
  }
  const float msAngle = o2::gpu::CAMath::Sqrt(ms2);
  const float cosTheta1half = o2::gpu::CAMath::Sqrt(1.f - o2::its::math_utils::Sq(0.5f * r1 * clampedOneOverR));
  const float cosTheta2half = o2::gpu::CAMath::Sqrt(1.f - o2::its::math_utils::Sq(0.5f * r2 * clampedOneOverR));
  const float x = (r2 * cosTheta1half) - (r1 * cosTheta2half);
  const float delta = o2::gpu::CAMath::Sqrt(1.f / (1.f - 0.25f * o2::its::math_utils::Sq(x * clampedOneOverR)) *
                                            (o2::its::math_utils::Sq((0.25f * r1 * r2 * o2::its::math_utils::Sq(clampedOneOverR) / cosTheta2half) + cosTheta1half) * o2::its::math_utils::Sq(res1) +
                                             o2::its::math_utils::Sq((0.25f * r1 * r2 * o2::its::math_utils::Sq(clampedOneOverR) / cosTheta1half) + cosTheta2half) * o2::its::math_utils::Sq(res2)));
  const float phiCut = o2::gpu::CAMath::Min(o2::gpu::CAMath::ASin(0.5f * x * clampedOneOverR) + 2.f * msAngle + delta, o2::constants::math::PI * 0.5f);
  return TransitionScatteringBendingPrep{msAngle, phiCut};
}
} // namespace

BOOST_AUTO_TEST_CASE(CylinderScatteringAngleMatchesFrozenITSFormula)
{
  // Bit-exact vs the frozen ITS expression (ITS/tracking/src/TimeFrame.cxx:347):
  // math_utils::MSangle(0.14f, trkParam.TrackletMinPt, trkParam.LayerxX0[iLayer]).
  const std::array<float, 4> xX0Values{0.f, -0.001f, 5.e-3f, 1.e-2f};
  const std::array<float, 3> trackletMinPtValues{0.1f, 0.3f, 2.5f};
  for (float xX0 : xX0Values) {
    for (float trackletMinPt : trackletMinPtValues) {
      const float reference = o2::its::math_utils::MSangle(0.14f, trackletMinPt, xX0);
      const float actual = layerMultipleScatteringAngle<TransitionPolicyTag::CylinderCylinder>(
        LayerScatteringInputs<TransitionPolicyTag::CylinderCylinder>{xX0}, trackletMinPt);
      BOOST_CHECK_EQUAL(actual, reference);
    }
  }
  // xX0 <= 0 behavior, explicit: legacy MSangle maps this to zero, not a
  // rejection; the typed operation must not add validation beyond it.
  BOOST_CHECK_EQUAL(layerMultipleScatteringAngle<TransitionPolicyTag::CylinderCylinder>(
                      LayerScatteringInputs<TransitionPolicyTag::CylinderCylinder>{0.f}, 0.3f),
                    0.f);
}

BOOST_AUTO_TEST_CASE(DiskScatteringAngleMatchesLegacyMftFormulaWithExplicitReferenceZ)
{
  // Bit-exact vs the legacy detail::mftLayerMSAngle(layer, params), except
  // the DiskDisk operation receives referenceCoordinate/layerRadius
  // explicitly instead of calling mftLayerZ()/LayerZCoordinate() internally.
  // mftLayerZ() is used here only to construct the *expected* legacy value,
  // exactly as this operation's caller (TrackerTraits::initialiseTimeFrame(),
  // via bindLegacyMFTReferenceCoordinates()) is required to do.
  TrackingParameters legacy;
  resetDetectorDefaults(legacy, o2::detectors::DetID::MFT);
  for (int layer : {0, 3, o2::mft::constants::mft::LayersNumber - 1}) {
    const float referenceZ = detail::mftLayerZ(layer);
    const float radius = legacy.LayerRadii[layer];
    const float xX0 = legacy.LayerxX0[layer];

    const float reference = detail::mftLayerMSAngle(layer, legacy);
    const float actual = layerMultipleScatteringAngle<TransitionPolicyTag::DiskDisk>(
      LayerScatteringInputs<TransitionPolicyTag::DiskDisk>{xX0, radius, referenceZ}, legacy.TrackletMinPt);
    BOOST_CHECK_EQUAL(actual, reference);
  }

  // xX0 == 0 behavior, explicit: the legacy formula has no special case for
  // it (sqrt(0 * cscLambda) == 0), and this operation must not add one.
  const float referenceZ = detail::mftLayerZ(0);
  const float zeroX0Actual = layerMultipleScatteringAngle<TransitionPolicyTag::DiskDisk>(
    LayerScatteringInputs<TransitionPolicyTag::DiskDisk>{0.f, legacy.LayerRadii[0], referenceZ}, legacy.TrackletMinPt);
  BOOST_CHECK_EQUAL(zeroX0Actual, 0.f);
}

BOOST_AUTO_TEST_CASE(DiskScatteringAngleNearZeroReferenceRadiusFallback)
{
  // Legacy fallback: |rRef| <= 1e-6 => tanlRef = 0 (detail::mftLayerMSAngle),
  // rather than dividing by a near-zero radius.
  TrackingParameters legacy;
  resetDetectorDefaults(legacy, o2::detectors::DetID::MFT);
  legacy.LayerRadii[0] = 1.e-9f; // below the legacy 1e-6 fallback threshold
  const float referenceZ = detail::mftLayerZ(0);

  const float reference = detail::mftLayerMSAngle(0, legacy);
  const float actual = layerMultipleScatteringAngle<TransitionPolicyTag::DiskDisk>(
    LayerScatteringInputs<TransitionPolicyTag::DiskDisk>{legacy.LayerxX0[0], legacy.LayerRadii[0], referenceZ},
    legacy.TrackletMinPt);
  BOOST_CHECK_EQUAL(actual, reference);

  // Cross-check the fallback actually engages: tanlRef == 0 (rRef below the
  // 1e-6 threshold) makes absTanl == 0, which is *not* > 1e-6 either, so
  // cscLambda takes the near-parallel-incidence sentinel 1e6f, not 1 -- i.e.
  // this input is genuinely exercising the near-zero-radius branch, not
  // merely reproducing an unrelated formula.
  const float expectedWithSentinelCscLambda = 0.0136f * (1.f / legacy.TrackletMinPt) * std::sqrt(legacy.LayerxX0[0] * 1.e6f);
  BOOST_CHECK_EQUAL(reference, expectedWithSentinelCscLambda);
}

BOOST_AUTO_TEST_CASE(LayerGeometryConfigViewChecksSpanSizeOnlyNotNumericValues)
{
  TrackingParameters legacy;
  legacy.LayerRadii = {0.f, -1.f, 3.91924f, 19.6213f, 24.5597f, 34.388f, 39.3329f}; // degenerate/negative radii
  legacy.LayerxX0 = {5.e-3f, 5.e-3f, 5.e-3f, 1.e-2f, 1.e-2f, 1.e-2f, 1.e-2f};       // valid
  const auto legacyMaterial = toMaterial(legacy.LayerxX0);

  const auto attachHitConfig = bindAttachHitPolicyConfig(gsl::span<const NominalSurfaceMaterial>(legacyMaterial), legacy);
  BOOST_CHECK(attachHitConfig.isValid(7));

  const auto geometryConfig = bindLayerGeometryConfig(legacy, attachHitConfig);
  // Legacy TimeFrame::initialise() never rejected degenerate/zero/negative
  // radii; this slice must not silently start doing so.
  BOOST_CHECK(geometryConfig.isValid(7));
  BOOST_CHECK(!geometryConfig.isValid(8)); // span-size check still applies

  // Negative xOverX0 must be rejected -- by AttachHitPolicyConfigView, the
  // single established contract for that data. LayerGeometryConfigView
  // borrows the same (rejected) span rather than independently re-validating
  // it, so it must not be read as a numeric-validity signal on its own.
  auto corrupted = legacy;
  corrupted.LayerxX0[3] = -1.f;
  const auto corruptedMaterial = toMaterial(corrupted.LayerxX0);
  const auto corruptedAttachHitConfig = bindAttachHitPolicyConfig(gsl::span<const NominalSurfaceMaterial>(corruptedMaterial), corrupted);
  BOOST_CHECK(!corruptedAttachHitConfig.isValid(7));
  const auto corruptedGeometryConfig = bindLayerGeometryConfig(corrupted, corruptedAttachHitConfig);
  BOOST_CHECK(corruptedGeometryConfig.isValid(7)); // size-only: still reports valid
}

BOOST_AUTO_TEST_CASE(BindLayerGeometryConfigBorrowsAttachHitLayerMaterialSpan)
{
  TrackingParameters legacy;
  const auto legacyMaterial = toMaterial(legacy.LayerxX0);
  const auto attachHitConfig = bindAttachHitPolicyConfig(gsl::span<const NominalSurfaceMaterial>(legacyMaterial), legacy);
  const auto geometryConfig = bindLayerGeometryConfig(legacy, attachHitConfig);
  BOOST_CHECK_EQUAL(geometryConfig.layerMaterial.data(), attachHitConfig.layerMaterial.data());
  BOOST_CHECK_EQUAL(geometryConfig.layerMaterial.size(), attachHitConfig.layerMaterial.size());
  BOOST_CHECK_EQUAL(geometryConfig.layerRadii.data(), legacy.LayerRadii.data());
}

BOOST_AUTO_TEST_CASE(ClampTransitionCurvatureMatchesExactLegacyExpressionPerFamily)
{
  // Integration review finding: the two legacy branches compare `0.5 *
  // oneOverR` (CylinderCylinder: double-promoted) against `0.5f * oneOverR`
  // (DiskDisk: float) before the same `>= 1.f / r2` clamp. Preserved
  // verbatim per family, not canonicalized -- see
  // ClampTransitionCurvatureFloatVersusDoubleDiscriminatorAttempt for why no
  // observable difference could be constructed, and note that preservation,
  // not the (negative) discriminator search, is what this test asserts.
  const std::array<std::pair<float, float>, 5> samples{{
    {0.001f, 50.f}, // clamp does not trigger
    {3.0f, 1.0f},   // clamp triggers
    {0.02f, 25.f},
    {0.5f, 0.9f},
    {0.0001f, 4.f},
  }};
  for (const auto& sample : samples) {
    const float oneOverR = sample.first;
    const float r2 = sample.second;

    const float cylinderActual = clampTransitionCurvature<TransitionPolicyTag::CylinderCylinder>(oneOverR, r2);
    const float cylinderReference = (0.5 * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
    BOOST_CHECK_EQUAL(cylinderActual, cylinderReference);

    const float diskActual = clampTransitionCurvature<TransitionPolicyTag::DiskDisk>(oneOverR, r2);
    const float diskReference = (0.5f * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
    BOOST_CHECK_EQUAL(diskActual, diskReference);
  }
}

BOOST_AUTO_TEST_CASE(ClampTransitionCurvatureFloatVersusDoubleDiscriminatorAttempt)
{
  // Attempted discriminator search. Multiplying/dividing an IEEE-754 value by
  // exactly 0.5 (a power of two) is *exact* in both binary32 and binary64 --
  // it only decrements the exponent, no mantissa rounding is required --
  // provided the result does not underflow into the subnormal range.
  // Consequently `0.5 * oneOverR` (double: `0.5` is an exact double literal
  // and `oneOverR` promotes to double losslessly) and `0.5f * oneOverR`
  // (float) are provably bit-identical once compared against the same
  // `1.f / r2`, for any `oneOverR` that is not itself subnormal. Physically,
  // `oneOverR == 0.001f*0.3f*|Bz|/TrackletMinPt` is many orders of magnitude
  // away from the float underflow boundary (~1e-4..1e-2 for realistic
  // Bz/TrackletMinPt), so no discriminating input exists in the physically
  // meaningful domain this code operates in. This sweeps representative and
  // extreme-but-normal values, including the clamp boundary itself, and
  // finds none; it does not probe genuinely subnormal `oneOverR` (below
  // ~1.18e-38), since those are far outside any value this code can produce
  // and would not demonstrate a real risk.
  const std::array<float, 9> oneOverRSamples{
    0.f, 1.e-6f, 1.e-4f, 1.e-3f, 1.e-2f, 0.1f, 1.f, 10.f, 1.e10f};
  const std::array<float, 7> r2Samples{
    0.5f, 1.f, 2.33959f, 5.f, 19.6213f, 39.3329f, 1.e6f};

  bool discriminatorFound = false;
  for (float oneOverR : oneOverRSamples) {
    for (float r2 : r2Samples) {
      const float viaFloatLiteral = (0.5f * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
      const float viaDoubleLiteral = (0.5 * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
      if (viaFloatLiteral != viaDoubleLiteral) {
        discriminatorFound = true;
      }
    }
  }
  BOOST_CHECK(!discriminatorFound); // documents the (negative) search result described above
}

BOOST_AUTO_TEST_CASE(CurvatureRatchetThreadsInIncreasingLegacyTransitionIdOrder)
{
  // The legacy oneOverR is a loop-carried variable, not reset per transition:
  // clampTransitionCurvature<Tag> must be called once per transition, in
  // increasing legacy transitionId order, threading its return value into the
  // next call. This proves the computation is genuinely order-sensitive (an
  // unproven iteration order, e.g. a policy-grouping span, cannot be
  // substituted without first proving it matches legacy transitionId order)
  // and pins the exact sequence increasing order must produce.
  constexpr float initialOneOverR = 3.f;
  const std::array<float, 3> r2InIncreasingTransitionIdOrder{1.f, 4.f, 0.5f};

  float running = initialOneOverR;
  std::array<float, 3> forwardResults{};
  for (int i = 0; i < 3; ++i) {
    running = clampTransitionCurvature<TransitionPolicyTag::CylinderCylinder>(running, r2InIncreasingTransitionIdOrder[i]);
    forwardResults[i] = running;
  }

  float manual = initialOneOverR;
  manual = (0.5 * manual >= 1.f / r2InIncreasingTransitionIdOrder[0]) ? (2.f / r2InIncreasingTransitionIdOrder[0]) - o2::constants::math::Almost0 : manual;
  BOOST_CHECK_EQUAL(forwardResults[0], manual);
  manual = (0.5 * manual >= 1.f / r2InIncreasingTransitionIdOrder[1]) ? (2.f / r2InIncreasingTransitionIdOrder[1]) - o2::constants::math::Almost0 : manual;
  BOOST_CHECK_EQUAL(forwardResults[1], manual);
  manual = (0.5 * manual >= 1.f / r2InIncreasingTransitionIdOrder[2]) ? (2.f / r2InIncreasingTransitionIdOrder[2]) - o2::constants::math::Almost0 : manual;
  BOOST_CHECK_EQUAL(forwardResults[2], manual);

  // A different processing order must not be assumed to reproduce the same
  // first step: proves the caller cannot substitute an unproven order.
  const float reversedFirstStep = clampTransitionCurvature<TransitionPolicyTag::CylinderCylinder>(
    initialOneOverR, r2InIncreasingTransitionIdOrder[2]);
  BOOST_CHECK_NE(reversedFirstStep, forwardResults[0]);
}

BOOST_AUTO_TEST_CASE(PrepareTransitionScatteringAndBendingMatchesFrozenFormulaForITSAndMFTShapedInputs)
{
  // ITS-shaped (cm-scale barrel radii from TrackingParameters defaults).
  {
    const std::array<float, 7> msAngles{1.e-3f, 1.1e-3f, 1.2e-3f, 2.e-3f, 2.1e-3f, 2.2e-3f, 2.3e-3f};
    constexpr int fromLayer = 0;
    constexpr int toLayer = 3; // half-open: sums layers 0,1,2 only
    constexpr float r1 = 2.33959f;
    constexpr float r2 = 19.6213f;
    constexpr float res1 = 5.e-4f;
    constexpr float res2 = 5.e-4f;
    const float oneOverR = clampTransitionCurvature<TransitionPolicyTag::CylinderCylinder>(
      0.001f * 0.3f * std::abs(Bz) / 0.3f, r2);
    const gsl::span<const float> msSpan(msAngles.data(), msAngles.size());
    const auto actual = prepareTransitionScatteringAndBending(msSpan, fromLayer, toLayer, r1, r2, oneOverR, res1, res2);
    const auto reference = referenceTransitionScatteringAndBending(msSpan, fromLayer, toLayer, r1, r2, oneOverR, res1, res2);
    BOOST_CHECK_EQUAL(actual.msAngle, reference.msAngle);
    BOOST_CHECK_EQUAL(actual.phiCut, reference.phiCut);

    // Half-open range: layer index `toLayer` itself must not contribute.
    const auto includingToLayer = referenceTransitionScatteringAndBending(msSpan, fromLayer, toLayer + 1, r1, r2, oneOverR, res1, res2);
    BOOST_CHECK_NE(actual.msAngle, includingToLayer.msAngle);
  }

  // MFT-shaped, deliberately skipped/non-adjacent transition (fromLayer=1,
  // toLayer=4: sums layers 1,2,3, skipping layer 4 itself as the endpoint).
  {
    TrackingParameters mft;
    resetDetectorDefaults(mft, o2::detectors::DetID::MFT);
    std::array<float, o2::mft::constants::mft::LayersNumber> msAngles{};
    for (int layer = 0; layer < o2::mft::constants::mft::LayersNumber; ++layer) {
      msAngles[layer] = layerMultipleScatteringAngle<TransitionPolicyTag::DiskDisk>(
        LayerScatteringInputs<TransitionPolicyTag::DiskDisk>{mft.LayerxX0[layer], mft.LayerRadii[layer], detail::mftLayerZ(layer)},
        mft.TrackletMinPt);
    }
    constexpr int fromLayer = 1;
    constexpr int toLayer = 4;
    const float r1 = mft.LayerRadii[fromLayer];
    const float r2 = mft.LayerRadii[toLayer];
    constexpr float res1 = 5.e-4f;
    constexpr float res2 = 6.e-4f;
    const float oneOverR = clampTransitionCurvature<TransitionPolicyTag::DiskDisk>(
      0.001f * 0.3f * std::abs(Bz) / mft.TrackletMinPt, r2);
    const gsl::span<const float> msSpan(msAngles.data(), msAngles.size());
    const auto actual = prepareTransitionScatteringAndBending(msSpan, fromLayer, toLayer, r1, r2, oneOverR, res1, res2);
    const auto reference = referenceTransitionScatteringAndBending(msSpan, fromLayer, toLayer, r1, r2, oneOverR, res1, res2);
    BOOST_CHECK_EQUAL(actual.msAngle, reference.msAngle);
    BOOST_CHECK_EQUAL(actual.phiCut, reference.phiCut);
  }
}

BOOST_AUTO_TEST_CASE(PrepareTransitionScatteringAndBendingZeroFieldAndDegenerateRadiusMatchLegacyFormula)
{
  const std::array<float, 3> msAngles{1.e-3f, 1.2e-3f, 1.4e-3f};
  const gsl::span<const float> msSpan(msAngles.data(), msAngles.size());

  // Zero field: oneOverR's initial value (before any clamp) is exactly 0,
  // matching the legacy `0.001f * 0.3f * std::abs(mBz) / trkParam.TrackletMinPt`.
  {
    const float zeroFieldOneOverR = 0.001f * 0.3f * std::abs(0.f) / 0.3f;
    BOOST_CHECK_EQUAL(zeroFieldOneOverR, 0.f);
    const float clamped = clampTransitionCurvature<TransitionPolicyTag::CylinderCylinder>(zeroFieldOneOverR, 5.f);
    BOOST_CHECK_EQUAL(clamped, 0.f); // 0.5*0 >= 1/5 is false: clamp does not trigger
    const auto actual = prepareTransitionScatteringAndBending(msSpan, 0, 2, 2.f, 5.f, clamped, 5.e-4f, 5.e-4f);
    const auto reference = referenceTransitionScatteringAndBending(msSpan, 0, 2, 2.f, 5.f, clamped, 5.e-4f, 5.e-4f);
    BOOST_CHECK_EQUAL(actual.msAngle, reference.msAngle);
    BOOST_CHECK_EQUAL(actual.phiCut, reference.phiCut);
  }

  // Degenerate radius (r2 == 0): legacy does not reject this -- it flows
  // through to whatever the floating-point expression produces. This test
  // asserts parity with that expression, not any particular finiteness.
  {
    const float oneOverR = clampTransitionCurvature<TransitionPolicyTag::CylinderCylinder>(0.01f, 0.f);
    const auto actual = prepareTransitionScatteringAndBending(msSpan, 0, 2, 2.f, 0.f, oneOverR, 5.e-4f, 5.e-4f);
    const auto reference = referenceTransitionScatteringAndBending(msSpan, 0, 2, 2.f, 0.f, oneOverR, 5.e-4f, 5.e-4f);
    // BOOST_CHECK_EQUAL on NaN is always false (NaN != NaN); compare the bit
    // pattern so a NaN-vs-NaN legacy-parity match is still recognized as a pass.
    BOOST_CHECK(std::memcmp(&actual.msAngle, &reference.msAngle, sizeof(float)) == 0);
    BOOST_CHECK(std::memcmp(&actual.phiCut, &reference.phiCut, sizeof(float)) == 0);
  }
}

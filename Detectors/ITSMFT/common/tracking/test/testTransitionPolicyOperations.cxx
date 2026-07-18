// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT TransitionPolicyOperations
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <limits>
#include <cmath>

#include <boost/test/unit_test.hpp>

#include <TGeoGlobalMagField.h>

#include "Field/MagneticField.h"
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/TransitionPolicyBinding.h"
#include "ITSMFTTracking/TransitionPolicyOperations.h"

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

  const auto attach = bindAttachHitPolicyConfig(legacy);
  BOOST_REQUIRE_EQUAL(attach.layerxX0.size(), 3u);
  BOOST_CHECK_CLOSE(attach.layerxX0[0], 0.011f, 1e-6);
  BOOST_CHECK_CLOSE(attach.layerxX0[1], 0.022f, 1e-6);
  BOOST_CHECK_CLOSE(attach.layerxX0[2], 0.033f, 1e-6);
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
  BOOST_CHECK(!bindAttachHitPolicyConfig(materialNaN).isValid(materialNaN.LayerxX0.size()));
  auto materialInf = legacy;
  materialInf.LayerxX0[4] = inf;
  BOOST_CHECK(!bindAttachHitPolicyConfig(materialInf).isValid(materialInf.LayerxX0.size()));
  auto invalidCorrection = legacy;
  invalidCorrection.CorrType = static_cast<o2::base::PropagatorF::MatCorrType>(99);
  BOOST_CHECK(!bindAttachHitPolicyConfig(invalidCorrection).isValid(invalidCorrection.LayerxX0.size()));
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
  BOOST_CHECK(attachHit<TransitionPolicyTag::CylinderCylinder>(policyState, hit, 0.f, NoMaterialCorrection, Bz, policyChi2, exact));
  BOOST_CHECK(legacyBarrelAttach(inlineState, hit, 0.f, NoMaterialCorrection, Bz, exact.maxChi2ClusterAttachment, inlineChi2));
  checkBarrelStateEqual(policyState, inlineState);
  BOOST_CHECK_EQUAL(policyChi2, inlineChi2);

  CylinderCylinderPolicyParams below = exact;
  below.maxChi2ClusterAttachment = std::nextafter(exactChi2, 0.f);
  policyState = initial;
  inlineState = initial;
  policyChi2 = inlineChi2 = 7.25f;
  BOOST_CHECK(!attachHit<TransitionPolicyTag::CylinderCylinder>(policyState, hit, 0.f, NoMaterialCorrection, Bz, policyChi2, below));
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
  BOOST_CHECK(!attachHit<TransitionPolicyTag::CylinderCylinder>(rotationState, rotationHit, 0.f, NoMaterialCorrection, Bz, chi2, params));
  checkBarrelStateEqual(rotationState, rotationBefore);
  BOOST_CHECK_EQUAL(chi2, 3.f);

  auto propagationState = makeBarrelState(0.f, 0.3f, 0.1f, 1.f, 0.f, 0.4f, 2000.f);
  const auto propagationHit = makeBarrelHit(5.f, 0.3f, 0.1f, 1.f);
  chi2 = 4.f;
  BOOST_CHECK(!attachHit<TransitionPolicyTag::CylinderCylinder>(propagationState, propagationHit, 0.f, NoMaterialCorrection, Bz, chi2, params));
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
  BOOST_CHECK(!attachHit<TransitionPolicyTag::CylinderCylinder>(updateState, singularHit, 0.f, NoMaterialCorrection, Bz, chi2, params));
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
  BOOST_CHECK(attachHit<TransitionPolicyTag::DiskDisk>(policyState, hit, 0.017f, NoMaterialCorrection, Bz, policyChi2, exact));
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
  BOOST_CHECK(!attachHit<TransitionPolicyTag::DiskDisk>(policyState, hit, 0.017f, NoMaterialCorrection, Bz, policyChi2, below));
  BOOST_CHECK(!legacyDiskAttach(inlineState, hit, 0.017f, Bz, below.maxChi2ClusterAttachment, inlineChi2));
  checkDiskStateEqual(policyState, before);
  checkDiskStateEqual(policyState, inlineState);
  BOOST_CHECK_EQUAL(policyChi2, 8.5f);
  BOOST_CHECK_EQUAL(policyChi2, inlineChi2);
}

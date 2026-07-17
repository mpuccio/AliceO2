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

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/TransitionPolicyBinding.h"
#include "ITSMFTTracking/TransitionPolicyOperations.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

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
}

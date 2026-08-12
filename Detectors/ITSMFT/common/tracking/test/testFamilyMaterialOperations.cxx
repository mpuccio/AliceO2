// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFTFamilyMaterialOperations
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

#include "ITSMFTTracking/BarrelSurfaceStateOperations.h"
#include "ITSMFTTracking/ForwardSurfaceStateOperations.h"
#include "ITSMFTTracking/detail/SurfaceKinematicStateLegacyAdapters.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

namespace
{
using namespace o2::itsmft::tracking;

constexpr float AbsTol = 5.e-5f;
constexpr float RelTol = 1.e-3f;

bool closeTo(float a, float b, float absTol = AbsTol, float relTol = RelTol)
{
  const float diff = std::fabs(a - b);
  return diff <= absTol || diff <= relTol * std::fabs(b);
}

template <typename T>
bool bitEqual(const T& lhs, const T& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
}

SurfaceKinematicState makeBarrelState(uint8_t absCharge = 1, float q2pt = 0.8f, o2::track::PID pid = o2::track::PID::Pion)
{
  SurfaceKinematicState state{};
  state.parameters[0] = 1.25f;  // Y
  state.parameters[1] = -0.75f; // Z
  state.parameters[2] = 0.2f;   // Snp
  state.parameters[3] = 0.35f;  // Tgl
  state.parameters[4] = q2pt;   // Q2Pt
  state.referenceCoordinate = 4.f;
  state.alpha = 0.3f;
  state.kind = SurfaceKind::Cylinder;
  state.flags = 0x5a;
  state.absCharge = absCharge;
  state.pid = pid;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.01f * (row + 1) : 0.0002f * (row + column + 1);
    }
  }
  return state;
}

SurfaceKinematicState makeForwardState(uint8_t absCharge = 1, float q2pt = 0.8f, o2::track::PID pid = o2::track::PID::Pion)
{
  SurfaceKinematicState state{};
  state.parameters[0] = 1.25f;  // X
  state.parameters[1] = -0.75f; // Y
  state.parameters[2] = 0.35f;  // Phi
  state.parameters[3] = -2.5f;  // Tanl
  state.parameters[4] = q2pt;   // Q2Pt
  state.referenceCoordinate = -45.f;
  state.alpha = 0.f;
  state.kind = SurfaceKind::Disk;
  state.flags = 0x5a;
  state.absCharge = absCharge;
  state.pid = pid;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.01f * (row + 1) : 0.0002f * (row + column + 1);
    }
  }
  return state;
}

// Rebuilds slot 4 so the physical momentum derived from the state's own
// (unchanged) slope equals targetMomentumGeV, for the given absCharge.
SurfaceKinematicState withMomentum(SurfaceKinematicState state, float targetMomentumGeV, uint8_t absCharge)
{
  const float t = state.parameters[3];
  const float pT = targetMomentumGeV / std::sqrt(1.f + t * t);
  state.absCharge = absCharge;
  state.parameters[4] = (absCharge == 0) ? (1.f / pT) : (static_cast<float>(absCharge) / pT);
  return state;
}

SurfaceKinematicState makeValidBarrel() { return makeBarrelState(); }
SurfaceKinematicState makeValidForward() { return makeForwardState(); }

using CorrectFn = material::MaterialOperationResult (*)(SurfaceKinematicState&, material::IntegratedMaterialBudget,
                                                        material::MaterialTraversalDirection);
using MakeValidFn = SurfaceKinematicState (*)();

struct Family {
  const char* name;
  CorrectFn correct;
  MakeValidFn makeValid;
};

const Family kFamilies[] = {
  {"barrel", &barrel::correctForMaterial, &makeValidBarrel},
  {"forward", &forward::correctForMaterial, &makeValidForward},
};

void checkPreflightFailure(const material::MaterialOperationResult& result, material::MaterialFailureReason expected)
{
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.failure == expected);
  BOOST_CHECK_EQUAL(result.momentumBeforeGeV, 0.f);
  BOOST_CHECK_EQUAL(result.momentumAfterGeV, 0.f);
  BOOST_CHECK_EQUAL(result.signedEnergyChangeGeV, 0.f);
  BOOST_CHECK_EQUAL(result.highlandTheta2Rad2, 0.f);
  BOOST_CHECK_EQUAL(result.relativeInverseMomentumVariance, 0.f);
  BOOST_CHECK_EQUAL(result.energyLossSubsteps, 0);
  BOOST_CHECK(result.flags == material::MaterialOperationFlags::None);
  BOOST_CHECK_EQUAL(result.reserved, 0);
}

o2::track::TrackParCovFwd makeForwardOracle(const SurfaceKinematicState& state)
{
  o2::track::SMatrix5 parameters{};
  o2::track::SMatrix55Sym covariance{};
  for (uint8_t i = 0; i < 5; ++i) {
    parameters(i) = state.parameters[i];
  }
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      covariance(row, column) = state.covariance[packedCovarianceIndex(row, column)];
    }
  }
  return {state.referenceCoordinate, parameters, covariance, 0.};
}

} // namespace

BOOST_AUTO_TEST_CASE(ValidOperationsSucceed)
{
  material::IntegratedMaterialBudget materialBudget{0.02f, 0.03f};
  for (const auto& family : kFamilies) {
    auto state = family.makeValid();
    auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    BOOST_CHECK_MESSAGE(result.ok(), family.name << " failed with reason " << static_cast<int>(result.failure));
  }
}

BOOST_AUTO_TEST_CASE(FamilyMismatchBothDirections)
{
  material::IntegratedMaterialBudget materialBudget{0.02f, 0.01f};
  auto barrelState = makeBarrelState();
  const auto barrelBefore = barrelState;
  auto resultAsForward = forward::correctForMaterial(barrelState, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
  checkPreflightFailure(resultAsForward, material::MaterialFailureReason::SourceSurfaceKindMismatch);
  BOOST_CHECK(bitEqual(barrelState, barrelBefore));

  auto forwardState = makeForwardState();
  const auto forwardBefore = forwardState;
  auto resultAsBarrel = barrel::correctForMaterial(forwardState, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
  checkPreflightFailure(resultAsBarrel, material::MaterialFailureReason::SourceSurfaceKindMismatch);
  BOOST_CHECK(bitEqual(forwardState, forwardBefore));
}

BOOST_AUTO_TEST_CASE(NonFiniteStateRejectedInEveryLocation)
{
  const float nan = std::numeric_limits<float>::quiet_NaN();
  material::IntegratedMaterialBudget materialBudget{0.02f, 0.01f};
  for (const auto& family : kFamilies) {
    for (uint8_t i = 0; i < 5; ++i) {
      auto state = family.makeValid();
      state.parameters[i] = nan;
      const auto before = state;
      auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
      checkPreflightFailure(result, material::MaterialFailureReason::NonFiniteState);
      BOOST_CHECK(bitEqual(state, before));
    }
    for (uint8_t i = 0; i < 15; ++i) {
      auto state = family.makeValid();
      state.covariance[i] = nan;
      const auto before = state;
      auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
      checkPreflightFailure(result, material::MaterialFailureReason::NonFiniteState);
      BOOST_CHECK(bitEqual(state, before));
    }
    {
      auto state = family.makeValid();
      state.referenceCoordinate = nan;
      const auto before = state;
      auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
      checkPreflightFailure(result, material::MaterialFailureReason::NonFiniteState);
      BOOST_CHECK(bitEqual(state, before));
    }
    {
      // alpha is checked for finiteness before the forward-specific "alpha ==
      // 0" kinematics check, so NaN alpha must surface NonFiniteState even for
      // forward states.
      auto state = family.makeValid();
      state.alpha = nan;
      const auto before = state;
      auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
      checkPreflightFailure(result, material::MaterialFailureReason::NonFiniteState);
      BOOST_CHECK(bitEqual(state, before));
    }
  }
}

BOOST_AUTO_TEST_CASE(BarrelSnpBoundaryRejection)
{
  material::IntegratedMaterialBudget materialBudget{0.f, 0.f};
  for (float snp : {1.f, -1.f, 1.0001f, -1.0001f}) {
    auto state = makeBarrelState();
    state.parameters[2] = snp;
    const auto before = state;
    auto result = barrel::correctForMaterial(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    checkPreflightFailure(result, material::MaterialFailureReason::InvalidStateKinematics);
    BOOST_CHECK(bitEqual(state, before));
  }
  auto state = makeBarrelState();
  state.parameters[2] = 0.999999f;
  auto result = barrel::correctForMaterial(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
  BOOST_CHECK(result.ok());
}

BOOST_AUTO_TEST_CASE(ForwardNonzeroAlphaRejected)
{
  material::IntegratedMaterialBudget materialBudget{0.f, 0.f};
  for (float alpha : {1.e-6f, -0.3f, 3.14f}) {
    auto state = makeForwardState();
    state.alpha = alpha;
    const auto before = state;
    auto result = forward::correctForMaterial(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    checkPreflightFailure(result, material::MaterialFailureReason::InvalidStateKinematics);
    BOOST_CHECK(bitEqual(state, before));
  }
}

BOOST_AUTO_TEST_CASE(ZeroSlot4Rejected)
{
  material::IntegratedMaterialBudget materialBudget{0.f, 0.f};
  for (const auto& family : kFamilies) {
    auto state = family.makeValid();
    state.parameters[4] = 0.f;
    const auto before = state;
    auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    checkPreflightFailure(result, material::MaterialFailureReason::InvalidStateKinematics);
    BOOST_CHECK(bitEqual(state, before));
  }
}

BOOST_AUTO_TEST_CASE(InvalidPidRejectedBeforeMassLookup)
{
  material::IntegratedMaterialBudget materialBudget{0.f, 0.f};
  for (const auto& family : kFamilies) {
    auto state = family.makeValid();
    state.pid = o2::track::PID(static_cast<o2::track::PID::ID>(255));
    const auto before = state;
    auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    checkPreflightFailure(result, material::MaterialFailureReason::InvalidPID);
    BOOST_CHECK(bitEqual(state, before));
  }
}

BOOST_AUTO_TEST_CASE(ChargedMasslessPidRejected)
{
  material::IntegratedMaterialBudget materialBudget{0.f, 0.f};
  for (const auto& family : kFamilies) {
    auto state = family.makeValid(); // absCharge defaults to 1
    state.pid = o2::track::PID::Photon;
    const auto before = state;
    auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    checkPreflightFailure(result, material::MaterialFailureReason::ChargedMasslessPID);
    BOOST_CHECK(bitEqual(state, before));
  }
}

BOOST_AUTO_TEST_CASE(NegativeCovarianceDiagonalRejected)
{
  material::IntegratedMaterialBudget materialBudget{0.01f, 0.01f};
  for (const auto& family : kFamilies) {
    for (uint8_t i = 0; i < 5; ++i) {
      auto state = family.makeValid();
      state.covariance[packedCovarianceIndex(i, i)] = -1.f;
      const auto before = state;
      auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
      checkPreflightFailure(result, material::MaterialFailureReason::InvalidCovariance);
      BOOST_CHECK(bitEqual(state, before));
    }
  }
}

BOOST_AUTO_TEST_CASE(PhysicalMomentumDerivation)
{
  material::IntegratedMaterialBudget materialBudget{0.f, 0.f}; // zero material: momentumAfter == momentumBefore exactly
  struct Case {
    uint8_t absCharge;
    float slot4;
  };
  const Case cases[] = {{0, 0.5f}, {0, -0.5f}, {1, 0.8f}, {1, -0.8f}, {2, 0.4f}, {3, 0.25f}};
  for (const auto& family : kFamilies) {
    for (const auto& c : cases) {
      auto state = family.makeValid();
      const float t = state.parameters[3];
      state.absCharge = c.absCharge;
      state.parameters[4] = c.slot4;
      const float absU = std::abs(c.slot4);
      const float expectedPt = (c.absCharge == 0) ? (1.f / absU) : (static_cast<float>(c.absCharge) / absU);
      const float expectedP = expectedPt * std::sqrt(1.f + t * t);
      auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
      BOOST_REQUIRE_MESSAGE(result.ok(), family.name << " absCharge=" << static_cast<int>(c.absCharge) << " slot4=" << c.slot4);
      BOOST_CHECK(closeTo(result.momentumBeforeGeV, expectedP));
      BOOST_CHECK(closeTo(result.momentumAfterGeV, expectedP));
    }
  }
}

BOOST_AUTO_TEST_CASE(ZeroMaterialIsByteIdentical)
{
  material::IntegratedMaterialBudget materialBudget{0.f, 0.f};
  for (const auto& family : kFamilies) {
    for (uint8_t absCharge : {0, 1, 2}) {
      auto state = family.makeValid();
      state.absCharge = absCharge;
      const auto before = state;
      auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
      BOOST_REQUIRE(result.ok());
      BOOST_CHECK(bitEqual(state, before));
    }
  }
}

BOOST_AUTO_TEST_CASE(NeutralNonzeroMaterialIsByteIdentical)
{
  material::IntegratedMaterialBudget materialBudget{0.2f, 5.f}; // would matter if charged
  for (const auto& family : kFamilies) {
    auto state = family.makeValid();
    state.absCharge = 0;
    const auto before = state;
    auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    BOOST_REQUIRE(result.ok());
    BOOST_CHECK(bitEqual(state, before));
    BOOST_CHECK_EQUAL(result.momentumAfterGeV, result.momentumBeforeGeV);
  }
}

BOOST_AUTO_TEST_CASE(ScalarKernelFailuresPropagateWithoutMutation)
{
  for (const auto& family : kFamilies) {
    // InvalidDirection
    {
      auto state = withMomentum(family.makeValid(), 1.f, 1);
      const auto before = state;
      auto badDirection = static_cast<material::MaterialTraversalDirection>(255);
      material::IntegratedMaterialBudget materialBudget{0.01f, 0.01f};
      auto result = family.correct(state, materialBudget, badDirection);
      BOOST_CHECK(result.failure == material::MaterialFailureReason::InvalidDirection);
      BOOST_CHECK(bitEqual(state, before));
    }
    // InvalidMaterial
    {
      auto state = withMomentum(family.makeValid(), 1.f, 1);
      const auto before = state;
      material::IntegratedMaterialBudget materialBudget{-1.f, 0.01f};
      auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
      BOOST_CHECK(result.failure == material::MaterialFailureReason::InvalidMaterial);
      BOOST_CHECK(bitEqual(state, before));
    }
    // MomentumBelowMinimum: input momentum > 0 but below the kernel's 0.01 GeV floor.
    {
      auto state = withMomentum(family.makeValid(), 0.001f, 1);
      state.pid = o2::track::PID::Pion;
      const auto before = state;
      material::IntegratedMaterialBudget materialBudget{0.f, 0.f};
      auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
      BOOST_CHECK(result.failure == material::MaterialFailureReason::MomentumBelowMinimum);
      BOOST_CHECK(bitEqual(state, before));
    }
    // StoppedInMaterial: grossly exceeds a 0.5 GeV/c proton's kinetic energy.
    {
      auto state = withMomentum(family.makeValid(), 0.5f, 1);
      state.pid = o2::track::PID::Proton;
      const auto before = state;
      material::IntegratedMaterialBudget materialBudget{0.f, 50.f};
      auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
      BOOST_CHECK(result.failure == material::MaterialFailureReason::StoppedInMaterial);
      BOOST_CHECK(bitEqual(state, before));
    }
    // ExcessiveScattering: absurdly thick, drives theta^2 past pi^2.
    {
      auto state = withMomentum(family.makeValid(), 0.1f, 1);
      state.pid = o2::track::PID::Pion;
      const auto before = state;
      material::IntegratedMaterialBudget materialBudget{500.f, 0.f};
      auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
      BOOST_CHECK(result.failure == material::MaterialFailureReason::ExcessiveScattering);
      BOOST_CHECK(bitEqual(state, before));
    }
    // NonFiniteResult: entry-derived total energy overflows float range.
    {
      auto state = withMomentum(family.makeValid(), std::numeric_limits<float>::max() / 4.f, 1);
      state.pid = o2::track::PID::Proton;
      const auto before = state;
      material::IntegratedMaterialBudget materialBudget{0.01f, 1.e10f};
      auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
      BOOST_CHECK(result.failure == material::MaterialFailureReason::NonFiniteResult);
      BOOST_CHECK(bitEqual(state, before));
    }
  }
}

BOOST_AUTO_TEST_CASE(LossAndGainUpdateQ2PtInCorrectDirectionPreservingPolarity)
{
  for (const auto& family : kFamilies) {
    auto state = withMomentum(family.makeValid(), 1.5f, 1);
    state.pid = o2::track::PID::Proton;
    material::IntegratedMaterialBudget materialBudget{0.f, 0.01f};
    const float kBefore = state.parameters[4];

    auto lossState = state;
    auto lossResult = family.correct(lossState, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    BOOST_REQUIRE(lossResult.ok());
    BOOST_CHECK_LT(lossResult.momentumAfterGeV, lossResult.momentumBeforeGeV);
    BOOST_CHECK_GT(std::abs(lossState.parameters[4]), std::abs(kBefore));
    BOOST_CHECK_EQUAL(std::signbit(lossState.parameters[4]), std::signbit(kBefore)); // polarity preserved

    auto gainState = state;
    auto gainResult = family.correct(gainState, materialBudget, material::MaterialTraversalDirection::OppositeMomentum);
    BOOST_REQUIRE(gainResult.ok());
    BOOST_CHECK_GT(gainResult.momentumAfterGeV, gainResult.momentumBeforeGeV);
    BOOST_CHECK_LT(std::abs(gainState.parameters[4]), std::abs(kBefore));
    BOOST_CHECK_EQUAL(std::signbit(gainState.parameters[4]), std::signbit(kBefore)); // polarity preserved
  }
}

BOOST_AUTO_TEST_CASE(CovarianceProjectionExactIndices)
{
  for (const auto& family : kFamilies) {
    auto state = withMomentum(family.makeValid(), 2.f, 1);
    state.pid = o2::track::PID::Kaon;
    material::IntegratedMaterialBudget materialBudget{0.03f, 0.02f};
    const auto before = state;
    const float tBefore = before.parameters[3];
    const float kBefore = before.parameters[4];
    const float pTBefore = static_cast<float>(before.absCharge) / std::abs(kBefore);
    const float pBefore = pTBefore * std::sqrt(1.f + tBefore * tBefore);

    auto scalarReference = material::calculateMaterialPhysics(pBefore, before.pid, before.absCharge,
                                                              material::MaterialTraversalDirection::AlongMomentum, materialBudget);
    BOOST_REQUIRE(scalarReference.ok());
    const float h = scalarReference.highlandTheta2Rad2;
    const float R = scalarReference.relativeInverseMomentumVariance;
    const float A = 1.f + tBefore * tBefore;

    auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    BOOST_REQUIRE(result.ok());

    const bool isBarrel = (before.kind == SurfaceKind::Cylinder);
    const float snp = before.parameters[2];
    const float c2 = 1.f - snp * snp;
    const float expectedIndex5Increment = isBarrel ? h * A * c2 : h * A;
    const float expectedIndex9Increment = h * A * A;
    const float expectedIndex13Increment = h * A * tBefore * kBefore;
    const float expectedIndex14Increment = h * (tBefore * kBefore) * (tBefore * kBefore) + kBefore * kBefore * R;

    BOOST_CHECK(closeTo(state.covariance[5] - before.covariance[5], expectedIndex5Increment));
    BOOST_CHECK(closeTo(state.covariance[9] - before.covariance[9], expectedIndex9Increment));
    BOOST_CHECK(closeTo(state.covariance[13] - before.covariance[13], expectedIndex13Increment));
    BOOST_CHECK(closeTo(state.covariance[14] - before.covariance[14], expectedIndex14Increment));

    for (uint8_t idx = 0; idx < 15; ++idx) {
      if (idx == 5 || idx == 9 || idx == 13 || idx == 14) {
        continue;
      }
      BOOST_CHECK_EQUAL(state.covariance[idx], before.covariance[idx]);
    }
    BOOST_CHECK_EQUAL(state.parameters[0], before.parameters[0]);
    BOOST_CHECK_EQUAL(state.parameters[1], before.parameters[1]);
    BOOST_CHECK_EQUAL(state.parameters[2], before.parameters[2]);
    BOOST_CHECK_EQUAL(state.parameters[3], before.parameters[3]);
  }
}

BOOST_AUTO_TEST_CASE(Slot13SignFollowsTAndK)
{
  material::IntegratedMaterialBudget materialBudget{0.05f, 0.f}; // MCS-only, guarantees h > 0
  for (const auto& family : kFamilies) {
    for (float tSign : {1.f, -1.f}) {
      for (float kSign : {1.f, -1.f}) {
        auto state = family.makeValid();
        state.parameters[3] = tSign * std::abs(state.parameters[3]);
        state.parameters[4] = kSign * std::abs(state.parameters[4]);
        const auto before = state;
        auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
        BOOST_REQUIRE(result.ok());
        const float increment13 = state.covariance[13] - before.covariance[13];
        if (tSign * kSign > 0.f) {
          BOOST_CHECK_GT(increment13, 0.f);
        } else {
          BOOST_CHECK_LT(increment13, 0.f);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(StragglingContributesToSlot14)
{
  for (const auto& family : kFamilies) {
    auto state = withMomentum(family.makeValid(), 2.f, 1);
    state.pid = o2::track::PID::Proton;
    material::IntegratedMaterialBudget materialBudget{0.f, 0.02f}; // energy-loss only: h == 0, R > 0
    const auto before = state;
    const float kBefore = before.parameters[4];
    auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    BOOST_REQUIRE(result.ok());
    BOOST_CHECK_EQUAL(result.highlandTheta2Rad2, 0.f);
    BOOST_CHECK_GT(result.relativeInverseMomentumVariance, 0.f);
    const float expectedIncrement14 = kBefore * kBefore * result.relativeInverseMomentumVariance;
    BOOST_CHECK(closeTo(state.covariance[14] - before.covariance[14], expectedIncrement14));
    BOOST_CHECK_EQUAL(state.covariance[13], before.covariance[13]); // h == 0: slot13 increment is zero
  }
}

BOOST_AUTO_TEST_CASE(PidAndAbsChargeAreIndependent)
{
  // PID::Electron has a nominal PID::getCharge() of 1, deliberately different
  // from absCharge here; the momentum derivation must use absCharge alone
  // (PID::getCharge() == absCharge is never required).
  material::IntegratedMaterialBudget materialBudget{0.f, 0.f};
  for (const auto& family : kFamilies) {
    auto state = family.makeValid();
    state.pid = o2::track::PID::Electron;
    state.absCharge = 3;
    const float t = state.parameters[3];
    const float u = state.parameters[4];
    const float expectedP = (3.f / std::abs(u)) * std::sqrt(1.f + t * t);
    auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    BOOST_REQUIRE(result.ok());
    BOOST_CHECK(closeTo(result.momentumBeforeGeV, expectedP));
  }
}

BOOST_AUTO_TEST_CASE(MetadataPreservedOnSuccess)
{
  material::IntegratedMaterialBudget materialBudget{0.02f, 0.01f};
  for (const auto& family : kFamilies) {
    auto state = withMomentum(family.makeValid(), 1.5f, 2);
    state.pid = o2::track::PID::Kaon;
    state.flags = 0x7b;
    const auto before = state;
    auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    BOOST_REQUIRE(result.ok());
    BOOST_CHECK_EQUAL(static_cast<uint8_t>(state.kind), static_cast<uint8_t>(before.kind));
    BOOST_CHECK_EQUAL(state.flags, before.flags);
    BOOST_CHECK_EQUAL(state.absCharge, before.absCharge);
    BOOST_CHECK_EQUAL(static_cast<uint8_t>(state.pid), static_cast<uint8_t>(before.pid));
    BOOST_CHECK_EQUAL(state.referenceCoordinate, before.referenceCoordinate);
    BOOST_CHECK_EQUAL(state.alpha, before.alpha);
  }
}

BOOST_AUTO_TEST_CASE(RepeatedCallsAreDeterministic)
{
  material::IntegratedMaterialBudget materialBudget{0.03f, 0.02f};
  for (const auto& family : kFamilies) {
    auto stateA = withMomentum(family.makeValid(), 1.3f, 1);
    stateA.pid = o2::track::PID::Kaon;
    auto stateB = stateA;
    auto resultA = family.correct(stateA, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    auto resultB = family.correct(stateB, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    BOOST_CHECK(bitEqual(resultA, resultB));
    BOOST_CHECK(bitEqual(stateA, stateB));

    auto failA = family.makeValid();
    failA.kind = SurfaceKind::Undefined;
    auto failB = failA;
    auto rA = family.correct(failA, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    auto rB = family.correct(failB, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    BOOST_CHECK(bitEqual(rA, rB));
    BOOST_CHECK(bitEqual(failA, failB));
  }
}

BOOST_AUTO_TEST_CASE(BarrelIntentionalUnitChargeCorrectionIsolated)
{
  // absCharge == 1 exercises the retained TrackParametrizationWithError<float>::
  // correctForMaterial() bug: it multiplies the Tgl/Q2Pt Jacobian factor (and
  // its self-product) by getCharge2Pt() only when getAbsCharge()^2 != 1,
  // silently omitting q/pT from covariance slots 13/14 for unit charge. This
  // operation's accepted Jacobian includes q/pT unconditionally.
  auto state = withMomentum(makeBarrelState(), 1.2f, 1);
  state.pid = o2::track::PID::Kaon;
  // xOverX0 is deliberately large (but still far below the ExcessiveScattering
  // threshold) so the intentional slot-13/14 divergence is well above float
  // characterization noise.
  material::IntegratedMaterialBudget materialBudget{20.f, 0.03f};
  const auto before = state;

  o2::track::TrackParCovF oracle;
  BOOST_REQUIRE(legacy::exportBarrelTrackParCov(before, oracle));
  const float legacyXrho = -materialBudget.arealDensityGPerCm2; // AlongMomentum == energy loss == negative xrho
  BOOST_REQUIRE(oracle.correctForMaterial(materialBudget.xOverX0, legacyXrho, false));

  auto result = barrel::correctForMaterial(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
  BOOST_REQUIRE(result.ok());

  // Snp/Tgl diagonal terms are unaffected by the unit-charge bug.
  BOOST_CHECK(closeTo(state.covariance[5], oracle.getSigmaSnp2()));
  BOOST_CHECK(closeTo(state.covariance[9], oracle.getSigmaTgl2()));

  // Slots 13/14 intentionally diverge from the retained oracle for absCharge == 1.
  BOOST_CHECK_GT(std::fabs(state.covariance[13] - oracle.getSigma1PtTgl()), 1.e-5f);

  const float t = before.parameters[3];
  const float k = before.parameters[4];
  const float A = 1.f + t * t;
  auto scalarReference = material::calculateMaterialPhysics(result.momentumBeforeGeV, before.pid, before.absCharge,
                                                            material::MaterialTraversalDirection::AlongMomentum, materialBudget);
  BOOST_REQUIRE(scalarReference.ok());
  const float expectedIncrement13 = scalarReference.highlandTheta2Rad2 * A * t * k;
  BOOST_CHECK(closeTo(state.covariance[13], before.covariance[13] + expectedIncrement13));
}

BOOST_AUTO_TEST_CASE(BarrelMatchesLegacyOracleForChargeAboveUnity)
{
  // absCharge == 2 does not hit the retained bug (getAbsCharge()^2 != 1), so
  // the retained oracle and this operation should agree on all four
  // covariance slots within float characterization tolerance.
  auto state = withMomentum(makeBarrelState(), 1.2f, 2);
  state.pid = o2::track::PID::Kaon;
  material::IntegratedMaterialBudget materialBudget{0.02f, 0.03f};
  const auto before = state;

  o2::track::TrackParCovF oracle;
  BOOST_REQUIRE(legacy::exportBarrelTrackParCov(before, oracle));
  const float legacyXrho = -materialBudget.arealDensityGPerCm2;
  BOOST_REQUIRE(oracle.correctForMaterial(materialBudget.xOverX0, legacyXrho, false));

  auto result = barrel::correctForMaterial(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
  BOOST_REQUIRE(result.ok());

  BOOST_CHECK(closeTo(state.covariance[5], oracle.getSigmaSnp2(), 5.e-4f, 3.e-3f));
  BOOST_CHECK(closeTo(state.covariance[9], oracle.getSigmaTgl2(), 5.e-4f, 3.e-3f));
  BOOST_CHECK(closeTo(state.covariance[13], oracle.getSigma1PtTgl(), 5.e-4f, 3.e-3f));
  BOOST_CHECK(closeTo(state.covariance[14], oracle.getSigma1Pt2(), 5.e-4f, 3.e-3f));
  BOOST_CHECK(closeTo(state.parameters[4], oracle.getQ2Pt(), 5.e-4f, 3.e-3f));
}

BOOST_AUTO_TEST_CASE(BarrelCovarianceLimitingActivates)
{
  auto state = withMomentum(makeBarrelState(), 0.05f, 1);
  state.pid = o2::track::PID::Pion;
  state.covariance[packedCovarianceIndex(2, 2)] = 0.99f;       // Snp diag, close to kCSnp2max
  state.covariance[packedCovarianceIndex(3, 3)] = 0.99f;       // Tgl diag, close to kCTgl2max
  state.covariance[packedCovarianceIndex(4, 3)] = 0.02f;       // slot 13 seed, to observe rescaling
  material::IntegratedMaterialBudget materialBudget{1.f, 0.f}; // MCS-only, large enough to blow past both thresholds

  const float t = state.parameters[3];
  const float k = state.parameters[4];
  const float pT = 1.f / std::abs(k);
  const float pBefore = pT * std::sqrt(1.f + t * t);
  auto scalarReference = material::calculateMaterialPhysics(pBefore, state.pid, state.absCharge,
                                                            material::MaterialTraversalDirection::AlongMomentum, materialBudget);
  BOOST_REQUIRE(scalarReference.ok());
  const float A = 1.f + t * t;
  const float snp = state.parameters[2];
  const float c2 = 1.f - snp * snp;
  const float rawIndex5 = 0.99f + scalarReference.highlandTheta2Rad2 * A * c2;
  const float rawIndex9 = 0.99f + scalarReference.highlandTheta2Rad2 * A * A;
  const float rawIndex13 = 0.02f + scalarReference.highlandTheta2Rad2 * A * t * k;
  BOOST_REQUIRE_GT(rawIndex5, o2::track::kCSnp2max);
  BOOST_REQUIRE_GT(rawIndex9, o2::track::kCTgl2max);

  auto result = barrel::correctForMaterial(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
  BOOST_REQUIRE(result.ok());

  BOOST_CHECK_EQUAL(state.covariance[packedCovarianceIndex(2, 2)], o2::track::kCSnp2max);
  BOOST_CHECK_EQUAL(state.covariance[packedCovarianceIndex(3, 3)], o2::track::kCTgl2max);
  // The clamp must have rescaled the dependent Q2Pt/Tgl cross term downward
  // from what the raw (unclamped) projection would have produced.
  BOOST_CHECK_LT(std::fabs(state.covariance[packedCovarianceIndex(4, 3)]), std::fabs(rawIndex13));
}

BOOST_AUTO_TEST_CASE(ForwardDoesNotInheritBarrelCovarianceLimiting)
{
  auto state = withMomentum(makeForwardState(), 0.05f, 1);
  state.pid = o2::track::PID::Pion;
  state.covariance[packedCovarianceIndex(2, 2)] = 0.99f;
  state.covariance[packedCovarianceIndex(3, 3)] = 0.99f;
  material::IntegratedMaterialBudget materialBudget{1.f, 0.f};

  auto result = forward::correctForMaterial(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_GT(state.covariance[packedCovarianceIndex(2, 2)], o2::track::kCSnp2max);
  BOOST_CHECK_GT(state.covariance[packedCovarianceIndex(3, 3)], o2::track::kCTgl2max);
  BOOST_CHECK(std::isfinite(state.covariance[packedCovarianceIndex(2, 2)]));
  BOOST_CHECK(std::isfinite(state.covariance[packedCovarianceIndex(3, 3)]));
}

BOOST_AUTO_TEST_CASE(ForwardMatchesLegacyMCSUltraRelativisticUnitCharge)
{
  auto state = withMomentum(makeForwardState(), 20.f, 1); // ultra-relativistic: beta^2 ~ 1
  state.pid = o2::track::PID::Pion;
  material::IntegratedMaterialBudget materialBudget{0.05f, 0.f}; // zero areal density: MCS only
  const auto before = state;

  auto oracle = makeForwardOracle(before);
  const float tanl = before.parameters[3];
  const float cscLambda = std::fabs(std::sqrt(1.f + tanl * tanl) / tanl);
  // TrackParCovFwd::addMCSEffect() internally multiplies its argument by
  // cscLambda; undo that here since our material budget is already
  // path-integrated.
  const float legacyXOverX0 = materialBudget.xOverX0 / cscLambda;
  oracle.addMCSEffect(legacyXOverX0);

  auto result = forward::correctForMaterial(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
  BOOST_REQUIRE(result.ok());

  BOOST_CHECK(closeTo(state.covariance[5], static_cast<float>(oracle.getCovariances()(2, 2)), 1.e-4f, 5.e-3f));
  BOOST_CHECK(closeTo(state.covariance[9], static_cast<float>(oracle.getCovariances()(3, 3)), 1.e-4f, 5.e-3f));
}

BOOST_AUTO_TEST_CASE(ForwardAnalyticReferenceChargeAboveUnityWithEnergyLoss)
{
  // This does not compare to any legacy forward object: TrackParCovFwd has no
  // charge/PID/energy-loss awareness. The reference here is the documented
  // projection formula, independently recomputed from the separately
  // characterized scalar kernel (see testMaterialPhysics.cxx).
  auto state = withMomentum(makeForwardState(), 1.5f, 3);
  state.pid = o2::track::PID::Proton;
  material::IntegratedMaterialBudget materialBudget{0.02f, 0.01f};
  const auto before = state;
  const float t = before.parameters[3];
  const float k = before.parameters[4];
  const float pT = 3.f / std::abs(k);
  const float pBefore = pT * std::sqrt(1.f + t * t);

  auto scalarReference = material::calculateMaterialPhysics(pBefore, before.pid, before.absCharge,
                                                            material::MaterialTraversalDirection::AlongMomentum, materialBudget);
  BOOST_REQUIRE(scalarReference.ok());
  const float A = 1.f + t * t;
  const float h = scalarReference.highlandTheta2Rad2;
  const float R = scalarReference.relativeInverseMomentumVariance;

  auto result = forward::correctForMaterial(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
  BOOST_REQUIRE(result.ok());

  BOOST_CHECK(closeTo(state.covariance[5], before.covariance[5] + h * A));
  BOOST_CHECK(closeTo(state.covariance[9], before.covariance[9] + h * A * A));
  BOOST_CHECK(closeTo(state.covariance[13], before.covariance[13] + h * A * t * k));
  BOOST_CHECK(closeTo(state.covariance[14], before.covariance[14] + h * (t * k) * (t * k) + k * k * R));
  BOOST_CHECK_LT(result.momentumAfterGeV, result.momentumBeforeGeV);
}

BOOST_AUTO_TEST_CASE(UnconditionalNoopLeavesOverLimitBarrelCovarianceUntouched)
{
  // Once the scalar kernel succeeds, absCharge == 0 or an exactly-{0,0}
  // materialBudget must return immediately, before projectCovariance (and
  // the barrel covariance-range limiting it applies) or the slot-4 update
  // ever run -- even when the source covariance diagonals already exceed
  // the retained checkCovariance limits, which must not be silently
  // clamped by an operation that has no material to apply.
  for (auto direction : {material::MaterialTraversalDirection::AlongMomentum, material::MaterialTraversalDirection::OppositeMomentum}) {
    {
      auto state = withMomentum(makeBarrelState(), 1.2f, 1);
      state.pid = o2::track::PID::Pion;
      state.covariance[packedCovarianceIndex(2, 2)] = 2.f * o2::track::kCSnp2max;
      state.covariance[packedCovarianceIndex(3, 3)] = 2.f * o2::track::kCTgl2max;
      const auto before = state;
      material::IntegratedMaterialBudget materialBudget{0.f, 0.f}; // zero material, charged
      auto result = barrel::correctForMaterial(state, materialBudget, direction);
      BOOST_REQUIRE(result.ok());
      BOOST_CHECK(bitEqual(state, before));
    }
    {
      auto state = withMomentum(makeBarrelState(), 1.2f, 0);
      state.pid = o2::track::PID::Pion;
      state.covariance[packedCovarianceIndex(2, 2)] = 2.f * o2::track::kCSnp2max;
      state.covariance[packedCovarianceIndex(3, 3)] = 2.f * o2::track::kCTgl2max;
      const auto before = state;
      material::IntegratedMaterialBudget materialBudget{0.5f, 3.f}; // nonzero, neutral (would matter if charged)
      auto result = barrel::correctForMaterial(state, materialBudget, direction);
      BOOST_REQUIRE(result.ok());
      BOOST_CHECK(bitEqual(state, before));
    }
  }
}

BOOST_AUTO_TEST_CASE(Slot4ZeroChangePreservesBitForBit)
{
  // MCS-only material (xOverX0 > 0, arealDensity == 0) leaves momentum
  // exactly unchanged even though covariance is projected; the equality
  // branch of the slot-4 update must recover kBefore bit-for-bit here.
  for (const auto& family : kFamilies) {
    auto state = withMomentum(family.makeValid(), 2.f, 1);
    state.pid = o2::track::PID::Kaon;
    material::IntegratedMaterialBudget materialBudget{0.05f, 0.f};
    const float kBefore = state.parameters[4];
    auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    BOOST_REQUIRE(result.ok());
    BOOST_CHECK_EQUAL(result.momentumBeforeGeV, result.momentumAfterGeV);
    BOOST_CHECK_EQUAL(state.parameters[4], kBefore);
  }
}

BOOST_AUTO_TEST_CASE(Slot4OrdinaryLossAndGainMatchExplicitProductFirstFormula)
{
  // For an ordinary (non-equal) momentum change, the retained/accepted
  // arithmetic is left-to-right: (kBefore * pBefore) / pAfter, not
  // kBefore * (pBefore / pAfter).
  for (const auto& family : kFamilies) {
    for (auto direction : {material::MaterialTraversalDirection::AlongMomentum, material::MaterialTraversalDirection::OppositeMomentum}) {
      auto state = withMomentum(family.makeValid(), 1.5f, 1);
      state.pid = o2::track::PID::Proton;
      material::IntegratedMaterialBudget materialBudget{0.f, 0.01f}; // energy-loss only: pBefore != pAfter
      const float kBefore = state.parameters[4];
      auto result = family.correct(state, materialBudget, direction);
      BOOST_REQUIRE(result.ok());
      BOOST_REQUIRE_NE(result.momentumBeforeGeV, result.momentumAfterGeV);
      const float expectedKAfter = (kBefore * result.momentumBeforeGeV) / result.momentumAfterGeV;
      BOOST_CHECK_EQUAL(state.parameters[4], expectedKAfter);
    }
  }
}

BOOST_AUTO_TEST_CASE(RatioFirstArithmeticUnderflowsWhileSelectedFormulationRemainsFinite)
{
  // Constructing a scalar-kernel-success physical fixture that reaches this
  // exact float-underflow boundary was not possible: any state whose slot 4
  // is small enough (or whose derived momentum is large enough) to approach
  // this regime causes the physical-momentum derivation or the scalar
  // kernel's own entry-energy computation (e0 = sqrt(p^2 + mass^2)) to
  // overflow first, which is rejected earlier as InvalidStateKinematics or
  // NonFiniteResult -- never reaching this arithmetic. This is therefore a
  // direct, scalar-kernel-independent characterization of why the ratio
  // must not be computed before multiplying by kBefore, using the exact
  // formulas the production code chooses between.
  const float kBefore = 1.e25f;
  const float pBefore = 1.e-25f;
  const float pAfter = 1.e25f;

  const float ratioFirst = kBefore * (pBefore / pAfter);
  const float productFirst = (kBefore * pBefore) / pAfter;

  BOOST_CHECK_EQUAL(ratioFirst, 0.f); // pBefore/pAfter underflows to exactly zero first
  BOOST_CHECK(std::isfinite(productFirst));
  BOOST_CHECK_GT(productFirst, 0.f); // (kBefore * pBefore) / pAfter remains finite and nonzero
}

BOOST_AUTO_TEST_CASE(Slot4ArithmeticRepeatedExecutionIsDeterministic)
{
  material::IntegratedMaterialBudget materialBudget{0.f, 0.02f};
  for (const auto& family : kFamilies) {
    for (auto direction : {material::MaterialTraversalDirection::AlongMomentum, material::MaterialTraversalDirection::OppositeMomentum}) {
      auto stateA = withMomentum(family.makeValid(), 1.4f, 1);
      stateA.pid = o2::track::PID::Proton;
      auto stateB = stateA;
      auto resultA = family.correct(stateA, materialBudget, direction);
      auto resultB = family.correct(stateB, materialBudget, direction);
      BOOST_CHECK(bitEqual(resultA, resultB));
      BOOST_CHECK(bitEqual(stateA, stateB));
    }
  }
}

BOOST_AUTO_TEST_CASE(PostProjectionValidationDoesNotFalsePositiveNearBoundary)
{
  // The three new post-projection checks (finite output, slot 4 valid plus
  // family kinematics, re-derivable momentum, non-negative covariance
  // diagonals) are defense-in-depth for the current projection formulas:
  // parameters[0..3] are never written by projectCovariance, so the
  // family-specific kinematics check (which only reads those) cannot
  // legitimately fail post-projection given it already passed preflight,
  // and every covariance increment (h*A*c2, h*A^2, h*A*t*k, h*(t*k)^2+k^2*R)
  // is provably non-negative or, for slot 13, bounded by quantities already
  // required finite by the scalar kernel's own success. A fixture with
  // large-but-legitimate t/k that would overflow the projection arithmetic
  // was found to overflow the scalar kernel's own e0 = sqrt(p^2 + mass^2)
  // computation first (or the physical-momentum derivation itself), so no
  // reachable public fixture drives NonFiniteResult, InvalidStateKinematics,
  // or InvalidCovariance out of the post-projection checks specifically
  // (as opposed to preflight or the scalar kernel). This test instead
  // confirms the checks do not spuriously reject a legitimate large-but-
  // valid boundary fixture.
  for (const auto& family : kFamilies) {
    auto state = withMomentum(family.makeValid(), 5.f, 1);
    state.parameters[3] = 1.e4f; // large but finite slope, well short of overflowing t*t
    state.pid = o2::track::PID::Pion;
    material::IntegratedMaterialBudget materialBudget{5.f, 0.01f}; // sizable but under the ExcessiveScattering cap
    auto result = family.correct(state, materialBudget, material::MaterialTraversalDirection::AlongMomentum);
    BOOST_CHECK_MESSAGE(result.ok(), family.name << " unexpectedly failed with reason " << static_cast<int>(result.failure));
  }
}

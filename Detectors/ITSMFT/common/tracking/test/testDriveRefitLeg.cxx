// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFTDriveRefitLeg
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

#include "ITSMFTTracking/BarrelSurfaceStateOperations.h"
#include "ITSMFTTracking/ForwardSurfaceStateOperations.h"
#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/TransitionPolicyOperations.h"

/// Stage-B refit-leg orchestration slice: focused coverage for
/// driveRefitLeg<Tag> (TransitionPolicyOperations.h), unwired like refitHit
/// (testRefitHit.cxx): there is no production call site anywhere in this
/// slice. Every scenario below is exercised for both CylinderCylinder and
/// DiskDisk through one shared templated test body, following the
/// checkChi2ConfigRejection<Tag> convention already established in
/// testRefitHit.cxx.
using namespace o2::itsmft::tracking;

namespace
{

template <typename T>
bool bitEqual(const T& lhs, const T& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
}

// SurfaceKinematicState has no inter-field padding (92 bytes, 4-byte aligned
// fields summing exactly to 92), so bitEqual is exact for it.
// SurfaceLinearizationReference, by contrast, has 3 trailing padding bytes
// after `family` (28+1=29 bytes of fields, 32-byte aligned size) whose
// content is unspecified after an ordinary (non-memcpy) assignment chain --
// two logically-identical objects reached through differently-shaped call
// paths (e.g. this test's independent oracle vs. the production
// orchestration) can legitimately disagree on those 3 bytes alone. Compare
// the documented fields instead of raw bytes for this type.
bool linRefEqual(const SurfaceLinearizationReference& lhs, const SurfaceLinearizationReference& rhs)
{
  for (int i = 0; i < 5; ++i) {
    if (lhs.parameters[i] != rhs.parameters[i]) {
      return false;
    }
  }
  return lhs.referenceCoordinate == rhs.referenceCoordinate && lhs.alpha == rhs.alpha && lhs.family == rhs.family;
}

// --- Per-Tag fixtures -----------------------------------------------------

template <TransitionPolicyTag Tag>
struct LegTraits;

template <>
struct LegTraits<TransitionPolicyTag::CylinderCylinder> {
  static constexpr float bz = 5.f;
  static constexpr SurfaceKind kind = SurfaceKind::Cylinder;

  static SurfaceKinematicState baseState()
  {
    SurfaceKinematicState state{};
    state.parameters[0] = 1.25f;
    state.parameters[1] = -0.75f;
    state.parameters[2] = 0.2f;
    state.parameters[3] = -0.35f;
    state.parameters[4] = 0.8f;
    state.referenceCoordinate = 4.f;
    state.alpha = 0.3f;
    state.family = StateFamily::Barrel;
    state.absCharge = 1;
    state.pid = o2::track::PID::Pion;
    for (uint8_t row = 0; row < 5; ++row) {
      for (uint8_t column = 0; column <= row; ++column) {
        state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.01f * (row + 1) : 0.0002f * (row + column + 1);
      }
    }
    return state;
  }

  static SurfaceLinearizationReference baseLinRef(const SurfaceKinematicState& state)
  {
    SurfaceLinearizationReference ref{};
    BOOST_REQUIRE(makeLinearizationReference(state, ref));
    return ref;
  }

  static NominalSurfaceMaterial material() { return NominalSurfaceMaterial{0.01f, 0.001f}; }
  static NominalSurfaceMaterial badMaterial() { return NominalSurfaceMaterial{std::numeric_limits<float>::quiet_NaN(), 0.f}; }

  // Same measurement frame for every SurfaceId: frameAngle equals the
  // fixture state's own alpha (no rotation needed), only cluster/surface
  // identity vary between slots.
  static SurfaceMeasurement measurement(SurfaceId id)
  {
    SurfaceMeasurement m{};
    m.frame.q = 2.5f;
    m.frame.frameAngle = 0.3f;
    m.frame.u = 0.8f;
    m.frame.v = -0.45f;
    m.covariance = {0.04f, 0.012f, 0.09f};
    m.cluster = ClusterRef{ClusterSourceId{0}, id.value()};
    m.surface = id;
    return m;
  }

  static SurfaceMeasurement hugeResidualMeasurement(SurfaceId id)
  {
    auto m = measurement(id);
    m.frame.u = 500.f;
    return m;
  }
};

template <>
struct LegTraits<TransitionPolicyTag::DiskDisk> {
  static constexpr float bz = 5.f;
  static constexpr SurfaceKind kind = SurfaceKind::Disk;

  static SurfaceKinematicState baseState()
  {
    SurfaceKinematicState state{};
    state.parameters[0] = 1.25f;
    state.parameters[1] = -0.75f;
    state.parameters[2] = 0.35f;
    state.parameters[3] = -2.5f;
    state.parameters[4] = 0.8f;
    state.referenceCoordinate = -45.f;
    state.family = StateFamily::Forward;
    state.absCharge = 1;
    state.pid = o2::track::PID::Pion;
    for (uint8_t row = 0; row < 5; ++row) {
      for (uint8_t column = 0; column <= row; ++column) {
        state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.01f * (row + 1) : 0.0002f * (row + column + 1);
      }
    }
    return state;
  }

  static SurfaceLinearizationReference baseLinRef(const SurfaceKinematicState& state)
  {
    SurfaceLinearizationReference ref{};
    BOOST_REQUIRE(makeLinearizationReference(state, ref));
    return ref;
  }

  static NominalSurfaceMaterial material() { return NominalSurfaceMaterial{0.02f, 0.002f}; }
  static NominalSurfaceMaterial badMaterial() { return NominalSurfaceMaterial{std::numeric_limits<float>::quiet_NaN(), 0.f}; }

  static SurfaceMeasurement measurement(SurfaceId id)
  {
    SurfaceMeasurement m{};
    m.global = {0.8f, -0.45f, -50.f};
    m.frame.q = -50.f;
    m.frame.u = 0.8f;
    m.frame.v = -0.45f;
    m.covariance = {0.04f, 0.f, 0.09f};
    m.cluster = ClusterRef{ClusterSourceId{0}, id.value()};
    m.surface = id;
    return m;
  }

  static SurfaceMeasurement hugeResidualMeasurement(SurfaceId id)
  {
    auto m = measurement(id);
    m.global.x += 500.f;
    m.frame.u += 500.f;
    return m;
  }
};

// --- Catalog construction ---------------------------------------------------

constexpr uint32_t kNoBadIndex = std::numeric_limits<uint32_t>::max();

template <TransitionPolicyTag Tag>
std::vector<SurfaceDescriptor> makeCatalogSurfaces(uint32_t n, uint32_t badMaterialIndex = kNoBadIndex)
{
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.reserve(n);
  for (uint32_t i = 0; i < n; ++i) {
    SurfaceDescriptor d{};
    d.id = SurfaceId{static_cast<uint16_t>(i)};
    d.detectorSurfaceIndex = static_cast<uint16_t>(i);
    d.kind = LegTraits<Tag>::kind;
    d.material = (i == badMaterialIndex) ? LegTraits<Tag>::badMaterial() : LegTraits<Tag>::material();
    surfaces.push_back(d);
  }
  return surfaces;
}

// --- Independent oracle: a plain re-transcription of driveRefitLeg<Tag>'s own
// documented contract (skip holes, validate surface/catalog association in
// the documented order, gate on a from-zero local count, call refitHit<Tag>
// directly), used as the byte-identical replay oracle for every success case
// below, including out-of-order SurfaceId sequences (processed in the exact
// supplied order, just like the production operation). -------------------

template <TransitionPolicyTag Tag>
struct LegReplayResult {
  bool ok{false};
  SurfaceKinematicState state{};
  SurfaceLinearizationReference linRef{};
  float chi2{0.f};
  uint32_t acceptedHitCount{0};
  OperationFailureReason reason{};
};

template <TransitionPolicyTag Tag>
LegReplayResult<Tag> replayDriveRefitLeg(SurfaceKinematicState state, SurfaceLinearizationReference linRef, float chi2,
                                         gsl::span<const SurfaceMeasurement> orderedSlots, SurfaceCatalogView catalog,
                                         float bz, material::MaterialTraversalDirection direction, bool shiftRef,
                                         float maxChi2)
{
  LegReplayResult<Tag> result{};
  if (!std::isfinite(chi2)) {
    result.reason = OperationFailureReason::NonFiniteInput;
    return result;
  }
  if (chi2 < 0.f) {
    result.reason = OperationFailureReason::PredictedChi2Failure;
    return result;
  }

  uint32_t count = 0;
  for (const auto& measurement : orderedSlots) {
    if (!measurement.cluster.isValid()) {
      continue;
    }
    if (!measurement.surface.isValid()) {
      result.reason = OperationFailureReason::InvalidSurfaceCatalogAssociation;
      return result;
    }
    if (!(catalog.nSurfaces == 0 || catalog.surfaces != nullptr)) {
      result.reason = OperationFailureReason::InvalidSurfaceCatalogAssociation;
      return result;
    }
    if (!(measurement.surface.value() < catalog.nSurfaces)) {
      result.reason = OperationFailureReason::InvalidSurfaceCatalogAssociation;
      return result;
    }
    const SurfaceDescriptor& descriptor = catalog.getSurface(measurement.surface);
    if (!(descriptor.id == measurement.surface)) {
      result.reason = OperationFailureReason::InvalidSurfaceCatalogAssociation;
      return result;
    }
    const bool gate = count >= kRefitLegChi2GateMinAcceptedHits;
    OperationFailureReason reason{};
    if (!refitHit<Tag>(state, linRef, measurement, descriptor.material, bz, direction, gate, maxChi2, chi2, shiftRef, reason)) {
      result.reason = reason;
      return result;
    }
    ++count;
  }

  result.ok = true;
  result.state = state;
  result.linRef = linRef;
  result.chi2 = chi2;
  result.acceptedHitCount = count;
  return result;
}

template <TransitionPolicyTag Tag>
void expectSuccessMatchesOracle(const SurfaceKinematicState& state, const SurfaceLinearizationReference& linRef, float chi2,
                                gsl::span<const SurfaceMeasurement> orderedSlots, SurfaceCatalogView catalog, float bz,
                                material::MaterialTraversalDirection direction, bool shiftRef, float maxChi2,
                                uint32_t sentinelAcceptedHitCount = 0)
{
  auto oracle = replayDriveRefitLeg<Tag>(state, linRef, chi2, orderedSlots, catalog, bz, direction, shiftRef, maxChi2);
  BOOST_REQUIRE(oracle.ok);

  auto actualState = state;
  auto actualLinRef = linRef;
  float actualChi2 = chi2;
  uint32_t actualCount = sentinelAcceptedHitCount;
  OperationFailureReason reason{};
  BOOST_REQUIRE(driveRefitLeg<Tag>(actualState, actualLinRef, actualChi2, actualCount, orderedSlots, catalog, bz,
                                   direction, shiftRef, maxChi2, reason));
  BOOST_CHECK(bitEqual(actualState, oracle.state));
  BOOST_CHECK(linRefEqual(actualLinRef, oracle.linRef));
  BOOST_CHECK_EQUAL(actualChi2, oracle.chi2);
  BOOST_CHECK_EQUAL(actualCount, oracle.acceptedHitCount);
}

template <TransitionPolicyTag Tag>
void expectFailureRollsBack(const SurfaceKinematicState& state, const SurfaceLinearizationReference& linRef, float chi2,
                            gsl::span<const SurfaceMeasurement> orderedSlots, SurfaceCatalogView catalog, float bz,
                            material::MaterialTraversalDirection direction, bool shiftRef, float maxChi2,
                            OperationFailureReason expectedReason, uint32_t sentinelAcceptedHitCount)
{
  auto actualState = state;
  auto actualLinRef = linRef;
  float actualChi2 = chi2;
  uint32_t actualCount = sentinelAcceptedHitCount;
  OperationFailureReason reason{};
  BOOST_CHECK(!driveRefitLeg<Tag>(actualState, actualLinRef, actualChi2, actualCount, orderedSlots, catalog, bz,
                                  direction, shiftRef, maxChi2, reason));
  BOOST_CHECK(reason == expectedReason);
  BOOST_CHECK(bitEqual(actualState, state));
  BOOST_CHECK(linRefEqual(actualLinRef, linRef));
  BOOST_CHECK_EQUAL(actualChi2, chi2);
  BOOST_CHECK_EQUAL(actualCount, sentinelAcceptedHitCount);
}

// Sequential, gate-always-disabled probe used only to discover a leg's
// natural per-hit chi2 contribution (predChi2 == updateChi2 for an accepted
// step, as already established by testRefitHit.cxx's own boundary tests) --
// never itself asserted against driveRefitLeg, only used to derive maxChi2
// thresholds for the gating-boundary test below.
template <TransitionPolicyTag Tag>
std::vector<float> rawSequentialChi2AfterEachHit(SurfaceKinematicState state, SurfaceLinearizationReference linRef,
                                                 gsl::span<const SurfaceMeasurement> measurements,
                                                 const NominalSurfaceMaterial& material, float bz,
                                                 material::MaterialTraversalDirection direction, bool shiftRef)
{
  std::vector<float> chi2AfterEach;
  float chi2 = 0.f;
  OperationFailureReason reason{};
  for (const auto& m : measurements) {
    BOOST_REQUIRE(refitHit<Tag>(state, linRef, m, material, bz, direction, false, 0.f, chi2, shiftRef, reason));
    chi2AfterEach.push_back(chi2);
  }
  return chi2AfterEach;
}

// --- Scenario bodies, shared between CylinderCylinder and DiskDisk --------

template <TransitionPolicyTag Tag>
void runOutputCountStartsFromZeroDespiteSentinel()
{
  using Traits = LegTraits<Tag>;
  const auto baseState = Traits::baseState();
  const auto baseLinRef = Traits::baseLinRef(baseState);
  std::array<SurfaceMeasurement, 1> slots{Traits::hugeResidualMeasurement(SurfaceId{0})};
  auto surfaces = makeCatalogSurfaces<Tag>(1);
  SurfaceCatalogView catalog{surfaces.data(), static_cast<uint32_t>(surfaces.size())};

  auto state = baseState;
  auto linRef = baseLinRef;
  float chi2 = 0.f;
  // A nonzero sentinel at/above kRefitLegChi2GateMinAcceptedHits: if the gate
  // were (incorrectly) seeded from this entry value instead of counting from
  // zero internally, this single huge-residual hit would be wrongly gated
  // against the tiny maxChi2 below and rejected.
  uint32_t acceptedHitCount = 999;
  OperationFailureReason reason{};
  BOOST_REQUIRE(driveRefitLeg<Tag>(state, linRef, chi2, acceptedHitCount, slots, catalog, Traits::bz,
                                   material::MaterialTraversalDirection::AlongMomentum, false, 1.e-6f, reason));
  BOOST_CHECK_EQUAL(acceptedHitCount, 1u);
}

template <TransitionPolicyTag Tag>
void runSentinelUnchangedOnFailure()
{
  using Traits = LegTraits<Tag>;
  const auto baseState = Traits::baseState();
  const auto baseLinRef = Traits::baseLinRef(baseState);
  std::array<SurfaceMeasurement, 1> slots{LegTraits<Tag>::measurement(SurfaceId{5})};
  SurfaceCatalogView emptyCatalog{nullptr, 0}; // out-of-range: value()=5 >= nSurfaces=0

  expectFailureRollsBack<Tag>(baseState, baseLinRef, 0.f, slots, emptyCatalog, Traits::bz,
                              material::MaterialTraversalDirection::AlongMomentum, false, 10.f,
                              OperationFailureReason::InvalidSurfaceCatalogAssociation, 555);
}

template <TransitionPolicyTag Tag>
void runHits1To3UngatedHit4Gated()
{
  using Traits = LegTraits<Tag>;
  const auto baseState = Traits::baseState();
  const auto baseLinRef = Traits::baseLinRef(baseState);
  const auto material = Traits::material();

  std::array<SurfaceMeasurement, 4> measurements{
    Traits::measurement(SurfaceId{0}), Traits::measurement(SurfaceId{1}),
    Traits::measurement(SurfaceId{2}), Traits::measurement(SurfaceId{3})};

  const auto chi2AfterEach = rawSequentialChi2AfterEachHit<Tag>(
    baseState, baseLinRef, measurements, material, Traits::bz, material::MaterialTraversalDirection::AlongMomentum, false);
  const float natural4 = chi2AfterEach[3] - chi2AfterEach[2];
  BOOST_REQUIRE_GT(natural4, 0.f);

  auto surfaces = makeCatalogSurfaces<Tag>(4);
  SurfaceCatalogView catalog{surfaces.data(), static_cast<uint32_t>(surfaces.size())};

  // Accepted: maxChi2 comfortably above hit 4's natural contribution.
  expectSuccessMatchesOracle<Tag>(baseState, baseLinRef, 0.f, measurements, catalog, Traits::bz,
                                  material::MaterialTraversalDirection::AlongMomentum, false, natural4 * 2.f + 1.f);

  // Rejected: maxChi2 strictly below hit 4's natural contribution. Hits 1-3
  // still succeed internally (their own natural chi2 may just as well exceed
  // this same threshold -- they are ungated regardless), only hit 4 is
  // rejected once the gate activates, and the whole leg rolls back.
  expectFailureRollsBack<Tag>(baseState, baseLinRef, 0.f, measurements, catalog, Traits::bz,
                              material::MaterialTraversalDirection::AlongMomentum, false, natural4 * 0.5f,
                              OperationFailureReason::PredictedChi2Failure, 0);
}

template <TransitionPolicyTag Tag>
void runNaNMaxChi2IgnoredUntilHit4()
{
  using Traits = LegTraits<Tag>;
  const auto baseState = Traits::baseState();
  const auto baseLinRef = Traits::baseLinRef(baseState);

  std::array<SurfaceMeasurement, 4> measurements{
    Traits::measurement(SurfaceId{0}), Traits::measurement(SurfaceId{1}),
    Traits::measurement(SurfaceId{2}), Traits::measurement(SurfaceId{3})};
  auto surfaces = makeCatalogSurfaces<Tag>(4);
  SurfaceCatalogView catalog{surfaces.data(), static_cast<uint32_t>(surfaces.size())};

  const float nanMaxChi2 = std::numeric_limits<float>::quiet_NaN();

  // Only 3 accepted hits: the gate never activates, so a poisoned NaN
  // maxChi2 is never even read and must be completely harmless.
  gsl::span<const SurfaceMeasurement> firstThree(measurements.data(), 3);
  expectSuccessMatchesOracle<Tag>(baseState, baseLinRef, 0.f, firstThree, catalog, Traits::bz,
                                  material::MaterialTraversalDirection::AlongMomentum, false, nanMaxChi2);

  // All 4: hit 4 activates the gate, and refitHit<Tag>'s own chi2/config
  // hardening rejects the non-finite maxChi2.
  expectFailureRollsBack<Tag>(baseState, baseLinRef, 0.f, gsl::span<const SurfaceMeasurement>(measurements), catalog,
                              Traits::bz, material::MaterialTraversalDirection::AlongMomentum, false, nanMaxChi2,
                              OperationFailureReason::NonFiniteInput, 0);
}

template <TransitionPolicyTag Tag>
void runZeroLengthLegSucceeds()
{
  using Traits = LegTraits<Tag>;
  const auto baseState = Traits::baseState();
  const auto baseLinRef = Traits::baseLinRef(baseState);
  gsl::span<const SurfaceMeasurement> empty;
  SurfaceCatalogView poisonedCatalog{nullptr, 123}; // must never be inspected: no present slot

  auto state = baseState;
  auto linRef = baseLinRef;
  float chi2 = 0.25f;
  uint32_t acceptedHitCount = 42;
  OperationFailureReason reason{};
  BOOST_REQUIRE(driveRefitLeg<Tag>(state, linRef, chi2, acceptedHitCount, empty, poisonedCatalog, Traits::bz,
                                   material::MaterialTraversalDirection::AlongMomentum, false,
                                   std::numeric_limits<float>::quiet_NaN(), reason));
  BOOST_CHECK(bitEqual(state, baseState));
  BOOST_CHECK(bitEqual(linRef, baseLinRef));
  BOOST_CHECK_EQUAL(chi2, 0.25f);
  BOOST_CHECK_EQUAL(acceptedHitCount, 0u);
}

template <TransitionPolicyTag Tag>
void runAllHoleLegSucceeds()
{
  using Traits = LegTraits<Tag>;
  const auto baseState = Traits::baseState();
  const auto baseLinRef = Traits::baseLinRef(baseState);

  std::array<SurfaceMeasurement, 3> holes{};
  holes[0].surface = SurfaceId{9999}; // garbage/out-of-range, must never be inspected
  holes[1].surface = SurfaceId{};     // invalid, must never be inspected
  holes[2].surface = SurfaceId{2};    // in-range-looking, still just a hole
  SurfaceCatalogView poisonedCatalog{nullptr, 5};

  auto state = baseState;
  auto linRef = baseLinRef;
  float chi2 = 0.1f;
  uint32_t acceptedHitCount = 7;
  OperationFailureReason reason{};
  BOOST_REQUIRE(driveRefitLeg<Tag>(state, linRef, chi2, acceptedHitCount, holes, poisonedCatalog, Traits::bz,
                                   material::MaterialTraversalDirection::OppositeMomentum, true,
                                   std::numeric_limits<float>::quiet_NaN(), reason));
  BOOST_CHECK(bitEqual(state, baseState));
  BOOST_CHECK(bitEqual(linRef, baseLinRef));
  BOOST_CHECK_EQUAL(chi2, 0.1f);
  BOOST_CHECK_EQUAL(acceptedHitCount, 0u);
}

template <TransitionPolicyTag Tag>
void runCatalogAssociationRejections()
{
  using Traits = LegTraits<Tag>;
  const auto baseState = Traits::baseState();
  const auto baseLinRef = Traits::baseLinRef(baseState);

  auto goodSurfaces = makeCatalogSurfaces<Tag>(4);
  SurfaceCatalogView goodCatalog{goodSurfaces.data(), static_cast<uint32_t>(goodSurfaces.size())};

  // Invalid measurement.surface (default-constructed SurfaceId).
  {
    auto m = Traits::measurement(SurfaceId{0});
    m.surface = SurfaceId{};
    std::array<SurfaceMeasurement, 1> slots{m};
    expectFailureRollsBack<Tag>(baseState, baseLinRef, 0.f, slots, goodCatalog, Traits::bz,
                                material::MaterialTraversalDirection::AlongMomentum, false, 10.f,
                                OperationFailureReason::InvalidSurfaceCatalogAssociation, 0);
  }

  // Out-of-range measurement.surface.value() (>= nSurfaces).
  {
    std::array<SurfaceMeasurement, 1> slots{Traits::measurement(SurfaceId{10})};
    expectFailureRollsBack<Tag>(baseState, baseLinRef, 0.f, slots, goodCatalog, Traits::bz,
                                material::MaterialTraversalDirection::AlongMomentum, false, 10.f,
                                OperationFailureReason::InvalidSurfaceCatalogAssociation, 0);
  }

  // Resolved descriptor's own id does not match measurement.surface.
  {
    auto mismatchedSurfaces = goodSurfaces;
    mismatchedSurfaces[1].id = SurfaceId{99};
    SurfaceCatalogView mismatchedCatalog{mismatchedSurfaces.data(), static_cast<uint32_t>(mismatchedSurfaces.size())};
    std::array<SurfaceMeasurement, 1> slots{Traits::measurement(SurfaceId{1})};
    expectFailureRollsBack<Tag>(baseState, baseLinRef, 0.f, slots, mismatchedCatalog, Traits::bz,
                                material::MaterialTraversalDirection::AlongMomentum, false, 10.f,
                                OperationFailureReason::InvalidSurfaceCatalogAssociation, 0);
  }

  // nullptr surfaces pointer with a nonzero nSurfaces -- must be rejected
  // before any getSurface() dereference (never a dangling/null access).
  {
    SurfaceCatalogView badCatalog{nullptr, 4};
    std::array<SurfaceMeasurement, 1> slots{Traits::measurement(SurfaceId{0})};
    expectFailureRollsBack<Tag>(baseState, baseLinRef, 0.f, slots, badCatalog, Traits::bz,
                                material::MaterialTraversalDirection::AlongMomentum, false, 10.f,
                                OperationFailureReason::InvalidSurfaceCatalogAssociation, 0);
  }
}

template <TransitionPolicyTag Tag>
void runHolePositions()
{
  using Traits = LegTraits<Tag>;
  const auto baseState = Traits::baseState();
  const auto baseLinRef = Traits::baseLinRef(baseState);
  auto surfaces = makeCatalogSurfaces<Tag>(3);
  SurfaceCatalogView catalog{surfaces.data(), static_cast<uint32_t>(surfaces.size())};

  SurfaceMeasurement hole{};
  hole.surface = SurfaceId{42}; // garbage, must never be inspected (cluster invalid)

  auto m0 = Traits::measurement(SurfaceId{0});
  auto m1 = Traits::measurement(SurfaceId{1});

  {
    std::array<SurfaceMeasurement, 3> slots{hole, m0, m1};
    expectSuccessMatchesOracle<Tag>(baseState, baseLinRef, 0.f, slots, catalog, Traits::bz,
                                    material::MaterialTraversalDirection::AlongMomentum, false, 1.e6f);
  }
  {
    std::array<SurfaceMeasurement, 3> slots{m0, hole, m1};
    expectSuccessMatchesOracle<Tag>(baseState, baseLinRef, 0.f, slots, catalog, Traits::bz,
                                    material::MaterialTraversalDirection::AlongMomentum, false, 1.e6f);
  }
  {
    std::array<SurfaceMeasurement, 3> slots{m0, m1, hole};
    expectSuccessMatchesOracle<Tag>(baseState, baseLinRef, 0.f, slots, catalog, Traits::bz,
                                    material::MaterialTraversalDirection::AlongMomentum, false, 1.e6f);
  }
}

template <TransitionPolicyTag Tag>
void runRollbackPositions()
{
  using Traits = LegTraits<Tag>;
  const auto baseState = Traits::baseState();
  const auto baseLinRef = Traits::baseLinRef(baseState);
  constexpr uint32_t kBadIndex = 2;
  auto surfaces = makeCatalogSurfaces<Tag>(3, kBadIndex);
  SurfaceCatalogView catalog{surfaces.data(), static_cast<uint32_t>(surfaces.size())};

  auto good0 = Traits::measurement(SurfaceId{0});
  auto good1 = Traits::measurement(SurfaceId{1});
  auto bad = Traits::measurement(SurfaceId{kBadIndex});

  {
    std::array<SurfaceMeasurement, 3> slots{bad, good0, good1};
    expectFailureRollsBack<Tag>(baseState, baseLinRef, 0.f, slots, catalog, Traits::bz,
                                material::MaterialTraversalDirection::AlongMomentum, false, 1.e6f,
                                OperationFailureReason::MaterialFailure, 0);
  }
  {
    std::array<SurfaceMeasurement, 3> slots{good0, bad, good1};
    expectFailureRollsBack<Tag>(baseState, baseLinRef, 0.f, slots, catalog, Traits::bz,
                                material::MaterialTraversalDirection::AlongMomentum, false, 1.e6f,
                                OperationFailureReason::MaterialFailure, 0);
  }
  {
    std::array<SurfaceMeasurement, 3> slots{good0, good1, bad};
    expectFailureRollsBack<Tag>(baseState, baseLinRef, 0.f, slots, catalog, Traits::bz,
                                material::MaterialTraversalDirection::AlongMomentum, false, 1.e6f,
                                OperationFailureReason::MaterialFailure, 0);
  }
}

template <TransitionPolicyTag Tag>
void runNonMonotonicOrderingMatchesOracle()
{
  using Traits = LegTraits<Tag>;
  const auto baseState = Traits::baseState();
  const auto baseLinRef = Traits::baseLinRef(baseState);
  auto surfaces = makeCatalogSurfaces<Tag>(3);
  SurfaceCatalogView catalog{surfaces.data(), static_cast<uint32_t>(surfaces.size())};

  std::array<SurfaceMeasurement, 3> slots{
    Traits::measurement(SurfaceId{2}), Traits::measurement(SurfaceId{0}), Traits::measurement(SurfaceId{1})};

  expectSuccessMatchesOracle<Tag>(baseState, baseLinRef, 0.f, slots, catalog, Traits::bz,
                                  material::MaterialTraversalDirection::OppositeMomentum, true, 1.e6f);
}

// A single-hit leg, matching testRefitHit.cxx's own RefitHit*ShiftReferenceOnOff
// scope exactly: shiftReferenceToMeasurement is applied only after that same
// hit's update, so it cannot affect this hit's own fitted state. (A
// multi-hit leg would not have this property: a shifted linRef legitimately
// changes the *next* hit's linRef-aware propagation/transport, so comparing
// final states on/off across more than one hit is not a valid invariant.)
template <TransitionPolicyTag Tag>
void runShiftReferenceOnOff()
{
  using Traits = LegTraits<Tag>;
  const auto baseState = Traits::baseState();
  const auto baseLinRef = Traits::baseLinRef(baseState);
  auto surfaces = makeCatalogSurfaces<Tag>(1);
  SurfaceCatalogView catalog{surfaces.data(), static_cast<uint32_t>(surfaces.size())};
  std::array<SurfaceMeasurement, 1> slots{Traits::measurement(SurfaceId{0})};

  auto onState = baseState;
  auto onRef = baseLinRef;
  float onChi2 = 0.f;
  uint32_t onCount = 0;
  OperationFailureReason reason{};
  BOOST_REQUIRE(driveRefitLeg<Tag>(onState, onRef, onChi2, onCount, slots, catalog, Traits::bz,
                                   material::MaterialTraversalDirection::AlongMomentum, true, 1.e6f, reason));
  BOOST_CHECK_EQUAL(onRef.parameters[0], slots[0].frame.u);
  BOOST_CHECK_EQUAL(onRef.parameters[1], slots[0].frame.v);

  auto offState = baseState;
  auto offRef = baseLinRef;
  float offChi2 = 0.f;
  uint32_t offCount = 0;
  BOOST_REQUIRE(driveRefitLeg<Tag>(offState, offRef, offChi2, offCount, slots, catalog, Traits::bz,
                                   material::MaterialTraversalDirection::AlongMomentum, false, 1.e6f, reason));
  BOOST_CHECK_NE(offRef.parameters[0], slots[0].frame.u);
  BOOST_CHECK(bitEqual(onState, offState)); // shift never touches the fitted state
}

template <TransitionPolicyTag Tag>
void runBothMaterialDirections()
{
  using Traits = LegTraits<Tag>;
  const auto baseState = Traits::baseState();
  const auto baseLinRef = Traits::baseLinRef(baseState);
  auto surfaces = makeCatalogSurfaces<Tag>(1);
  SurfaceCatalogView catalog{surfaces.data(), static_cast<uint32_t>(surfaces.size())};
  std::array<SurfaceMeasurement, 1> slots{Traits::measurement(SurfaceId{0})};

  auto alongState = baseState;
  auto alongRef = baseLinRef;
  float alongChi2 = 0.f;
  uint32_t alongCount = 0;
  OperationFailureReason reason{};
  BOOST_REQUIRE(driveRefitLeg<Tag>(alongState, alongRef, alongChi2, alongCount, slots, catalog, Traits::bz,
                                   material::MaterialTraversalDirection::AlongMomentum, false, 1.e6f, reason));

  auto oppositeState = baseState;
  auto oppositeRef = baseLinRef;
  float oppositeChi2 = 0.f;
  uint32_t oppositeCount = 0;
  BOOST_REQUIRE(driveRefitLeg<Tag>(oppositeState, oppositeRef, oppositeChi2, oppositeCount, slots, catalog, Traits::bz,
                                   material::MaterialTraversalDirection::OppositeMomentum, false, 1.e6f, reason));

  // Same isolation rationale as RefitHitBarrelDirectionAffectsQ2PtOppositely
  // in testRefitHit.cxx: compare the two full results against each other.
  BOOST_CHECK_GT(alongState.parameters[4], oppositeState.parameters[4]);
}

template <TransitionPolicyTag Tag>
void runByteIdenticalOracleReplay()
{
  using Traits = LegTraits<Tag>;
  const auto baseState = Traits::baseState();
  const auto baseLinRef = Traits::baseLinRef(baseState);
  auto surfaces = makeCatalogSurfaces<Tag>(3);
  SurfaceCatalogView catalog{surfaces.data(), static_cast<uint32_t>(surfaces.size())};
  std::array<SurfaceMeasurement, 3> slots{
    Traits::measurement(SurfaceId{0}), Traits::measurement(SurfaceId{1}), Traits::measurement(SurfaceId{2})};

  expectSuccessMatchesOracle<Tag>(baseState, baseLinRef, 0.5f, slots, catalog, Traits::bz,
                                  material::MaterialTraversalDirection::OppositeMomentum, true, 1.e6f);
}

} // namespace

// ===========================================================================
// CylinderCylinder
// ===========================================================================

BOOST_AUTO_TEST_CASE(DriveRefitLegBarrelOutputCountStartsFromZeroDespiteSentinel) { runOutputCountStartsFromZeroDespiteSentinel<TransitionPolicyTag::CylinderCylinder>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegBarrelSentinelUnchangedOnFailure) { runSentinelUnchangedOnFailure<TransitionPolicyTag::CylinderCylinder>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegBarrelHits1To3UngatedHit4Gated) { runHits1To3UngatedHit4Gated<TransitionPolicyTag::CylinderCylinder>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegBarrelNaNMaxChi2IgnoredUntilHit4) { runNaNMaxChi2IgnoredUntilHit4<TransitionPolicyTag::CylinderCylinder>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegBarrelZeroLengthLegSucceeds) { runZeroLengthLegSucceeds<TransitionPolicyTag::CylinderCylinder>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegBarrelAllHoleLegSucceeds) { runAllHoleLegSucceeds<TransitionPolicyTag::CylinderCylinder>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegBarrelCatalogAssociationRejections) { runCatalogAssociationRejections<TransitionPolicyTag::CylinderCylinder>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegBarrelHolePositions) { runHolePositions<TransitionPolicyTag::CylinderCylinder>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegBarrelRollbackPositions) { runRollbackPositions<TransitionPolicyTag::CylinderCylinder>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegBarrelNonMonotonicOrderingMatchesOracle) { runNonMonotonicOrderingMatchesOracle<TransitionPolicyTag::CylinderCylinder>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegBarrelShiftReferenceOnOff) { runShiftReferenceOnOff<TransitionPolicyTag::CylinderCylinder>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegBarrelBothMaterialDirections) { runBothMaterialDirections<TransitionPolicyTag::CylinderCylinder>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegBarrelByteIdenticalOracleReplay) { runByteIdenticalOracleReplay<TransitionPolicyTag::CylinderCylinder>(); }

// ===========================================================================
// DiskDisk
// ===========================================================================

BOOST_AUTO_TEST_CASE(DriveRefitLegDiskOutputCountStartsFromZeroDespiteSentinel) { runOutputCountStartsFromZeroDespiteSentinel<TransitionPolicyTag::DiskDisk>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegDiskSentinelUnchangedOnFailure) { runSentinelUnchangedOnFailure<TransitionPolicyTag::DiskDisk>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegDiskHits1To3UngatedHit4Gated) { runHits1To3UngatedHit4Gated<TransitionPolicyTag::DiskDisk>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegDiskNaNMaxChi2IgnoredUntilHit4) { runNaNMaxChi2IgnoredUntilHit4<TransitionPolicyTag::DiskDisk>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegDiskZeroLengthLegSucceeds) { runZeroLengthLegSucceeds<TransitionPolicyTag::DiskDisk>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegDiskAllHoleLegSucceeds) { runAllHoleLegSucceeds<TransitionPolicyTag::DiskDisk>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegDiskCatalogAssociationRejections) { runCatalogAssociationRejections<TransitionPolicyTag::DiskDisk>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegDiskHolePositions) { runHolePositions<TransitionPolicyTag::DiskDisk>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegDiskRollbackPositions) { runRollbackPositions<TransitionPolicyTag::DiskDisk>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegDiskNonMonotonicOrderingMatchesOracle) { runNonMonotonicOrderingMatchesOracle<TransitionPolicyTag::DiskDisk>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegDiskShiftReferenceOnOff) { runShiftReferenceOnOff<TransitionPolicyTag::DiskDisk>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegDiskBothMaterialDirections) { runBothMaterialDirections<TransitionPolicyTag::DiskDisk>(); }
BOOST_AUTO_TEST_CASE(DriveRefitLegDiskByteIdenticalOracleReplay) { runByteIdenticalOracleReplay<TransitionPolicyTag::DiskDisk>(); }

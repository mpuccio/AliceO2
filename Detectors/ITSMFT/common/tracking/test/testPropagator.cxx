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

#define BOOST_TEST_MODULE ITSMFTPropagator
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>

#include "CommonConstants/MathConstants.h"
#include "ITSMFTTracking/RefitDriver.h"
#include "ITSMFTTracking/detail/SurfaceStateOperations.h"
#include "ITSMFTTracking/Propagator.h"

#if __has_include("ITSMFTTracking/BarrelSurfaceStateOperations.h") || __has_include("ITSMFTTracking/ForwardSurfaceStateOperations.h")
#error "coordinate-family state operations must remain private to Propagator"
#endif

using namespace o2::itsmft::tracking;

namespace
{

template <typename T>
bool bitEqual(const T& lhs, const T& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
}

// --- Barrel fixtures (same convention as testRefitHit.cxx's barrelState()) --

SurfaceTrackState barrelState(uint8_t absCharge = 1, o2::track::PID pid = o2::track::PID::Pion)
{
  SurfaceTrackState state{};
  state.parameters[0] = 1.25f;
  state.parameters[1] = -0.75f;
  state.parameters[2] = 0.2f;
  state.parameters[3] = -0.35f;
  state.parameters[4] = 0.8f;
  state.referenceCoordinate = 4.f;
  state.alpha = 0.3f;
  state.kind = SurfaceKind::Cylinder;
  state.absCharge = absCharge;
  state.pid = pid;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.01f * (row + 1) : 0.0002f * (row + column + 1);
    }
  }
  return state;
}

SurfaceTrackParameters barrelLinRef(const SurfaceTrackState& state)
{
  return SurfaceTrackParameters{state};
}

SurfaceMeasurement barrelMeasurement()
{
  SurfaceMeasurement measurement{};
  measurement.frame.q = 2.5f;
  measurement.frame.frameAngle = 0.3f; // same alpha as barrelState(): no rotation needed
  measurement.frame.u = 0.8f;
  measurement.frame.v = -0.45f;
  measurement.covariance = {0.04f, 0.012f, 0.09f};
  return measurement;
}

constexpr float BarrelBz = 5.f;

SurfaceDescriptor cylinderDescriptor(NominalSurfaceMaterial material)
{
  SurfaceDescriptor descriptor{};
  descriptor.kind = SurfaceKind::Cylinder;
  descriptor.referenceCoordinate = 2.5f;
  descriptor.material = material;
  return descriptor;
}

// --- Disk fixtures (same convention as testRefitHit.cxx's diskState()) -----

SurfaceTrackState diskState(uint8_t absCharge = 1, o2::track::PID pid = o2::track::PID::Pion)
{
  SurfaceTrackState state{};
  state.parameters[0] = 1.25f;
  state.parameters[1] = -0.75f;
  state.parameters[2] = 0.35f;
  state.parameters[3] = -2.5f;
  state.parameters[4] = 0.8f;
  state.referenceCoordinate = -45.f;
  state.kind = SurfaceKind::Disk;
  state.absCharge = absCharge;
  state.pid = pid;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.01f * (row + 1) : 0.0002f * (row + column + 1);
    }
  }
  return state;
}

SurfaceTrackParameters diskLinRef(const SurfaceTrackState& state)
{
  return SurfaceTrackParameters{state};
}

SurfaceMeasurement diskMeasurement()
{
  SurfaceMeasurement measurement{};
  measurement.frame = {-50.f, 0.8f, -0.45f, 0.f};
  measurement.frame.q = -50.f;
  measurement.frame.u = 0.8f;
  measurement.frame.v = -0.45f;
  measurement.covariance = {0.04f, 0.f, 0.09f};
  return measurement;
}

constexpr float DiskBz = 5.f;

SurfaceDescriptor diskDescriptor(NominalSurfaceMaterial material)
{
  SurfaceDescriptor descriptor{};
  descriptor.kind = SurfaceKind::Disk;
  descriptor.referenceCoordinate = -50.f;
  descriptor.material = material;
  return descriptor;
}

} // namespace

// --- 1/2: same-family propagate-to-measurement succeeds ---------------------

BOOST_AUTO_TEST_CASE(CylinderToCylinderPropagateAndUpdateSucceeds)
{
  auto state = barrelState();
  auto linRef = barrelLinRef(state);
  const auto measurement = barrelMeasurement();
  const auto descriptor = cylinderDescriptor(NominalSurfaceMaterial{0.f, 0.f});
  float chi2 = 0.f;
  OperationFailureReason reason{};

  BOOST_REQUIRE(Propagator::propagateToMeasurement(state, linRef, descriptor, measurement, BarrelBz,
                                                   material::MaterialTraversalDirection::AlongMomentum,
                                                   false, 0.f, chi2, false, reason));
  BOOST_CHECK_EQUAL(static_cast<int>(state.kind), static_cast<int>(SurfaceKind::Cylinder));
  BOOST_CHECK_EQUAL(state.referenceCoordinate, measurement.frame.q);
  BOOST_CHECK(std::isfinite(chi2));
  BOOST_CHECK_GE(chi2, 0.f);
}

BOOST_AUTO_TEST_CASE(DiskToDiskPropagateAndUpdateSucceeds)
{
  auto state = diskState();
  auto linRef = diskLinRef(state);
  const auto measurement = diskMeasurement();
  const auto descriptor = diskDescriptor(NominalSurfaceMaterial{0.f, 0.f});
  float chi2 = 0.f;
  OperationFailureReason reason{};

  BOOST_REQUIRE(Propagator::propagateToMeasurement(state, linRef, descriptor, measurement, DiskBz,
                                                   material::MaterialTraversalDirection::AlongMomentum,
                                                   false, 0.f, chi2, false, reason));
  BOOST_CHECK_EQUAL(static_cast<int>(state.kind), static_cast<int>(SurfaceKind::Disk));
  BOOST_CHECK_EQUAL(state.referenceCoordinate, measurement.frame.q);
  BOOST_CHECK(std::isfinite(chi2));
  BOOST_CHECK_GE(chi2, 0.f);
}

BOOST_AUTO_TEST_CASE(AcceptedForwardPropagationSelectsFieldAndLowFieldPaths)
{
  auto fieldOn = diskState();
  auto lowPositive = diskState();
  auto lowNegative = diskState();
  OperationFailureReason reason{};

  BOOST_REQUIRE(Propagator::propagateToReference(fieldOn, -50.f, 5.f, reason));
  BOOST_REQUIRE(Propagator::propagateToReference(lowPositive, -50.f, 0.01f, reason));
  BOOST_REQUIRE(Propagator::propagateToReference(lowNegative, -50.f, -0.01f, reason));
  BOOST_CHECK(bitEqual(lowPositive, lowNegative));
  BOOST_CHECK(!bitEqual(fieldOn, lowPositive));
}

BOOST_AUTO_TEST_CASE(PropagatorSelectsCompatibilityFromStateKind)
{
  auto cylinderReference = barrelState();
  auto cylinderCandidate = cylinderReference;
  auto diskReference = diskState();
  auto diskCandidate = diskReference;
  float chi2 = -1.f;
  OperationFailureReason reason{};

  BOOST_REQUIRE(Propagator::stateChi2(cylinderReference, cylinderCandidate, chi2, reason));
  BOOST_CHECK_EQUAL(chi2, 0.f);
  BOOST_REQUIRE(Propagator::stateChi2(diskReference, diskCandidate, chi2, reason));
  BOOST_CHECK_EQUAL(chi2, 0.f);
  BOOST_CHECK(!Propagator::stateChi2(cylinderReference, diskCandidate, chi2, reason));
  BOOST_CHECK(reason == OperationFailureReason::SourceSurfaceKindMismatch);
}

// --- 3: compatible family never converts -- exact agreement with a direct
// detail::barrel::rotate/propagate/correctForMaterial/predictedChi2/update replay ---

BOOST_AUTO_TEST_CASE(CompatibleFamilyMatchesDirectBarrelPrimitiveReplay)
{
  auto viaPropagator = barrelState();
  auto viaPropagatorRef = barrelLinRef(viaPropagator);
  auto viaDirect = viaPropagator;
  auto viaDirectRef = viaPropagatorRef;
  const auto measurement = barrelMeasurement();
  const auto material = NominalSurfaceMaterial{0.01f, 0.001f};
  const auto descriptor = cylinderDescriptor(material);
  float chi2Propagator = 0.f;
  float chi2Direct = 0.f;
  OperationFailureReason reason{};

  BOOST_REQUIRE(Propagator::propagateToMeasurement(viaPropagator, viaPropagatorRef, descriptor, measurement, BarrelBz,
                                                   material::MaterialTraversalDirection::OppositeMomentum,
                                                   false, 0.f, chi2Propagator, true, reason));

  BOOST_REQUIRE(detail::barrel::rotate(viaDirect, viaDirectRef, measurement.frame.frameAngle, BarrelBz, reason));
  BOOST_REQUIRE(detail::barrel::propagate(viaDirect, viaDirectRef, measurement.frame.q, BarrelBz, reason));
  const auto materialResult = detail::barrel::correctForMaterial(
    viaDirect, material::IntegratedMaterialBudget{material.xOverX0, material.arealDensityGPerCm2},
    material::MaterialTraversalDirection::OppositeMomentum);
  BOOST_REQUIRE(materialResult.ok());
  float predChi2 = 0.f;
  BOOST_REQUIRE(detail::barrel::predictedChi2(viaDirect, measurement, predChi2, reason));
  float updateChi2 = 0.f;
  BOOST_REQUIRE(detail::barrel::update(viaDirect, measurement, updateChi2, reason));
  chi2Direct = updateChi2;
  BOOST_REQUIRE(detail::barrel::shiftReferenceToMeasurement(viaDirectRef, measurement, reason));

  BOOST_CHECK(bitEqual(viaPropagator, viaDirect));
  BOOST_CHECK(bitEqual(viaPropagatorRef, viaDirectRef));
  BOOST_CHECK_EQUAL(chi2Propagator, chi2Direct);
}

// --- 4: incompatible family converts, then propagates -----------------------

BOOST_AUTO_TEST_CASE(BarrelStateConvertsToForwardThenPropagatesToDiskMeasurement)
{
  auto state = barrelState();
  auto linRef = barrelLinRef(state);
  const auto poisonState = state;

  // A disk far enough along z that the converted (Forward) state can reach it.
  SurfaceMeasurement measurement{};
  measurement.frame.q = -10.f;
  measurement.frame.u = 5.f;
  measurement.frame.v = -5.f;
  measurement.covariance = {10.f, 0.f, 10.f}; // loose: the point is not expected to land exactly here
  const auto descriptor = diskDescriptor(NominalSurfaceMaterial{0.f, 0.f});
  float chi2 = 0.f;
  OperationFailureReason reason{};

  const bool ok = Propagator::propagateToMeasurement(state, linRef, descriptor, measurement, BarrelBz,
                                                     material::MaterialTraversalDirection::AlongMomentum,
                                                     false, 0.f, chi2, false, reason);
  BOOST_REQUIRE(ok);
  BOOST_CHECK_EQUAL(static_cast<int>(state.kind), static_cast<int>(SurfaceKind::Disk));
  BOOST_CHECK_EQUAL(state.referenceCoordinate, measurement.frame.q);
  BOOST_CHECK_EQUAL(state.absCharge, poisonState.absCharge);
  BOOST_CHECK(state.pid == poisonState.pid);
  for (float value : state.parameters) {
    BOOST_CHECK(std::isfinite(value));
  }
  for (float value : state.covariance) {
    BOOST_CHECK(std::isfinite(value));
  }
}

BOOST_AUTO_TEST_CASE(KindConversionRelinearizesAtConvertedState)
{
  auto nominalState = barrelState();
  auto nominalRef = barrelLinRef(nominalState);
  auto perturbedState = nominalState;
  auto perturbedRef = nominalRef;
  perturbedRef.parameters[0] += 0.1f;
  perturbedRef.parameters[1] -= 0.2f;
  perturbedRef.parameters[2] += 0.01f;
  perturbedRef.parameters[3] -= 0.02f;
  perturbedRef.parameters[4] += 0.001f;

  SurfaceMeasurement measurement{};
  measurement.frame.q = -10.f;
  measurement.frame.u = 5.f;
  measurement.frame.v = -5.f;
  measurement.covariance = {10.f, 0.f, 10.f};
  const auto descriptor = diskDescriptor(NominalSurfaceMaterial{0.f, 0.f});
  float nominalChi2 = 0.f;
  float perturbedChi2 = 0.f;
  OperationFailureReason reason{};

  BOOST_REQUIRE(Propagator::propagateToMeasurement(nominalState, nominalRef, descriptor, measurement, BarrelBz,
                                                   material::MaterialTraversalDirection::AlongMomentum,
                                                   false, 0.f, nominalChi2, false, reason));
  BOOST_REQUIRE(Propagator::propagateToMeasurement(perturbedState, perturbedRef, descriptor, measurement, BarrelBz,
                                                   material::MaterialTraversalDirection::AlongMomentum,
                                                   false, 0.f, perturbedChi2, false, reason));

  BOOST_CHECK(bitEqual(perturbedState, nominalState));
  BOOST_CHECK(bitEqual(perturbedRef, nominalRef));
  BOOST_CHECK_EQUAL(perturbedChi2, nominalChi2);
}

BOOST_AUTO_TEST_CASE(ReverseKindConversionRelinearizesAtConvertedState)
{
  auto nominalState = diskState();
  auto nominalRef = diskLinRef(nominalState);
  auto perturbedState = nominalState;
  auto perturbedRef = nominalRef;
  perturbedRef.parameters[0] += 0.1f;
  perturbedRef.parameters[1] -= 0.2f;
  perturbedRef.parameters[2] += 0.01f;
  perturbedRef.parameters[3] -= 0.02f;
  perturbedRef.parameters[4] += 0.001f;

  const auto measurement = barrelMeasurement();
  const auto descriptor = cylinderDescriptor(NominalSurfaceMaterial{0.f, 0.f});
  float nominalChi2 = 0.f;
  float perturbedChi2 = 0.f;
  OperationFailureReason reason{};

  BOOST_REQUIRE(Propagator::propagateToMeasurement(nominalState, nominalRef, descriptor, measurement, DiskBz,
                                                   material::MaterialTraversalDirection::AlongMomentum,
                                                   false, 0.f, nominalChi2, false, reason));
  BOOST_REQUIRE(Propagator::propagateToMeasurement(perturbedState, perturbedRef, descriptor, measurement, DiskBz,
                                                   material::MaterialTraversalDirection::AlongMomentum,
                                                   false, 0.f, perturbedChi2, false, reason));

  BOOST_CHECK(bitEqual(perturbedState, nominalState));
  BOOST_CHECK(bitEqual(perturbedRef, nominalRef));
  BOOST_CHECK_EQUAL(perturbedChi2, nominalChi2);
}

BOOST_AUTO_TEST_CASE(ConvertFamilyPreservesChargeAndPID)
{
  auto state = barrelState(2, o2::track::PID::Kaon);
  OperationFailureReason reason{};
  BOOST_REQUIRE(Propagator::convertKind(state, SurfaceKind::Disk, reason));
  BOOST_CHECK_EQUAL(static_cast<int>(state.kind), static_cast<int>(SurfaceKind::Disk));
  BOOST_CHECK_EQUAL(state.absCharge, uint8_t{2});
  BOOST_CHECK(state.pid == o2::track::PID::Kaon);
}

BOOST_AUTO_TEST_CASE(ConvertFamilySameFamilyIsNoOpSuccess)
{
  auto state = barrelState();
  const auto before = state;
  OperationFailureReason reason{};
  BOOST_REQUIRE(Propagator::convertKind(state, SurfaceKind::Cylinder, reason));
  BOOST_CHECK(bitEqual(state, before));
}

// --- 5: degenerate conversion fails, transactionally ------------------------

BOOST_AUTO_TEST_CASE(ForwardToBarrelConversionFailsAtOriginTransactionally)
{
  auto state = diskState();
  state.parameters[0] = 0.f; // X
  state.parameters[1] = 0.f; // Y: R == 0, alpha undefined
  const auto poison = state;
  OperationFailureReason reason{};

  BOOST_CHECK(!Propagator::convertKind(state, SurfaceKind::Cylinder, reason));
  BOOST_CHECK_EQUAL(static_cast<int>(reason), static_cast<int>(OperationFailureReason::SurfaceKindConversionFailure));
  BOOST_CHECK(bitEqual(state, poison));
}

// --- Zero-material and nonzero-material (MatLUT/nominal-material) paths -----

BOOST_AUTO_TEST_CASE(ZeroMaterialPathSucceeds)
{
  auto state = barrelState();
  auto linRef = barrelLinRef(state);
  const auto measurement = barrelMeasurement();
  const auto descriptor = cylinderDescriptor(NominalSurfaceMaterial{0.f, 0.f});
  float chi2 = 0.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(Propagator::propagateToMeasurement(state, linRef, descriptor, measurement, BarrelBz,
                                                   material::MaterialTraversalDirection::AlongMomentum,
                                                   false, 0.f, chi2, false, reason));
}

BOOST_AUTO_TEST_CASE(NonzeroNominalMaterialChangesResultRelativeToZeroMaterial)
{
  auto zeroState = barrelState();
  auto zeroRef = barrelLinRef(zeroState);
  auto materialState = barrelState();
  auto materialRef = barrelLinRef(materialState);
  const auto measurement = barrelMeasurement();
  const auto zeroDescriptor = cylinderDescriptor(NominalSurfaceMaterial{0.f, 0.f});
  const auto materialDescriptor = cylinderDescriptor(NominalSurfaceMaterial{0.05f, 0.01f});
  float zeroChi2 = 0.f;
  float materialChi2 = 0.f;
  OperationFailureReason reason{};

  BOOST_REQUIRE(Propagator::propagateToMeasurement(zeroState, zeroRef, zeroDescriptor, measurement, BarrelBz,
                                                   material::MaterialTraversalDirection::OppositeMomentum,
                                                   false, 0.f, zeroChi2, false, reason));
  BOOST_REQUIRE(Propagator::propagateToMeasurement(materialState, materialRef, materialDescriptor, measurement, BarrelBz,
                                                   material::MaterialTraversalDirection::OppositeMomentum,
                                                   false, 0.f, materialChi2, false, reason));

  // The material budget is read from the target SurfaceDescriptor (the
  // "MatLUT" mechanism, task requirement 6) -- not equal, not a parallel
  // model producing a byte-identical result either.
  BOOST_CHECK(!bitEqual(zeroState, materialState));
}

// --- Holes are skipped by the native refit driver ----------------------------

BOOST_AUTO_TEST_CASE(RefitDriverSkipsHoleSlots)
{
  auto state = barrelState();
  auto linRef = barrelLinRef(state);
  const auto measurement = barrelMeasurement();

  std::array<SurfaceDescriptor, 1> surfaces{cylinderDescriptor(NominalSurfaceMaterial{0.f, 0.f})};
  SurfaceCatalogView catalog{surfaces.data(), static_cast<uint32_t>(surfaces.size())};

  const detail::RefitMeasurementSlot present{measurement, LayerId{0}, true};
  const detail::RefitMeasurementSlot hole{};

  std::array<detail::RefitMeasurementSlot, 3> slots{hole, present, hole};
  float chi2 = 0.f;
  uint32_t acceptedHitCount = 999;
  OperationFailureReason reason{};

  BOOST_REQUIRE(detail::driveRefitLeg(state, linRef, chi2, acceptedHitCount, slots, catalog, BarrelBz,
                                      material::MaterialTraversalDirection::AlongMomentum, false, 100.f, reason));
  BOOST_CHECK_EQUAL(acceptedHitCount, 1u);
}

// --- 10/11: chi2-gate failure and atomicity ----------------------------------

BOOST_AUTO_TEST_CASE(Chi2GateRejectsOversizedPredictedChi2Transactionally)
{
  auto state = barrelState();
  auto linRef = barrelLinRef(state);
  const auto poisonState = state;
  const auto poisonRef = linRef;
  auto measurement = barrelMeasurement();
  measurement.frame.u += 5.f; // far outlier vs the state's predicted local Y
  const auto descriptor = cylinderDescriptor(NominalSurfaceMaterial{0.f, 0.f});
  float chi2 = 0.f;
  const float poisonChi2 = chi2;
  OperationFailureReason reason{};

  const bool ok = Propagator::propagateToMeasurement(state, linRef, descriptor, measurement, BarrelBz,
                                                     material::MaterialTraversalDirection::AlongMomentum,
                                                     true, 1.e-6f, chi2, false, reason);
  BOOST_CHECK(!ok);
  BOOST_CHECK_EQUAL(static_cast<int>(reason), static_cast<int>(OperationFailureReason::PredictedChi2Failure));
  BOOST_CHECK(bitEqual(state, poisonState));
  BOOST_CHECK(bitEqual(linRef, poisonRef));
  BOOST_CHECK_EQUAL(chi2, poisonChi2);
}

BOOST_AUTO_TEST_CASE(UnrecognizedTargetSurfaceKindFails)
{
  auto state = barrelState();
  auto linRef = barrelLinRef(state);
  const auto poisonState = state;
  const auto measurement = barrelMeasurement();
  SurfaceDescriptor descriptor = cylinderDescriptor(NominalSurfaceMaterial{0.f, 0.f});
  // SurfaceKind currently only has Cylinder/Disk (both recognized); this
  // proves the routing guard itself, not a reachable production input.
  descriptor.kind = static_cast<SurfaceKind>(0xFFu);
  float chi2 = 0.f;
  OperationFailureReason reason{};

  const bool ok = Propagator::propagateToMeasurement(state, linRef, descriptor, measurement, BarrelBz,
                                                     material::MaterialTraversalDirection::AlongMomentum,
                                                     false, 0.f, chi2, false, reason);
  BOOST_CHECK(!ok);
  BOOST_CHECK_EQUAL(static_cast<int>(reason), static_cast<int>(OperationFailureReason::SurfaceKindConversionFailure));
  BOOST_CHECK(bitEqual(state, poisonState));
}

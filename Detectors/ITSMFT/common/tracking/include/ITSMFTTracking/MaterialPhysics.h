// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_MATERIALPHYSICS_H_
#define ALICEO2_ITSMFT_TRACKING_MATERIALPHYSICS_H_

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include "ReconstructionDataFormats/PID.h"

// This header and its implementation are host-only. Nothing here has been
// exercised on GPU, and OPENCL_HOST-style host-side compilation of GPU
// translation units is not evidence of CUDA/HIP support.

namespace o2::itsmft::tracking::material
{

// Semantic direction of material traversal relative to the particle's
// momentum. This is independent of the sign of any propagation/covariance
// convention used by a caller's track state.
enum class MaterialTraversalDirection : uint8_t {
  AlongMomentum = 0,
  OppositeMomentum = 1
};

// Unsigned, already path-integrated material budget crossed by a state.
// Both fields are non-negative by construction; direction is carried
// separately by MaterialTraversalDirection.
struct IntegratedMaterialBudget {
  float xOverX0;             ///< thickness in units of radiation length
  float arealDensityGPerCm2; ///< crossed length*density, g/cm^2
};

// Reasons a scalar material-physics operation can fail to produce a usable
// result. SourceFamilyMismatch, NonFiniteState, InvalidStateKinematics, and
// InvalidCovariance are reserved for later composite family operations that
// carry a full track state; calculateMaterialPhysics() never emits them.
enum class MaterialFailureReason : uint8_t {
  None = 0,
  SourceFamilyMismatch = 1,
  NonFiniteState = 2,
  InvalidStateKinematics = 3,
  InvalidPID = 4,
  ChargedMasslessPID = 5,
  InvalidDirection = 6,
  InvalidMaterial = 7,
  StoppedInMaterial = 8,
  MomentumBelowMinimum = 9,
  ExcessiveScattering = 10,
  InvalidCovariance = 11,
  NonFiniteResult = 12
};

enum class MaterialOperationFlags : uint8_t {
  None = 0,
  SubstepCountClamped = 1
};

// Result of one scalar material-physics evaluation.
//
// Diagnostics on failure: momentumBeforeGeV always echoes the raw
// momentumGeV argument as supplied by the caller (even if that argument
// itself was the cause of failure, e.g. non-finite). failure identifies the
// reason. Every other field (momentumAfterGeV, signedEnergyChangeGeV,
// highlandTheta2Rad2, relativeInverseMomentumVariance,
// energyLossSubsteps, flags) is set to its deterministic zero/None value
// and carries no physical meaning on failure. reserved is always zero.
struct MaterialOperationResult {
  float momentumBeforeGeV;
  float momentumAfterGeV;
  float signedEnergyChangeGeV;
  float highlandTheta2Rad2;
  float relativeInverseMomentumVariance;
  uint8_t energyLossSubsteps;
  MaterialOperationFlags flags;
  MaterialFailureReason failure;
  uint8_t reserved;

  bool ok() const noexcept { return failure == MaterialFailureReason::None; }
};

// Representation properties asserted here describe the current in-memory
// layout of these types. They are not yet a durable serialized or device ABI
// commitment; a future slice may need to revisit them explicitly before
// making such a commitment.
#define O2_ITSMFT_MATERIAL_ASSERT_LAYOUT(Type, Size, Alignment) \
  static_assert(std::is_standard_layout_v<Type>);               \
  static_assert(std::is_trivially_copyable_v<Type>);            \
  static_assert(sizeof(Type) == Size);                          \
  static_assert(alignof(Type) == Alignment)

O2_ITSMFT_MATERIAL_ASSERT_LAYOUT(IntegratedMaterialBudget, 8, 4);
O2_ITSMFT_MATERIAL_ASSERT_LAYOUT(MaterialOperationResult, 24, 4);

#undef O2_ITSMFT_MATERIAL_ASSERT_LAYOUT

static_assert(offsetof(MaterialOperationResult, momentumBeforeGeV) == 0);
static_assert(offsetof(MaterialOperationResult, momentumAfterGeV) == 4);
static_assert(offsetof(MaterialOperationResult, signedEnergyChangeGeV) == 8);
static_assert(offsetof(MaterialOperationResult, highlandTheta2Rad2) == 12);
static_assert(offsetof(MaterialOperationResult, relativeInverseMomentumVariance) == 16);
static_assert(offsetof(MaterialOperationResult, energyLossSubsteps) == 20);
static_assert(offsetof(MaterialOperationResult, flags) == 21);
static_assert(offsetof(MaterialOperationResult, failure) == 22);
static_assert(offsetof(MaterialOperationResult, reserved) == 23);

// Detector-neutral, PID/absCharge-aware scalar material-physics kernel.
//
// pid supplies the mass; absCharge supplies the physical |q| used for the
// q^2 scaling of energy loss and multiple scattering. There is no required
// equality between absCharge and PID::getCharge() -- callers may model
// exotic or non-standard charge states by supplying absCharge independently
// of the nominal PID charge. For absCharge == 0 the state is treated as
// neutral: after the common validation below, the operation always
// succeeds with an unchanged momentum and zero energy-loss/scattering
// contributions, even for nonzero material and even for a valid massive
// neutral PID.
//
// Validation precedence (first failing check wins):
//   1. direction is a recognized MaterialTraversalDirection value, else
//      InvalidDirection.
//   2. both material fields are finite and non-negative, else
//      InvalidMaterial.
//   3. momentumGeV is finite and strictly positive, else
//      MomentumBelowMinimum.
//   4. pid.getID() is below PID::NIDsTot, else InvalidPID (checked before
//      any mass-array access).
//   5. absCharge != 0 with zero PID mass, else ChargedMasslessPID.
//
// For a charged, massive state, the entry-derived total energy (e0) and
// beta^2 are explicitly checked to be finite (and beta^2 additionally
// checked to be strictly positive) before any energy-loss or multiple-
// scattering arithmetic runs, since both quantities feed both branches;
// this is a deterministic NonFiniteResult, not a reliance on later
// inf-inf/inf-over-inf/NaN propagation to eventually surface the problem.
//
// For a charged, massive state, momentumGeV is the physical input momentum
// already selected by the caller (no covariance projection is performed by
// this kernel). Energy loss is evaluated with the same capped-substep
// Bethe-Bloch algorithm as
// o2::track::TrackParametrizationWithError::correctForMaterial(): the
// requested substep count is
//   1 + floor(|dE_full| / eKin * o2::track::ELoss2EKinThreshInv)
// (computed in a manner that never converts a non-finite or out-of-range
// float to int), capped at o2::track::MaxELossIter (50); flags carries
// SubstepCountClamped iff the requested count exceeds 50. Regardless of
// clamping, the complete arealDensityGPerCm2 is processed (the capped
// substep count only changes the granularity, not the total crossed
// material), and Bethe-Bloch (o2::track::BetheBlochSolidOpt<float>()) is
// recomputed from the current momentum at every substep.
// MaterialTraversalDirection::AlongMomentum subtracts energy per substep;
// OppositeMomentum adds it; signedEnergyChangeGeV is always
// Eafter - Ebefore. A particle whose energy would fall to or below its rest
// mass fails with StoppedInMaterial; a particle that completes with
// momentum below 0.01 GeV/c fails with MomentumBelowMinimum.
//
// highlandTheta2Rad2 and relativeInverseMomentumVariance use the
// simplified O2 Highland variance (no logarithmic correction) and the
// pre-material momentum, energy, and beta -- never the post-energy-loss
// values. Both scale with absCharge^2 (relativeInverseMomentumVariance via
// the already-charge-scaled total energy change). highlandTheta2Rad2 > pi^2
// fails with ExcessiveScattering.
//
// This kernel does not include or construct TrackParCovF, TrackParCovFwd,
// TrackFwd, SurfaceKinematicState, detector geometry, or any
// ITS/MFT/topology/policy/Propagator header.
MaterialOperationResult calculateMaterialPhysics(
  float momentumGeV,
  o2::track::PID pid,
  uint8_t absCharge,
  MaterialTraversalDirection direction,
  IntegratedMaterialBudget material) noexcept;

} // namespace o2::itsmft::tracking::material

#endif // ALICEO2_ITSMFT_TRACKING_MATERIALPHYSICS_H_

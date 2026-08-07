// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEKINEMATICSTATE_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEKINEMATICSTATE_H_

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>

#include "GPUCommonDef.h"
#include "ITSMFTTracking/StateFamily.h"
#include "ReconstructionDataFormats/PID.h"

namespace o2::itsmft::tracking
{

// The interpretation of parameters and covariance is selected by a typed view:
// Barrel:  (Y, Z, Snp, Tgl, Q2Pt), referenceCoordinate is local X, alpha is frame angle.
// Forward: (X, Y, Phi, Tanl, InvQPt), referenceCoordinate is global Z, alpha is unused (zero).
struct SurfaceKinematicState {
  float parameters[5]{};
  float covariance[15]{};
  float referenceCoordinate{0.f};
  float alpha{0.f};
  StateFamily family{StateFamily::Invalid};
  uint8_t flags{0}; // Reserved: no flag semantics are defined in this slice.
  uint8_t absCharge{0};
  o2::track::PID pid{o2::track::PID::Pion};

  // This only recognizes the two supported family tags. It does not validate
  // parameter ranges, covariance, PID, alpha, or family-specific metadata.
  GPUhdi() constexpr bool hasRecognizedFamily() const noexcept { return family == StateFamily::Barrel || family == StateFamily::Forward; }
};

static_assert(std::is_standard_layout_v<SurfaceKinematicState>);
static_assert(std::is_trivially_copyable_v<SurfaceKinematicState>);
static_assert(sizeof(SurfaceKinematicState) == 92);
static_assert(alignof(SurfaceKinematicState) == 4);
static_assert(offsetof(SurfaceKinematicState, parameters) == 0);
static_assert(offsetof(SurfaceKinematicState, covariance) == 20);
static_assert(offsetof(SurfaceKinematicState, referenceCoordinate) == 80);
static_assert(offsetof(SurfaceKinematicState, alpha) == 84);
static_assert(offsetof(SurfaceKinematicState, family) == 88);
static_assert(offsetof(SurfaceKinematicState, flags) == 89);
static_assert(offsetof(SurfaceKinematicState, absCharge) == 90);
static_assert(offsetof(SurfaceKinematicState, pid) == 91);

// Transient covariance-free trajectory used as a Kalman linearization point;
// it is transported alongside exactly one fitted SurfaceKinematicState.
//
// Deliberately excludes covariance, chi2, cluster bookkeeping, and
// PID/absCharge: a linearization reference always accompanies exactly one
// SurfaceKinematicState (the fitted state being linearized), which already
// carries the particle hypothesis. Operations that need absCharge/pid
// (e.g. curvature) read them from that paired state, never from a
// duplicated copy here. This also means a SurfaceLinearizationReference
// has no meaning on its own; it is only ever valid in the context of the
// SurfaceKinematicState it was constructed from or is paired with.
//
// Barrel: (Y, Z, Snp, Tgl, Q2Pt), referenceCoordinate is local X, alpha is
// frame angle -- identical parameter/frame convention to SurfaceKinematicState's
// own Barrel interpretation.
// Forward: (X, Y, Phi, Tanl, InvQPt), referenceCoordinate is global Z,
// alpha is unused (zero) -- identical to SurfaceKinematicState's Forward
// interpretation.
//
// The layout assertions describe the current in-memory representation, not a
// serialized ABI.
struct SurfaceLinearizationReference {
  float parameters[5]{};
  float referenceCoordinate{0.f};
  float alpha{0.f};
  StateFamily family{StateFamily::Invalid};

  GPUhdi() constexpr bool hasRecognizedFamily() const noexcept { return family == StateFamily::Barrel || family == StateFamily::Forward; }
};

static_assert(std::is_standard_layout_v<SurfaceLinearizationReference>);
static_assert(std::is_trivially_copyable_v<SurfaceLinearizationReference>);
static_assert(sizeof(SurfaceLinearizationReference) == 32);
static_assert(alignof(SurfaceLinearizationReference) == 4);
static_assert(offsetof(SurfaceLinearizationReference, parameters) == 0);
static_assert(offsetof(SurfaceLinearizationReference, referenceCoordinate) == 20);
static_assert(offsetof(SurfaceLinearizationReference, alpha) == 24);
static_assert(offsetof(SurfaceLinearizationReference, family) == 28);

#ifndef GPUCA_GPUCODE

// Host-only validating factory. Copies parameters, referenceCoordinate, alpha,
// and family from state; rejects unrecognized families or non-finite fields
// and leaves out unchanged on failure.
bool makeLinearizationReference(const SurfaceKinematicState& state, SurfaceLinearizationReference& out) noexcept;

#endif // GPUCA_GPUCODE

GPUhdi() constexpr uint8_t packedCovarianceIndex(uint8_t row, uint8_t column) noexcept
{
  return row >= column ? row * (row + 1) / 2 + column : column * (column + 1) / 2 + row;
}

// Covariance sanitizer used after successful covariance mutations. It takes
// absolute values of diagonals, clamps each to maxDiagonal[i] while rescaling
// its row/column by sqrt(maxDiagonal[i]/diagonal), then clamps each
// off-diagonal to the pairwise Cauchy-Schwarz bound
// sqrt(diagonal_i * diagonal_j). Pairwise bounds are necessary but do not
// guarantee positive semi-definiteness of the full matrix; a 3x3 correlation
// block additionally requires 1 + 2*rho01*rho02*rho12 - rho01^2 - rho02^2 -
// rho12^2 >= 0.
GPUhdi() void sanitizeCovariance(SurfaceKinematicState& state, const float (&maxDiagonal)[5]) noexcept
{
  auto& c = state.covariance;
  for (uint8_t i = 0; i < 5; ++i) {
    const uint8_t diagIndex = packedCovarianceIndex(i, i);
    c[diagIndex] = c[diagIndex] < 0.f ? -c[diagIndex] : c[diagIndex];
    if (c[diagIndex] > maxDiagonal[i]) {
      const float scale = std::sqrt(maxDiagonal[i] / c[diagIndex]);
      c[diagIndex] = maxDiagonal[i];
      for (uint8_t j = 0; j < 5; ++j) {
        if (j == i) {
          continue;
        }
        c[packedCovarianceIndex(i, j)] *= scale;
      }
    }
  }
  for (uint8_t i = 0; i < 5; ++i) {
    for (uint8_t j = 0; j < i; ++j) {
      const float bound = std::sqrt(c[packedCovarianceIndex(i, i)] * c[packedCovarianceIndex(j, j)]);
      const uint8_t offIndex = packedCovarianceIndex(i, j);
      if (c[offIndex] > bound) {
        c[offIndex] = bound;
      } else if (c[offIndex] < -bound) {
        c[offIndex] = -bound;
      }
    }
  }
}

class BarrelStateView;
class ForwardStateView;

// Host validation factories. On failure, `view` is left unchanged.
bool makeBarrelStateView(SurfaceKinematicState& state, BarrelStateView& view) noexcept;
bool makeForwardStateView(SurfaceKinematicState& state, ForwardStateView& view) noexcept;

class BarrelStateView
{
 public:
  // An invalid default view has a null pointer. Dereference only after a host
  // validation factory has returned true.
  BarrelStateView() = default;

  GPUhdi() bool isValid() const noexcept { return mState != nullptr; }

  GPUhdi() float getY() const noexcept { return mState->parameters[0]; }
  GPUhdi() float getZ() const noexcept { return mState->parameters[1]; }
  GPUhdi() float getSnp() const noexcept { return mState->parameters[2]; }
  GPUhdi() float getTgl() const noexcept { return mState->parameters[3]; }
  GPUhdi() float getQ2Pt() const noexcept { return mState->parameters[4]; }
  GPUhdi() float getCovariance(uint8_t row, uint8_t column) const noexcept { return mState->covariance[packedCovarianceIndex(row, column)]; }
  GPUhdi() float getReferenceX() const noexcept { return mState->referenceCoordinate; }
  GPUhdi() float getAlpha() const noexcept { return mState->alpha; }

  GPUhdi() void setY(float value) noexcept { mState->parameters[0] = value; }
  GPUhdi() void setZ(float value) noexcept { mState->parameters[1] = value; }
  GPUhdi() void setSnp(float value) noexcept { mState->parameters[2] = value; }
  GPUhdi() void setTgl(float value) noexcept { mState->parameters[3] = value; }
  GPUhdi() void setQ2Pt(float value) noexcept { mState->parameters[4] = value; }
  GPUhdi() void setCovariance(uint8_t row, uint8_t column, float value) noexcept { mState->covariance[packedCovarianceIndex(row, column)] = value; }

 private:
  explicit BarrelStateView(SurfaceKinematicState* state) noexcept : mState{state} {}

  SurfaceKinematicState* mState{nullptr};

  friend bool makeBarrelStateView(SurfaceKinematicState&, BarrelStateView&) noexcept;
};

class ForwardStateView
{
 public:
  // An invalid default view has a null pointer. Dereference only after a host
  // validation factory has returned true.
  ForwardStateView() = default;

  GPUhdi() bool isValid() const noexcept { return mState != nullptr; }

  GPUhdi() float getX() const noexcept { return mState->parameters[0]; }
  GPUhdi() float getY() const noexcept { return mState->parameters[1]; }
  GPUhdi() float getPhi() const noexcept { return mState->parameters[2]; }
  GPUhdi() float getTanl() const noexcept { return mState->parameters[3]; }
  GPUhdi() float getInvQPt() const noexcept { return mState->parameters[4]; }
  // Forward slot 4 is signed q/pT; getInvQPt()/setInvQPt() are aliases.
  GPUhdi() float getQ2Pt() const noexcept { return mState->parameters[4]; }
  GPUhdi() float getCovariance(uint8_t row, uint8_t column) const noexcept { return mState->covariance[packedCovarianceIndex(row, column)]; }
  GPUhdi() float getReferenceZ() const noexcept { return mState->referenceCoordinate; }

  GPUhdi() void setX(float value) noexcept { mState->parameters[0] = value; }
  GPUhdi() void setY(float value) noexcept { mState->parameters[1] = value; }
  GPUhdi() void setPhi(float value) noexcept { mState->parameters[2] = value; }
  GPUhdi() void setTanl(float value) noexcept { mState->parameters[3] = value; }
  GPUhdi() void setInvQPt(float value) noexcept { mState->parameters[4] = value; }
  GPUhdi() void setQ2Pt(float value) noexcept { mState->parameters[4] = value; }
  GPUhdi() void setCovariance(uint8_t row, uint8_t column, float value) noexcept { mState->covariance[packedCovarianceIndex(row, column)] = value; }

 private:
  explicit ForwardStateView(SurfaceKinematicState* state) noexcept : mState{state} {}

  SurfaceKinematicState* mState{nullptr};

  friend bool makeForwardStateView(SurfaceKinematicState&, ForwardStateView&) noexcept;
};

inline bool makeBarrelStateView(SurfaceKinematicState& state, BarrelStateView& view) noexcept
{
  if (state.family != StateFamily::Barrel) {
    return false;
  }
  view = BarrelStateView{&state};
  return true;
}

inline bool makeForwardStateView(SurfaceKinematicState& state, ForwardStateView& view) noexcept
{
  if (state.family != StateFamily::Forward) {
    return false;
  }
  view = ForwardStateView{&state};
  return true;
}

static_assert(std::is_standard_layout_v<BarrelStateView>);
static_assert(std::is_trivially_copyable_v<BarrelStateView>);
static_assert(sizeof(BarrelStateView) == sizeof(SurfaceKinematicState*));
static_assert(alignof(BarrelStateView) == alignof(SurfaceKinematicState*));
static_assert(std::is_standard_layout_v<ForwardStateView>);
static_assert(std::is_trivially_copyable_v<ForwardStateView>);
static_assert(sizeof(ForwardStateView) == sizeof(SurfaceKinematicState*));
static_assert(alignof(ForwardStateView) == alignof(SurfaceKinematicState*));

} // namespace o2::itsmft::tracking

#endif // ALICEO2_ITSMFT_TRACKING_SURFACEKINEMATICSTATE_H_

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

GPUhdi() constexpr uint8_t packedCovarianceIndex(uint8_t row, uint8_t column) noexcept
{
  return row >= column ? row * (row + 1) / 2 + column : column * (column + 1) / 2 + row;
}

// Detector-neutral covariance-validity invariant, applied unconditionally
// after every successful covariance-mutating state operation (propagate,
// rotate, measurement update -- barrel and forward alike). It knows nothing
// about materials, detectors, or families beyond the five per-slot upper
// bounds the caller supplies. Operates directly on the packed
// lower-triangular storage, so symmetry is preserved by construction; it
// never constructs or depends on any legacy track-parametrization type.
//
// DECLARED INVARIANT (audited against o2::track::TrackParametrizationWithError
// <float>::checkCovariance(), TrackParametrizationWithError.cxx -- read in
// full, not excerpted): every diagonal is non-negative, no diagonal exceeds
// maxDiagonal[i], and no individual off-diagonal pair's Pearson correlation
// exceeds 1 in magnitude. This is the same three-part invariant legacy's own
// checkCovariance() doc comment claims ("forces the diagonal elements... to
// be positive and abs of correlation coefficients to be <1. In case the
// diagonal element is bigger than the maximal allowed value, it is set to
// the limit...").
//
// PASS 1 (diagonal abs + range clamp): byte-faithful to legacy's actual
// code. For each diagonal, take its absolute value, and if it still exceeds
// maxDiagonal[i], clamp it to that bound and rescale every off-diagonal
// entry sharing that row/column by sqrt(maxDiagonal[i]/diagonal). Diagonals
// are processed in slot order (0..4), matching legacy's own sequential
// order, so cumulative rescaling of an entry shared by two out-of-range
// diagonals is deterministic.
//
// PASS 2 (pairwise correlation clamp): legacy's own comment claims this
// ("abs of correlation coefficients to be <1") but legacy's own CODE does
// NOT implement it as a general check -- only pass 1 exists there, so a
// pair whose correlation exceeds 1 while BOTH diagonals remain within range
// is never touched by legacy's real implementation (confirmed by reading
// the complete function body: there is no separate correlation-coefficient
// pass at all). This is exactly the failure mode ADR 0008's covariance-
// fault-localization investigation traced: a large-step propagate can leave
// |correlation(Y,Q2Pt)| > 1 (Cauchy-Schwarz violated) while c(Y,Y) and
// c(Q2Pt,Q2Pt) both stay well inside their maxDiagonal ceilings, so pass 1
// alone never repairs it. Pass 2 completes legacy's own DOCUMENTED (not
// coded) intent using the value that intent already names -- 1, the
// mathematical Cauchy-Schwarz bound every valid covariance matrix must
// satisfy, not a chosen/invented threshold -- by clamping each off-diagonal
// entry's magnitude to sqrt(diagonal_i * diagonal_j) (using the
// already-pass-1-clamped diagonals), independently per pair, in a single
// non-iterative sweep.
//
// EXPLICIT LIMITATION, not silently omitted: pairwise correlation bounds
// are necessary but NOT sufficient for the matrix to be positive
// semi-definite (a joint/multivariate property of 3 or more parameters
// simultaneously, e.g. a valid 3x3 correlation submatrix additionally
// requires 1 + 2*rho01*rho02*rho12 - rho01^2 - rho02^2 - rho12^2 >= 0).
// Verified empirically against the real captured ITS legB reproducer state
// (candidate "13,6,6,5,4,9,5", hit 5): clamping every one of that state's
// three simultaneously-violated pairs ((Y,Snp), (Y,Q2Pt), (Snp,Q2Pt)) to
// exactly touch their pairwise Cauchy-Schwarz bound measurably shrinks the
// magnitude of the negative diagonal the subsequent measurement update's
// naive Kalman subtraction still manufactures (Q2Pt-Q2Pt: -0.0328 -> only
// -0.0067), but does not eliminate it -- proving pass 1+2 together are not
// a full PSD guarantee. The diagonal-abs step (pass 1) remains necessary
// AND, for the non-negative-diagonal part of the declared invariant,
// sufficient; achieving true full-matrix PSD validity would require an
// eigenvalue-based projection this function deliberately does not attempt
// -- that is a separate physics/model decision, not implemented here or
// anywhere in the frozen legacy code this migration reproduces.
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
  // The Stage-B forward slot carries signed q/pT for arbitrary absCharge
  // (not merely the unit-charge inverse pT the legacy name suggests).
  // getQ2Pt()/setQ2Pt() are the accepted post-migration names; getInvQPt()/
  // setInvQPt() remain as documented migration aliases over the same slot,
  // with no layout or ABI change.
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

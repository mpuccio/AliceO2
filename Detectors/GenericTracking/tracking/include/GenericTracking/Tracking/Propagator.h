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

/// \file Propagator.h
/// \brief Exact helix propagation between sensor planes in a uniform Bz.
///
/// A Propagator is configuration and nothing else: field, tolerances and the
/// mass hypothesis.  It is a plain struct so it copies into device constant memory as
/// one object, and every operation on it is a free function.  There are no
/// virtual calls: a field-map propagator later is a different type with the
/// same operations, chosen at compile time.

#ifndef GTRK_TRACKING_PROPAGATOR_H_
#define GTRK_TRACKING_PROPAGATOR_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/Trig.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "GPUCommonMath.h"
#include "MathUtils/detail/trigonometric.h"

#include <cstdint>

namespace o2::gtrk
{

/// Field, tolerances and the mass hypothesis. Copied by value into kernels.
struct Propagator {
  /// Longitudinal field in kGauss, O2's unit: 20 is 2 T.
  Scalar bzKiloGauss = 20.F;

  /// Smallest |normal . direction| at which the bound representation is still
  /// meaningful: at 0.01, about 89.4 degrees, a track crosses a 125 um barrel
  /// sensor over 1.25 cm, still through it; primaries from the diamond's tail
  /// reach 0.04 at the far end of barrel layer 0.
  Scalar minAbsCosIncidence = 0.01F;

  /// The mass hypothesis [GeV], positive: beta for the scattering, and the
  /// energy loss of the fitting stages.
  Scalar massGeV = static_cast<Scalar>(o2::constants::physics::MassPionCharged);

  /// Transverse curvature in 1/cm.
  ///
  /// kappa = (q/pT) * Bz * B2C, with q/pT = (q/p) / sin(theta). B2C carries
  /// its own sign, which is how the charge convention arrives without being
  /// restated here.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() ExactScalar curvature(ExactScalar qOverP, ExactScalar sinTheta) const
  {
    return (qOverP / sinTheta) * static_cast<ExactScalar>(bzKiloGauss) * static_cast<ExactScalar>(constants::math::B2C);
  }
};

/// What a propagation did and how far it went.
struct PropagationResult {
  PropagationStatus status = PropagationStatus::kNotFinite;
  Scalar pathLength = 0.F; ///< signed 3D path length [cm]
};

/// The uniform propagator a step is taken with.
///
/// `propagate`, `jacobianAnalytic` and `curvature` only ever see a uniform
/// field. Code above them asks for each step's propagator through these two
/// and hands on what they return. For a `Propagator` that is the propagator
/// itself: both return their argument by reference and look at nothing else,
/// so the default model pays nothing for them. A
/// field-map propagator is a different type and returns a `Propagator` by
/// value, with the step's own Bz.
///
/// The reference they return is the argument's, so a temporary propagator
/// would leave it dangling; the rvalue overloads are deleted for that.
///
/// The step from the chain's measurement \a from to its measurement \a to.
/// It takes indices rather than measurements so that a chain which reads its
/// rows from columns (`RowChain`) is not read for a field that does not look.
template <typename Geometry, typename Chain>
[[nodiscard]] GPUhdi() GPUAlwaysInline() const Propagator& stepPropagator(const Propagator& prop, const Geometry& /*geometry*/, const Chain& /*measurements*/,
                                                                          int /*from*/, int /*to*/)
{
  // The deleted rvalue overload below covers it; clang-tidy does not see that
  // for a template.
  return prop; // NOLINT(bugprone-return-const-ref-from-parameter)
}
template <typename Geometry, typename Chain>
const Propagator& stepPropagator(const Propagator&& prop, const Geometry& geometry, const Chain& measurements, int from, int to) = delete;

/// A step from \a point -- a bound state or a measurement -- on \a surface
/// towards a point not known yet: a layer crossing to predict, a candidate
/// cluster to try, the perigee. The point is passed whole, not its
/// coordinates, so the default model reads nothing for it either.
template <typename Point>
[[nodiscard]] GPUhdi() GPUAlwaysInline() const Propagator& pointPropagator(const Propagator& prop, const TrackingSurface& /*surface*/, const Point& /*point*/)
{
  // As for stepPropagator.
  return prop; // NOLINT(bugprone-return-const-ref-from-parameter)
}
template <typename Point>
const Propagator& pointPropagator(const Propagator&& prop, const TrackingSurface& surface, const Point& point) = delete;

namespace detail
{

/// The reference direction already evaluated by the FP64 plane solve.
/// Shared with the Jacobian to avoid evaluating it a second time.
struct PropagationTrig {
  ExactScalar sinPhi = 0.;
  ExactScalar cosPhi = 0.;
  ExactScalar sinTheta = 0.;
  ExactScalar cosTheta = 0.;
  /// The direction rotated to the arrival, as the incidence test had it: the
  /// covariance transport reuses it on every backend alike.
  Scalar sinPhiArrival = 0.F;
  Scalar cosPhiArrival = 0.F;
};

/// Global position of a bound state on its own surface.
GPUhdi() GPUAlwaysInline() void boundToGlobal(const TrackingSurface& surface, ExactScalar loc0, ExactScalar loc1, ExactScalar global[3])
{
  for (int i = 0; i < 3; ++i) {
    global[i] =
      static_cast<ExactScalar>(surface.origin[i]) + loc0 * static_cast<ExactScalar>(surface.axis0[i]) + loc1 * static_cast<ExactScalar>(surface.axis1[i]);
  }
}

/// Displacement of a bound state from the target surface origin.
///
/// Keeping the origin difference separate avoids forming a large global point
/// only to subtract another large global point in the intersection and final
/// projection. Both origins are floats, so their difference is exact, and what
/// the solve then carries is a distance between two sensors rather than a
/// distance from the interaction point.
GPUhdi() GPUAlwaysInline() void boundDisplacement(const TrackingSurface& from, const TrackingSurface& to, ExactScalar loc0, ExactScalar loc1,
                                                  ExactScalar delta[3])
{
  for (int i = 0; i < 3; ++i) {
    delta[i] = (static_cast<ExactScalar>(from.origin[i]) - static_cast<ExactScalar>(to.origin[i])) + loc0 * static_cast<ExactScalar>(from.axis0[i]) +
               loc1 * static_cast<ExactScalar>(from.axis1[i]);
  }
}

/// Where a crossing displaced from the target by \a fromTarget + (dxy, dz)
/// lands in the target plane, as local coordinates.
GPUhdi() GPUAlwaysInline() void crossingLocal(const TrackingSurface& to, const ExactScalar fromTarget[3], const ExactScalar dxy[2], ExactScalar dz,
                                              ExactScalar loc[2])
{
  const ExactScalar d[3] = {fromTarget[0] + dxy[0], fromTarget[1] + dxy[1], fromTarget[2] + dz};
  loc[0] = 0.;
  loc[1] = 0.;
  for (int i = 0; i < 3; ++i) {
    loc[0] += d[i] * static_cast<ExactScalar>(to.axis0[i]);
    loc[1] += d[i] * static_cast<ExactScalar>(to.axis1[i]);
  }
}

GPUhdi() GPUAlwaysInline() ExactScalar offCentreSquared(const TrackingSurface& to, const ExactScalar fromTarget[3], const ExactScalar dxy[2], ExactScalar dz)
{
  ExactScalar loc[2];
  crossingLocal(to, fromTarget, dxy, dz, loc);
  return (loc[0] * loc[0]) + (loc[1] * loc[1]);
}

/// A crossing as the other crossing rules see it (`CrossingRule`): its
/// squared distance from the target's centre, as `offCentreSquared`, and
/// whether it lies inside the sensor's bounds.
struct CrossingPlace {
  ExactScalar offCentre2 = 0.;
  bool inside = false;
};
GPUhdi() GPUAlwaysInline() CrossingPlace crossingPlace(const TrackingSurface& to, const ExactScalar fromTarget[3], const ExactScalar dxy[2], ExactScalar dz)
{
  ExactScalar loc[2];
  crossingLocal(to, fromTarget, dxy, dz, loc);
  return {.offCentre2 = (loc[0] * loc[0]) + (loc[1] * loc[1]), .inside = to.bounds.isInside(static_cast<Scalar>(loc[0]), static_cast<Scalar>(loc[1]))};
}

} // namespace detail

/// Advance \a state from the plane it is bound to onto \a to.
///
/// The intersection and the path-length solve run in ExactScalar and narrow on
/// return. The budget is 1 um across a detector reaching |z| = 300 cm, where a
/// float ulp is already 3.05e-5 cm, and the incidence cut admits crossings at
/// |n . d| = 0.01, which turns an error in the distance to the target plane
/// into a hundred times that along it. Ten nanometres of headroom in the solve,
/// 1e-6 cm, is below a float's resolution at this range, with or without an
/// error-free transform, so the solve is where the FP64 goes.
///
/// Everything outside the solve is Scalar, because it either ends up in a float
/// record or only has to clear a tolerance: the final direction, the incidence
/// test, the bounds test. Inputs are taken as finite -- they come from the fit,
/// not from a file -- so only results that can legitimately run away, the
/// second root of a straight track above all, are tested.
///
/// On kOk and kOutOfBounds the state is updated; on every other status it is
/// left untouched, so a caller that ignores the return value cannot silently
/// continue with a meaningless state. kOutOfBounds updates because the bound
/// is reported rather than enforced: a fit may want the intersection, a
/// search will not.
///
/// \a Rule is which of a barrel plane's two crossings is taken
/// (`CrossingRule`); the default takes the first root ahead, and the others
/// are their own instantiations.
template <CrossingRule Rule = CrossingRule::kFirstAhead>
GPUdi() PropagationResult propagate(const Propagator& prop, BoundTrackState& state, const TrackingSurface& from, const TrackingSurface& to,
                                    PropagationDirection dir, detail::PropagationTrig* referenceTrig = nullptr)
{
  using Math = o2::gpu::GPUCommonMath;
  PropagationResult res;

  const auto phi0 = static_cast<ExactScalar>(state.phi());
  const auto theta0 = static_cast<ExactScalar>(state.theta());

  // Evaluate the stored direction in Scalar, then retain FP64 throughout the
  // plane solve. The independent helix-closure sweep bounds the direction
  // rounding together with the solve/storage error.
  ExactScalar sinPhi = 0.;
  ExactScalar cosPhi = 0.;
  ExactScalar sinTheta = 0.;
  ExactScalar cosTheta = 0.;
  directionSincosOf(static_cast<Scalar>(phi0), sinPhi, cosPhi);
  directionSincosOf(static_cast<Scalar>(theta0), sinTheta, cosTheta);

  // A track exactly along the beam has no transverse motion and no bound
  // representation on any plane here.
  if (Math::Abs(sinTheta) < ExactScalar{1e-9}) {
    res.status = PropagationStatus::kNotFinite;
    return res;
  }

  const ExactScalar kappa = prop.curvature(static_cast<ExactScalar>(state.qOverP()), sinTheta);

  ExactScalar fromTarget[3] = {};
  detail::boundDisplacement(from, to, static_cast<ExactScalar>(state.loc0()), static_cast<ExactScalar>(state.loc1()), fromTarget);

  // Transverse direction and its left normal; the pair spans the plane the
  // circle lives in.
  const ExactScalar dT[2] = {cosPhi, sinPhi};
  const ExactScalar perp[2] = {-sinPhi, cosPhi};

  const auto nx = static_cast<ExactScalar>(to.normal[0]);
  const auto ny = static_cast<ExactScalar>(to.normal[1]);
  const auto nz = static_cast<ExactScalar>(to.normal[2]);

  ExactScalar path = 0.;   // 3D path length
  ExactScalar psi = 0.;    // turn angle
  ExactScalar dxy[2] = {}; // transverse displacement from r0
  // The turn's sine and cosine, kept because the chord already contains them.
  ExactScalar turnSin = 0.;
  ExactScalar turnCos = 1.;

  if (Math::Abs(nz) > ExactScalar{0.5}) {
    // A disc: the normal is along z, so z is linear in path length and the
    // solve is a division. Nothing about the helix enters it.
    if (Math::Abs(cosTheta) < static_cast<ExactScalar>(prop.minAbsCosIncidence)) {
      res.status = PropagationStatus::kGrazing;
      return res;
    }
    path = -fromTarget[2] / cosTheta;
    const ExactScalar ell = path * sinTheta;
    psi = kappa * ell;
    const auto sA = sinc<ExactScalar>(psi);
    const auto sB = versine<ExactScalar>(psi);
    dxy[0] = ell * (dT[0] * sA + perp[0] * sB);
    dxy[1] = ell * (dT[1] * sA + perp[1] * sB);
    // sin(psi) = psi sinc(psi) and cos(psi) = 1 - psi versine(psi), from the
    // series that is already evaluated rather than from a second sincos.
    turnSin = psi * sA;
    turnCos = 1. - psi * sB;
  } else {
    // A barrel sensor: the normal has no z component, so the plane projects to
    // a *line* in the transverse plane while the trajectory projects to a
    // circle. Circle-line intersection is a quadratic, and the whole solve is
    // closed form.
    //
    // Writing the circle as |delta|^2 = 2 delta . e about the start point --
    // rather than about its centre, which runs away as the curvature vanishes
    // -- gives
    //
    //     A mu^2 - 2 B mu + C = 0,   B = dT . n,  C = a (A a - 2 dT . t)
    //
    // in the offset mu along the line, with a the signed distance to the plane
    // and A = kappa |n|^2. Both coefficients stay finite as kappa -> 0, and
    // taking the root in the form C / (B + sign(B) sqrt(D)) removes the
    // cancellation *and* the straight-line special case: at kappa = 0 it
    // collapses exactly to mu = -a (dT.t) / (dT.n).
    const ExactScalar tx = -ny;
    const ExactScalar ty = nx;

    // The stored normal is a float triplet, so its norm is one only to within
    // a float's last bit. Carrying |n|^2 costs one divide and makes the solve
    // exact for the plane as stored, instead of scaling that last bit by the
    // hundredfold a grazing crossing applies to the distance to the plane.
    const ExactScalar normal2 = nx * nx + ny * ny;
    const ExactScalar a = (nx * fromTarget[0] + ny * fromTarget[1]) / normal2;
    const ExactScalar quadratic = kappa * normal2;

    const ExactScalar bCoef = dT[0] * nx + dT[1] * ny;
    const ExactScalar tangent = dT[0] * tx + dT[1] * ty;
    const ExactScalar cCoef = a * (quadratic * a - ExactScalar{2} * tangent);

    // sqrt(D) is the transverse incidence at the crossing, so D cancels only
    // where the circle is tangent to the line -- and there the crossing does
    // not clear minAbsCosIncidence anyway. Inside the accepted band the
    // cancellation is bounded by that cut, 1e-4 of the terms, which leaves a
    // double twelve digits.
    const ExactScalar disc = bCoef * bCoef - quadratic * cCoef;
    if (disc < ExactScalar{0}) {
      // The circle never reaches the line: the track curls before the target
      // radius. A fact about the track, not a failure.
      res.status = PropagationStatus::kNoIntersection;
      return res;
    }
    const ExactScalar den = bCoef + (bCoef >= ExactScalar{0} ? ExactScalar{1} : ExactScalar{-1}) * o2::math_utils::detail::sqrt(disc);
    if (Math::Abs(den) < ExactScalar{1e-30}) {
      res.status = PropagationStatus::kNoIntersection;
      return res;
    }

    // A line cuts the circle twice, and which crossing is the sensor's is a
    // question about the geometry, not about the arithmetic. C / den is the
    // root of smaller magnitude: the stable form, and the right crossing
    // whenever the turn is small. Once the target's normal is far enough from
    // radial at the crossing -- a tilted stave, a cluster away from the module
    // centre -- the two swap, and at a 30 to 100 degree turn between adjacent
    // layers, which is what a 0.3 GeV track makes, the small root becomes a
    // crossing of the plane's infinite extension tens of cm from the sensor.
    //
    // So the small root is tried first and kept whenever it is usable; the
    // second root, den / A, is reached only where the first is not usable at
    // all. As kappa -> 0 the second runs away and stops being finite, which is how
    // the straight line keeps having one solution without a special case; it
    // is also why it is not divided out before it is needed.
    //
    // The second root is taken only when it lands *nearer the target's centre*
    // than the first, which is the diagnosis stated as a condition: the first
    // root is the extension and the second is the sensor. Without that guard
    // it also answers a question the wrong-direction status was right to
    // decline -- a track pointing away from the target reaches it after
    // sweeping most of a circle, 20 to 60 m of path, and no linearisation of
    // that is worth anything.
    //
    // Not addressed by the default: the small root can also be the wrong
    // crossing while still running the right way, in which case nothing is
    // reported and the update simply sees a large residual. A looper's inward
    // refit from its way back onto an outer barrel layer of its way out does
    // exactly that, and fails a step or two later. The other rules (`CrossingRule`) look at the second root even when
    // the first lies ahead: kNearerCentre takes it where it lies ahead and
    // nearer the centre, kInsideFirst only where the first is outside the
    // sensor and the second inside or, neither inside, nearer the centre.
    // Selecting on the bounds moves the chosen root at the sensor edge, where
    // testCovarianceTransport takes a derivative by finite differences; the
    // ties go to the first root.
    // The default: the first root, unless it is not usable and the second
    // lands nearer the centre. The other rules, each its own instantiation:
    // the second where the rule moves off the first.
    [[maybe_unused]] ExactScalar firstOffCentre = 0.;
    [[maybe_unused]] detail::CrossingPlace first;
    [[maybe_unused]] bool firstAhead = false;
    for (int root = 0; root < 2; ++root) {
      const ExactScalar mu = root == 0 ? cCoef / den : den / quadratic;
      const ExactScalar tryDxy[2] = {-a * nx + mu * tx, -a * ny + mu * ty};

      // Recover the turn from the chord: sin(psi) = kappa (delta . dT) and
      // cos(psi) = 1 - kappa (delta . perp), both exact and both bounded. It
      // is the principal turn, so this holds for |turn| <= pi over one step,
      // well above a turn between adjacent layers.
      const ExactScalar alongD = (tryDxy[0] * dT[0]) + (tryDxy[1] * dT[1]);
      const ExactScalar alongP = (tryDxy[0] * perp[0]) + (tryDxy[1] * perp[1]);
      const ExactScalar trySin = kappa * alongD;
      const ExactScalar tryCos = ExactScalar{1} - kappa * alongP;
      const auto tryPsi = atan2Of(trySin, tryCos);

      // ell = psi / kappa is 0/0 for a straight track; ell = (delta.dT)/sinc(psi)
      // is the same number and is regular there. The quotient is the one to
      // take near zero only: sinc vanishes at a half turn, where psi / kappa
      // is the well-conditioned form.
      const ExactScalar ell = Math::Abs(tryPsi) < ExactScalar{0.01} ? alongD / sinc<ExactScalar>(tryPsi) : tryPsi / kappa;
      const ExactScalar tryPath = ell / sinTheta;

      const bool usable = isFiniteScalar(static_cast<Scalar>(tryPath)) && tryPath * static_cast<ExactScalar>(static_cast<int>(dir)) >= ExactScalar{0};
      bool take = root == 0;
      if constexpr (Rule == CrossingRule::kFirstAhead) {
        if (root != 0) {
          take = usable && detail::offCentreSquared(to, fromTarget, tryDxy, tryPath * cosTheta) < firstOffCentre;
        }
      } else {
        const auto place = detail::crossingPlace(to, fromTarget, tryDxy, tryPath * cosTheta);
        if (root == 0) {
          first = place;
          firstAhead = usable;
        } else if constexpr (Rule == CrossingRule::kInsideFirst) {
          // The first, if ahead, is outside: the second if it is inside or,
          // neither inside, nearer the centre; with the first behind, as the
          // default.
          take = usable && (place.offCentre2 < first.offCentre2 || (firstAhead && place.inside));
        } else {
          take = usable && place.offCentre2 < first.offCentre2;
        }
      }
      if (take) {
        dxy[0] = tryDxy[0];
        dxy[1] = tryDxy[1];
        psi = tryPsi;
        path = tryPath;
        turnSin = trySin;
        turnCos = tryCos;
      }
      if constexpr (Rule == CrossingRule::kFirstAhead) {
        if (root == 0) {
          if (usable) {
            break; // the overwhelming majority: nothing below is reached or paid for
          }
          firstOffCentre = detail::offCentreSquared(to, fromTarget, tryDxy, tryPath * cosTheta);
        }
      } else if constexpr (Rule == CrossingRule::kInsideFirst) {
        if (root == 0 && usable && first.inside) {
          break;
        }
      }
    }
  }

  // Asked as a float, because a float is what gets stored: a double too large
  // to narrow is exactly the case this has to reject.
  if (!isFiniteScalar(static_cast<Scalar>(path))) {
    res.status = PropagationStatus::kNotFinite;
    return res;
  }

  const auto wanted = static_cast<ExactScalar>(static_cast<int>(dir));
  if (path * wanted < ExactScalar{0}) {
    res.status = PropagationStatus::kWrongDirection;
    res.pathLength = static_cast<Scalar>(path);
    return res;
  }

  // The direction only rotates in azimuth; theta and |p| are invariant in a
  // uniform longitudinal field, which is why neither is recomputed. The turn's
  // sine and cosine came out of the chord, so this is a rotation rather than a
  // second sincos -- and it is asked in Scalar, because the azimuth is stored
  // as one and the incidence only has to clear a tolerance.
  ExactScalar phi1 = phi0 + psi;
  math_utils::detail::bringToPMPi<ExactScalar>(phi1);

  const auto sinTurn = static_cast<Scalar>(turnSin);
  const auto cosTurn = static_cast<Scalar>(turnCos);
  const Scalar sinPhi1 = static_cast<Scalar>(sinPhi) * cosTurn + static_cast<Scalar>(cosPhi) * sinTurn;
  const Scalar cosPhi1 = static_cast<Scalar>(cosPhi) * cosTurn - static_cast<Scalar>(sinPhi) * sinTurn;
  const auto sinThetaF = static_cast<Scalar>(sinTheta);
  const Scalar dir1[3] = {sinThetaF * cosPhi1, sinThetaF * sinPhi1, static_cast<Scalar>(cosTheta)};

  const Scalar cosInc = cosIncidence(to, dir1);
  countWork(WorkCounter::kPropagations);
  if (Math::Abs(cosInc) <= 0.1F) {
    countWork(WorkCounter::kPropagationsGrazing);
  }
  if (Math::Abs(cosInc) < prop.minAbsCosIncidence) {
    res.status = PropagationStatus::kGrazing;
    res.pathLength = static_cast<Scalar>(path);
    return res;
  }

  // Bind to the target plane.  The out-of-plane component is not computed: it
  // is zero by construction, and checking it at run time would be checking the
  // algebra rather than the data.  testPropagator.cxx checks it instead.
  const ExactScalar d0 = fromTarget[0] + dxy[0];
  const ExactScalar d1 = fromTarget[1] + dxy[1];
  const ExactScalar d2 = fromTarget[2] + path * cosTheta;
  const ExactScalar loc0 = d0 * static_cast<ExactScalar>(to.axis0[0]) + d1 * static_cast<ExactScalar>(to.axis0[1]) + d2 * static_cast<ExactScalar>(to.axis0[2]);
  const ExactScalar loc1 = d0 * static_cast<ExactScalar>(to.axis1[0]) + d1 * static_cast<ExactScalar>(to.axis1[1]) + d2 * static_cast<ExactScalar>(to.axis1[2]);

  // What is about to be stored, asked as floats like the path. A disc's path
  // is finite whatever the azimuth or the curvature, so a NaN in either shows
  // up here first: in the displacement, and so in the bound coordinates, and
  // the incidence test above does not reject a NaN. Every comparison on the
  // way down evaluates false for one, so without this test a NaN state would
  // come back as kOutOfBounds with the NaNs written.
  if (!isFiniteScalar(static_cast<Scalar>(loc0)) || !isFiniteScalar(static_cast<Scalar>(loc1)) || !isFiniteScalar(static_cast<Scalar>(phi1))) {
    res.status = PropagationStatus::kNotFinite;
    return res;
  }

  state.par[kLoc0] = static_cast<Scalar>(loc0);
  state.par[kLoc1] = static_cast<Scalar>(loc1);
  state.par[kPhi] = static_cast<Scalar>(phi1);
  // theta and q/p are untouched on purpose.
  state.surface = to.id;
  state.surfacePhi = to.surfacePhi;

  res.pathLength = static_cast<Scalar>(path);
  res.status = to.bounds.isInside(state.par[kLoc0], state.par[kLoc1]) ? PropagationStatus::kOk : PropagationStatus::kOutOfBounds;
  if (referenceTrig != nullptr) {
    *referenceTrig = {.sinPhi = sinPhi, .cosPhi = cosPhi, .sinTheta = sinTheta, .cosTheta = cosTheta, .sinPhiArrival = sinPhi1, .cosPhiArrival = cosPhi1};
  }
  return res;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_PROPAGATOR_H_

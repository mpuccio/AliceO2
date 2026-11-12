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

/// \file Follow.h
/// \brief A track followed beyond its rows, outward and then inward, as a CKF
///        follows it: the crossing and the update,
///        one step's search, the observer, the stage (`FollowStage`), the
///        passes after it, and the host's orchestration (`Follow`).
///
/// After an iteration's refit each of its tracks is followed, and, with
/// `followEarlier`, every earlier iteration's track first. The rows proposed
/// then pass the bidding between the tracks followed in their order, the
/// kink steps' witnesses, the combination, the refit under the road fit's
/// chi2 cut, and the marks: the bid is over one iteration's tracks, and the
/// marks freeze the rows kept, so that later iterations do not build on them.

#ifndef GTRK_TRACKING_FOLLOW_H_
#define GTRK_TRACKING_FOLLOW_H_

#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Tracking/CovarianceTransport.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"
#include "GenericTracking/Tracking/MeasurementUpdate.h"
#include "GenericTracking/Tracking/Profile.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/SensorCells.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackFitter.h"
#include "GenericTracking/Tracking/TrackRefit.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/Trig.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include "GPUCommonMath.h"
#include "MathUtils/detail/trigonometric.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numbers>
#include <span>
#include <type_traits>

namespace o2::gtrk
{

/// A track's intersection with a layer's nominal cylinder or disc plane, in
/// the two coordinates the measurement grid indexes. The search is wider
/// than the difference between that nominal surface and an individual sensor;
/// the candidate is propagated to its actual sensor before it is chosen.
struct LayerCrossingPrediction {
  Scalar phi = 0.F;
  Scalar second = 0.F;
  PropagationStatus status = PropagationStatus::kNotFinite;
};

namespace follow
{

GPUhdi() GPUAlwaysInline() ExactScalar directedTurn(ExactScalar turn, ExactScalar curvature, PropagationDirection direction)
{
  const auto wanted = static_cast<ExactScalar>(static_cast<int>(direction));
  if ((turn / curvature) * wanted < 0.) {
    turn += wanted * (curvature > 0. ? 2 * std::numbers::pi_v<ExactScalar> : -(2 * std::numbers::pi_v<ExactScalar>));
  }
  return turn;
}

/// Whether a crossing \a turn radians along the helix is past half a turn.
/// Such a crossing is a looper coming round again -- a track at eta ~ 0 reaches
/// a forward disc after about ten turns -- and not the next layer of the same
/// arc; its transported covariance is also past what float can factorise.
GPUhdi() GPUAlwaysInline() bool beyondHalfTurn(ExactScalar turn) { return o2::gpu::CAMath::Abs(turn) > std::numbers::pi_v<ExactScalar>; }

GPUhdi() GPUAlwaysInline() bool nearerPath(ExactScalar candidate, ExactScalar chosen) { return exactlyEqual(chosen, ExactScalar{0}) || o2::gpu::CAMath::Abs(candidate) < o2::gpu::CAMath::Abs(chosen); }

} // namespace follow

/// Predict where \a state crosses \a layer, without choosing a sensor. The
/// nominal radius or z is copied into IterationParameters on the host, so the
/// runtime-indexed path reads no detector constexpr table on a device.
template <typename Prop, DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() LayerCrossingPrediction predictLayerCrossing(const Prop& prop, const IterationParameters<Det>& parameters, const TrackingSurface& from, const BoundTrackState& state, LayerID layer, PropagationDirection direction)
{
  using Math = o2::gpu::CAMath;
  LayerCrossingPrediction result;
  const auto& propagator = pointPropagator(prop, from, state);

  Scalar startScalar[3]{};
  surfaceToGlobal(from, state.loc0(), state.loc1(), startScalar);
  const ExactScalar start[3] = {static_cast<ExactScalar>(startScalar[0]), static_cast<ExactScalar>(startScalar[1]), static_cast<ExactScalar>(startScalar[2])};
  const auto phi0 = static_cast<ExactScalar>(state.phi());
  const auto theta = static_cast<ExactScalar>(state.theta());
  ExactScalar sinPhi = 0.;
  ExactScalar cosPhi = 0.;
  ExactScalar sinTheta = 0.;
  ExactScalar cosTheta = 0.;
  sincosOf(phi0, sinPhi, cosPhi);
  sincosOf(theta, sinTheta, cosTheta);
  if (Math::Abs(sinTheta) < 1.e-9) {
    result.status = PropagationStatus::kNotFinite;
    return result;
  }

  const ExactScalar curvature = propagator.curvature(static_cast<ExactScalar>(state.qOverP()), sinTheta);
  const auto directionValue = static_cast<ExactScalar>(static_cast<int>(direction));
  ExactScalar path = 0.;
  ExactScalar x = start[0];
  ExactScalar y = start[1];
  ExactScalar z = start[2];

  if (Detector<Det>::hasDiscs && !Detector<Det>::isBarrelLayer(layer)) {
    if (Math::Abs(cosTheta) < static_cast<ExactScalar>(propagator.minAbsCosIncidence)) {
      result.status = PropagationStatus::kGrazing;
      return result;
    }
    path = (static_cast<ExactScalar>(parameters.layerSpans[layer].position) - start[2]) / cosTheta;
    if (path * directionValue < 0.) {
      result.status = PropagationStatus::kWrongDirection;
      return result;
    }
    const ExactScalar turn = curvature * path * sinTheta;
    if (follow::beyondHalfTurn(turn)) {
      result.status = PropagationStatus::kTooManyTurns;
      return result;
    }
    const ExactScalar transverse = path * sinTheta;
    x += transverse * ((cosPhi * sinc<ExactScalar>(turn)) - (sinPhi * versine<ExactScalar>(turn)));
    y += transverse * ((sinPhi * sinc<ExactScalar>(turn)) + (cosPhi * versine<ExactScalar>(turn)));
    z = static_cast<ExactScalar>(parameters.layerSpans[layer].position);
  } else {
    const auto radius = static_cast<ExactScalar>(parameters.layerSpans[layer].position);
    const ExactScalar dot = (start[0] * cosPhi) + (start[1] * sinPhi);
    const ExactScalar radial = (start[0] * start[0]) + (start[1] * start[1]) - (radius * radius);
    if (Math::Abs(curvature) < 1.e-12) {
      const ExactScalar discriminant = (dot * dot) - radial;
      if (discriminant < 0.) {
        result.status = PropagationStatus::kNoIntersection;
        return result;
      }
      const auto root = static_cast<ExactScalar>(Math::Sqrt(static_cast<Scalar>(discriminant)));
      const ExactScalar first = -dot - root;
      const ExactScalar second = -dot + root;
      ExactScalar transverse = 0.;
      if (first * directionValue >= 0.) {
        transverse = first;
      }
      if (second * directionValue >= 0. && follow::nearerPath(second, transverse)) {
        transverse = second;
      }
      if (exactlyEqual(transverse, ExactScalar{0})) {
        result.status = PropagationStatus::kWrongDirection;
        return result;
      }
      path = transverse / sinTheta;
      x += transverse * cosPhi;
      y += transverse * sinPhi;
      z += path * cosTheta;
    } else {
      const ExactScalar inverse = 1. / curvature;
      const ExactScalar centreX = start[0] - (sinPhi * inverse);
      const ExactScalar centreY = start[1] + (cosPhi * inverse);
      const auto centreDistance = static_cast<ExactScalar>(Math::Sqrt(static_cast<Scalar>((centreX * centreX) + (centreY * centreY))));
      const ExactScalar circleRadius = Math::Abs(inverse);
      if (centreDistance > radius + circleRadius || centreDistance < Math::Abs(radius - circleRadius) || !(centreDistance > 0.)) {
        result.status = PropagationStatus::kNoIntersection;
        return result;
      }
      const ExactScalar along = ((radius * radius) - (circleRadius * circleRadius) + (centreDistance * centreDistance)) / (2. * centreDistance);
      const ExactScalar heightSquared = (radius * radius) - (along * along);
      if (heightSquared < 0.) {
        result.status = PropagationStatus::kNoIntersection;
        return result;
      }
      const auto height = static_cast<ExactScalar>(Math::Sqrt(static_cast<Scalar>(heightSquared)));
      const ExactScalar unitX = centreX / centreDistance;
      const ExactScalar unitY = centreY / centreDistance;
      ExactScalar chosenPath = 0.;
      ExactScalar chosenX = 0.;
      ExactScalar chosenY = 0.;
      for (int side = -1; side <= 1; side += 2) {
        const auto sideValue = static_cast<ExactScalar>(side);
        const ExactScalar candidateX = (along * unitX) - (sideValue * height * unitY);
        const ExactScalar candidateY = (along * unitY) + (sideValue * height * unitX);
        const ExactScalar deltaX = candidateX - start[0];
        const ExactScalar deltaY = candidateY - start[1];
        const ExactScalar alongDirection = (deltaX * cosPhi) + (deltaY * sinPhi);
        const ExactScalar alongNormal = (-deltaX * sinPhi) + (deltaY * cosPhi);
        ExactScalar turn = atan2Of(curvature * alongDirection, 1. - (curvature * alongNormal));
        turn = follow::directedTurn(turn, curvature, direction);
        const ExactScalar candidatePath = turn / (curvature * sinTheta);
        if (candidatePath * directionValue >= 0. && follow::nearerPath(candidatePath, chosenPath)) {
          chosenPath = candidatePath;
          chosenX = candidateX;
          chosenY = candidateY;
        }
      }
      if (exactlyEqual(chosenPath, ExactScalar{0})) {
        result.status = PropagationStatus::kWrongDirection;
        return result;
      }
      // The nearer of the two circle crossings: if it is past half a turn, so
      // is the other.
      if (follow::beyondHalfTurn(curvature * chosenPath * sinTheta)) {
        result.status = PropagationStatus::kTooManyTurns;
        return result;
      }
      path = chosenPath;
      x = chosenX;
      y = chosenY;
      z += path * cosTheta;
    }
  }

  result.phi = o2::math_utils::detail::to02Pi<Scalar>(Math::ATan2(static_cast<Scalar>(y), static_cast<Scalar>(x)));
  result.second = Detector<Det>::isBarrelLayer(layer) ? static_cast<Scalar>(z) : Math::Sqrt(static_cast<Scalar>((x * x) + (y * y)));
  result.status = isFiniteScalar(result.phi) && isFiniteScalar(result.second) ? PropagationStatus::kOk : PropagationStatus::kNotFinite;
  return result;
}

/// Rows a track's follow would add in each direction. Outward rows are in
/// increasing road order; inward rows are in the order they are found,
/// outermost to innermost, and are reversed only when the combined RowChain
/// is read. With the end decided again (`followEnd`), \a replacedEnd is the
/// track's outermost row where the outward rows take its place; it goes only
/// with them, if outward rows are kept after the passes.
template <std::size_t K>
struct FollowProposal {
  std::uint32_t outward[K]{};
  std::uint32_t inward[K]{};
  std::uint32_t replacedEnd = std::numeric_limits<std::uint32_t>::max();
  std::uint8_t outwardCount = 0U;
  std::uint8_t inwardCount = 0U;
};

namespace follow
{

/// A track's kink steps in one follow: how many each direction took -- the budget, `maxKinks`, is the
/// follow's, outward and inward together -- and, per
/// direction, the rows that stand on the next one, the witness that
/// confirmed them (bit i for row i). A row the beam line confirmed stands on
/// none.
struct KinkSteps {
  std::uint64_t outwardWitnessed = 0U;
  std::uint64_t inwardWitnessed = 0U;
  std::uint8_t outward = 0U;
  std::uint8_t inward = 0U;
};

/// How many of a direction's first \a count rows to keep so that every kept
/// kink step keeps its witness: the rows from the last one standing on a row
/// not kept are dropped with it.
[[nodiscard]] GPUhdi() GPUAlwaysInline() constexpr std::uint8_t keptWithWitnesses(std::uint8_t count, std::uint64_t witnessed)
{
  while (count != 0U && ((witnessed >> (count - 1U)) & 1U) != 0U) {
    --count;
  }
  return count;
}

/// \a input, its covariance \a inputFactors (CovarianceFactors), carried to
/// the plane \a to with the material's noise and energy loss between, as
/// `detail::filterStep` carries a fit's state (`detail::transportStep`):
/// \a output and \a outputFactors. Unlike the fit, it takes the field at
/// \a input (`pointPropagator`) and refuses a step beyond half a turn.
/// `cov` of neither state is read or written. False if the step is not usable.
/// The follow's window carries the state to a layer's nominal surface with
/// it (`predictedWindow`), and `updateCandidate` is this plus the update.
/// With \a Scaled, \a noiseScale multiplies the step's theta0: the kink
/// allowance's step; without, the nominal step,
/// which never multiplies (`detail::scaledTheta0`). The body, inlined where
/// it is called; the follow calls `transportCandidate`.
template <bool Scaled, typename Prop, DetectorLayout Det, typename Material>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool transportCandidateInline(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const TrackingSurface& to, PropagationDirection direction, const BoundTrackState& input, const CovarianceFactors& inputFactors, BoundTrackState& output, CovarianceFactors& outputFactors, Scalar noiseScale)
{
  const auto& from = geometry.planeOf(input);
  const auto& propagator = pointPropagator(prop, from, input);
  Scalar start[3]{};
  surfaceToGlobal(from, input.loc0(), input.loc1(), start);
  SolenoidJacobian jacobian;
  PropagationResult propagated;
  if (jacobianAnalytic(propagator, input, from, to, direction, jacobian, &output, &propagated) != PropagationStatus::kOk) {
    return false;
  }
  // The turn is kappa times the transverse path, (q/p) Bz B2C times the 3D
  // path: curvature() at sin(theta) = 1 against the 3D length.
  if (beyondHalfTurn(propagator.curvature(static_cast<ExactScalar>(input.qOverP()), 1.) * static_cast<ExactScalar>(propagated.pathLength))) {
    return false;
  }
  Scalar end[3]{};
  surfaceToGlobal(to, output.loc0(), output.loc1(), end);
  const auto step = detail::stepMaterial(material, start, end);
  const Scalar momentum = 1.F / o2::gpu::CAMath::Abs(input.qOverP());
  const Scalar betaP = betaMomentum(momentum, propagator.massGeV);
  // As filterStep: the nearer half's noise before the transport, at the
  // direction the step leaves with.
  const auto before = scatteringNoise(detail::scaledTheta0<Scaled>(highlandTheta0(step.nearerX2X0, betaP), noiseScale), input.theta());
  outputFactors = inputFactors;
  return detail::transportStep<Scaled>(propagator, step, betaP, before, direction, jacobian, output, outputFactors, noiseScale) == PropagationStatus::kOk;
}

/// `transportCandidateInline`, never inlined on the device: the one copy
/// every search shares. Only
/// `FollowStage::passingOn` takes the body inline (`updateCandidate`'s
/// \a Inlined).
template <bool Scaled, typename Prop, DetectorLayout Det, typename Material>
[[nodiscard]] GPUhdnii() bool transportCandidate(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const TrackingSurface& to, PropagationDirection direction, const BoundTrackState& input, const CovarianceFactors& inputFactors, BoundTrackState& output, CovarianceFactors& outputFactors, Scalar noiseScale)
{
  return transportCandidateInline<Scaled>(prop, geometry, material, to, direction, input, inputFactors, output, outputFactors, noiseScale);
}

/// \a input, its covariance \a inputFactors (CovarianceFactors), carried to
/// \a measurement (`transportCandidate`; with \a Scaled the step's scattering
/// scaled by \a noiseScale, its theta0, so the variances by its square: the
/// kink allowance's step) and updated with it:
/// \a output and \a outputFactors, and the update's \a chi2 -- and, with
/// \a Determinant, its innovation covariance's determinant into
/// \a determinant (the follow's EndDensity). `cov` of neither state is read
/// or written. False if the step or the update is not usable. With
/// \a Inlined the transport is inlined (`transportCandidateInline`), else
/// called.
template <bool Scaled, bool Determinant = false, bool Inlined = false, typename Prop, DetectorLayout Det, typename Material>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool updateCandidate(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, const BoundTrackState& input, const CovarianceFactors& inputFactors, BoundTrackState& output, CovarianceFactors& outputFactors, Scalar& chi2, Scalar noiseScale, [[maybe_unused]] Scalar* determinant = nullptr)
{
  bool carried = false;
  if constexpr (Inlined) {
    carried = transportCandidateInline<Scaled>(prop, geometry, material, geometry.planeOf(measurement), direction, input, inputFactors, output, outputFactors, noiseScale);
  } else {
    carried = transportCandidate<Scaled>(prop, geometry, material, geometry.planeOf(measurement), direction, input, inputFactors, output, outputFactors, noiseScale);
  }
  if (!carried) {
    return false;
  }
  const auto update = updateWithMeasurement<Determinant>(output, outputFactors, measurement, determinant);
  chi2 = update.chi2;
  // A negative chi2 would otherwise win the best-candidate comparison: the
  // same test the fit passes make (`updateUsable`).
  return updateUsable(update, outputFactors);
}

/// The half-widths of the planes `predictedWindow` builds: no bound [cm].
inline constexpr Scalar kOpenBound = 1.e4F;

/// The window of one step of a follow, narrowing \a window about the nominal crossing \a nominal. The
/// state is
/// carried to the plane touching the layer there, with the material's noise
/// (`transportCandidate`, with \a Scaled the step's scattering scaled by
/// \a noiseScale), and the window's half-widths are:
///  - the gate's sigma of the prediction's spread, its covariance plus the
///    pair windows' resolution term, projected onto phi and the second
///    coordinate;
///  - plus the layer's envelope about its nominal position, its
///    measurements' radii on a barrel layer or z on a disc, times the track's
///    slope there: a sensor anywhere in it, tilted, overlapping or curved, is
///    crossed within that;
///  - at most \a window's as given.
/// Where the transport fails, \a window stands, and it returns false.
/// Which candidate is taken is still the update's chi2 on its own sensor.
/// \a reach, where given, gets the envelope's term alone, azimuth then the
/// second coordinate: how far from the nominal crossing the track itself
/// may cross a sensor of the layer (the follow's hole test).
template <bool Scaled, typename Prop, DetectorLayout Det, typename Material>
GPUhdi() bool predictedWindow(const Prop& prop, const IterationParameters<Det>& parameters, const TrackingGeometry<Det>& geometry, const Material& material, const BoundTrackState& state, const CovarianceFactors& factors, LayerID target, PropagationDirection direction, const LayerCrossingPrediction& nominal, CandidateWindow& window, Scalar noiseScale, Scalar* reach = nullptr)
{
  using Math = o2::gpu::CAMath;
  const bool barrel = ((parameters.barrelLayers >> target) & 1U) != 0U;
  // The plane touching the layer at the nominal crossing: on a barrel layer
  // loc0 along r phi and loc1 along z, on a disc loc0 along r and loc1
  // along r phi.
  const Scalar position = parameters.layerSpans[target].position;
  const Scalar radius = barrel ? position : nominal.second;
  if (!(radius > 0.F)) {
    return false;
  }
  Scalar sinPhi = 0.F;
  Scalar cosPhi = 0.F;
  sincosOf(nominal.phi, sinPhi, cosPhi);
  TrackingSurface plane;
  plane.origin[0] = radius * cosPhi;
  plane.origin[1] = radius * sinPhi;
  plane.origin[2] = barrel ? nominal.second : position;
  const Scalar radial[3] = {cosPhi, sinPhi, 0.F};
  const Scalar azimuthal[3] = {-sinPhi, cosPhi, 0.F};
  const Scalar along[3] = {0.F, 0.F, 1.F};
  for (int i = 0; i < 3; ++i) {
    plane.axis0[i] = barrel ? azimuthal[i] : radial[i];
    plane.axis1[i] = barrel ? along[i] : azimuthal[i];
    plane.normal[i] = barrel ? radial[i] : along[i];
  }
  plane.bounds = {.halfX = kOpenBound, .halfY = kOpenBound};
  plane.layer = target;
  BoundTrackState predicted;
  CovarianceFactors predictedFactors;
  if (!transportCandidate<Scaled>(prop, geometry, material, plane, direction, state, factors, predicted, predictedFactors, noiseScale)) {
    return false;
  }
  covarianceFromFactors(predictedFactors, predicted.cov);
  const Scalar variance0 = predicted.variance(kLoc0) + parameters.followPositionVariance;
  const Scalar variance1 = predicted.variance(kLoc1) + parameters.followPositionVariance;
  if (!(variance0 > 0.F) || !(variance1 > 0.F)) {
    return false;
  }
  // The window spans what the update's gate admits on one axis,
  // sqrt(clusterChi2Max), not the pair cut's nSigma: a row the gate takes is
  // never outside the window.
  const Scalar nSigma = Math::Sqrt(parameters.clusterChi2Max);
  const Scalar phiSpread = nSigma * Math::Sqrt(barrel ? variance0 : variance1) / radius;
  const Scalar secondSpread = nSigma * Math::Sqrt(barrel ? variance1 : variance0);
  // How far the crossing moves per unit of the layer's position, along the
  // track there: alpha is its transverse direction from the radial one. The
  // transverse path is 1 / cos(alpha) per unit of radius on a barrel layer,
  // tan(theta) per unit of z on a disc.
  Scalar sinAlpha = 0.F;
  Scalar cosAlpha = 0.F;
  sincosOf(predicted.phi() - nominal.phi, sinAlpha, cosAlpha);
  Scalar sinTheta = 0.F;
  Scalar cosTheta = 0.F;
  sincosOf(predicted.theta(), sinTheta, cosTheta);
  const Scalar transverse = barrel ? 1.F / Math::Abs(cosAlpha) : Math::Abs(sinTheta / cosTheta);
  const Scalar phiSlope = transverse * Math::Abs(sinAlpha) / radius;
  const Scalar secondSlope = transverse * (barrel ? Math::Abs(cosTheta / sinTheta) : Math::Abs(cosAlpha));
  const auto& span = parameters.layerSpans[target];
  const Scalar envelope = Math::Max(span.positionMax - position, position - span.positionMin);
  // The window given first: a NaN or infinite slope leaves it.
  window.phiHalfWidth = Math::Min(window.phiHalfWidth, phiSpread + (envelope * phiSlope));
  window.secondHalfWidth = Math::Min(window.secondHalfWidth, secondSpread + (envelope * secondSlope));
  if (reach != nullptr) {
    reach[0] = envelope * phiSlope;
    reach[1] = envelope * secondSlope;
  }
  return true;
}

/// The nominal crossing \a prediction of \a layer as a point [cm], into
/// \a point: on a barrel layer its nominal radius at the crossing's azimuth
/// and z, on a disc its radius and azimuth at the disc's z.
template <DetectorLayout Det>
GPUhdi() GPUAlwaysInline() void crossingPoint(const IterationParameters<Det>& parameters, LayerID layer, const LayerCrossingPrediction& prediction, Scalar* point)
{
  const bool barrel = ((parameters.barrelLayers >> layer) & 1U) != 0U;
  const Scalar position = parameters.layerSpans[layer].position;
  const Scalar r = barrel ? position : prediction.second;
  Scalar sinPhi = 0.F;
  Scalar cosPhi = 0.F;
  sincosOf(prediction.phi, sinPhi, cosPhi);
  point[0] = r * cosPhi;
  point[1] = r * sinPhi;
  point[2] = barrel ? prediction.second : position;
}

/// The centre of a step's window on \a layer, in the rows' coordinates: the
/// nominal crossing \a prediction, which is the detector frame's, as its
/// point's azimuth and second coordinate about the beam line, as the
/// measurements have theirs (`gridPoint`).
template <DetectorLayout Det>
GPUhdi() GPUAlwaysInline() void windowCentre(const IterationParameters<Det>& parameters, LayerID layer, const LayerCrossingPrediction& prediction, Scalar& phi, Scalar& second)
{
  Scalar point[3]{};
  crossingPoint(parameters, layer, prediction, point);
  const Scalar x = point[0] - parameters.beamX;
  const Scalar y = point[1] - parameters.beamY;
  phi = azimuth(x, y);
  second = secondCoordinate<Det>(layer, o2::gpu::CAMath::Hypot(x, y), point[2]);
}

/// What a follow's observer sees (`FollowStage::follow`); the default sees
/// nothing and compiles to nothing:
///  - each direction's start, and each step: the layer it leaves, its state,
///    its covariance's factors (the state's own `cov` is not kept up to
///    date) and the rows the direction has taken;
///  - each target the step searches, with its window;
///  - each row the search visits, with its outcome (`visit`) and the update's
///    chi2 (0 where the update was not made);
///  - a ranking by confirmation (the lookahead), with the rows passing;
///  - each hole test, with whether the track crosses a sensor;
///  - the kink pass, with the kinks left; each row it visits, with its
///    outcome under the scaled step; each row passing it, with what confirms
///    it: the next step's row (`ClusterList::kNoCluster` for the beam line)
///    and its chi2, or nothing;
///  - the row a step takes, its chi2 and whether a kink step took it;
///  - the last row dropped as nothing confirms it (the unconfirmed drop), refused at the
///    layout's end (EndDensity) or against the beam line (BeamEnd), and the
///    end.
struct NoObserver {
  GPUhdi() GPUAlwaysInline() void begin(PropagationDirection /*direction*/) const {}
  GPUhdi() GPUAlwaysInline() void step(LayerID /*layer*/, const BoundTrackState& /*state*/, const CovarianceFactors& /*factors*/, std::uint32_t /*rows*/) const {}
  GPUhdi() GPUAlwaysInline() void target(LayerID /*layer*/, const CandidateWindow& /*window*/) const {}
  GPUhdi() GPUAlwaysInline() void candidate(LayerID /*layer*/, std::uint32_t /*row*/, StepOutcome /*outcome*/, Scalar /*chi2*/) const {}
  GPUhdi() GPUAlwaysInline() void lookahead(std::uint32_t /*passing*/) const {}
  GPUhdi() GPUAlwaysInline() void hole(LayerID /*layer*/, bool /*crossed*/) const {}
  GPUhdi() GPUAlwaysInline() void kink(std::uint32_t /*kinksLeft*/) const {}
  GPUhdi() GPUAlwaysInline() void scaled(LayerID /*layer*/, std::uint32_t /*row*/, StepOutcome /*outcome*/, Scalar /*chi2*/) const {}
  GPUhdi() GPUAlwaysInline() void confirmation(std::uint32_t /*row*/, bool /*confirmed*/, std::uint32_t /*next*/, Scalar /*chi2*/) const {}
  /// the lookahead's confirmation of a nominal row (`confirm`): the row, its chi2, and
  /// the next step's row confirming it and that row's chi2 (kNoCluster: none).
  GPUhdi() GPUAlwaysInline() void lookaheadConfirmation(std::uint32_t /*row*/, Scalar /*chi2*/, std::uint32_t /*next*/, Scalar /*nextChi2*/) const {}
  GPUhdi() GPUAlwaysInline() void taken(std::uint32_t /*row*/, Scalar /*chi2*/, bool /*kinked*/) const {}
  GPUhdi() GPUAlwaysInline() void unconfirmed(std::uint32_t /*row*/) const {}
  /// The last row refused at the layout's end (EndDensity's layout-end refusal).
  GPUhdi() GPUAlwaysInline() void refusedEnd(std::uint32_t /*row*/) const {}
  /// The inward last row refused against the beam line (BeamEnd).
  GPUhdi() GPUAlwaysInline() void refusedBeamEnd(std::uint32_t /*row*/) const {}
  GPUhdi() GPUAlwaysInline() void end(FollowEnd /*why*/) const {}
};

/// What a step's search reads: the iteration, the frame and its grid, the
/// sensors, the material, the field and the used mask.
template <DetectorLayout Det, typename Material>
struct StepContext {
  const IterationParameters<Det>* parameters = nullptr;
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  const Material* material = nullptr;
  const Propagator* propagator = nullptr;
  const std::uint8_t* used = nullptr;
};

/// A step's row: the row, its update's chi2, and the state updated with it.
struct Child {
  std::uint32_t row = std::numeric_limits<std::uint32_t>::max();
  Scalar chi2 = 0.F;
  BoundTrackState state;
  CovarianceFactors factors;
};

/// A step's row as EndDensity reads it: with its update's innovation covariance's determinant
/// [cm^4]. A type of its own, so that the follow without it keeps `Child`.
struct DensityChild : Child {
  Scalar determinant = 0.F;
};

/// The follow's row type with or without \a EndDensity.
template <bool EndDensity>
using ChildOf = std::conditional_t<EndDensity, DensityChild, Child>;

/// Where \a child's update writes its determinant (`updateCandidate`):
/// nowhere for a `Child`.
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar* determinantOf(Child& /*child*/) { return nullptr; }
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar* determinantOf(DensityChild& child) { return &child.determinant; }

/// One search of a step: each row of \a layer the grid holds in \a window, in
/// the grid's order, through the step's tests in the order `StepOutcome`
/// names them: free (with \a Own, \a own counts as free: the track's own end
/// decided again), inside the window's spatial bounds, its segment from
/// (\a r, \a z) -- where the track stands -- crossing none of \a holes (0: no
/// test, LinkHoles.h), and the update of \a state (its covariance \a factors)
/// with it under the gate: the nominal update, or with \a Scaled the step's
/// scattering scaled by \a noiseScale (`updateCandidate`). With
/// \a EndDensity a candidate holds its update's determinant.
/// \a report(row, outcome, chi2) sees every row; \a counter counts the
/// updates made.
///
/// The lanes of \a group walk the window in step
/// (`MeasurementGrid::visitCandidates`), each its share of the rows, and
/// every lane calls \a step(candidate) at every step of the walk:
/// \a candidate its row passing the gate, with the state updated with it and
/// its factors, or none there (its row `kNoCluster`). So \a step may call
/// the group's collectives. One lane has nothing to keep in step: it calls
/// \a step for each row passing the gate alone, right after its report, in
/// the grid's order. \a Inlined: the transport of each row inlined
/// (`updateCandidate`). The body, inlined where it is called; the follow
/// calls `visit`.
template <bool Scaled, bool Own = false, bool EndDensity = false, bool Inlined = false, DetectorLayout Det, typename Material, typename Group, typename Report, typename Step>
GPUhdi() GPUAlwaysInline() void visitInline(const StepContext<Det, Material>& context, const BoundTrackState& state, const CovarianceFactors& factors, Scalar r, Scalar z, LayerID layer, PropagationDirection direction, const CandidateWindow& window, std::uint64_t holes, Scalar noiseScale, WorkCounter counter, const Group& group, const Report& report, const Step& step, [[maybe_unused]] std::uint32_t own = std::numeric_limits<std::uint32_t>::max())
{
  const auto& parameters = *context.parameters;
  const auto& frame = *context.frame;
  const FollowGridScope followCounts;
  // The tests before the update, each reported where it refuses the row.
  const auto reaches = [&](std::uint32_t row) {
    if (context.used[row] != kFree && (!Own || row != own)) {
      report(row, StepOutcome::kUsed, 0.F);
      return false;
    }
    const auto measurement = frame.measurement(row);
    if (!spatiallyCompatible(window, measurement.phi, secondCoordinate<Det>(layer, measurement.r, measurement.z))) {
      report(row, StepOutcome::kOutsideWindow, 0.F);
      return false;
    }
    if (holes != 0U && segmentLeavesHole(holes, parameters.barrelLayers, parameters.layerSpans, r, z, measurement.r, measurement.z)) {
      report(row, StepOutcome::kHole, 0.F);
      return false;
    }
    return true;
  };
  // The update into \a candidate, its row set where it passes the gate.
  const auto update = [&](std::uint32_t row, ChildOf<EndDensity>& candidate) {
    countWork(counter);
    Scalar chi2 = 0.F;
    const bool usable = updateCandidate<Scaled, EndDensity, Inlined>(*context.propagator, *context.geometry, *context.material, measurementAt(frame, row), direction, state, factors, candidate.state, candidate.factors, chi2, noiseScale, determinantOf(candidate));
    if (!usable) {
      report(row, StepOutcome::kNoUpdate, 0.F);
      return;
    }
    const bool passes = chi2 <= parameters.clusterChi2Max;
    report(row, passes ? StepOutcome::kPasses : StepOutcome::kChi2, chi2);
    if (passes) {
      candidate.row = row;
      candidate.chi2 = chi2;
    }
  };
  if constexpr (Group::kLanes == 1U) {
    // One lane: the grid's walk of one thread, the candidate of the row
    // alone, made where it is updated.
    frame.grid().visitCandidates(layer, window, frame.grid().geometry().allBlocks(), [&](std::uint32_t row) {
      if (!reaches(row)) {
        return;
      }
      ChildOf<EndDensity> candidate;
      update(row, candidate);
      if (candidate.row != std::numeric_limits<std::uint32_t>::max()) {
        step(candidate);
      }
    });
  } else {
    // A group: one candidate a lane for the window, its row reset at every
    // step; the state and its factors are written by the update before they
    // are read.
    ChildOf<EndDensity> candidate;
    frame.grid().visitCandidates(
      layer, window, frame.grid().geometry().allBlocks(),
      [&](std::uint32_t row, bool has) {
        candidate.row = std::numeric_limits<std::uint32_t>::max();
        if (has && reaches(row)) {
          update(row, candidate);
        }
        step(candidate);
      },
      group);
  }
}

/// `visitInline`, never inlined on the device: one copy for each call site,
/// which the kernels then do not absorb.
template <bool Scaled, bool Own = false, bool EndDensity = false, DetectorLayout Det, typename Material, typename Group, typename Report, typename Step>
GPUhdnii() void visit(const StepContext<Det, Material>& context, const BoundTrackState& state, const CovarianceFactors& factors, Scalar r, Scalar z, LayerID layer, PropagationDirection direction, const CandidateWindow& window, std::uint64_t holes, Scalar noiseScale, WorkCounter counter, const Group& group, const Report& report, const Step& step, [[maybe_unused]] std::uint32_t own = std::numeric_limits<std::uint32_t>::max())
{
  visitInline<Scaled, Own, EndDensity>(context, state, factors, r, z, layer, direction, window, holes, noiseScale, counter, group, report, step, own);
}

/// \a candidate of lane \a from of \a group, on every lane: its row, chi2,
/// parameters, surface and factors (and with \a EndDensity its determinant).
/// Its covariance is \a origin's -- the state
/// the step's updates started from: `transportCandidate` and the update carry
/// the input's `cov` unchanged (the factors are the covariance), so every
/// candidate of a step holds the same one and it is not sent.
template <typename Group, typename Row>
[[nodiscard]] GPUhdi() GPUAlwaysInline() decltype(auto) shared(const Group& group, const Row& candidate, std::uint32_t from, const BoundTrackState& origin)
{
  if constexpr (Group::kLanes == 1U) {
    // The candidate itself: one lane has it already.
    static_cast<void>(group);
    static_cast<void>(from);
    static_cast<void>(origin);
    return (candidate);
  } else {
    struct Wire {
      std::uint32_t row;
      Scalar chi2;
      Scalar par[kNBoundParams];
      SensorID surface;
      Scalar surfacePhi;
      CovarianceFactors factors;
    };
    Wire wire{
      .row = candidate.row,
      .chi2 = candidate.chi2,
      .par = {},
      .surface = candidate.state.surface,
      .surfacePhi = candidate.state.surfacePhi,
      .factors = candidate.factors,
    };
    for (int i = 0; i < kNBoundParams; ++i) {
      wire.par[i] = candidate.state.par[i];
    }
    wire = group.broadcast(wire, from);
    Row out;
    if constexpr (std::is_same_v<Row, DensityChild>) {
      out.determinant = group.broadcast(candidate.determinant, from);
    }
    out.row = wire.row;
    out.chi2 = wire.chi2;
    for (int i = 0; i < kNBoundParams; ++i) {
      out.state.par[i] = wire.par[i];
    }
    for (int i = 0; i < kBoundCovSize; ++i) {
      out.state.cov[i] = origin.cov[i];
    }
    out.state.surface = wire.surface;
    out.state.surfacePhi = wire.surfacePhi;
    out.factors = wire.factors;
    return out;
  }
}

/// What `FollowStage::passingOn` sees on a layer: the rows passing the nominal gate, whether the row asked
/// about is among them and its chi2, and the smallest chi2 of the others.
struct Passing {
  std::uint32_t count = 0U;
  bool own = false;
  Scalar ownChi2 = 0.F;
  Scalar otherChi2 = std::numeric_limits<Scalar>::max();
};

/// The witness of a confirmation: the row the next step takes, by the
/// nominal pass's choice -- the smallest chi2, the row breaking ties -- or
/// none (`kNoCluster`).
struct Witness {
  std::uint32_t row = std::numeric_limits<std::uint32_t>::max();
  Scalar chi2 = 0.F;
};

/// Whether \a row with \a chi2 beats \a best: the nominal pass's choice, as
/// the one-lane walk makes it row by row. Over the usable, gated chi2 (no
/// NaN) a strict total order on the rows, so the choice does not depend on
/// the order the rows come in.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool beats(std::uint32_t row, Scalar chi2, const Witness& best) { return best.row == std::numeric_limits<std::uint32_t>::max() || chi2 < best.chi2 || (!(best.chi2 < chi2) && row < best.row); }

/// The best of the lanes' \a own witnesses (`beats`), on every lane of
/// \a group.
template <typename Group>
[[nodiscard]] GPUhdi() GPUAlwaysInline() Witness bestOf(const Group& group, const Witness& own)
{
  if constexpr (Group::kLanes == 1U) {
    return own;
  } else {
    Witness best;
    for (auto lanes = group.ballot(own.row != std::numeric_limits<std::uint32_t>::max()); lanes != 0U; lanes &= lanes - 1U) {
      const auto other = group.broadcast(own, o2::gpu::CAMath::Ctz64(lanes));
      if (beats(other.row, other.chi2, best)) {
        best = other;
      }
    }
    return best;
  }
}

/// The most links leaving or entering one layer of \a Det: the targets a
/// step orders.
template <DetectorLayout Det>
inline constexpr std::size_t maxDegree = [] {
  std::size_t most = 1U;
  for (std::size_t layer = 0; layer < static_cast<std::size_t>(Detector<Det>::nLayers); ++layer) {
    most = std::max(most, outLinks<Det>.count(static_cast<LayerID>(layer)));
    most = std::max(most, inLinks<Det>.count(static_cast<LayerID>(layer)));
  }
  return most;
}();

/// One direction of a track's follow: where it stands, its covariance's
/// factors, the rows it took, the last one's chi2 (the unconfirmed drop removes it), the holes it
/// crossed, whether its latest step, one that took no row, crossed a
/// sensor without a row, and the track's own row its searches take as free
/// (the end decided again, `followEnd`; none otherwise).
template <std::size_t K>
struct Branch {
  BoundTrackState state;
  CovarianceFactors factors;
  std::uint32_t rows[K]{};
  Scalar lastChi2 = 0.F;
  std::uint32_t count = 0U;
  std::uint32_t holes = 0U;
  std::uint32_t own = std::numeric_limits<std::uint32_t>::max();
  bool crossed = false;
};

/// A direction as EndDensity reads it: with its last row's update's
/// determinant [cm^4], and whether a kink step took it -- its update's V the
/// scaled one, its confirmation the witness. A type of its own, so that the
/// follow without it keeps `Branch`.
template <std::size_t K>
struct DensityBranch : Branch<K> {
  Scalar lastDeterminant = 0.F;
  bool lastKink = false;
};

/// The follow's direction type with or without \a EndDensity.
template <std::size_t K, bool EndDensity>
using BranchOf = std::conditional_t<EndDensity, DensityBranch<K>, Branch<K>>;

/// What BeamEnd reads of an inward walk: the states its last two steps started from -- entry 0 the
/// latest -- their factors, whether a kink step took the row, and how many
/// rows the walk took. A row the unconfirmed drop or the layout's end removes makes the one
/// before the last: entry 1.
struct BeamEndHistory {
  BoundTrackState state[2];
  CovarianceFactors factors[2];
  bool kink[2] = {true, true};
  std::uint32_t taken = 0U;

  /// A step from \a from (its factors \a fromFactors) taking a row, by a
  /// kink step where \a kinked.
  GPUhdi() GPUAlwaysInline() void remember(const BoundTrackState& from, const CovarianceFactors& fromFactors, bool kinked)
  {
    state[1] = state[0];
    factors[1] = factors[0];
    kink[1] = kink[0];
    state[0] = from;
    factors[0] = fromFactors;
    kink[0] = kinked;
    ++taken;
  }
};

/// The follow without BeamEnd keeps nothing.
struct NoBeamEndHistory {
};

/// The inward walk's history with or without \a BeamEnd.
template <bool BeamEnd>
using BeamEndHistoryOf = std::conditional_t<BeamEnd, BeamEndHistory, NoBeamEndHistory>;

/// A layer a step may reach: its nominal crossing and its distance from the
/// track, squared.
struct Target {
  Scalar distance2 = 0.F;
  LayerID layer = 0;
  LayerCrossingPrediction prediction;
};

} // namespace follow

/// One track's follow: outward from its outward
/// state, then inward from its inner state, into its proposal. It reads the
/// grid and the used mask, and writes only its own proposal and kink steps.
///
/// A step, in either direction:
///  - **Navigation.** The layers every link of the layout reaches from the
///    one the track stands on (enters it, inward), skip links included, in
///    the order the track crosses them: by the distance of the nominal
///    crossing from the track. The first with a free row within the gate
///    gives the step's row. A layer the state cannot be carried to is
///    passed.
///  - **The window.** The prediction's own (`predictedWindow`) and the
///    track's time box, which every row taken narrows, outward then inward.
///  - **Holes.** A layer without a row is a hole only where the track crosses
///    one of its sensors within its bounds (SensorCells.h); a gap costs
///    nothing. Past the iteration's `followHoles` the direction ends; each
///    direction has its own count.
///  - **The lookahead.** Where more than one row passes, each is taken and the next step
///    made from there: the rows it confirms rank first, by their chi2 plus
///    the confirming row's; the others by their chi2; then the row.
///  - **Kinks.** Only where the nominal pass ends the step, the same walk
///    again with each target's scattering scaled by the iteration's kink
///    allowance, its own hole count, under the kink budget the directions
///    share. A row it admits is taken only if the next step, nominal, has a
///    row from it -- or, past the innermost layer inward with the beam
///    constraint, the beam line -- by its chi2 plus that row's; its witness
///    bit says the next row stands on it (`KeepWitnesses`).
///  - **The unconfirmed drop.** A direction whose last step crossed a sensor without a row
///    drops its last row, which nothing confirms -- unless its chi2 is below
///    the iteration's `followKeepChi2` (0: never), or, inward, unless that step
///    saw a row near: one within its window the update refused above the
///    gate, or one the kink pass admitted and could not confirm (inward the
///    sensor usually held the particle's row, refused; outward the row near
///    is a noise row too often). The row
///    seen is not taken.
///
/// Inward starts from the inner state updated with the beam line where the
/// iteration has the beam constraint (\a BeamLine, its own instantiation),
/// and does not run where the iteration shares innermost rows (\a inward).
/// The directions share one capacity, the longest chain less the track's
/// rows.
///
/// With \a DecidesEnd, the iteration's `followEnd`, the track's outermost row is decided again (`decideEnd`): from
/// the road fit's outward pass without it (`follow::EndFit`, a pass of its
/// own before the stage, \a endStarts), the walk outward, that row free to
/// it. Its first row the end: the walk's rows beyond it are the
/// outward follow. Another: the walk's rows replace the end. None: the end
/// stays and the follow is as without. An iteration without `follow` walks
/// one row, the end alone, and nothing inward.
///
/// With \a Greedy, the iteration's `followGreedy`, a step more than one row passes takes the smallest chi2, the
/// row breaking ties (`follow::beats`), instead of the lookahead's confirmation.
///
/// With \a EndDensity, where the iteration has `followEndMaxSqrtDet` or
/// `followKeepDensity`, a
/// direction's end is decided on its last row's own update, V = H P H^T + R
/// and the chi2, where a nominal step took it:
///  - **inward at the layout's end** -- the layer it stands on has no
///    inward link -- the last row goes where sqrt|V| is above
///    `followEndMaxSqrtDet`: a prediction carried that far admits strangers;
///  - **outward, where the unconfirmed drop would remove it**, the last row stays where its
///    Gaussian density g = exp(-chi2 / 2) / (2 pi sqrt|V|) is above
///    `followKeepDensity`: the keep's chi2 with a penalty for a wide
///    prediction.
///
/// With \a BeamEnd, the iteration's `followBeamEnd` (only with \a BeamLine), the inward walk's last
/// row, where a nominal step took it, is tested against that step's
/// prediction updated with the beam line (`refusesBeamEnd`): refused above
/// the gate. Its own step saw only the rows outside; the beam line is the
/// evidence beyond the last row inward, where the strangers of the vertex
/// detector -- QED pairs in Pb-Pb -- pass the step alone.
template <DetectorLayout Det, typename Material, bool BeamLine = false, bool DecidesEnd = false, bool Greedy = false, bool EndDensity = false, bool BeamEnd = false>
class FollowStage
{
 public:
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  static_assert(kMaxClusters <= 64U, "a kink step's witness is a bit of a 64-bit mask");
  static_assert(BeamLine || !BeamEnd, "the beam-line end needs the iteration's beam constraint");
  /// Whether it decides the track's end again (`DecidesEnd`).
  static constexpr bool kDecidesEnd = DecidesEnd;
  using TrackView = typename BoundedColumns<TrackOf<kMaxClusters>::template RecordT>::ConstView;
  using Proposal = FollowProposal<kMaxClusters>;
  using Branch = follow::BranchOf<kMaxClusters, EndDensity>;
  using Child = follow::ChildOf<EndDensity>;

  FollowStage(const IterationParameters<Det>* parameters, const TimeFrame<Det>* frame, const TrackingGeometry<Det>* geometry, const Material* material, const Propagator& propagator, TrackView tracks, const std::uint8_t* used, const std::uint32_t* maxCandidateWidthPS, SensorCellsView cells, Proposal* proposals, follow::KinkSteps* kinkSteps, bool inward, std::uint32_t trackOffset, const BoundTrackState* endStarts = nullptr)
    : mParameters(parameters),
      mFrame(frame),
      mGeometry(geometry),
      mMaterial(material),
      mPropagator(propagator),
      mTracks(tracks),
      mUsed(used),
      mMaxCandidateWidthPS(maxCandidateWidthPS),
      mCells(cells),
      mProposals(proposals),
      mKinkSteps(kinkSteps),
      mInward(inward),
      mTrackOffset(trackOffset),
      mEndStarts(endStarts)
  {
  }

  GPUhdi() void operator()(std::uint32_t track) const { together(track, SoloGroup{}); }

  /// `operator()` for the lanes of \a group, which follow \a track together:
  /// every lane makes the same decisions, the
  /// first writes the proposal.
  template <typename Group>
  GPUhdi() void together(std::uint32_t track, const Group& group) const
  {
    const std::uint32_t input = mTrackOffset + track;
    Proposal proposal;
    for (std::size_t i = 0; i < kMaxClusters; ++i) {
      proposal.outward[i] = ClusterList<kMaxClusters>::kNoCluster;
      proposal.inward[i] = ClusterList<kMaxClusters>::kNoCluster;
    }
    follow::KinkSteps steps;
    if (mTracks.status[input] == PropagationStatus::kOk) {
      follow(mTracks.state[input], mTracks.outwardState[input], mTracks.clusters[input], mTracks.count[input], proposal, steps, follow::NoObserver{}, mEndStarts != nullptr ? mEndStarts + track : nullptr, group);
    }
    if (group.lane() == 0U) {
      mProposals[track] = proposal;
      if (mKinkSteps != nullptr) {
        mKinkSteps[track] = steps;
      }
    }
  }

  /// A track's follow from \a inner and \a outer, its inner and outward
  /// states, its \a count rows \a rows giving the time box, into
  /// \a proposal's rows and \a steps' kink steps (both empty before),
  /// \a observe seeing every decision: what `operator()` does for a track,
  /// for the bench's replay too. With \a DecidesEnd, \a endStart is the
  /// track's state without its outermost row (`follow::EndFit`). The lanes of
  /// \a group follow it together.
  template <typename Observer, typename Group = SoloGroup>
  GPUhdi() void follow(const BoundTrackState& inner, const BoundTrackState& outer, const ClusterList<kMaxClusters>& rows, std::uint32_t count, Proposal& proposal, follow::KinkSteps& steps, const Observer& observe, [[maybe_unused]] const BoundTrackState* endStart = nullptr, const Group& group = {}) const
  {
    // A track as long as a chain can be is not followed.
    if (count >= kMaxClusters) {
      return;
    }
    if constexpr (DecidesEnd) {
      if (endStart != nullptr && endStart->surface != kInvalidSensorID && decideEnd(inner, *endStart, rows, count, proposal, steps, observe, group)) {
        return;
      }
      if (!mParameters->follow) {
        return;
      }
    }
    const std::uint32_t capacity = static_cast<std::uint32_t>(kMaxClusters) - count;
    auto box = trackTimeBox(mFrame->measurements(), rows, static_cast<std::uint8_t>(count));
    std::uint32_t kinks = 0U;
    proposal.outwardCount = walk(outer, PropagationDirection::kForward, proposal.outward, capacity, box, kinks, steps.outwardWitnessed, observe, std::numeric_limits<std::uint32_t>::max(), group);
    walkInward(inner, rows.index[0], capacity - proposal.outwardCount, box, kinks, proposal, steps, observe, group);
  }

  /// The free rows of \a layer passing the nominal gate from \a state, its
  /// covariance \a factors, in \a direction within the time \a box: what a
  /// nominal step from there would see on that layer, \a own among them or
  /// not (`follow::TrimAmbiguous`). None where no link of the layout leads from
  /// the state's layer to \a layer, or the window cannot be made.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() follow::Passing passingOn(const BoundTrackState& state, const CovarianceFactors& factors, LayerID layer, std::uint32_t own, const TimeEstPS& box, PropagationDirection direction) const
  {
    follow::Target targets[follow::maxDegree<Det>];
    const std::uint32_t nTargets = targetsOf(state, direction, targets);
    for (std::uint32_t t = 0; t < nTargets; ++t) {
      const auto& target = targets[t];
      if (target.layer != layer) {
        continue;
      }
      CandidateWindow window;
      Scalar reach[2]{};
      if (!targetWindow<false>(state, factors, box, target, direction, 1.F, window, reach)) {
        return {};
      }
      follow::Passing passing;
      // The search inlined with its transports: called, it slows
      // `TrimAmbiguous`'s kernel on a device.
      follow::visitInline<false, false, false, true>(
        context(), state, factors, 0.F, 0.F, target.layer, direction, window, 0U, 1.F, WorkCounter::kFollowTrimUpdates, SoloGroup{}, [](std::uint32_t /*row*/, StepOutcome /*outcome*/, Scalar /*chi2*/) {},
        [&](const follow::Child& candidate) {
          ++passing.count;
          if (candidate.row == own) {
            passing.own = true;
            passing.ownChi2 = candidate.chi2;
          } else {
            passing.otherChi2 = o2::gpu::CAMath::Min(passing.otherChi2, candidate.chi2);
          }
        });
      return passing;
    }
    return {};
  }

  /// The walk inward after the outward one, from \a inner (`inwardStart`)
  /// over \a capacity rows at most, where the stage follows inward and rows
  /// are left: \a kinks the outward walk's so far, then both counts into
  /// \a steps.
  template <typename Observer, typename Group>
  GPUhdi() GPUAlwaysInline() void walkInward(const BoundTrackState& inner, std::uint32_t innermostRow, std::uint32_t capacity, TimeEstPS& box, std::uint32_t& kinks, Proposal& proposal, follow::KinkSteps& steps, const Observer& observe, const Group& group) const
  {
    steps.outward = static_cast<std::uint8_t>(kinks);
    if (mInward && capacity != 0U) {
      proposal.inwardCount = walk(inwardStart(inner, innermostRow), PropagationDirection::kBackward, proposal.inward, capacity, box, kinks, steps.inwardWitnessed, observe, std::numeric_limits<std::uint32_t>::max(), group);
    }
    steps.inward = static_cast<std::uint8_t>(kinks - steps.outward);
  }

  /// The track's outermost row decided again (\a DecidesEnd): from \a start,
  /// \a rows' first \a count - 1 fitted outward (`follow::EndFit`), the
  /// walk outward with the end free to it -- the follow's capacity, or one
  /// row without `follow`. Whether the walk took a row: then \a proposal and
  /// \a steps hold the follow, inward too with `follow`; else they are
  /// untouched.
  template <typename Observer, typename Group>
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool decideEnd(const BoundTrackState& inner, const BoundTrackState& start, const ClusterList<kMaxClusters>& rows, std::uint32_t count, Proposal& proposal, follow::KinkSteps& steps, const Observer& observe, const Group& group) const
  {
    const std::uint32_t end = rows.index[count - 1U];
    const std::uint32_t capacity = mParameters->follow ? static_cast<std::uint32_t>(kMaxClusters) - (count - 1U) : 1U;
    auto box = trackTimeBox(mFrame->measurements(), rows, static_cast<std::uint8_t>(count - 1U));
    std::uint32_t kinks = 0U;
    std::uint64_t witnessed = 0U;
    std::uint32_t walked[kMaxClusters];
    const std::uint8_t taken = walk(start, PropagationDirection::kForward, walked, capacity, box, kinks, witnessed, observe, end, group);
    if (taken == 0U) {
      return false;
    }
    // The end taken again: the rows beyond it, and their witnesses, are the
    // follow outward; another row: the walk's rows replace the end.
    const bool replaced = walked[0] != end;
    const std::uint32_t first = replaced ? 0U : 1U;
    if (replaced) {
      countWork(WorkCounter::kFollowEndReplaced);
      proposal.replacedEnd = end;
    }
    for (std::uint32_t i = first; i < taken; ++i) {
      proposal.outward[i - first] = walked[i];
    }
    proposal.outwardCount = static_cast<std::uint8_t>(taken - first);
    steps.outwardWitnessed = witnessed >> first;
    const std::uint32_t held = count - (replaced ? 1U : 0U) + proposal.outwardCount;
    walkInward(inner, rows.index[0], mParameters->follow && held < kMaxClusters ? static_cast<std::uint32_t>(kMaxClusters) - held : 0U, box, kinks, proposal, steps, observe, group);
    return true;
  }

  /// Where the inward follow starts: \a inner, with \a BeamLine updated with
  /// the beam line, \a innermost the track's innermost row -- the nominal
  /// update, or, where the iteration has a kink allowance and the nominal
  /// update is above the gate, the update with the scattering inside scaled
  /// by it, as `clearsBeamLine` admitted the road.
  [[nodiscard]] GPUhdi() BoundTrackState inwardStart(const BoundTrackState& inner, std::uint32_t innermost) const
  {
    if constexpr (BeamLine) {
      BoundTrackState constrained = inner;
      const auto& plane = mGeometry->planeOf(inner);
      const LayerID layer = mGeometry->surface(inner.surface).layer;
      const auto point = pairPointAt(mFrame->measurements(), innermost);
      const auto update = updateWithBeamLine(mPropagator, *mParameters, plane, layer, point, constrained);
      if (update.ok && !(update.chi2 <= mParameters->beamChi2Max) && mParameters->kinkScale2 > 1.F) {
        constrained = inner;
        static_cast<void>(updateWithBeamLine(mPropagator, *mParameters, plane, layer, point, constrained, mParameters->kinkScale2));
      }
      return constrained;
    } else {
      static_cast<void>(innermost);
      return inner;
    }
  }

  /// One direction of a follow from \a start: at most \a capacity rows into
  /// \a rows, \a box the track's time box, \a kinks the kink steps the
  /// follow took so far, \a witnessed a bit per row standing on the next
  /// one, \a own a row of the track its searches take as free (with
  /// \a DecidesEnd); how many rows. The lanes of \a group walk it together.
  template <typename Observer, typename Group = SoloGroup>
  GPUhdi() GPUAlwaysInline() std::uint8_t walk(const BoundTrackState& start, PropagationDirection direction, std::uint32_t* rows, std::uint32_t capacity, TimeEstPS& box, std::uint32_t& kinks, std::uint64_t& witnessed, const Observer& observe, std::uint32_t own = std::numeric_limits<std::uint32_t>::max(), const Group& group = {}) const
  {
    observe.begin(direction);
    Branch branch;
    branch.own = own;
    branch.state = start;
    branch.factors = factorCovariance(start.cov);
    auto end = FollowEnd::kCapacity;
    // Whether the latest step saw a row near where it took none (the unconfirmed drop, below).
    bool nearMiss = false;
    [[maybe_unused]] follow::BeamEndHistoryOf<BeamEnd> history; // NOLINT(misc-const-correctness): written in a dependent call the check does not see
    while (branch.count < capacity) {
      const LayerID current = mGeometry->surface(branch.state.surface).layer;
      observe.step(current, branch.state, branch.factors, branch.count);
      follow::Target targets[follow::maxDegree<Det>];
      const std::uint32_t nTargets = targetsOf(branch.state, direction, targets);
      if (nTargets == 0U) {
        end = FollowEnd::kNoTarget;
        break;
      }
      const std::uint32_t holesBefore = branch.holes;
      branch.crossed = false;
      nearMiss = false;
      Child child;
      if (nominal(branch, box, direction, targets, nTargets, child, nearMiss, observe, group)) {
        observe.taken(child.row, child.chi2, false);
        if constexpr (BeamEnd) {
          if (direction == PropagationDirection::kBackward) {
            history.remember(branch.state, branch.factors, false);
          }
        }
        take(branch, box, child);
        continue;
      }
      end = branch.holes > mParameters->followHoles ? FollowEnd::kHoles : FollowEnd::kNoRow;
      // The kink pass: no target had a row within the holes.
      if (mParameters->kinkScale2 > 1.F && kinks < mParameters->maxKinks) {
        observe.kink(mParameters->maxKinks - kinks);
        countWork(WorkCounter::kFollowKinkSteps);
        std::uint32_t kinkHoles = holesBefore;
        bool onNext = false;
        if (kinkStep(branch, box, direction, targets, nTargets, child, kinkHoles, onNext, nearMiss, observe, group)) {
          ++kinks;
          countWork(WorkCounter::kFollowKinksTaken);
          witnessed |= onNext ? std::uint64_t{1} << branch.count : std::uint64_t{0};
          observe.taken(child.row, child.chi2, true);
          if constexpr (BeamEnd) {
            if (direction == PropagationDirection::kBackward) {
              history.remember(branch.state, branch.factors, true);
            }
          }
          take(branch, box, child);
          if constexpr (EndDensity) {
            branch.lastKink = true;
          }
          branch.holes = kinkHoles;
          continue;
        }
      }
      break;
    }
    // The unconfirmed drop: the last step crossed a sensor
    // without a row; the last row, which nothing confirms, goes -- unless,
    // inward, the step saw a row near (inward that row is almost always the
    // particle's own, outward a noise row often enough to add fakes, so
    // the outward unconfirmed drop stands).
    // A last row under the iteration's followKeepChi2 stays (the outward unconfirmed drop
    // drops mostly the particle's own rows, the wrong ones at the higher
    // chi2).
    const bool evidence = (direction == PropagationDirection::kBackward && nearMiss) || branch.lastChi2 < mParameters->followKeepChi2;
    // With EndDensity also an
    // outward last row of a high density stays (`densityKept`).
    if (branch.crossed && !evidence && branch.count != 0U && !densityKept(branch, direction)) {
      --branch.count;
      countWork(WorkCounter::kFollowUnconfirmed);
      observe.unconfirmed(branch.rows[branch.count]);
    }
    if constexpr (EndDensity) {
      // At every exit: a row taken as the last the capacity allows leaves
      // the loop before the search for targets. Never where the unconfirmed drop removed a
      // row: a layer without an inward link has no targets, so the walk
      // crossed nothing there. A refused row a kink row stood on takes that
      // row with it (`KeepWitnesses`).
      if (refusesLayoutEnd(branch, direction)) {
        --branch.count;
        countWork(WorkCounter::kFollowEndRefused);
        observe.refusedEnd(branch.rows[branch.count]);
      }
    }
    if constexpr (BeamEnd) {
      // After the unconfirmed drop and the layout's end, on the row that stands last; a
      // refused row a kink row stood on takes that row with it
      // (`KeepWitnesses`).
      if (refusesBeamEnd(branch, direction, history)) {
        --branch.count;
        countWork(WorkCounter::kFollowBeamEndRefused);
        observe.refusedBeamEnd(branch.rows[branch.count]);
      }
    }
    observe.end(end);
    for (std::uint32_t i = 0; i < branch.count; ++i) {
      rows[i] = branch.rows[i];
    }
    return static_cast<std::uint8_t>(branch.count);
  }

 private:
  /// EndDensity's density keep, asked where the unconfirmed drop would remove \a branch's last row:
  /// whether, outward and taken by a nominal step, it has its density g
  /// above the iteration's `followKeepDensity` -- in logs, chi2 + log|V|
  /// below -2 log(2 pi key). A determinant that is not positive never keeps.
  /// Without EndDensity, false.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool densityKept(const Branch& branch, PropagationDirection direction) const
  {
    if constexpr (EndDensity) {
      using Math = o2::gpu::CAMath;
      const Scalar key = mParameters->followKeepDensity;
      const bool kept = direction == PropagationDirection::kForward && key > 0.F && !branch.lastKink && branch.lastDeterminant > 0.F && branch.lastChi2 + Math::Log(branch.lastDeterminant) < -2.F * Math::Log(2.F * std::numbers::pi_v<Scalar> * key);
      countWork(WorkCounter::kFollowDensityKept, kept ? 1U : 0U);
      return kept;
    } else {
      static_cast<void>(branch);
      static_cast<void>(direction);
      return false;
    }
  }

  /// EndDensity's layout-end refusal: whether \a branch, inward, stands at the layout's
  /// end -- no inward link from its layer, not a search whose every
  /// prediction failed -- on a row a nominal step took with sqrt|V| above the
  /// iteration's `followEndMaxSqrtDet`.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool refusesLayoutEnd(const Branch& branch, PropagationDirection direction) const
  {
    const Scalar key = mParameters->followEndMaxSqrtDet;
    if (direction != PropagationDirection::kBackward || !(key > 0.F) || branch.count == 0U || branch.lastKink) {
      return false;
    }
    const LayerID current = mGeometry->surface(branch.state.surface).layer;
    return mParameters->allInLinkOffsets[current] == mParameters->allInLinkOffsets[current + 1U] && branch.lastDeterminant > key * key;
  }

  /// BeamEnd's rule: whether \a branch, inward, ends on a row a nominal step
  /// took (from \a history's state) whose chi2 against that step's
  /// prediction updated with the beam line is above the gate. The
  /// prediction carried as the step carried it, the update at the row's
  /// sensor from the row's pair point (`updateWithBeamLine`), the residual's
  /// covariance in double as the update's. A prediction or an update that
  /// fails, or a determinant that is not positive, refuses nothing.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool refusesBeamEnd(const Branch& branch, PropagationDirection direction, const follow::BeamEndHistory& history) const
  {
    if (direction != PropagationDirection::kBackward || branch.count == 0U) {
      return false;
    }
    // The rows the unconfirmed drop and the layout's end dropped: at most one of the two acts.
    const std::uint32_t dropped = history.taken - branch.count;
    if (dropped > 1U || history.kink[dropped]) {
      return false;
    }
    const std::uint32_t row = branch.rows[branch.count - 1U];
    const auto measurement = measurementAt(*mFrame, row);
    const auto& plane = mGeometry->planeOf(measurement);
    BoundTrackState predicted;
    CovarianceFactors carried;
    if (!follow::transportCandidate<false>(mPropagator, *mGeometry, *mMaterial, plane, direction, history.state[dropped], history.factors[dropped], predicted, carried, 1.F)) {
      return false;
    }
    covarianceFromFactors(carried, predicted.cov);
    const LayerID layer = mGeometry->surface(measurement.surface).layer;
    if (!updateWithBeamLine(mPropagator, *mParameters, plane, layer, pairPointAt(mFrame->measurements(), row), predicted).ok) {
      return false;
    }
    const auto v00 = static_cast<ExactScalar>(predicted.variance(kLoc0)) + static_cast<ExactScalar>(measurement.varLoc0);
    const auto v11 = static_cast<ExactScalar>(predicted.variance(kLoc1)) + static_cast<ExactScalar>(measurement.varLoc1);
    const auto v01 = static_cast<ExactScalar>(predicted.covAt(kLoc0, kLoc1)) + static_cast<ExactScalar>(measurement.covLoc01);
    const ExactScalar det = (v00 * v11) - (v01 * v01);
    const auto r0 = static_cast<ExactScalar>(measurement.loc0) - static_cast<ExactScalar>(predicted.loc0());
    const auto r1 = static_cast<ExactScalar>(measurement.loc1) - static_cast<ExactScalar>(predicted.loc1());
    return det > 0. && ((r0 * r0 * v11) - (2. * r0 * r1 * v01) + (r1 * r1 * v00)) / det > static_cast<ExactScalar>(mParameters->clusterChi2Max);
  }

  /// What a step's search reads (`follow::visit`).
  [[nodiscard]] GPUhdi() GPUAlwaysInline() follow::StepContext<Det, Material> context() const { return {.parameters = mParameters, .frame = mFrame, .geometry = mGeometry, .material = mMaterial, .propagator = &mPropagator, .used = mUsed}; }

  /// The window of a step from \a state (its covariance \a factors) to
  /// \a target: the track's time \a box, narrowed to the prediction
  /// (`follow::predictedWindow`, with \a Scaled the step's scattering scaled
  /// by \a scale), its \a reach beside. False where the prediction fails.
  template <bool Scaled>
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool targetWindow(const BoundTrackState& state, const CovarianceFactors& factors, const TimeEstPS& box, const follow::Target& target, PropagationDirection direction, Scalar scale, CandidateWindow& window, Scalar (&reach)[2]) const
  {
    Scalar phiCentre = 0.F;
    Scalar secondCentre = 0.F;
    follow::windowCentre(*mParameters, target.layer, target.prediction, phiCentre, secondCentre);
    window = {
      .earliestPS = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(mMaxCandidateWidthPS[target.layer]) + 1,
      .latestPS = static_cast<std::int64_t>(box.upper()) - 1,
      .afterPS = box.lower(),
      .phiCentre = phiCentre,
      .phiHalfWidth = std::numbers::pi_v<Scalar>,
      .secondCentre = secondCentre,
      .secondHalfWidth = follow::kOpenBound,
    };
    return follow::predictedWindow<Scaled>(mPropagator, *mParameters, *mGeometry, *mMaterial, state, factors, target.layer, direction, target.prediction, window, scale, reach);
  }

  /// The layers a step from \a state may reach in \a direction, into
  /// \a targets, nearest crossing first; how many.
  [[nodiscard]] GPUhdi() std::uint32_t targetsOf(const BoundTrackState& state, PropagationDirection direction, follow::Target* targets) const
  {
    const bool forward = direction == PropagationDirection::kForward;
    const LayerID current = mGeometry->surface(state.surface).layer;
    const auto& from = mGeometry->planeOf(state);
    Scalar here[3]{};
    surfaceToGlobal(from, state.loc0(), state.loc1(), here);
    const auto* offsets = forward ? mParameters->allOutLinkOffsets : mParameters->allInLinkOffsets;
    const auto* indices = forward ? mParameters->allOutLinkIndices : mParameters->allInLinkIndices;
    std::uint32_t nTargets = 0U;
    for (auto at = offsets[current]; at < offsets[current + 1U]; ++at) {
      const auto link = indices[at];
      const LayerID layer = forward ? mParameters->linkTargets[link] : mParameters->linkSources[link];
      const auto prediction = predictLayerCrossing(mPropagator, *mParameters, from, state, layer, direction);
      if (prediction.status != PropagationStatus::kOk) {
        continue;
      }
      Scalar point[3]{};
      follow::crossingPoint(*mParameters, layer, prediction, point);
      const Scalar dx = point[0] - here[0];
      const Scalar dy = point[1] - here[1];
      const Scalar dz = point[2] - here[2];
      const follow::Target target{.distance2 = (dx * dx) + (dy * dy) + (dz * dz), .layer = layer, .prediction = prediction};
      // With one link there is nothing to sort; compiling the loop out also
      // spares GCC's -Warray-bounds a subscript it cannot prove unreachable.
      // The link tables are `outLinks<Det>`/`inLinks<Det>`'s (Iteration.h), so
      // a layer's links never outnumber `maxDegree<Det>`, the targets' size;
      // the analyzer, not seeing the tables' origin, takes them as any.
      // NOLINTBEGIN(clang-analyzer-security.ArrayBound)
      if constexpr (follow::maxDegree<Det> > 1U) {
        std::uint32_t slot = nTargets++;
        while (slot > 0U && (target.distance2 < targets[slot - 1U].distance2 || (!(targets[slot - 1U].distance2 < target.distance2) && target.layer < targets[slot - 1U].layer))) {
          targets[slot] = targets[slot - 1U];
          --slot;
        }
        targets[slot] = target;
      } else {
        targets[nTargets++] = target;
      }
      // NOLINTEND(clang-analyzer-security.ArrayBound)
    }
    return nTargets;
  }

  /// The nominal pass's walk over \a branch's \a nTargets \a targets: per
  /// target its window (`predictedWindow`; a target the state cannot be
  /// carried to is passed, the track taken not to cross it), and
  /// \a search(target, window), whether it found a row there; a target
  /// without one crossed in a sensor is a hole, counted on the branch, past
  /// the iteration's holes the pass ends. Whether a search found a row.
  template <typename Observer, typename Search>
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool overTargets(Branch& branch, const TimeEstPS& box, PropagationDirection direction, const follow::Target* targets, std::uint32_t nTargets, const Observer& observe, const Search& search) const
  {
    for (std::uint32_t t = 0; t < nTargets; ++t) {
      const auto& target = targets[t]; // NOLINT(clang-analyzer-security.ArrayBound): nTargets <= maxDegree<Det> (`targetsOf`)
      CandidateWindow window;
      Scalar reach[2]{};
      if (!targetWindow<false>(branch.state, branch.factors, box, target, direction, 1.F, window, reach)) {
        continue;
      }
      observe.target(target.layer, window);
      if (search(target, window)) {
        return true;
      }
      const bool crossed = crossesSensor(branch.state, target, reach, direction);
      observe.hole(target.layer, crossed);
      branch.crossed = branch.crossed || crossed;
      if (crossed && ++branch.holes > mParameters->followHoles) {
        return false;
      }
    }
    return false;
  }

  /// The nominal pass of a step of \a branch over its \a nTargets
  /// \a targets (`overTargets`): the first target with a row in the gate
  /// gives \a child -- the row, or, with more than one, the lookahead's best confirmed
  /// (`confirm`). Whether it found a row; its updates count as the follow's.
  /// \a nearMiss is set where a target held a row within its window that
  /// the update refused above the gate (the unconfirmed drop reads it).
  ///
  /// The lanes of \a group search a window in step: after each step of the
  /// walk the rows that passed, in lane order -- the grid's order -- are
  /// counted for the whole window; the first stays on its lane until a second
  /// comes, then both are confirmed, and every later one as it comes, each
  /// confirmation by the whole group from the row's state on every lane.
  template <typename Observer, typename Group>
  [[nodiscard]] GPUhdi() bool nominal(Branch& branch, const TimeEstPS& box, PropagationDirection direction, const follow::Target* targets, std::uint32_t nTargets, Child& child, bool& nearMiss, const Observer& observe, const Group& group) const
  {
    return overTargets(branch, box, direction, targets, nTargets, observe, [&](const follow::Target& target, const CandidateWindow& window) {
      std::uint32_t passing = 0U;
      std::uint32_t firstLane = 0U; // \a child holds the first passing row on this lane
      Confirmed confirmed;
      bool gated = false; // this lane saw a row refused above the gate
      follow::visit<false, DecidesEnd, EndDensity>(
        context(), branch.state, branch.factors, 0.F, 0.F, target.layer, direction, window, 0U, 1.F, WorkCounter::kFollowUpdates, group,
        [&](std::uint32_t row, StepOutcome outcome, Scalar chi2) {
          gated = gated || outcome == StepOutcome::kChi2;
          observe.candidate(target.layer, row, outcome, chi2);
        },
        [&](const Child& candidate) {
          for (auto lanes = group.ballot(candidate.row != ClusterList<kMaxClusters>::kNoCluster); lanes != 0U; lanes &= lanes - 1U) {
            const auto from = static_cast<std::uint32_t>(o2::gpu::CAMath::Ctz64(lanes));
            if constexpr (Greedy) {
              // Every lane holds the best so far, each passing row shared as
              // it comes.
              const auto one = follow::shared(group, candidate, from, branch.state);
              if (++passing == 1U || follow::beats(one.row, one.chi2, follow::Witness{.row = child.row, .chi2 = child.chi2})) {
                child = one;
              }
              continue;
            }
            if (++passing == 1U) {
              firstLane = from;
              if (group.lane() == from) {
                child = candidate;
              }
              continue;
            }
            if (passing == 2U) {
              confirm(branch, box, direction, follow::shared(group, child, firstLane, branch.state), confirmed, observe, group);
            }
            confirm(branch, box, direction, follow::shared(group, candidate, from, branch.state), confirmed, observe, group);
          }
        },
        branch.own);
      // Any lane's: a collective, every lane calls it.
      nearMiss = nearMiss || group.ballot(gated) != 0U;
      if (passing == 0U) {
        return false;
      }
      if constexpr (Greedy) {
        return true;
      }
      if (passing == 1U) {
        if constexpr (Group::kLanes > 1U) {
          child = follow::shared(group, child, firstLane, branch.state);
        }
      } else {
        countWork(WorkCounter::kFollowLookaheads);
        observe.lookahead(passing);
        child = confirmed.child;
      }
      return true;
    });
  }

  /// The nominal pass of a step of \a branch as a confirmation sees it
  /// (`overTargets`): the first target with a row in the gate gives
  /// \a witness, its smallest chi2, the row breaking ties. Whether it found
  /// one; its updates count on \a counter. The lanes of \a group search a
  /// window each its own rows, and agree on the best after it.
  template <typename Group>
  [[nodiscard]] GPUhdi() bool witnessOf(Branch& branch, const TimeEstPS& box, PropagationDirection direction, const follow::Target* targets, std::uint32_t nTargets, follow::Witness& witness, WorkCounter counter, const Group& group) const
  {
    return overTargets(branch, box, direction, targets, nTargets, follow::NoObserver{}, [&](const follow::Target& target, const CandidateWindow& window) {
      follow::Witness own;
      follow::visit<false, DecidesEnd>(
        context(), branch.state, branch.factors, 0.F, 0.F, target.layer, direction, window, 0U, 1.F, counter, group, [](std::uint32_t /*row*/, StepOutcome /*outcome*/, Scalar /*chi2*/) {},
        [&](const follow::Child& candidate) {
          if (candidate.row != ClusterList<kMaxClusters>::kNoCluster && follow::beats(candidate.row, candidate.chi2, own)) {
            own = {.row = candidate.row, .chi2 = candidate.chi2};
          }
        },
        branch.own);
      witness = follow::bestOf(group, own);
      return witness.row != ClusterList<kMaxClusters>::kNoCluster;
    });
  }

  /// The lookahead: of the rows of a window passing
  /// the gate, where more than one does, the one the next step confirms, by
  /// its chi2 plus the confirming row's; where none is confirmed, the
  /// smallest chi2; the row breaks ties. `nominal` hands it each passing row
  /// in the window's order from the second on, the first with it.
  struct Confirmed {
    Child child; // its row kNoCluster
    bool confirmed = false;
    Scalar key = 0.F;
  };

  /// \a taken, a row passing the gate from \a branch, against \a best; the
  /// lanes of \a group together. \a observe sees the confirmation.
  template <typename Observer, typename Group>
  GPUhdi() void confirm(const Branch& branch, const TimeEstPS& box, PropagationDirection direction, const Child& taken, Confirmed& best, const Observer& observe, const Group& group) const
  {
    std::uint32_t next = ClusterList<kMaxClusters>::kNoCluster;
    Scalar witness = 0.F;
    const bool confirmed = confirmation<false>(branch, box, direction, taken, next, witness, WorkCounter::kFollowUpdates, group);
    const Scalar key = confirmed ? taken.chi2 + witness : taken.chi2;
    observe.lookaheadConfirmation(taken.row, taken.chi2, next, witness);
    const bool better = best.child.row == ClusterList<kMaxClusters>::kNoCluster || (confirmed != best.confirmed ? confirmed : follow::beats(taken.row, key, follow::Witness{.row = best.child.row, .chi2 = best.key}));
    if (better) {
      best.child = taken;
      best.confirmed = confirmed;
      best.key = key;
    }
  }

  /// The kink pass of a step of \a branch over its \a nTargets \a targets:
  /// the nominal pass's walk with each target's scattering scaled by the
  /// kink allowance, its holes counted on \a holes; the first target with a
  /// confirmed row gives \a child, the smallest chi2 plus its confirmation's,
  /// and \a onNext whether a row of the next step confirmed it (else the
  /// beam line). Whether it found one; \a nearMiss is set where a row passed
  /// the scaled gate, confirmed or not (the unconfirmed drop reads it). The lanes of \a group
  /// search a window in step and confirm its passing rows together, in the
  /// grid's order.
  template <typename Observer, typename Group>
  [[nodiscard]] GPUhdi() bool kinkStep(const Branch& branch, const TimeEstPS& box, PropagationDirection direction, const follow::Target* targets, std::uint32_t nTargets, Child& child, std::uint32_t& holes, bool& onNext, bool& nearMiss, const Observer& observe, const Group& group) const
  {
    const Scalar scale = o2::gpu::CAMath::Sqrt(mParameters->kinkScale2);
    for (std::uint32_t t = 0; t < nTargets; ++t) {
      const auto& target = targets[t];
      CandidateWindow window;
      Scalar reach[2]{};
      if (!targetWindow<true>(branch.state, branch.factors, box, target, direction, scale, window, reach)) {
        continue;
      }
      Scalar bestScore = 0.F;
      child.row = ClusterList<kMaxClusters>::kNoCluster;
      follow::visit<true, DecidesEnd, EndDensity>(
        context(), branch.state, branch.factors, 0.F, 0.F, target.layer, direction, window, 0U, scale, WorkCounter::kFollowKinkUpdates, group, [&](std::uint32_t row, StepOutcome outcome, Scalar chi2) { observe.scaled(target.layer, row, outcome, chi2); },
        [&](const Child& candidate) {
          for (auto lanes = group.ballot(candidate.row != ClusterList<kMaxClusters>::kNoCluster); lanes != 0U; lanes &= lanes - 1U) {
            nearMiss = true;
            const Child taken = follow::shared(group, candidate, static_cast<std::uint32_t>(o2::gpu::CAMath::Ctz64(lanes)), branch.state);
            std::uint32_t next = ClusterList<kMaxClusters>::kNoCluster;
            Scalar witness = 0.F;
            const bool confirmed = confirmation<true>(branch, box, direction, taken, next, witness, WorkCounter::kFollowWitnessUpdates, group);
            observe.confirmation(taken.row, confirmed, next, witness);
            if (!confirmed) {
              continue;
            }
            const Scalar score = taken.chi2 + witness;
            if (follow::beats(taken.row, score, follow::Witness{.row = child.row, .chi2 = bestScore})) {
              child = taken;
              bestScore = score;
              onNext = next != ClusterList<kMaxClusters>::kNoCluster;
            }
          }
        },
        branch.own);
      if (child.row != ClusterList<kMaxClusters>::kNoCluster) {
        return true;
      }
      if (crossesSensor(branch.state, target, reach, direction) && ++holes > mParameters->followHoles) {
        return false;
      }
    }
    return false;
  }

  /// What confirms \a taken, a row passing the gate from \a branch: the
  /// next step's row by the nominal pass from past it, into \a next and
  /// \a chi2, its updates counted on \a counter; with \a BeamLineConfirms (a
  /// kink step's row), where no link leads on inward and the iteration has
  /// the beam constraint, the beam line's update under its gate (\a next
  /// none). Whether anything does: a second kink step does not confirm the
  /// first. The lanes of \a group together.
  template <bool BeamLineConfirms, typename Group>
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool confirmation(const Branch& branch, const TimeEstPS& box, PropagationDirection direction, const Child& taken, std::uint32_t& next, Scalar& chi2, WorkCounter counter, const Group& group) const
  {
    Branch probe = branch;
    auto probeBox = box;
    take(probe, probeBox, taken);
    follow::Target targets[follow::maxDegree<Det>];
    const std::uint32_t nTargets = targetsOf(probe.state, direction, targets);
    if (nTargets == 0U) {
      if constexpr (BeamLine && BeamLineConfirms) {
        if (direction == PropagationDirection::kBackward) {
          BoundTrackState constrained = probe.state;
          covarianceFromFactors(probe.factors, constrained.cov);
          const LayerID layer = mGeometry->surface(constrained.surface).layer;
          const auto beam = updateWithBeamLine(mPropagator, *mParameters, mGeometry->planeOf(constrained), layer, pairPointAt(mFrame->measurements(), probe.rows[probe.count - 1U]), constrained);
          chi2 = beam.chi2;
          return beam.ok && beam.chi2 <= mParameters->beamChi2Max;
        }
      }
      return false;
    }
    follow::Witness witness;
    if (!witnessOf(probe, probeBox, direction, targets, nTargets, witness, counter, group)) {
      return false;
    }
    next = witness.row;
    chi2 = witness.chi2;
    return true;
  }

  /// \a branch past \a child's row; \a box narrowed by it.
  GPUhdi() GPUAlwaysInline() void take(Branch& branch, TimeEstPS& box, const Child& child) const
  {
    branch.rows[branch.count++] = child.row;
    branch.state = child.state;
    branch.factors = child.factors;
    branch.lastChi2 = child.chi2;
    if constexpr (EndDensity) {
      branch.lastDeterminant = child.determinant;
      branch.lastKink = false;
    }
    branch.crossed = false;
    box += timeBoxAt(mFrame->measurements(), child.row);
  }

  /// Whether \a state crosses a sensor of \a target's layer within its
  /// bounds in \a direction: carried to each sensor listed near its nominal
  /// crossing, as far from it as the layer's envelope \a reach lets the
  /// crossing move; a curved sensor is tested on its tangent plane at the
  /// crossing's angle.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool crossesSensor(const BoundTrackState& state, const follow::Target& target, const Scalar* reach, PropagationDirection direction) const
  {
    countWork(WorkCounter::kFollowHoleTests);
    const auto& from = mGeometry->planeOf(state);
    Scalar point[3]{};
    follow::crossingPoint(*mParameters, target.layer, target.prediction, point);
    const auto& propagator = pointPropagator(mPropagator, from, state);
    const ExactScalar curvature = propagator.curvature(static_cast<ExactScalar>(state.qOverP()), 1.);
    const bool crossed = mCells.anySensor(target.layer, target.prediction.phi, reach[0], target.prediction.second, reach[1], [&](std::uint32_t id) {
      const TrackingSurface& surface = mGeometry->surface(static_cast<SensorID>(id));
      TrackingSurface plane = surface;
      if (surface.radius > 0.F) {
        Scalar along0 = 0.F;
        Scalar along1 = 0.F;
        for (int i = 0; i < 3; ++i) {
          along0 += (point[i] - surface.origin[i]) * surface.axis0[i];
          along1 += (point[i] - surface.origin[i]) * surface.axis1[i];
        }
        plane = tangentPlane(surface, o2::gpu::CAMath::ATan2(along1, along0));
      }
      BoundTrackState at = state;
      const auto result = propagate(propagator, at, from, plane, direction);
      return result.status == PropagationStatus::kOk && !follow::beyondHalfTurn(curvature * static_cast<ExactScalar>(result.pathLength));
    });
    countWork(WorkCounter::kFollowHoles, crossed ? 1U : 0U);
    return crossed;
  }

  const IterationParameters<Det>* mParameters = nullptr;
  const TimeFrame<Det>* mFrame = nullptr;
  const TrackingGeometry<Det>* mGeometry = nullptr;
  const Material* mMaterial = nullptr;
  Propagator mPropagator;
  TrackView mTracks;
  const std::uint8_t* mUsed = nullptr;
  const std::uint32_t* mMaxCandidateWidthPS = nullptr;
  SensorCellsView mCells;
  Proposal* mProposals = nullptr;
  follow::KinkSteps* mKinkSteps = nullptr;
  bool mInward = true;
  std::uint32_t mTrackOffset = 0U;
  const BoundTrackState* mEndStarts = nullptr;
};

/// Whether the iteration drops tracks after the follow: it has a length
/// after the follow, or it follows and has a pT floor.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool dropsAfterFollow(const IterationParameters<Det>& parameters)
{
  bool floors = false;
  for (const auto floor : parameters.minPtByClusters) {
    floors = floors || floor > 0.F;
  }
  return parameters.minClustersAfterFollow != 0U || (parameters.follow && floors);
}

/// Whether a track of \a count clusters, \a state its refit as the follow
/// left it, stays after the follow: at least the iteration's length after
/// it, and, where the iteration follows, at or above its pT floor for that
/// count.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool keptAfterFollow(const IterationParameters<Det>& parameters, std::uint32_t count, const BoundTrackState& state)
{
  return count >= parameters.minClustersAfterFollow && (!parameters.follow || clearsPtFloor(parameters, count, state));
}

namespace follow
{

/// The road fit's outward pass of each track's rows less its outermost,
/// a pass of its own before the follow
/// decides the end again (`FollowStage`'s DecidesEnd): \a starts gets the
/// state, or one on no surface where the track has fewer than four rows --
/// three left, a degree of freedom -- or is as long as a chain can be, or the
/// fit refuses them. Not with a kink allowance: the road fit takes its last
/// step with it (`fitOutward`'s LastStepAllowance), this one does not.
template <DetectorLayout Det, typename Material>
struct EndFit {
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using TrackView = typename BoundedColumns<TrackOf<kMaxClusters>::template RecordT>::ConstView;

  const IterationParameters<Det>* parameters = nullptr;
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  const Material* material = nullptr;
  Propagator propagator;
  TrackView tracks;
  BoundTrackState* starts = nullptr;
  std::uint32_t trackOffset = 0U;

  GPUhdi() void operator()(std::uint32_t element) const
  {
    const std::uint32_t track = trackOffset + element;
    starts[element].surface = kInvalidSensorID;
    const std::uint32_t count = tracks.count[track];
    if (tracks.status[track] != PropagationStatus::kOk || count < 4U || count >= kMaxClusters) {
      return;
    }
    countWork(WorkCounter::kFollowEndFits);
    const auto fit = fitOutward(propagator, *geometry, *material, RowChain<Det>{.frame = frame, .rows = tracks.clusters[track].index}, static_cast<int>(count - 1U), parameters->roadChi2MaxByNdf, kinkBudget(*parameters));
    if (fit.status == RoadFitStatus::kAccepted) {
      starts[element] = fit.state;
    }
  }
};

/// Each direction cut back so that every kink step kept keeps the row that
/// witnessed it (`keptWithWitnesses`): after the bidding, before the refit.
template <std::size_t K>
struct KeepWitnesses {
  FollowProposal<K>* proposals = nullptr;
  const KinkSteps* steps = nullptr;

  GPUhdi() void operator()(std::uint32_t track) const
  {
    auto& proposal = proposals[track];
    const std::uint32_t before = static_cast<std::uint32_t>(proposal.outwardCount) + proposal.inwardCount;
    proposal.outwardCount = keptWithWitnesses(proposal.outwardCount, steps[track].outwardWitnessed);
    proposal.inwardCount = keptWithWitnesses(proposal.inwardCount, steps[track].inwardWitnessed);
    countWork(WorkCounter::kFollowWitnessTrimmed, before - proposal.outwardCount - proposal.inwardCount);
  }
};

template <std::size_t K>
struct Bid {
  const FollowProposal<K>* proposals = nullptr;
  std::uint32_t* slots = nullptr;
  std::uint32_t trackOffset = 0U;
  std::uint32_t trackCount = 0U;

  GPUhdi() void operator()(std::uint32_t track) const
  {
    const std::uint32_t priority = trackCount - (trackOffset + track);
    const auto& proposal = proposals[track];
    for (std::uint32_t i = 0; i < proposal.outwardCount; ++i) {
      o2::gpu::CAMath::AtomicMax(slots + proposal.outward[i], priority);
    }
    for (std::uint32_t i = 0; i < proposal.inwardCount; ++i) {
      o2::gpu::CAMath::AtomicMax(slots + proposal.inward[i], priority);
    }
  }
};

/// Keep each direction only through the first row on which another track has
/// the higher output-order priority. A later row won by this track is still
/// dropped: one bidding round does not retry a truncated continuation.
/// \a earlier counts the earlier iterations' follow apart.
template <std::size_t K>
struct Resolve {
  FollowProposal<K>* proposals = nullptr;
  const std::uint32_t* slots = nullptr;
  std::uint32_t trackOffset = 0U;
  std::uint32_t trackCount = 0U;
  bool earlier = false;

  GPUhdi() void operator()(std::uint32_t track) const
  {
    const std::uint32_t priority = trackCount - (trackOffset + track);
    auto& proposal = proposals[track];
    const std::uint32_t proposed = static_cast<std::uint32_t>(proposal.outwardCount) + proposal.inwardCount;
    proposal.outwardCount = kept(proposal.outward, proposal.outwardCount, priority);
    proposal.inwardCount = kept(proposal.inward, proposal.inwardCount, priority);
    countWork(earlier ? WorkCounter::kFollowEarlierProposed : WorkCounter::kFollowProposed, proposed);
    countWork(earlier ? WorkCounter::kFollowEarlierOutbid : WorkCounter::kFollowOutbid, proposed - proposal.outwardCount - proposal.inwardCount);
  }

 private:
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint8_t kept(const std::uint32_t* rows, std::uint8_t count, std::uint32_t priority) const
  {
    for (std::uint8_t i = 0U; i < count; ++i) {
      if (slots[rows[i]] != priority) {
        return i;
      }
    }
    return count;
  }
};

template <std::size_t K>
struct Combine {
  const ClusterList<K>* original = nullptr;
  const std::uint8_t* originalCount = nullptr;
  const FollowProposal<K>* proposals = nullptr;
  ClusterList<K>* rows = nullptr;
  std::uint8_t* counts = nullptr;
  std::uint32_t trackOffset = 0U;

  GPUhdi() void operator()(std::uint32_t element) const
  {
    const std::uint32_t track = trackOffset + element;
    const auto& proposal = proposals[element];
    ClusterList<K> combined;
    for (std::size_t i = 0; i < K; ++i) {
      combined.index[i] = ClusterList<K>::kNoCluster;
    }
    for (std::uint32_t i = 0; i < proposal.inwardCount; ++i) {
      combined.index[i] = proposal.inward[proposal.inwardCount - i - 1U];
    }
    // A replaced end goes with the outward rows that replace it.
    const bool replaced = proposal.replacedEnd != ClusterList<K>::kNoCluster && proposal.outwardCount != 0U;
    const std::uint32_t kept = originalCount[track] - (replaced ? 1U : 0U);
    for (std::uint32_t i = 0; i < kept; ++i) {
      combined.index[static_cast<std::uint32_t>(proposal.inwardCount) + i] = original[track].index[i];
    }
    for (std::uint32_t i = 0; i < proposal.outwardCount; ++i) {
      combined.index[static_cast<std::uint32_t>(proposal.inwardCount) + kept + i] = proposal.outward[i];
    }
    rows[element] = combined;
    counts[element] = static_cast<std::uint8_t>(kept + proposal.inwardCount + proposal.outwardCount);
  }
};

/// Whether \a fit, a track refitted after its follow, passes the road's
/// chi2 cut (`Refit`, `TrimAmbiguous`).
template <DetectorLayout Det, typename Fit>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool passesRoadChi2(const IterationParameters<Det>& parameters, const Fit& fit)
{
  return fit.status == PropagationStatus::kOk && fit.chi2 <= parameters.roadChi2MaxByNdf[fit.ndf];
}

/// The extended refit of each followed track, all or nothing. \a Rule is
/// the fit's choice of crossing (`CrossingRule`), the iteration's
/// `crossingRule`. With \a Prefix, the iteration's `followPrefix`, a refit the cut
/// refuses or that fails is tried again on the follow's rows' prefixes in
/// the walks' order -- the outward rows first, then the inward ones --
/// the longest first, and the first the refit accepts is kept; a prefix
/// that would keep a kink step without its witness (`keptWithWitnesses`)
/// is not tried.
template <DetectorLayout Det, typename Material, CrossingRule Rule = CrossingRule::kFirstAhead, bool Prefix = false>
struct Refit {
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using TrackView = typename BoundedColumns<TrackOf<kMaxClusters>::template RecordT>::View;

  const IterationParameters<Det>* parameters = nullptr;
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  const Material* material = nullptr;
  Propagator propagator;
  TrackView tracks;
  FollowProposal<kMaxClusters>* proposals = nullptr;
  const ClusterList<kMaxClusters>* rows = nullptr;
  const std::uint8_t* counts = nullptr;
  std::uint32_t trackOffset = 0U;
  bool earlier = false;
  /// Each proposal's kink steps, where the iteration takes them; a prefix
  /// keeps them with their witnesses.
  const KinkSteps* steps = nullptr;

  GPUhdi() void operator()(std::uint32_t element) const
  {
    const std::uint32_t track = trackOffset + element;
    auto& proposal = proposals[element];
    const std::uint32_t added = static_cast<std::uint32_t>(proposal.outwardCount) + proposal.inwardCount;
    if (added == 0U) {
      return;
    }
    const std::uint8_t count = counts[element];
    if (!accept(track, rows[element], count, proposal.outwardCount, proposal.inwardCount)) {
      if constexpr (Prefix) {
        if (keepPrefix(element, track, proposal, count)) {
          return;
        }
      }
      countWork(earlier ? WorkCounter::kFollowEarlierRefitRejected : WorkCounter::kFollowRefitRejected, added);
      proposal.outwardCount = 0U;
      proposal.inwardCount = 0U;
    }
  }

 private:
  /// The track with \a list's rows (\a n of them, \a outward and \a inward
  /// of them the follow's) refitted: whether the fit passes the road's chi2
  /// cut; then it is written as the track, its added counts grown -- added
  /// to what an earlier follow added: an earlier iteration's track is
  /// followed again with `followEarlier`.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool accept(std::uint32_t track, const ClusterList<kMaxClusters>& list, std::uint8_t n, std::uint8_t outward, std::uint8_t inward) const
  {
    const auto fit = fitMeasurementsReference<Rule>(propagator, *geometry, *material, RowChain<Det>{.frame = frame, .rows = list.index}, static_cast<int>(n), parameters->fitReference, kinkBudget(*parameters));
    if (!passesRoadChi2(*parameters, fit)) {
      return false;
    }
    writeFit(tracks, track, fit, list, n, trackTimeBox(frame->measurements(), list, n));
    tracks.addedOutward[track] = static_cast<std::uint8_t>(tracks.addedOutward[track] + outward);
    tracks.addedInward[track] = static_cast<std::uint8_t>(tracks.addedInward[track] + inward);
    return true;
  }

  /// `Prefix`'s second chance: the combined rows are the inward ones
  /// innermost first, the track's, the outward ones; a prefix of \a outward
  /// outward and \a inward inward rows is the track's rows with the first
  /// of each. Whether one was kept, its counts then the proposal's.
  GPUhdi() GPUAlwaysInline() bool keepPrefix(std::uint32_t element, std::uint32_t track, FollowProposal<kMaxClusters>& proposal, std::uint8_t count) const
  {
    const auto outwardAdded = static_cast<std::uint32_t>(proposal.outwardCount);
    const auto inwardAdded = static_cast<std::uint32_t>(proposal.inwardCount);
    const std::uint32_t kept = count - outwardAdded - inwardAdded;
    const auto& combined = rows[element];
    for (std::uint32_t length = outwardAdded + inwardAdded; length-- > 1U;) {
      const auto outward = length < outwardAdded ? length : outwardAdded;
      const auto inward = length - outward;
      if (steps != nullptr && (keptWithWitnesses(static_cast<std::uint8_t>(outward), steps[element].outwardWitnessed) != outward || keptWithWitnesses(static_cast<std::uint8_t>(inward), steps[element].inwardWitnessed) != inward)) {
        continue;
      }
      ClusterList<kMaxClusters> prefix{};
      std::uint32_t n = 0U;
      for (std::uint32_t i = inward; i-- > 0U;) {
        prefix.index[n++] = proposal.inward[i];
      }
      for (std::uint32_t i = 0; i < kept; ++i) {
        prefix.index[n++] = combined.index[inwardAdded + i];
      }
      for (std::uint32_t i = 0; i < outward; ++i) {
        prefix.index[n++] = proposal.outward[i];
      }
      if (!accept(track, prefix, static_cast<std::uint8_t>(n), static_cast<std::uint8_t>(outward), static_cast<std::uint8_t>(inward))) {
        continue;
      }
      countWork(earlier ? WorkCounter::kFollowEarlierRefitRejected : WorkCounter::kFollowRefitRejected, outwardAdded + inwardAdded - length);
      proposal.outwardCount = static_cast<std::uint8_t>(outward);
      proposal.inwardCount = static_cast<std::uint8_t>(inward);
      return true;
    }
    return false;
  }
};

/// The trim of an ambiguous last row, after the refit, where the iteration has a `followTrimMargin`: a
/// track whose last row is one its follow added outward, neither a kink
/// step's row nor a kink row's witness, is refitted without it. Where that
/// fit passes the road's chi2 cut, the track one row shorter keeps the cuts
/// after the follow (`keptAfterFollow`), and from the fit's outward state two
/// rows or more of that row's layer pass the nominal gate (`passingOn`: free
/// rows, the row itself among them -- nothing is marked before `Mark`) with
/// the row itself not among them, or the smallest other chi2 less the row's
/// below the margin (a margin above the gate takes every such
/// row), the row is dropped and the track written with that fit; else the
/// track stays as the refit left it. The row dropped is not marked: it is free for the
/// iterations after. Not with `followEnd`, whose replaced end goes with the
/// outward rows.
template <DetectorLayout Det, typename Material, CrossingRule Rule = CrossingRule::kFirstAhead>
struct TrimAmbiguous {
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using TrackView = typename BoundedColumns<TrackOf<kMaxClusters>::template RecordT>::View;

  /// The search of the step's rows, as a nominal step makes it.
  FollowStage<Det, Material> search;
  const IterationParameters<Det>* parameters = nullptr;
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  const Material* material = nullptr;
  Propagator propagator;
  TrackView tracks;
  FollowProposal<kMaxClusters>* proposals = nullptr;
  /// Each proposal's kink steps, where the iteration takes them.
  const KinkSteps* steps = nullptr;
  std::uint32_t trackOffset = 0U;

  GPUhdi() void operator()(std::uint32_t element) const
  {
    const std::uint32_t track = trackOffset + element;
    auto& proposal = proposals[element];
    const std::uint32_t outward = proposal.outwardCount;
    const std::uint32_t count = tracks.count[track];
    if (outward == 0U || tracks.status[track] != PropagationStatus::kOk) {
      return;
    }
    // A kink row stands on the row after it (bit i for row i): the last row
    // is that row's witness where the row before it is one.
    if (steps != nullptr) {
      const std::uint64_t witnessed = steps[element].outwardWitnessed;
      const std::uint32_t last = outward - 1U;
      if (((witnessed >> last) & 1U) != 0U || (last != 0U && ((witnessed >> (last - 1U)) & 1U) != 0U)) {
        return;
      }
    }
    // The outward rows are the track's last (`Combine`).
    ClusterList<kMaxClusters> rows = tracks.clusters[track];
    const std::uint32_t row = rows.index[count - 1U];
    rows.index[count - 1U] = ClusterList<kMaxClusters>::kNoCluster;
    const auto shorter = static_cast<std::uint8_t>(count - 1U);
    countWork(WorkCounter::kFollowTrimFits);
    const auto fit = fitMeasurementsReference<Rule>(propagator, *geometry, *material, RowChain<Det>{.frame = frame, .rows = rows.index}, static_cast<int>(shorter), parameters->fitReference, kinkBudget(*parameters));
    if (!passesRoadChi2(*parameters, fit) || !keptAfterFollow(*parameters, shorter, fit.state)) {
      return;
    }
    const auto box = trackTimeBox(frame->measurements(), rows, shorter);
    const LayerID layer = geometry->surface(frame->measurements().surface[row]).layer;
    const auto passing = search.passingOn(fit.outwardState, factorCovariance(fit.outwardState.cov), layer, row, box, PropagationDirection::kForward);
    if (passing.count < 2U || (passing.own && !(passing.otherChi2 - passing.ownChi2 < parameters->followTrimMargin))) {
      return;
    }
    writeFit(tracks, track, fit, rows, shorter, box);
    tracks.addedOutward[track] = static_cast<std::uint8_t>(tracks.addedOutward[track] - 1U);
    proposal.outwardCount = static_cast<std::uint8_t>(outward - 1U);
    countWork(WorkCounter::kFollowTrimmed);
  }
};

template <std::size_t K>
struct Mark {
  const FollowProposal<K>* proposals = nullptr;
  std::uint8_t* used = nullptr;

  GPUhdi() void operator()(std::uint32_t track) const
  {
    const auto& proposal = proposals[track];
    // A replaced end, its replacement kept through the refit, is free again:
    // no other track proposed it, so no other one marks it.
    if (proposal.replacedEnd != ClusterList<K>::kNoCluster && proposal.outwardCount != 0U) {
      used[proposal.replacedEnd] = kFree;
    }
    for (std::uint32_t i = 0; i < proposal.outwardCount; ++i) {
      used[proposal.outward[i]] = kExclusive;
    }
    for (std::uint32_t i = 0; i < proposal.inwardCount; ++i) {
      used[proposal.inward[i]] = kExclusive;
    }
  }
};

} // namespace follow

template <DetectorLayout Det>
class Follow;

namespace follow
{

/// What the follow's passes read and write (`runSequence`): the iteration's
/// parameters as the passes read them and as the host decides on them (one
/// object on the host, the device's copy and the host's on a device), the
/// track model, the tracks from \a trackOffset to \a trackCount, the used
/// mask (\a rowCount rows) and the passes' storage, one element per track
/// followed -- the kink steps only where the iteration takes them, the end
/// starts only where it decides the end again.
template <DetectorLayout Det, typename Material>
struct Inputs {
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using TrackView = typename BoundedColumns<TrackOf<kMaxClusters>::template RecordT>::View;
  using ConstTrackView = typename BoundedColumns<TrackOf<kMaxClusters>::template RecordT>::ConstView;

  const IterationParameters<Det>* parameters = nullptr;
  const IterationParameters<Det>* choices = nullptr;
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  const Material* material = nullptr;
  Propagator propagator;
  TrackView tracks;
  ConstTrackView constTracks;
  std::uint32_t trackOffset = 0U;
  std::uint32_t trackCount = 0U;
  bool inward = false;
  bool earlier = false;
  std::uint8_t* used = nullptr;
  std::size_t rowCount = 0U;
  const std::uint32_t* maxCandidateWidthPS = nullptr;
  SensorCellsView cells;
  FollowProposal<kMaxClusters>* proposals = nullptr;
  std::uint32_t* bids = nullptr;
  ClusterList<kMaxClusters>* rows = nullptr;
  std::uint8_t* counts = nullptr;
  KinkSteps* steps = nullptr;
  BoundTrackState* endStarts = nullptr;
};

/// The host's launches of the follow's passes: each in its profile zone,
/// none of its own for the earlier iterations' tracks (`Follow::run`).
template <typename ForEach>
struct HostLaunch {
  ForEach* forEach;
  bool earlier = false;

  [[nodiscard]] profile::Phase phase(profile::Phase own) const { return earlier ? profile::Phase::kNone : own; }

  template <typename Body>
  void pass(profile::Phase own, std::uint32_t count, const Body& body) const
  {
    const profile::Zone zone{phase(own)};
    (*forEach)(count, body);
  }

  template <typename Stage>
  void propose(profile::Phase own, std::uint32_t count, const Stage& stage) const
  {
    pass(own, count, stage);
  }

  /// The bids zeroed in parallel, as the selection's slots are.
  void clearBids(std::uint32_t* bids, std::size_t rows) const
  {
    const profile::Zone zone{phase(profile::Phase::kFollowPropose)};
    (*forEach)(static_cast<std::uint32_t>(rows), [bids](std::uint32_t row) { bids[row] = 0U; });
  }
};

/// The follow over \a in's tracks, the host's and a device's alike: the bids
/// cleared, the end fit where the iteration decides the end again (on its
/// own tracks only), the proposals, then bid, resolve, the kink steps'
/// witnesses, combine, the refit, the trim of the ambiguous last rows
/// (`followTrimMargin`, not with the end decided again, on the iteration's
/// own tracks) and the mark. Each choice of the iteration is its own
/// instantiation. \a launch runs a pass (`pass(phase, count, body)`), the
/// proposals (`propose(phase, count, stage)`: a device's lanes) and
/// `clearBids(bids, rows)`.
template <DetectorLayout Det, typename Material, typename Launch>
void runSequence(const Launch& launch, const Inputs<Det, Material>& in)
{
  using profile::Phase;
  constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  const IterationParameters<Det>& choices = *in.choices;
  const std::uint32_t elements = in.trackCount - in.trackOffset;
  launch.clearBids(in.bids, in.rowCount);
  // The beam line starts the inward follow where the iteration has the
  // constraint, and the end is decided
  // again where it has followEnd, on its own tracks
  // only.
  const bool end = choices.followEnd && !in.earlier;
  if (end) {
    launch.pass(Phase::kFollowPropose, elements,
                EndFit<Det, Material>{
                  .parameters = in.parameters,
                  .frame = in.frame,
                  .geometry = in.geometry,
                  .material = in.material,
                  .propagator = in.propagator,
                  .tracks = in.constTracks,
                  .starts = in.endStarts,
                  .trackOffset = in.trackOffset,
                });
  }
  // The iteration's choice among the passing rows, its end decisions on the
  // last row's density and against the beam line (with the beam line only).
  withBools(
    [&]<bool BeamLine, bool DecidesEnd, bool Greedy, bool EndDensity, bool BeamEnd>(std::bool_constant<BeamLine>, std::bool_constant<DecidesEnd>, std::bool_constant<Greedy>, std::bool_constant<EndDensity>, std::bool_constant<BeamEnd>) {
      if constexpr (BeamLine || !BeamEnd) {
        const FollowStage<Det, Material, BeamLine, DecidesEnd, Greedy, EndDensity, BeamEnd> stage{in.parameters, in.frame, in.geometry, in.material, in.propagator, in.constTracks, in.used, in.maxCandidateWidthPS, in.cells, in.proposals, in.steps, in.inward, in.trackOffset, DecidesEnd ? in.endStarts : nullptr};
        launch.propose(Phase::kFollowPropose, elements, stage);
      }
    },
    choices.beamConstraint, end, choices.followGreedy, Follow<Det>::endDensity(choices), choices.followBeamEnd);
  launch.pass(Phase::kFollowResolve, elements, Bid<kMaxClusters>{in.proposals, in.bids, in.trackOffset, in.trackCount});
  launch.pass(Phase::kFollowResolve, elements, Resolve<kMaxClusters>{in.proposals, in.bids, in.trackOffset, in.trackCount, in.earlier});
  if (in.steps != nullptr) {
    launch.pass(Phase::kFollowResolve, elements, KeepWitnesses<kMaxClusters>{.proposals = in.proposals, .steps = in.steps});
  }
  launch.pass(Phase::kFollowResolve, elements, Combine<kMaxClusters>{in.constTracks.clusters, in.constTracks.count, in.proposals, in.rows, in.counts, in.trackOffset});
  // The iteration's choice of crossing and of
  // the prefix.
  withCrossingRule(choices.crossingRule, [&]<CrossingRule Rule>(std::integral_constant<CrossingRule, Rule>) { withBools([&]<bool Prefix>(std::bool_constant<Prefix>) { launch.pass(Phase::kFollowRefit, elements, Refit<Det, Material, Rule, Prefix>{in.parameters, in.frame, in.geometry, in.material, in.propagator, in.tracks, in.proposals, in.rows, in.counts, in.trackOffset, in.earlier, in.steps}); }, choices.followPrefix); });
  // The ambiguous last rows trimmed, on the iteration's own
  // tracks, before the marks.
  if (choices.followTrimMargin > 0.F && !end && !in.earlier) {
    const FollowStage<Det, Material> search{in.parameters, in.frame, in.geometry, in.material, in.propagator, in.constTracks, in.used, in.maxCandidateWidthPS, in.cells, nullptr, nullptr, false, in.trackOffset};
    withCrossingRule(choices.crossingRule, [&]<CrossingRule Rule>(std::integral_constant<CrossingRule, Rule>) {
      launch.pass(Phase::kFollowRefit, elements,
                  TrimAmbiguous<Det, Material, Rule>{
                    .search = search,
                    .parameters = in.parameters,
                    .frame = in.frame,
                    .geometry = in.geometry,
                    .material = in.material,
                    .propagator = in.propagator,
                    .tracks = in.tracks,
                    .proposals = in.proposals,
                    .steps = in.steps,
                    .trackOffset = in.trackOffset,
                  });
    });
  }
  launch.pass(Phase::kFollowMark, elements, Mark<kMaxClusters>{in.proposals, in.used});
}

} // namespace follow

/// Host orchestration and bounded storage of the follow: propose, then bid,
/// resolve, keep the kink steps' witnesses, combine, refit, trim the
/// ambiguous last rows (`followTrimMargin`) and mark. The
/// tracker runs it on an iteration's own tracks and, with `followEarlier`,
/// on the earlier iterations' tracks first; a device runs the same sequence
/// (`follow::runSequence`) with its own launches.
template <DetectorLayout Det>
class Follow
{
 public:
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  static constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);
  using Proposal = FollowProposal<kMaxClusters>;
  using Tracks = typename TrackRefit<Det>::Tracks;

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool = std::move(pool);
    clearStorage();
  }

  /// The arrays back to the pool.
  void clearStorage() noexcept { deepVectorClear(mMemoryPool.get(), mProposals, mBids, mRows, mCounts, mKinkSteps, mEndStarts); }

  /// Follow the tracks from \a trackOffset on, inward too with \a inward,
  /// each layer's candidates bound by \a maxCandidateWidthPS, the frame's
  /// widest time box there, the sensors
  /// near a crossing from \a cells. \a earlier: the earlier iterations'
  /// tracks (`followEarlier`), counted and timed apart -- within the
  /// caller's `kFollowEarlier` zone, no phase of its own.
  template <typename Material, typename ForEach>
  [[nodiscard]] StageStatus run(ForEach forEach, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator, Tracks& tracks, std::uint32_t trackOffset, bool inward, bool earlier, std::span<std::uint8_t> used, std::span<const std::uint32_t, kLayerCount> maxCandidateWidthPS, SensorCellsView cells)
  {
    using profile::Phase;
    return guardOutOfMemory([&] {
      const auto trackCount = static_cast<std::uint32_t>(tracks.size());
      const auto elements = trackCount - trackOffset;
      const bool end = parameters.followEnd && !earlier;
      {
        const profile::Zone zone{earlier ? Phase::kNone : Phase::kFollowPropose};
        // Kept from iteration to iteration and filled by the passes, not
        // reallocated.
        resizeBoundedVector(mProposals, elements, mMemoryPool.get());
        resizeBoundedVector(mBids, used.size(), mMemoryPool.get());
        resizeBoundedVector(mRows, elements, mMemoryPool.get());
        resizeBoundedVector(mCounts, elements, mMemoryPool.get());
        resizeBoundedVector(mKinkSteps, takesKinkSteps(parameters) ? elements : 0U, mMemoryPool.get());
        if (end) {
          resizeBoundedVector(mEndStarts, elements, mMemoryPool.get());
        }
      }
      follow::runSequence(follow::HostLaunch<ForEach>{.forEach = &forEach, .earlier = earlier}, follow::Inputs<Det, Material>{
                                                                                                  .parameters = &parameters,
                                                                                                  .choices = &parameters,
                                                                                                  .frame = &frame,
                                                                                                  .geometry = &geometry,
                                                                                                  .material = &material,
                                                                                                  .propagator = propagator,
                                                                                                  .tracks = tracks.view(),
                                                                                                  .constTracks = static_cast<const Tracks&>(tracks).view(),
                                                                                                  .trackOffset = trackOffset,
                                                                                                  .trackCount = trackCount,
                                                                                                  .inward = inward,
                                                                                                  .earlier = earlier,
                                                                                                  .used = used.data(),
                                                                                                  .rowCount = used.size(),
                                                                                                  .maxCandidateWidthPS = maxCandidateWidthPS.data(),
                                                                                                  .cells = cells,
                                                                                                  .proposals = mProposals.data(),
                                                                                                  .bids = mBids.data(),
                                                                                                  .rows = mRows.data(),
                                                                                                  .counts = mCounts.data(),
                                                                                                  .steps = mKinkSteps.empty() ? nullptr : mKinkSteps.data(),
                                                                                                  .endStarts = end ? mEndStarts.data() : nullptr,
                                                                                                });
      return StageStatus::kComplete;
    });
  }

  /// Each track's proposal after the passes: the rows it took.
  [[nodiscard]] std::span<const Proposal> proposals() const { return mProposals; }

  /// Each proposal's kink steps, as proposed (before the bidding); empty
  /// without them (`takesKinkSteps`).
  [[nodiscard]] std::span<const follow::KinkSteps> kinkSteps() const { return mKinkSteps; }

  /// Whether an iteration's follow takes kink steps: it has a kink allowance.
  [[nodiscard]] static bool takesKinkSteps(const IterationParameters<Det>& parameters) { return parameters.kinkScale2 > 1.F && parameters.maxKinks != 0U; }

  /// Whether an iteration's follow decides its ends on the last row's
  /// density (`FollowStage`'s EndDensity): either of its keys is set.
  [[nodiscard]] static bool endDensity(const IterationParameters<Det>& parameters) { return parameters.followEndMaxSqrtDet > 0.F || parameters.followKeepDensity > 0.F; }

  /// Proposal \a first + i taken from proposal kept[i], for every i, as the
  /// tracker drops short tracks in parallel: \a kept ascends
  /// and kept[i] >= first + i, so the moves go through a scratch.
  /// \a forEach(count, body) runs body(i) for i in [0, count).
  template <typename Index, typename ForEach>
  void keepProposals(std::span<const Index> kept, std::size_t first, const ForEach& forEach)
  {
    const auto count = static_cast<Index>(kept.size());
    const auto move = [&](auto& values) {
      using Value = typename std::remove_reference_t<decltype(values)>::value_type;
      BoundedVector<Value> scratch(std::pmr::polymorphic_allocator<Value>{mMemoryPool.get()});
      resizeBoundedVector(scratch, count);
      forEach(count, [&](Index at) { scratch[at] = values[kept[at]]; });
      forEach(count, [&](Index at) { values[first + at] = scratch[at]; });
    };
    move(mProposals);
    if (!mKinkSteps.empty()) {
      move(mKinkSteps);
    }
  }
  void truncateProposals(std::size_t count)
  {
    mProposals.resize(count);
    if (!mKinkSteps.empty()) {
      mKinkSteps.resize(count);
    }
  }

 private:
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  BoundedVector<Proposal> mProposals;
  BoundedVector<std::uint32_t> mBids;
  BoundedVector<ClusterList<kMaxClusters>> mRows;
  BoundedVector<std::uint8_t> mCounts;
  BoundedVector<follow::KinkSteps> mKinkSteps;
  BoundedVector<BoundTrackState> mEndStarts; ///< with followEnd: each track's state without its outermost row
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_FOLLOW_H_

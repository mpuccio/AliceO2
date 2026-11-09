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

/// \file TrackFitter.h
/// \brief Independent outward and inward Kalman filters over a chain, each step linearised about its filtered state.

#ifndef GTRK_TRACKING_TRACKFITTER_H_
#define GTRK_TRACKING_TRACKFITTER_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/CovarianceTransport.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/MeasurementUpdate.h"
#include "GenericTracking/Tracking/Seed.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include "GPUCommonMath.h"

#include <cstdint>

namespace o2::gtrk
{

/// How many steps of a filter pass may be kinks: a
/// step whose spatial increment is beyond \a gate, with a measurement after
/// it in the pass, is taken again from the surface before with the track's
/// direction freed -- a break point, two degrees of freedom. None by default.
/// The last step, with nothing after it to check a freed direction, never
/// is; the road fit's takes the kink allowance there instead
/// (`fitOutward`'s LastStepAllowance).
struct KinkBudget {
  int max = 0;
  Scalar gate = 0.F; ///< `chiSquareQuantile(2, nSigma)`
};

/// An iteration's kink budget, from its parameters' `maxKinks` and
/// `kinkChi2Min`: the refit's and the extension's.
template <typename Parameters>
[[nodiscard]] GPUhdi() GPUAlwaysInline() KinkBudget kinkBudget(const Parameters& parameters)
{
  return KinkBudget{.max = static_cast<int>(parameters.maxKinks), .gate = parameters.kinkChi2Min};
}

/// Whether an iteration's road fit takes its kink allowance on the last step
/// (`fitOutward`'s LastStepAllowance, against the road-end fakes): every iteration with an allowance. The one rule the host's and the
/// device's road-fit stages are instantiated by and `fitRoad` follows.
template <typename Parameters>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool takesLastStepAllowance(const Parameters& parameters)
{
  return parameters.kinkScale2 > 1.F;
}

/// Result of the two independent filter passes.
///
/// `state` is the inward result on the innermost measurement surface. The
/// outward state is the starting point for extension; its chi2 remains
/// available for QA.
struct TrackFitResult {
  BoundTrackState state{};
  BoundTrackState outwardState{};
  Scalar chi2 = 0.F;
  Scalar outwardChi2 = 0.F;
  /// Two coordinates per measurement against the five parameters, two fewer
  /// per kink.
  std::int16_t ndf = 0;
  /// kOk, or what the failed step's propagation gave (`failedPass`,
  /// `failedMeasurement` say where).
  PropagationStatus status = PropagationStatus::kOk;
  FitPass failedPass = FitPass::kNone;
  /// The measurement the failed step was propagating to.
  std::int16_t failedMeasurement = -1;
  /// The kinks each pass took. `ndf` is the inward pass's, two fewer per
  /// kink; the outward pass has `2 N - 5 - 2 outwardKinks`.
  std::uint8_t kinks = 0;
  std::uint8_t outwardKinks = 0;
  /// The measurement the inward pass's last kink updated with; -1 without.
  std::int16_t kinkMeasurement = -1;
};

/// A fit's per-measurement chi2 increments, for QA: called once per update
/// with the pass's direction, the measurement's index in the chain and its
/// increment (two dof). The default does nothing.
struct NoFitObserver {
  GPUhdi() GPUAlwaysInline() void operator()(PropagationDirection /*direction*/, int /*index*/, Scalar /*chi2*/) const {}
};

/// The mass a fit step's energy loss takes: the propagator's hypothesis,
/// the one its scattering takes. A QA oracle's propagator type overloads it
/// to set the two apart (SpeciesOracle.h).
GPUhdi() GPUAlwaysInline() Scalar lossMassGeV(const Propagator& prop) { return prop.massGeV; }

namespace detail
{

/// A step's material in two halves about its midpoint, and its `rho * length`
/// for the energy loss.
///
/// The Highland noise of the half nearer the surface the step leaves goes into
/// the covariance *before* the transport, so that its kink becomes position
/// spread at the next surface through `F Q F^T`; the other half's goes in
/// after it, before the update. Added after the transport, all of it, every
/// prediction would miss `theta0^2 L^2` of the material near the start and the
/// outward state would carry one step's unconsumed `theta0^2` in its direction.
/// With a layer's material
/// attributed from half-gap to half-gap, the two halves are the stave just
/// left and the stave about to be reached.
struct StepMaterial {
  Scalar nearerX2X0 = 0.F;  ///< x/X0 from the start to the midpoint
  Scalar fartherX2X0 = 0.F; ///< x/X0 from the midpoint to the end
  Scalar xRho = 0.F;        ///< [g/cm^2] over the whole step
};

template <typename Material>
GPUdi() StepMaterial stepMaterial(const Material& material, const Scalar start[3], const Scalar end[3])
{
  const Scalar mid[3] = {0.5F * (start[0] + end[0]), 0.5F * (start[1] + end[1]), 0.5F * (start[2] + end[2])};
  const auto nearer = material.getMatBudget(start[0], start[1], start[2], mid[0], mid[1], mid[2]);
  const auto farther = material.getMatBudget(mid[0], mid[1], mid[2], end[0], end[1], end[2]);
  return {.nearerX2X0 = nearer.meanX2X0, .fartherX2X0 = farther.meanX2X0, .xRho = (nearer.meanRho * nearer.length) + (farther.meanRho * farther.length)};
}

/// The mass hypothesis's energy loss over a step's \a xRho [g/cm^2], O2's
/// way: taken at \a state's momentum and applied to it, the straggling's q/p
/// variance returned in \a varQOverP for the step's noise (StepNoise). Lost
/// along the momentum, regained against it. False if the particle stopped.
/// The mass is `lossMassGeV(prop)`, the propagator's hypothesis.
template <typename Prop>
GPUhdi() GPUAlwaysInline() bool applyEnergyLoss(const Prop& prop, Scalar xRho, PropagationDirection direction, BoundTrackState& state, Scalar& varQOverP)
{
  varQOverP = 0.F;
  const auto loss = energyLoss(state.qOverP(), lossMassGeV(prop), direction == PropagationDirection::kForward ? -xRho : xRho);
  if (!loss.ok) {
    return false;
  }
  state.par[kQOverP] *= loss.momentumScale;
  varQOverP = loss.varQOverP;
  return true;
}

/// \a theta0 times \a scale with \a Scaled, \a theta0 alone without: the
/// nominal step never multiplies.
template <bool Scaled>
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar scaledTheta0(Scalar theta0, [[maybe_unused]] Scalar scale)
{
  if constexpr (Scaled) {
    return scale * theta0;
  } else {
    return theta0;
  }
}

/// A step's second half, once \a state is on the plane the step reaches and
/// \a jacobian carries the covariance there: the farther half's scattering
/// of \a step at the direction \a state arrives with, the energy loss over
/// \a step at its momentum, and the covariance \a factors transported
/// with \a before, the nearer half's noise, and the farther half's
/// (StepMaterial). With \a Scaled, \a noiseScale multiplies the farther
/// half's theta0 (`scaledTheta0`). The fit's step (`filterStep`) and the
/// follow's (`transportCandidate`) both end so; each takes the nearer half's
/// noise before its propagation.
template <bool Scaled, typename Prop>
GPUdi() PropagationStatus transportStep(const Prop& prop, const StepMaterial& step, Scalar betaP, const ScatteringNoise& before, PropagationDirection direction, const SolenoidJacobian& jacobian, BoundTrackState& state, CovarianceFactors& factors, Scalar noiseScale)
{
  const auto after = scatteringNoise(scaledTheta0<Scaled>(highlandTheta0(step.fartherX2X0, betaP), noiseScale), state.theta());
  StepNoise noise{.phiBefore = before.phi, .thetaBefore = before.theta, .phiAfter = after.phi, .thetaAfter = after.theta};
  if (!applyEnergyLoss(prop, step.xRho, direction, state, noise.qOverPAfter)) {
    return PropagationStatus::kStopped;
  }
  if (!transportFactors(jacobian, noise, factors)) {
    return PropagationStatus::kNotFinite;
  }
  return PropagationStatus::kOk;
}

/// One step of a filter pass: \a state carried to \a measurement's plane
/// and updated with it, linearised about itself -- the filtered state the step
/// starts from, as ACTS's Kalman filter (a fixed whole-chain reference, never
/// corrected by the measurements, falls short of the outer planes of a low-momentum track near
/// its turning point). The covariance is \a factors, the pass's
/// (CovarianceFactors): `state.cov` is neither read nor written. \a Scaled,
/// the step's scattering theta0 is multiplied by \a noiseScale, both halves:
/// the kink allowance's last step (`allowanceStep`); the nominal step is its
/// own instantiation. \a Rule is the propagation's choice of crossing
/// (`CrossingRule`), and with it every pass below. Never inlined on the
/// device: one copy for each instantiation, not one in every fit of every
/// stage variant.
template <bool Scaled = false, CrossingRule Rule = CrossingRule::kFirstAhead, DetectorLayout Det, typename Material, typename Prop>
GPUdnii() PropagationStatus filterStep(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, CovarianceFactors& factors, Scalar& chi2, [[maybe_unused]] Scalar noiseScale = 1.F)
{
  const auto& from = geometry.planeOf(state);
  const auto& to = geometry.planeOf(measurement);

  Scalar start[3]{};
  surfaceToGlobal(from, state.loc0(), state.loc1(), start);
  const Scalar leaving = state.theta();
  // An intersection just outside the sensor's active area is still the point
  // to update at: the cluster is there, so the track was. Bounds are for
  // searches, not fits (`jacobianAnalytic` takes kOutOfBounds).
  SolenoidJacobian jacobian;
  const auto jacobianStatus = jacobianAnalytic<Scalar, Rule>(prop, state, from, to, direction, jacobian, &state);
  if (jacobianStatus != PropagationStatus::kOk) {
    return jacobianStatus;
  }

  Scalar end[3]{};
  surfaceToGlobal(to, state.loc0(), state.loc1(), end);

  const auto step = stepMaterial(material, start, end);
  const Scalar momentum = 1.F / o2::gpu::GPUCommonMath::Abs(state.qOverP());
  const Scalar betaP = betaMomentum(momentum, prop.massGeV);
  // The nearer half before the transport (StepMaterial), at the direction on
  // the surface the step leaves; the farther half after it.
  const auto before = scatteringNoise(scaledTheta0<Scaled>(highlandTheta0(step.nearerX2X0, betaP), noiseScale), leaving);
  const auto transported = transportStep<Scaled>(prop, step, betaP, before, direction, jacobian, state, factors, noiseScale);
  if (transported != PropagationStatus::kOk) {
    return transported;
  }

  const auto update = updateWithMeasurement(state, factors, measurement);
  if (!updateUsable(update, factors)) {
    return PropagationStatus::kNotFinite;
  }
  chi2 += update.chi2;
  return PropagationStatus::kOk;
}

/// \a budget for a pass over \a n measurements: no more kinks than leave the
/// fit a degree of freedom, `2 n - 5 - 2 k >= 1`, so `n - 3` -- one for four
/// measurements, two for five. Each kink frees the direction, two parameters;
/// beyond that the chain no longer constrains the track, and the chi2 cut by
/// degrees of freedom has no entry for it.
GPUhdi() GPUAlwaysInline() KinkBudget kinkBudgetFor(KinkBudget budget, int n)
{
  const int allowed = ((2 * n) - kNBoundParams - 1) / 2;
  if (budget.max > allowed) {
    budget.max = allowed > 0 ? allowed : 0;
  }
  return budget;
}

/// `freeDirection` on a pass's factors: the covariance formed, its direction
/// freed, factored again. A kink takes it once; the direction it frees is
/// what made the covariance hard to hold, so the factoring loses nothing.
GPUhdi() GPUAlwaysInline() void freeDirection(BoundTrackState& state, CovarianceFactors& factors)
{
  covarianceFromFactors(factors, state.cov);
  o2::gtrk::freeDirection(state);
  factors = factorCovariance(state.cov);
}

/// A filter step that may be a kink: taken as `filterStep`; if kinks are left,
/// \a followed -- a measurement comes after this one in the pass -- and its
/// increment is beyond the gate, taken again from the surface before with the
/// direction freed, and counted in \a kinks.
template <CrossingRule Rule = CrossingRule::kFirstAhead, DetectorLayout Det, typename Material, typename Prop>
GPUdi() PropagationStatus kinkStep(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, CovarianceFactors& factors, Scalar& chi2, const KinkBudget& budget, bool followed, int& kinks)
{
  if (kinks >= budget.max || !followed) {
    return filterStep<false, Rule>(prop, geometry, material, measurement, direction, state, factors, chi2);
  }
  const BoundTrackState keptState = state;
  const CovarianceFactors keptFactors = factors;
  const Scalar keptChi2 = chi2;
  auto status = filterStep<false, Rule>(prop, geometry, material, measurement, direction, state, factors, chi2);
  if (status != PropagationStatus::kOk || !(chi2 - keptChi2 > budget.gate)) {
    return status;
  }
  state = keptState;
  factors = keptFactors;
  chi2 = keptChi2;
  freeDirection(state, factors);
  status = filterStep<false, Rule>(prop, geometry, material, measurement, direction, state, factors, chi2);
  kinks += status == PropagationStatus::kOk ? 1 : 0;
  return status;
}

/// The road fit's last step with the kink allowance (`fitOutward`'s
/// LastStepAllowance, against the road-end fakes): taken as
/// `filterStep`; if kinks are left and its increment is beyond the gate,
/// taken again from the surface before with its scattering's theta0 scaled by
/// \a noiseScale, as the extension's kink step, and counted in \a kinks;
/// \a refused if the increment is still beyond the gate. It frees no
/// parameter.
template <CrossingRule Rule = CrossingRule::kFirstAhead, DetectorLayout Det, typename Material>
GPUdi() PropagationStatus allowanceStep(const Propagator& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, CovarianceFactors& factors, Scalar& chi2, const KinkBudget& budget, Scalar noiseScale, int& kinks, bool& refused)
{
  refused = false;
  if (kinks >= budget.max) {
    return filterStep<false, Rule>(prop, geometry, material, measurement, direction, state, factors, chi2);
  }
  const BoundTrackState keptState = state;
  const CovarianceFactors keptFactors = factors;
  const Scalar keptChi2 = chi2;
  auto status = filterStep<false, Rule>(prop, geometry, material, measurement, direction, state, factors, chi2);
  if (status != PropagationStatus::kOk || !(chi2 - keptChi2 > budget.gate)) {
    return status;
  }
  state = keptState;
  factors = keptFactors;
  chi2 = keptChi2;
  status = filterStep<true, Rule>(prop, geometry, material, measurement, direction, state, factors, chi2, noiseScale);
  if (status == PropagationStatus::kOk) {
    ++kinks;
    refused = !(chi2 - keptChi2 <= budget.gate);
  }
  return status;
}

/// A pass's start on its first measurement: the
/// limit of its update under a position prior infinitely wide. The position
/// becomes the measurement's, with the measurement's covariance and none with
/// the other parameters, which keep the seed's prior. The seed's position is
/// not a prior the measurement could be tested against: the measurement
/// determines the position, its chi2 is zero, and the pass's `2 N - 5`
/// degrees of freedom already count it. As
/// `updateWithMeasurement`, with the factors of \a state's covariance.
GPUhdi() GPUAlwaysInline() MeasurementUpdateResult startOnMeasurement(BoundTrackState& state, CovarianceFactors& factors, const Measurement& measurement)
{
  MeasurementUpdateResult result;
  state.decorrelate(kLoc0, kLoc1);
  state.par[kLoc0] = measurement.loc0;
  state.par[kLoc1] = measurement.loc1;
  state.variance(kLoc0) = measurement.varLoc0;
  state.variance(kLoc1) = measurement.varLoc1;
  state.covAt(kLoc0, kLoc1) = measurement.covLoc01;
  factors = factorCovariance(state.cov);
  // Positive definite as the update's innovations would be: the measurement's
  // two whitened variances (`updateWithMeasurement`).
  result.positiveDefinite = measurement.varLoc0 > 0.F && measurement.varLoc1 - (measurement.covLoc01 * measurement.covLoc01 / measurement.varLoc0) > 0.F;
  return result;
}

/// A filter pass from \a seed over the chain. Its covariance is carried as
/// factors from the first update to the last (CovarianceFactors) and written
/// to \a result's `cov` once, at the end.
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Prop, DetectorLayout Det, typename Material, typename Chain, typename Observer>
GPUdi() PropagationStatus runPass(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, PropagationDirection direction, BoundTrackState seed, BoundTrackState& result, Scalar& chi2, std::int16_t& failedMeasurement, const KinkBudget& passBudget, int& kinks, std::int16_t& kinkMeasurement, const Observer& observe)
{
  const auto budget = kinkBudgetFor(passBudget, n);
  const int first = direction == PropagationDirection::kForward ? 0 : n - 1;
  const int increment = direction == PropagationDirection::kForward ? 1 : -1;

  BoundTrackState state = seed;
  CovarianceFactors factors;
  const auto firstUpdate = startOnMeasurement(state, factors, measurements[first]);
  if (!updateUsable(firstUpdate, factors)) {
    failedMeasurement = static_cast<std::int16_t>(first);
    return PropagationStatus::kNotFinite;
  }
  chi2 = firstUpdate.chi2;
  observe(direction, first, chi2);
  for (int step = 1; step < n; ++step) {
    const int index = first + step * increment;
    const Scalar before = chi2;
    const int kinksBefore = kinks;
    const auto& stepProp = stepPropagator(prop, geometry, measurements, index - increment, index);
    const auto status = kinkStep<Rule>(stepProp, geometry, material, measurements[index], direction, state, factors, chi2, budget, step < n - 1, kinks);
    if (status != PropagationStatus::kOk) {
      failedMeasurement = static_cast<std::int16_t>(index);
      return status;
    }
    if (kinks != kinksBefore) {
      kinkMeasurement = static_cast<std::int16_t>(index);
    }
    observe(direction, index, chi2 - before);
  }
  covarianceFromFactors(factors, state.cov);
  result = state;
  return PropagationStatus::kOk;
}

/// The inward pass of `fitMeasurements`: its seed, the outward pass's end
/// state with a diffuse covariance, then the filter from there in.
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Prop, DetectorLayout Det, typename Material, typename Chain, typename Observer>
GPUdi() void inwardPass(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, Scalar qOverPSigma, const KinkBudget& kinkBudget, const Observer& observe, TrackFitResult& fit)
{
  // From the outward pass's end state, its covariance made diffuse with the
  // chain's q/p width, as ACTS's reverse filter starts from the forward
  // one's: the whole-chain helix carried there unfitted misses
  // the outermost rows at low momentum.
  BoundTrackState outerSeed = fit.outwardState;
  setDiffuseCovariance(outerSeed, qOverPSigma);

  int kinks = 0;
  const auto status = runPass<Rule>(prop, geometry, material, measurements, n, PropagationDirection::kBackward, outerSeed, fit.state, fit.chi2, fit.failedMeasurement, kinkBudget, kinks, fit.kinkMeasurement, observe);
  if (status != PropagationStatus::kOk) {
    fit.status = status;
    fit.failedPass = FitPass::kInward;
    return;
  }
  fit.kinks = static_cast<std::uint8_t>(kinks);
  fit.ndf = static_cast<std::int16_t>(fit.ndf - (2 * kinks));
}

} // namespace detail

/// Fit an ordered chain containing one measurement per logical layer.
///
/// The outward filter starts on the first surface from a helix fitted
/// through every measurement of the chain (`makeChainSeed`) with a
/// diffuse covariance; the inward one on the last from the outward filter's
/// end state, its covariance made diffuse again, as ACTS's reverse filter.
/// Each step linearises about the filtered state it starts from
/// (`filterStep`). The inward filter does not inherit the outward
/// covariance, so no measurement's weight is counted twice in the returned
/// state.
///
/// The measurements' times are boxes, which the fit does not take: the
/// track's time is their intersection, set by the refit (TrackRefit.h).
///
/// \a measurements is anything indexed by `[i]` giving a `Measurement`, as
/// for `makeChainSeed`. Each pass may take \a kinkBudget's kinks. \a observe
/// sees every update's chi2 increments (`NoFitObserver`).
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Prop, DetectorLayout Det, typename Material, typename Chain, typename Observer = NoFitObserver>
GPUdi() TrackFitResult fitMeasurements(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, const KinkBudget& kinkBudget = {}, const Observer& observe = {})
{
  TrackFitResult fit{.ndf = static_cast<std::int16_t>((2 * n) - kNBoundParams)};
  const auto seed = makeChainSeed(prop, geometry, measurements, n);

  int outwardKinks = 0;
  std::int16_t outwardKinkMeasurement = -1;
  auto status = detail::runPass<Rule>(prop, geometry, material, measurements, n, PropagationDirection::kForward, seed.state, fit.outwardState, fit.outwardChi2, fit.failedMeasurement, kinkBudget, outwardKinks, outwardKinkMeasurement, observe);
  fit.outwardKinks = static_cast<std::uint8_t>(outwardKinks);
  if (status != PropagationStatus::kOk) {
    fit.status = status;
    fit.failedPass = FitPass::kOutward;
    return fit;
  }
  detail::inwardPass<Rule>(prop, geometry, material, measurements, n, seed.qOverPSigma, kinkBudget, observe, fit);
  return fit;
}

/// `fitMeasurements` given its outward pass's end state, \a outwardState
/// -- as the road fit keeps it for an accepted
/// road (`RoadFitResult::state`): only the inward pass runs. The outward
/// pass's chi2 and kinks are not known here and stay zero.
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Prop, DetectorLayout Det, typename Material, typename Chain, typename Observer = NoFitObserver>
GPUdi() TrackFitResult fitInward(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, const BoundTrackState& outwardState, const KinkBudget& kinkBudget = {}, const Observer& observe = {})
{
  TrackFitResult fit{.ndf = static_cast<std::int16_t>((2 * n) - kNBoundParams)};
  const auto seed = makeChainSeed(prop, geometry, measurements, n);
  fit.outwardState = outwardState;
  detail::inwardPass<Rule>(prop, geometry, material, measurements, n, seed.qOverPSigma, kinkBudget, observe, fit);
  return fit;
}

/// The road fit's candidate pass: the outward filter
/// alone, as in `fitMeasurements`.
struct RoadFitResult {
  Scalar chi2 = 0.F;
  RoadFitStatus status = RoadFitStatus::kAccepted;
  std::uint8_t kinks = 0;
  /// The outward pass's state at the outermost measurement: the pT floor
  /// reads it, and `fitMeasurements`' outward pass ends on it, so the refit
  /// of the road takes it instead of running that pass again.
  BoundTrackState state;
};

/// Its chi2 is the whole road's; the inward pass only gives the state on the
/// innermost surface, which the refit of an accepted track provides. The cut
/// is \a chi2MaxByNdf read at the road's degrees of freedom, `2 N - 5` less
/// two per kink. A filter's chi2 only grows as measurements are added, and the
/// cut only falls as kinks are taken, so the pass stops as soon as the chi2
/// exceeds the cut so far: the decision the complete pass would reach.
///
/// With \a LastStepAllowance, for an iteration with a kink allowance, the
/// last step is `allowanceStep` with the allowance's theta0 factor
/// \a lastStepScale (against the road-end fakes): refused
/// beyond the gate, it frees no degree of freedom, and the result's kinks
/// count it. A compile-time choice: without it the last step is a `kinkStep`
/// like the others.
template <bool LastStepAllowance = false, CrossingRule Rule = CrossingRule::kFirstAhead, typename Prop, DetectorLayout Det, typename Material, typename Chain>
GPUdi() RoadFitResult fitOutward(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, const Scalar* chi2MaxByNdf, const KinkBudget& kinkBudget = {}, [[maybe_unused]] Scalar lastStepScale = 1.F)
{
  const auto budget = detail::kinkBudgetFor(kinkBudget, n);
  BoundTrackState state = makeChainSeed(prop, geometry, measurements, n).state;
  CovarianceFactors factors;
  RoadFitResult result;
  // The start and the sum as `runPass` keeps them, so that the kink gate,
  // which reads its increments, rounds as in the refit's outward pass.
  const auto firstUpdate = detail::startOnMeasurement(state, factors, measurements[0]);
  if (!updateUsable(firstUpdate, factors)) {
    // Not a measurement of this road (a negative or NaN chi2 would pass the
    // cut below and rank first): counted with the propagation failures.
    countWork(WorkCounter::kRoadPropagation);
    result.status = RoadFitStatus::kPropagationFailed;
    return result;
  }
  result.chi2 = firstUpdate.chi2;
  const int ndf = (2 * n) - kNBoundParams;
  int kinks = 0;
  // The allowance's last step, which frees no degree of freedom; counted by its
  // instantiation alone.
  [[maybe_unused]] int allowed = 0; // NOLINT(misc-const-correctness)
  countWork(WorkCounter::kRoadStepsPossible, static_cast<std::uint64_t>(n - 1));
  // A road is at least a triplet, so `ndf` >= 1, and `kinkBudgetFor` keeps
  // `ndf - 2 kinks` there; the analyzer, not knowing n >= 3, takes n = 2.
  // NOLINTBEGIN(clang-analyzer-security.ArrayBound)
  for (int index = 1; index < n; ++index) {
    if (!(result.chi2 <= chi2MaxByNdf[ndf - (2 * kinks)])) {
      countWork(WorkCounter::kRoadChi2);
      result.status = RoadFitStatus::kRejected;
      return result;
    }
    countWork(WorkCounter::kRoadSteps);
    const auto& stepProp = stepPropagator(prop, geometry, measurements, index - 1, index);
    if constexpr (LastStepAllowance) {
      if (index == n - 1) {
        int taken = kinks;
        bool refused = false;
        const auto status = detail::allowanceStep<Rule>(stepProp, geometry, material, measurements[index], PropagationDirection::kForward, state, factors, result.chi2, budget, lastStepScale, taken, refused);
        if (status != PropagationStatus::kOk) {
          countWork(WorkCounter::kRoadPropagation);
          result.status = RoadFitStatus::kPropagationFailed;
          return result;
        }
        if (refused) {
          countWork(WorkCounter::kRoadChi2);
          result.status = RoadFitStatus::kRejected;
          return result;
        }
        allowed = taken - kinks;
        continue;
      }
    }
    const auto status = detail::kinkStep<Rule>(stepProp, geometry, material, measurements[index], PropagationDirection::kForward, state, factors, result.chi2, budget, index < n - 1, kinks);
    if (status != PropagationStatus::kOk) {
      countWork(WorkCounter::kRoadPropagation);
      result.status = RoadFitStatus::kPropagationFailed;
      return result;
    }
  }
  result.kinks = static_cast<std::uint8_t>(kinks + allowed);
  covarianceFromFactors(factors, state.cov);
  result.state = state;
  if (!(result.chi2 <= chi2MaxByNdf[ndf - (2 * kinks)])) {
    // NOLINTEND(clang-analyzer-security.ArrayBound)
    countWork(WorkCounter::kRoadChi2);
    result.status = RoadFitStatus::kRejected;
  }
  return result;
}

/// The outward road fit as an iteration of \a parameters takes it: its kink
/// budget, and with \a LastStepAllowance its kink allowance on the last step
/// (`fitOutward`'s), with the allowance's theta0 factor. The stages make the
/// choice at compile time (`RoadFitStage`).
template <bool LastStepAllowance, CrossingRule Rule = CrossingRule::kFirstAhead, typename Parameters, typename Prop, DetectorLayout Det, typename Material, typename Chain>
[[nodiscard]] GPUdi() RoadFitResult fitRoad(const Parameters& parameters, const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n)
{
  if constexpr (LastStepAllowance) {
    return fitOutward<true, Rule>(prop, geometry, material, measurements, n, parameters.roadChi2MaxByNdf, kinkBudget(parameters), o2::gpu::CAMath::Sqrt(parameters.kinkScale2));
  } else {
    return fitOutward<false, Rule>(prop, geometry, material, measurements, n, parameters.roadChi2MaxByNdf, kinkBudget(parameters));
  }
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRACKFITTER_H_

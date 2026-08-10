# Fast local triplet estimator

- Status: accepted fit direction; observation and cell-payload placement
  superseded by [global measurements and the triplet fit factor](0021-global-measurements-and-triplet-state.md)
- Date: 2026-08-10
- Predecessor: [momentum-dependent material experiment](0019-momentum-dependent-local-triplet-material.md)

## Goal

Neighbour finding needs a fast, approximately calibrated estimate of a
triplet's momentum and direction parameters and their covariance. It does not
need a self-consistent material fit or a precision energy-loss correction.

## Decision

The local fit is evaluated once with a fixed physical space-angle variance at
the middle observation. The caller derives that variance from the existing
`TrackletMinPt` scattering envelope prepared for the candidate graph. This is
conservative over the accepted momentum domain and avoids feeding the noisy
three-hit momentum estimate back into its own weight.

The fit therefore remains the closed-form operation

```text
fit(observations, sigmaMS2_at_TrackletMinPt) -> local parameters and covariance
```

with no PID, material traversal, convergence loop, or family dispatch. The
fitted curvature still provides the triplet momentum point estimate. Its
uncertainty includes full hit covariance and the supplied scattering process
variance.

The fixed-point implementation is retained at `ae27032797` for possible
future fits over longer triplet topologies, where the momentum lever arm may
justify self-consistent material weighting.

## Energy loss

Mean energy loss is not included in this local compatibility estimator. Over
one three-hit span it is a small deterministic change in momentum, while the
dominant material effect for direction compatibility is the stochastic
multiple-scattering kick. Applying energy loss here would require a particle
hypothesis and propagation direction, increase cost, and create stopping or
material-domain failures that must not become implicit neighbour rejections.

Energy loss remains part of the later state propagation/refit. If a future
longer-topology triplet fit transports momentum between several material
surfaces, it may add energy loss explicitly there.

## Global observation boundary

The successor design stores index-aligned global and native measurement
records. The cell retains the Equation-19 `Psi`, `rho`, and `H` factor rather
than independent segment marginals, so shared-hit covariance can be assembled
exactly when neighbouring triplets are compared.

## Remaining integration requirement

The factor contract is defined in the successor note. Its multi-triplet
Equation-19 assembly, MS model, candidate cut-flow, and physics calibration
must still be validated before changing neighbour selection.

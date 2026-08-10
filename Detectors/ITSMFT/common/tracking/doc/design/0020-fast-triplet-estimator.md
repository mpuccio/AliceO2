# Fast local triplet estimator

- Status: accepted standalone direction
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

`SurfaceMeasurement` already stores global `(x,y,z)`. Its three covariance
numbers intentionally remain in the measured two-dimensional surface frame;
that compact representation is sufficient because the normal coordinate is
exact at this stage. The descriptor-selected leaf lifts it once into the full
rank-deficient global covariance consumed by the neutral fit.

Replacing the stored covariance with a packed global `3x3` covariance would
increase every measurement by three floats and make the native cylinder
Kalman update project it back to `(u,v)`. It would move, rather than eliminate,
surface-frame work. Before neighbour integration, benchmark two narrower
placements:

1. cache the neutral global observation once per measurement; or
2. fit each accepted cell once and store only the compact triplet parameters
   and covariance needed by neighbour compatibility.

The second placement avoids repeating both observation construction and the
fit for every candidate neighbour and is the preferred direction. Surface
measurement and material-normal construction remain leaf operations; the fit
and neighbour comparison stay family-neutral.

## Remaining integration requirement

The current result exposes curvature, curvature variance, fit quality, and
the angular residual covariance. Before it replaces neighbour compatibility,
the compact cell payload must explicitly define the middle-reference
direction parameters and their covariance/cross-covariance with curvature.
That parameter contract, its comparison chi-square, and candidate-level cost
must be validated before changing a production selection.

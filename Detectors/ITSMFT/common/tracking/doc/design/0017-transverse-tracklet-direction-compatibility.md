# Transverse tracklet-direction compatibility

- Status: implemented candidate
- Date: 2026-08-10
- Production commit: `1e2b552ac6`
- Tests and guards: `f47b7aef18`
- Predecessor: [shared cell-direction compatibility](0016-cell-direction-compatibility.md)

## Decision

The combination of two tracklets now applies the configured
`TrackletMinPt` in the transverse plane without fitting the triplet. The
physical curvature interval is profiled first; measurement uncertainty and
multiple scattering then normalize only the observed bend outside that
interval. Cylinder and disk use the same algorithm and the common
`NSigmaCut` acceptance contract.

For ordered transverse observations `p0`, `p1`, and `p2`, let

```text
d01 = |p0-p1|,  d12 = |p1-p2|
kappaMax = |B2C Bz| / TrackletMinPt
kappa = min(kappaMax, 2/d01, 2/d12)
deltaPhiMax = asin(kappa d01 / 2) + asin(kappa d12 / 2).
```

The two `asin` terms are the half central angles subtended by the two chords
on one circle. The chord-domain bounds keep both terms defined without a
family guard. With the stored, inward-oriented tracklet chord azimuths,

```text
deltaPhi = remainder(phi01 - phi12, 2 pi)
e = max(0, |deltaPhi| - deltaPhiMax).
```

The unknown charge sign and any curvature magnitude allowed by
`TrackletMinPt` are nuisance quantities inside the symmetric interval
`[-deltaPhiMax,+deltaPhiMax]`; therefore they contribute no penalty. The
decision is

```text
chi2Phi = e^2 / Var(deltaPhi)
chi2Phi < NSigmaCut^2.
```

This is a cheap analytic gate. It does not estimate candidate momentum,
invoke propagation, or introduce a transverse resolution parameter.

## Uncertainty propagation

Each neutral transverse observation contains `x`, `y`, `Var(x)`,
`Cov(x,y)`, and `Var(y)`. For

```text
phi01 = atan2(y0-y1, x0-x1)
phi12 = atan2(y1-y2, x1-x2),
gij = (-(yi-yj)/dij^2, (xi-xj)/dij^2),
```

the derivatives of `deltaPhi=phi01-phi12` are

```text
J0 =  g01
J1 = -g01 - g12
J2 =  g12.
```

Thus the shared middle hit and its correlation between both chords are
retained exactly:

```text
Var_measurement(deltaPhi) = sum_i Ji Ci Ji^T.
```

The outgoing transition begins at the shared middle hit. Its already
prepared `TrackletMinPt` RMS scattering angle is an equivalent thin angular
kick in the transverse direction, so

```text
Var(deltaPhi) = Var_measurement(deltaPhi) + sigmaTheta^2.
```

There is no uncertainty floor. Inputs, covariances, chord lengths, and final
variance must be finite and valid; failures reject transactionally.

## Measurement leaves

The real decoder/state conventions define the only family-specific work.

- Disk measurements expose global `x/y`, and their packed covariance is
  `(Cxx,Cxy,Cyy)`. The transverse observation copies those values.
- Cylinder measurements expose tracking-frame `(q,u)` at frame angle
  `alpha`, with `q` normal to the sensor and treated as exact. The leaf uses

  ```text
  x = q cos(alpha) - u sin(alpha)
  y = q sin(alpha) + u cos(alpha)
  d(x,y)/du = (-sin(alpha), cos(alpha))
  ```

  to rotate `Cuu` into the rank-one global `x/y` covariance. The complete
  measured `(u,v)` covariance is still validated before projection.

The rank-one cylinder projection is accepted within machine roundoff only;
materially non-positive-semidefinite inputs still fail closed.

## Orchestration boundary

`TrackerTraits` resolves the same three measurements and descriptors used by
cell construction, invokes the descriptor-selected transverse-observation
leaf three times, and makes one family-neutral compatibility call. It passes
the two stored tracklet azimuths, `Bz`, `TrackletMinPt`, the outgoing
transition process noise, and `NSigmaCut`. There is no cell-direction
`SurfaceKind`, detector-ID, or source-ID dispatch.

The existing longitudinal `(r,z)` collinearity gate follows the transverse
gate and consumes the same process-noise record. Graph, holes, `(6,7)`,
timing, neighbour/road logic, scalar policy, and legacy workflows are
unchanged.

## Outlook

The analytic gate establishes the cheap `TrackletMinPt`, `deltaPhi`,
measurement-covariance, and MS contract. The detector-neutral triplet fit
remains the next R&D step if its added parameter estimate and covariance are
worth the measured cost. It is not required to restore candidates below the
configured transverse-momentum threshold.

# Shared cell-direction compatibility

- Status: implemented candidate
- Date: 2026-08-10
- Production commit: `ec9deeb178`
- Tests and guards: `58d87eb8e2`
- Targeted cylinder observation correction: `1e2c8634be`
- Targeted correction tests: `4435cfe08f`
- Predecessor: [tracklet and cell harmonization](0015-tracklet-cell-harmonization.md)

## Decision

Cylinder and disk cell construction use one uncertainty-normalized
collinearity test in the physical `(r,z)` plane. For three ordered
observations,

```text
K = (z0 - z1) (r1 - r2) - (z1 - z2) (r0 - r1)
```

The decision is

```text
K^2 / (Var_measurement(K) + Var_process(K)) < NSigmaCut^2
```

There is no slope division, fixed direction resolution, family-specific
threshold, or uncertainty floor. `Tracklet::tanLambda` remains `dz/dr`, but
cell direction no longer compares stored tracklet slopes.

The shared algorithm accepts three neutral observations containing `r`, `z`,
`Var(r)`, `Cov(r,z)`, and `Var(z)`. Its derivatives for observation
`i = 0,1,2` are

```text
dK/d(r,z)[0] = (-(z1-z2),        r1-r2)
dK/d(r,z)[1] = ( (z0-z1)+(z1-z2), -((r0-r1)+(r1-r2)))
dK/d(r,z)[2] = (-(z0-z1),        r0-r1)
```

`Var_measurement(K)` is the sum of each derivative pair contracted with that
observation's full symmetric `(r,z)` covariance.

The process term models the outgoing transition's already-prepared multiple
scattering as an equivalent thin angular kick at the middle observation. For
the incoming and outgoing segment vectors `a=(r1-r0,z1-z0)` and
`b=(r2-r1,z2-z1)`, with `L=|b|` and normal `n=(-b_z,b_r)/L`, the induced
outer-point covariance is

```text
Q_process = L^2 sigmaTheta^2 n n^T
Var_process(K) = J2 Q_process J2^T
               = (a dot b)^2 sigmaTheta^2.
```

There is no slope denominator and no spatial uncertainty floor. The neutral
`sigmaTheta^2` input is the square of the outgoing transition RMS scattering
angle already computed from nominal material and `TrackletMinPt` during
transition preparation. A skipped outgoing transition is represented by its
quadrature-summed angle at the middle observation; this is conservative
because it gives downstream thin scatterers the maximum outgoing lever arm.
The cell loop neither estimates momentum nor repeats family material math.

Inputs must be finite, angular variance must be non-negative, each measurement
covariance must be positive semidefinite, `r` must be positive, and the final
variance must be finite and positive. Failure is transactional and rejects the
candidate.

## Measurement models

The real decoder and state-operation conventions determine the leaf
conversion.

- Disk decoding maps the two measured coordinates to global `x/y`; its
  covariance is therefore `(Cxx,Cxy,Cyy)` at this boundary. The leaf computes
  `r = hypot(x,y)` and

  ```text
  Var(r) = (x^2 Cxx + 2xy Cxy + y^2 Cyy) / (x^2 + y^2).
  ```

  Disk `z` is the descriptor's fixed reference coordinate, expressed as
  `Var(z)=Cov(r,z)=0`.

- Cylinder decoding stores covariance in tracking-frame `(u,v)=(y,z)`, with
  `q` normal to the sensor plane. The tracking frame differs from global
  coordinates by a rotation around Z, so the leaf computes

  ```text
  r = hypot(q,u),  z = v,  dr/du = u/r
  Var(r) = (dr/du)^2 Cuu
  Cov(r,z) = (dr/du) Cuv
  Var(z) = Cvv.
  ```

  The normal coordinate `q` is treated as exact. The descriptor radius is a
  nominal layer reference, not the hit-level physical radius, and therefore
  is not used as the observation coordinate.

## Code boundary

`TrackerTraits` resolves the three measurements and descriptors, asks the
descriptor-selected observation leaf to construct each neutral observation,
binds the already-prepared outgoing transition angle as neutral process noise,
and makes one shared compatibility call. It contains no family branch for the
direction decision. Existing surface selection remains only for cell-seed,
propagation, and transition-material leaves.

`CellFinding.h` and `TrackletFinding.h` were consolidated into
`CandidateFinding.h`, with their implementation in `CandidateFinding.cxx`.
That boundary contains the cohesive candidate pipeline: tracklet projection
and search, neutral cell-direction observations, the shared statistic, cell
seed construction, attachment, extension, and transition preparation.

The common `CellDeltaTanLambdaSigma` / `cellDeltaTanLambdaSigma` parameter and
its configuration projection were removed. Legacy ITS production tracking
retains its own parameter and implementation outside this common-CA slice.

## Transverse extension

This correction stopped at longitudinal `(r,z)` compatibility. The subsequent
[transverse tracklet-direction compatibility](0017-transverse-tracklet-direction-compatibility.md)
adds the cheap analytic `TrackletMinPt`, `deltaPhi`, measurement-covariance,
and MS contract. The full detector-neutral triplet fit remains a later R&D
migration.

## Deferred deletion inventory

The existing MFT-only common-tracker helpers and publication compatibility
facilities remain on the previously recorded later-deletion inventory. This
change does not add a new MFT adapter, policy, parameter, or compatibility
path. Legacy workflow removal, neighbour/road harmonization, holes, and graph
changes remain outside this slice.

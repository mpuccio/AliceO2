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
K^2 / Var(K) < NSigmaCut^2
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

`Var(K)` is the sum of each derivative pair contracted with that
observation's full symmetric `(r,z)` covariance. Inputs must be finite, each
covariance must be positive semidefinite, `r` must be positive, and the final
variance must be finite and positive. Failure is transactional and rejects
the candidate.

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
and makes one shared compatibility call. It contains no family branch for
the direction decision. Existing surface selection remains only for cell-seed
and propagation leaves.

`CellFinding.h` and `TrackletFinding.h` were consolidated into
`CandidateFinding.h`, with their implementation in `CandidateFinding.cxx`.
That boundary contains the cohesive candidate pipeline: tracklet projection
and search, neutral cell-direction observations, the shared statistic, cell
seed construction, attachment, extension, and transition preparation.

The common `CellDeltaTanLambdaSigma` / `cellDeltaTanLambdaSigma` parameter and
its configuration projection were removed. Legacy ITS production tracking
retains its own parameter and implementation outside this common-CA slice.

## Deferred deletion inventory

The existing MFT-only common-tracker helpers and publication compatibility
facilities remain on the previously recorded later-deletion inventory. This
change does not add a new MFT adapter, policy, parameter, or compatibility
path. Legacy workflow removal, neighbour/road harmonization, holes, and graph
changes remain outside this slice.

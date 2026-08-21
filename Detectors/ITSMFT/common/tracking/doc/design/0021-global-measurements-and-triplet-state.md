# Global measurements and the triplet fit factor

- Status: implementation design
- Date: 2026-08-10
- Predecessor: [fast local triplet estimator](0020-fast-triplet-estimator.md)

## Event measurement split

`MultiSourceFrame`, owned by `TimeFrame`, stores two index-aligned arrays per
surface. A `(SurfaceId, SurfaceMeasurementIndex)` resolves the same hit in
both arrays.

`GlobalMeasurement` contains tracking-candidate geometry and hit identity:

```text
(x, y, z, r)
packed global Cov(x,y,z)
sensor, cluster, shape, source ROF, surface, flags
```

`r = hypot(x,y)` is a validated cache. It is constructed from `(x,y)` during
decoding and is never accepted as an independent input coordinate.

`SurfaceMeasurement` contains only the native measurement used by a surface
state update:

```text
(q, u, v, frameAngle)
packed Cov(u,v)
```

The decoder constructs both records atomically. For cylinders it rotates the
measured `(u,v)` covariance into global coordinates with exact `q`; for disks
it copies the established global `(x,y)` covariance and expresses exact `z`
with zero variance. Global covariance validity is checked once at loading.

Tracklet search, cell direction compatibility, and triplet-factor construction
consume `GlobalMeasurement`. Publication and label lookup use hit identity from
`GlobalMeasurement`. Family-specific covariance projection therefore ends at
decoding rather than being repeated in candidate loops.

`CellSeed` is the CA cell and geometric triplet. It contains the three hit
references, two tracklet references, time/level metadata, and the linearized
triplet-fit factor described below. It does not contain a
`SurfaceKinematicState`, `q/pT`, or a fit `chi2`.

`TrackerTraits::buildTrackSeed` is the single fallible transition from a
selected `CellSeed` to the state-bearing `TrackSeed`. It resolves the three
global/native measurements, anchors one state at the outer measurement, fills
the native parameter/covariance values, and runs one common attachment loop
over the middle and inner measurements. One constant-solenoid construction
derives curvature, `q/pT`, averaged `tanLambda`, and the outer direction from
the three global positions. The last-mile cylinder/disk branch stores that
direction as `snp` or `phi` and applies the corresponding angular covariance.
Hit attachment and refit then consume only native measurements.

The former `buildCellSeed` facade and the family-local
`detail::barrel::buildSeed`/`detail::forward::buildSeed` leaves are retired;
there is no second seed-construction path alongside this boundary.

Rank-deficient global covariances can acquire a negative projected variance
at float roundoff. The neutral `(r,z)` projection validates the global matrix,
then repairs only violations within a term-scaled float roundoff tolerance by
clamping to the nearest positive-semidefinite boundary. Larger violations
still fail closed. No uncertainty floor is introduced.

## Triplet fit factor

Cell construction produces only the linearized factor needed by Eq. (19) of
the General Triplet Track Fit:

```text
Psi = (thetaTilde, phiTilde)
rho = (rhoTheta, rhoPhi)
H   = d(Psi + rho*kappa) / d(x0,y0,z0,...,x2,y2,z2)
```

`H` is evaluated at the geometric reference
`kappaRef = -Psi_phi/rho_phi`. The 2x9 matrix is stored as two xyz gradients
for each of the three hits. Float conversion is checked and fails closed.
The earlier local momentum and segment-error estimator is not part of this
path: its outputs are not consumed by production. `CellSeed` stores only the
geometry/Jacobian factor used by adjacent-cell compatibility; track-state
construction remains independent in `TrackerTraits::buildTrackSeed`.

The cell already stores the three local cluster indices and their positional
hit mask. `getClusterReference(i)` exposes their unambiguous
`(surfacePosition, clusterIndex)` binding; the factor does not duplicate
those references. Factor hit slot `i` uses the same ordering.

For adjacent triplets, neighbour finding resolves the global hit covariance
through these references and assembles the four-kink system

```text
C = D_MS + H V_hit H^T,
chi2 = min_kappa (Psi + rho*kappa)^T C^-1 (Psi + rho*kappa).
```

Here `Psi`, `rho`, and `H` concatenate the `(theta,phi)` rows of both
triplets. Shared hits naturally populate the off-diagonal blocks of `C`; their
uncertainty is neither discarded nor counted as independent twice. The MS
term remains outside the stored factor because it depends on the material
and momentum hypothesis used by the compatibility operation.

## Material and reference convention

For each triplet the caller supplies its second transition's space-angle
variance evaluated at `TrackletMinPt`. It enters the angular covariance as
`Var(theta)=deltaTheta^2` and
`Var(phi)=deltaTheta^2/sin(theta)^2`. There is no energy-loss correction or
momentum/material iteration. `tanLambda` is the longitudinal displacement
divided by the transverse arc length of the corresponding fitted segment.

The factor is evaluated once when a cell is accepted. It depends only on the
three global observations; magnetic field and material do not enter factor
construction. Its 88-byte float payload stores `Psi`, `rho`, and all 18
elements of `H`. Invalid factor construction rejects the cell. The local fit
then fails closed if its resolved covariance or material term is invalid.
Neither triplet linearization nor state construction is repeated per
prospective neighbour within a road step.

Neighbour finding visits the global scheduled-cell span once. It verifies the
two shared hit references, resolves the four unique global observations,
assembles their joint covariance, solves the symmetric positive-definite 4x4
system by Cholesky decomposition, and applies the existing common
`maxChi2ClusterAttachment` threshold. There is no surface-kind, detector, or
source dispatch in this path. Cylinder/disk distinctions end when the decoder
constructs `GlobalMeasurement`.

The host characterization test evaluates 20,000 factors. `CellSeed` no longer
embeds the native state and `chi2`; the 88-byte factor is its only fit payload.
`TrackSeed` remains the first state-bearing road object. The normalized
per-hit storage uses a global record plus a native record, keeping candidate
geometry independent from state-update coordinates.

## Closed-form adjacent fit

With `W=C^-1`, the common-curvature solution is

```text
kappa = -(rho^T W Psi) / (rho^T W rho)
Var(kappa) = 1 / (rho^T W rho)
chi2 = Psi^T W Psi - (rho^T W Psi)^2 / (rho^T W rho).
```

No slope denominator is formed. Non-finite observations, covariance outside
the float-storage PSD tolerance, non-positive factor covariance, and invalid
solutions fail closed. The validation and population characterization are in
[adjacent triplet neighbour fit](../validation/adjacent-triplet-neighbour-fit.md).

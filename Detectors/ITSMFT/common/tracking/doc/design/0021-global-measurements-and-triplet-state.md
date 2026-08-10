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

Tracklet search, cell direction compatibility, and triplet fitting consume
only `GlobalMeasurement`. Seed construction, hit attachment, and native refit
consume `SurfaceMeasurement`. Publication and label lookup use hit identity
from `GlobalMeasurement`. Family-specific covariance projection therefore
ends at decoding rather than being repeated in candidate loops.

Cylindrical seed initialization remains a geometric exception to the native
update rule: its established three-point formula takes the inner and middle
global positions and anchors the state in the outer native frame. Rebuilding
those two points from float surface-frame values changes the seed. The cell
leaf therefore receives the global and native records explicitly; subsequent
propagation and updates consume only the native record.

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
path: its outputs were discarded by production, so retaining it would make
every cell pay for covariance contractions and automatic-differentiation
propagation that neighbour finding does not consume. That prototype remains
available in the repository history.

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
elements of `H`. An invalid factor does not reject the accepted cell during
cell construction, but it fails closed if the cell is later considered for a
neighbour. Neither global observation construction nor triplet linearization
is repeated per prospective neighbour.

Neighbour finding visits the global scheduled-cell span once. It verifies the
two shared hit references, resolves the four unique global observations,
assembles their joint covariance, solves the symmetric positive-definite 4x4
system by Cholesky decomposition, and applies the existing common
`maxChi2ClusterAttachment` threshold. There is no surface-kind, detector, or
source dispatch in this path. Cylinder/disk distinctions end when the decoder
constructs `GlobalMeasurement`.

The host characterization test evaluates 20,000 factors and reports about
0.24 microseconds per factor on the validation machine. `CellSeed` grows from
132 to 220 bytes for the 88-byte factor. The normalized per-hit storage grows
from one 72-byte mixed record to a 72-byte global record plus a 28-byte native
record. This deliberately spends 28 bytes per hit to remove repeated global
projection and keep the candidate geometry independent from state-update
coordinates.

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

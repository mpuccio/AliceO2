# Global measurements and triplet-factor validation

- Date: 2026-08-10
- Package: `daily-20260717-0700-local1`
- Branch: `codex/itsmft-triplet-rnd-scratch`
- Active parent: `056f955cb6`
- Iterative prototype retained at: `ae27032797`
- Verdict: safe structural candidate; neighbour selection not activated
- Subsequent activation: [adjacent-triplet neighbour fit](adjacent-triplet-neighbour-fit.md)

## Implemented boundary

`TimeFrame` now owns index-aligned global and native measurement arrays for
every surface. `GlobalMeasurement` carries `(x,y,z)`, cached
`r = hypot(x,y)`, the packed global covariance, and hit identity. The cache is
recomputed and validated at loading. `SurfaceMeasurement` carries only the
native `(q,u,v,frameAngle)` coordinates and packed `(u,v)` covariance.

Tracklet search, transverse and longitudinal cell compatibility, and the fast
triplet estimator consume the global record. Native seed/update/refit leaves
consume the surface record. The cylinder seed keeps its established inner and
middle global points while using the outer native frame as its anchor. No
detector/source dispatch was added to candidate geometry or the triplet fit.

Each accepted cell receives an 88-byte Equation-19 factor containing
`Psi = (thetaTilde,phiTilde)`, `rho = (rhoTheta,rhoPhi)`, and the full 2x9 hit
Jacobian `H`. The factor's three hit slots bind to the cell's existing
positional mask and local cluster-index array; no identity is duplicated.
The local fit still exposes segment directions, signed `q/pT`, `pT`, and their
covariances as diagnostics, but those marginals are not persisted in the
cell.

The factor deliberately excludes the MS covariance. The future neighbour
operation must construct `D_MS` from its material/momentum hypothesis and
resolve `V_hit` from the referenced global measurements before assembling
`C = D_MS + H V_hit H^T`. This retains the off-diagonal covariance induced by
the two hits shared by adjacent triplets.

The estimator is non-iterative. It uses the middle-transition space-angle
variance evaluated at `TrackletMinPt`, includes full hit covariance, and does
not apply energy loss. Energy loss remains in propagation/refit. An invalid
triplet estimate leaves the payload invalid but does not reject an otherwise
accepted cell. Neighbour finding does not read the new payload in this slice.

## Numerical boundary correction

An initial replay produced 185 rather than 184 ITS tracks. Stage diagnostics
found unchanged tracklets and compatibility decisions but one fewer cylinder
cell: the seed formula had inadvertently reconstructed inner/middle points
from native frame values. Restoring the original global geometric inputs
returned the count to 184.

One same-label outer-hit substitution remained. The absent cell used cluster
1946, whose valid rank-deficient global covariance projected to
`Var(r) = -7.3864182144965215e-15` solely through float cancellation; the
physical projected value is zero and `Cov(r,z)` is zero. The strict projected
PSD check rejected it before cell construction. The projection now clamps
only violations within a term-scaled float-roundoff tolerance to the nearest
PSD boundary. A materially non-PSD or non-finite input still fails closed.
This restored the complete ITS reference population and hash.

## Cost and representation

The 20,000-repetition host characterization reports 489.3 ns per triplet fit
on the validation machine. The sizes in the pinned build are:

| Type/storage | Bytes |
|---|---:|
| `GlobalMeasurement` | 72 |
| `SurfaceMeasurement` | 28 |
| previous mixed per-hit record | 72 |
| current combined per-hit storage | 100 |
| `TripletFitFactor` | 88 |
| `CellSeed` before/after payload | 132 / 220 |
| `TrackSeed` | 248 |

The per-hit cost is +28 bytes and each accepted cell costs +88 bytes. The fit
is performed once per accepted cell rather than once per prospective
neighbour.

## Validation campaign

Strict preflight resolved the exact package, found a valid AliEn token from
2026-08-07 05:11:10 UTC through 2026-09-05 09:33:43 UTC, confirmed
`ALICEO2_CCDB_NOTOKENCHECK` was unset, and found every required Geant4
dataset. The fixed fixture's 43-entry checksum manifest passed before and
after the campaign.

The complete serial suite passed:

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 94
Total Test time (real) = 93.44 sec
```

`git diff --check` is clean and pinned `git clang-format --diff HEAD` reports
no changes.

Final replay characterization is identical to the parent candidate:

| Product | Tracks | Content hash | Matched / reconstructable | Efficiency | Fake | Clone |
|---|---:|---|---:|---:|---:|---:|
| ITS standalone | 184 | `932e38144a0575f495b29b532012f468` | 118 / 142 | 0.830986 | 15 (0.0815217) | 1 (0.00704225) |
| ITS combined | 184 | `932e38144a0575f495b29b532012f468` | 118 / 142 | 0.830986 | 15 (0.0815217) | 1 (0.00704225) |
| MFT standalone | 69 | `f6dee3f7a5f7def6b55900dbac734ef0` | 57 / 109 | 0.522936 | 2 (0.0289855) | 0 |
| MFT combined | 100 | `98f9730d6fca9e738c3f20afc66d296d` | 78 / 109 | 0.715596 | 3 (0.03) | 6 (0.0550459) |

The persisted-product comparator reports field-level equality for the
preceding/current ITS pairs and the pre-correction/current MFT pairs,
including cluster references, ROFs, seed patterns, and labels. The unchanged
pre-correction MFT products had already matched the preceding candidate. The
final MFT comparisons cover 3,036 standalone and 4,400 combined forward float
values with maximum absolute and relative delta zero.

The standalone/combined MFT population difference is unchanged by this
migration. It remains the already-characterized composition effect of their
different application-level scalar policies, not a result of the global/native
split or triplet payload.

## Equation-19 factor refinement

After the campaign above, the unconsumed cell payload was changed from the
two segment marginals to `Psi`, `rho`, and `H`. Three focused tests pass for
the factor numerics, cell/reference representation, and cylinder/disk cell
orchestration. Recontracting the stored float `H` with the full non-diagonal
hit covariances reproduces all three measurement-covariance components of the
double-precision local fit within `2e-6` relative tolerance. The measured
host cost is 482.4 ns per fit, compared with 489.3 ns before the payload
change. Neighbour selection still does not read the payload, so this
refinement introduces no physics selection.

After rebuilding the complete common-tracking test subtree, the full pinned
serial `itsmft` suite passed again: 94/94 tests in 116.57 seconds.

## Outlook

The next step is to assemble the joint Equation-19 system for adjacent
triplets from their factors, referenced global hit covariances, and an
explicit MS model. It can then fit their common curvature and evaluate the
compatibility chi-square without an independence approximation. Candidate
cut-flow and ITS/MFT physics calibration are required before activation. No
neighbour/road or hole-policy change belongs to this slice.

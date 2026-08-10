# Shared cell-direction compatibility validation

- Status: safe targeted-correction candidate; unified physics sign-off required
- Date: 2026-08-10
- Harmonization commits: `ec9deeb178`, `58d87eb8e2`, `9fdeaaf318`
- Targeted cylinder correction: `1e2c8634be`
- Targeted correction tests: `4435cfe08f`
- Design: [shared cell-direction compatibility](../design/0016-cell-direction-compatibility.md)
- Rejected nominal-radius population inventory: [cell-direction population delta](cell-direction-population-delta.md)
- Pinned O2 package: `daily-20260717-0700-local1`
- Build: `/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`
- Fresh correction artifacts: `/private/tmp/o2-cell-direction-cylinder-correction-20260810-agvm4P`

## Correction and covariance contract

The first harmonized implementation treated the nominal cylinder descriptor
radius as the exact hit radius. Real ITS sensors and hits do not lie exactly
at that nominal layer radius. The old stored `tanLambda=dz/dr` instead used
the per-hit global radius, so combining nominal radii with micrometre-scale Z
uncertainties created artificial pulls as large as `chi2=136644.5` and reduced
the ITS fixture from 212 to 17 tracks.

The targeted correction changes only the cylinder observation leaf. ITS
decoding supplies the sensor tracking frame `(q,u,v)`, where `(u,v)=(y,z)`
are measured and `q` is normal to the sensor plane. That frame differs from
global coordinates by a rotation around Z, which preserves radius. The leaf
therefore constructs

```text
r = hypot(q,u),  z = v,  dr/du = u/r
Var(r) = (dr/du)^2 Cuu
Cov(r,z) = (dr/du) Cuv
Var(z) = Cvv.
```

`q` is the exact surface-normal coordinate. No uncertainty floor, fixed
direction scale, family cut, detector dispatch, or generic family branch was
added. The shared compatibility algorithm remains

```text
K = (z0-z1)(r1-r2) - (z1-z2)(r0-r1)
chi2 = K^2 / Var(K)
chi2 < NSigmaCut^2.
```

Disk observation construction and the common `NSigmaCut=5` contract are
unchanged.

## Tests and source boundary

The focused cylinder test now proves a nontrivial tracking-frame projection:
for `(q,u)=(3,4)`, `(Cuu,Cuv,Cvv)=(0.04,0.006,0.09)`, the leaf returns
`r=5`, `Var(r)=0.0256`, `Cov(r,z)=0.0048`, and `Var(z)=0.09`. Zero radius and
non-PSD cylinder covariance fail without mutating the output.

The combined-composition synthetic helix fixtures were also corrected to
populate their cylinder tracking frames. Their former all-zero `q/u/v`
values were not valid decoded cylinder measurements; the tests continue to
exercise ROF ownership and explicit ITS-then-MFT orchestration.

The final clean serial campaign passed:

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 93
Total Test time (real) = 93.77 sec
```

The source guards still prove that `CellFinding.h`, the retired disk-specific
compatibility function, a cell-direction `SurfaceKind` switch, and generic
detector/source dispatch are absent. `TrackerTraits` still constructs three
neutral observations and makes one shared compatibility call.

Strict preflight passed with a valid AliEn token, no
`ALICEO2_CCDB_NOTOKENCHECK` bypass, the exact pinned package, and all Geant4
datasets. The fixture checksum manifest passed 43/43 before and 43/43 after
the campaign. Replays used condition cap `1784207296000`.

## Stage-level effect

Temporary diagnostic logging was removed before the final clean rebuild.
Tracklet formation remains unchanged; the first population decision is cell
creation.

| ITS stage | Prior fixed slope pull | Nominal-radius collinearity | Corrected hit-radius collinearity |
|---|---:|---:|---:|
| Accepted tracklets | 1,886 | 1,886 | 1,886 |
| Accepted cells | 1,455 | 136 | 770 |
| Neighbour links | 1,049 | 29 | 364 |
| Accepted/published tracks | 212 | 17 | 174 |

The correction recovers 634 cells, 335 neighbour links, and 157 tracks from
the rejected nominal-radius candidate. The remaining 38-track difference
from the former fixed slope pull is the direct and downstream effect of the
requested measurement-uncertainty-normalized collinearity selection. It is
not recoverable by correcting the radius mapping alone.

## Numerical representative cases

The common acceptance threshold is `chi2 < 25`. `oldPull` characterizes the
retired `abs(delta tanLambda)/0.007` decision and is not used in production.
Each observation below is `(r,z,Var(r),Cov(r,z),Var(z))`.

| Candidate refs | Old slope pull | Corrected observations | `K` | `Var(K)` | chi2 | Decision |
|---|---:|---|---:|---:|---:|---|
| `117,327,515` | `0.0435625` | `(2.235743,-4.604550,4.43e-11,0,5.12e-7)`, `(3.203698,-6.598718,2.84e-8,0,3.79e-7)`, `(3.864410,-7.959709,7.86e-9,0,3.80e-7)` | `-1.94458e-4` | `1.93502e-6` | `0.01954` | accept |
| `639,812,1010` | `0.0458317` | `(19.739592,-7.349940,6.96e-10,0,3.79e-7)`, `(24.676400,-9.196802,1.90e-9,0,3.24e-6)`, `(34.507762,-12.877872,5.65e-11,0,3.79e-7)` | `0.0155680` | `7.53472e-4` | `0.32166` | accept |
| `470,639,812` | `0.346856` | `(3.795571,-1.423988,5.90e-9,0,5.00e-7)`, `(19.739592,-7.349940,6.96e-10,0,3.79e-7)`, `(24.676400,-9.196802,1.90e-9,0,3.24e-6)` | `0.191120` | `1.00240e-3` | `36.4397` | reject |

The first row is the clearest proof of the defect correction: the same
candidate changed from an artificial `chi2=136644.5` under nominal layer
radii to `0.01954` using its physical hit radii. The last row illustrates a
remaining intentional selection difference: its coordinates are valid, but
the measured triplet is incompatible with exact `(r,z)` collinearity at five
measurement sigmas. Any future recovery of such candidates would require an
explicit physical extension of `Var(K)` for curvature or process uncertainty,
not a larger family-specific magic cut.

## Population characterization

Standalone and combined ITS outputs are field-identical for the correction.
The clean standalone repeat is also TSV-identical, including ordered
references, ROFs, labels, fake bits, chi2, and hit counts.

| ITS candidate | Tracks / hash | Matched / 142 | Efficiency | Fake / rate | Clone | Hit distribution |
|---|---:|---:|---:|---:|---:|---|
| Prior fixed slope pull | 212 / `46913a67a7e2fe7462e29df0db264fa8` | 135 | `0.950704` | 16 / `0.0754717` | 0 | 46x4, 16x5, 29x6, 121x7 |
| Rejected nominal-radius candidate | 17 / `b194927a8306cfac57ada54fc55c0c6f` | 10 | `0.0704225` | 3 / `0.176471` | 0 | 13x4, 3x5, 1x7 |
| Corrected hit-radius candidate | 174 / `5103bc6679a801339828eb7e684e4ba0` | 112 | `0.788732` | 13 / `0.0747126` | 0 | 135x4, 20x5, 10x6, 9x7 |

Against the prior 212 ordered reference sequences, the corrected population
has 39 exact matches, 173 removals, and 135 additions. The low exact-overlap
count reflects changed road availability and candidate competition; the
stage counts localize the initiating difference to cell creation. The
complete 212-to-17 appendix remains useful only as characterization of the
rejected nominal-radius implementation and is not the corrected population.

MFT is field-identical to the already reviewed disk-collinearity products:

| MFT composition | Tracks / hash | Matched / 109 | Efficiency | Fake / rate | Clone / rate |
|---|---:|---:|---:|---:|---:|
| Standalone | 69 / `f6dee3f7a5f7def6b55900dbac734ef0` | 57 | `0.522936` | 2 / `0.0289855` | 0 / 0 |
| Combined | 100 / `98f9730d6fca9e738c3f20afc66d296d` | 78 | `0.715596` | 3 / `0.03` | 6 / `0.0550459` |

## Verdict

**Safe targeted-correction candidate; unified physics sign-off required.**
The catastrophic 212-to-17 loss was caused by an incorrect cylinder
observation coordinate, not by the common five-sigma threshold. The
correction restores the physical hit radius and propagates every available
tracking-frame covariance term without adding tuning or family dispatch.
The deterministic 212-to-174 difference is now the visible candidate-physics
effect of measurement-only `(r,z)` collinearity. It must not be tuned away in
this slice; whether curvature or process uncertainty belongs in the common
variance model is a separate physics decision.

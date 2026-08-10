# Shared cell-direction compatibility validation

- Status: TrackletMinPt process-noise correction validated; transverse step deferred
- Date: 2026-08-10
- Harmonization commits: `ec9deeb178`, `58d87eb8e2`, `9fdeaaf318`
- Cylinder observation correction: `1e2c8634be`, `4435cfe08f`, `09f821f1ce`
- Design: [shared cell-direction compatibility](../design/0016-cell-direction-compatibility.md)
- Pinned O2 package: `daily-20260717-0700-local1`
- Build: `/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`
- Fresh artifacts: `/private/tmp/o2-cell-direction-ms-20260810-JlTxV5`

## Targeted correction

The measurement-only statistic was physically incomplete because it treated a
charged particle as exactly straight between three measurements. The corrected
common statistic remains

```text
K = (z0-z1)(r1-r2) - (z1-z2)(r0-r1)
chi2 = K^2 / Var(K)
chi2 < NSigmaCut^2,
```

but now

```text
Var(K) = Var_measurement(K) + (a dot b)^2 sigmaTheta^2,
a = (r1-r0,z1-z0), b = (r2-r1,z2-z1).
```

This is the contraction of the analytic K derivative with the spatial
covariance induced by a middle-point angular kick,
`Q=L^2 sigmaTheta^2 n n^T`. The neutral `sigmaTheta^2` is the square of the
already-prepared outgoing-transition RMS multiple-scattering angle. That angle
uses nominal surface material and `TrackletMinPt`; no candidate momentum fit,
new parameter, uncertainty floor, or family-specific direction cut was added.
For a skipped outgoing transition, its quadrature-summed scattering is placed
at the middle observation, conservatively giving downstream material the
maximum lever arm.

`TrackerTraits` still constructs three descriptor-selected neutral
observations and makes one family-neutral compatibility call. It reads the
precomputed angle by transition slot and contains no cell-direction
`SurfaceKind`, detector-ID, or source-ID dispatch. The stashed triplet fitter,
transverse curvature, `deltaPhi`, graph, holes, `(6,7)` disconnection,
neighbour/road logic, timing, and scalar policy are unchanged.

## Focused and full tests

The analytic test uses non-axis-aligned segments with `a dot b=18` and
`sigmaTheta^2=0.01`, proving that the process contribution is
`18^2*0.01=3.24` and that the total variance changes from `6.55` to `9.79`.
Negative and non-finite angular variance fail transactionally. Existing
measurement covariance, rotation, disk/cylinder observation, collinearity,
common-threshold, and source-boundary tests remain active.

A new orchestration test runs the same non-collinear cylinder candidate at
two values of `TrackletMinPt`: it is accepted at `0.3 GeV/c` and rejected at
`1 GeV/c`. This proves that the cell loop consumes the precomputed
TrackletMinPt scattering term without a fit or family branch.

The final pinned serial campaign passed:

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 93
Total Test time (real) = 89.31 sec
```

`git diff --check` and `git clang-format --diff` are clean. Strict preflight
passed with a valid AliEn token, no `ALICEO2_CCDB_NOTOKENCHECK` bypass, the
exact pinned package, and every Geant4 dataset. The fixture checksum manifest
passed 43/43 before and 43/43 after the campaign. All temporary diagnostic
instrumentation was removed before the final clean build and tests.

## Stage-level cut flow

Tracklet formation is unchanged. The process term acts first at cell
construction.

| ITS stage | Fixed slope pull | Measurement-only K | K plus TrackletMinPt MS |
|---|---:|---:|---:|
| Accepted tracklets | 1,886 | 1,886 | 1,886 |
| Accepted cells | 1,455 | 770 | 1,489 |
| Neighbour links | 1,049 | 364 | 1,050 |
| Accepted/published tracks | 212 | 174 | 212 |

The correction recovers 719 cells, 686 neighbour links, and 38 tracks from
the measurement-only candidate. Relative to the former fixed pull it exposes
34 additional cells (`+2.34%`), one additional neighbour link, and the same
published count. This is recovery through physical process uncertainty, not
through restoration of the fixed `0.007` scale.

| Standalone MFT stage | Measurement-only K | K plus TrackletMinPt MS |
|---|---:|---:|
| Accepted tracklets | 2,723 | 2,723 |
| Accepted cells | 601 | 604 |
| Neighbour links | 314 | 314 |
| Accepted/published tracks | 69 | 69 |

The three added disk cells do not create a new neighbour or output track.

## Numerical representative

For the previously rejected ITS candidate with references `470,639,812`,
the corrected hit observations give

```text
K = 0.1911203284
Var_measurement(K) = 0.00100239554
measurement-only chi2 = 36.4397  (reject at 25)
```

Its middle surface has `x/X0=0.01`. With `TrackletMinPt=0.3 GeV/c` and the
established pion-mass transition preparation,

```text
sigmaTheta = 0.00412722 rad
a dot b = 89.6570 cm^2
Var_process(K) = 0.136925
Var_total(K) = 0.137928
corrected chi2 = 0.264827  (accept)
```

The large change is expected: micrometre-scale measurement errors alone
cannot represent milliradian multiple scattering over centimetre-scale
lever arms.

## Physics populations

Standalone and combined ITS products are field-identical. A clean standalone
repeat is also field-identical, including tracks, flattened references, ROFs,
labels, and writer sidecars.

| ITS candidate | Tracks / hash | Matched / 142 | Efficiency | Fake / rate | Clone | Hit distribution |
|---|---:|---:|---:|---:|---:|---|
| Fixed slope pull | 212 / `46913a67a7e2fe7462e29df0db264fa8` | 135 | `0.950704` | 16 / `0.0754717` | 0 | 46x4, 16x5, 29x6, 121x7 |
| Measurement-only K | 174 / `5103bc6679a801339828eb7e684e4ba0` | 112 | `0.788732` | 13 / `0.0747126` | 0 | 135x4, 20x5, 10x6, 9x7 |
| K plus TrackletMinPt MS | 212 / `3b5b9df964b54a5d59616eff73bb6cb4` | 135 | `0.950704` | 16 / `0.0754717` | 0 | 44x4, 19x5, 28x6, 121x7 |

The corrected and fixed-pull products share 207 exact ordered reference
sequences. Five fixed-pull sequences disappear and five alternatives are
published; therefore the restored count and aggregate metrics are not a
hidden compatibility reconstruction.

| Change | ROF | Hits | Label / fake | Ordered references |
|---|---:|---:|---|---|
| removed | 0 | 6 | `0x3000000a8` / 0 | `1187,934,762,581,412,187` |
| added | 0 | 7 | `0x3000000a8` / 0 | `1461,1187,934,762,581,412,187` |
| removed | 0 | 4 | `0x30000010a` / 0 | `647,469,286,46` |
| added | 0 | 5 | `0x30000010a` / 0 | `1318,1031,821,647,469` |
| removed | 1 | 7 | `0x48000000b` / 0 | `1867,1781,1691,1639,1589,1558,1531` |
| added | 1 | 5 | `0x48000000b` / 0 | `1867,1781,1691,1639,1589` |
| removed | 1 | 4 | `0x500000060` / 0 | `1628,1582,1550,1521` |
| added | 1 | 5 | `0x500000060` / 0 | `1841,1754,1675,1628,1582` |
| removed | 2 | 4 | `0x80000006800000a3` / 1 | `2136,2062,2011,1951` |
| added | 3 | 4 | `0x8000000780000033` / 1 | `2623,2586,2562,2527` |

MFT outputs are unchanged field-by-field from the measurement-only candidate:

| MFT composition | Tracks / hash | Matched / 109 | Efficiency | Fake / rate | Clone / rate |
|---|---:|---:|---:|---:|---:|
| Standalone | 69 / `f6dee3f7a5f7def6b55900dbac734ef0` | 57 | `0.522936` | 2 / `0.0289855` | 0 / 0 |
| Combined | 100 / `98f9730d6fca9e738c3f20afc66d296d` | 78 | `0.715596` | 3 / `0.03` | 6 / `0.0550459` |

## Runtime characterization

This fixture is too short for a controlled benchmark, so wall time is not
used as an acceptance metric. The clean ITS runs consumed `5.89 s` and
`5.84 s` user CPU, compared with `5.87 s` for the measurement-only replay;
maximum resident set size stayed at approximately `2.65 GB`. The corrected
standalone MFT and combined runs consumed `4.46 s` and `8.37 s` user CPU.
No material memory increase or timing regression is resolved at this scale.
The hot-loop addition is one transition-angle load, one square outside the
shared statistic, one two-dimensional dot product, and one variance term; it
does not invoke propagation or fitting.

## Verdict and outlook

**Safe candidate for the completed longitudinal correction.** The requested
compatibility is recovered: ITS returns from 174 to 212 tracks with the former
aggregate efficiency/fake/clone values, while retaining a distinct,
deterministic candidate population; MFT publication is unchanged. The first
two planned steps are complete.

The next step, intentionally not started here, is a cheap analytic transverse
triplet primitive: signed XY curvature or a normalized transverse
collinearity residual, with the same `TrackletMinPt` contract and propagated
measurement/MS uncertainty. Only after its physics and timing characterization
should the broader stashed detector-neutral triplet fitting design be
introduced. Neighbour/road and hole harmonization remain out of scope.

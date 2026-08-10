# Transverse tracklet-direction compatibility validation

- Status: safe candidate
- Date: 2026-08-10
- Production commit: `1e2b552ac6`
- Tests and guards: `f47b7aef18`
- Design: [transverse tracklet-direction compatibility](../design/0017-transverse-tracklet-direction-compatibility.md)
- Pinned O2 package: `daily-20260717-0700-local1`
- Build: `/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`
- Fixture: `/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717`
- Fresh artifacts: `/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/transverse-tracklet-direction-20260810-Sb909v`

## Scope and acceptance

The new cell gate profiles the physical bend allowed by `TrackletMinPt`, then
tests only the excess using the analytic three-hit azimuth covariance and the
outgoing-transition multiple-scattering variance:

```text
kappaMax = |B2C Bz| / TrackletMinPt
deltaPhiMax = asin(kappa d01/2) + asin(kappa d12/2)
e = max(0, |deltaPhi| - deltaPhiMax)
chi2Phi = e^2 / (Var_measurement(deltaPhi) + sigmaTheta^2)
accept when chi2Phi < NSigmaCut^2.
```

The chord-domain clamp uses one admissible curvature for both chords. The
algorithm reads the two stored tracklet azimuths and does not run a triplet
fit. Disk and cylinder differ only in descriptor-selected construction of a
neutral global `(x,y)` observation and covariance.

## Test evidence

Focused tests cover:

- exact non-axis-aligned analytic derivatives and covariance contraction,
  including non-zero `Cov(x,y)` terms;
- the exact circular-chord `TrackletMinPt` envelope;
- a candidate rejected without MS and accepted when the outgoing angular
  variance is supplied;
- global-Z rotation invariance of points and covariances;
- disk global covariance and cylinder tangential-covariance projection;
- non-finite, coincident-point, and non-PSD failure without output mutation;
- one common `TrackerTraits` call path for both surface families; and
- source guards excluding family, detector, or source dispatch from the
  shared algorithm and cell loop.

The three focused test executables pass. The complete pinned serial campaign
also passes:

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 93
Total Test time (real) = 93.81 sec
```

`git diff --check` is clean and `git clang-format --diff` reports no changes.
Strict preflight passed with a valid AliEn token, no
`ALICEO2_CCDB_NOTOKENCHECK` bypass, the exact pinned package, and every
Geant4 dataset. Fixture checksums passed 43/43 before and 43/43 after. All
temporary environment-gated diagnostics were removed before the final clean
build and tests.

## Stage-level population

Tracklet formation is unchanged. The transverse decision first acts in cell
construction.

| ITS stage | Longitudinal correction | Plus transverse gate | Change |
|---|---:|---:|---:|
| Accepted tracklets | 1,886 | 1,886 | 0 |
| Accepted cells | 1,489 | 1,211 | -278 |
| Neighbour links | 1,050 | 876 | -174 |
| Accepted/published tracks | 212 | 184 | -28 |

| Standalone MFT stage | Longitudinal correction | Plus transverse gate | Change |
|---|---:|---:|---:|
| Accepted tracklets | 2,723 | 2,723 | 0 |
| Accepted cells | 604 | 559 | -45 |
| Neighbour links | 314 | 312 | -2 |
| Accepted/published tracks | 69 | 69 | 0 |

The MFT cell reduction has no publication effect in this fixture. The
standalone and combined MFT products remain exactly equal to their prior
candidate products.

## ITS population diagnosis

Identity is equal ROF plus equal ordered cluster-reference sequence. Against
the preceding 212-track candidate, the new 184-track product has:

```text
exactly retained = 176
removed          = 36
added            = 8
net               = -28
```

Every removed sequence has a valid MC association. All 36 associated MC
particles are below `TrackletMinPt=0.3 GeV/c`; none is at or above it:

```text
MC pT range = 0.103378668 .. 0.232226312 GeV/c
MC pT mean  = 0.162678509 GeV/c
removed fake/non-fake = 2 / 34
```

The eight additions are alternative low-pT candidates created by changed
competition (`0.129739359 .. 0.232226312 GeV/c`, one fake). This explains why
36 identity removals produce a net loss of 28. The fixture efficiency
denominator has no `pT >= 0.3 GeV/c` restriction, so its aggregate efficiency
decreases when the configured threshold is finally applied at tracklet
combination; this is not evidence of above-threshold inefficiency.

The complete removed/added inventories, including hit count, ROF, raw label,
fake flag, reconstructed and MC `pT`, and ordered references, are retained as
`its-population-removed.tsv` and `its-population-added.tsv` in the artifact
directory.

## Numerical representatives

The highest-MC-`pT` removed sequence is

```text
ROF 2, label 0x6800000a0, MC pT = 0.232226312 GeV/c
references = 2393,2260,2156,2075,2021,1964,1891
```

For its near-boundary cell `1964,2021,2075`, the decision is

```text
|deltaPhi|       = 0.0570027828 rad
deltaPhiMax      = 0.0421481082 rad
excess           = 0.0148546746 rad
sigmaTheta       = 0.00282521034 rad
Var(deltaPhi)    = 8.57417793e-6 rad^2
chi2Phi          = 25.7355703
NSigmaCut^2      = 25
```

It is rejected narrowly, not by a hard raw-angle cut. The MS term contributes
`7.9818e-6 rad^2`, about 93% of the total variance, and almost rescues a true
particle whose transverse momentum is still 23% below threshold. Lower-pT
examples are decisively separated; for the `0x3000000a8` sequence at
`MC pT=0.139254838 GeV/c`, representative cells have `chi2Phi` from 104 to
249 after the same MS treatment.

## Published physics products

ITS standalone and combined outputs agree. MFT retains its two
composition-specific prior products.

| Product | Prior tracks / hash | Current tracks / hash | Current matched / denominator | Efficiency | Fake / rate | Clone / rate |
|---|---|---|---:|---:|---:|---:|
| ITS standalone/combined | 212 / `3b5b9df964b54a5d59616eff73bb6cb4` | 184 / `932e38144a0575f495b29b532012f468` | 118 / 142 | `0.830986` | 15 / `0.0815217` | 1 / `0.00704225` |
| MFT standalone | 69 / `f6dee3f7a5f7def6b55900dbac734ef0` | unchanged | 57 / 109 | `0.522936` | 2 / `0.0289855` | 0 / 0 |
| MFT combined | 100 / `98f9730d6fca9e738c3f20afc66d296d` | unchanged | 78 / 109 | `0.715596` | 3 / `0.03` | 6 / `0.0550459` |

A clean ITS repeat is metric-file identical to the first standalone replay,
including the 184-track content hash. The combined ITS product has the same
count, hash, population metrics, and distributions as standalone. The MFT
hashes are exactly the preceding candidate hashes, so no MFT determinism
repeat was required.

## Runtime characterization

This fixture is too short for a controlled benchmark; the values below are
trend evidence only.

| Replay | Prior user CPU | Current user CPU | Current max RSS |
|---|---:|---:|---:|
| ITS standalone | 5.84-5.89 s | 6.59 s; repeat 5.95 s | 2.65 GB |
| MFT standalone | 4.46 s | 4.57 s | 1.29 GB |
| Combined | 8.37 s | 8.51 s | 5.41 GB |

No material runtime or memory regression is resolved. The hot-loop work is
three leaf projections, chord lengths, one profiled curvature envelope, and
one analytic covariance contraction; it adds no propagation or fit.

## Verdict and outlook

**Safe candidate.** The transverse gate does what the configuration says:
all historical ITS identity removals are associated with MC particles below
`TrackletMinPt`, MS is quantitatively active, no above-threshold associated
track is removed, and MFT publication is unchanged. The result is
deterministic and the complete pinned campaign passes.

The next optional step is the detector-neutral triplet fit, which can provide
an explicit momentum estimate and additional correlated selections. It
should be evaluated for incremental physics value and cost; it is not needed
to repair this analytic gate. Neighbour/road harmonization and hole logic
remain untouched.

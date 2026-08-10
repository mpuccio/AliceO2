# Adjacent-triplet Equation-19 truth calibration

- Date: 2026-08-10
- Verdict: **accept**
- Recommendation: **keep `maxChi2ClusterAttachment = 60` for this candidate**
- Comparison parent: `1df60514c1`
- Candidate production: `aaf140ed60`
- Candidate tests and guards: `37492f29d1`
- Candidate design and initial validation: `255ad1713e`
- Pinned package: `daily-20260717-0700-local1`
- Build: `/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`
- Fixture: `/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717`
- Calibration artifacts: `/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/adjacent-triplet-chi2-calibration-20260810-MJP7Z1`

## Question and conclusion

The inherited value 60 is not harsh for the adjacent-triplet Equation-19
quality. It is a very permissive provisional bound. At the current value,
every truth-associated neighbour candidate in this fixture passes: 718/718
for ITS and 277/277 for each MFT composition. All nine rejected ITS
candidates and all 19 rejected MFT candidates are unmatched at the four-hit
level.

The three standalone and six net combined MFT track losses are therefore not
direct over-rejection of their own adjacent triplets. Every cell link forming
those tracks has chi2 below 0.05. Their publication changes are indirect CA
level and road-competition effects after unmatched links are removed.

For four angular residuals minimized over one common curvature, the nominal
Equation-19 quality has three degrees of freedom. A chi2 of 60 has a nominal
upper-tail probability of about `5.88e-13`; it is not a statistically
calibrated three-degree-of-freedom cut. The observed truth distributions are
also substantially narrower than an ideal chi2(3), especially for MFT. This
fixture cannot distinguish how much comes from the conservative
`TrackletMinPt` MS hypothesis, candidate preselection, and the available hit
covariances. It therefore supports keeping 60 as a non-over-rejecting bound,
but not calling 60 an optimal or calibrated value.

A lower common threshold must be a separately reviewed, multi-fixture physics
change. Threshold 10 improves the MFT aggregate metrics here, but rejects
eight truth-associated ITS links and loses one matched ITS track. Threshold
20 retains every observed truth link but gives mixed publication effects.
No threshold was changed in production.

## Phase-A preservation and gates

The clean candidate was preserved before calibration in three commits:

| Commit | Scope |
|---|---|
| `aaf140ed60` | production global-measurement, stored-factor, and adjacent-triplet neighbour migration |
| `37492f29d1` | focused tests and source guards |
| `255ad1713e` | design and initial validation records |

`git diff --check 056f955cb6..255ad1713e` and pinned
`git clang-format --diff 056f955cb6 255ad1713e` passed. Strict preflight
resolved the requested package exactly, found a valid token with
`ALICEO2_CCDB_NOTOKENCHECK` unset, and validated every Geant4 dataset. The
post-commit serial gate passed 94/94 tests:

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 94
Total Test time (real) = 91.60 sec
```

No CUDA or HIP claim is made; no device compiler was exercised.

## Diagnostic method

Production code was not instrumented. An external overlay was built from the
already compiled `libO2ITSMFTTracking` objects with only `TrackerTraits.cxx`
recompiled through a forced-include symbol rename. The renamed calls capture
the four resolved `GlobalMeasurement`s and the public Equation-19 inputs and
result. The diagnostic copy of the library is loaded only in the tracker
process through `DYLD_LIBRARY_PATH`.

For every candidate it records:

- workflow leg, source ROF, four external cluster references, surfaces, and
  global positions;
- the unmodified chi2, curvature, curvature variance, and both MS angular
  variances;
- `Psi`, `rho`, the fitted residual `Psi + rho kappa`, marginal pulls, a
  Cholesky-whitened residual, and the complete 4x4 covariance;
- first-label four-hit truth association, MC label and pT, and the number of
  correctly labelled detector hits; and
- exact parent/current track containment plus same-label parent/current
  publication counts.

At threshold 60 the overlay logged exactly 898 ITS, 372 standalone-MFT, and
1,271 combined candidates. The standalone-MFT overlay output matches an
uninjected replay field by field, including 2,904 floating values with zero
maximum difference plus references, ROFs, seed patterns, and labels. The
diagnostic threshold scan changes only the returned comparison value outside
the production build; the logged chi2 is always the unscaled Equation-19
value.

## Candidate-level truth and pull distributions

Truth-associated means that the first correct source-0 label of all four hits
is the same. Unmatched includes mixed, fake, invalid, and unset four-hit
associations. This is deliberately stricter than assigning a label to a
published track.

| Leg | Class | Candidates | chi2 median / p95 / p99 / max | absolute whitened pull median / p95 / p99 / max |
|---|---|---:|---|---|
| ITS standalone | truth | 718 | 0.693 / 4.439 / 9.388 / 18.949 | 0.264 / 1.152 / 1.926 / 3.846 |
| ITS standalone | unmatched | 180 | 1.886 / 29.864 / 86.505 / 145.467 | 0.379 / 2.800 / 6.105 / 8.780 |
| MFT standalone | truth | 277 | 0.0103 / 0.0928 / 0.8799 / 4.980 | 0.0323 / 0.1606 / 0.3576 / 2.199 |
| MFT standalone | unmatched | 95 | 11.129 / 128.890 / 166.732 / 173.649 | 0.525 / 6.325 / 11.028 / 13.029 |

Combined ITS is numerically identical to standalone. Combined MFT has the
same 277 truth candidates and one additional unmatched accepted candidate;
its truth maximum is 4.98004 and unmatched maximum is 173.645. The full
quantiles are in `pull-summary.tsv`.

Numerical boundary representatives are:

| Leg/class | ROF | References | chi2 | residual `(theta0,phi0,theta1,phi1)` | whitened residual |
|---|---:|---|---:|---|---|
| ITS truth maximum | 1 | `1589,1639,1691,1781` | 18.9485 | `(0.01306,0.00664,0.00537,-0.00955)` | `(3.114,1.583,1.392,-2.192)` |
| ITS first rejection, unmatched | 0 | `441,449,783,789` | 62.5184 | `(-0.00091,-0.02449,0.00445,0.03856)` | `(-0.309,-5.217,1.076,5.835)` |
| MFT truth maximum | 0 | `59,87,118,158` | 4.98009 | `(0.00026,0.08109,0.01687,-0.05446)` | `(0.023,0.328,2.199,-0.186)` |
| MFT first rejection, unmatched | 0 | `152,193,212,237` | 60.0479 | `(0.04667,0.02534,-0.00096,0.01710)` | `(3.600,-1.721,6.631,0.404)` |

The residual components have different coordinate scales, so the whitened
values, not raw residual magnitudes, are the comparable pull diagnostic.

## Threshold scan

The scan used the common true-chi2 thresholds 10, 20, 40, 60, 80, 120, and
180. Each point was replayed twice for standalone ITS, standalone MFT, and
combined ITS+MFT. Fake-track counts are invariant across this scan: 15 for
ITS, two for standalone MFT, and three for combined MFT.

| Threshold | ITS tracks / matched | MFT standalone tracks / matched | MFT combined tracks / matched / clones | Accepted ITS truth / unmatched links | Accepted MFT truth / unmatched links |
|---:|---:|---:|---:|---:|---:|
| 10 | 183 / 117 | 69 / 57 | 98 / 77 / 5 | 710 / 153 | 277 / 45 |
| 20 | 184 / 118 | 67 / 55 | 96 / 75 / 5 | 718 / 168 | 277 / 59 |
| 40 | 184 / 118 | 66 / 54 | 94 / 75 / 3 | 718 / 171 | 277 / 72 |
| 60 | 184 / 118 | 66 / 54 | 94 / 75 / 3 | 718 / 171 | 277 / 76 |
| 80 | 184 / 118 | 65 / 53 | 93 / 74 / 3 | 718 / 175 | 277 / 83 |
| 120 | 184 / 118 | 65 / 53 | 93 / 74 / 3 | 718 / 179 | 277 / 88 |
| 180 | 184 / 118 | 65 / 53 | 93 / 74 / 3 | 718 / 180 | 277 / 95 |

The combined MFT unmatched-link counts are exactly one above standalone at
every threshold. ITS outputs are field-identical from 40 through 180. MFT
outputs at 40 and 60 are field-identical even though four additional
unmatched links pass at 60; outputs at 80, 120, and 180 are likewise
field-identical. Published populations therefore respond non-monotonically
to added links, as expected from CA level and road competition. Output count
is not a valid one-dimensional threshold objective.

Threshold 10 reproduces the parent standalone-MFT count and matched count,
but not its ordered content hash, and it directly rejects eight ITS truth
links. It is not evidence for restoring the old path by choosing 10.
Threshold 20 accepts all observed truth links and rejects more unmatched
links than 60, but it changes ITS content at fixed aggregate metrics, recovers
only one of three standalone-MFT truth references, and increases combined
MFT clones from three to five. This fixture does not establish a common
optimum below 60.

## MFT publication changes at 60

### Standalone

All links forming the three lost tracks pass by orders of magnitude. The
current 66-track population also replaces one eight-hit track by its exact
seven-hit same-label suffix.

| ROF | Label | MC pT | Correct hits | Parent references | Direct cell chi2 | Classification |
|---:|---:|---:|---:|---|---|---|
| 0 | 12884902077 | 0.05756 | 8 | `247,229,203,186,142` | 0.04948; 0.00347 | unique truth publication loss; direct links accepted |
| 0 | 12884902109 | 0.13789 | 6 | `565,544,511,479,423` | 0.01247; 0.00117 | unique truth publication loss; direct links accepted |
| 2 | 30064771095 | 0.86428 | 9 | `899,877,860,845,834` | 0.00655; 0.00223 | unique truth publication loss; direct links accepted |

The substitution is label 10737418505, seed `510 -> 254`, from
`545,512,480,424,416,383,358,319` to the exact suffix
`512,480,424,416,383,358,319`; its retained direct link has chi2 0.00727.

### Combined

The net 100 to 94 count is not six simple identity removals. Exact sequence
comparison finds seven parent removals and one same-label one-hit shortening,
for a net reduction of six. Excluding that substitution, the six removed
tracks are:

| ROF | Label | MC pT | Correct hits | Parent references | Direct cell chi2 | Classification |
|---:|---:|---:|---:|---|---|---|
| 0 | 17179869229 | 0.30529 | 8 | `551,525,490,460` | 0.02854 | clone removed; truth particle remains reconstructed |
| 0 | 12884902109 | 0.13789 | 6 | `565,544,511,479,423` | 0.01247; 0.00117 | unique truth publication loss |
| 0 | 12884902077 | 0.05756 | 8 | `247,229,203,186,142` | 0.04948; 0.00347 | unique truth publication loss |
| 2 | 30064771095 | 0.86428 | 9 | `899,877,860,845,834` | 0.00655; 0.00223 | clone removed; truth particle remains reconstructed |
| 2 | 23622320206 | 0.57607 | 8 | `765,754,746,732` | 0.03320 | clone removed; truth particle remains reconstructed |
| 4 | 40802189479 | 0.30466 | 8 | `1217,1207,1187,1180` | 0.00051 | unique truth publication loss |

Thus three combined removals are clone reduction and three remove the only
published track for a truth particle, exactly explaining matched 78 to 75
and clones six to three. The same label-10737418505 shortening as standalone
accounts for the seventh exact removal and sole addition. Every listed
direct link also passes threshold 10; none can be attributed to Equation-19
over-rejection.

Rejected MFT links are all unmatched, including mixed branches sharing hits
with the affected neighbourhoods. Examples include
`423,481,513,545` at chi2 123.09 and `423,481,513,544` at chi2 166.73.
Removing such links changes cell levels and subsequent road competition. It
does not make each resulting publication loss intrinsically desirable: the
clone removals are beneficial on this fixture, while the three unique-truth
losses are an indirect cost and require broader physics statistics.

## ITS comparison

At threshold 60, ITS rejects nine unmatched candidates and no truth
candidate. Relative to the parent 184-track product, 183 ordered sequences
are retained exactly and one fake track is replaced by another fake track
sharing two references. Standalone and combined ITS remain identical.

Threshold 10 is the only scanned value that rejects truth-associated ITS
links: 8/718, including the truth maximum at chi2 18.95, and it decreases the
published matched population from 118 to 117. This is the concrete reason
not to choose the MFT-favourable threshold 10 as a common value from this
fixture.

## Replay integrity and determinism

There are 42 successful scan replays: seven thresholds, three workflows, and
two executions per point. Every successful replay records 43/43 fixture
checks before and after. Candidate TSV files are byte-identical between
repeats. Persisted products match field by field for all 28 output legs
(standalone ITS, standalone MFT, combined ITS, and combined MFT at seven
thresholds), including references, ROFs, seed patterns, labels, and floating
track fields.

One initial ITS threshold-180 repeat failed strict CCDB retrieval of
`ITS/Calib/ClusterDictionary/1547590800000`. It is retained separately as
`its-t180-repeat`, has a recorded post-failure 43/43 fixture check, and is not
used. `its-t180-repeat-retry1` passed without a token-check bypass and is the
repeat used in all comparisons.

The parent products were not regenerated in a second worktree. They are the
fresh pinned products recorded by parent commit `1df60514c1` itself under
`transverse-tracklet-direction-20260810-Sb909v`, with hashes
`932e38144a0575f495b29b532012f468` (ITS),
`f6dee3f7a5f7def6b55900dbac734ef0` (standalone MFT), and
`98f9730d6fca9e738c3f20afc66d296d` (combined MFT). The current threshold-60
hashes are `e6da9d94faed581d5bff044993698e30`,
`32555b198d9b094f3f3600ec619cd2e2`, and
`96f4c632b7e0111501a63660774480ef`, respectively.

## Durable evidence

The artifact root contains:

- `threshold-scan.tsv`: link acceptance and published metrics for every scan
  point;
- `pull-summary.tsv`: truth/unmatched chi2 and whitened-pull quantiles;
- `mft-loss-classification.tsv`: complete standalone and combined loss and
  substitution classification;
- `its-t60/enriched.tsv`, `mft-t60/enriched.tsv`, and the two
  `combined-t60/*-enriched.tsv` files: per-candidate truth and final-track
  associations;
- every replay product, raw candidate table, standard metric JSON,
  provenance, before/after fixture checks, and repeat comparator log;
- `tools/`: the external overlay source and binary, replay/extraction scripts,
  and ROOT truth/population helpers; and
- `manifest.sha256`: SHA-256 manifest for the durable campaign files.

## Verdict

**Accept the adjacent-triplet candidate and keep 60 in this slice.** The
selection does not over-reject any truth-associated adjacent triplet on this
fixture; its direct rejections are desirable unmatched-tail suppression.
The MFT publication losses are a mixture of three combined clone removals and
indirect unique-truth losses caused by downstream competition, not evidence
that 60 is too small.

This is not an endorsement of 60 as a statistically meaningful Equation-19
cut. It is demonstrably loose. A separate review may test a lower common
threshold, with multiple fixtures spanning pT, eta, occupancy, and material,
and must calibrate the covariance/pulls rather than optimize one output
count. No energy loss, material iteration, detector-specific threshold,
longer-topology fitting, neighbour/road redesign, or hole-policy change is
started here.

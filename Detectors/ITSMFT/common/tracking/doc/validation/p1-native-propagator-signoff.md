# P1 native-Propagator physics sign-off characterization

Status: **INSUFFICIENT EVIDENCE** (2026-08-06)

This is a measurement record, not a physics acceptance. It characterizes the
accepted structural migration's native refit candidate against the last
replayable common-CA commit before native refit activation. No tracking
physics, cut, covariance policy, workflow default, or frozen ITS code changed
in this campaign.

## Provenance and reproducibility

| Item | Native candidate | Frozen reference |
|---|---|---|
| Commit | `38a798dfd9e890bf6683cd65e83563d1c157b9fa` | `5eb255fd62ad181d42e2a8aa3a51828c0b89a930` |
| Relationship | Current accepted M7f integration head | Direct parent of `0aaa1a9c8f`, which activates native refit |
| Worktree | `/Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement` | `/Users/mpuccio/alice/run3/O2-worktrees/m5c-operation-local-binding` |
| Durable build | `/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement` | `/Users/mpuccio/alice/run3/O2-worktree-builds/m5c-operation-local-binding` |
| Package, resolved package | `daily-20260717-0700-local1` | `daily-20260717-0700-local1` |
| Fixture | `/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717` | Same read-only fixture |
| Run protocol | `RUNNUMBER=303000`, `TIMESTAMP=1784207296000`, static diamond `(0,0,0)`, `pvRes=0.05`, Sync | Identical |
| MFT threading | `MFTCATrackerParam.nThreads=1` | Identical |

The fixture manifest was verified with `shasum -a 256 -c checksums.sha256`:
**43/43** before the six replays and **43/43** afterward. The native durable
build completed the required serial `ctest -L itsmft --output-on-failure -j1`:
**98/98 passed**, with no `Not Run` entries. Both builds then built
`O2lib-ITSMFTTracking`, `o2-its-ca-tracker-workflow`,
`o2-mft-ca-tracker-workflow`, and
`o2-itsmft-combined-ca-tracker-workflow`; the frozen build was already
up-to-date.

The replays use the checked-in common-CA scripts
`gate3-slice3-its-ca-validation/replay_tracking_its_common_ca.sh` and
`gate0-baseline/replay_tracking_common_ca.sh`. A representative native ITS
command is:

```sh
FIXTURE_DIR=/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717 \
REPLAY_DIR=/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/p1-native-propagator-signoff-20260806/native/its \
TIMESTAMP=1784207296000 RUNNUMBER=303000 \
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
./.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
bash Detectors/ITSMFT/common/tracking/test/gate3-slice3-its-ca-validation/replay_tracking_its_common_ca.sh
```

The frozen leg substitutes only its worktree/build and writes to
`.../frozen/its`; MFT additionally sets `MFT_CA_NTHREADS=1`. Combined legs
pipe the two cluster readers into the combined tracker with the same timestamp
and these exact key-values:

```text
ITSMFTCombinedCATrackerParam.enabled=true;
ITSCommonCATrackerParam.useDiamond=true;
ITSCommonCATrackerParam.diamondPos[0]=0;
ITSCommonCATrackerParam.diamondPos[1]=0;
ITSCommonCATrackerParam.diamondPos[2]=0;
ITSCommonCATrackerParam.pvRes=0.05;
MFTCATrackerParam.nThreads=1
```

## Persisted-output observations

The existing metric extractors were run on every standalone and combined
output. These are observations, not acceptance thresholds.

| Application | Frozen refit | Native candidate | Delta |
|---|---:|---:|---:|
| ITS tracks / hash | 203 / `ee7f7c794d60f2362fd2564258b7887e` | 212 / `46913a67a7e2fe7462e29df0db264fa8` | +9 |
| ITS reconstructable / matched / efficiency | 142 / 135 / 0.950704 | 142 / 135 / 0.950704 | 0 / 0 / 0 |
| ITS fake tracks / fake rate / clones | 14 / 0.0689655 / 0 | 16 / 0.0754717 / 0 | +2 / +0.0065062 / 0 |
| ITS mean chi2 | 6.80670 | 8.29868 | +1.49198 |
| MFT tracks / hash | 70 / `24737e73b7146bf3bd35a90a2517c527` | 68 / `8106b08571ca593c6b76ff72b761a680` | -2 |
| MFT reconstructable / matched / efficiency | 109 / 60 / 0.550459 | 109 / 58 / 0.532110 | 0 / -2 / -0.018349 |
| MFT fake tracks / fake rate / clones | 1 / 0.0142857 / 0 | 1 / 0.0147059 / 0 | 0 / +0.0004202 / 0 |
| MFT mean chi2 | 0.229033 | 0.203916 | -0.025117 |

Native combined output reproduces its standalone legs exactly: ITS has the
same 212/hash and MFT the same 68/hash. The existing field comparator reports
the ITS products match field-for-field; for MFT it reports 2,992 forward
float-projection values with maximum absolute and relative delta zero. Frozen
combined ITS also matches its standalone product field-for-field. Frozen MFT
has the same count/MC summary but a different exact-state content hash
(`1a858f59a722891799c1dacfdaafdf71`): the established comparator accepts it
as the historical forward float projection of the standalone product (3,080
values, maximum absolute `9.77768e-05`, maximum relative `5.73662e-08`). It
is not used as a byte-identity reference for the native combined leg.

Comparing frozen to native writer products fails first on the track-array
length, then on cluster-reference, ROF, and MC-label products for both ITS and
MFT. Therefore the observed difference is not the known undefined
`MFTTrack.mInvQPtSeed` byte artifact.

## Deterministic track association

`test/p1-native-propagator-signoff/compare_native_frozen_tracks.C` is the
analysis-only helper added for this campaign. It reads persisted track,
cluster-reference, ROF, MC-label, state, and covariance products. Pairs must
be in the same output ROF and share at least one external cluster reference.
It greedily chooses the maximum ordered tuple: same MC raw label, shared
reference count, Jaccard overlap, native index, frozen index. The last two
fields make ties deterministic.

| Association measure | ITS | MFT |
|---|---:|---:|
| Associated pairs | 202 | 68 |
| Native-only / frozen-only | 10 / 1 | 0 / 2 |
| Same MC raw label among pairs | 199 | 68 |
| Exact cluster-reference set | 141 | 68 |
| Mean shared references / mean Jaccard | 5.90099 / 0.920026 | 6.45588 / 1.0 |
| Native-only mean holes / frozen-only mean holes | 2.7 / 3.0 | n/a / 4.5 |

For associated tracks, the macro records absolute state and covariance-diagonal
deltas. The compact JSON artifacts retain all components in their persisted
parameter order. Notable maxima are ITS `|delta q2pt|=21.6943`,
`|delta chi2|=56.8955`, and final covariance diagonal `|delta|=54.6733`; MFT
`|delta invQPt|=24.0805`, `|delta chi2|=0.552517`, and final covariance
diagonal `|delta|=91814`. These are characterization values at each fitter's
written reference state, not pull distributions and not pass/fail limits.

Both applications have one surface family in this fixture (ITS barrel,
MFT disk). The persisted analysis can relate differences to output ROF,
cluster references, MC identity, hole count, state, and covariance. It cannot
recover per-step material crossed, pre-publication timestamps, or intermediate
CA/refit rejection populations from a written accepted-track file.

## Reused native-refit and covariance evidence

The current M7f focused tests were rerun successfully:

```text
testMFTNormalizedRefit.cxx
testITSNativeVsLegacyRefitCharacterizationHarness.cxx
testCovarianceSanitization.cxx
```

The rerun ITS harness writes eight scenarios to the artifact root. In all of
them production `fitTrackSeedLegs` exactly matches the native driver. The
single-hit, chi2-gate, and MinPt scenarios are rejected by both native and
frozen paths with respectively `LegAcceptanceFailure`,
`PredictedChi2Failure`, and `MinPtFailure`; the comparable multi-hit,
hole, material, shift-reference, and repeat-refit scenarios preserve pattern
and timestamp. This confirms the existing M5a/M5d local characterization is
still representative of the current native implementation.

M5d's earlier event-level diagnostic found `InvalidCovariance` eliminated
after sanitization, with remaining observed failures ordinary propagation,
predicted-chi2, whole-leg acceptance, or momentum-minimum outcomes. The
current sanitizer test confirms the implemented diagonal and pairwise
correlation policy at M7f. It does **not** prove that every covariance is
positive semidefinite, and this campaign did not re-instrument event-level
rejection categories. Consequently the earlier observation is supporting
context, not a fresh proof that every current residual loss is ordinary.

## Evidence versus inference

**Evidence**: both replays are reproducible on the same fixture and conditions;
the native candidate anchors and standalone/combined isolation hold; accepted
tracks have substantial reference-level continuity (especially MFT); and the
persisted differences include counts, references, ROFs, labels, states, and
covariances rather than only an undefined writer byte. The current unit/harness
tests find no reproduced invalid-covariance regression.

**Inference**: the MFT delta is consistent with two frozen-only selections
under a different final fit, while ITS changes mostly preserve the same
reconstructable MC denominator and matching count but alter selections and
fit states. This is consistent with the intended replacement of independent
fitters. It is not proof that the numerical/covariance/material behavior is
physics-acceptable.

## Verdict and required review questions

**INSUFFICIENT EVIDENCE.** No implementation defect is proven, so this task
does not propose a code correction. It also does not recommend acceptance:
there is no pre-agreed tolerance or physics-quality threshold, and the required
event-level CA cut-flow plus refit-rejection categories are not observable from
the current persisted products.

Physics review should decide whether to authorize a diagnostic-only next
campaign that exposes, without changing decisions, per-event counts for
tracklets/cells/neighbours/roads/refit attempts/accepted tracks and the
propagation/material/invalid-covariance/predicted-chi2/whole-leg/MinPt/adapter
rejection taxonomy. It should also define the physics quantities and
thresholds for matched residuals and for the ITS `+9` / MFT `-2` selection
deltas. Any material, covariance, propagation, MinPt, chi2, or hole-policy
change requires a separate reviewed decision after that evidence exists.

## Durable artifact index

All generated replay outputs, DPL logs/configuration snapshots, metric JSON,
association JSON, and the rerun harness JSON are under:

`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/p1-native-propagator-signoff-20260806`

Key machine-readable summaries are `native/its_metrics.json`,
`native/mft_metrics.json`, `frozen/its_metrics.json`, `frozen/mft_metrics.json`,
`its_track_association.json`, `mft_track_association.json`, and
`itsNativeVsLegacyRefitCharacterization.json`.

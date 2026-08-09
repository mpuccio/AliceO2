# Production pair-list and combined-plan validation

- Status: structurally validated review candidate; combined MFT physics sign-off pending
- Date: 2026-08-09
- Design: [flat pair-list graph and one combined plan](../design/0013-flat-pair-list-combined-plan.md)
- Package: `daily-20260717-0700-local1`
- Build: `O2-worktree-builds/itsmft-header-graph-cleanup`
- Fresh replay artifacts: `/private/tmp/o2-itsmft-flat-single-policy-20260809`

## Corrected contract

- Production graph authority is one ordered surface list plus one flat list of
  immediate forward pairs.
- The combined workflow has one graph, binding, workspace, `Tracker`,
  `TrackerTraits`, schedule, parameter record, and tracker run.
- The absent `(6,7)` pair is the ITS/MFT component boundary. There is no
  subgraph or source-plan API and no cross-component transition.
- The sole combined `TrackingParameters` record derives all scalar selection
  values from ITS. Its surface arrays concatenate seven ITS values and ten MFT
  values in global order; the start and seeding masks are `0x1ffff`.
- `SurfaceKind` selects cylinder/disk state and operation leaves only. It is not
  a tuning key. Source IDs remain loading/publication provenance only.
- Shared-workspace ROF timing resolves each surface's owner view and local
  layer rather than consulting a first/global source view.

Commit `a44649ac0c` restored the old 68-track combined MFT result with
family-keyed parameter records, but that architecture is superseded. The
production correction removes `parametersByKind`, `parametersForKind`, the
TimeFrame commit/storage plumbing, and every equivalent per-kind parameter
selection. The generic core again has one parameter record per iteration.

## Focused tests and guards

The focused suite verifies the following contracts:

- compile-time and source guards reject `parametersByKind`,
  `parametersForKind`, `mTrkParamsByKind`, and an
  `array<TrackingParameters, 2>` in the generic core;
- every combined scalar field equals the ITS configuration, while
  `AddTimeError`, geometry, material, resolution, systematic-error, and column
  extent arrays contain the ITS prefix followed by the MFT suffix;
- all 17 surfaces are eligible to seed, the graph has 15 transitions, and no
  transition crosses from cylinder rank 6 to disk rank 7;
- tracklet/cell operation partitions are derived from each descriptor's
  `SurfaceKind`, with guards against source ID, detector ID, or tuning-policy
  dispatch;
- a shared-workspace run with a 40-BC ITS ROF and an 80-BC MFT ROF reproduces
  each component's standalone timestamp interval. The longer MFT interval
  makes the test discriminating against the old first/global ITS clamp.

All focused tests pass. The complete pinned build covered the common tracking
library, ITS, MFT, and combined tracking/workflow libraries and executables,
the affected writers and all ITSMFT tests, plus Framework AnalysisSupport and
CCDBSupport. The complete serial suite passed with zero failures and zero
`Not Run`:

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 93
```

## Preflight and fixture integrity

Strict preflight resolved the exact pinned package, found a valid AliEn token,
confirmed `ALICEO2_CCDB_NOTOKENCHECK` was unset, and found every required Geant4
dataset. No authentication bypass was used.

The 43-file fixture manifest passed before and after all fresh replays:

```text
/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/
  pp-20ev-run303000-seed20260716-daily20260717/checksums.sha256
43/43 checksums passed before replay
43/43 checksums passed after replay
```

## Replay results

| output | tracks | content hash | comparison |
|---|---:|---|---|
| ITS standalone | 212 | `46913a67a7e2fe7462e29df0db264fa8` | frozen baseline |
| ITS combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` | field-for-field identical to standalone |
| MFT standalone | 68 | `8106b08571ca593c6b76ff72b761a680` | frozen baseline; 2,992 projected values, zero delta |
| MFT combined | 102 | `07e5addd74f374d2a4a70da69d0147aa` | unified-selection candidate |

The current 102-track MFT result is field-for-field identical to the earlier
deterministic 102-track result from `625778dcc3`: all 4,488 compared forward
values have zero delta, and references, ROFs, seed patterns, and labels match.
The intentionally undefined `MFTTrack.mInvQPtSeed` byte remains excluded from
binary comparison.

### Exact ordered-reference association

Using the ordered cluster-reference sequence as the track identity:

- 66 of the 68 standalone tracks occur verbatim in combined output;
- all 66 retain the same ROF, seed pattern, and raw MC label;
- 36 sequences occur only in combined output and two only in standalone;
- pairing two sequence changes as substitutions leaves 34 genuinely new
  combined tracks and no unpaired removal.

The two substitutions are exact:

1. ROF 0, label `0x280000898`: standalone
   `(410,386,351,323,292)` becomes the ten-hit extension
   `(566,543,514,476,427,410,386,351,323,292)`.
2. ROF 4, label `0x900000804`: standalone
   `(1154,1144,1123,1112,1095,1082)` becomes
   `(1144,1123,1112,1095,1082)`; reference 1154 occurs in no combined track.

The 34 genuinely new sequences consist of 32 four-hit tracks, one six-hit
track `(169,158,118,87,59,25)`, and one ten-hit track
`(771,762,751,741,722,717,707,699,691,683)`. The complete deterministic list,
including ROF, seed pattern, label, and chi2 for every raw addition and removal,
is recorded in `mft-association.log` beside the fresh replay artifacts.

Hit multiplicities are:

| hits | standalone | combined | raw additions | raw removals |
|---:|---:|---:|---:|---:|
| 4 | 0 | 32 | 32 | 0 |
| 5 | 27 | 27 | 1 | 1 |
| 6 | 14 | 14 | 1 | 1 |
| 7 | 13 | 13 | 0 | 0 |
| 8 | 1 | 1 | 0 | 0 |
| 9 | 9 | 9 | 0 | 0 |
| 10 | 4 | 6 | 2 | 0 |

Track counts by populated ROF are standalone `[36,4,12,3,13]` and combined
`[50,7,22,4,19]`. Raw additions distribute as `[15,3,10,1,7]`; the two raw
removals are in ROFs 0 and 4. All 2,304 persisted MFT ROF records retain the
same BC, orbit, RO-frame, and flags. Entry ranges all change because their
cumulative track population changes.

All 68 standalone and 102 combined tracks have valid labels. Standalone has
one fake label; combined has two. The 66 retained sequences preserve their raw
labels exactly and include the same one fake. All 36 raw additions are valid
and contain one fake; both raw removals are valid and non-fake. The two
substitutions preserve their respective labels. Aggregate reconstruction
metrics change from 58 to 80 matched tracks, one to two fakes, and zero to six
clones.

The writer does not persist `CommonTrack.timestamp`, so replay ROOT output
cannot associate a timestamp with each ordered MFT reference sequence. The
focused 40-BC/80-BC test supplies that evidence directly: both combined
components match their corresponding standalone timestamp interval, and the
MFT interval remains wider than ITS. This proves that the owner/local timing
lookup fixes the old cross-component clamp independently of the 102-track
population change.

### Localization to ITS scalar selection

The persisted output supports the following attribution. It does not expose
intermediate rejected candidates, so effects that share a stage are reported
as residual rather than assigned speculatively.

| field | ITS combined value | standalone MFT value | localized effect |
|---|---:|---:|---|
| `MinTrackLength` | 4 | 5 | Directly admits all 32 four-hit additions once their roads exist. |
| `TrackletMinAbsX` | 0 (disabled) | 0.05 | One genuinely new ten-hit track has fitted `abs(x)=0.00396024`; this is consistent with the relaxed cut, but persisted fitted X cannot prove the earlier candidate passed because of it. |
| `PVres` | 0.05 | 0.01 | Changes projection/search windows; individual final tracks cannot be uniquely assigned to this change. |
| `ColBins` | 256 | 64 | Changes disk LUT quantization and candidate formation; individual output differences cannot be isolated from persisted data. |
| `CorrType` | `NONE` | `LUT` | The current common path does not provide a persisted per-track discriminator for this setting; no independent population delta is claimed. |
| `StartLayerMask` | `0x1ffff` | local `0x3ff` | Both make every MFT surface start-eligible. With `(6,7)` absent, this does not explain cross-component traversal. |
| `IndexRowMin/Max` | 0/0 | -20/20 | Disk binding maps 0/0 to the same -20/20 fallback, so there is no effective difference. |
| `MinPt` entries | four zeroes | six zeroes | No threshold is applied for the MFT hit-count slots in this replay. |

The remaining operation thresholds are equal in this campaign and therefore
do not directly explain the population delta: `TrackletMinPt=0.3`,
`RowBins=128`, `NSigmaCut=5`, `CellDeltaTanLambdaSigma=0.007`,
`CellRoadRCut=0.05`, `MaxChi2ClusterAttachment=60`, `MaxChi2NDF=30`, sharing
disabled with `SharedMaxClusters=0`, and reseeding at six clusters. All 66
retained tracks have changed fitted parameter tuples, which localizes a broader
effective selection/refit change but cannot separate LUT, projection-window,
and candidate-competition effects from persisted output alone.

## Static checks and limitation

`git diff --check` and pinned-environment `git clang-format --diff` are clean.
Source guards reject family-keyed tracking parameters and source/detector-based
operation dispatch while retaining descriptor-selected cylinder/disk leaves.

`CMAKE_CUDA_COMPILER` is `NOTFOUND` in this build. Host and device-guard
compilation pass, but this validation makes no real CUDA or HIP compiler claim.

The 102-track combined MFT output is recorded only as a candidate for later
unified physics sign-off. This correction does not declare it an accepted
baseline and does not redesign detector tuning, hole policy, or graph topology.

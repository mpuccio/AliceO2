# Coordinate-cut retirement validation

- Status: structurally validated review candidate; unified MFT physics sign-off pending
- Date: 2026-08-09
- Design: [coordinate-cut retirement](../design/0014-coordinate-cut-retirement.md)
- Base: `90773323ba`
- Production commit: `bea08be655`
- Test commit: `7287cfeb56`
- Package: `daily-20260717-0700-local1`
- Build: `O2-worktree-builds/triplet-tracking-rnd-scratch`
- Fresh artifacts: `/private/tmp/o2-itsmft-coordinate-cut-retirement-20260809-pI5Oqm`
- Historical characterization: `/private/tmp/o2-itsmft-flat-single-policy-20260809`

## Validated contract

`TrackletMinAbsX` was an orientation-dependent MFT compatibility selection,
not a forward-coordinate validity guard. It is removed without replacement.
`CellRoadRCut` was a legacy coarse-road search optimization whose common-port
three-cell inequalities did not preserve the legacy call site or mathematics.
The transplanted selection and its cylinder no-op are removed without a hidden
constant or renamed disk parameter.

Both names are absent from common production configuration, projection,
device parameters, Tracker/TrackerTraits orchestration, cell construction,
and refit. The source guard scans common `include` and `src`; the only common
source occurrence of either name is the negative test's forbidden-token list.
No detector/source dispatch was introduced. Cylinder and disk enter the same
four public traversal stages, with descriptor-selected coordinate leaves.

Forward seed construction still rejects unordered or coincident Z surfaces
and zero transverse separation. Disk tracklet projection still protects its
layer/vertex-Z denominators. Focused leaf tests exercise these boundaries.
The flat pair list, hole policy, one combined scalar record, missing `(6,7)`
pair, owner/local ROF timing, and current publication ownership are unchanged.
Legacy production MFT `MFTTrackingParam::ROADclsRCut` is outside this slice and
was not modified.

The later-deletion inventory remains the one in the design note:
`MFTFwdTrackHelpers.{h,cxx}` can disappear after normalized identity checking
and direct native-refit invocation move to their generic/application
boundaries; `DetectorRefitSupport.h` may then be narrowed. No legacy workflow
removal is included here.

## Build and tests

The fresh build is `RelWithDebInfo` with `BUILD_TESTING=ON`. Configuration
needed the established exclusion
`O2_ROOT_MACRO_EXCLUSION_LIST=Detectors/ITSMFT/common/tracking/test/p1-native-propagator-signoff/compare_native_frozen_tracks.C`
because that standalone validation macro is intentionally not a CTest macro.
All 89 ITS/MFT-labelled test targets, the ITS, MFT, and combined replay
executables, and `FrameworkAnalysisSupport`/`FrameworkCCDBSupport` were built.

The eight directly affected tests passed first. The complete required serial
suite then passed with no failure or `Not Run`:

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 93
Total Test time (real) = 131.09 sec
```

`git diff --check` and pinned-environment `git clang-format --diff` are clean.
`CMAKE_CUDA_COMPILER` is `NOTFOUND`; host and device-guard compilation pass,
but this campaign makes no real CUDA/HIP compiler claim.

## Preflight and fixture integrity

Strict preflight resolved the exact pinned package, found a valid AliEn token,
confirmed `ALICEO2_CCDB_NOTOKENCHECK` was unset, and found every required
Geant4 dataset. Replays used CCDB condition cap `1784207296000` and no token
bypass.

The durable fixture's 43-entry manifest passed before replay and again after
all primary and repeat replays:

```text
/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/
  pp-20ev-run303000-seed20260716-daily20260717/checksums.sha256
43/43 checksums passed before replay
43/43 checksums passed after all replays
```

The artifact directory retains two failed initial ITS attempts for
reproducibility: `its-standalone` exposed the not-yet-built DPL plugin
libraries, and `its-standalone-rerun` was denied IPC socket creation by the
sandbox. After building both plugin libraries and granting normal DPL
IPC/network access, `its-standalone-success`, `mft-standalone`, and `combined`
completed. Neither failed attempt modified the fixture or produced a usable
physics output.

## Replay summary

| output | historical tracks/hash | current tracks/hash | matched | fake | clone |
|---|---|---|---:|---:|---:|
| ITS standalone | 212 / `46913a67a7e2fe7462e29df0db264fa8` | 212 / `46913a67a7e2fe7462e29df0db264fa8` | 135 | 16 | 0 |
| ITS combined | 212 / `46913a67a7e2fe7462e29df0db264fa8` | 212 / `46913a67a7e2fe7462e29df0db264fa8` | 135 | 16 | 0 |
| MFT standalone | 68 / `8106b08571ca593c6b76ff72b761a680` | 69 / `f6dee3f7a5f7def6b55900dbac734ef0` | 57 | 2 | 0 |
| MFT combined | 102 / `07e5addd74f374d2a4a70da69d0147aa` | 100 / `4dacc3058740cce83656a3bcb71def95` | 78 | 3 | 6 |

The ITS standalone product matches its historical product field-for-field,
the combined ITS product matches its historical product field-for-field, and
current standalone and combined ITS match each other field-for-field. ITS
efficiency remains `135/142 = 0.950704`, fake rate `16/212 = 0.0754717`, and
clone rate zero.

Independent repeat MFT standalone and combined replays match their respective
first fresh products field-for-field, including ordered references, ROFs, seed
patterns, labels, and initialized track fields. The deliberately undefined
`MFTTrack.mInvQPtSeed` byte remains excluded by the established comparator.

MFT standalone efficiency changes from `58/109 = 0.532110` to
`57/109 = 0.522936`; fake rate changes from `1/68 = 0.0147059` to
`2/69 = 0.0289855`; clones remain zero. Combined efficiency changes from
`80/109 = 0.733945` to `78/109 = 0.715596`; fake rate changes from
`2/102 = 0.0196078` to `3/100 = 0.0300000`; clones remain `6/109 = 0.0550459`.
These are characterization results, not preservation failures.

## Exact MFT association and attribution

Ordered external cluster references define track identity. All retained
tracks preserve ROF, seed pattern, and raw MC label. All 2,304 persisted ROF
records preserve BC, orbit, RO-frame, and flags; cumulative entry ranges
change when track populations change. Every current label is valid:
standalone has 69/69 valid with two fake labels, and combined has 100/100 valid
with three fake labels.

### Standalone: 68 to 69

There are 65 exact retained sequences, four additions, and three removals.
All 65 retained sequences are also field-for-field identical. The changes are:

| result | ROF | seed | label | ordered references | interpretation |
|---|---:|---:|---|---|---|
| add 10-hit | 0 | `0x3ff` | `0x280000898` | `566,543,514,476,427,410,386,351,323,292` | replaces the removed five-hit suffix `410,386,351,323,292`; removal of the compatibility X selection allows the longer candidate already seen when that selection was disabled |
| add 10-hit | 2 | `0x3ff` | `0x6800000b1` | `771,762,751,741,722,717,707,699,691,683` | direct X-cut case: fitted `x=-0.003969064 cm`, below the retired `0.05 cm` threshold |
| add 6-hit | 2 | `0x3f0` | `0x680000385` | `874,875,843,840,821,815` | road-candidate extension; the same extension appears in the isolated combined road-cut delta |
| add 9-hit | 4 | `0x3fe` | fake `0x8000000980000059` | `1220,1211,1192,1182,1158,1150,1133,1121,1103` | road-candidate substitution for the removed valid sequence using reference `1212`; the same substitution appears in combined |
| remove 7-hit | 0 | `0x1fc` | `0x400000136` | `550,524,489,458,438,395,366` | road-candidate competition; identically removed in combined |

The removed nine-hit substitute is
`1220,1212,1192,1182,1158,1150,1133,1121,1103`, seed `0x3fe`, valid label
`0x980000059`. The five-hit removal is the suffix paired with the first row.
The full records, including chi2, are in `mft-standalone-change.log`.

### Combined: 102 to 100

The historical combined configuration already disabled `TrackletMinAbsX`
through the one ITS-derived scalar record. Therefore this comparison isolates
`CellRoadRCut` removal: 97 sequences are retained field-for-field, three are
added, and five are removed as the enlarged cell-candidate set changes road
extension and competition.

The three additions are:

1. ROF 0, seed `0xfc`, label `0x4000010b4`, six-hit extension
   `180,175,173,162,119,88` of removed four-hit
   `173,162,119,88`.
2. ROF 2, seed `0x3f0`, label `0x680000385`, six-hit extension
   `874,875,843,840,821,815` of removed four-hit
   `874,875,843,840`.
3. ROF 4, seed `0x3fe`, fake label `0x8000000980000059`, nine-hit
   `1220,1211,1192,1182,1158,1150,1133,1121,1103`, substituting for the
   removed valid-label sequence with reference `1212` in place of `1211`.

The other removals are ROF-0 seven-hit
`550,524,489,458,438,395,366` and ROF-2 four-hit
`842,839,820,815`. Thus deleting the old pre-cut does not monotonically add
published tracks: it admits cells before road construction, after which the
existing candidate extension, conflict resolution, and publication ordering
can replace or suppress shorter competitors. No physical selection was
silently restored to force the historical count.

### Current standalone versus combined

This is structural evidence only. Of 69 standalone sequences, 68 occur in
combined output with identical ROF, seed pattern, and raw label. Combined has
32 additional sequences (29 four-hit, one five-hit, and two six-hit), while
standalone alone has the six-hit sequence
`1154,1144,1123,1112,1095,1082`; combined publishes its five-hit suffix
`1144,1123,1112,1095,1082` with the same label. Population equality is not an
acceptance condition.

## Timestamp evidence and review boundary

The ROOT writer does not persist `CommonTrack.timestamp`, so ordered-reference
tracks cannot be assigned a timestamp from replay output. The full suite's
`CombinedComponentsUseOwnROFTimingInOneCombinedPass` test supplies direct
nonzero-track evidence: a 40-BC ITS source and an 80-BC MFT source each match
its standalone timestamp interval, and MFT remains wider. This slice did not
alter that owner/local timing implementation.

The MFT changes are deterministic and localized to removal of the obsolete X
compatibility selection and the transplanted road optimization. They remain a
candidate for later unified physics sign-off. This slice does not redesign
holes, the flat pair graph, combined scalar policy, or legacy workflows.

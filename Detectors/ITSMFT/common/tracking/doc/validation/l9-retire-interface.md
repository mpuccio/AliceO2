# L9 TrackingInterface retirement validation

Date: 2026-08-07
Branch: `codex/itsmft-postm7-l9-retire-interface`
Base: `3d15feac7a` (integrated L8)
Production commit: `4d6b2114ee`
Tests/guard commit: `ce7664691d`
Package: `daily-20260717-0700-local1`
Source: `/Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement`
Build: `/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement`
Replay descendant: `f36a6e5b4c` (includes the subsequent L10 commits)

## Scope and result

L9 deletes the standalone `TrackingInterface` and makes ITS and MFT workflows
compose the generic components directly:

```text
workflow-owned timing/publication context
  -> Loader::load(TimeFrame&, inputs)
  -> Tracker::run(TimeFrame&, TrackerTraits&)
  -> workflow publication adapter
```

The workflows retain ownership of `TimeFrame` lifetime, raw ROFs, timing-table
backing storage, publication validity, typed sidecars, typed output conversion,
and DPL lifecycle. No replacement facade or lifecycle-owner class was added.
The combined workflow remains the reference for the same direct composition.

The initial replay attempt was blocked by unavailable live CCDB authentication.
A retry at the current descendant completed successfully and provides replay
evidence for that descendant. Because the branch now includes L10 commits, the
retry is not presented as an isolated L9-only replay.

## Provenance and user work

The durable build cache contains:

```text
CMAKE_HOME_DIRECTORY:INTERNAL=/Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement
```

The pinned user stash remains untouched:

```text
stash@{0}: On codex/itsmft-integration: user WIP: TripletTrackingRAndD.md
```

The fixture is:

```text
/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717
```

with `RUNNUMBER=303000`, `TIMESTAMP=1784207296000`, static ITS diamond
position `(0,0,0)`, `pvRes=0.05`, synchronous tracking, and one MFT CA thread.

Durable artifacts are in:

```text
/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l9-retire-interface-20260807/
```

The pre-replay and post-replay checksum commands were:

```zsh
cd /Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717
shasum -a 256 -c checksums.sha256
```

Both durable logs report **43/43 `OK`**, with zero checksum failures.

## Production and ownership checks

The deleted production/test surface is:

```text
Detectors/ITSMFT/common/tracking/include/ITSMFTTracking/TrackingInterface.h
Detectors/ITSMFT/common/tracking/src/TrackingInterface.cxx
testTrackingInterfaceLoadFailureContract.cxx
testTrackingModeOwnership.cxx
testM6e1StandaloneMFTMigration.cxx
testITSCommonCANThreads.cxx
```

The standalone workflow sources now construct and retain their own `TimeFrame`,
`TrackerTraits`, `Tracker`, and operation adapter, and invoke the generic
loader/tracker APIs directly. The generic library has no interface-only
failure sentinel (`kDroppedTimeFrameResult`/`isDroppedTimeFrame`). The focused
`testL9StandaloneCompositionGuard` checks the deleted names, direct calls,
frame ownership, reset path, and absence of a replacement lifecycle facade.

## Build and test validation

The affected common library, standalone workflows, combined workflow, writers,
and focused tests were built in the durable build. The final serial command
was:

```zsh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
ctest --test-dir /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
  -L itsmft --output-on-failure -j1
```

The L9 head result was **97/97 passed, 0 failed, 0 `Not Run`**. The current
descendant replay/build result is **96/96 passed, 0 failed, 0 `Not Run`**; the
one removed registration belongs to the subsequent L10 cleanup. The focused
standalone-composition guard and existing failure, reset, source-isolation,
publication, and combined-composition tests all passed.

The final source checks are:

```zsh
git diff --check 3d15feac7a HEAD
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
git clang-format --diff 3d15feac7a HEAD
```

## Replay protocol and evidence

The native standalone scripts were run with the pinned fixture and durable
build:

```zsh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
FIXTURE_DIR=/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717 \
REPLAY_DIR=/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l9-retire-interface-20260807/its-standalone \
TIMESTAMP=1784207296000 RUNNUMBER=303000 DIAMOND_POS_X=0 DIAMOND_POS_Y=0 DIAMOND_POS_Z=0 PV_RES=0.05 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
bash Detectors/ITSMFT/common/tracking/test/gate3-slice3-its-ca-validation/replay_tracking_its_common_ca.sh

O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
FIXTURE_DIR=/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717 \
REPLAY_DIR=/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l9-retire-interface-20260807/mft-standalone \
TIMESTAMP=1784207296000 RUNNUMBER=303000 MFT_CA_NTHREADS=1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
bash Detectors/ITSMFT/common/tracking/test/gate0-baseline/replay_tracking_common_ca.sh
```

The combined invocation used the two cluster readers piped into
`o2-itsmft-combined-ca-tracker-workflow` with:

```text
ITSMFTCombinedCATrackerParam.enabled=true
ITSCommonCATrackerParam.useDiamond=true
ITSCommonCATrackerParam.diamondPos[0..2]=0
ITSCommonCATrackerParam.pvRes=0.05
MFTCATrackerParam.nThreads=1
```

The first-attempt logs reached the DPL CCDB backend and then failed before
conditions were delivered to the tracker:

```text
TGrid::Connect returned nullptr. May be due to missing alien token
internal-dpl-ccdb-backend: exit 128
```

Those first-attempt files are retained under
`l9-retire-interface-20260807/`; they are not used as replay evidence.

### Successful retry

The successful retry artifacts are under:

```text
/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l9-retire-interface-20260807-rerun/
```

The retry used the exact standalone commands above with each `REPLAY_DIR`
changed to the corresponding `l9-retire-interface-20260807-rerun/` directory,
and the same two-reader combined pipeline and configuration values recorded
above. It ran against current HEAD `f36a6e5b4c`.

The retry verified:

| Leg | Tracks | Content hash | Standalone/combined | L8-parent initialized content |
| --- | ---: | --- | --- | --- |
| ITS standalone | 212 | `46913a67a7e2fe7462e29df0db264fa8` | — | field match |
| ITS combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` | field match | field match |
| MFT standalone | 68 | `8106b08571ca593c6b76ff72b761a680` | — | field match; 2992 projected floats, max abs/rel 0 |
| MFT combined | 68 | `8106b08571ca593c6b76ff72b761a680` | field match | field match; 2992 projected floats, max abs/rel 0 |

The standard metric extractors reported the same input counts and MC/fake/
clone summaries for standalone and combined legs. The established
field-level comparison returned zero for ITS standalone versus combined, MFT
standalone versus combined, and each current leg versus the retained L8
parent. ROOT file bytes are not used because ROOT metadata is run-specific;
the known undefined `MFTTrack.mInvQPtSeed` byte artifact remains excluded.

The retry checksum log reports **43/43 `OK`** after replay, matching the
43/43 pre-replay result.

For context only, the successful L8 parent artifacts recorded the accepted
native values ITS `212` / `46913a67a7e2fe7462e29df0db264fa8` and MFT `68` /
`8106b08571ca593c6b76ff72b761a680`, with standalone/combined parity. Those
artifacts are not substituted for a current L9 replay.

No real GPU/device build was run: `nvcc` and `hipcc` are unavailable in the
pinned environment.

## Exact L10 boundary

L10 is still bounded to the existing operation seam. It may move
`completeAccepted()` and `resetAdapterState()` into application publication
handling and retain only the call-scoped refit operation if more than one real
implementation remains; otherwise that seam should be deleted. It must not
introduce callbacks, a central detector dispatch, or another lifecycle owner.
Any remaining public compatibility header or adapter operation requires a
named live owner and caller before deletion. L10 must repeat the full replay
gate and native-refit ordering/holes/rejection checks before changing this
boundary.

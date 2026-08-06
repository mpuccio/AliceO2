# L2 dead scratch and fixed-source forwarders validation

- Branch: `codex/itsmft-postm7-l2-dead-forwarders`
- Base: `6bb9471bcf055948c0c9e9519dba0dfc6d521348`
- Production commit: `c3c9b4f3eb`
- Tests/guard commit: `2b3cd68f05`
- Package: `daily-20260717-0700-local1`
- Worktree: `/Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement`
- Build: `/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement`

## Scope and deletion evidence

This behavior-preserving slice removed only:

- `SurfaceTrackingScratch::mPValphaX` and its reset, allocator, swap, and
  staged-allocator bookkeeping;
- the forwarding `SurfaceTrackingScratch::resetScratch()` spelling, retaining
  `reset()` as the sole scratch operation; and
- `MultiSourceTimeFrameLoader::loadITSAndMFT()` and
  `resetITSAndMFTEvent()`.

The real-scratch loader tests now construct ordered
`LoadTargetImplSurface`/`AtomicLoadBinding` values and call `loadEvent()`. They
retain coverage for multi-source atomicity, non-dense-source rejection,
source-qualified measurements, independent scratch backfills, retry after a
malformed replacement, and reset of event-derived state while configured
capacity remains owned by the caller.

The focused guard is the `L2DeadScratchAndFixedSourceForwarderGuard` case in
`testMultiSourceTimeFrameLoader.cxx`. It recursively scans common tracking
`include`, `src`, and `test` files, including CMake files, for all four deleted
names. The literals are split in the guard source so the guard cannot exempt
itself. A repository search after the test migration returned no occurrence in
common tracking include/src/test.

No graph, `TimeFrame`, loader ownership, workflow, output, sidecar, Propagator,
or physics behavior was changed. L3 remains the next authorized slice.

## Build and test validation

The durable build was reconfigured in place and the affected common library,
common tests, and ITS/MFT/combined workflow libraries and executables were
built with:

```sh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
sh -c 'build=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement; targets=$(ninja -C "$build" -t targets all | awk "/^O2test-itsmft-tracking-[^:]+: phony$/{sub(/:$/,\"\",\$1); print \$1}"); ninja -C "$build" -j12 O2lib-ITSMFTTracking $targets'
```

The workflow consumers were then rebuilt explicitly:

```sh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
ninja -C /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement -j8 \
  O2lib-ITSMFTCombinedCAWorkflow O2exe-itsmft-combined-ca-tracker-workflow \
  O2lib-ITSCAWorkflow O2exe-its-ca-tracker-workflow \
  O2lib-MFTTracking O2exe-mft-ca-tracker-workflow
```

The first post-edit test run used stale test objects and produced ABI-shaped
bus/segmentation failures in tests that instantiate `SurfaceTrackingScratch`.
The library object was newer than those dependent objects. Cleaning the
common tracking test target graphs and rebuilding their dependencies in the
same durable build corrected the issue; the affected tests then passed. A
direct CTest invocation without the O2 environment also failed in ROOT setup
(`O2_ROOT`/`ROOT_DYN_PATH` were unset). The final run used the repository O2
environment:

```sh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
ctest --test-dir /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
  -L itsmft --output-on-failure -j1
```

Result: **100% tests passed, 0 tests failed out of 98**. Every registered test
executed; there were no `Not Run` entries. The focused loader and normalized
source tests also passed independently.

Formatting checks passed:

```sh
git diff --check 6bb9471bcf055948c0c9e9519dba0dfc6d521348 HEAD
git clang-format --diff 6bb9471bcf055948c0c9e9519dba0dfc6d521348 HEAD
```

## Replay gate

Fixture: `/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717`.

```sh
cd /Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717
shasum -a 256 -c checksums.sha256
```

All 43 entries passed before replay and all 43 passed after replay. Fresh
outputs and logs are under:

`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l2-dead-forwarders-20260806/`

The accepted replay scripts were run with run `303000`, condition-not-after
`1784207296000`, static ITS diamond position `(0,0,0)`, `pvRes=0.05`, and one
MFT thread. The successful replay products were:

```sh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
FIXTURE_DIR=/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717 \
REPLAY_DIR=/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l2-dead-forwarders-20260806/its-standalone-env \
TIMESTAMP=1784207296000 RUNNUMBER=303000 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
Detectors/ITSMFT/common/tracking/test/gate3-slice3-its-ca-validation/replay_tracking_its_common_ca.sh

O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
FIXTURE_DIR=/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717 \
REPLAY_DIR=/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l2-dead-forwarders-20260806/mft-standalone-env \
TIMESTAMP=1784207296000 RUNNUMBER=303000 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- bash \
Detectors/ITSMFT/common/tracking/test/gate0-baseline/replay_tracking_common_ca.sh
```

| Leg | Tracks | Content hash |
| --- | ---: | --- |
| ITS standalone | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT standalone | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT combined | 68 | `8106b08571ca593c6b76ff72b761a680` |

The field-level comparator reported initialized writer, CommonTrack,
cluster-reference, ROF, label, and sidecar content equal to the accepted L1
outputs for ITS and MFT. It also reported standalone/combined field equality
for both legs. The MFT comparison covered 2,992 forward-state/covariance
values with maximum absolute and relative projection deltas of zero. Only the
known undefined `MFTTrack.mInvQPtSeed` byte artifact remains excluded by the
established comparison protocol.

No CUDA or HIP compiler was available in the pinned environment (`nvcc` and
`hipcc` were absent), so no device build result is claimed.

The user stash `stash@{0}` (`user WIP: TripletTrackingRAndD.md`) remained
untouched and was neither restored nor committed.

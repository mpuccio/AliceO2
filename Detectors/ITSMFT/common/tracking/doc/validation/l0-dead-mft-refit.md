# L0 dead common typed-MFT refit/export validation

Status: **implemented and fully validated; 98/98 ITS/MFT tests passed**

L0 deletes the unused common typed-MFT refit/export path and migrates
`testMFTNormalizedRefit.cxx` to the live generic native-refit result path. No
tracking algorithm, workflow default, publication adapter, writer, or physics
code was changed.

## Provenance

| Item | Value |
| --- | --- |
| Branch | `codex/itsmft-postm7-l0-dead-mft-refit` |
| Base | `8f033406dd075628c2d2a4e381b6a7243ef9a552` |
| Package | `daily-20260717-0700-local1` |
| Worktree | `/Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement` |
| Durable build | `/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement` |
| Fixture | `/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717` |
| Run protocol | `RUNNUMBER=303000`, `TIMESTAMP=1784207296000`, static diamond `(0,0,0)`, `pvRes=0.05`, synchronous, MFT threads `1` |

The build cache was confirmed with:

```text
CMAKE_HOME_DIRECTORY:INTERNAL=/Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement
```

The named user stash remained untouched:
`stash@{0}: On codex/itsmft-integration: user WIP: TripletTrackingRAndD.md`.

## Deletion proof

The following common-tracking files and their CMake source entry were deleted:

- `include/ITSMFTTracking/MFTAdapterRefit.h`;
- `include/ITSMFTTracking/MFTCATrack.h`;
- `src/MFTAdapterRefit.cxx`.

Repository searches over common tracking and the MFT/ITS application paths found
no production or test reference to the deleted identifiers. The live
`MFTCATrackWriterSpec` publication adapter and MFT workflow/configuration code
were retained; they are unrelated writer/publication code, not the deleted
common typed refit path. `testM7eAdapterBoundaryGuard` checks exact identifier
boundaries so the live `MFTCATrackerParam` name is not treated as a deleted
`MFTCATrack` type.

## Build and tests

The existing durable build was reconfigured and the affected targets built:

```sh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
./.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
ninja -C /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement -j4 \
  O2lib-ITSMFTTracking O2lib-ITSMFTTrackingParams \
  O2test-itsmft-tracking-mft-normalized-refit \
  O2test-itsmft-tracking-m7e-adapter-boundary-guard \
  O2test-itsmft-tracking-m7f-runtime-core-guard \
  O2exe-its-ca-tracker-workflow O2exe-mft-ca-tracker-workflow \
  O2exe-itsmft-combined-ca-tracker-workflow
```

Focused tests passed: the migrated generic MFT refit test and both source
guards passed 3/3.

An initial serial run executed all 98 registered ITS/MFT tests with 95 passes
and the following three reported failures:

| Test | Observed failure |
| --- | --- |
| `testTimeFrameDetectorLayouts.cxx` | segmentation violation at the unchanged `TrackerTraits::initialiseTimeFrame` line 409 |
| `testCATrackerFailureContract.cxx` | subprocess abort: `CA tracker exceeded memory limit: New set maximum below current used (newMax: 531, used: 532)` |
| `testM7dNontemplatedTrackerGuard.cxx` | segmentation violation at the unchanged `TrackerTraits::initialiseTimeFrame` line 409 |

Those results were not a valid parent/feature classification because the exact
test binaries had not been rebuilt after the feature build state changed. The
exact feature targets were then rebuilt with:

```sh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
run-in-o2-env.zsh -- ninja -C /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement -j4 \
  O2test-itsmft-tracking-timeframe-detector-layouts \
  O2test-itsmft-tracking-catracker-failure-contract \
  O2test-itsmft-tracking-m7d-nontemplated-tracker-guard
```

The same target command was run in the independently configured parent build,
replacing both build paths with
`/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement-l0-parent`.
Each rebuilt executable was then run directly with detailed output:

```sh
O2_BUILD_DIR=<feature-build> O2_PACKAGE=daily-20260717-0700-local1 \
run-in-o2-env.zsh -- <absolute-test-executable> \
  --log_level=all --report_level=detailed --catch_system_errors=yes
```

All three feature executables returned exit status 0. The detached parent at
`10a0dcd1c9^ = 8f033406dd075628c2d2a4e381b6a7243ef9a552` was independently
configured and built with the same package, build type, macro-audit cache
exclusion, executable arguments, and environment. All three parent
executables also returned exit status 0.

| Test | Initial full-suite symptom | Rebuilt feature | Rebuilt parent | Classification |
| --- | --- | ---: | ---: | --- |
| `testTimeFrameDetectorLayouts` | segfault at `TrackerTraits::initialiseTimeFrame:409` | 0 | 0 | stale-build/environment issue |
| `testCATrackerFailureContract` | memory-limit abort (`newMax: 531, used: 532`) | 0 | 0 | stale-build/environment issue |
| `testM7dNontemplatedTrackerGuard` | segfault at `TrackerTraits::initialiseTimeFrame:409` | 0 | 0 | stale-build/environment issue |

The identical targeted CTest filter passed 3/3 in both builds after the exact
rebuild. Detailed feature and parent logs are retained at:

`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l0-dead-mft-parent-validation-20260806/`

The parent evidence therefore rules out an L0 regression and no production
correction is warranted. After the exact feature binaries were rebuilt, the
complete required command was rerun:

```sh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
run-in-o2-env.zsh -- ctest \
  --test-dir /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
  -L itsmft --output-on-failure -j1
```

Final result: **100% tests passed, 0 tests failed out of 98**. Every registered
test executed; there were no `Not Run` entries.

The exact parent/reference direct-run logs and the feature logs used for the
failure classification remain at:

`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l0-dead-mft-parent-validation-20260806/`

No L0 production source changed in `TrackerTraits`, `Tracker`, or CA runtime
code. The three reported failures were stale-build/environment issues, not L0
regressions.

The durable build was configured with this cache-only exclusion because the
repository already contains an unregistered P1 diagnostic macro that makes the
root-macro audit fail during configuration:

```text
-DO2_ROOT_MACRO_EXCLUSION_LIST:STRING=Detectors/ITSMFT/common/tracking/test/p1-native-propagator-signoff/compare_native_frozen_tracks.C
```

No source CMake exclusion was added, and all 98 registered `itsmft` tests still
executed.

## Replay evidence

Fixture checksums passed all 43 entries before replay and all 43 entries after
replay with:

```sh
cd /Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717
shasum -a 256 -c checksums.sha256
```

Replay commands reused the checked-in scripts:

```sh
FIXTURE_DIR=<fixture> REPLAY_DIR=<artifact>/its-standalone \
TIMESTAMP=1784207296000 RUNNUMBER=303000 O2_BUILD_DIR=<build> \
O2_PACKAGE=daily-20260717-0700-local1 \
run-in-o2-env.zsh -- bash \
Detectors/ITSMFT/common/tracking/test/gate3-slice3-its-ca-validation/replay_tracking_its_common_ca.sh

FIXTURE_DIR=<fixture> REPLAY_DIR=<artifact>/mft-standalone \
TIMESTAMP=1784207296000 RUNNUMBER=303000 MFT_CA_NTHREADS=1 O2_BUILD_DIR=<build> \
O2_PACKAGE=daily-20260717-0700-local1 \
run-in-o2-env.zsh -- bash \
Detectors/ITSMFT/common/tracking/test/gate0-baseline/replay_tracking_common_ca.sh
```

The combined replay used both cluster readers and
`o2-itsmft-combined-ca-tracker-workflow` with the fixed diamond, `pvRes=0.05`,
and `MFTCATrackerParam.nThreads=1` configuration. All replay logs and ROOT
outputs are retained outside Git at:

`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l0-dead-mft-refit-20260806/`

| Leg | Tracks | Content hash |
| --- | ---: | --- |
| ITS standalone | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT standalone | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT combined | 68 | `8106b08571ca593c6b76ff72b761a680` |

The existing field-level comparator reported standalone/combined equality for
both legs, including zero MFT forward-projection delta. It also reported
field-level equality between each L0 standalone output and the durable M7f
parent output. ROOT file bytes themselves were not compared because ROOT file
metadata is non-semantic; the known undefined `MFTTrack.mInvQPtSeed` artifact
was not observed as an initialized-content difference.

The replay evidence supports unchanged MFT publication behavior. The targeted
failure investigation is resolved as a stale-build/environment issue, and the
complete serial L0 test gate is now green.

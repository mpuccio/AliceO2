# L0 dead common typed-MFT refit/export validation

Status: **implemented; validation blocked by three pre-existing integration-test failures**

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

The required serial suite executed all 98 registered ITS/MFT tests. 95 passed;
the following three failed on the integration base's existing runtime tests,
outside the L0 diff:

| Test | Observed failure |
| --- | --- |
| `testTimeFrameDetectorLayouts.cxx` | segmentation violation at the unchanged `TrackerTraits::initialiseTimeFrame` line 409 |
| `testCATrackerFailureContract.cxx` | subprocess abort: `CA tracker exceeded memory limit: New set maximum below current used (newMax: 531, used: 532)` |
| `testM7dNontemplatedTrackerGuard.cxx` | segmentation violation at the unchanged `TrackerTraits::initialiseTimeFrame` line 409 |

No L0 production source changed in `TrackerTraits`, `Tracker`, or CA runtime
code. These failures prevent claiming a green full-suite gate and require a
separate baseline/runtime investigation; they were not corrected in L0.

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

The replay evidence supports unchanged MFT publication behavior, but the full
L0 validation remains blocked until the three unrelated baseline tests are
resolved or explicitly waived by review.

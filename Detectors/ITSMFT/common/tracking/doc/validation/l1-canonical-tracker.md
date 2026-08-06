# L1 canonical `Tracker` validation

Status: **complete; behavior-preserving canonical-file migration validated**

L1 makes `Tracker.h` and `Tracker.cxx` the canonical declaration and
definition of the non-templated common tracker. `CATracker.{h,cxx}` and the
forwarding relationship are deleted. No tracking algorithm, workflow default,
publication adapter, writer, or physics behavior was changed.

## Provenance

| Item | Value |
| --- | --- |
| Branch | `codex/itsmft-postm7-l1-canonical-tracker` |
| L1 base | `08583c2047f3280004e4c0699b7e4f8426ef3a8d` |
| Production commit | `6dcd72c4cf` |
| Test/guard commit | `bdb976ed1c` |
| Package | `daily-20260717-0700-local1` |
| Worktree | `/Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement` |
| Durable build | `/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement` |
| Fixture | `/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717` |
| Protocol | run `303000`, timestamp `1784207296000`, static diamond `(0,0,0)`, `pvRes=0.05`, synchronous, MFT threads `1` |

The build cache resolves to the L1 worktree. The named user stash remains
untouched:
`stash@{0}: On codex/itsmft-integration: user WIP: TripletTrackingRAndD.md`.

## Ownership and deletion proof

- `include/ITSMFTTracking/Tracker.h` now contains the `Tracker` declaration,
  result types, and outcome contract; it is no longer a forwarding header.
- `src/Tracker.cxx` contains the existing orchestration definition, moved
  without algorithmic changes.
- CMake, ITS/MFT workflow includes, and common tests use `Tracker.h`.
- The failure-contract fixture is named `testTrackerFailureContract.cxx`.
- The M7f guard checks that the old header/source paths, old includes, the
  standalone `CATracker` identifier, and old CMake entries are absent from
  common production include/src and CMake paths.

Outside the focused guard's intentional negative-check literals, the remaining
`CATracker` spellings are application workflow/configuration vocabulary such as
`CATrackerSpec` and `CATrackerParam`; no common-tracking `CATracker` file,
include, class, alias, declaration, or CMake entry remains.

## Build and tests

The existing build was reconfigured in place, retaining the cache-only root
macro exclusion already required by the unregistered P1 diagnostic macro:

```sh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
cmake -S /Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement \
  -B /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
  -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_TESTING=ON \
  -DCMAKE_INSTALL_PREFIX=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement/stage \
  -DO2_ROOT_MACRO_EXCLUSION_LIST=Detectors/ITSMFT/common/tracking/test/p1-native-propagator-signoff/compare_native_frozen_tracks.C
```

Affected common libraries, ITS/MFT/combined workflow targets, and all 111
registered common tracking test targets rebuilt successfully. The focused
canonical-file guard, M7e adapter guard, and renamed Tracker failure-contract
test all passed.

The required complete suite was run with:

```sh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
run-in-o2-env.zsh -- ctest \
  --test-dir /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
  -L itsmft --output-on-failure -j1
```

Result: **100% tests passed, 0 tests failed out of 98**. Every registered test
executed and there were no `Not Run` entries.

## Replay gate

The fixture checksum command passed all 43 entries before replay and all 43
entries after replay:

```sh
cd /Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717
shasum -a 256 -c checksums.sha256
```

Standalone replays used the checked-in scripts and the pinned package. Fresh
outputs, logs, and metrics are retained under:

`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l1-canonical-tracker-20260806/`

The combined replay used both cluster readers and the explicit combined opt-in:

```sh
o2-its-cluster-reader-workflow --with-mc --input-dir <fixture> --its-cluster-infile o2clus_its.root \
  | o2-mft-cluster-reader-workflow --with-mc --input-dir <fixture> --mft-cluster-infile mftclusters.root \
  | o2-itsmft-combined-ca-tracker-workflow -b --run --condition-not-after 1784207296000 \
      --configKeyValues 'ITSMFTCombinedCATrackerParam.enabled=true;ITSCommonCATrackerParam.useDiamond=true;ITSCommonCATrackerParam.diamondPos[0]=0;ITSCommonCATrackerParam.diamondPos[1]=0;ITSCommonCATrackerParam.diamondPos[2]=0;ITSCommonCATrackerParam.pvRes=0.05;MFTCATrackerParam.nThreads=1'
```

The successful combined output is in `combined-final/`. A preceding attempt
that omitted the required combined opt-in stopped in workflow preflight before
any tracking device ran; its log is retained in `combined-rerun/` and is not
part of the acceptance result.

| Leg | Tracks | Content hash |
| --- | ---: | --- |
| ITS standalone | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT standalone | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT combined | 68 | `8106b08571ca593c6b76ff72b761a680` |

The existing field-level comparator reported initialized writer, CommonTrack,
cluster-reference, ROF, label, and sidecar content matching the L0 output for
both standalone legs. It also reported standalone/combined equality for both
legs. The MFT comparison covered 2,992 forward-state/covariance values with
zero projection delta. Only the previously approved undefined
`MFTTrack.mInvQPtSeed` byte artifact remains excluded by the comparison
protocol; no initialized-content difference was observed.

No CUDA or HIP device build was claimed; the pinned environment has no usable
`nvcc` or `hipcc` toolchain.

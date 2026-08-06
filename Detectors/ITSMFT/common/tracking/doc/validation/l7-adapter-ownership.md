# L7 adapter ownership validation

Status: complete, behavior-preserving structural migration.

L7 moves detector timing and publication compatibility ownership to the ITS,
MFT, and combined workflow edges. The generic `TimeFrame`, loader, tracker,
and tracker-traits headers retain only borrowed runtime views and generic
results; they do not own detector ROF tables, clocks, validity, typed
sidecars, or typed output state.

## Provenance

| Item | Value |
| --- | --- |
| Branch | `codex/itsmft-postm7-l7-adapter-ownership` |
| L6 base | `ed79914bfe` |
| Production commit | `90d12e7f74` |
| Test/guard commit | `6f0e4fe37a` |
| Package | `daily-20260717-0700-local1` |
| Source | `/Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement` |
| Durable build | `/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement` |
| CMake source cache | `/Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement` |
| Fixture | `pp-20ev-run303000-seed20260716-daily20260717` |
| Run / condition timestamp | `303000` / `1784207296000` |
| Diamond / PV resolution | `(0,0,0)` / `0.05` |
| MFT threads | `1` |
| User stash | `stash@{0}` remains the untouched `TripletTrackingRAndD.md` WIP stash |

Durable replay products and metric summaries are under:

`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l7-adapter-ownership-20260806/`

## Ownership result

Standalone ITS and MFT workflow contexts now own their fixed ROF tables,
timing construction, publication clock views, and detector compatibility
sidecars. The combined DPL task owns the corresponding combined contexts.
`TrackingInterface` and `SurfacePlanTrackingParticipant` borrow runtime ROF
views and publication adapters; they do not allocate, reset, invalidate, or
decide publication validity for detector state. The borrowed views are valid
only for the event invocation and are invalidated by the owning workflow on
load failure, drop, structural failure, exception, and replacement.

The generic boundary guard is `testL7AdapterOwnership`. It checks the public
generic headers for fixed ROF-table, clock/validity, workflow-mask, and
sidecar ownership tokens, and checks that the workflow contexts contain the
required adapter-owned state. It also exercises timing/mask view semantics
and detector-local sidecar clearing. The remaining participant/interface
forwarders are intentionally narrow borrowed-context bridges and are the
L8/L9 cleanup boundary, not independent owners.

## Build and tests

Affected targets were built with the existing durable build:

```zsh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
ninja -C /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement -j2 \
  libO2ITSMFTTracking.dylib \
  Detectors/ITSMFT/common/tracking/test/all \
  o2-its-ca-tracker-workflow o2-mft-ca-tracker-workflow \
  o2-itsmft-combined-ca-tracker-workflow
```

The required complete suite was run serially:

```zsh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
ctest --test-dir /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
  -L itsmft --output-on-failure -j1
```

Result: **102/102 registered tests passed, 0 failed, and no `Not Run`
entries**. The focused L7 ownership and standalone/combined workflow tests
also passed.

## Replay protocol and results

The fixture was checked with the repository manifest before the replay and
again after all replay and extraction work:

```zsh
cd /Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717
shasum -a 256 -c checksums.sha256
```

Both checks produced **43/43 `OK`** entries. The standalone commands were:

```zsh
FIXTURE_DIR=/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717 \
REPLAY_DIR=/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l7-adapter-ownership-20260806/its-standalone \
TIMESTAMP=1784207296000 RUNNUMBER=303000 \
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- bash \
Detectors/ITSMFT/common/tracking/test/gate3-slice3-its-ca-validation/replay_tracking_its_common_ca.sh

FIXTURE_DIR=/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717 \
REPLAY_DIR=/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l7-adapter-ownership-20260806/mft-standalone \
TIMESTAMP=1784207296000 RUNNUMBER=303000 MFT_CA_NTHREADS=1 \
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- bash \
Detectors/ITSMFT/common/tracking/test/gate0-baseline/replay_tracking_common_ca.sh
```

The combined replay used the same fixture, timestamp, diamond, PV, and MFT
thread configuration. Its workflow log and products are in the `combined`
artifact directory; the required digitizer configuration was copied there
from the fixture before launching the pipeline.

Metric extraction produced:

| Leg | Input clusters | Input ROFs | Tracks | Candidate hash |
| --- | ---: | ---: | ---: | --- |
| ITS standalone | 4103 | 2304 | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS combined | 4103 | 2304 | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT standalone | 1249 | 2304 | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT combined | 1249 | 2304 | 68 | `8106b08571ca593c6b76ff72b761a680` |

The established `compare_common_ca_output.C` was run for current standalone
versus combined ITS and MFT, and current L7 versus the L6 standalone rerun
for both detectors. All four returned zero and reported field-by-field
equality. The MFT comparisons covered 2,992 forward float values with maximum
absolute and relative difference zero. This checks initialized writer fields,
CommonTrack content, references, ROFs, labels, and sidecars through their
persisted representations. ROOT file bytes are not compared because ROOT
metadata is run-specific; the known undefined `MFTTrack.mInvQPtSeed` byte
artifact remains the only excluded field.

No CUDA or HIP validation is claimed: neither `nvcc` nor `hipcc` is available
in the pinned environment.

## Evidence and boundary

Measured facts are the test result, checksum result, candidate metrics, and
field-level comparisons above. The ownership conclusion is an architectural
inference from the source guard and the focused lifecycle tests: detector
backing storage and publication compatibility now live at workflow/adapter
edges, while the generic core consumes only event-lifetime views and generic
results. No tracking equation, refit policy, covariance policy, ordering,
hole behavior, output format, or workflow default was changed.

The next authorized work is L8/L9 only: retire the remaining borrowed
participant/interface composition and direct workflow composition of Loader
plus Tracker. That work may delete the temporary borrowed-context forwarders,
but must not move raw ROFs or publication state into `TimeFrame` and must
preserve the L7 ownership and replay gates. L8/L9 are not started by this
commit.

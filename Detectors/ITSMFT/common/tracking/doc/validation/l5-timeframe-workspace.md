# L5 TimeFrame workspace and event-reset validation

Status: historical validation record; L5 complete and L6 subsequently completed.

L5 makes `TimeFrame` the sole owner of generic mutable event state and CA
workspace. The frame remains passive: graph construction and configuration are
performed by `Tracker`, while Loader and adapters invoke narrow frame storage
operations. Raw ROFs, timing backing storage, publication validity, typed
sidecars, writers, and workflow lifecycle remain outside `TimeFrame`.

## Provenance

| Item | Value |
| --- | --- |
| Branch | `codex/itsmft-postm7-l5-timeframe-workspace` |
| L4 base | `fdb1d5eeb5` |
| Production commit | `a673237e62` |
| Test/guard commit | `12a06424d0` |
| Runtime ROF repair | `ff7dbbef3e` |
| Lifecycle test migration | `5bcceaf349` |
| Package | `daily-20260717-0700-local1` |
| Build | `/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement` |
| Source | `/Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement` |
| CMake source cache | The build cache points to the source path above. |
| Fixture | `pp-20ev-run303000-seed20260716-daily20260717` |
| Run / timestamp | `303000` / `1784207296000` |
| Diamond / PV resolution | `0,0,0` / `0.05` |
| User stash | `stash@{0}` remains `user WIP: TripletTrackingRAndD.md`; it was not restored or modified. |

The validation artifact directory is:

`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l5-timeframe-workspace-20260806/`

## Ownership and reset evidence

During atomic configuration, `TimeFrame` creates one private workspace entry
for each source-qualified binding. `Tracker`, `TrackingInterface`,
`SurfacePlanTrackingParticipant`, and the direct loader borrow the corresponding
workspace. No participant or interface owns a second generic workspace,
configuration, or reset authority. `SurfaceTrackingScratch` remains an
implementation type for kernels and loader staging, but its live configured
instances are frame-owned. L6 removed the loader-target indirection.

`TimeFrame::resetEvent()` is the single generic event reset. It clears
normalized measurements, runtime ROF views, CA working storage, CommonTracks,
references, labels, and event-derived support state. It preserves graphs,
per-iteration parameters, source-qualified partitions, allocator identity, and
reserved capacities. Failed configuration and failed staging are non-mutating;
successful replacement resets the prior event once before installing the new
validated event view.

The focused `l5-timeframe-workspace` test covers allocator and capacity reuse,
configuration preservation, event-data clearing, reset counting, and failed
configuration preservation. Existing failure-contract and composition tests
cover successful replacement, recoverable drops, structural failures, and
exception cleanup. The source guard checks that no participant/interface owns a
generic scratch member or independent reset path and that the frame owns the
workspace collection.

## Build and test commands

The affected targets were built in the pinned ALICE environment with the
existing durable build:

```zsh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
ninja -C /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement -j2 \
  libO2ITSMFTTracking.dylib \
  Detectors/ITSMFT/common/tracking/test/all \
  Detectors/ITSMFT/ITS/tracking/test/all \
  o2-itsmft-combined-ca-tracker-workflow
```

The authoritative full suite was run serially in the same environment:

```zsh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
ctest --test-dir /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
  -L itsmft --output-on-failure -j1
```

Result: **101/101 tests passed, 0 failed, and no `Not Run` entries**. Focused
workspace, failure-contract, MFT, ITS, and combined-composition tests were
also rebuilt and passed. A preliminary direct-shell CTest invocation lacked
the pinned O2 runtime environment and was not used as the validation result;
the environment-wrapped run above is authoritative.

## Fixture and replay evidence

Before and after the replay campaign, the fixture was checked with:

```zsh
cd /Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717
shasum -a 256 -c checksums.sha256
```

All **43/43** entries were `OK` before and after.

The standalone ITS and MFT replay scripts and the combined workflow were run
with the package, run, timestamp, diamond, PV-resolution, and one-thread
settings in the provenance table. Their logs, configs, ROOT outputs, and
metrics are in the artifact directory above. The resulting candidate counts
and hashes were:

| Leg | Tracks | Candidate hash |
| --- | ---: | --- |
| ITS standalone | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT standalone | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT combined | 68 | `8106b08571ca593c6b76ff72b761a680` |

The existing `gate4-b42b-validation/compare_common_ca_output.C` comparator
reported matching initialized output content for standalone ITS and MFT
against the L4 parent, and for combined legs against both their standalone
products and the L4 combined parent. The MFT comparison covered 2992 forward
float values with maximum absolute and relative differences of zero. This
comparison includes CommonTracks, references, ROFs, labels, and initialized
writer fields; it excludes only the established undefined
`MFTTrack.mInvQPtSeed` byte artifact.

No device validation was claimed: neither `nvcc` nor `hipcc` was available in
the pinned environment.

## Exact L6 boundary

L6 may make Loader act directly on `TimeFrame` and remove the current private
load-target hierarchy once an equivalent staging/commit contract is proven.
It must preserve atomic failed staging, successful replacement, and the current
event/configuration separation. L6 must not move raw ROFs, timing construction
or backing storage, publication validity, typed sidecars, writers, or workflow
event ownership into `TimeFrame`; those remain later workflow/adapter concerns.

# L8 Tracker orchestration validation

- Branch: `codex/itsmft-postm7-l8-retire-engine`
- L8 base: `a59e038a76` (`codex/itsmft-integration` at branch creation)
- Production migration: `1d65906aed`
- Test/guard migration: `3a58a2fb82`
- Replay-parity correction: `15109d23d0`
- Package: `daily-20260717-0700-local1`
- Source: `/Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement`
- Durable build: `/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement`
- Build cache: `CMAKE_HOME_DIRECTORY` points to the source path above.
- User work: `stash@{0}` remains `user WIP: TripletTrackingRAndD.md`; it was not
  restored, modified, dropped, or committed.

## Result

L8 is **complete**. `Tracker` is the only generic orchestration component.
`TimeFrame` remains passive and owns the configured graphs, partitions,
workspace, event state, and generic results. `TrackerTraits` remains the
borrowed CPU/GPU kernel seam. The combined task now composes workflow-owned
timing/publication context, direct atomic loading, direct `Tracker::run()` for
the explicit ITS/MFT order, and workflow publication adapters.

The standalone `TrackingInterface` remains intentionally as the narrow L9
bridge. L8 does not retire it and does not move raw ROFs, timing backing
storage, publication validity, sidecars, writers, or typed output into the
generic core.

## Structural change

The production deletion removed:

- `TrackingEngine.{h,cxx}`;
- `TrackingParticipant.h`;
- `ParticipantId.h`; and
- `SurfacePlanTrackingParticipant.{h,cxx}`.

`Tracker::run(TimeFrame&, TrackerTraits&)` now performs the existing sequence:
per-iteration binding/graph initialization, tracklets, cells, neighbours,
roads/refit, generic result commit, and failure reset/classification. The
tracker stores only its non-owning operation adapter and source identifier;
it stores no frame, graph, workspace, event, timing, sidecar, or output state.

The guard `testL8TrackerOrchestrationGuard` proves that the deleted engine,
participant, and schedule vocabulary is absent from common tracking and the
combined workflow, that `Tracker` has no frame/event ownership members, and
that the combined source directly contains Tracker construction and `run()`
calls. Existing direct-composition and failure-contract tests retain schedule,
source-isolation, CPU-traits, and reset coverage.

## Build and test commands

The durable build was reconfigured in place only as needed. Affected targets
were built with:

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

The final gate was:

```zsh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
ctest --test-dir /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
-L itsmft --output-on-failure -j1
```

Result: **100/100 passed, 0 failed, 0 Not Run**. The suite included the
common library, ITS, MFT, combined workflow, source guards, reset/failure
contracts, and compiled validation macros.

## Replay provenance and evidence

Fixture:
`/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717`.
The pinned replay used `RUNNUMBER=303000`,
`TIMESTAMP=1784207296000`, static ITS diamond position `(0,0,0)`,
`pvRes=0.05`, synchronous tracking, and one MFT common-CA thread.

The fixture manifest was checked before replay with:

```zsh
cd /Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717
shasum -a 256 -c checksums.sha256
```

The recorded final check is
`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l8-retire-engine-20260806/checksums-after.log`:
**43/43 `OK`**, with no failures. The initial pre-replay check also returned
43/43 `OK`.

Durable replay artifacts are under:

`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l8-retire-engine-20260806/`

with `its-standalone-rerun/`, `mft-standalone/`, and `combined/` containing
the logs, DPL configurations, ROOT products, and extracted metrics. The first
ITS replay attempt produced 103 tracks because the refactor had omitted the
old frame-to-traits magnetic-field handoff. That discrepancy was localized to
the L8 change, corrected by `traits.setBz(frame.getBz())` in `Tracker::run()`,
and the rebuilt replay then reproduced the accepted result. The failed first
attempt is retained as `its-standalone/` for auditability; acceptance uses the
corrected `its-standalone-rerun/` artifact.

| Leg | Input clusters | Input ROFs | Tracks | Candidate hash |
| --- | ---: | ---: | ---: | --- |
| ITS standalone, corrected | 4103 | 2304 | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS combined | 4103 | 2304 | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT standalone | 1249 | 2304 | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT combined | 1249 | 2304 | 68 | `8106b08571ca593c6b76ff72b761a680` |

The established `gate4-b42b-validation/compare_common_ca_output.C` returned
zero for corrected standalone versus combined ITS and MFT. It also returned
zero for corrected L8 standalone versus the L7 parent artifacts. The MFT
comparison covered 2,992 forward-projected float values with maximum absolute
and relative differences of zero. These comparisons found equality for the
initialized writer fields, CommonTracks, references, ROFs, labels, and
sidecars. ROOT file bytes were not used because ROOT metadata is run-specific;
the known undefined `MFTTrack.mInvQPtSeed` byte artifact remains excluded.

## Device validation and next boundary

No GPU validation is claimed. `nvcc` and `hipcc` are unavailable in the pinned
environment, so no real device build was possible.

L9 is now unambiguous and bounded: remove the standalone `TrackingInterface`
and its float sentinel by composing the standalone workflow edges directly as
workflow-owned timing/publication context -> `Loader::load(TimeFrame&, ...)`
-> `Tracker::run(TimeFrame&, TrackerTraits&)` -> publication adapter. L9 must
not create a facade/coordinator and must preserve raw-ROF ownership, timing
backing storage, sidecars, writers, and the current physics/output behavior.

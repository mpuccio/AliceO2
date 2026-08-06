# L4 TimeFrame configuration validation

Status: implementation complete; L5 is the next authorized slice.

L4 makes `TimeFrame` the sole owner of static tracking configuration. The
workflow/application supplies declarations; `Tracker::initialize()` builds
and validates graphs locally and commits the complete replacement through
`TimeFrame::commitConfiguration()`. The frame remains passive and
`TimeFrame::wipe()` clears event data without discarding configured graphs,
parameters, source-qualified bindings, allocator identity, or capacity
requirements.

## Evidence

The focused test is
`O2test-itsmft-tracking-l4-timeframe-configuration` and covers:

- failed configuration preserving a prior valid configuration and allocator;
- source-qualified bindings with non-identity ordered positions;
- distinct per-iteration seeding/hole declarations;
- successful replacement of the complete configuration; and
- event reset preserving static configuration and capacity.

The ownership scan checks the canonical Tracker, TrackingInterface,
SurfacePlanTrackingParticipant, and their implementations for independent
graph, parameter, and pool members. `TrackerTraits` retains only its borrowed
kernel allocator handle pending L5 workspace ownership migration.

Reproducible commands, using the durable build and pinned package:

```zsh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
ninja -C /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement

O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement/stage/tests/o2-test-itsmft-tracking-l4-timeframe-configuration
```

The full serial ITS/MFT CTest result was **100% passed, 0 failed out of 100**;
all registered tests executed and there were no `Not Run` entries. No physics,
output, or workflow-default change is authorized in L4; replay comparison is
against the L3 parent.

```zsh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
ctest --test-dir /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
  -L itsmft --output-on-failure -j1
```

## Replay and output evidence

The 43-entry fixture manifest passed before and after replay with:

```zsh
cd /Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717
shasum -a 256 -c checksums.sha256
```

The fresh artifacts are under
`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l4-timeframe-configuration-20260806/`.
The standalone replays used the checked-in ITS and MFT scripts with run
`303000`, condition-not-after `1784207296000`, the fixed diamond `(0,0,0)`,
`pvRes=0.05`, and one MFT thread. The combined replay used both cluster
readers and `ITSMFTCombinedCATrackerParam.enabled=true` with the same settings.
The first combined invocation omitted the copied HBFUtils configuration and
failed during workflow setup; the replay was rerun after copying the fixture
configuration as required by the established protocol. The successful run
exited all eight devices with status zero.

| Leg | Tracks | Content hash |
| --- | ---: | --- |
| ITS standalone | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT standalone | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT combined | 68 | `8106b08571ca593c6b76ff72b761a680` |

The established `compare_common_ca_output.C` reported field-by-field equality
for standalone versus combined ITS and MFT products. Its MFT comparison covered
2,992 forward-state/covariance values with maximum absolute and relative delta
zero. The same comparator reported initialized-content equality for both
standalone and combined L4 products against the accepted L3 parent products,
including CommonTracks, cluster references, ROFs, labels, sidecars, and writer
fields. The comparison excludes only the known undefined
`MFTTrack.mInvQPtSeed` byte artifact; no initialized-content difference was
observed.

No CUDA or HIP compiler was available in the pinned environment (`nvcc` and
`hipcc` were absent), so no device-build result is claimed. The named user
stash `stash@{0}` (`user WIP: TripletTrackingRAndD.md`) remained untouched.

## Boundary retained for L5

Physical `SurfaceTrackingScratch`, its event reset, loader targets, and
adapter-side timing/publication state remain outside `TimeFrame`. L5 may move
workspace/reset ownership only after proving atomic load/reset behavior and
preserving raw-ROF and publication ownership at the workflow edge.

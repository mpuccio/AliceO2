# L6 direct `TimeFrame` loader validation

Status: historical L6 validation. P1 later removed the one-method
`MultiSourceTimeFrameLoader` façade, and the TimeFrame ownership cleanup later
superseded this transaction contract with direct reset-and-fill loading.

L6 makes `MultiSourceTimeFrameLoader::load(TimeFrame&, ...)` the direct,
non-owning loading component. It stages all decoded normalized input and
workspace backfills before one `TimeFrame` event commit. The live frame is not
changed by a decode, source-qualification, allocator, or capacity failure.
`TimeFrame::resetTimeFrame()` is invoked once by a successful replacement. Runtime
ROF views are copied into frame-owned event state, while raw ROFs and the table
storage behind those views remain owned by the adapter/workflow.

## Provenance

| Item | Value |
| --- | --- |
| Branch | `codex/itsmft-postm7-l6-direct-loader` |
| L5 base | `3c4e48e3e7` |
| Package | `daily-20260717-0700-local1` |
| Source | `/Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement` |
| Durable build | `/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement` |
| CMake source cache | `CMAKE_HOME_DIRECTORY:INTERNAL=/Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement` |
| Fixture | `pp-20ev-run303000-seed20260716-daily20260717` |
| Run / condition timestamp | `303000` / `1784207296000` |
| Diamond / PV resolution | `(0,0,0)` / `0.05` |
| MFT threads | `1` |
| User stash | `stash@{0}` remains the untouched `TripletTrackingRAndD.md` WIP stash |

The durable external artifacts are under:

`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l6-direct-loader-20260806/`

The final replay outputs are in `its-standalone-rerun/`,
`mft-standalone-rerun/`, and `combined-rerun2/`. The first ITS attempt in
`its-standalone/` was retained as diagnostic evidence: it used a stale
standalone path before the adapter was migrated to the direct loader and
produced an empty writer tree. It was not used for acceptance; the fresh
rerun after that scoped correction is the result below.

## Direct-loader contract and guards

`MultiSourceTimeFrameLoader::load()` now performs the following sequence:

1. reject an unconfigured frame;
2. decode and validate every source into a local `MultiSourceFrame`;
3. validate dense source count, source-qualified binding order, allocator
   identity, and plan-sized workspace capacity;
4. backfill one local workspace per source using the existing normalized
   decoder path; and
5. call the one `TimeFrame::commitLoadedEvent()` operation.

The commit resets old event data once and swaps only validated event-produced
workspace data. Standalone `TrackingInterface::loadTimeFrame()` and the
combined workflow both construct `ClusterSourceInput` at their adapter edge
and call this same operation. No production or test source in the scoped
tracking/workflow paths names `LoadTarget`, `loadTarget()`,
`AtomicLoadBinding`, or `loadEvent()`. The focused
`testMultiSourceTimeFrameLoader` guard scans the common include/src/test paths
for those deleted seams.

The direct-loader tests cover three-source atomic installation, partial source
failure without mutation, retry/replacement, unconfigured frames, invalid
source qualification, configured workspace ownership, reset-count behavior,
and combined ITS/MFT isolation. Existing lifecycle and adapter failure tests
cover recoverable drops, structural failures, exception cleanup, and runtime
ROF-view invalidation.

## Build and test commands

The affected library, tests, standalone workflows, and combined workflow were
built in the pinned environment with the existing build directory:

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

The focused loader, normalized-source, failure-contract, composition, ITS,
and MFT migration tests passed 6/6. The complete registered suite was then
run serially:

```zsh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
ctest --test-dir /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
  -L itsmft --output-on-failure -j1
```

Final result: **101/101 registered tests passed, 0 failed, and no `Not Run`
entries**.

## Replay protocol and results

The checked-in replay scripts were used with the provenance above:

```zsh
FIXTURE_DIR=<fixture> REPLAY_DIR=<artifact>/its-standalone-rerun \
TIMESTAMP=1784207296000 RUNNUMBER=303000 O2_BUILD_DIR=<build> \
O2_PACKAGE=daily-20260717-0700-local1 run-in-o2-env.zsh -- bash \
Detectors/ITSMFT/common/tracking/test/gate3-slice3-its-ca-validation/replay_tracking_its_common_ca.sh

FIXTURE_DIR=<fixture> REPLAY_DIR=<artifact>/mft-standalone-rerun \
TIMESTAMP=1784207296000 RUNNUMBER=303000 MFT_CA_NTHREADS=1 O2_BUILD_DIR=<build> \
O2_PACKAGE=daily-20260717-0700-local1 run-in-o2-env.zsh -- bash \
Detectors/ITSMFT/common/tracking/test/gate0-baseline/replay_tracking_common_ca.sh

o2-its-cluster-reader-workflow --with-mc --input-dir <fixture> --its-cluster-infile o2clus_its.root \
  | o2-mft-cluster-reader-workflow --with-mc --input-dir <fixture> --mft-cluster-infile mftclusters.root \
  | o2-itsmft-combined-ca-tracker-workflow -b --run \
      --condition-not-after 1784207296000 \
      --configKeyValues 'ITSMFTCombinedCATrackerParam.enabled=true;ITSCommonCATrackerParam.useDiamond=true;ITSCommonCATrackerParam.diamondPos[0]=0;ITSCommonCATrackerParam.diamondPos[1]=0;ITSCommonCATrackerParam.diamondPos[2]=0;ITSCommonCATrackerParam.pvRes=0.05;MFTCATrackerParam.nThreads=1'
```

The fixture manifest passed **43/43** before replay and **43/43** afterward.
The metric extractors reported:

| Leg | Input clusters | Input ROFs | Tracks | Candidate hash |
| --- | ---: | ---: | ---: | --- |
| ITS standalone | 4103 | 2304 | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS combined | 4103 | 2304 | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT standalone | 1249 | 2304 | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT combined | 1249 | 2304 | 68 | `8106b08571ca593c6b76ff72b761a680` |

The established `compare_common_ca_output.C` comparator reported field-level
equality for current standalone versus combined ITS and MFT products, and for
each current L6 output versus its L5 parent output. The MFT comparisons each
covered 2,992 forward float values with maximum absolute and relative delta
zero. This checks initialized CommonTracks, references, ROFs, labels, and
writer fields; only the known undefined `MFTTrack.mInvQPtSeed` byte artifact is
excluded. ROOT file bytes themselves are not a semantic comparison because
ROOT metadata varies between runs.

No CUDA or HIP validation is claimed: neither `nvcc` nor `hipcc` is available
in the pinned environment.

## Exact L7 boundary

L7 may relocate fixed ITS/MFT ROF table construction and backing storage,
publication clocks/validity, typed sidecars, and writer compatibility fully to
workflow/application owners. L6 leaves those owners and raw-ROF lifetime
unchanged. L7 must not add a generic timing manager, move raw ROFs into
`TimeFrame`, or alter output adapters, workflow reset decisions, tracking
physics, or frozen legacy ITS code.

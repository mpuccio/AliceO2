# Header and graph cleanup: Campaign A validation

- Status: accepted
- Date: 2026-08-07
- Base: `86113f2a14`
- Accepted head: `5bfcfae910`
- Package: `daily-20260717-0700-local1`
- Integration build: `O2-worktree-builds/itsmft-header-graph-cleanup`
- Frozen build: `O2-worktree-builds/itsmft-header-graph-cleanup-frozen-86113f2a14`

## Integrated slices

| Slice | Commit | Result |
|---|---|---|
| P1 scheduling and binding | `847c7179ce` | traversal validation and deterministic schedules moved into `SurfacePlanBinding` |
| P2 kernel parameters | `89c1194d72` | one device-portable `TrackingKernelParameters` record keyed by `SurfaceKind` |
| P3 stage operations | `49b2be6f9c` | explicit cylinder/disk tracklet and cell leaves in `detail/` |
| P4 final retirement | `5bfcfae910` | policy tags, headers, dispatch names, and compatibility surface removed |

Each slice passed `git diff --check`, pinned-environment
`git clang-format --diff`, changed-header compilation, focused Ninja targets,
and source-name-focused tests before integration. The final source guard found
no retired policy identifier, dispatch wrapper, compatibility alias, or one of
the five retired filenames under `Detectors/ITSMFT`.

## Build and test gate

Both the integration head and detached base were built with the common
tracking library, ITS, MFT, and combined workflow libraries, the CA writer,
and the Framework analysis and CCDB support libraries. The integration build
also built every directory contributing an `itsmft`-labelled test.

The required serial label run completed with no missing executable and no
`Not Run` entry:

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 93
Total Test time (real) = 119.45 sec
```

The strict preflight passed without a token-check bypass. It resolved the
requested package exactly, found a valid AliEn token, confirmed
`ALICEO2_CCDB_NOTOKENCHECK` was unset, and validated all reported Geant4 data
sets.

No CUDA or HIP compiler is available in the pinned environment. The campaign
therefore claims the exercised host/device-guard compilation only and makes no
real GPU-build claim.

## Frozen replay gate

Fresh feature and detached-base replays were created outside all worktrees at
`/private/tmp/o2-itsmft-header-graph-cleanup-gate-a-20260807`. The fixed
fixture manifest passed 43/43 checksums before and after the six replays.

| Leg | Input clusters | Input ROFs | Tracks | Content hash |
|---|---:|---:|---:|---|
| ITS feature standalone | 4103 | 2304 | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS feature combined | 4103 | 2304 | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS frozen standalone | 4103 | 2304 | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS frozen combined | 4103 | 2304 | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT feature standalone | 1249 | 2304 | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT feature combined | 1249 | 2304 | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT frozen standalone | 1249 | 2304 | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT frozen combined | 1249 | 2304 | 68 | `8106b08571ca593c6b76ff72b761a680` |

The complete extractor JSON was byte-identical for feature versus frozen and
standalone versus combined outputs. Eight persisted-product comparisons then
covered both detector legs across both dimensions. Every comparison reported
field-by-field equality for initialized writer products, ROFs, cluster-index
vectors, labels, references, `CommonTrack` content, and sidecars. Each MFT
comparison covered 2,992 projected float values with maximum absolute and
relative delta zero. The comparator retains only the previously documented
exclusion for the undefined `MFTTrack.mInvQPtSeed` byte.

The named user stash remained the unchanged object
`6a90bbcd7e187673a7eeaedc2f8df07c471c09b4` throughout the gate.

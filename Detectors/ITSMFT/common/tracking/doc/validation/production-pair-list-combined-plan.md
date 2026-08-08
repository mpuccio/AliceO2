# Production pair-list and combined-plan validation

- Status: review candidate after standalone/combined parity restoration (host build; no CUDA compiler available)
- Date: 2026-08-08
- Design: [flat pair-list graph and one combined plan](../design/0013-flat-pair-list-combined-plan.md)
- Package: `daily-20260717-0700-local1`
- Integration build: `O2-worktree-builds/itsmft-header-graph-cleanup`

## Implemented contract

- Production graph authority is one ordered surface list plus one flat list of
  immediate forward pairs.
- The combined workflow has one graph, binding, workspace, tracker, traits
  object, schedule, and tracker run. Family parameters are keyed only by
  `SurfaceKind`; they do not form source plans.
- The missing `(6,7)` pair defines the ITS/MFT component boundary without a
  subgraph or source-plan API.
- The default combined start and seeding masks are `0x1ffff`.
- Source IDs remain measurement/error/refit/publication provenance only.

## Focused evidence

The pair-list parity, graph builder, binding, loader, normalized-source,
workspace, detector-layout, tracklet, and combined-composition tests pass in
the pinned build. The combined composition fixture produces nonzero ITS and
MFT tracks and matches both standalone fixture counts. It also checks the
17-surface masks before execution and one shared 17-surface/15-transition/
13-cell workspace. A mixed-ROF-length regression verifies that timestamp
clamping for one disconnected component depends only on surfaces hit by that
track.

The focused combined-composition executable passes all 15 cases after the
family-semantics correction. The complete serial label suite was rebuilt and
passes with zero failures and zero `Not Run`:

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 93
```

The campaign build matrix covers the common tracking library, ITS and MFT
tracking/workflow libraries, the combined workflow, affected writers and
tests, the ITS tracking interface, and the Framework analysis/CCDB support
libraries. Changed-header onboarding/self-containment checks pass.

## Preflight and fixture integrity

The strict simulation/CCDB preflight passed in the pinned environment with a
valid token and all required Geant4 datasets; no token-check bypass was used.
The fixture manifest contains 43 files and verifies before and after replay:

```text
/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/
  pp-20ev-run303000-seed20260716-daily20260717/checksums.sha256
43/43 checksums passed
```

The rejected 102-track artifacts are under
`/private/tmp/o2-itsmft-flat-combined-gate-20260808/feature/combined-fixed2`.
The final parity-correction replay is under
`/private/tmp/o2-itsmft-mft-parity-fix-20260808/combined-final`.

## Replay results

Standalone output is unchanged from frozen commit `86113f2a14`:

| workflow | tracks | content hash | persisted output |
|---|---:|---|---|
| ITS standalone | 212 | `46913a67a7e2fe7462e29df0db264fa8` | field-for-field match |
| MFT standalone | 68 | `8106b08571ca593c6b76ff72b761a680` | field-for-field match; 2,992 projected values, zero delta |

The rejected combined plan produced:

| component | tracks | content hash | comparison |
|---|---:|---|---|
| ITS | 212 | `46913a67a7e2fe7462e29df0db264fa8` | field-for-field match with frozen combined ITS |
| MFT | 102 | `07e5addd74f374d2a4a70da69d0147aa` | deterministic but not accepted as a physics baseline |

The 102-track result came from applying ITS scalar defaults to disk leaves.
Compared by ordered cluster-reference sequence with standalone/frozen MFT, 66
tracks were retained, 36 sequences appeared, and two disappeared. Of the 36
appearing sequences, 32 have four clusters and are directly admitted by
`MinTrackLength=4` instead of MFT's established value 5. The remaining sequence
changes are exact: one old five-hit sequence is the suffix of a new ten-hit
extension, and one old six-hit sequence is replaced by the same suffix without
its first reference. Those substitutions do not change the count. The final
net +2 comes from one unrelated six-hit sequence and one unrelated ten-hit
sequence. The latter is the only raw addition with fitted `|x| < 0.05`,
consistent with the combined record also changing MFT's
`TrackletMinAbsX=0.05` to zero, although persisted fitted X alone does not
identify the earlier cluster/candidate rejection stage. Additions occur in
every populated ROF (`15,3,10,1,7` raw additions across ROFs 0--4), so this is
not a timestamp- or single-ROF-local population.

The chosen implementation direction is to restore standalone/combined equivalence, not to
nominate 102 as a physics candidate baseline. The correction preserves the
one plan and `StartLayerMask=0x1ffff`, but selects the established scalar and
LUT parameters by `SurfaceKind`. A fresh combined replay now gives:

| component | tracks | content hash | comparison |
|---|---:|---|---|
| ITS | 212 | `46913a67a7e2fe7462e29df0db264fa8` | field-for-field match with frozen combined ITS |
| MFT | 68 | `8106b08571ca593c6b76ff72b761a680` | field-for-field match with frozen combined and standalone MFT; 2,992 projected values, zero delta |

The ROF correction is independent of the track-count regression. Before the
shared-workspace change, each detector workspace's global ROF view was local
and correct. In one workspace the legacy global accessor would retain the
first source's ITS view, use ITS timing for MFT, and index a seven-layer view
with global MFT positions 7--16. The per-surface lookup resolves each hit's
owner view and local layer. The regression now uses an 80-BC MFT ROF against a
40-BC ITS ROF and compares both combined components with their standalone
timestamps. This is discriminating: the old global clamp would shorten MFT to
the ITS half-ROF, while the corrected path matches standalone MFT.

All initialized writer products, ROFs, labels, references, `CommonTrack`
records, and compatibility sidecars match within the applicable comparison.
The only excluded byte remains the previously documented undefined
`MFTTrack.mInvQPtSeed` byte.

## Static checks and limitation

`git diff --check` and pinned-environment `git clang-format --diff HEAD` are
clean. Source guards require the flat one-plan API and reject the retired
subgraph/source-plan and traversal-operation adapter surfaces.

`CMAKE_CUDA_COMPILER` is `NOTFOUND` in this build. Host/device-guard compilation
passes, but this validation makes no real CUDA or HIP build claim.

# Production pair-list and combined-plan validation

- Status: accepted (host build; no CUDA compiler available)
- Date: 2026-08-08
- Design: [flat pair-list graph and one combined plan](../design/0013-flat-pair-list-combined-plan.md)
- Package: `daily-20260717-0700-local1`
- Integration build: `O2-worktree-builds/itsmft-header-graph-cleanup`

## Implemented contract

- Production graph authority is one ordered surface list plus one flat list of
  immediate forward pairs.
- The combined workflow has one graph, binding, parameter record, workspace,
  tracker, traits object, schedule, and tracker run.
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

The focused combined-composition executable passes all 15 cases. The complete
serial label suite passes 93/93 with zero failures and zero `Not Run`:

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

Fresh artifacts are under
`/private/tmp/o2-itsmft-flat-combined-gate-20260808/feature`. The final
timestamp-corrected combined replay is `combined-fixed2`.

## Replay results

Standalone output is unchanged from frozen commit `86113f2a14`:

| workflow | tracks | content hash | persisted output |
|---|---:|---|---|
| ITS standalone | 212 | `46913a67a7e2fe7462e29df0db264fa8` | field-for-field match |
| MFT standalone | 68 | `8106b08571ca593c6b76ff72b761a680` | field-for-field match; 2,992 projected values, zero delta |

The combined plan has the following deterministic result:

| component | tracks | content hash | comparison |
|---|---:|---|---|
| ITS | 212 | `46913a67a7e2fe7462e29df0db264fa8` | field-for-field match with frozen combined ITS |
| MFT | 102 | `07e5addd74f374d2a4a70da69d0147aa` | two fresh runs match field-for-field; 4,488 projected values, zero delta |

The MFT combined count is intentionally not the old 68-track baseline. The
requested `StartLayerMask=0x1ffff` makes every MFT surface a valid road start,
and the one combined parameter record inherits ITS scalar defaults. Those are
explicit semantic changes to the combined default, so requiring old
standalone/combined track parity would contradict the chosen design. ITS has
no such content change.

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

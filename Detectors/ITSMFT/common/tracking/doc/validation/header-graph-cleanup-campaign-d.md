# Header and graph cleanup: Campaign D validation

- Status: accepted
- Date: 2026-08-07
- Frozen base: `86113f2a14`
- Accepted implementation head: `5317d16d1e`
- Package: `daily-20260717-0700-local1`
- Integration build: `O2-worktree-builds/itsmft-header-graph-cleanup`
- Frozen build: `O2-worktree-builds/itsmft-header-graph-cleanup-frozen-86113f2a14`

## Comment cleanup and inventory enforcement

Three comment-only slices removed migration history and obvious narration while
retaining ownership, lifetime, device, failure, numerical, and physics
invariants. A coordinator review removed the remaining milestone narration.
The inventory/source guard enforces the 36-public/12-detail endpoint, absence
of retired paths and policy symbols, workflow isolation from detail-only
workspace and compatibility headers, and detector-neutral generic plan
authority.

| Slice | Commit |
| --- | --- |
| inventory/source guard | `e061b871bd` |
| workflows, decoding, timing, and output comments | `c440354c0c` |
| Tracker, Traits, TimeFrame, binding, and Loader comments | `be9cb7ef3a` |
| core, device, and graph comments | `847b1efcbc` |
| coordinator comment review | `23258666a8` |

Each slice passed `git diff --check`, pinned-environment
`git clang-format --diff`, an affected Ninja build, and focused tests after
its separate integration. The final whole-branch diff and format audit from
`86113f2a14` was clean.

## Test-only pair-list evidence

Commit `5317d16d1e` adds only test and validation files. No production header,
source, graph representation, `SurfaceGraphBuilder`, `SurfacePlanBinding`,
`TimeFrame`, tracker behavior, default, or hole policy changed.

The prototype derives expanded transitions and skipped-surface witnesses,
transition and cell IDs, successor CSR, scheduled cells, and road starts from
ordered active surfaces, immediate pair declarations, an independent hole
policy, and an independent seeding mask. Raw-byte or explicit value parity
with the production graph passed for ITS, MFT, disconnected combined ITS+MFT,
non-monotonic IDs, empty topology, zero/one/multiple holes, hole-mask
rejection, and varied seeding masks. Tests also demonstrate that active
selection, adjacency, hole policy, and seeding eligibility remain independent
inputs.

The focused integration suite passed 4/4. Counts, ordering evidence, and the
measured pair-declaration versus current-topology byte counts are recorded in
[the pair-list prototype validation](campaign-d-pair-list-prototype.md). The
prototype is evidence only and makes no production representation
recommendation.

## Complete Gate D validation

The integration build covered common tracking, ITS tracking and workflow,
MFT workflow, combined workflow, CA writer, tracking interface, Framework
analysis/CCDB support, and every directory contributing an `itsmft`-labelled
test.

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 94
No Not Run entries.
```

The strict preflight passed with a valid AliEn token,
`ALICEO2_CCDB_NOTOKENCHECK` unset, the exact pinned package, and every reported
Geant4 data set present. The integration configuration reports
`CMAKE_CUDA_COMPILER=NOTFOUND` and no HIP compiler, so no real GPU build is
claimed; host/device guards are covered by the compiled suite.

## Frozen replay gate

Fresh outputs were created under
`/private/tmp/o2-itsmft-header-graph-cleanup-gate-d-20260807`. The fixture
manifest passed all 43 checksums before and after replay.

| Leg | Tracks | Content hash |
| --- | ---: | --- |
| ITS feature standalone and combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS frozen standalone and combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT feature standalone and combined | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT frozen standalone and combined | 68 | `8106b08571ca593c6b76ff72b761a680` |

All eight extractor JSON comparisons were byte-identical across the
feature/frozen and standalone/combined dimensions. Eight persisted-product
comparisons reported field-by-field equality for initialized writer output,
ROFs, cluster indices, labels, references, `CommonTrack` content, and
sidecars. Each MFT comparison covered 2,992 projected values with maximum
absolute and relative delta zero. Only the established undefined
`MFTTrack.mInvQPtSeed` byte remains excluded.

The named user stash remained the unchanged object
`6a90bbcd7e187673a7eeaedc2f8df07c471c09b4`. Gate D is the campaign stop line:
production pair-list authority, mixed cylinder/disk edges, and dynamic hole
traversal remain deferred to a later user decision and physics campaign.

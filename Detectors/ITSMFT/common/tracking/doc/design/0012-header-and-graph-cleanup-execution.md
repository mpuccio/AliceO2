# Header and graph cleanup execution

- Status: Campaign D accepted; execution complete at the test-only prototype stop line
- Date: 2026-08-07
- Integration branch: `codex/itsmft-header-graph-cleanup`
- Base: `86113f2a14`
- Architecture source: [header intentionality and surface-graph audit](0011-header-intentionality-and-surface-graph-audit.md)
- Follow-on decision: [production flat pair-list and one combined plan](0013-flat-pair-list-combined-plan.md)
- Pinned O2 package: `daily-20260717-0700-local1`
- Integration build: `O2-worktree-builds/itsmft-header-graph-cleanup`
- Physics status: behavior-preserving structural cleanup only
- Named user stash at campaign start: `6a90bbcd7e187673a7eeaedc2f8df07c471c09b4`

## Scope and stop line

This campaign executes the corrected audit in dependency order. It retires
the temporary transition-policy layer, consolidates accidental header
fragments, moves implementation-only facilities to `detail/`, cleans stale
comments, and adds a test-only pair-list graph derivation prototype. It does
not change the production graph representation, hole behavior, physics
formulas, cuts, defaults, ordering, output bytes, workflow ownership, or
detector compatibility.

The target inventory is 36 public and 12 detail headers when the downstream
`SeedAnchor` evidence gate is empty. If downstream evidence is unavailable or
non-empty, `SeedAnchor` remains and the endpoint is 49 headers.

## Execution order

### Campaign A: transition-policy retirement

1. P1 moves cycle/rank validation and deterministic binding schedules into
   `SurfacePlanBinding`, removes the redundant expected-policy declaration,
   tag grouping, dispatch, and test-only no-binding traversal fallback.
2. P2 replaces policy parameter records with one device-portable
   `TrackingKernelParameters` record keyed by `SurfaceKind`, bound once during
   `TrackerTraits::initialiseTimeFrame`. Unwired material-correction preflight
   code is deleted only after repository and downstream evidence.
3. P3 replaces policy operations with structurally parallel cylinder/disk
   leaves in `detail/TrackletFinding.h` and `detail/CellFinding.h`, selected
   once from endpoint kind. Existing state, propagation, and native-refit
   arithmetic stays in its current owner.
4. P4 converts `IndexTableConfiguration` to `SurfaceKind`, deletes all policy
   files/names/tests, and installs a zero-policy source guard.

Campaign gate A runs before any header consolidation.

Campaign A passed its complete build, test, authentication, fixture, and
frozen-replay gate. The exact evidence is recorded in
[Campaign A validation](../validation/header-graph-cleanup-campaign-a.md).

### Campaign B: device/public consolidation

- Merge `SurfaceMeasurementIndex` into `SurfaceId` and `LayerMask` into
  `LayerMask`.
- Merge `SurfaceLinearizationReference` into `SurfaceKinematicState` and
  `RefitLegAssembly` into `RefitDriver`. Delete `SeedAnchor` only if the
  repository and downstream search is available and empty.
- Merge catalog view into descriptor, static descriptor into surface spec,
  and ITS/MFT static specs into `ITSMFTDetectorDefinitions`.

There are no forwarding headers or compatibility aliases. Each focused
commit is validated and integrated separately before campaign gate B.

Campaign B passed its complete build, test, authentication, fixture, and
frozen-replay gate. The exact evidence is recorded in
[Campaign B validation](../validation/header-graph-cleanup-campaign-b.md).

### Campaign C: host boundaries and detail relocation

Merge clock publication and ROF timing validation into their owners. Build
one host-only cluster-decoding boundary while retaining live `IOUtils` APIs,
then consolidate the Loader declarations and failures. Move the audited
implementation-only headers, including the TimeFrame-owned workspace, to
`detail/` without changing ownership or atomic replacement. Decoding precedes
Loader work. `SurfacePlanBinding` remains a detail implementation type.

Campaign C passed its complete build, test, authentication, fixture, and
frozen-replay gate. The exact evidence is recorded in
[Campaign C validation](../validation/header-graph-cleanup-campaign-c.md).

### Comment and inventory campaign

Delete migration history and obvious narration while preserving ownership,
lifetime, device, numerical, failure, and physics invariants. Add a source
guard for the final inventory, retired paths/symbols, detail-only workflow
includes, and generic-core detector/layer-count authority.

The comment-only slices and the final inventory/source guard were integrated
without behavioral changes. The enforced endpoint is exactly 36 public and 12
`detail/` headers.

### Campaign D: test-only pair-list prototype

Add a non-production helper that accepts ordered active components,
immediate pairs, separate hole policy, and separate seeding mask. It derives
expanded pairs and witnesses, stable transition and cell IDs, successor CSR,
scheduled cells, and road starts. Tests compare every derived sequence and
byte representation with the current graph for ITS, MFT, disconnected
combined topology, non-monotonic IDs, empty topology, hole variants,
hole-mask rejection, and seeding variants. The four inputs—active selection,
hole policy, seeding eligibility, and adjacency—remain independent.

Campaign gate D is the stop point. Production pair-list authority, dynamic
hole traversal, and mixed cylinder/disk edges require a later decision and
physics campaign.

Campaign D passed its complete build, test, authentication, fixture, and
frozen-replay gate. The pair-list evidence and full gate record are in
[the prototype validation](../validation/campaign-d-pair-list-prototype.md)
and [Campaign D validation](../validation/header-graph-cleanup-campaign-d.md).

## Validation contract

Every implementation slice runs `git diff --check`, pinned-environment
`git clang-format --diff`, changed-header self-containment, affected Ninja
targets, and source-name-focused CTest filters.

Every campaign gate additionally rebuilds common tracking, ITS, MFT,
combined workflow, affected writer/test targets, and Framework analysis and
CCDB support; runs serial `ctest -L itsmft` with no failure or `Not Run`; runs
strict simulation/CCDB preflight; verifies all 43 fixture checksums before
and after fresh replay; and compares against the frozen `86113f2a14` build.
Required standalone results are 212 ITS tracks with hash
`46913a67a7e2fe7462e29df0db264fa8` and 68 MFT tracks with hash
`8106b08571ca593c6b76ff72b761a680`. Standalone/combined fields, initialized
writers, ROFs, labels, references, `GenericTrack`s, sidecars, and MFT projected
values must match. Only the documented undefined `MFTTrack.mInvQPtSeed` byte
is excluded.

Strict authentication or another external prerequisite failure stops the
campaign at that gate. Real CUDA/HIP absence is recorded as a limitation;
host/device-guard compilation is not reported as a GPU build.

## Coordination record

The coordinator alone modifies this reusable integration checkout. Isolated
Luna/high agents use temporary branches/worktrees and return focused commits;
no more than three run concurrently, and all `TrackerTraits`, binding,
kernel-operation, decoding, and Loader changes are serialized. Temporary
branches remain until final acceptance. Replay directories and fixtures stay
outside every worktree. The named user stash is never applied, modified,
dropped, or committed.

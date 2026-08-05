# Design note 0005: M7c runtime topology and ROF ownership

Status: M7c implementation and validation record
Date: 2026-08-05
Branch: `codex/itsmft-m7c-runtime-topology-rof`
Base: `dbaf26f52d` (integrated M7b equivalent)
Package: `daily-20260717-0700-local1`
Durable build: `/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement`

Companion plan: [GenericTrackingEngineMigration.md](../GenericTrackingEngineMigration.md)
Preceding design: [M7 runtime-plan Tracker/TrackerTraits migration](0003-m7-runtime-plan-tracker-migration.md)
Preceding implementation record: [M7b runtime-count authority](0004-m7b-runtime-count-authority.md)
Architecture decision: [ADR 0007](../decisions/0007-generic-tracking-engine-boundary.md)

This note records the completed M7c structural slice. It removes the
layer-indexed topology and detector-specific ROF-table ownership from the
common tracking path while deliberately leaving `Tracker<NLayers>` and
`TrackerTraits<NLayers>` for M7d. The slice preserves the native Propagator
path, CA policy and ordering, hole behavior, output adapters, raw-ROF
ownership, and workflow defaults.

## 1. Decision

The common tracking path has exactly one traversal topology and one ROF
boundary:

| Responsibility | Common-core authority | Owner of storage/lifetime |
|---|---|---|
| Surface descriptors and sparse transitions/cells | `DetectorLayoutView::topology`, a `SparseTrackingTopologyView` | `DetectorLayoutSet`/application plan |
| Which surfaces and topology entries a participant owns | `SurfacePlanBinding::getOrderedSurfaces()`, `getGlobalTransitions()`, and `getGlobalCells()` | The application participant |
| Compact event workspace | `SurfaceTrackingScratch::getNOwnedSurfaces()`, `getNTransitions()`, and `getNCells()` | The participant's scratch |
| Timing, ROF assignment, overlap, vertex lookup, and masks | `RuntimeROFViews` and its non-owning views | Adapter-owned fixed tables, configured per event |
| Raw ROFs, publication validity, clocks, and reset sequencing | DPL workflow state and the atomic loader | Combined workflow task |

`SurfaceTrackingScratch` no longer stores either a layer-indexed topology or
an ITS/MFT pair of ROF tables. It borrows the sparse topology during plan
initialisation and one `RuntimeROFViews` value for the current event. The
runtime views contain pointers and runtime counts; they do not own tables or
raw ROFs.

The adapter edge may still construct the frozen fixed-capacity
`ROFOverlapTable<N>`, `ROFVertexLookupTable<N>`, and `ROFMaskTable<N>` values
for an individual ITS or MFT leg. `TrackingInterface` and
`SurfacePlanTrackingParticipant` retain those owners and immediately publish
their views to the scratch. This is an explicitly named compatibility edge,
not a second core representation or a coordinator. Timing arithmetic is
shared through `ROFTimingLayer`; it is not duplicated in the core.

## 2. Sparse topology migration

The old common `TrackingTopology<NLayers>` header, source-level storage, and
test were deleted. Traversal now follows this sequence:

1. An application constructs a `DetectorLayoutSet` with static descriptors
   and a `SparseTrackingTopology` for each iteration.
2. `SurfacePlanBinding` validates ownership and retains the ordered surface,
   transition, and cell IDs. It also maps global topology IDs to compact
   scratch slots.
3. `TrackerTraits` validates and binds the sparse layout directly. It resolves
   transition endpoints through the binding's surface-position map and fills
   compact transition/cell scratch entries in binding order.
4. CA stages read the sparse transitions, cells, neighbours, scheduled cells,
   and road starts through the existing layout/binding views. No numeric
   `SurfaceId` range is interpreted as a traversal order.
5. `TrackSeed::SurfaceMask` records compact plan positions. A mask bit is not
   a global or detector-layer identity.

There is no second sparse topology in `SurfaceTrackingScratch`, and no
runtime-plan wrapper was introduced. Fixed device-portable value types retain
their established maximum capacities; their populated prefix is supplied by
the runtime binding/workspace counts.

The focused sparse test uses a non-identity topology with transitions
`5 -> 2 -> 7`, verifies the sparse IDs and cell membership, and the combined
composition tests retain the real ITS/MFT schedule and source-qualified
exports. Existing CA failure-contract tests continue to exercise fail-closed
checks for invalid topology references and structural tracking failures.

## 3. Runtime ROF view boundary

`ROFViews.h` defines the detector-neutral view layer:

- `ROFTimingLayer` contains the established timing arithmetic, including the
  diamond timing error envelope;
- `RuntimeROFOverlapView` exposes overlap, clock, and timestamp lookup;
- `RuntimeROFVertexLookupView` exposes vertex ranges and compatibility checks;
- `RuntimeROFMaskView` exposes enabled-ROF lookup; and
- `RuntimeROFViews` groups the four non-owning views, including the UPC mask.

The fixed ITS/MFT builders remain at the adapter edge because their frozen
capacity and table layouts are application compatibility details. The common
core sees only the view group. `SurfaceTrackingScratch::reset()` clears the
views and the UPC selection, so a view from an old event cannot be used after
load failure, dropped-frame reset, structural tracking failure, or event
replacement. A successful adapter load installs the new views only after the
event data is ready.

The event lifecycle is intentionally split as follows:

| Event transition | Combined DPL owner | Core/adapter consequence |
|---|---|---|
| Begin/replacement | Clears publication sidecars and invalidates clocks | The prior publication cannot be reused |
| Atomic load failure | Calls engine event reset and invalidates publication | `TimeFrame`, participant scratch, masks, and views are empty |
| Dropped time frame | Runs the existing reset path and invalidates publication | No timing or sidecar state is published |
| Structural tracking failure | Engine resets all participants and wipes the frame; workflow invalidates publication | No partial tracks, stale masks, or stale clocks escape |
| Successful replacement | Commits the new load, configures adapter views, then executes the ordered schedule | Publication is marked valid only after both legs succeed |

Raw ROF records remain workflow-owned. The live frozen ITS raw-ROF methods and
MFT raw-ROF `IOUtils`/workflow ownership remain outside the common CA path. The
dead common-scratch `loadROFrameData()`, `resetROFrameData()`, and
`prepareROFrameData()` family remains deleted; no live common caller was
found that would justify restoring it.

## 4. Deletions and source guards

M7c deletes or removes the following from the common path:

- `TrackingTopology.h` and its obsolete `testTrackingTopology.cxx` coverage;
- dual topology and ROF members in `SurfaceTrackingScratch`;
- templated scratch topology/ROF accessors and `checkSupportedNLayers()`;
- `IndexTableUtils<N>` compatibility alias, with callers using
  `IndexTableUtilsCore`; and
- layer-indexed topology setup from common tracker and onboarding paths.

`testM7cRuntimeTopologyRofGuard` scans common tracking `include/` and `src/`
and rejects `TrackingTopology<`, `ROFOverlapTable<`,
`ROFVertexLookupTable<`, `ROFMaskTable<`, `checkSupportedNLayers`, and
`IndexTableUtils<`. Its only production exclusions are the four named
adapter-edge files that own frozen ROF builders:
`TrackingInterface.{h,cxx}` and `SurfacePlanTrackingParticipant.{h,cxx}`.
Frozen legacy ITS sources are outside the scanned common tree and remain
unchanged. The guard also checks that the deleted common topology header does
not return.

The M7b count guard remains active. Its deferred category is now explicitly
`M7d Tracker/TrackerTraits boundary`; it no longer classifies a topology or
ROF bridge. Every remaining code-level `NLayers` use is classified as fixed
device ABI/capacity, private operation implementation, adapter edge, or that
M7d boundary.

## 5. Ranked simplification inventory

### Safe during M7c

1. Delete the common layer-indexed topology and obsolete topology test.
2. Delete the scratch's dual table members, template dispatchers, and
   `checkSupportedNLayers()`.
3. Delete `IndexTableUtils<N>` rather than preserving a source-compatibility
   alias.
4. Keep one runtime ROF view path and one scratch reset contract; do not add a
   second timing representation.
5. Convert traversal and seed-mask tests to sparse IDs, binding positions,
   and runtime views while retaining ITS/MFT schedule and export coverage.

### Blocked until M7d completes

1. Make `Tracker` and `TrackerTraits` non-templated and delete their common
   7/10 explicit instantiations.
2. Remove the `CATracker`/`TrackerTraits` detector aliases and the remaining
   `stateFamilyFromNLayers`/typed core bridge, resolving representation from
   actual surface descriptors and `SurfaceKinematicState` metadata.
3. Move the fixed ROF-builder setup out of the templated application seams
   once those seams can construct the single runtime core directly. The
   temporary fixed builders are permitted only at those adapter edges until
   their callers are migrated.
4. Reduce the remaining typed `TrackingInterface`,
   `SurfacePlanTrackingParticipant`, and `DetectorTraits` seams after the
   non-templated core has a single production caller. This is not a reason to
   reintroduce a common coordinator or a second tracker implementation.

### Deferred physics/algorithm decisions

The following are intentionally not M7c deletions: native Propagator and
covariance sanitization choices; tracklet/cell/road cuts and policy grouping;
holes and candidate ordering; cylinder/disk leaf arithmetic; refit acceptance
and output conversion; CommonTrack and writer semantics; and frozen legacy ITS
workflow behavior. Removing or simplifying any of these requires a separate
physics or algorithm sign-off and an A/B replay, even if a runtime-plan API
would make the code shorter.

## 6. Exact M7d boundary

M7d starts with the existing runtime plan and this M7c ownership boundary. It
may change only the common `Tracker`/`TrackerTraits`/`CATracker` template
composition and the explicit instantiations needed to remove that boundary.
It must not add a wrapper around `Tracker<NLayers>`, revive a layer topology or
ROF table in the core, change `TransitionPolicyTag` visibility, or alter
Propagator, covariance, CA, holes, output, raw-ROF, or workflow semantics.

The first M7d proof is a real deletion: no common production
`Tracker<7>`, `Tracker<10>`, `TrackerTraits<7>`, or `TrackerTraits<10>`
instantiation remains, and one shared non-templated algorithm body serves both
application plans. Adapter-only templates may remain only when they own
detector configuration/output compatibility and pass their own explicit
source guard.

## 7. Validation record

The durable build was reused from the requested worktree and package. The
affected common libraries, ITS/MFT application code, combined workflow, and
tests rebuilt successfully. The final serial selector executed all 95
registered `itsmft` tests with no `Not Run` entries; all 95 passed. The final
replay record below is the authoritative M7c result and retains the accepted
candidate hashes:

- fixture checksum manifest: 43/43 before and 43/43 after;
- standalone ITS: 212 tracks,
  `46913a67a7e2fe7462e29df0db264fa8`;
- standalone MFT: 68 tracks,
  `8106b08571ca593c6b76ff72b761a680`;
- combined ITS and MFT legs match their standalone products field-by-field;
- initialized writer content matches the M7b parent for both standalone and
  combined products, excluding only the known undefined
  `MFTTrack.mInvQPtSeed` byte artifact. The MFT comparison covered 2,992
  float-projected values with maximum absolute and relative delta zero;
- the fixture was unchanged: `shasum -a 256 -c checksums.sha256` passed 43/43
  both before and after replay; and
- no GPU/device result is claimed: `nvcc`, `hipcc`, `nvidia-smi`, and
  `rocminfo` are absent in the pinned environment.

Replay outputs and metric JSONs are in
`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/m7c-runtime-topology-rof/`.
The combined run is in `combined-retry/`; the initial out-of-environment
attempt in `combined/` produced no workflow output and is retained separately
as an execution trace, not as validation evidence. The source/dependency
guards passed, and both `git diff --check dbaf26f52d HEAD` and
`git clang-format --diff dbaf26f52d HEAD` are clean.

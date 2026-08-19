# Post-M7 intentionality and cleanup audit

- Status: audit complete; implementation not started
- Date: 2026-08-06
- Audited revision: `51a49e5d4a` (`codex/itsmft-integration`, after M7f and P1)
- Follow-up: revised after CPU/GPU traits and graph-initialization review
- Scope: `Detectors/ITSMFT/common/tracking/` production headers and sources,
  their production construction sites, and tests that constrain their shape

This is an architecture and ownership audit. It authorizes no C++, CMake,
test, workflow, algorithm, configuration, or physics change. In particular,
items below labelled “safe” still require their stated build and replay gate
when implemented.

## 1. Executive conclusion

M7 achieved one runtime-plan CA core, but the route used to reach it left two
application compositions and several migration-era seams around that core.
The core algorithm is not duplicated: `Tracker` and `TrackerTraits` are the
single production implementation. The duplication is primarily lifecycle,
adapter composition, and compatibility ownership:

- standalone ITS/MFT use `TrackingInterface<7/10>`, which owns loading,
  timing tables, scratch, tracker, publication compatibility, reset policy,
  and float-sentinel translation;
- combined ITS+MFT uses `SurfacePlanTrackingParticipant<7/10>`,
  `MultiSourceTimeFrameLoader`, and `TrackingEngine`, while the DPL task owns
  event timing/publication validity;
- `Tracker` resets its scratch and the shared `TimeFrame` on failure even
  when `TrackingEngine` is the whole-event reset owner, producing a real
  duplicate reset path;
- `TrackingOperationAdapter` mixes native-refit execution, accepted-result
  compatibility staging, and adapter reset;
- several public headers and APIs now have no production caller and survive
  only because migration-era tests name them.

The initial version of this audit incorrectly recommended merging
`TrackerTraits` into `Tracker`. The follow-up CPU/GPU review disproves that
recommendation. The frozen-but-live ITS architecture has a CPU
`TrackerTraits<NLayers>` implementation and a `TrackerTraitsGPU<NLayers>`
override of the same kernel entry points. The common class deliberately
retains the corresponding virtual kernel methods, `getName()`, and `isGPU()`.
No common GPU backend has been ported yet, but the seam is intentional, not
residue:

- `Tracker` is the target-independent algorithm orchestrator;
- `TrackerTraits` is the target-architecture kernel implementation; and
- a future common GPU traits implementation must be able to replace the CPU
  traits without replacing the orchestration.

The next follow-up also corrects graph/workspace ownership. `Tracker` should
not become a new aggregate root merely because engine/interface wrappers are
deleted. `TimeFrame` is the reusable tracking entity. It owns persistent
tracking configuration and allocator-backed state together with current-event
measurements/results; `Loader` and `Tracker` are components that act on a
borrowed `TimeFrame&`. `TrackerTraits` supplies the borrowed architecture
backend used by `Tracker`.

This recommendation depends on an observed lifecycle fact, not on the name
alone: standalone `TrackingInterface` and the combined DPL task each retain a
`TimeFrame` member across `run()` calls. Therefore graph/configuration can be
installed once and preserved while an event reset clears only measurements,
workspace contents, and results. If an implementation instead reconstructed
or recopied the graph for every event, it would fail this ownership target;
that would be a regression, not an alternative interpretation of the design.

The intended end state does not need a new wrapper or manager. It should:

1. delete demonstrably dead typed/refit and forwarding residue;
2. make the existing loader directly target `TimeFrame` and remove fixed
   ITS+MFT wrapper entry points and public scratch-target plumbing;
3. give `TimeFrame` one fully constructed immutable
   `std::vector<SurfaceGraph>`, parameters, partitions, allocator-backed
   workspaces, normalized input, and generic results;
4. move generic tracking execution and failure orchestration from
   `TrackingEngine`/`TrackingInterface` into `Tracker`, which borrows the
   configured `TimeFrame` and a selected traits backend, while leaving raw-ROF
   and publication lifecycle in the DPL owner;
5. delete `TrackingEngine`, `TrackingInterface`, and the participant
   coordinator vocabulary once `Tracker` directly consumes the explicit
   graph partitions/schedule;
6. narrow the operation seam to refit only and move publication compatibility
   fully to application/workflow code; and
7. preserve `TrackerTraits` as the CPU/GPU kernel strategy seam, making its
   backend contract more explicit rather than merging it away.

`TrackingEngine` and `TrackingInterface` should be retired. The engine is a
short schedule/reset loop around participant-owned `Tracker` instances; the
interface is a standalone composition and lifecycle holder. Their useful
tracking responsibilities belong in `Tracker`, but workflow responsibilities
do not: raw input ownership, timing construction, publication validity, and
writer adaptation remain at the DPL edge. `TrackerTraits` remains distinct as
the kernel backend. Neither Loader nor Tracker owns `TimeFrame` data.

The follow-up initialization audit also found a real model simplification
opportunity. `DetectorLayout` owns only a `SparseTrackingTopology`,
`DetectorLayoutSet` is mostly an iteration-vector wrapper, and
`DetectorLayoutConfigurationKey` is no longer a key. In combined tracking,
one authoritative topology is copied into separate ITS and MFT
`DetectorLayoutSet`s. The recommended vocabulary is `SurfaceGraph`, not
`DetectorGraph`: the graph connects surface nodes and may span detector
identities. Section 6.5 defines the proposed graph/plan model without
authorizing a rename implementation.

The first implementation slice was unambiguous: remove the unused common
typed-MFT refit/export path and migrate its fixture to the already-live
generic native-refit result path. L0 now records that deletion and its absence
guard; Section 8.1 remains the exact boundary and exit criterion.

## 2. Audit method and invariants

The audit traced declarations, definitions, constructors, members, mutation
sites, reset/load/publication calls, explicit instantiations, CMake entries,
production callers in ITS/MFT workflows, and tests. A repository search was
also used to distinguish common-library APIs from same-named frozen ITS/MFT
utilities.

Recommendations preserve these invariants:

- the generic core names no detector identity, source-position convention,
  DPL/workflow/writer type, or typed detector output;
- `TransitionPolicyTag` stays a private compatibility implementation detail;
- `SurfaceKind` remains representation metadata, never schedule or topology;
- the future `SurfaceGraphView` (today `DetectorLayoutView`), graph
  partitions, and TimeFrame-internal workspace remain the runtime
  authorities;
- raw ROFs and publication lifecycle remain workflow-owned;
- fixed capacities remain only where device-safe value representation needs
  them;
- cylinder/disk differences remain narrow descriptor/state operations;
- frozen `Detectors/ITSMFT/ITS/tracking/` is excluded explicitly, never by a
  blanket exclusion of common tracking; and
- no recommendation in this note changes Propagator, covariance, material,
  holes, cuts, ordering, or acceptance.

## 3. Ownership and lifecycle map

### 3.1 Current standalone path

```text
ITS or MFT DPL task
  owns raw ROFs and typed output staging
  |
  +-- TrackingInterface<7 or 10>
        owns TimeFrame + SurfaceTrackingScratch
        owns static plan + SurfacePlanBinding
        owns TrackerTraits + Tracker
        owns fixed ROF tables/masks + publication clock
        owns DetectorPublicationAdapter + compatibility sidecars
        |
        +-- configure timing/beam/masks
        +-- scratch.loadNormalizedSource(frame, ...)
        +-- Tracker::clustersToTracks(adapter = interface)
        |     +-- TrackerTraits: tracklets -> cells -> neighbours -> roads
        |     +-- adapter refit -> generic GenericTrack
        |     +-- adapter accepted-result compatibility
        +-- translate TrackingResult to float/drop sentinel
        +-- workflow consumes GenericTrack and sidecars
```

The DPL task owns the actual event boundary, but the interface owns most event
lifecycle decisions. On a load or tracking failure, the interface resets its
clock, sidecars, scratch, and frame. `Tracker` also resets scratch and frame
for tracking failures. The ownership is therefore operationally split.

### 3.2 Current combined path

```text
CombinedCATrackerSpec DPL task (whole-event owner)
  owns raw ROFs, publication clocks/validity, TimeFrame, layout, schedule
  |
  +-- SurfacePlanTrackingParticipant<7>  (ITS leg)
  +-- SurfacePlanTrackingParticipant<10> (MFT leg)
  |     each owns scratch, binding, TrackerTraits, Tracker,
  |     fixed ROF tables, publication adapter/sidecar, load target
  |
  +-- MultiSourceTimeFrameLoader::loadEvent()
  |     stages normalized frame and both scratch backfills atomically
  |
  +-- TrackingEngine::executeEvent(schedule)
        runs participants in explicit order
        resets every participant, then TimeFrame, on failure
```

The combined workflow correctly owns raw ROFs, clocks, validity, and event
publication. However, a participant calls `Tracker::clustersToTracks()`, and
that method clears its scratch and the shared `TimeFrame` before propagating a
failure. The engine then clears all participant state and the frame again.
This is duplicate lifecycle ownership and means the leaf tracker mutates
whole-event state before the engine applies its all-or-nothing policy.

### 3.3 Target lifecycle

```text
DPL task (standalone or combined)
  owns raw ROFs, timing construction, publication validity and event policy
  |
  +-- TimeFrame entity
  |     owns vector<SurfaceGraph> + per-iteration partitions/parameters
  |     owns allocator-backed workspaces
  |     owns normalized measurements + GenericTrack results
  |
  +-- Loader::load(TimeFrame&, normalized source inputs/views)
  |
  +-- Tracker::run(TimeFrame&, TrackerTraits&)
        orchestrates initialise -> tracklets -> cells -> neighbours -> roads
        sequences the generic tracking transaction
        invokes TimeFrame's event reset on tracking failure
        |
        +-- application publication adapter consumes successful generic result
```

Only the number of graph partitions stored by `TimeFrame` differs between
standalone and combined. A standalone frame has one; combined ITS+MFT has two
explicit source-qualified partitions. Neither needs an engine, participant
coordinator, interface, or float sentinel around `Tracker`.

This does not pull workflow ownership into `Tracker`. Raw compact clusters,
raw ROFs, CCDB-derived timing construction, publication clocks/validity,
typed output conversion, and writers remain outside. The tracker consumes
normalized event input and runtime timing/mask views.

## 4. Dependency map of the central cluster

Solid arrows below mean construction/ownership or a required call. “Adapter”
arrows are detector compatibility dependencies that should move outward.

```text
Combined DPL task -----------> TrackingEngine -----> TrackingParticipant
       |                                                |
       +--> MultiSourceTimeFrameLoader                  v
       |                                    SurfacePlanTrackingParticipant<7/10>
       |                                                |
       |                                                +--> Tracker
       |                                                       |
       |                                                       v
       |                                                  TrackerTraits
       |
Standalone DPL task --------> TrackingInterface<7/10>
                                  |
                                  +--> another Tracker + TrackerTraits composition

Target:

Standalone/combined DPL task
       |
       +--> TimeFrame
       |      owns vector<SurfaceGraph>, partitions, parameters, pool,
       |      workspaces, normalized event data, and GenericTrack results
       |
       +--> MultiSourceTimeFrameLoader::load(TimeFrame&, inputs)
       |
       +--> Tracker::run(TimeFrame&, TrackerTraits&)
              +--> CPU or GPU kernel backend
              +--> graph/partition order read from TimeFrame
              +--> generic event transaction over borrowed frame state
              +--successful result--> application publication adapter

TimeFrame <----- normalized measurements, GenericTrack results, and
                TimeFrame-owned graph partitions/bindings
ROFViews <----- borrowed by core; built and owned at adapter edge
TrackingOperationAdapter <----- refit + accepted compatibility + reset seam
```

The highest-cost dependencies are not algorithmic. They are the extra
engine/participant/interface composition around `Tracker`, the inclusion of
`TrackingInterface.h` by `SurfacePlanTrackingParticipant.h` solely to reuse
sidecar-owner mixins, and the public detector-specialized
`DetectorPublicationAdapter.h` beneath the nominally common include root.

## 5. Responsibility audit of the focal modules

| Module | Actual single responsibility and data relationship | Classification | Duplication and disposition | Risk |
|---|---|---|---|---|
| `Tracker.h` | Currently a forwarding include for `CATracker.h`, despite `Tracker` being the intended public class name. | Canonical API name obscured by residue | Move the actual declaration here, then delete the forwarding relationship. | Behavior-preserving header migration; compile all public-header users. |
| `CATracker.h` / `Tracker` | Target-independent orchestrator for configured iterations and CA stages. It currently borrows traits/frame/scratch/layout but also owns parameter and pool handles. | Generic core component in a stale file name | Move it to canonical `Tracker.{h,cxx}` and make it act on `TimeFrame&` plus `TrackerTraits&`. Parameters, graphs, partitions, pool, workspace, and event data move to TimeFrame; Tracker retains no persistent tracking entity. | Reset and multi-partition migration are exception/order-sensitive and need full failure/replay gates. |
| `TrackerTraits` | CPU kernel implementation and architecture strategy seam for initialisation, tracklets, cells, neighbours, roads, and associated backend caches. Its virtual API mirrors the live frozen ITS `TrackerTraitsGPU` override pattern. | Generic kernel backend | Keep distinct from `Tracker`. Make the CPU/GPU substitution contract explicit; do not merge kernels into the orchestrator. Reduce only non-kernel compatibility state. | A future common GPU port requires real device ABI/build validation; CPU behavior remains replay-gated. |
| `TrackingEngine` | Executes an explicit participant schedule and repeats whole-event reset policy around participant-owned trackers. It owns no data or kernel behavior. | Migration-era orchestration layer | Retire after `Tracker` reads the explicit graph-partition order from TimeFrame and applies the generic transaction. Do not rename it or move DPL concerns into Tracker. | Combined ordering, failure classification, all-or-nothing reset, and output isolation must be pinned before deletion. |
| `TrackingParticipant` | Dynamic wrapper around one plan-bound tracker composition: identity, surfaces, track/reset/export. | Migration-era application seam | Retire with the engine once partitions are immutable plan data consumed directly by `Tracker`. Preserve source-qualified partition order as data, not polymorphic coordinator objects. | Combined heterogeneous schedule and future-backend tests; avoid closing the plan over ITS/MFT types. |
| `TrackingInterface<7/10>` | Standalone application composition and event coordinator. Owns frame, scratch, tracker, plan/binding, decoder, ROF state, clocks, sidecars, and publication adapters. | Application adapter plus compatibility residue | Delete. Move generic configuration/state/reset mechanics into `TimeFrame`, orchestration into `Tracker`, and raw loading/timing/publication responsibilities outward to the workflow. Do not recreate an interface facade. | High lifecycle/output compatibility risk; stage under standalone and combined replay gates. |
| `SurfacePlanTrackingParticipant<7/10>` | Combined plan-bound leg. Owns scratch/binding/tracker composition, fixed ROF compatibility tables, load target, sidecars, and tracked flag; borrows plan/frame. | Migration-era application wrapper | Delete after graph/partition/workspace data moves into TimeFrame and timing/publication ownership moves to workflows. It need not become a non-templated permanent participant. | Bounded adapter migration; exact schedule, writer, and sidecar replay required. |
| `TrackingOperationAdapter` | Supplies seed refit, accepted-result completion, and adapter reset. It borrows generic candidate/scratch/measurement views. | Mixed core/adapter seam | Three responsibilities. Narrow to the one operation the algorithm needs (refit); run publication conversion/reset after generic success at participant/workflow edge. Delete if refit can be directly owned without detector dispatch. | Call order and rejection classification are physics-sensitive; preserve exact boundaries. |
| `DetectorPublicationAdapter<7/10>` | Specializes accepted-result publication compatibility for ITS shared-cluster flags and MFT sidecars. | Detector application adapter | Responsibility is real but location/visibility is wrong. Move to ITS/MFT application code and give each workflow direct ownership. No generic dispatcher. | Writer/sidecar compatibility gate. |
| `TimeFrame` | Reusable tracking entity. It currently owns normalized multi-source event data, vertices/labels, beam/Bz, pool, and generic results while scratch/graphs/parameters are owned elsewhere. | Generic aggregate root | Expand intentionally to own immutable graph vector, per-iteration parameters/partitions, allocator-backed workspaces, normalized inputs, and generic outputs. One reset primitive clears event state but preserves configuration/capacity. | Broad ownership migration; allocator, reset, CPU/GPU, and replay gates. |
| `SurfaceTrackingScratch` | Owns runtime-sized mutable CA workspace and per-source legacy-compatible measurement backfill. Borrows allocator/pool. | Split state container created for participant composition | Move its storage behind TimeFrame, either merged directly or retained as a private partition-workspace implementation detail. It should not remain a separately configured public entity. | Allocator identity, reset, atomic load, and sparse-plan tests. |
| `SurfacePlanBinding` | Immutable source-qualified projection from a global topology to one graph partition's ordered surfaces and compact transition/cell slots. | Generic graph-partition data | Keep the information, not necessarily the class/name. Make it per iteration and TimeFrame-owned; prefer `SurfaceGraphPartition` or `SurfaceGraphBinding`. | Sparse/non-identity and multi-iteration topology tests are mandatory. |
| `MultiSourceTimeFrameLoader` | Provides atomic normalized-frame plus scratch-backfill staging/commit. | Generic component acting on TimeFrame | Keep the transaction, remove virtual scratch-target hierarchy and fixed wrappers, and expose a direct `load(TimeFrame&, inputs)` boundary. | Atomicity/allocator/failure-retry tests and combined replay. |
| `ROFViews` | Defines non-owning runtime timing, overlap, vertex-lookup, and mask views consumed by common tracking. | Generic runtime input view | Keep. Fixed detector tables belong at adapter/workflow edge and should not migrate into scratch/frame. | Timing/diamond/mask parity tests. |

## 6. Duplicate lifecycle and representation findings

### 6.1 Initialization and construction

`TrackingInterface` and `SurfacePlanTrackingParticipant` each build the same
`TrackerTraits`/`Tracker`/scratch/binding/pool composition. Both also select
ITS/MFT fixed timing tables and publication compatibility through `NLayers`.
The single state composition should become `TimeFrame`, configured with one
immutable vector of iteration graphs, parameters, partitions, allocator, and
workspaces. `Loader` and `Tracker` act on it; `TrackerTraits` supplies the
selected backend. Adapter-owned timing/publication objects stay outside.

The follow-up found duplicate runtime topology ownership in the combined
path. `buildCombinedLayout()` creates one authoritative combined
`DetectorLayout`; `ownDetectorPlan()` then copies that layout, including the
owning sparse-topology vectors, into both `mITSPlan` and `mMFTPlan`. The two
plans differ only in their configuration keys while presenting the same
global graph. This should become one owned graph with two borrowed
partitions.

There is also a plan/binding lifetime mismatch hidden by current one-iteration
acceptance workflows. `DetectorLayoutSet` can own a different topology per
iteration, but standalone and combined construction build one
`SurfacePlanBinding` from `getLayoutView(0)` and reuse its transition/cell IDs
and schedules for every iteration. MFT async parameters can vary start masks
and holes by iteration. A general multi-iteration contract must therefore
either prove graph topology invariant or own a binding/partition per
iteration. The latter matches the data model and avoids a hidden restriction.

### 6.2 Loading

Three loading shapes remain:

1. `TrackingInterface::processTimeFrame()` loads one source directly through
   `SurfaceTrackingScratch::loadNormalizedSource()`;
2. `MultiSourceTimeFrameLoader::loadEvent()` performs the production atomic
   multi-source transaction;
3. `loadITSAndMFT()` is a fixed-position wrapper used only by tests, along
   with `resetITSAndMFTEvent()`.

The nested virtual `LoadTarget` hierarchy has one implementation and every
production target ultimately mutates TimeFrame plus a
`SurfaceTrackingScratch`. Once scratch is internal entity storage,
`Loader::load(TimeFrame&, inputs)` can resolve source-qualified partitions
from the configured frame, stage all normalized/workspace state privately,
and commit atomically. This deletes virtual dispatch, participant
`mLoadTarget`, public scratch targeting, friendship-driven forwarding, and
one lifetime constraint.

### 6.3 Reset and failure handling

Reset is currently implemented at four levels:

- `SurfaceTrackingScratch::reset()` and its forwarding alias
  `resetScratch()`;
- `resetTimeFrameEvent(frame, scratch)`;
- `TrackingInterface::resetEvent()` plus operation-adapter compatibility
  reset;
- `TrackingEngine::resetEvent()` plus participant `eventReset()`.

`TimeFrame` should own the reset operation because it owns all generic state.
`Tracker` decides that a tracking failure requires reset and invokes that
single entity operation; Loader staging failure remains non-mutating by
contract. The DPL owner then invalidates workflow publication state; it does
not duplicate generic tracking reset. This cannot be a casual cleanup:
exception paths, dropped-TF classification, retryability, and publication
invalidation must be pinned first.

### 6.4 Tracking and publication

There is one CA algorithm body, but two adapter entry paths. Publication is
also split: generic `GenericTrack` append occurs in the tracker path,
`TrackingOperationAdapter::completeAccepted()` stages compatibility, the
participant exposes a publication export, and the workflow builds typed
outputs. The intended split is simpler:

- `Tracker` produces generic accepted results and `GenericTrack` references;
- a detector application adapter converts the successful generic result;
- the workflow owns validity, clocks, raw ROFs, and final publication.

Adapter publication must never be called from a failed generic transaction.
Moving it outward requires preserving serial accepted-result order and the
current “final” boundary.

### 6.5 Surface graph model and naming

The current names obscure a simpler structure:

| Current type | What it actually represents | Finding |
|---|---|---|
| `StaticSurfaceDescriptor` / `SurfaceDescriptor` | Compile-time and runtime payload for a graph node. | Two representations are justified by constexpr authoring versus device/runtime ABI. Projection should remain one-way and construction-time. |
| `SparseTrackingTopology` | The actual graph: surface-to-surface transitions, three-surface cells, seeding nodes, and adjacency offsets. | This is the object whose public concept should be “surface graph.” |
| `DetectorLayout` | One owning sparse topology plus validation error; it does not own surface layout/geometry. | Forwarding wrapper. Merge with/rename the owning graph rather than retain both names. |
| `DetectorLayoutView` | Device-facing descriptors, kind masks, and sparse-topology view for one iteration. | This is naturally `SurfaceGraphView`. |
| `DetectorLayoutSet` | Shared catalog/masks plus a vector of iteration layouts and a stale configuration-key record. | After graph ownership is consolidated, its only live role is vector ownership; that does not justify a public plan class. |
| `DetectorLayoutConfigurationKey` | Original ordered surfaces and build inputs. It is never used as a currency/cache key; production reads only `orderedSurfaces`. | Delete the “key” abstraction. Store live traversal data directly in the plan/graph. Do not retain build inputs solely for equality. |
| `SurfacePlanBinding` | Source-qualified ordered graph partition plus global-to-compact mappings and precomputed transition/cell schedules. | Necessary information for partitioned tracking, but it must be associated with each iteration graph and can be named as graph partition/binding data. |

The recommended vocabulary is:

- **`SurfaceGraph`**, not `DetectorGraph`, for one iteration's immutable
  surface nodes and edges. `DetectorGraph` would reintroduce the false “one
  detector, one graph” assumption and is awkward for combined or future
  systems;
- **`SurfaceGraphView`** for the device-facing POD view;
- **`std::vector<SurfaceGraph>`** for the ordered iteration graphs, owned
  directly by `TimeFrame`; no public `TrackingPlan` wrapper; and
- **`SurfaceGraphPartition`** (or, if the mapping action remains central,
  `SurfaceGraphBinding`) for one source-qualified ordered subset and compact
  slot map.

This is a consolidation, not a new wrapper hierarchy. The implementation target
deletes `DetectorLayout`, `DetectorLayoutSet`, and
`DetectorLayoutConfigurationKey` vocabulary as their responsibilities move
into the irreducible graph and partition concepts, with iteration order
represented by the vector itself.
`SparseTrackingTopology` may remain a private storage component inside
`SurfaceGraph`; it should not remain a competing public graph model.

Parameters and partitions still have to be paired with the correct graph.
That invariant belongs to `TimeFrame` construction/configuration, which can
validate the caller-supplied spans and store a private per-iteration record
containing a graph, parameters, and ordered partitions. Such a private record
is entity storage, not a public `TrackingPlan` abstraction.

### 6.6 Initialization audit and target phases

Current construction permits many partially valid states:

1. resolve parameters;
2. construct and distribute a pool separately to `TimeFrame`, scratch,
   traits, and tracker;
3. construct traits, then tracker, then bind scratch/frame;
4. project a static catalog and build a layout set;
5. build one binding from iteration 0;
6. size scratch from that binding; and
7. separately set field, thread arena, ROF views, and publication adapters.

The order differs by path. Standalone installs the pool before
`scratch.adoptPlan()`. Combined calls `adoptPlan()` first, creating/resizing
containers with the scratch's current allocator state, and only later calls
`setMemoryPool()`, which clears and rebinds those containers. This is at least
redundant initialization and makes allocator correctness dependent on setter
order.

`TrackerTraits::initialiseTimeFrame()` then repeats both static and event work
for every event/iteration: graph grouping and policy selection, binding and
material validation, operation-function binding, event measurement-span/ROF
validation, index-table configuration, workspace clearing/allocation, and
cluster LUT preparation. Static graph/partition validation and operation
selection need not be redone per event; normalized measurement, ROF, LUT, and
workspace work genuinely does.

The target has three explicit phases:

| Phase | Owner | Work allowed |
|---|---|---|
| Configure once | `TimeFrame` | Consume parameters and one immutable `std::vector<SurfaceGraph>`; validate graph/parameter counts and every iteration partition; establish one allocator/pool; size internal workspaces from maximum required graph extents. No event data. |
| Load event | `MultiSourceTimeFrameLoader` acting on `TimeFrame&` | Decode normalized source-qualified measurements, consume workflow/adapter-built runtime timing and mask views, and atomically commit normalized input plus internal partition backfills. Raw ROFs and timing construction remain workflow-owned. |
| Execute event/iteration | `Tracker` acting on `TimeFrame&` through `TrackerTraits&` | Select the configured graph/partition, bind event measurement/ROF views, clear/reuse internal workspace, build/reuse LUTs according to pass flags, execute kernels in graph order, and commit generic results or invoke the frame's one reset operation. |

A `TimeFrame` constructor or one fallible `configure()` operation should
establish entity configuration atomically; the present sequence of `adopt*`
and `set*` calls should not be replaced by a different setter sequence.
The configured frame is then reused across DPL `run()` calls. Its event reset
must retain graphs, parameters, partitions, allocator identity, and reserved
capacity while clearing all event-derived state.
Backend selection remains explicit at execution: Tracker borrows one CPU or
GPU traits implementation whose lifetime covers the call; neither component
owns the frame.

## 7. Rest-of-public-include-tree audit

Every public common-tracking header is covered below. Grouping means the
listed files share one disposition; it does not imply a new module.

| Files/modules | Current owner and recommendation |
|---|---|
| `BarrelSurfaceStateOperations.h`, `ForwardSurfaceStateOperations.h`, `SurfaceStateOperationResult.h`, `MaterialPhysics.h`, `Propagator.h`, `RefitDriver.h`, `RefitLegAssembly.h`, `SurfaceLinearizationReference.h` | Narrow descriptor/state leaf operations and native refit. Generic core; retain. Any formula, covariance, material, leg, or ordering consolidation needs physics approval. |
| `NativeCylinderCylinderRefitDriver.h` | Pre-activation comparison driver with no production caller; tests and P1-era characterization still name it. Archive unique evidence, then delete or fold only genuinely unique test helpers into the live driver. |
| `Cell.h`, `SeedAnchor.h`, `LayerMask.h`, `LayerMask.h`, `SurfaceMeasurementIndex.h` | Device-safe CA value types and fixed-capacity metadata. Retain fixed capacity with runtime active prefixes. Do not remove live `SeedMetadataBase<N>` used by `CellSeed`. |
| `GenericTrack.h`, `GenericTrackOutputAdapter.h`, `TimeFrame.h`/`Tracker.cxx` MC-label sidecar | Generic result and adapter conversion utilities. `TimeFrame` owns labels parallel to generic tracks; `Tracker` computes them once, and both detector publication paths consume them. Keep `GenericTrack`; move detector-shaped output helpers outward if their callers are only workflows. |
| `ClusterSource.h`, `DecodedCluster.h`, `ClusterDecoding.h`, `ClusterDecoder.h`, `MultiSourceFrame.h`, `MultiSourceLoading.h`, `TimeFrameLoadFailure.h` | Generic normalized multi-source input and transactional decoding. Retain. Narrow includes where implementation-only decoder dependencies leak into public headers. |
| `ClockTimingPublicationView.h`, `SurfaceTiming.h`, `ROFTimingUniformity.h`, `ROFViews.h` | Runtime timing views and validation. Retain generic views; timing construction and publication clocks stay workflows/adapters. |
| `Configuration.h`, `TrackingConfigParam.h`, `ConfigKeyValuesPreflight.h`, `IndexTableConfiguration.h`, `IndexTableUtils.h` | Algorithm parameters/configuration and fixed-capacity index tables. Retain live runtime-prefix storage. Detector-named configuration registration is adapter compatibility, not a core routing authority; move only under a separately gated configuration migration. |
| `DetectorLayout.h`, `DetectorLayoutBuilder.h`, `DetectorLayoutSet.h`, `SparseTrackingTopology.h`, `SurfaceCatalogView.h`, `SurfaceDescriptor.h`, `IdTypes.h`, `StaticSurfaceDescriptor.h`, `ITSMFTDetectorDefinitions.h`, `SurfaceSpec.h` | Authoritative descriptors and graph data, but too many public ownership names. Consolidate `DetectorLayout` plus public sparse topology into `SurfaceGraph`/`SurfaceGraphView`; replace `DetectorLayoutSet` with a TimeFrame-owned `std::vector<SurfaceGraph>`; delete the stale configuration-key concept. `ITSMFTDetectorDefinitions` remains application data and should move outward without duplicating graph construction. |
| `ITSSurfaceSpec.h`, `MFTSurfaceSpec.h`, `ITSMFTDetectorDefinitions.h` | ITS/MFT application data. Valid compatibility owners, but not generic core concepts. Move to detector application include locations in a bounded include/API migration. |
| `SurfaceKinematicState.h`, `SurfaceKind.h` | Generic state representation. Retain; `SurfaceKind` must not become dispatch policy. |
| `SurfaceMeasurement.h` | Generic normalized measurement. Retain. |
| `SurfaceKinematicStateLegacyAdapters.h`, `SurfaceMeasurementAdapters.h`, `TrackingFrameInfoAdapters.h` | Compatibility conversions. The first two remain live at adapter/decoder leaves. `TrackingFrameInfoAdapters` and common `loadClusterTrackingFrameInfo` have no production caller found and are test-preserved deletion candidates. |
| `IOUtils.h` | Loader boundary plus decoder covariance/systematic helpers. The common `getClusterLayer`, iterator decoding overload, and `convertCompactClusters` API were deleted after a whole-repository audit found no common or downstream caller; detector-local ITS/MFT variants are unaffected. |
| `ITSSharedClusterCompatibility.h`, `MFTPublicationCompatibility.h`, `DetectorPublicationAdapter.h`, `DetectorTrackingOperationAdapterSupport.h`, `MFTFwdTrackHelpers.h` | Detector application/refit/publication compatibility. Responsibilities can remain, but ownership should move outside generic core headers as participant/interface consolidation proceeds. Do not replace with a central detector dispatcher. |
| Former common typed-MFT refit/export files | Typed MFT compatibility path with no production consumer; only `testMFTNormalizedRefit` and guards retained it. **Deleted by L0.** |
| `ParticipantId.h`, `TrackingParticipant.h`, `TrackingEngine.h` | Migration-era heterogeneous schedule contract around tracker instances. Retire after source-qualified graph partitions and schedule order become immutable TimeFrame data consumed by `Tracker`. |
| `SurfacePlanTrackingParticipant.h`, `TrackingInterface.h`, `TrackingOperationAdapter.h` | Application-composition cluster described in Sections 5–6. Delete participant/interface after the TimeFrame/Loader/Tracker composition is established; narrow the operation seam to architecture/refit work only. |
| `TimeFrame.h`, `SurfaceTrackingScratch.h`, `detail/SurfacePlanBinding.h`, `IOUtils.h` | One intended entity, currently split workspace, immutable graph partition, and atomic loader component. Fold scratch/partition ownership behind TimeFrame; keep Loader separate and acting on the entity. |
| `CATracker.h`, `TrackerTraits.h`, `Tracker.h` | Orchestrator plus CPU/GPU kernel-backend seam and one reversed forwarding relationship. Make `Tracker.h` canonical, delete the stale `CATracker` file name, and retain `TrackerTraits` distinctly. |
| `detail/TransitionPolicy.h`, `detail/TransitionPolicyBinding.h`, `detail/TransitionPolicyDispatch.h`, `detail/TransitionPolicyOperations.h`, `detail/TransitionPolicyState.h` | Private compatibility implementation for descriptor-family leaf selection. Keep private. Structural collapse is not safe if it changes operation order or arithmetic. |

## 8. Ranked cleanup inventory

### 8.1 Safe behavior-preserving cleanup

“Safe” means no intended behavior change, not that validation may be skipped.

| Rank | Exact files/action | Current callers | Replacement owner | Required gate | Deletion criterion |
|---:|---|---|---|---|---|
| 1 | Delete the three unused common typed-MFT refit/export files and their CMake entry. Migrate `testMFTNormalizedRefit.cxx` away from typed export. | No production caller; only `testMFTNormalizedRefit`, M7 guards, and stale CMake comments. | Existing generic `MFTFwdTrackHelpers`/`RefitDriver` result path; MFT workflow publication already consumes `GenericTrack`/adapter compatibility. | Focused native-refit and MFT publication tests, full serial `itsmft` suite, 43/43 fixture checks, 212/68 hashes, combined-leg and writer parity. | **Completed by L0:** repository search finds no common typed-MFT refit/export symbol or CMake entry; no production target loses a symbol. |
| 2 | Move the `Tracker` declaration/definition from `CATracker.{h,cxx}` to canonical `Tracker.{h,cxx}` and delete the forwarding/stale file names. | Common tracker users include either name; `TrackingInterface.h` is the forwarding header's production consumer found. | `Tracker.{h,cxx}`. | Build all common/ITS/MFT/combined targets; public-header dependency guard. | One canonical Tracker header/source and no `CATracker` implementation filename. |
| 3 | Delete scratch member `mPValphaX` and allocator/reset/swap bookkeeping. | No read or semantic write found; only initialization, clearing, allocator-match, and swap. | None. | Scratch allocator/swap/reset tests, sanitizer-capable focused tests, full suite and replay. | No field/reference remains and staged/live allocator matching still passes. |
| 4 | Use one scratch reset name; delete forwarding `resetScratch()` after changing direct callers to `reset()`. | Participant, loader wrapper, and migration-era tests. | `SurfaceTrackingScratch::reset()`. | Reset ordering, retry, dropped/structural failure tests. | One public scratch reset method and no semantic change in ordering. |
| 5 | Delete `loadITSAndMFT()` and `resetITSAndMFTEvent()` after migrating their tests to `loadEvent()` plus the current generic reset path. | No production caller found; fixed-source behavior is test-only. | Generic `loadEvent()`; reset ultimately belongs to `TimeFrame`. | Atomic load success/failure/retry and combined isolation tests. | No fixed ITS=0/MFT=1 loader API under common tracking. |
| 6 | Correct migration-era comments and test names that claim templated `Tracker`/`TrackerTraits`, M2/M6 temporary ownership, or an unfinished M7e seam. | Documentation and source comments/tests only. | Current responsibility wording. | Source guard plus documentation link/headings validation. | No stale architectural claim; behavioral assertions retained under intent-based names. |

The first cleanup slice is rank 1 only. Combining later ranks would make
review and regression localization worse.

### 8.2 Requires a bounded migration

| Rank | Exact files/action | Current callers | Replacement owner | Required gate | Deletion criterion |
|---:|---|---|---|---|---|
| 1 | Replace `MultiSourceTimeFrameLoader::{LoadTarget,LoadTargetImplSurface}` and participant `mLoadTarget` with direct atomic `load(TimeFrame&, inputs)` staging. | Combined workflow, participant, loader tests. | Loader stages a temporary TimeFrame-compatible event state and commits through TimeFrame's private/entity boundary. | Allocator identity, partial-stage failure, retry, dropped TF, combined source isolation, full replay. | No virtual scratch target, friendship-only public forwarding, or participant load-target member. |
| 2 | Consolidate `DetectorLayout`, `DetectorLayoutView`, `SparseTrackingTopology`, and `DetectorLayoutSet` into `SurfaceGraph`/`SurfaceGraphView` plus a TimeFrame-owned `std::vector<SurfaceGraph>`; delete `DetectorLayoutConfigurationKey` and always-true `rebuilt`. | Tracker/traits, standalone/combined construction, many graph fixtures. | TimeFrame owns each iteration graph once; sparse adjacency may remain private graph storage. | Structural graph parity, non-contiguous order, holes/seeding, CPU/device POD layout, full suite/replay. | No duplicate public layout/topology model, no retained fake key or plan wrapper, and combined owns one graph rather than two copies. |
| 3 | Replace the partial `adopt*`/`set*` construction sequence with one fallible TimeFrame entity-configuration boundary establishing graphs, parameters, partitions, pool, and workspace capacities in order. | Interface, participant, combined workflow, construction tests. | `TimeFrame` configure-once state; workflow still supplies adapter timing/publication inputs per event. | Allocator identity, construction failure, retry, destruction order, memory-limit, CPU/device backend tests. | TimeFrame cannot be loaded/executed partially configured; graph sizing never precedes allocator establishment. |
| 4 | Make graph partition/binding mandatory and per iteration; delete identity/source-0 fallback mapping and iteration-0 binding reuse. | Direct unit fixtures, standalone and combined construction. | `TimeFrame` stores one validated ordered `SurfaceGraphPartition` list per graph iteration. | Multi-iteration differing-hole/start-mask fixture, non-identity/sparse plan, invalid-partition failures, full replay. | No nullable binding, synthesized numeric traversal, or topology IDs borrowed from the wrong iteration. |
| 5 | Move fixed ROF tables, masks, and publication helpers out of participant/interface wrappers. | Standalone and combined workflows. | ITS/MFT workflow/application setup; TimeFrame/core continues borrowing runtime `ROFViews` only for the event. | Empty/first/last/diamond timing, mask, load-failure replacement, writer and sidecar parity. | TimeFrame/Tracker/traits own no fixed detector table or detector publication sidecar. |
| 6 | Make `Tracker::run(TimeFrame&, TrackerTraits&)` orchestrate TimeFrame's explicit ordered graph partitions and invoke its all-or-nothing generic reset; delete `TrackingEngine.{h,cxx}`, `TrackingParticipant.h`, and `ParticipantId.h`. | Combined workflow and engine/participant tests. | Tracker is a component; TimeFrame owns state and partition order; no dynamic coordinator wrapper. | Exact partition order, source isolation, success/recoverable/structural/exception reset-count tests, combined replay. | No engine/participant class and no equivalent renamed schedule executor; Tracker owns no frame state. |
| 7 | Migrate standalone and combined wrappers to Loader + TimeFrame + Tracker; delete `TrackingInterface.{h,cxx}`, `SurfacePlanTrackingParticipant.{h,cxx}`, aliases, mixins, and float sentinel. | ITS/MFT/combined CA workflow tasks and wrapper-heavy tests. | TimeFrame owns generic state; Loader/Tracker act on it; workflow owns raw ROFs/clocks/publication. | Standalone lifecycle/config/output tests, dropped-TF behavior, combined parity, exact writer/replay baseline. | One entity/component composition path and no interface/participant wrapper remains. |
| 8 | Narrow `TrackingOperationAdapter`: move `completeAccepted()` and `resetAdapterState()` to application publication handling; retain only the exact refit operation if still needed. | `TrackerTraits`, interface, participant. | TimeFrame owns generic accepted state; application adapter consumes only successful ordered results. | Candidate ordering, holes, all refit rejection paths, shared/pattern compatibility, failure-stale-state and replay gates. | Core operation seam neither stages publication nor owns lifecycle; delete seam entirely if one generic refit implementation remains. |
| 9 | Make the `TrackerTraits` architecture contract explicit and port tests away from assuming CPU is the only backend; retain virtual kernel entry points. | Tracker, direct kernel fixtures, future common GPU integration. | Tracker borrows CPU or GPU traits; traits implement backend kernels and act on TimeFrame/backend-local views without owning the entity. | CPU replay plus a real pinned CUDA/HIP build when a common GPU traits backend exists. | Backend substitution changes no TimeFrame/Tracker API and introduces no detector/layer specialization. |
| 10 | Remove `TimeFrame`'s no-op virtual device-propagator hook/member if full build configurations confirm no common-derived device frame. | A migration-era ownership test; frozen ITS has a separate live override in a different class. | None in common CPU runtime; common GPU traits must declare its actual device state ownership. | Whole-repository inheritance/API audit and actual device build when pinned toolchain exists. | Common `TimeFrame` contains no placeholder device pointer; GPU ownership is explicit. |
| 11 | Delete common `TrackingFrameInfoAdapters`, `loadClusterTrackingFrameInfo`, and obsolete compact-cluster conversion APIs. | Tests only for the former; a whole-repository audit found no common or downstream caller for `getClusterLayer`, the iterator decoding overload, or `convertCompactClusters`. | `SurfaceMeasurement` decoder path for production; detector-local frozen utilities remain untouched. | Whole-repository symbol/header audit, downstream build, covariance/systematics tests. | No downstream consumer and no lost compatibility contract. |
| 12 | Retire `NativeCylinderCylinderRefitDriver.h` after preserving the comparison evidence it uniquely supplies. | Dedicated test, guard, and historical design/P1 harness. | Live `RefitDriver`; durable validation docs/artifacts for historical evidence. | Demonstrate all production numerical/refit properties remain covered; no replay delta. | No production or active validation need for the unwired driver. |

### 8.3 Requires explicit physics or algorithm approval

| Candidate | Why it is not safe cleanup | Decision evidence required |
|---|---|---|
| Collapse cylinder/disk private transition-operation templates into one formula or branch. | It can change projection, rotation, reference coordinate, evaluation order, or policy grouping. Similar orchestration does not imply identical mathematics. | Focused descriptor-pair numerical study, cut-flow comparison, covariance/parameter deltas, and physics approval. |
| Change `Propagator`, native refit leg assembly, covariance sanitization, material treatment, MinPt/chi2 gates, holes, road ordering, or candidate sorting while merging classes. | These are the accepted candidate-physics path, not architecture residue. | Separate physics campaign and explicit approval; never bundle with cleanup. |
| Shrink fixed `TrackSeed`, `CellSeed`, surface-mask, or index-table capacities, or replace them with dynamic device values. | Changes device ABI, capacity/failure behavior, and potentially candidate populations. | CPU/GPU capacity characterization and algorithm/device sign-off. |
| Remove private transition-policy grouping because descriptors appear sufficient. | Grouping currently fixes schedule and hot-loop operation selection; removing it can alter order and floating-point behavior. | Order proof, performance/device evidence, replay and physics approval. |
| Generalize arbitrary mixed cylinder/disk plans. | Current supported descriptor combinations encode an algorithm validity boundary. | New algorithm specification and detector validation; outside cleanup. |

### 8.4 Frozen legacy and out of scope

| Area | Treatment |
|---|---|
| `Detectors/ITSMFT/ITS/tracking/` and its GPU `TimeFrame`/tracker utilities | Frozen legacy reference. Similar names such as `mPValphaX`, `convertCompactClusters`, and `setDevicePropagator` do not make common dead-code deletion applicable there. |
| Legacy ITS workflows | No modification. They remain a validation/reference application. |
| Raw ROF ownership, DPL defaults, writers, sidecars, typed persisted outputs | Workflow/application-owned and behavior-frozen for cleanup. Ownership may move outward, content may not change. |
| `SeedMetadataBase<N>` used by `CellSeed` | Live fixed device value support; retain. |
| `TimeFrame::mNTotalLowPtVertices` and other support state without a whole-repository liveness proof | Retain until a dedicated audit. |
| Triplet R&D, mixed-detector validity, and ALICE3 integration | Future work. Cleanup should lower their entry cost but not implement them. |

## 9. Tests and dependency cost

Useful behavior coverage to retain includes sparse/non-identity plan order,
source-qualified measurements, runtime counts, ROF timing/masks, atomic load
and retry, non-templated core, refit/covariance behavior, publication
sidecars, failure reset, and standalone/combined output parity.

Several tests instead pin temporary type shape or names:

- `testM6dMFTMigration`, `testM6e1StandaloneMFTMigration`, and
  `testM6e2ITSWorkspaceMigration` should be renamed or folded into current
  adapter/lifecycle contracts after `TrackingInterface` retirement;
- `testGate4B31OwnershipContract` explicitly preserves the inert common
  `TimeFrame::setDevicePropagator()` hook and should be narrowed to live
  ownership assertions before that hook is deleted;
- `testMultiSourceTimeFrameLoader` preserves fixed `loadITSAndMFT()` and
  `resetITSAndMFTEvent()` exemptions; migrate the same atomicity assertions to
  generic bindings rather than weakening coverage;
- `testM7eAdapterBoundaryGuard` now proves the dead typed-MFT compatibility
  path is absent;
- `testNativeCylinderCylinderRefitDriver` is a historical comparison harness,
  not production-path coverage; preserve its evidence before retirement.

The final M7 runtime-core guard remains valuable. Cleanup should tighten it
around current responsibilities rather than accumulate named exemptions for
deleted migration bridges.

Public include cost can be reduced without a new pimpl/service layer:

- make `Tracker.h` canonical and delete the stale `CATracker` file name;
- stop including `TrackingInterface.h` from the participant merely for
  sidecar-owner mixins;
- move detector-specific specs/publication helpers out of generic public
  includes as their callers migrate;
- keep private transition-policy headers private to implementation-facing
  headers; and
- keep `TrackerTraits` public only to the degree needed by CPU/GPU backend
  substitution; move non-kernel compatibility details out of its header.

## 10. Explicit class recommendation

| Class | Recommendation | Reason |
|---|---|---|
| `TrackingEngine` | **Delete after schedule/result sequencing moves into `Tracker` and reset mechanics into `TimeFrame`.** | It is a stateless loop around participant-owned trackers, not an independent kernel or data owner. Retaining it creates two orchestrators. |
| `TimeFrame` | **Become the sole generic tracking entity.** | It owns configuration, allocator-backed workspace, normalized event input, and generic results. Graphs, parameters, partitions, and scratch state should not be spread across wrappers or Tracker. |
| `Tracker` | **Remain as a non-owning generic tracking component.** | It sequences iterations and CA stages by reading TimeFrame's graphs/partitions and invoking a borrowed traits backend. It retains no graph, parameter, workspace, pool, event, publication, or detector state. |
| `TrackerTraits` | **Remain distinct as the architecture kernel backend.** | CPU and GPU implementations need the same orchestration with different kernels/storage. The virtual kernel seam is intentional and mirrors the live ITS CPU/GPU architecture. |
| `TrackingInterface` | **Delete after entity/component ownership is established.** | Generic data/configuration belongs in TimeFrame, loading in Loader, execution in Tracker, and raw timing/publication in DPL/application adapters. No replacement interface facade is justified. |
| `SurfacePlanTrackingParticipant` | **Delete with `TrackingParticipant`.** | Its useful graph-partition/workspace content moves into TimeFrame. Its timing/publication/load-target composition belongs elsewhere. |
| `TrackingOperationAdapter` | **Narrow to refit-only, then reassess deletion.** | Refit is an algorithm operation; publication completion and reset are application lifecycle. Do not replace it with a callback framework. |

## 11. Ordered cleanup slices

Each slice deletes its bridge before the next one starts; no long-lived dual
implementation is proposed.

1. **C0 — dead typed-MFT refit/export retirement.** Delete the three rank-1
   files and migrate the typed fixture to generic native-refit assertions.
2. **C1 — zero-risk residue.** Canonicalize `Tracker.{h,cxx}`, then delete
   the `CATracker` file name; delete `mPValphaX` and the duplicate scratch
   reset name in separate reviewable commits.
3. **C2 — loader simplification.** Replace the one-implementation virtual
   load target with direct `load(TimeFrame&, inputs)` staging; remove fixed
   ITS+MFT wrappers and public scratch-target plumbing.
4. **C3 — surface-graph consolidation.** Introduce the `SurfaceGraph` name by
   replacing, not wrapping, `DetectorLayout`/public sparse-topology ownership;
   replace `DetectorLayoutSet` with a TimeFrame-owned graph vector; delete the
   stale key, plan wrapper, and combined duplicate graph copies.
5. **C4 — one TimeFrame entity.** Move parameters, per-iteration partitions,
   allocator, and scratch/workspace capacities behind one fallible TimeFrame
   configuration boundary; eliminate partial setter ordering and public
   scratch ownership.
6. **C5 — binding correctness.** Make partition data mandatory and
   per-iteration; remove identity/source-0 and iteration-0 reuse fallbacks.
7. **C6 — application ownership.** Move ROF fixed tables and detector
   publication compatibility out of the participant/interface, preserving
   runtime `ROFViews` in core.
8. **C7 — one tracking component.** Make Tracker act on `TimeFrame&` through
   `TrackerTraits&`, reading frame-owned partition order and invoking the
   frame's generic reset; delete `TrackingEngine` and `TrackingParticipant`
   without adding a renamed executor.
9. **C8 — wrapper retirement.** Migrate standalone/combined workflows to
   Loader + TimeFrame + Tracker; delete `TrackingInterface`,
   `SurfacePlanTrackingParticipant`, and the float sentinel.
10. **C9 — operation seam narrowing.** Move accepted-result compatibility and
   reset outward; retain only a necessary refit operation.
11. **C10 — backend contract.** Keep CPU/GPU TrackerTraits substitution
   explicit while moving non-kernel compatibility state out of the backend
   seam; run a real device build only when the common backend exists.
12. **C11 — public compatibility tail.** Resolve dead common IO adapters,
    inert device hooks, old comparison driver, and detector-spec header
    placement one item at a time.

No slice above authorizes a physics simplification. Section 8.3 candidates
must remain separately gated even if a class merge makes them tempting.

## 12. First-slice implementation contract

The next task can be bounded exactly as follows:

**Input scope**

- the three deleted common typed-MFT refit/export files
- common tracking CMake source list and stale comments
- `testMFTNormalizedRefit.cxx`
- M7 boundary/guard tests that currently whitelist the typed path

**Required result**

- no production or test target includes or links the typed MFT refit/export;
- meaningful refit tests assert the existing generic
  `SurfaceKinematicState`/`TrackingCandidate` result and unchanged failure
  behavior rather than recreating a typed MFT output record;
- detector publication/output code is unchanged;
- repository guards change from “typed path isolated” to “typed path absent.”

**Validation**

- build common tracking, focused refit tests, ITS/MFT/combined workflows;
- run the complete serial `itsmft` suite;
- verify 43/43 fixture checksums before and after;
- reproduce ITS 212 / `46913a67a7e2fe7462e29df0db264fa8` and
  MFT 68 / `8106b08571ca593c6b76ff72b761a680`;
- require combined legs and initialized writer content to match the accepted
  P1/M7f candidate, retaining only the known undefined
  `MFTTrack.mInvQPtSeed` byte exclusion.

**Exit criterion**

A repository search finds no common typed-MFT refit/export symbol or typed
refit overload, and no replacement compatibility wrapper has been introduced.
This deletion is independent of reset, participant, interface, or physics
changes. L0 satisfies this criterion.

## 13. Decision record

Measured structure:

- there is one production CA body;
- there are two application compositions;
- reset policy is duplicated across tracker and event owners;
- the load target is one-implementation type erasure;
- common `TrackerTraits` currently has one CPU implementation, but its virtual
  kernel API intentionally mirrors the live frozen ITS CPU/GPU backend seam;
- combined tracking copies one authoritative graph into two layout-set
  owners;
- a single binding built from iteration 0 is reused with a plan that can own
  different per-iteration graphs;
- `DetectorLayoutConfigurationKey` is not used as a key and retains build
  inputs primarily to recover ordered surfaces;
- standalone and combined initialize allocator/plan/workspace in different
  orders, with combined sizing scratch before rebinding its pool;
- current generic state is split between `TimeFrame`, participant-owned
  scratch, layout-set owners, bindings, and Tracker-owned parameter/pool
  handles;
- typed MFT refit/export files have no production caller.

Inference:

- deleting dead typed compatibility first is lower risk than beginning with
  lifecycle consolidation;
- preserving `Tracker` as orchestrator and `TrackerTraits` as CPU/GPU kernel
  backend is the intentional architecture;
- making `TimeFrame` the sole reusable tracking entity removes split state
  ownership; Loader and Tracker can remain non-owning components acting on
  it;
- retiring engine, participant, and interface wrappers removes more
  complexity than narrowing them, provided TimeFrame stores explicit
  detector-neutral graph partitions and workflow publication state stays
  outside;
- `SurfaceGraph` is more accurate than `DetectorGraph` because graph nodes
  are surfaces and one graph may span multiple detector identities;
- consolidating layout/topology names and configure-once initialization is a
  bounded structural migration, not safe textual cleanup.

No physics inference or acceptance decision is made by this audit.

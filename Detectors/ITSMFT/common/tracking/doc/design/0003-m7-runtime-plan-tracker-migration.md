# Design note 0003: runtime-plan Tracker/TrackerTraits migration

Status: M7a design and M7b/M7c/M7d implementation records; M7d validation
recorded in [design note 0006](0006-m7d-nontemplated-tracker-core.md)
Date: 2026-08-05
Scope: `Detectors/ITSMFT/common/tracking`
Companion plan: [GenericTrackingEngineMigration.md](../GenericTrackingEngineMigration.md)
Architecture decision: [ADR 0007](../decisions/0007-generic-tracking-engine-boundary.md)
Related designs: [descriptor-driven operations](0001-descriptor-driven-operation-boundary.md),
[M6 workspace migration](0002-m6-generic-workspace-migration.md)

This note is the Gate 4 M7a design slice. It records the remaining
`NLayers` dependencies, assigns ownership, and makes the next production
slices bounded. It deliberately does not change C++ or begin
`Tracker`/`TrackerTraits` de-templating.

## 1. Decision summary

The generic tracker will not receive a new `RuntimeTrackingPlan` wrapper.
The runtime plan is the existing composition of four concrete objects:

| Role | Existing owner | Runtime information it is authoritative for |
|---|---|---|
| Static surface and sparse topology plan | `DetectorLayoutSet` and its `DetectorLayoutView` | ordered `SurfaceDescriptor`s, active-surface count, transitions, cells, and the plan's sparse topology |
| Application-to-plan binding | `SurfacePlanBinding` | the ordered owned-surface positions, source-qualified binding, and global transition/cell to compact-scratch slots |
| Mutable plan-sized event workspace | `SurfaceTrackingScratch` | per-owned-surface data, transition/cell storage, index tables, and event-local CA working state; its `getNOwnedSurfaces()`, `getNTransitions()`, and `getNCells()` are capacity facts |
| Shared event result and normalized input | `TimeFrame` | normalized measurements, vertices, `CommonTrack`s, and `TrackClusterReference`s |

The ordered span supplied to `SurfacePlanBinding::build()` is retained by the
`DetectorLayoutSet` configuration key and exposed by the application
participant; `SurfacePlanBinding` validates that span and owns the global-ID
to positional/scratch-slot mapping. That validated ordered span is the
loop-order source. The `DetectorLayoutSet` iteration's `activeCount` and the
ordered-surface span are the active count. The scratch counts and the sparse
topology counts are the transition/cell capacity sources.
`TrackingParameters::NLayers` is only an adapter-configuration compatibility
field during the migration; it must not be used as the generic core's loop or
capacity authority.

`TrackSeed::MaxSurfaces == MaxLayoutSurfaces` remains the sole fixed whole-track
seed capacity. Its `SurfaceMask` says which positional slots are active. A
runtime plan may use fewer positions, but it may not make the device value
type's array shorter at runtime. `CellSeed` remains a three-cluster cell value
because three measurements are the CA cell definition, not a detector layer
count. `SeedMetadataBase<N>` therefore remains live for `CellSeed` only.

The first architectural proof will be deletion, not renaming: the common
production declarations and explicit instantiations of `Tracker<NLayers>` and
`TrackerTraits<NLayers>`, together with the `CATracker` alias, will disappear.
An application adapter may remain templated temporarily, but it will contain
non-templated `Tracker`/`TrackerTraits` objects. There will be no
`Tracker<7>` or `Tracker<10>` hidden inside a new wrapper.

No new ADR is needed for M7a. This note operationalizes ADR 0007 decisions 2,
6, 7, 8, and 10; it does not change them.

## 2. Boundary and ownership invariants

The ownership model after M6g is the starting point:

```text
combined DPL task
  owns event lifecycle, raw ROFs, timing/publication context, and schedule
    -> MultiSourceTimeFrameLoader (atomic load)
    -> TrackingEngine (ordered track/reset over participants)
      -> ITS application participant
      -> MFT application participant
    -> CommonTrack/output adapters and publication
```

The DPL task remains the event owner. In particular, it creates, invalidates,
and resets publication clocks and sidecars at event boundaries and on every
load/tracking failure path. `TrackingEngine` resets participant-local state and
the shared `TimeFrame` according to its existing contract; it is not a
combined ITS+MFT coordinator and does not own publication timing.

The following remain hard boundaries for every M7 implementation slice:

- `TrackingEngine`, `TrackingParticipant`, `TimeFrame`,
  `SurfaceTrackingScratch`, `SurfacePlanBinding`, `TrackSeed`, sparse topology,
  and surface operations remain detector-neutral.
- No generic core header or source gains ITS/MFT identity, source 0/1
  conventions, DPL, workflow, writer, or output-track dependencies.
- `SurfaceKind` describes the representation carried by a
  `SurfaceKinematicState`. It is not a schedule, topology, detector, or
  operation-selection policy.
- `TransitionPolicyTag` remains a contained private compatibility detail. M7
  does not expose it, store a replacement `SurfaceKindPair`, or create a
  public family-dispatch taxonomy.
- Raw ROF input and the output sidecars remain workflow/adapter concerns.
  The dead common `loadROFrameData()`/`resetROFrameData()`/
  `prepareROFrameData()` family stays deleted; live frozen ITS raw-ROF code
  and MFT raw-ROF workflow ownership are not part of this migration.
- The native `Propagator` path, covariance policy, CA choices, hole behavior,
  ordering, output format, and workflow defaults remain unchanged unless a
  later slice obtains a separate physics decision.

## 3. Current `NLayers` audit

The audit covered the common tracking include/src/test tree, the ITS/MFT
adapter seams, refit helpers, explicit instantiations, and the frozen ITS
GPU-facing hierarchy. The classifications below use the requested numbers:

1. compile-time/device-facing for now;
2. runtime count from the plan/workspace;
3. adapter compatibility outside the core;
4. dead or redundant, to delete in a later slice;
5. separately gated physics/algorithm decision.

### 3.1 Core orchestration and plan storage

| Current dependency | Evidence and actual responsibility | Classification | M7 disposition |
|---|---|---:|---|
| `Tracker<NLayers>` in `CATracker.h/.cxx` | Owns the iteration loop and forwards to traits; its only explicit production instantiations are 7 and 10. The `CATracker` alias and `TrackerITS`/`TrackerMFT` aliases are naming/ABI residue. | 2 for orchestration; 4 for aliases | Make `Tracker` a single non-templated orchestrator taking a non-templated traits object. Delete `CATracker` and the detector aliases in the same slice. |
| `TrackerTraits<NLayers>` in `TrackerTraits.h/.cxx` | Contains the shared tracklet, cell, neighbour, road, acceptance, and marking flow, but still uses `NLayers` for arrays, topology views, refit signatures, explicit detector branches, and 7/10 instantiations. | 2 for loop/capacity data; 3 for current output/refit hooks; 4 for compile-time aliases/instantiations | Convert the body once. Read active positions from the plan/binding and workspace; pass operation-local data to adapter-owned refit/output code. Do not make a second runtime and templated implementation. |
| `TrackerTraits` `mLayerMaterial`, `mLayerMeasurements`, `candidateReachableLayers`, and similar `std::array<T, NLayers>` temporaries | These arrays are indexed in the current legacy ordered-surface space. Their contents are per-plan data, not device ABI. | 2 | Use a runtime span/vector sized to the binding's ordered surfaces. Keep a fixed `MaxLayoutSurfaces` scratch only where a device value type requires it. |
| `TrackingParameters::NLayers`, `NeighboursPerRoad()`, `CellsPerRoad()`, and `TrackletsPerRoad()` | Configuration and formulas still encode the old layer-count vocabulary. The current tracker validates it and uses it as a hot-loop bound in places. | 3 now; 4 after adapter migration | Validate the configured value against `DetectorLayoutIterationConfiguration::activeCount` at the application edge. Core algorithms derive road/cell counts from the runtime schedule/topology. Remove the field from the core-facing configuration only after all adapters consume the plan. |
| `stateFamilyFromNLayers<NLayers>()` and `if constexpr` family branches in traits | These are hidden 7/10 selectors. The state family itself is representation metadata. | 3 at adapter edges; 4 in core | Resolve representation from the actual `SurfaceDescriptor` and the state's `SurfaceKind`. Do not use the mapping as a new central dispatch type. |
| `SurfaceTrackingScratch::adoptPlan()` and runtime vector storage | Already owns runtime-sized surface/transition/cell containers. It still exposes templated dispatchers for dual legacy ROF/topology objects. | 2 for the runtime storage; 3 for legacy ROF objects | Retain the existing object and counts. Remove only the templated compatibility accessors after their application-owned views exist. |
| `TrackingTopology<NLayers>` | Layer-indexed topology with `N`-derived `MaxTransitions`/`MaxCells`, dual-stored in `SurfaceTrackingScratch`; its view is device-readable but represents the old layer graph. | 4 for the common core; 1 only for the transitional device view | Make `SparseTrackingTopologyView` the plan topology. If a device copy needs fixed storage, size it from `MaxLayoutSurfaces`/`MaxLayoutTransitions`/`MaxLayoutCellTopologies`, never from ITS=7 or MFT=10. Delete the layer-topology storage once all callers consume sparse IDs. |

The important distinction is between a runtime capacity and a fixed storage
limit. A plan count is data. A fixed maximum used to make a host/device value
copyable is not a detector choice.

### 3.2 Index tables and ROF/timing helpers

| Current dependency | Evidence and actual responsibility | Classification | M7 disposition |
|---|---|---:|---|
| `IndexTableUtilsCore`'s `std::array<float, MaxLayoutSurfaces>` extent storage | `getColBinIndex()` and related functions are `GPUhdi()`. The current fixed maximum is needed for device-portable per-surface lookup data; the populated prefix is runtime-sized. | 1 for fixed capacity; 2 for active prefix | Keep the fixed maximum and pass a runtime active-surface span/count. It must not regain an `NLayers` template. |
| `template<int> using IndexTableUtils = IndexTableUtilsCore` | The alias is now a source-compatibility mask around a non-templated type. It does not provide a real capacity guarantee. | 4 | Delete the alias after all call sites use `IndexTableUtilsCore`; this is a safe structural deletion, not a physics change. |
| `bindIndexTableConfiguration<Tag, NLayers>()` and `indexTableConfigurationsMatch<NLayers>()` | Validation checks `params.NLayers`, fills `NLayers` extents, loops over `NLayers`, and has explicit 7/10/tag instantiations. The tag is currently private implementation dispatch. | 2 for active extents; 3 for parameter decoding; private tag remains temporarily | Bind a runtime span of surface extents and an active count from the plan. Keep tag containment unchanged; do not replace it with a public pair taxonomy. |
| `ROFOverlapTable<NLayers>`, `ROFVertexLookupTable<NLayers>`, `ROFMaskTable<NLayers>` in `SurfaceTrackingScratch` | These are dual-typed frozen `ITStracking` timing/mask structures selected by `if constexpr` and used through templated scratch helpers. They carry per-ROF timing, not detector-neutral CA topology. | 3 | The application/workflow owns timing construction and supplies a runtime view to the core. The DPL task remains responsible for raw ROFs, validity, and reset. Do not redesign timing semantics in M7. |
| `SurfacePlanTrackingParticipant<NLayers>::configureRofTables()` | Constructs the three templated ROF tables and enables all ROFs for one application leg. This is adapter setup, not generic tracking orchestration. | 3 | Keep the application setup temporarily at the edge; remove its dependence on `Tracker<NLayers>` once the core accepts runtime views. No new common coordinator is allowed. |
| `TrackingInterface<NLayers>` ROF configuration and explicit 7/10 instantiations | This facade binds detector identity, source decoding, ROF timing, beam configuration, sidecars, and the public single-detector workflow. | 3 | It may remain a detector-specific application facade. Its template is not permission for the generic tracker to remain templated. Any future interface cleanup is separately scoped from core `Tracker`/`TrackerTraits`. |

### 3.3 Seeds, cells, and output representation

| Current dependency | Evidence and actual responsibility | Classification | M7 disposition |
|---|---|---:|---|
| `CellSeed` / `SeedMetadataBase<ClustersPerCell>` | A cell always carries three cluster slots and two tracklet references. The base is still live and is not a whole-track compatibility bridge. | 1 | Retain `SeedMetadataBase<N>` for `CellSeed`; do not broaden M6f's deletion. Its 16-bit `LayerMask` is a cell-local/legacy-pattern field, not the generic plan mask. |
| `TrackSeed` | Non-templated, `GPUhd()`, trivially copyable, with `MaxSurfaces == MaxLayoutSurfaces` cluster slots and a 32-bit `SurfaceMask`. Active positions come from the mask and plan order. | 1 for fixed device capacity; 2 for active count | Retain exactly this representation. Do not reintroduce `TrackSeedTpl` or make the array heap-backed. |
| `CATrackTypeHelper<NLayers>`/`CATrackType` and `bounded_vector<CATrackType<NLayers>>` | Resolves to `o2::its::TrackITSExt`, holds typed acceptance/output staging, and is consumed by `DetectorTraits` and output-side compatibility. It is not a generic tracking result. | 3 | Move typed output construction and shared-cluster/pattern compatibility to ITS/MFT adapters. Core acceptance publishes `CommonTrack` and references. Delete the alias and the generic typed vector after adapter parity is proven. |
| `LayerMask` views on `TrackSeed` and `CommonTrack` | `LayerMask` is 16-bit and cannot represent all 32 plan positions. `TrackSeed` already keeps `SurfaceMask` authoritative and exposes a compatibility view for old hole/pattern code. | 1 for compatibility value types; 4 for core use as a plan mask | Core schedule, topology, and whole-track hit selection use `SurfaceMask`/`SurfaceId`. Keep a `LayerMask` view only at a legacy output/algorithm edge until those consumers are migrated. |

The old ITS `TrackSeed<NLayers>` and related `ITStracking` track/cell types are
outside the common-CA production tree and belong to frozen legacy ITS
workflows. They are not evidence that common `TrackSeed` is incomplete and
are not touched by M7.

### 3.4 Refit and operation dispatch

| Current dependency | Evidence and actual responsibility | Classification | M7 disposition |
|---|---|---:|---|
| `LayerMeasurementSpans<NLayers>` | `std::array<gsl::span<const SurfaceMeasurement>, NLayers>` is used by `TrackerTraits`, `DetectorTraits`, `RefitLegAssembly`, and `NativeRefitDriver`. It is host-only operation input, not a device ABI. | 2 | Replace it with a span-of-spans whose length is the binding's ordered-surface count. Preserve traversal order and holes exactly. |
| `fitTrackSeedLegs<NLayers>()`, `assembleRefitLegSlots<NLayers>()`, and `nativeRefitTrackCylinderCylinder<NLayers>()` | The refit walkers use `NLayers` for leg endpoints, stack buffers, and the MinPt slot. The active algorithm is already the native `Propagator` path. | 2 for traversal and buffer bounds; 5 for any formula/acceptance change | Make leg endpoints and measurement slots runtime plan positions. Preserve the existing three-leg sequence, acceptance gates, hole handling, covariance sanitization, and `Propagator` calls. Do not unify formulas or alter MinPt semantics in M7. |
| `DetectorTraits<NLayers>::refitSeed()` | Selects ITS/MFT output type, refit helper, pattern copy, transient pattern clearing, polarity, and beam configuration through `if constexpr`. | 3 | Move this set of hooks to the concrete application adapters. The generic core receives/produces `SurfaceKinematicState`, `TrackSeed`, and `CommonTrack`; output adapters retain `TrackITSExt`/MFT-specific representation. |
| `TransitionPolicyTag` templates in private operations | The current tracklet/cell/road leaves use the tag as a contained implementation selector. It is not an API or stored plan property. | 1/3 temporary implementation detail | Leave containment intact through the runtime-plan migration. At the operation boundary inspect actual endpoint descriptors and state representation; delete the tag only in a separate, replay-gated cleanup when one implementation is available. |
| `passesCellRoadPrecut<DiskDisk>`, legacy MFT reference-Z lookup, and other special cuts | These preserve known MFT/legacy behavior. They are not evidence that orchestration needs two trackers. | 5 | Keep unchanged. A descriptor-driven replacement or removal requires its own physics/algorithm sign-off and A/B evidence. |

The existing [descriptor-driven operation design](0001-descriptor-driven-operation-boundary.md)
already identifies the narrow leaf differences. M7 only changes how the
shared orchestration supplies runtime plan positions and chooses those leaves;
it does not change a leaf's numerical behavior.

### 3.5 Adapter seams, explicit instantiations, tests, and GPU scope

| Area | Remaining dependency | Classification | Ownership result |
|---|---|---:|---|
| `SurfacePlanTrackingParticipant<NLayers>` | Its static assertion accepts only 7/10; it owns a leg's scratch, binding, parameters, and compatibility sidecar, and currently embeds `TrackerTraits<NLayers>`/`Tracker<NLayers>`. | 3 | It is a temporary application participant seam, not a coordinator. It may remain templated only until it embeds the non-templated core and its detector-specific sidecar setup moves to the adapter edge. It must never own the combined schedule, clocks, or event reset sequencing. |
| `TrackingLoadPolicy<DetId, NLayers>` and `TrackingLoadPolicyN` | Beam-position and source/configuration setup uses detector identity and compile-time counts. | 3 | Keep at the ITS/MFT application boundary; remove the `N` alias when its callers are direct. The core sees an already-configured `TimeFrame`. |
| Explicit instantiations | Common production currently instantiates `Tracker<7/10>`, `TrackerTraits<7/10>`, `DetectorTraits<7/10>`, `TrackingInterface<7/10>`, `SurfacePlanTrackingParticipant<7/10>`, scratch helper families, and index-table binders for 7/10 and both private tags. | 3 while callers migrate; 4 after migration | Treat each as a deletion checklist, not an ABI to preserve. The common core's M7 exit guard must show no `Tracker<7>`/`Tracker<10>` or `TrackerTraits<...>` instantiation. Adapter-only instantiations may remain until their own bridge exit criteria are met. |
| Common tests and fixtures | Orchestration, load-failure, ROF, refit, topology, and adapter tests use `template<int NLayers>` fixtures and explicit 7/10 cases. | 3 in tests only | Keep parameterized fixtures while they cover both plans, but migrate assertions to plan/binding counts and `TrackSeed`/`SurfaceMask`. Tests must not preserve a production compatibility type merely because a fixture uses it. |
| Common GPU/device-facing types | `SurfaceId`, `SurfaceMask`, `SurfaceDescriptor`, `SurfaceKinematicState`, `SurfaceMeasurement`, `SparseTrackingTopologyView`, `CellSeed`, `TrackSeed`, and `IndexTableUtilsCore` are `GPUhdi()`/trivially-copyable value/view types. | 1 for fixed representation and ABI | Keep fixed widths and `MaxLayoutSurfaces` capacity. Runtime counts are scalar/view fields. No GPU validation is claimed for M7a; a device build is required only in a later implementation slice when the toolchain exists. |
| Frozen ITS GPU hierarchy | `ITStracking` retains `TimeFrame<NLayers>`, `TrackingTopology<NLayers>`, `TrackITSInternal<NLayers>`, follower/refit types, and CUDA kernels with fixed ITS layer arrays. | 3, frozen workflow exclusion | Exclude narrowly: this is not common-CA production and frozen ITS workflows must not change. The M7 guard scans common tracking core paths and records this exclusion explicitly. |

## 4. Runtime plan API and ownership target

The non-templated core will use these existing values directly:

```text
DetectorLayoutView layout = layoutSet.getLayoutView(iteration)
orderedSurfaces = layoutSet.getConfigurationKey().orderedSurfaces
                 (validated by SurfacePlanBinding)
activeCount = orderedSurfaces.size()
transitionView = layout.topology / binding transition span
cellView = layout.topology / binding cell span
scratch capacities = scratch.getNOwnedSurfaces(), getNTransitions(), getNCells()
```

There is no new owner object around this composition. The following rules make
the source of every bound explicit:

1. A loop over surfaces walks the binding's ordered surface span. It never
   loops to `TrackingParameters::NLayers` or a detector constant.
2. A loop over transitions/cells walks the binding's compact span or the
   `SparseTrackingTopologyView` count. It never assumes a dense local
   `NLayers` topology.
3. A per-surface temporary is a runtime span or a fixed
   `MaxLayoutSurfaces` device buffer with an explicit populated count.
4. A seed's active count is `TrackSeed::getSurfaceMask().count()`; its slot
   interpretation is the adopted ordered plan, not the numeric `SurfaceId`.
5. A road/cell/hole decision reads the sparse topology and operation-local
   policy parameters. `StartLayerMask` and `HoleLayerMask` remain compatibility
   inputs until their plan-position migration is complete; they do not become
   a new generic taxonomy.
6. A `TrackingParameters::NLayers` mismatch with the plan is a structural
   adapter error, not permission to resize or reinterpret the core's plan.

This model handles a plan with any valid number and ordering of surfaces within
the established `MaxLayoutSurfaces` limit. ITS=7 and MFT=10 become two
application plan instances, not two core class instantiations.

## 5. One orchestration for cylinders and disks

The shared tracklet/cell/road/refit loops consume a `SurfaceTransition` and
resolve its `from`/`to` IDs through the current `DetectorLayoutView`. At the
narrow operation boundary, the operation receives the actual
`SurfaceDescriptor`s and the current `SurfaceKinematicState`. Cylinder/disk
differences are confined to the already identified operations: projection
coordinate, propagation target, local rotation, material conversion, and
representation conversion.

The implementation shape is:

```text
one runtime-plan loop
  -> resolve from/to SurfaceDescriptor
  -> one operation-local projection/propagation/update/refit call
  -> continue with the same topology, candidate, hole, and acceptance flow
```

This is not a new central dispatch taxonomy. In particular:

- no `DetectorId` switch selects the core algorithm;
- no `SurfaceKind` value selects scheduling or topology;
- no persisted `SurfaceKindPair`/equivalent tag is introduced;
- no second tracker or type-erased operation bridge is introduced;
- the current private `TransitionPolicyTag` machinery may bind a narrow
  operation-local implementation temporarily, but the tag remains private and
  is not the runtime plan.

Any future mixed-surface operation must inspect the actual descriptors at the
operation boundary. If a formula changes, that change is a separate
physics/algorithm decision even if the surrounding orchestration is already
shared.

## 6. Staged implementation plan

Every production slice below has the same exact replay contract, called **R**:

- build the affected common library, ITS and MFT application/workflow targets,
  tests, and the combined workflow with the durable M6 build;
- run `ctest --test-dir <durable-build> -L itsmft --output-on-failure -j1`;
  every registered test must execute and pass;
- verify all 43 fixture checksums before and after;
- run standalone ITS and MFT common-CA replays and require 212 /
  `46913a67a7e2fe7462e29df0db264fa8` and 68 /
  `8106b08571ca593c6b76ff72b761a680` respectively;
- run the combined replay and require each combined detector leg to be
  byte-identical to its standalone leg;
- compare initialized writer leaves, sidecars, `CommonTrack`s, cluster
  references, ROFs, and labels with the accepted M6f output. The known
  uninitialized `MFTTrack.mInvQPtSeed` artifact remains excluded exactly as in
  the M6f/M6g records; it is not a defined output value;
- run the source/dependency guard and `git diff --check`; run a real device
  build only when `nvcc`/`hipcc` and the pinned toolchain exist.

The current M7a slice is documentation-only and therefore runs none of R.

### M7b — make the existing runtime plan the only count authority

**Status: complete; implementation details and final validation are recorded
in [design note 0004](0004-m7b-runtime-count-authority.md).**

**Target API and ownership:** retain the existing `DetectorLayoutSet`,
`SurfacePlanBinding`, and `SurfaceTrackingScratch` APIs; add focused contract
coverage rather than a wrapper. Convert common hot-loop bounds and temporary
array sizes that are semantically plan-sized to the binding ordered span,
layout active count, or scratch capacity.

**Temporary bridge:** `Tracker<NLayers>`, `TrackerTraits<NLayers>`, and dual
ROF/topology tables remain only at the explicitly classified M7c/M7d
compatibility boundary. The former `LayerMeasurementSpans<NLayers>` host alias
has been removed. No caller may use `NLayers` as a loop bound when a plan count
is available; `params.NLayers` is checked only at the adapter edge.

**Deletion criterion:** a source guard can identify every remaining
`NLayers` use in common production as either a type/ABI boundary, a fixed
device capacity, a private operation implementation, or an adapter seam; no
unclassified hot-loop bound remains.

**Validation:** R, plus the sparse/non-identity plan test and classified
source guard in design note 0004. The implementation is behavior-preserving;
the accepted candidate anchors and the known undefined
`MFTTrack.mInvQPtSeed` comparison exception remain unchanged.

### M7c — remove the layer-topology/ROF compatibility from the core path

**Target API and ownership:** use sparse topology IDs/views for CA traversal.
Keep ROF timing construction and raw-ROF ownership in the application/workflow;
the core receives runtime views/counts and does not own detector-specific
`ROF*Table<N>` objects. Preserve the current timing values and mask semantics.

**Temporary bridge:** an application adapter may construct the existing frozen
ROF table and expose its view while callers are converted. This is a
short-lived edge conversion, not a new common class or a second topology.

**Deletion criterion:** `SurfaceTrackingScratch` no longer needs dual ITS/MFT
members, `checkSupportedNLayers()`, or templated scratch dispatchers for
ROF/topology access; common traversal consumes `SparseTrackingTopologyView`
and runtime ROF views directly. The `IndexTableUtils<N>` alias is deleted once
its last caller is converted.

**Validation:** R, plus failure-lifecycle tests for stale ROF masks/timing
after load failure, dropped TF, structural tracking failure, and successful
replacement. Behavior-preserving. No raw-ROF workflow or publication state is
deleted.

### M7d — make `TrackerTraits` and `Tracker` non-templated

**Status: complete (2026-08-05);** production, tests/guards, and the
validation record are in [design note 0006](0006-m7d-nontemplated-tracker-core.md).

**Target API and ownership:** one non-templated `TrackerTraits` owns the shared
CA orchestration and one non-templated `Tracker` owns iteration execution.
They consume the existing runtime plan/workspace composition and fixed
`TrackSeed`; they contain no detector identity, source convention, writer,
DPL, workflow, or typed output track. `CATracker`, `TrackerITS`, and
`TrackerMFT` aliases and the 7/10 core instantiations are deleted.

The application participant may remain `SurfacePlanTrackingParticipant<7/10>`
while it owns sidecars and adapter configuration, but its members are plain
`TrackerTraits` and `Tracker`. Typed refit/output work crosses the core only
through a call-scoped operation reference; it is not stored by the core and is
the exact M7e deletion boundary. No compatibility class or second algorithm
implementation was introduced.

**Deletion criterion:** satisfied. The common production tree has no
`Tracker<7>`, `Tracker<10>`, `TrackerTraits<7>`, or `TrackerTraits<10>`
declaration, definition, or instantiation, and the source guard rejects the
retired template/alias spellings. This is the first deletion proving the
migration is real rather than a wrapper around the old tracker.

**Validation:** the M7d record applies R, the non-templated-core dependency
guard, the non-7/non-10 sparse-plan execution test, adapter participant
coverage, and failure/reset tests. The slice is behavior-preserving.

### M7e — move refit/output hooks to application adapters

**Target API and ownership:** the core accepts runtime measurement spans and
returns generic `SurfaceKinematicState`/`TrackSeed`/`CommonTrack` data. ITS and
MFT application adapters retain typed `TrackITSExt`/MFT conversion,
publication, pattern, polarity, and writer decisions. The native `Propagator`
path remains the only production refit path for this common tracker.

**Temporary bridge:** a concrete adapter may call the existing
`DetectorTraits<7/10>` implementation while it is being moved, but the bridge
must be a call-site migration with one implementation, not an alias or a
newly templated core. It is deleted in the same bounded adapter migration
before M7e closes.

**Deletion criterion:** `DetectorTraits<NLayers>`, `CATrackType<NLayers>`,
`LayerMeasurementSpans<NLayers>`, and the common core's typed accepted-track
vector have no production consumers. `CommonTrackOutputAdapter` and the
workflow-owned sidecars retain their existing output contract.

**Validation:** R, plus bit-level refit-leg/hole tests using a runtime span and
writer comparisons for all initialized fields. Behavior-preserving. Any
formula, covariance, material, acceptance, or MinPt change is out of scope and
requires classification 5 approval.

### M7f — delete redundant compile-time bridges and close the guard

**Target API and ownership:** delete the remaining common aliases, legacy
layer-index conversion helpers, redundant fixed-layer topology, explicit 7/10
core instantiations, and comments that describe the old ownership. Keep only
fixed device capacities and adapter-specific templates that have a separate
owner.

**Temporary bridge:** none in the core. Frozen ITS `ITStracking` code remains
outside the guard's scan scope and is not migrated.

**Deletion criterion:** the common core has exactly one runtime plan/workspace/
seed model; source scans show no production use of `Tracker<NLayers>`,
`TrackerTraits<NLayers>`, `TrackSeedTpl`, a detector switch, or an
`NLayers`-selected core behavior. The guard reports its narrow adapter and
frozen-legacy exclusions explicitly.

**Validation:** R, full source/dependency guard, and a real GPU build/replay
only if the toolchain is available. Behavior-preserving structural cleanup.

## 7. Deletion and simplification inventory

### Safe during M7

These are structural or already redundant once their direct callers are
migrated:

| Item | Why safe | Required proof |
|---|---|---|
| `CATracker`, `TrackerITS`, and `TrackerMFT` aliases | They add no behavior and only preserve the old template vocabulary. | No common or adapter call site relies on the alias after M7d. |
| `IndexTableUtils<N>` alias | It already names one `IndexTableUtilsCore` object and supplies no capacity. | All callers use `IndexTableUtilsCore`; index-table tests and R pass. |
| Explicit 7/10 instantiations of the non-templated core | They are generated ABI for a class whose behavior becomes one runtime implementation. | Link/build all affected targets and source guard. |
| `SurfaceTrackingScratch` templated dispatchers for old ROF/topology objects | They are forwarding selectors over dual members, not algorithm logic. | Runtime ROF/topology view tests and stale-state failure tests. |
| `CATrackTypeHelper`/typed accepted-track staging in generic core | It is a detector-output representation, not CA state. | `CommonTrack` and writer field parity against M6f. |
| `SurfaceId` to legacy-layer bridge arrays and layer-index loops | They are compatibility translation after the plan already has ordered surfaces. | Sparse/non-identity plan tests, seed/hole/refit order parity. |
| Duplicate comments and helper names that claim `NLayers` is generic storage | They become actively misleading after M7d. | Documentation/source guard review. |

### Blocked until M7 completion

These should not be deleted before the preceding runtime-plan slices provide
their replacement and the full replay gate:

| Item | Blocking dependency |
|---|---|
| `Tracker<NLayers>`/`TrackerTraits<NLayers>` declarations | M7b plan-count conversion and M7d's single runtime implementation. |
| Dual `ROFOverlapTable<N>`/`ROFVertexLookupTable<N>`/`ROFMaskTable<N>` storage | M7c application-owned timing views and reset tests. |
| `TrackingTopology<N>` layer graph | Sparse topology consumers and device-view parity. |
| `DetectorTraits<N>` and `TrackingLoadPolicyN` | M7e application adapter extraction and output/refit parity. |
| `LayerMeasurementSpans<N>` and runtime refit slot assembly | M7e runtime-span refit tests. |
| Common `SurfacePlanTrackingParticipant<N>` template | Non-templated core first; then sidecar/configuration ownership can be moved to concrete adapters. |

### Deferred physics or algorithm choices

These are not simplifications to smuggle into M7:

- changing the native `Propagator` equations, covariance sanitization,
  material formulas, propagation direction, refit acceptance, MinPt indexing,
  hole semantics, or candidate ordering;
- removing or changing the MFT road pre-cut, legacy reference-Z lookup, or
  curvature-literal behavior;
- replacing `TransitionPolicyTag` with a different public dispatch taxonomy;
- changing `TrackSeed`/`SurfaceMask`/`LayerMask` widths or fixed capacities, or
  introducing dynamic device allocation;
- making arbitrary mixed cylinder/disk plans legal when current validation
  rejects them;
- unifying independent refit arithmetic or deleting a legacy path based only on
  numerical similarity; this needs an A/B report and a separate decision;
- removing frozen ITS tracker, vertexer, ROF, or CUDA code;
- deleting `TimeFrame::mNTotalLowPtVertices` or other support state before a
  dedicated whole-repository liveness audit confirms it is dead.

## 8. Guard and test contract for implementation

The M7 guard should scan
`Detectors/ITSMFT/common/tracking/{include,src}` and distinguish three narrow
sets instead of blanket-exempting the whole tree:

1. **Core scan:** `Tracker`, `TrackerTraits`, `CATracker`,
   `SurfaceTrackingScratch`, `IndexTableUtils`, sparse topology, refit
   orchestration, `TrackingEngine`, and `TrackingParticipant`. This scan must
   reject `NLayers`-selected behavior, `DetID`/source-0/1 switches, DPL,
   workflow, writer, and typed detector output dependencies once M7f closes.
2. **Adapter transition scan:** `TrackingInterface`,
   `SurfacePlanTrackingParticipant`, configuration/load policy, and concrete
   ITS/MFT output adapters. These may contain compatibility templates only
   until their stated M7 deletion criterion, and the guard must name each
   path explicitly.
3. **Frozen legacy exclusion:** `Detectors/ITSMFT/ITS/tracking` and other
   detector-specific frozen legacy code are outside this common-CA guard. The
   exclusion must be path-specific; it must not exempt common tracking
   production files.

Focused tests required by the implementation slices are:

- a runtime-plan contract test with non-7/non-10 active counts and sparse
  surface IDs, checking ordered traversal, source-qualified measurements,
  transition/cell slots, and `TrackSeed` masks;
- a compile/dependency test proving the non-templated core has no detector,
  DPL, workflow, writer, or public `TransitionPolicyTag` dependency;
- fixed-capacity tests for `CellSeed`, `TrackSeed`, `SurfaceMask`, and device
  views, including the distinction between active count and maximum capacity;
- failure lifecycle tests proving load failure, dropped TF, structural
  tracking failure, and successful replacement cannot publish stale clocks,
  validity flags, sidecars, or normalized data;
- runtime-span refit tests covering forward/reverse traversal, holes,
  acceptance gates, and the existing candidate baseline;
- schedule/order tests proving combined ITS/MFT source-qualified exports stay
  unchanged while the core sees only an explicit participant schedule.

## 9. M7a/M7b/M7c/M7d exit record

M7a was complete when this note and the migration-plan navigation were
committed as documentation only. M7b, M7c, and M7d are complete only with
their corresponding production/test/documentation commits and replay gates.
M7d's remaining boundary is deliberately narrow: M7e must remove the
call-scoped typed refit/output operation seam and relocate `DetectorTraits<N>`
conversion and sidecar sealing to the application adapters. The accepted
physics anchors remain ITS 212 / `46913a67a7e2fe7462e29df0db264fa8` and MFT
68 / `8106b08571ca593c6b76ff72b761a680`.

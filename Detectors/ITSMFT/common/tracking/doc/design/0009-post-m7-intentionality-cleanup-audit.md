# Post-M7 intentionality and cleanup audit

- Status: audit complete; implementation not started
- Date: 2026-08-06
- Audited revision: `51a49e5d4a` (`codex/itsmft-integration`, after M7f and P1)
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

The intended end state does not need a new wrapper or manager. It should:

1. delete demonstrably dead typed/refit and forwarding residue;
2. make the existing loader directly target `SurfaceTrackingScratch` and
   remove fixed ITS+MFT wrapper entry points;
3. move reset ownership entirely to the event owner;
4. evolve the existing participant in place into the one plan-bound
   application composition used by standalone and combined workflows;
5. delete `TrackingInterface`;
6. narrow the operation seam to refit only and move publication compatibility
   fully to application/workflow code; and
7. merge the implementation currently named `TrackerTraits` into `Tracker`
   once tests no longer construct it independently.

`TrackingEngine` and `Tracker` should remain distinct. The former owns an
ordered multi-participant event transaction; the latter owns one
participant's CA algorithm. `TrackingInterface` should be deleted, not merged
into either. `TrackerTraits` has no traits role after M7d and should not remain
a public peer indefinitely.

The first implementation slice is unambiguous: remove the unused typed MFT
refit/export path (`MFTAdapterRefit.{h,cxx}` and `MFTCATrack.h`) and migrate its
fixture to the already-live generic native-refit result path. Section 8.1
defines the exact boundary and exit criterion.

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
- `StateFamily` remains representation metadata, never schedule or topology;
- `DetectorLayoutView`, sparse topology, `SurfacePlanBinding`, and
  `SurfaceTrackingScratch` remain the runtime authorities;
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
        |     +-- adapter refit -> generic CommonTrack
        |     +-- adapter accepted-result compatibility
        +-- translate TrackingResult to float/drop sentinel
        +-- workflow consumes CommonTrack and sidecars
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
  +-- atomic loader -> TimeFrame + participant scratches
  +-- one or more plan-bound participants in explicit schedule
  +-- TrackingEngine -> whole-event success/failure/reset
        |
        +-- Tracker -> one participant's CA/refit/CommonTrack algorithm
        +-- application publication adapter consumes successful generic result
```

Only the schedule cardinality differs between standalone and combined. A
standalone task uses one participant; it does not need a second orchestration
class or a float sentinel.

## 4. Dependency map of the central cluster

Solid arrows below mean construction/ownership or a required call. “Adapter”
arrows are detector compatibility dependencies that should move outward.

```text
Combined DPL task ---------------------> TrackingEngine
       |                                      |
       |                                      v
       +--> MultiSourceTimeFrameLoader   TrackingParticipant
       |          |                           ^
       |          v                           |
       |   SurfaceTrackingScratch <--- SurfacePlanTrackingParticipant<7/10>
       |                                      |
       |                                      +--> Tracker --> TrackerTraits
       |                                      |      |            |
       |                                      |      +------------+
       |                                      |       TimeFrame / binding /
       |                                      |       scratch / sparse layout
       |                                      |
       |                                      +--adapter--> ROF fixed tables
       |                                      +--adapter--> DetectorPublicationAdapter
       |
Standalone DPL task --> TrackingInterface<7/10>
                              |
                              +--> same Tracker + TrackerTraits + scratch
                              +--> same ROF tables + publication adapter
                              +--> loading/reset/clock/sentinel policy

TimeFrame <----- normalized measurements and CommonTrack results
SurfacePlanBinding <----- immutable plan/source/compact-slot mapping
ROFViews <----- borrowed by core; built and owned at adapter edge
TrackingOperationAdapter <----- refit + accepted compatibility + reset seam
```

The highest-cost dependencies are not algorithmic. They are the inclusion of
`TrackingInterface.h` by `SurfacePlanTrackingParticipant.h` solely to reuse
sidecar-owner mixins, and the public detector-specialized
`DetectorPublicationAdapter.h` beneath the nominally common include root.

## 5. Responsibility audit of the focal modules

| Module | Actual single responsibility and data relationship | Classification | Duplication and disposition | Risk |
|---|---|---|---|---|
| `Tracker.h` | Includes `CATracker.h`; owns and implements nothing. | Compatibility residue | Delete and include `CATracker.h` directly at its sole production consumer. | Safe include cleanup; compile all public-header users. |
| `CATracker.h` / `Tracker` | Orchestrates configured CA iterations for one bound plan. Borrows traits, scratch, frame, binding/layout; owns parameter copies and a shared pool handle. | Generic core | Keep the one-participant algorithm role. Remove float sentinel after standalone migration. Move event reset out. Eventually absorb `TrackerTraits`. | Reset migration is exception-sensitive and needs failure-contract tests plus replay. |
| `TrackerTraits` | Implements tracklet, cell, neighbour, road, refit/acceptance, traversal caches, and operation binding. Borrows plan/binding/scratch/frame; owns iteration caches and accepted-candidate staging. | Generic core implementation | Unique algorithm body, but not “traits” and not an independent production strategy. Narrow, then merge into `Tracker`; do not replace it with another interface. | Broad build/test surface; physics replay required even for mechanical merge. |
| `TrackingEngine` | Executes an explicit heterogeneous participant schedule and owns whole-event failure reset. It owns no event data. | Generic event core | Keep distinct and stateless. Its responsibility is unique and above one-leg tracking. Shorten migration-era comments. | Low if API unchanged. |
| `TrackingParticipant` | Minimal dynamic seam for a plan-bound application leg: identity, owned surfaces, track/reset/export. Borrows the shared frame. | Generic application seam | Keep while combined scheduling is heterogeneous. Narrow publication export if its surface list remains unused. | ABI/caller audit; combined contract tests. |
| `TrackingInterface<7/10>` | Standalone application composition and event coordinator. Owns frame, scratch, tracker, plan/binding, decoder, ROF state, clocks, sidecars, and publication adapters. | Application adapter plus compatibility residue | Duplicates participant composition and workflow lifecycle. Migrate standalone tasks to one participant plus loader/engine, then delete; do not merge into `Tracker` or `TrackingEngine`. | High lifecycle/output compatibility risk; stage under standalone and combined replay gates. |
| `SurfacePlanTrackingParticipant<7/10>` | Combined plan-bound leg. Owns scratch/binding/tracker composition, fixed ROF compatibility tables, load target, sidecars, and tracked flag; borrows plan/frame. | Application adapter | Coherent participant role, but detector count, loading target, timing tables, and publication compatibility make it too broad. Evolve in place into one runtime-plan participant used by standalone and combined; do not add a parallel participant. | Bounded adapter migration; exact writer/sidecar replay required. |
| `TrackingOperationAdapter` | Supplies seed refit, accepted-result completion, and adapter reset. It borrows generic candidate/scratch/measurement views. | Mixed core/adapter seam | Three responsibilities. Narrow to the one operation the algorithm needs (refit); run publication conversion/reset after generic success at participant/workflow edge. Delete if refit can be directly owned without detector dispatch. | Call order and rejection classification are physics-sensitive; preserve exact boundaries. |
| `DetectorPublicationAdapter<7/10>` | Specializes accepted-result publication compatibility for ITS shared-cluster flags and MFT sidecars. | Detector application adapter | Responsibility is real but location/visibility is wrong. Move to ITS/MFT application code and give each workflow direct ownership. No generic dispatcher. | Writer/sidecar compatibility gate. |
| `TimeFrame` | Owns normalized multi-source event data, vertices/labels, beam/Bz, and generic `CommonTrack` results. | Generic event data owner | Keep distinct from scratch. Audit and remove inert device-propagator virtual hook separately. | Public ABI/device tests for hook removal. |
| `SurfaceTrackingScratch` | Owns runtime-sized mutable CA workspace and per-source legacy-compatible measurement backfill. Borrows allocator/pool. | Generic participant workspace | Keep distinct from `TimeFrame`. Delete dead `mPValphaX`; collapse `resetScratch()` alias; make binding mandatory before removing fallback mappings. | Allocator identity, reset, and sparse-plan tests. |
| `SurfacePlanBinding` | Immutable source-qualified projection from layout topology to one participant's ordered surfaces and compact transition/cell slots. | Generic plan binding | Keep. Apparent copied IDs are the participant-scoped mapping, not a duplicate topology. | Low; sparse/non-identity tests remain mandatory. |
| `MultiSourceTimeFrameLoader` | Provides atomic normalized-frame plus scratch-backfill staging/commit. | Generic load transaction | Keep transaction, remove virtual `LoadTarget` hierarchy and fixed `loadITSAndMFT`/reset wrappers. Bind source directly to scratch. | Atomicity/allocator/failure-retry tests and combined replay. |
| `ROFViews` | Defines non-owning runtime timing, overlap, vertex-lookup, and mask views consumed by common tracking. | Generic runtime input view | Keep. Fixed detector tables belong at adapter/workflow edge and should not migrate into scratch/frame. | Timing/diamond/mask parity tests. |

## 6. Duplicate lifecycle and representation findings

### 6.1 Initialization and construction

`TrackingInterface` and `SurfacePlanTrackingParticipant` each build the same
`TrackerTraits`/`Tracker`/scratch/binding/pool composition. Both also select
ITS/MFT fixed timing tables and publication compatibility through `NLayers`.
The combined participant should become the single composition, configured by
an application-supplied plan and explicit adapter-owned timing/publication
objects. Standalone workflows can schedule one instance.

No duplicate runtime topology was found. `DetectorLayoutView` owns global
sparse topology; `SurfacePlanBinding` owns the necessary participant-scoped
projection; scratch owns mutable compact-slot data. These should not be
merged.

### 6.2 Loading

Three loading shapes remain:

1. `TrackingInterface::processTimeFrame()` loads one source directly through
   `SurfaceTrackingScratch::loadNormalizedSource()`;
2. `MultiSourceTimeFrameLoader::loadEvent()` performs the production atomic
   multi-source transaction;
3. `loadITSAndMFT()` is a fixed-position wrapper used only by tests, along
   with `resetITSAndMFTEvent()`.

The nested virtual `LoadTarget` hierarchy has one implementation and every
production target is a `SurfaceTrackingScratch`. Replacing
`AtomicLoadBinding::{source, LoadTarget&}` with
`{source, SurfaceTrackingScratch&}` preserves the transaction while deleting
virtual dispatch, participant `mLoadTarget`, friendship-driven forwarding,
and one lifetime constraint.

### 6.3 Reset and failure handling

Reset is currently implemented at four levels:

- `SurfaceTrackingScratch::reset()` and its forwarding alias
  `resetScratch()`;
- `resetTimeFrameEvent(frame, scratch)`;
- `TrackingInterface::resetEvent()` plus operation-adapter compatibility
  reset;
- `TrackingEngine::resetEvent()` plus participant `eventReset()`.

The semantic owner should be the DPL event path via `TrackingEngine`.
Participant reset clears participant-local scratch and adapter state;
`TrackingEngine` wipes shared `TimeFrame` exactly once. `Tracker` should
report typed failure without clearing either. This cannot be a casual
cleanup: exception paths, dropped-TF classification, retryability, and
publication invalidation must be pinned first.

### 6.4 Tracking and publication

There is one CA algorithm body, but two adapter entry paths. Publication is
also split: generic `CommonTrack` append occurs in the tracker path,
`TrackingOperationAdapter::completeAccepted()` stages compatibility, the
participant exposes a publication export, and the workflow builds typed
outputs. The intended split is simpler:

- `Tracker` produces generic accepted results and `CommonTrack` references;
- a detector application adapter converts the successful generic result;
- the workflow owns validity, clocks, raw ROFs, and final publication.

Adapter publication must never be called from a failed generic transaction.
Moving it outward requires preserving serial accepted-result order and the
current “final” boundary.

## 7. Rest-of-public-include-tree audit

Every public common-tracking header is covered below. Grouping means the
listed files share one disposition; it does not imply a new module.

| Files/modules | Current owner and recommendation |
|---|---|
| `BarrelSurfaceStateOperations.h`, `ForwardSurfaceStateOperations.h`, `SurfaceStateOperationResult.h`, `MaterialPhysics.h`, `Propagator.h`, `NativeRefitDriver.h`, `RefitLegAssembly.h`, `SurfaceLinearizationReference.h` | Narrow descriptor/state leaf operations and native refit. Generic core; retain. Any formula, covariance, material, leg, or ordering consolidation needs physics approval. |
| `NativeCylinderCylinderRefitDriver.h` | Pre-activation comparison driver with no production caller; tests and P1-era characterization still name it. Archive unique evidence, then delete or fold only genuinely unique test helpers into the live driver. |
| `Cell.h`, `SeedAnchor.h`, `LayerMask.h`, `SurfaceMask.h`, `SurfaceMeasurementIndex.h` | Device-safe CA value types and fixed-capacity metadata. Retain fixed capacity with runtime active prefixes. Do not remove live `SeedMetadataBase<N>` used by `CellSeed`. |
| `CommonTrack.h`, `CommonTrackShadow.h`, `CommonTrackOutputAdapter.h`, `MCLabelAccumulator.h` | Generic result and adapter conversion utilities. Keep `CommonTrack`; move detector-shaped output helpers outward if their callers are only workflows. Audit `CommonTrackShadow` after P1 tooling no longer needs it. |
| `ClusterSource.h`, `DecodedCluster.h`, `ClusterDecoding.h`, `ClusterDecoder.h`, `MultiSourceFrame.h`, `MultiSourceLoading.h`, `TimeFrameLoadFailure.h` | Generic normalized multi-source input and transactional decoding. Retain. Narrow includes where implementation-only decoder dependencies leak into public headers. |
| `ClockTimingPublicationView.h`, `SurfaceTiming.h`, `ROFTimingUniformity.h`, `ROFViews.h` | Runtime timing views and validation. Retain generic views; timing construction and publication clocks stay workflows/adapters. |
| `Configuration.h`, `TrackingConfigParam.h`, `ConfigKeyValuesPreflight.h`, `IndexTableConfiguration.h`, `IndexTableUtils.h` | Algorithm parameters/configuration and fixed-capacity index tables. Retain live runtime-prefix storage. Detector-named configuration registration is adapter compatibility, not a core routing authority; move only under a separately gated configuration migration. |
| `DetectorLayout.h`, `DetectorLayoutBuilder.h`, `DetectorLayoutSet.h`, `SparseTrackingTopology.h`, `SurfaceCatalogView.h`, `SurfaceDescriptor.h`, `SurfaceId.h`, `StaticSurfaceDescriptor.h`, `StaticDetectorCatalogs.h`, `SurfaceSpec.h` | Authoritative generic descriptors and sparse plan/topology. Retain. `StaticDetectorCatalogs` is application plan data under a common namespace; consider moving ITS/MFT tables outward, but do not duplicate the builder. |
| `ITSSurfaceSpec.h`, `MFTSurfaceSpec.h`, `NominalSurfaceMaterialDefaults.h` | ITS/MFT application data. Valid compatibility owners, but not generic core concepts. Move to detector application include locations in a bounded include/API migration. |
| `SurfaceKinematicState.h`, `StateFamily.h` | Generic state representation. Retain; `StateFamily` must not become dispatch policy. |
| `SurfaceMeasurement.h` | Generic normalized measurement. Retain. |
| `SurfaceKinematicStateLegacyAdapters.h`, `SurfaceMeasurementAdapters.h`, `TrackingFrameInfoAdapters.h` | Compatibility conversions. The first two remain live at adapter/decoder leaves. `TrackingFrameInfoAdapters` and common `loadClusterTrackingFrameInfo` have no production caller found and are test-preserved deletion candidates. |
| `IOUtils.h` | Mixed decoder covariance/systematic helpers plus obsolete public conversions. Split only by deleting proven-dead declarations; avoid a new utility namespace hierarchy. Common `convertCompactClusters` needs a whole-repository/public-API check before deletion despite no current common production caller. |
| `ITSSharedClusterCompatibility.h`, `MFTPublicationCompatibility.h`, `DetectorPublicationAdapter.h`, `DetectorTrackingOperationAdapterSupport.h`, `MFTFwdTrackHelpers.h` | Detector application/refit/publication compatibility. Responsibilities can remain, but ownership should move outside generic core headers as participant/interface consolidation proceeds. Do not replace with a central detector dispatcher. |
| `MFTAdapterRefit.h`, `MFTCATrack.h` | Typed MFT compatibility path with no production consumer; only `testMFTNormalizedRefit` and guards retain it. Delete first. |
| `ParticipantId.h`, `TrackingParticipant.h`, `TrackingEngine.h` | Generic heterogeneous schedule contract. Retain, narrow comments and any unused export fields after caller audit. |
| `SurfacePlanTrackingParticipant.h`, `TrackingInterface.h`, `TrackingOperationAdapter.h` | Application-composition cluster described in Sections 5–6. Consolidate by evolving the participant, deleting the interface, and narrowing the operation seam. |
| `TimeFrame.h`, `SurfaceTrackingScratch.h`, `detail/SurfacePlanBinding.h`, `MultiSourceTimeFrameLoader.h` | Distinct event owner, participant workspace, immutable participant projection, and atomic loader. Retain responsibilities; simplify dead fields, forwarding aliases, and loader type erasure. |
| `CATracker.h`, `TrackerTraits.h`, `Tracker.h` | One CA implementation plus forwarding residue. Delete `Tracker.h`; eventually merge traits implementation into `Tracker`. |
| `detail/TransitionPolicy.h`, `detail/TransitionPolicyBinding.h`, `detail/TransitionPolicyDispatch.h`, `detail/TransitionPolicyOperations.h`, `detail/TransitionPolicyState.h` | Private compatibility implementation for descriptor-family leaf selection. Keep private. Structural collapse is not safe if it changes operation order or arithmetic. |

## 8. Ranked cleanup inventory

### 8.1 Safe behavior-preserving cleanup

“Safe” means no intended behavior change, not that validation may be skipped.

| Rank | Exact files/action | Current callers | Replacement owner | Required gate | Deletion criterion |
|---:|---|---|---|---|---|
| 1 | Delete `include/ITSMFTTracking/MFTAdapterRefit.h`, `include/ITSMFTTracking/MFTCATrack.h`, `src/MFTAdapterRefit.cxx`, and its CMake entry. Migrate `testMFTNormalizedRefit.cxx` away from typed export. | No production caller; only `testMFTNormalizedRefit`, M7 guards, and stale CMake comments. | Existing generic `MFTFwdTrackHelpers`/`NativeRefitDriver` result path; MFT workflow publication already consumes `CommonTrack`/adapter compatibility. | Focused native-refit and MFT publication tests, full serial `itsmft` suite, 43/43 fixture checks, 212/68 hashes, combined-leg and writer parity. | Repository search finds no common `MFTCATrack` type or `MFTAdapterRefit`; no production target loses a symbol. |
| 2 | Delete forwarding-only `Tracker.h`; include `CATracker.h` directly and remove stale references. | `TrackingInterface.h` is the sole production include found. | `CATracker.h`. | Build all common/ITS/MFT/combined targets; header dependency guard. | No include/reference to `ITSMFTTracking/Tracker.h`. |
| 3 | Delete scratch member `mPValphaX` and allocator/reset/swap bookkeeping. | No read or semantic write found; only initialization, clearing, allocator-match, and swap. | None. | Scratch allocator/swap/reset tests, sanitizer-capable focused tests, full suite and replay. | No field/reference remains and staged/live allocator matching still passes. |
| 4 | Use one scratch reset name; delete forwarding `resetScratch()` after changing direct callers to `reset()`. | Participant, loader wrapper, and migration-era tests. | `SurfaceTrackingScratch::reset()`. | Reset ordering, retry, dropped/structural failure tests. | One public scratch reset method and no semantic change in ordering. |
| 5 | Delete `loadITSAndMFT()` and `resetITSAndMFTEvent()` after migrating their tests to `loadEvent()` plus engine reset. | No production caller found; fixed-source behavior is test-only. | Generic `loadEvent()` and `TrackingEngine::resetEvent()`. | Atomic load success/failure/retry and combined isolation tests. | No fixed ITS=0/MFT=1 loader API under common tracking. |
| 6 | Correct migration-era comments and test names that claim templated `Tracker`/`TrackerTraits`, M2/M6 temporary ownership, or an unfinished M7e seam. | Documentation and source comments/tests only. | Current responsibility wording. | Source guard plus documentation link/headings validation. | No stale architectural claim; behavioral assertions retained under intent-based names. |

The first cleanup slice is rank 1 only. Combining later ranks would make
review and regression localization worse.

### 8.2 Requires a bounded migration

| Rank | Exact files/action | Current callers | Replacement owner | Required gate | Deletion criterion |
|---:|---|---|---|---|---|
| 1 | Replace `MultiSourceTimeFrameLoader::{LoadTarget,LoadTargetImplSurface}` and participant `mLoadTarget` with atomic bindings that borrow `SurfaceTrackingScratch&`. | Combined workflow, participant, loader tests. | `MultiSourceTimeFrameLoader::loadEvent()` directly stages one scratch per binding. | Allocator identity, partial-stage failure, retry, dropped TF, combined source isolation, full replay. | No virtual load target, friendship-only forwarding, or participant load-target member. |
| 2 | Make `SurfacePlanBinding` mandatory in `TrackerTraits`; delete identity/source-0 fallback mapping. | Direct unit fixtures and tracker construction; production already binds. | Caller-owned validated binding. | Non-identity/sparse plan tests, invalid-binding failure tests, full replay. | No nullable binding branch or synthesized numeric-surface traversal in core. |
| 3 | Remove reset from `Tracker::clustersToTracks()` and make it return/throw typed failure without mutating event ownership. | `TrackingInterface`, participant, engine. | `TrackingEngine` as whole-event reset owner; participant clears only local state. | Exact success/recoverable/structural/exception reset-count tests, stale sidecar tests, standalone and combined replay. | Each failure path clears participant state once and shared frame once; tracker contains no `wipe`/event reset. |
| 4 | Move fixed ROF tables, masks, and publication helpers out of `SurfacePlanTrackingParticipant` and `TrackingInterface`. | Standalone and combined workflows. | ITS/MFT workflow/application setup; core continues borrowing `ROFViews`. | Empty/first/last/diamond timing, mask, load-failure replacement, writer and sidecar parity. | Participant owns no fixed detector table or detector publication sidecar. |
| 5 | Evolve `SurfacePlanTrackingParticipant<7/10>` in place into one runtime-plan participant; remove detector-count static assertion, aliases, mixins, and explicit instantiations. | Combined workflow; future standalone after next rank. | Existing participant configured by runtime plan and a narrow application refit operation. | Valid neither-7-nor-10 plan test, ITS/MFT schedule/isolation, all failure classes, full replay. | One non-templated participant with no detector identity or sidecar ownership. |
| 6 | Migrate standalone ITS/MFT DPL tasks from `TrackingInterface<7/10>` to loader + one participant + `TrackingEngine`; delete `TrackingInterface.{h,cxx}` and float drop sentinel. | ITS and MFT CA workflow tasks and interface-heavy tests. | Workflow owns raw ROFs/clocks/publication; shared participant/engine own tracking. | Standalone lifecycle/config/output contract tests, dropped-TF behavior, combined parity, exact writer/replay baseline. | No `TrackingInterface`, duplicate composition, or float-sentinel translation remains. |
| 7 | Narrow `TrackingOperationAdapter`: move `completeAccepted()` and `resetAdapterState()` to application publication handling; retain only the exact refit operation if still needed. | `TrackerTraits`, interface, participant. | Tracker owns generic acceptance; application adapter consumes only successful ordered generic results. | Candidate ordering, holes, all refit rejection paths, shared/pattern compatibility, failure-stale-state and replay gates. | Core operation seam neither stages publication nor owns lifecycle; delete seam entirely if one generic refit implementation remains. |
| 8 | Merge `TrackerTraits` implementation into `Tracker`, or make it translation-unit-private as an intermediate step; remove virtual methods and public header. | Tracker and many direct unit fixtures; no alternate production implementation. | `Tracker` as the sole one-participant CA algorithm class. | Migrate tests to public runtime-plan behavior, build CPU/device configurations, full suite and replay. | One public algorithm class, no `TrackerTraits` type, no second implementation body. |
| 9 | Remove `TimeFrame`'s no-op virtual device-propagator hook/member if full build configurations confirm no common-derived device frame. | A migration-era ownership test; frozen ITS has a separate live override in a different class. | None in common CPU runtime; a future device owner must be explicit. | Whole-repository inheritance/API audit and actual device build when pinned toolchain exists. | Common `TimeFrame` is non-polymorphic and contains no unused device pointer. |
| 10 | Delete common `TrackingFrameInfoAdapters` and `loadClusterTrackingFrameInfo`; separately decide common `convertCompactClusters`. | Tests only for the former; no current common production caller found. Public API risk remains. | `SurfaceMeasurement` decoder path for production; detector-local frozen utilities remain untouched. | Whole-repository symbol/header audit, downstream build, covariance/systematics tests. | No downstream consumer and no lost compatibility contract. |
| 11 | Retire `NativeCylinderCylinderRefitDriver.h` after preserving the comparison evidence it uniquely supplies. | Dedicated test, guard, and historical design/P1 harness. | Live `NativeRefitDriver`; durable validation docs/artifacts for historical evidence. | Demonstrate all production numerical/refit properties remain covered; no replay delta. | No production or active validation need for the unwired driver. |

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
- `testM7eAdapterBoundaryGuard` currently proves typed MFT compatibility is
  isolated in `MFTAdapterRefit`; after rank 1 it should prove absence;
- `testNativeCylinderCylinderRefitDriver` is a historical comparison harness,
  not production-path coverage; preserve its evidence before retirement.

The final M7 runtime-core guard remains valuable. Cleanup should tighten it
around current responsibilities rather than accumulate named exemptions for
deleted migration bridges.

Public include cost can be reduced without a new pimpl/service layer:

- delete `Tracker.h`;
- stop including `TrackingInterface.h` from the participant merely for
  sidecar-owner mixins;
- move detector-specific specs/publication helpers out of generic public
  includes as their callers migrate;
- keep private transition-policy headers private to implementation-facing
  headers; and
- merge `TrackerTraits` only after direct test construction is removed.

## 10. Explicit class recommendation

| Class | Recommendation | Reason |
|---|---|---|
| `TrackingEngine` | **Remain distinct and narrow.** | It is the stateless whole-event schedule/failure boundary. Merging it into a participant algorithm would reintroduce coordinator-shaped detector/application concerns. |
| `Tracker` | **Remain as the single one-participant CA algorithm class; narrow its lifecycle role.** | It uniquely sequences iterations and CA stages. It should not wipe whole-event state or translate workflow sentinels. |
| `TrackerTraits` | **Merge into `Tracker` after a bounded fixture migration.** | It is neither a traits policy nor an alternate implementation; it is the bulk of the same algorithm. Keeping it public creates two apparent core trackers and heavy private-policy exposure. |
| `TrackingInterface` | **Delete after standalone migration.** | It duplicates participant composition and owns workflow lifecycle that belongs to the DPL task. Merging it into the engine would contaminate the generic boundary. |
| `SurfacePlanTrackingParticipant` | **Retain the role, evolve the existing class in place, and make it runtime-plan/non-templated.** | A plan-bound leg is the necessary bridge between an application and `TrackingParticipant`; its detector timing/publication/load-target holdings are not necessary. |
| `TrackingOperationAdapter` | **Narrow to refit-only, then reassess deletion.** | Refit is an algorithm operation; publication completion and reset are application lifecycle. Do not replace it with a callback framework. |

## 11. Ordered cleanup slices

Each slice deletes its bridge before the next one starts; no long-lived dual
implementation is proposed.

1. **C0 — dead typed-MFT refit/export retirement.** Delete the three rank-1
   files and migrate the typed fixture to generic native-refit assertions.
2. **C1 — zero-risk residue.** Delete `Tracker.h`, `mPValphaX`, and the
   duplicate scratch reset name in separate reviewable commits.
3. **C2 — loader simplification.** Replace the one-implementation virtual
   load target with direct scratch bindings; remove fixed ITS+MFT wrappers.
4. **C3 — binding and reset authority.** Make binding mandatory, then move
   reset out of `Tracker` under exact failure-count tests.
5. **C4 — application ownership.** Move ROF fixed tables and detector
   publication compatibility out of the participant/interface, preserving
   runtime `ROFViews` in core.
6. **C5 — one participant composition.** De-template the existing
   participant in place and use it for both combined and standalone tasks.
7. **C6 — standalone interface retirement.** Delete `TrackingInterface` and
   float sentinel once one-participant engine execution is live.
8. **C7 — operation seam narrowing.** Move accepted-result compatibility and
   reset outward; retain only a necessary refit operation.
9. **C8 — one visible tracker.** Merge `TrackerTraits` into `Tracker` and
   remove direct traits fixtures/public header.
10. **C9 — public compatibility tail.** Resolve dead common IO adapters,
    inert device hooks, old comparison driver, and detector-spec header
    placement one item at a time.

No slice above authorizes a physics simplification. Section 8.3 candidates
must remain separately gated even if a class merge makes them tempting.

## 12. First-slice implementation contract

The next task can be bounded exactly as follows:

**Input scope**

- `include/ITSMFTTracking/MFTAdapterRefit.h`
- `include/ITSMFTTracking/MFTCATrack.h`
- `src/MFTAdapterRefit.cxx`
- common tracking CMake source list and stale comments
- `testMFTNormalizedRefit.cxx`
- M7 boundary/guard tests that currently whitelist the typed path

**Required result**

- no production or test target includes or links the typed MFT refit/export;
- meaningful refit tests assert the existing generic
  `SurfaceKinematicState`/`TrackingCandidate` result and unchanged failure
  behavior rather than recreating `MFTCATrack`;
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

A repository search finds no `MFTAdapterRefit`, common `MFTCATrack` class, or
typed refit overload, and no replacement compatibility wrapper has been
introduced. This deletion is independent of reset, participant, interface,
or physics changes.

## 13. Decision record

Measured structure:

- there is one production CA body;
- there are two application compositions;
- reset policy is duplicated across tracker and event owners;
- the load target is one-implementation type erasure;
- `TrackerTraits` has one production implementation and no current traits
  selection role;
- typed MFT refit/export files have no production caller.

Inference:

- deleting dead typed compatibility first is lower risk than beginning with
  lifecycle consolidation;
- converging standalone onto the existing participant/engine path removes
  more complexity than a new facade would;
- retaining `TrackingEngine` separately protects the generic event boundary;
- eventually merging `TrackerTraits` into `Tracker` makes the core easier to
  understand without changing its algorithm.

No physics inference or acceptance decision is made by this audit.

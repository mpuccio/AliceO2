# Generic Tracking-Engine Migration Plan (post-Gate-4)

Status: Living plan
Decision anchor: [ADR 0007](decisions/0007-generic-tracking-engine-boundary.md)
Architecture: [Architecture.md](Architecture.md)
Status/decision log: [AgentCoordination.md](AgentCoordination.md)
M7a design/audit: [runtime-plan Tracker/TrackerTraits migration](design/0003-m7-runtime-plan-tracker-migration.md)
M7b count authority: [runtime-count authority](design/0004-m7b-runtime-count-authority.md)
M7c topology/ROF ownership: [runtime topology and ROF ownership](design/0005-m7c-runtime-topology-rof.md)
M7d non-templated core: [runtime-plan Tracker/TrackerTraits implementation](design/0006-m7d-nontemplated-tracker-core.md)
M7e adapter refit/output: [typed refit and output at the adapter boundary](design/0007-m7e-adapter-refit-output.md)
M7f final cleanup: [redundant runtime-core bridge retirement](design/0008-m7f-final-runtime-core-cleanup.md)
Post-M7 cleanup audit: [intentionality, ownership, and duplication audit](design/0009-post-m7-intentionality-cleanup-audit.md) — complete; bounded first slice is dead typed-MFT refit/export retirement.
Post-M7 cleanup implementation: [sequential Luna execution plan](design/0010-post-m7-cleanup-implementation-plan.md) — L9+L10 implementation complete; final replay gate open. [L10 operation-tail validation](validation/l10-operation-tail.md) records the refit-only seam, deleted compatibility paths, and final campaign requirements; [L9 TrackingInterface-retirement validation](validation/l9-retire-interface.md) records direct standalone composition; [L8 Tracker-orchestration validation](validation/l8-retire-engine.md) records engine/participant retirement and direct combined composition.
Header/graph follow-up: [header intentionality and SurfaceGraph simplification audit](design/0011-header-intentionality-and-surface-graph-audit.md) — documentation-only inventory and staged design; implementation not started.
Atomic loader cleanup: [direct TimeFrame-source loading](validation/p1-direct-atomic-loader.md).

L6 makes `TimeFrame` the sole live owner of generic configuration, workspace,
and event data, with `loadTimeFrameSources()` as the non-owning atomic load
operation. Standalone and combined adapters provide normalized inputs and
runtime ROF views; raw ROFs, timing backing storage, publication validity,
sidecars, and writers remain outside the frame. L7 relocated those remaining
compatibility owners to workflow edges; L8/L9 are limited to retiring the
borrowed participant/interface composition.

The post-M7 implementation now supersedes the original engine/participant
composition: `TimeFrame` owns reusable generic state, `Loader` and `Tracker`
act on it, `TrackerTraits` supplies the CPU/GPU kernels, and detector adapters
remain at the edges. Every milestone below is replay-gated; the only
intentional behavior change (native ITS refit) is explicitly fenced behind
its own decision (M5).

## Milestones

### M1 — TrackingEngine / TrackingParticipant contract

- **Goal / bounded API change**: commit the detector-neutral boundary as new
  headers plus a contract test: `TrackingEngine::executeEvent(TimeFrame&,
  schedule)` running `track()` over an ordered participant schedule against
  an event whose atomic loading has already been committed by the caller,
  returning a generic event result; `TrackingEngine::resetEvent(TimeFrame&,
  schedule)`, the same all-participant-reset-then-wipe-`TimeFrame`-once
  operation `executeEvent()` applies internally on a tracking failure, public
  so a caller whose own atomic load failed can reach it directly; the minimal
  `TrackingParticipant` interface (identity, owned surfaces, `track()`,
  participant-local event reset, publication export — deliberately no
  loading operation, see below); an explicit participant schedule as plan
  data (ADR 0007 decisions 3, 5, 6). Additive only — no existing header
  changes.
  - **Atomic loading vs. tracking execution**: these are two distinct,
    ordered steps, not one. Atomic event loading — every source's decoding
    and legacy backfill staged and committed together, or the live event
    left completely untouched — is a stronger whole-event transactional
    contract than a sequential per-participant `load()` inside
    `executeEvent()` could honor (an M1 draft of this contract had exactly
    that flaw; corrected before any M2 code depended on it). `TrackingEngine`
    therefore never loads anything: it only runs `track()` over an
    already-loaded event. Today's atomic loader is
    `MultiSourceTimeFrameLoader::loadITSAndMFT()`; M2 generalizes it into a
    participant-count-agnostic event loader that a caller runs, once, before
    calling `executeEvent()` — see M2 below.
- **Temporary bridge**: none (additive).
- **Acceptance/replay gate**: contract test pins the interface, including
  that a non-Success `track()` outcome or exception resets every scheduled
  participant and wipes the shared `TimeFrame` exactly once, and that
  `resetEvent()` applies that same contract directly; review confirms the
  interface names none of ITS/MFT/`NLayers`/`LegacyTrackerScratch`/
  source 0/1/output types, and no loading/event-origin operation.
- **Deletion/exit criterion**: M2 implements against the contract without
  amending it.
- **Dependency**: none.
- **Classification**: behavior-preserving cleanup.

### M2 — Engine seam and generic participant orchestration

- **Goal / bounded API change**: the concrete `TrackingEngine` runs `track()`
  over an ordered participant collection; a legacy CA participant
  implementation wraps today's `Tracker<NLayers>`/`TrackerTraits<NLayers>`/
  scratch/binding/sidecar composition unchanged; `CombinedTimeFrameCoordinator`
  builds two participants and delegates to the engine.
  `MultiSourceTimeFrameLoader` gains a participant-count-agnostic atomic
  event loader, with `loadITSAndMFT()`/`resetITSAndMFTEvent()` kept as
  forwarding bridges; the coordinator's `process()` runs that loader first,
  calls `TrackingEngine::executeEvent()` only once it reports success, and
  calls `TrackingEngine::resetEvent()` directly (never `executeEvent()`) on a
  load failure. No CA physics change; ITS+MFT is the first engine
  configuration.
- **Temporary bridge**: the coordinator's public API and the C4 DPL workflow
  stay untouched; the engine runs inside the coordinator.
- **Acceptance/replay gate**: existing coordinator and DPL contract tests
  green unchanged; combined, ITS-only, and MFT-only canonical replays
  byte-identical to the recorded baselines (see
  [Validation baseline](#validation-baseline)).
- **Deletion/exit criterion**: the coordinator holds no direct
  tracker/traits/scratch members — only the engine and its participants.
- **Dependency**: M1.
- **Classification**: behavior-preserving cleanup.

### M3 — Switch the C4 workflow to the engine; delete the coordinator

- **Goal / bounded API change**: `CombinedCATrackerSpec` constructs the engine
  and participants directly; `CombinedTimeFrameCoordinator.{h,cxx}` is
  deleted and its tests are ported to the engine (topology-copy identity,
  all-or-nothing failure reset, publication invalidation).
- **Temporary bridge**: none — this milestone removes M2's bridge.
- **Acceptance/replay gate**: combined replay byte-identical; DPL contract
  test green; the Gate 4 serial test selector (engine-ported equivalent)
  passes.
- **Deletion/exit criterion**: no reference to `CombinedTimeFrameCoordinator`
  remains (grep-verified); its contract coverage is demonstrably moved, not
  dropped.
- **Dependency**: M2.
- **Classification**: behavior-preserving cleanup.

### M4 — Contain TransitionPolicyTag

- **Goal / bounded API change**: the tag becomes nameable only in private
  legacy code (`src/`/detail headers). Remove only derivable stored/public
  exposure: the stored `SurfaceTransition::policyTag` field (derived from
  endpoint `SurfaceDescriptor::kind`s at layout build/validation, per ADR 0007
  decision 8), the public `TrackerTraits::processNeighbours<Tag>` surface, the
  `getPolicyBindingCount(TransitionPolicyTag)` accessor (replaced by a
  family-keyed accessor), and `TransitionPolicyDispatch.h` from the public
  include directory. No persisted replacement tag is introduced.
- **Temporary bridge**: the templated hot-loop machinery
  (`TransitionPolicyOperations`, grouping/dispatch) survives unchanged as
  private legacy implementation until M5 replaces it.
- **Acceptance/replay gate**: all canonical replays byte-identical; a guard
  test asserts no public header mentions `TransitionPolicyTag`.
- **Deletion/exit criterion**: tag confined to private legacy code; stored
  topology carries no tag field.
- **Dependency**: M3 (coordinator exposure); accessor renames may overlap M2.
- **Classification**: behavior-preserving cleanup.

### M5 — Descriptor-driven operations design and ITS refit A/B harness

- **Goal / bounded API change**: separately design (no production change) the
  single tracklet/cell/road/refit flow whose only family variation is a
  narrow operation boundary — projection, propagation target, rotation,
  material treatment, representation conversion — selected by inspecting the
  actual from/to `SurfaceDescriptor`s (ADR 0007 decision 10). Build the ITS
  legacy-vs-native refit A/B harness comparing the frozen
  `o2::its::track::refitTrackSeed` path against the native
  `driveRefitLeg`/`nativeRefitTrackCylinderCylinder` path on the canonical
  fixtures, reporting per-track parameter/chi2 deltas against declared
  tolerances.
- **Temporary bridge**: production refit stayed on the legacy path behind
  `DetectorTraits::refitSeed` until the separate activation decision below.
- **Acceptance/replay gate**: approved design note plus a recorded harness
  report; production replays untouched by construction.
- **Deletion/exit criterion**: a separate ADR records the native-refit
  activation decision (with its A/B evidence); only then may descriptor-driven
  implementation slices change production physics. Satisfied by [ADR
  0008](decisions/0008-native-refit-activation.md) (M5d): both
  `DetectorTraits::refitSeed` branches now call the shared `Propagator`
  (`Propagator.h`/`RefitDriver.h`); see that ADR for the exact legacy
  dependencies removed and the validation record.
- **Dependency**: M4.
- **Classification**: the design and harness (this milestone) are
  behavior-preserving; M5d's activation is the **intentional behavior
  change** ADR 0008 records.

### M6 — SurfaceTrackingScratch and legacy container removal

Split into a design/audit slice (M6a) plus six bounded, replay-gated
implementation slices (M6b–M6g), mirroring M5's a/b/c/d granularity. Full
field-by-field ownership audit, generic-workspace design, deletion order, and
per-slice acceptance gates are in [design note
0002](design/0002-m6-generic-workspace-migration.md); this section is the
summary anchor other milestones link against.

**Permanence rule** (design note 0002 §2, corrected): the ITS/MFT
*application-adapter role* survives M6 permanently (ADR 0007 decision 2).
`LegacyCATrackingParticipant` was renamed at M6f once its
`LegacyTrackerScratch`/`DetectorTraversalBinding` internals were gone. M6g
completed the remaining disposition: `ITSMFTLegacyParticipantSet` failed the
immutable-plan-builder test because it owned publication/timing state and
per-event reset/coordinator behavior, so it was deleted and its construction
logic was inlined into the combined DPL task. No `Legacy…`-named participant
set or coordinator remains in common-CA production code.

- **M6a — design/audit** (closed): read every `LegacyTrackerScratch<NLayers>`
  field, `DetectorTraversalBinding` responsibility, `LegacyCATrackingParticipant`
  responsibility, and `ITSMFTLegacyParticipantSet` responsibility in full;
  classified each as generic CA working state, `TimeFrame`-owned, adapter-private
  compatibility state, temporary migration bridge, or dead/redundant (design
  note 0002 §3). No production code change. Verdict: GO on the staged sequence
  below.
- **M6b — `SurfacePlanBinding`** (first bounded implementation slice): add the
  detector-neutral successor to `DetectorTraversalBinding` alongside it, unused
  by production — drops the two detector-switched lines
  (ITS/MFT allow-list; `expectedKind`/`expectedPolicy` derived by the caller
  instead of a `build()`-internal switch), keeps every other check unchanged
  (design note 0002 §3.2, §7, §9).
- **M6c — `SurfaceTrackingScratch`/`TrackSeed`**: add the non-templated,
  plan-sized workspace (sized from `ownedSurfaces().size()` and the topology's
  own `nTransitions`/`nCells`, per design note 0002 §4) and the fixed-capacity,
  GPU-portable `TrackSeed` replacement for `TrackSeedTpl<NLayers>`, unused by
  production.
- **M6d — wire MFT** onto `SurfaceTrackingScratch`/`SurfacePlanBinding` (MFT
  first — its refit is already fully normalized). Temporary bridge: an
  accessor seam both scratch types implement while ITS stays on the legacy
  type.
- **M6e — wire ITS** onto `SurfaceTrackingScratch`/`SurfacePlanBinding`;
  retire the legacy `mTracks`/`mTracksLabel` per-detector result staging so
  detector adapters build output from `TimeFrame::getCommonTracks()` plus the
  publication export alone (already populated in parallel today via
  `AcceptedTrackShadowPublisher` — this slice removes the now-redundant
  legacy-typed copy, not new publication logic). Writer-level output
  verification required, not only internal track-count/hash, since this
  slice changes an output-construction path.
- **M6f — delete** `LegacyTrackerScratch<NLayers>`, `DetectorTraversalBinding`,
  and the `NLayers`-templated `TrackSeedTpl` whole-track representation. Keep
  `SeedMetadataBase<N>` only for the still-live `CellSeed`; it is not a dead
  M6f bridge. Remove any other now-unreferenced `NLayers`-templated legacy CA
  container.
  **Rename** `LegacyCATrackingParticipant<NLayers>` (and its ITS/MFT aliases
  and extern-template instantiations) to a narrowly-scoped, non-`Legacy`
  ITS/MFT participant name — its own reason for the old name no longer applies
  once those deletions land (design note 0002 §3.3, §9). `ITSMFTLegacyParticipantSet`
  was disposed of at M6g by inlining its remaining application composition
  into the combined DPL task.
- **M6g — evaluate, relocate, and dispose of `ITSMFTLegacyParticipantSet`**:
  relocate its event-owned publication/timing bridge and per-event
  `clearPublicationSidecars()`/`configureRofTables()` calls into the
  combined DPL task's `process()`; the remaining layout, bindings,
  participants, and explicit schedule were then inlined into the task's
  application composition. No renamed holder was retained.
- **Acceptance/replay gate** (M6d–M6g, cumulative): per-detector and combined
  replays byte-identical to the M5d-era candidate baseline (or matching
  separately approved deltas); lifecycle/pool-destruction-order contracts
  re-pinned on the native scratch; memory/runtime changes recorded; a
  grep-guard test extended through M6f/M6g (design note 0002 §11) asserts zero
  remaining references in common-CA production include/src to
  `LegacyTrackerScratch`, `DetectorTraversalBinding`,
  `LegacyCATrackingParticipant`, `TrackSeedTpl`, and
  `ITSMFTLegacyParticipantSet`. The guard excludes only frozen legacy ITS code
  outside the common tree; no M6g participant-set exemption remains.
- **Deletion/exit criterion**: no production instantiation of
  `LegacyTrackerScratch<NLayers>` (grep-verified); native scratch has executed
  production traffic for both participants (ADR 0007 decision 9); no old M6f
  bridge name remains under `Detectors/ITSMFT/common/tracking/{include,src}`
  (grep-verified), including the deleted `ITSMFTLegacyParticipantSet`.
- **Dependency**: M5 implementation.
- **Classification**: M6a is a documentation-only audit; M6b/M6c are additive
  and behavior-preserving; M6d–M6g are behavior-preserving cleanup
  (replay-gated) — any residual output delta requires separate approval under
  M5's decision.
- **Flagged separately, not mandatory M6 scope, ranked by confidence** (design
  note 0002 §12): Rank 1, the already-dead
  `loadROFrameData()`/`resetROFrameData()`/`prepareROFrameData()` family (zero
  production callers today, independent of M6); Rank 2, the
  double-`TimeFrame::wipe()` on a single participant's recoverable-drop path
  (idempotent, not a correctness bug); Rank 3, verifying whether
  `mNTotalLowPtVertices` has any remaining reader (needs a wider grep than
  this audit performed).

**M6f validation record (2026-08-05)**: completed on
`codex/itsmft-m6f-retire-legacy-workspace`, based on the integrated M6e3 head
`ebcf6fc7608f2eb7562c248da008d799e3975515`, with
`daily-20260717-0700-local1` and the durable build specified by the gate. The
production guard scans only `Detectors/ITSMFT/common/tracking/{include,src}`;
it found no `LegacyTrackerScratch`, `DetectorTraversalBinding`,
`LegacyCATrackingParticipant`, `TrackSeedTpl`, `ScratchT`, or `BindingT` token.
The guard explicitly excludes frozen legacy ITS code outside that common tree
and permits the unchanged `ITSMFTLegacyParticipantSet` M6g adapter only.
`SeedMetadataBase<N>` remains solely because live `CellSeed` still uses it.

The common `SurfaceTrackingScratch` ROF-method audit found no production
definition or caller for `loadROFrameData()`, `resetROFrameData()`, or
`prepareROFrameData()`; that dead common-scratch family was removed. The
detector-specific raw-ROF APIs in ITS/MFT tracking, raw-ROF ownership, output
sidecars, and timing/publication compatibility state were not removed.

The durable affected-library/workflow/test rebuild completed, and serial
`ctest -L itsmft --output-on-failure -j1` executed all 95 registered tests:
95/95 passed, with no `Not Run` tests. On the 43-file fixture
`pp-20ev-run303000-seed20260716-daily20260717`, the checksum manifest passed
43/43 before and after replay. Standalone and combined replays both produced
ITS 212 tracks (`46913a67a7e2fe7462e29df0db264fa8`) and MFT 68 tracks
(`8106b08571ca593c6b76ff72b761a680`); the current standalone and combined
legs are field-identical, including the approved 2,992-value MFT float
projection. The M6f output was also compared with the accepted M6e3 parent:
all initialized writer leaves, ROFs, labels, cluster references, and the MFT
float-projected fields matched. The sole persisted difference was
`MFTTrack.mInvQPtSeed`, whose `TrackMFT` member is uninitialized in both the
parent and M6f production path (the parent contains run-dependent bytes such
as `0x60000002b`, while the M6f process contains `0x4`). It is not a defined
physics or GenericTrack value and was not rewritten or normalized in M6f, so
the writer/sidecar path remains unchanged; the discrepancy is recorded rather
than misreported as byte equality.

**M6g validation record (2026-08-05)**: completed on
`codex/itsmft-m6g-retire-participant-set`, based on accepted M6f commit
`65c0202650ad7d156cbade06b4a517483383d0bd`, with
`daily-20260717-0700-local1` and the same durable build. The affected libraries,
combined DPL workflow/executable, and changed tests rebuilt successfully. The
serial selector executed all 94 registered ITS/MFT tests: **94/94 passed**, with
no `Not Run` tests.

The 43-file fixture checksum manifest passed **43/43 before and after** the
replays. Standalone and combined common-CA replays both produced the accepted
candidates: ITS 212 / `46913a67a7e2fe7462e29df0db264fa8`, MFT 68 /
`8106b08571ca593c6b76ff72b761a680`. The combined ITS and MFT writer products
matched their standalone counterparts field-by-field; the MFT comparison covered
all 2,992 approved float-projected values with zero absolute/relative delta.
The combined output also matched the retained accepted pre-M6f M6e3 parent
artifact for all initialized ITS/MFT writer leaves, ROFs, labels, cluster
references, and sidecars. As in M6f, the only excluded parent difference was the
single uninitialized `MFTTrack.mInvQPtSeed` leaf; no defined writer or GenericTrack
content differed.

The M6g ownership guard and direct-composition tests cover load failure, dropped
TF, structural tracking failure, successful replacement, publication invalidation,
sidecar clearing, schedule order, and source-qualified exports. The common ROF
method audit found zero production definitions or callers for
`loadROFrameData()`, `resetROFrameData()`, or `prepareROFrameData()` in the common
tracking/workflow paths; those dead common methods remain deleted. The retained
exception is deliberately narrow: frozen ITS detector-specific `TimeFrame` raw-ROF
methods and MFT raw-ROF `IOUtils`/workflow ownership remain live outside common CA,
as do output sidecars and workflow-owned timing/publication state. No raw-ROF
workflow ownership or compatibility state was removed. `nvcc` and `hipcc` were
absent from the pinned environment, so no device build or GPU replay was claimed.

### M7a — runtime-plan `Tracker`/`TrackerTraits` design and ownership audit

**Status: closed as a documentation-only slice (2026-08-05).** The complete
remaining `NLayers` audit, ownership classification, runtime-plan contract,
cylinder/disk operation boundary, staged implementation plan, and ranked
deletion inventory are in [design note 0003](design/0003-m7-runtime-plan-tracker-migration.md).

The design selects the existing `DetectorLayoutSet`/`DetectorLayoutView`,
`SurfacePlanBinding`, and `SurfaceTrackingScratch` composition as the runtime
plan. It introduces no wrapper around `Tracker<NLayers>` and no new dispatch
taxonomy. The first implementation proof is the deletion of the
`Tracker<NLayers>`/`TrackerTraits<NLayers>` declarations, their 7/10 core
instantiations, and the `CATracker` aliases; adapter-only templates may remain
temporarily under their stated exit criteria.

M7a changed documentation only. It did not run replay, CTest, or a device
build. The next bounded implementation slice is M7b: make the existing plan
and workspace counts the only common-core loop/capacity authority, with the
full R replay/dependency gate specified in design note 0003. Frozen
`Detectors/ITSMFT/ITS/tracking` code remains outside the narrow common-CA
guard, and no production `ITSMFTLegacyParticipantSet` exemption is carried
forward from M6g.

### M7b — make the existing runtime plan the only count authority

**Status: complete (2026-08-05).** The production migration, focused guard,
and validation record are in
[design note 0004](design/0004-m7b-runtime-count-authority.md). The existing
`DetectorLayoutSet`/`SurfacePlanBinding`/`SurfaceTrackingScratch` composition
is now the source of ordered traversal, active-surface, transition/cell,
measurement, seed-mask, and index-table-prefix counts. No runtime-plan wrapper
or parallel implementation was added.

The M7b source guard classifies every remaining common production `NLayers`
use as fixed device ABI/capacity, temporary private operation implementation,
adapter edge, or explicitly deferred M7d Tracker/TrackerTraits compatibility.
It rejects new
plan-sized host arrays and shared-core `.NLayers` count authorities. The
sparse/non-identity plan test proves binding order, source-qualified
measurements, compact transition/cell slots, and positional `TrackSeed`
masks. `Tracker`/`TrackerTraits` de-templating remains the next M7d slice; no
implementation from that slice is included here.

The durable build reused package `daily-20260717-0700-local1`; the complete
serial selector executed all 95 registered ITS/MFT tests and passed 95/95 with
no `Not Run` entries. The 43-file fixture checksum manifest passed 43/43
before and after replay. Standalone and combined common-CA replays produced
ITS 212 / `46913a67a7e2fe7462e29df0db264fa8` and MFT 68 /
`8106b08571ca593c6b76ff72b761a680`; the combined legs matched their
standalone products field-by-field. The MFT comparison covered 2,992
float-projected values with zero maximum absolute/relative delta. Parent
writer comparison preserved all initialized content and excluded only the
known undefined `MFTTrack.mInvQPtSeed` leaf. No GPU result is claimed because
the pinned Darwin host has no `nvcc`/`hipcc` toolchain or device.

### M7c — remove layer topology and ROF compatibility from the common path

**Status: complete (2026-08-05).** The production migration, ownership guard,
focused sparse/ROF tests, and validation record are in [design note
0005](design/0005-m7c-runtime-topology-rof.md).

Common traversal now consumes the existing `DetectorLayoutView` sparse
topology and `SurfacePlanBinding` IDs/compact positions. `TrackingTopology.h`,
the old topology test, scratch topology members, templated scratch dispatchers,
`checkSupportedNLayers()`, and the `IndexTableUtils<N>` alias are gone.
`TrackSeed::SurfaceMask` remains in compact binding-position space.

`SurfaceTrackingScratch` receives one non-owning `RuntimeROFViews` value. The
frozen fixed-capacity ROF builders are explicitly confined to the two
application-edge seams (`TrackingInterface` and
`SurfacePlanTrackingParticipant`); raw ROFs, publication clocks, validity,
sidecars, and event reset sequencing remain workflow-owned. Frozen legacy ITS
sources are outside the common guard and remain unchanged.

No `Tracker`/`TrackerTraits` de-templating is included. The exact M7c build,
serial suite, source guards, and replay evidence are recorded in design note
0005 and in the final validation record for this plan.

**M7c validation record (2026-08-05)**: the requested durable build completed
for the affected libraries, adapters, combined workflow, and tests. The final
serial selector executed all 95 registered ITS/MFT tests and passed 95/95 with
no `Not Run` entries. The fixture checksum manifest passed 43/43 before and
after replay. Fresh standalone and combined replays produced ITS 212 /
`46913a67a7e2fe7462e29df0db264fa8` and MFT 68 /
`8106b08571ca593c6b76ff72b761a680`; combined legs matched their standalone
writer products field-by-field. M7c standalone and combined writer products
matched the M7b parent initialized content; only the known undefined
`MFTTrack.mInvQPtSeed` leaf was excluded, with 2,992 MFT float-projected
values compared at zero maximum absolute/relative delta. No GPU result is
claimed because the pinned environment has no `nvcc`, `hipcc`, `nvidia-smi`, or
`rocminfo`.

### M7d — make Tracker and TrackerTraits non-templated

**Status: complete (2026-08-05).** The production migration, focused tests,
source/dependency guards, and validation record are in [design note
0006](design/0006-m7d-nontemplated-tracker-core.md). `TrackerTraits` and
`Tracker` are now single non-templated runtime-plan classes. Their loop and
workspace authorities remain the M7b/M7c plan, binding, sparse-topology, and
scratch composition; no hidden 7/10 implementation or compatibility alias
remains in common production code.

`SurfacePlanTrackingParticipant<7/10>` and
`ITSMFTTrackingInterface<7/10>` remain application-edge templates with plain
core members. Typed refit/output conversion and sidecar operations cross the
core only through a call-scoped operation reference and remain the exact
narrow M7e boundary. The common core still contains no detector identity,
source convention, DPL/workflow/writer dependency, or public policy taxonomy.

The M7d implementation is behavior-preserving by construction; it does not
alter Propagator/refit arithmetic, covariance policy, CA choices, holes,
ordering, output semantics, raw-ROF ownership, or workflow defaults.

### M7e — move typed refit/output hooks to application adapters

**Status: complete (2026-08-05).** The production migration, adapter tests,
ownership guards, and validation record are in [design note 0007](design/0007-m7e-adapter-refit-output.md).
`DetectorTraits<NLayers>`, `CATrackType<NLayers>`,
`LayerMeasurementSpans<NLayers>`, `AcceptedTrackShadowPublisher`, and the
typed MFT helper hook are deleted from generic-core ownership. One narrow
`TrackingOperationAdapter` remained at this milestone for generic seed refit.
The subsequent [fine-comb cleanup](validation/fine-comb-header-boundaries.md)
replaced it with one typed function pointer and direct ITS/MFT refit functions;
no adapter object or detector-ID refit dispatcher remains. ITS/MFT edges consume the same
`TrackSeed`/`GenericTrack` result path, with ITS shared-status staging owned by
the ITS publication adapter.

The durable build passed all 97 registered serial ITS/MFT tests. Fixture
checksums passed 43/43 before and after. Standalone and combined replay legs
produced ITS 212 / `46913a67a7e2fe7462e29df0db264fa8` and MFT 68 /
`8106b08571ca593c6b76ff72b761a680`; combined legs matched standalone products
and the M7d parent matched in all initialized content. Only the undefined
`MFTTrack.mInvQPtSeed` byte artifact remains excluded. No GPU result is
claimed because the pinned environment has no CUDA/HIP/device tools.

### M7f — delete redundant compile-time bridges and close the runtime-core guard

**Status: complete (2026-08-06).** The production deletion, focused test
migrations, authoritative source/dependency guard, and validation record are
in [design note 0008](design/0008-m7f-final-runtime-core-cleanup.md). The
common road-start loop now derives its extent from
`SurfaceTrackingScratch::getNOwnedSurfaces()`, and sparse transition/cell
endpoints resolve through `SurfacePlanBinding` positions. The persistent
SurfaceId-to-slot bridge, `IndexTableUtilsN`, redundant `<NLayers>` refit
helpers, and dead typed native-output export helper are deleted.

The guard classifies every live layer-count token in common production
include/src by exact file and role. Only named adapter-edge configuration,
static descriptor, participant/interface, publication, and frozen-compatibility
uses remain; fixed device capacities retain their maximum storage with runtime
counts/masks. The frozen legacy ITS tree is excluded only by its explicit
`Detectors/ITSMFT/ITS/tracking` path and is unchanged.

M7f is behavior-preserving by construction. The full durable build passed, the
serial ITS/MFT selector passed 98/98 with no `Not Run` entries, fixture checksums
passed 43/43 before and after, and standalone/combined replay anchors remained
ITS 212 / `46913a67a7e2fe7462e29df0db264fa8` and MFT 68 /
`8106b08571ca593c6b76ff72b761a680`. Combined legs matched standalone writer
content and M7e initialized output matched, excluding only the known undefined
`MFTTrack.mInvQPtSeed` byte. No GPU result is claimed because the pinned host
has no CUDA/HIP/device tools.

## Pass-owned traversal schedules

`SurfacePlanBinding` has been deleted. `Tracker` now derives each pass's
active surfaces, selected edges/cells, compact slots, and road schedules into
`TraversalWorkspace`; `SurfaceGraph` remains immutable static configuration.
The focused migration record is [workspace-owned traversal plan](validation/workspace-owned-traversal-plan.md).

## Not safe to delete yet

| Artifact | Why it must stay | Removal gate |
|---|---|---|
| `LegacyTrackerScratch<NLayers>` (and its `NLayers` container family) | Retired in M6f after both participants executed production traffic on `SurfaceTrackingScratch`; no longer safe-to-delete work remains | **Deleted M6f**; see [design note 0002](design/0002-m6-generic-workspace-migration.md) §9–§10 and the M6f validation record below |
| `DetectorTraversalBinding` (`include/ITSMFTTracking/detail/DetectorTraversalBinding.h`) | Its plan-sizing/slot-assignment role is now carried by `SurfacePlanBinding` for both participants | **Deleted M6f** after M6d/M6e wiring and replay gates |
| `loadITSAndMFT()`/`resetITSAndMFTEvent()` fixed source-0/1 mapping | Current combined load/reset contract of the C4 workflow | M2 generalizes, M3 deletes the fixed-position entry points |
| `TransitionPolicyTag` machinery (dispatch, grouping, templated operations) | Only existing hot-loop implementation of the CA stages | contained at M4, replaced by M5 implementation |
| Policy/legacy compatibility code (`kDroppedTimeFrameResult` sentinel, `mLayerMaterial`/`LegacyMaterialMismatch`, `DiskDiskReferenceCoordinateView`, `passesCellRoadPrecut<DiskDisk>`) | Pins byte-identical replay parity against the frozen legacy implementations | respective M4–M6 slices, each under its replay gate |
| Output sidecars (`ITSSharedClusterCompatibility`, `MFTPublicationCompatibility`) | Legacy output conversion still requires per-detector compatibility state | M6, when adapters convert from `GenericTrack` alone |
| `SurfacePlanTrackingParticipant<NLayers>` | Narrow plan-driven ITS/MFT participant adapter; it owns detector/output sidecars but no retired workspace/binding bridge or combined event loop | **M7f:** retained only as the named application-edge compatibility seam; reduce it after workflow-edge ROF/configuration consumers are independently runtime-view based |
| `ITSMFTLegacyParticipantSet` | Coordinator-shaped holder of combined application construction and event-owned publication/reset state | **Deleted M6g**; construction is inlined into the combined DPL task and publication/timing ownership remains workflow-local |
| Common `TrackingTopology<NLayers>` | The common layer-indexed topology duplicated the sparse plan topology | **Deleted M7c**; frozen legacy ITS topology is outside the common-CA guard and remains unchanged |
| Common `ROF*Table<NLayers>` scratch ownership | The core now receives non-owning runtime timing/overlap/mask views | **Deleted M7c**; fixed-capacity builders remain only at explicitly named adapter edges until their own template seam is reduced |
| Call-scoped operation seam | The generic core still needs one operation-local refit/completion boundary while typed MFT refit and sidecars are adapter-owned | Deferred post-M7 ownership decision; do not remove without a stable generic accepted-result owner and replay-gated sidecar migration, see [design note 0007](design/0007-m7e-adapter-refit-output.md) |

## Validation baseline

Durable records this plan's replay gates compare against — linked, not
copied:

- **Gate 4 acceptance evidence** (combined workflow, coordinator tests,
  serial 8/8 rerun): the Gate 4 row of
  [AgentCoordination.md, Gate status](AgentCoordination.md#14-gate-status);
  workflow replay artifacts under
  `/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/gate4-c4-combined-workflow/`.
- **Canonical single-detector replay anchors**: ITS 203 tracks, state hash
  `ee7f7c794d60f2362fd2564258b7887e`; MFT 70 tracks, hash
  `24737e73b7146bf3bd35a90a2517c527` — as recorded in
  [ADR 0006](decisions/0006-provider-aggregation-cleanup.md) (validation
  gates) and the artifact records it links.
- **MFT float-projected output contract**: MFT's float-based
  `TrackingFrameInfo`-projected publication path is the approved output
  contract carried unchanged through the combined staging (Gate 4 row above);
  byte-identical evidence under
  `/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/stageb-mft-normalized-refit/`.

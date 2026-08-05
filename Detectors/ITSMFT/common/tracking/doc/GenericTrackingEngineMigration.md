# Generic Tracking-Engine Migration Plan (post-Gate-4)

Status: Living plan
Decision anchor: [ADR 0007](decisions/0007-generic-tracking-engine-boundary.md)
Architecture: [Architecture.md](Architecture.md)
Status/decision log: [AgentCoordination.md](AgentCoordination.md)

This plan turns the accepted Gate 4 implementation into the ADR 0007 end
state: one concrete `TrackingEngine::executeEvent()` over ordered
`TrackingParticipant`s, permanent `TimeFrame`/`CommonTrack` ownership,
detector adapters at the edges, and one algorithmic flow over
`SurfaceDescriptor`/`SurfaceKinematicState`. Every milestone below is
replay-gated; the only intentional behavior change (native ITS refit) is
explicitly fenced behind its own decision (M5).

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
  (`Propagator.h`/`NativeRefitDriver.h`); see that ADR for the exact legacy
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
physics or CommonTrack value and was not rewritten or normalized in M6f, so
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
single uninitialized `MFTTrack.mInvQPtSeed` leaf; no defined writer or CommonTrack
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

## Not safe to delete yet

| Artifact | Why it must stay | Removal gate |
|---|---|---|
| `LegacyTrackerScratch<NLayers>` (and its `NLayers` container family) | Retired in M6f after both participants executed production traffic on `SurfaceTrackingScratch`; no longer safe-to-delete work remains | **Deleted M6f**; see [design note 0002](design/0002-m6-generic-workspace-migration.md) §9–§10 and the M6f validation record below |
| `DetectorTraversalBinding` (`include/ITSMFTTracking/detail/DetectorTraversalBinding.h`) | Its plan-sizing/slot-assignment role is now carried by `SurfacePlanBinding` for both participants | **Deleted M6f** after M6d/M6e wiring and replay gates |
| `loadITSAndMFT()`/`resetITSAndMFTEvent()` fixed source-0/1 mapping | Current combined load/reset contract of the C4 workflow | M2 generalizes, M3 deletes the fixed-position entry points |
| `TransitionPolicyTag` machinery (dispatch, grouping, templated operations) | Only existing hot-loop implementation of the CA stages | contained at M4, replaced by M5 implementation |
| Policy/legacy compatibility code (`kDroppedTimeFrameResult` sentinel, `mLayerMaterial`/`LegacyMaterialMismatch`, `mSurfaceToLegacyLayer`, `DiskDiskReferenceCoordinateView`, `passesCellRoadPrecut<DiskDisk>`) | Pins byte-identical replay parity against the frozen legacy implementations | respective M4–M6 slices, each under its replay gate |
| Output sidecars (`ITSSharedClusterCompatibility`, `MFTPublicationCompatibility`) | Legacy output conversion still requires per-detector compatibility state | M6, when adapters convert from `CommonTrack` alone |
| `SurfacePlanTrackingParticipant<NLayers>` | Narrow plan-driven ITS/MFT participant adapter; it owns detector/output sidecars but no retired workspace/binding bridge | Remains through M6 and owns no event-loop coordination |
| `ITSMFTLegacyParticipantSet` | Coordinator-shaped holder of combined application construction and event-owned publication/reset state | **Deleted M6g**; construction is inlined into the combined DPL task and publication/timing ownership remains workflow-local |

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

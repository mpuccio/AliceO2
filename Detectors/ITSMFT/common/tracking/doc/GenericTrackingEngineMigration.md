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
  headers plus a contract test: `TrackingEngine::executeEvent(TimeFrame&, ...)`
  returning a generic event result; the minimal `TrackingParticipant`
  interface (identity, owned surfaces, load, track, event reset, publication
  export); an explicit participant schedule as plan data (ADR 0007 decisions
  3, 5, 6). Additive only — no existing header changes.
- **Temporary bridge**: none (additive).
- **Acceptance/replay gate**: contract test pins the interface; review
  confirms the interface names none of ITS/MFT/`NLayers`/
  `LegacyTrackerScratch`/source 0/1/output types.
- **Deletion/exit criterion**: M2 implements against the contract without
  amending it.
- **Dependency**: none.
- **Classification**: behavior-preserving cleanup.

### M2 — Engine seam and generic participant orchestration

- **Goal / bounded API change**: the concrete `TrackingEngine` executes an
  ordered participant collection; a legacy CA participant implementation wraps
  today's `Tracker<NLayers>`/`TrackerTraits<NLayers>`/scratch/binding/sidecar
  composition unchanged; `CombinedTimeFrameCoordinator` builds two
  participants and delegates `process()` to the engine.
  `MultiSourceTimeFrameLoader` gains a participant-count-agnostic load/reset,
  with `loadITSAndMFT()`/`resetITSAndMFTEvent()` kept as forwarding bridges.
  No CA physics change; ITS+MFT is the first engine configuration.
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
- **Temporary bridge**: production refit stays on the legacy path behind
  `DetectorTraits::refitSeed` until a separate activation decision.
- **Acceptance/replay gate**: approved design note plus a recorded harness
  report; production replays untouched by construction.
- **Deletion/exit criterion**: a separate ADR records the native-refit
  activation decision (with its A/B evidence); only then may descriptor-driven
  implementation slices change production physics.
- **Dependency**: M4.
- **Classification**: the design and harness are behavior-preserving;
  activating native refit (or any unified-flow output delta) is an
  **intentional behavior change** requiring that separate decision.

### M6 — SurfaceTrackingScratch and legacy container removal

- **Goal / bounded API change**: introduce the native, runtime-sized,
  non-templated `SurfaceTrackingScratch` (sized from the plan/binding compact
  slot counts; transients keyed by `SurfaceId`/compact transition and cell
  slots; results published to `TimeFrame::getCommonTracks()`). Participants
  migrate one at a time (MFT first — its refit is already fully normalized);
  legacy output conversion moves entirely into detector adapters reading
  `CommonTrack` plus the publication export. Then delete
  `LegacyTrackerScratch<NLayers>` and the remaining `NLayers`-templated legacy
  CA containers in the core.
- **Temporary bridge**: an accessor seam both scratch types implement during
  the per-participant switch.
- **Acceptance/replay gate**: per-detector and combined replays byte-identical
  (or matching separately approved M5 deltas, whichever decision stands);
  lifecycle/pool-destruction-order contracts re-pinned on the native scratch;
  memory/runtime changes recorded.
- **Deletion/exit criterion**: no production instantiation of
  `LegacyTrackerScratch<NLayers>` (grep-verified); native scratch has executed
  production traffic (ADR 0007 decision 9).
- **Dependency**: M5 implementation.
- **Classification**: behavior-preserving cleanup (replay-gated); any residual
  output delta requires separate approval under M5's decision.

## Not safe to delete yet

| Artifact | Why it must stay | Removal gate |
|---|---|---|
| `LegacyTrackerScratch<NLayers>` (and its `NLayers` container family) | Sole production owner of CA transients and legacy result staging | M6, after `SurfaceTrackingScratch` executes production traffic |
| `DetectorTraversalBinding` legacy-slot translation | Hot loops still index compact legacy scratch slots by global id | M6 (the plan-sizing role may survive; the legacy-slot translation is the temporary part) |
| `loadITSAndMFT()`/`resetITSAndMFTEvent()` fixed source-0/1 mapping | Current combined load/reset contract of the C4 workflow | M2 generalizes, M3 deletes the fixed-position entry points |
| `TransitionPolicyTag` machinery (dispatch, grouping, templated operations) | Only existing hot-loop implementation of the CA stages | contained at M4, replaced by M5 implementation |
| Policy/legacy compatibility code (`kDroppedTimeFrameResult` sentinel, `mLayerMaterial`/`LegacyMaterialMismatch`, `mSurfaceToLegacyLayer`, `DiskDiskReferenceCoordinateView`, `passesCellRoadPrecut<DiskDisk>`) | Pins byte-identical replay parity against the frozen legacy implementations | respective M4–M6 slices, each under its replay gate |
| Output sidecars (`ITSSharedClusterCompatibility`, `MFTPublicationCompatibility`) | Legacy output conversion still requires per-detector compatibility state | M6, when adapters convert from `CommonTrack` alone |

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

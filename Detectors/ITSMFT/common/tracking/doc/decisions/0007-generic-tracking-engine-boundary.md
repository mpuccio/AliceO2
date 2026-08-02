# ADR 0007: generic tracking-engine boundary (post-Gate-4 target)

Status: Accepted
Date: 2026-08-02

## Context

Gate 4 is accepted with complete evidence (see
[AgentCoordination.md, Gate status](../AgentCoordination.md#14-gate-status)):
ITS and MFT run as disconnected applications of one common CA tracker over one
shared `TimeFrame`, reached in production through the opt-in
`o2-itsmft-combined-ca-tracker-workflow`. That acceptance was deliberately
built on migration artifacts, not on the long-term architecture:
`CombinedTimeFrameCoordinator` hardcodes the two-detector shape,
`TransitionPolicyTag` still appears in public core headers and stored topology
data, and `LegacyTrackerScratch<NLayers>` still owns every per-layer CA
container. The long-term target is a detector-neutral solenoidal-field
tracking core that future adapters (for example ALICE 3 trackers) can supply
with plans and inputs without the core learning their identity.

This ADR fixes the end-state boundary those artifacts migrate toward. It
refines — without reversing — D007/D010 in the
[decision log](../AgentCoordination.md#13-decision-log) and
[Architecture.md](../Architecture.md) Sec 10.1: the tag-dispatch *mechanics*
accepted there remain valid hot-loop implementation, but the tag itself is now
classified as temporary. The milestone plan implementing this ADR is
[GenericTrackingEngineMigration.md](../GenericTrackingEngineMigration.md).

## Decision

1. **`TimeFrame` and `CommonTrack` are permanent.** The non-templated
   `TimeFrame` is the detector-neutral per-event data owner; `CommonTrack` (with
   `TrackClusterReference`) is the generic algorithm result. Neither is a
   migration artifact.

2. **Detectors are adapters.** ITS, MFT, and future detectors such as ALICE 3
   supply decoding, timing/ROF policy, static surface plans, DPL origins,
   writers, and legacy output conversion. The core never permanently knows
   ITS, MFT, layer counts 7/10, source 0/1, or current output track types.

3. **One concrete `TrackingEngine` with `executeEvent()`.** The engine executes
   one event over an ordered collection of participants and applies the
   whole-event all-or-nothing failure contract proven in Gate 4. There is
   exactly one such class; a second public "TrackingExecutor" abstraction (an
   engine interface plus a concrete executor) must not be introduced.

4. **`CombinedTimeFrameCoordinator` is a temporary Gate 4 wrapper.** It is
   removed once the combined workflow runs on the engine and the parity/replay
   gates pass.

5. **`TrackingParticipant` is the detector-leg boundary.** Its minimal
   interface (identity, owned surfaces, load, track, event reset, publication
   export) must not mention ITS, MFT, `NLayers`, `LegacyTrackerScratch`,
   source 0/1, or current output types. All of those live inside concrete
   participant implementations or adapters.

6. **Participant execution order is explicit plan/schedule data.** The engine
   consumes an ordered schedule supplied with the plan. Ordering is *not* an
   inherent consequence of static-catalog concatenation order, even though the
   first ITS+MFT configuration happens to coincide with it.

7. **`TransitionPolicyTag` is a temporary legacy implementation detail.** It is
   not a future public/core abstraction; no new public or core API may expose
   it.

8. **No persisted tag replacement.** In the current same-family-only topology,
   tag information is fully derivable from the endpoint
   `SurfaceDescriptor::kind`s. For future mixed surfaces, operations inspect
   the actual from/to `SurfaceDescriptor`s directly at the operation boundary;
   a persisted "SurfaceKindPair" (or equivalent) tag must not be introduced.

9. **`LegacyTrackerScratch<NLayers>` is a temporary bridge.** It cannot be
   removed before a native, runtime-sized `SurfaceTrackingScratch` replacement
   executes production traffic.

10. **One algorithmic flow for cylinders and disks.** Cylinder and disk
    orchestration must converge to one tracklet/cell/road/refit flow over
    `SurfaceDescriptor` and `SurfaceKinematicState`. Family differences belong
    only in narrow, explicit surface-operation boundaries (projection,
    propagation target, rotation, material treatment, representation
    conversion) — never in duplicated orchestration or per-detector pipelines.
    Mixed-surface scenarios follow from that model.

11. **Intentional behavior changes require separate decisions.** Any
    physics/output behavior change — in particular activating the native
    cylinder-cylinder refit in place of the frozen legacy ITS refit — requires
    its own recorded decision backed by A/B validation. Migration slices under
    this ADR are otherwise replay-gated and behavior-preserving.

## Non-goals

No big-bang rewrite; no tracking-physics change inside this migration charter;
no GPU migration decisions; no change to workflow defaults or user-facing
configuration keys.

## Consequences

- The engine seam can land immediately inside the coordinator without touching
  CA physics, then the coordinator can be hollowed and deleted.
- Containing the tag (decisions 7–8) removes derivable stored/public exposure
  now, while the templated hot-loop machinery remains private until the
  descriptor-driven operation boundary (decision 10) replaces it.
- The scratch replacement (decision 9) is the last milestone and the deletion
  gate for the remaining `NLayers`-templated core containers.
- Milestones, bridges, acceptance gates, and the not-safe-to-delete inventory
  are maintained in
  [GenericTrackingEngineMigration.md](../GenericTrackingEngineMigration.md),
  not duplicated here.

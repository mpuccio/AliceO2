# Design note 0006: M7d non-templated runtime-plan tracker core

Status: M7d implementation and validation record
Date: 2026-08-05
Branch: `codex/itsmft-m7d-nontemplated-tracker-core`
Base: `9649b71bc3` (integrated M7c equivalent)
Package: `daily-20260717-0700-local1`
Durable build: `/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement`

Companion plan: [GenericTrackingEngineMigration.md](../GenericTrackingEngineMigration.md)
Preceding design: [M7 runtime-plan migration](0003-m7-runtime-plan-tracker-migration.md)
Preceding slice: [M7c runtime topology and ROF ownership](0005-m7c-runtime-topology-rof.md)
Architecture decision: [ADR 0007](../decisions/0007-generic-tracking-engine-boundary.md)

This note records the M7d structural migration. It makes the shared
`TrackerTraits` and `Tracker` implementations non-templated while preserving
the existing runtime plan, sparse topology, ROF views, CA policy, native
refit, failure contract, and output adapters. It does not begin M7e typed
refit/output relocation or M7f final cleanup.

## 1. Decision and ownership result

There is now exactly one common orchestration implementation:

| Responsibility | M7d owner | Authority |
|---|---|---|
| Shared CA orchestration | non-templated `TrackerTraits` | adopted `DetectorLayoutView`, `SurfacePlanBinding`, and `SurfaceTrackingScratch` |
| Iteration execution and failure classification | non-templated `Tracker` | runtime scratch/plan plus the caller-owned operation invocation |
| Whole-track device-safe value | `TrackSeed` | fixed `MaxLayoutSurfaces` capacity and positional `LayerMask` |
| Detector/refit/output conversion | `SurfacePlanTrackingParticipant<N>` and `ITSMFTTrackingInterface<N>` | application-owned `DetectorTraits<N>` and publication sidecars |
| Event lifecycle and raw ROFs | DPL workflow/application task | existing loader, reset, timing, and publication ownership |

`Tracker` and `TrackerTraits` no longer have a layer-count template argument,
and the common library no longer emits 7- or 10-surface core instantiations.
`CATracker`, `TrackerITS`, `TrackerMFT`, `CATrackerITS`, and `CATrackerMFT`
compatibility aliases were deleted. The first proof is therefore a deletion
of the old core type vocabulary, not a wrapper around hidden 7/10 instances.

The application participants remain templated only because they still own
detector-specific configuration, fixed-capacity ROF builders, typed refit
conversion, and output sidecars. Each participant contains plain
`TrackerTraits` and `Tracker` objects. Neither participant owns a combined
schedule or event lifecycle.

## 2. Call-scoped M7e boundary

M7d leaves a narrow operation seam for the typed work that is intentionally
not moved until M7e. `TrackingOperationAdapter` is passed by reference to
`Tracker::clustersToTracks()` and down the road/accept/mark call for one
tracking invocation. The non-templated core does not store the adapter, own a
typed output track, or retain publication state. Failure cleanup calls the
same invocation adapter's `clearPublicationState()` before the existing
scratch/`TimeFrame` reset.

The seam is deliberately call-scoped and has no alternative implementation:

- refit converts a `TrackSeed` into a detector-neutral `TrackingCandidate`;
- acceptance publishes the candidate's `CommonTrack` shadow through the
  adapter-owned output path;
- polarity and final shared-cluster sealing remain operation-local hooks; and
- typed conversions are centralized in
  `DetectorTrackingOperationAdapterSupport.h`, at the application edge.

M7e must remove this temporary seam rather than grow it. Its bounded target is
to move `DetectorTraits<N>` refit/export calls, typed candidate conversion,
`DetectorTrackingOperationAdapterSupport.h`, and sidecar sealing fully into
the two application adapters. The generic body may retain only runtime
measurement spans, `SurfaceKinematicState`, `TrackSeed`, and `CommonTrack`
data. No new coordinator, type-erased core owner, detector switch, or second
tracking implementation is permitted.

The narrow remaining leaf exception is the existing descriptor-driven MFT
forward-state helper used by a disk operation. M7e may relocate that helper
if the typed refit/output extraction requires it, but it must remain a leaf
operation selected from actual descriptors/states rather than a detector or
layer-count dispatch.

## 3. Runtime-plan invariants

The M7b/M7c authorities remain unchanged:

- ordered traversal uses the validated `SurfacePlanBinding` positions;
- active-surface and temporary workspace extents use
  `SurfaceTrackingScratch`/sparse-topology counts;
- transition and cell slots use the binding's compact mappings;
- source-qualified measurements resolve through the binding and normalized
  `TimeFrame` view; and
- `TrackSeed::LayerMask` records compact plan positions, not numeric
  `SurfaceId` or ITS/MFT layer identity.

`TrackingParameters::NLayers` remains only an adapter-edge validation against
the adopted active-surface count. Fixed device capacities are not reduced:
`TrackSeed` retains `MaxLayoutSurfaces`, `CellSeed` retains its three
measurement slots, and fixed GPU index-table storage retains its maximum
capacity with a runtime populated prefix.

Cylinder/disk differences continue to be selected by the existing private
descriptor-driven operation binding. M7d does not expose `TransitionPolicyTag`,
introduce `SurfaceKind`/`SurfaceKindPair` routing, or add detector identity,
source conventions, DPL, workflow, writer, or typed output dependencies to
the core headers.

## 4. Tests and guards

The M7d guard scans common tracking production include/source files after
removing comments and preprocessor lines. It rejects declarations,
definitions, and instantiations containing `Tracker<`, `TrackerTraits<`,
`CATracker`, `TrackerITS`, or `TrackerMFT`. It separately checks the
`Tracker.h` and `TrackerTraits.h` source text for `DetID`, DPL/workflow/writer
tokens, ITS/MFT identity, and `SurfaceKindPair`.

The guard has no frozen-legacy exemption inside `common/tracking`; frozen ITS
sources are outside its scanned tree. The only intentional M7d boundary is
the adapter-only participant/interface templates and their named support
header, which are covered by the existing M7b classification guard rather
than treated as generic core.

The focused runtime test constructs a four-surface non-identity plan with
ordered positions `{1, 4, 7, 10}` and executes the single non-templated core.
It checks plan order, compact ownership, active `TrackSeed` mask positions,
and sparse-policy traversal. The migrated failure contract passes the
call-scoped adapter through successful, recoverable-drop, and structural
paths, and the combined composition test continues to cover both adapter
participants, schedule order, and source-qualified exports.

## 5. Ranked deletion and simplification inventory

### Safe in M7d

1. Delete the common `Tracker<NLayers>` and `TrackerTraits<NLayers>`
   declarations, definitions, aliases, and explicit 7/10 instantiations.
2. Delete the `CATracker`/`TrackerITS`/`TrackerMFT` compatibility aliases.
3. Keep one generic candidate representation in the shared orchestration and
   pass typed operation work through the call-scoped edge seam.
4. Remove the layer-count-to-`SurfaceKind` bridge, since policy/representation
   selection already comes from the validated plan and actual surface states.

### Blocked until M7e or later

1. Remove `DetectorTraits<N>` from refit/output construction and delete
   `DetectorTrackingOperationAdapterSupport.h` after its two adapter callers
   are converted.
2. Remove the `TrackingOperationAdapter` seam and any generic accepted
   candidate staging that exists only to service typed publication.
3. Reduce `SurfacePlanTrackingParticipant<N>` and
   `ITSMFTTrackingInterface<N>` templates only after detector configuration,
   fixed ROF builders, typed output, and sidecar ownership have independent
   parity coverage.
4. Audit the remaining MFT forward-state helper include in the shared source
   and move it to the narrow operation leaf if M7e's extraction permits.

### Deferred physics/algorithm choices

These are deliberately not M7d work: changing Propagator equations or
covariance sanitization; changing CA cuts, holes, policy grouping, candidate
ordering, timing arithmetic, or source qualification; changing fixed device
capacities; making currently rejected mixed cylinder/disk plans legal; or
deleting a legacy numerical path based only on code similarity. Each requires
its own physics/algorithm decision and replay evidence.

## 6. Validation record

The durable build and package named above were reused. The affected common
library, adapter/workflow, and test targets rebuilt after the migration. The
M7b runtime-count guard, M7d non-templated-core guard, failure-contract test,
and combined-composition test passed after the final call-scoped seam change.

The complete serial `itsmft` suite, fixture checksum replay, standalone and
combined candidate replays, initialized writer comparison against M7c, and
device-toolchain check completed as follows. The accepted candidate
requirements remain:

| Replay leg | Required result |
|---|---:|
| ITS standalone and combined | 212 / `46913a67a7e2fe7462e29df0db264fa8` |
| MFT standalone and combined | 68 / `8106b08571ca593c6b76ff72b761a680` |

The pinned environment ran `ctest -L itsmft --output-on-failure -j1` with
96/96 registered tests passing and no `Not Run` entries. The fixture manifest
passed 43/43 both before and after the replays. Both standalone legs and both
combined legs produced the required counts and hashes. The persisted ITS
and MFT products matched their corresponding standalone products
field-by-field, including ROFs, cluster indices, seed patterns, labels, and
all initialized writer leaves; the combined products also matched the M7c
parent field-by-field. MFT comparison excluded only
`MFTTrack.mInvQPtSeed`, the known undefined byte artifact; the 2,992
float-projected MFT values had zero maximum absolute and relative delta.
The available environment had no `nvcc`, `hipcc`, `nvidia-smi`, or `rocminfo`,
so no GPU/device validation is claimed.

No M7d validation claim changes the known exclusion for the undefined
`MFTTrack.mInvQPtSeed` byte artifact. No ADR is added: this slice implements
the ownership decision already established by ADR 0007 and the M7a design.

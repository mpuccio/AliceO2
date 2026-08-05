# Design note 0004: M7b runtime-count authority

Status: M7b implementation and validation record
Date: 2026-08-05
Scope: `Detectors/ITSMFT/common/tracking`
Companion plan: [GenericTrackingEngineMigration.md](../GenericTrackingEngineMigration.md)
Preceding design: [M7 runtime-plan Tracker/TrackerTraits migration](0003-m7-runtime-plan-tracker-migration.md)
Architecture decision: [ADR 0007](../decisions/0007-generic-tracking-engine-boundary.md)

This note closes the count-authority slice of Gate 4 M7. It does not
de-template `Tracker` or `TrackerTraits`, or change tracking physics. The
layer-topology/ROF compatibility removal is recorded in [design note
0005](0005-m7c-runtime-topology-rof.md); Tracker/TrackerTraits de-templating
remains the separately gated M7d slice described in design note 0003.

## 1. Decision and scope

The existing runtime-plan composition is the only source of semantically
runtime-sized common-core bounds:

| Need | Authority | M7b rule |
|---|---|---|
| Traversal order and extent | `SurfacePlanBinding::getOrderedSurfaces()` | Iterate the validated ordered positions. Never reconstruct order from numeric `SurfaceId` values. |
| Active surface workspace | `SurfaceTrackingScratch::getNOwnedSurfaces()` | Size host caches, per-surface working vectors, and shared hot-loop bounds from the adopted plan. |
| Transition and cell capacity | `SurfaceTrackingScratch::getNTransitions()` / `getNCells()` plus the binding's compact slot maps | Use runtime topology counts and binding positions; preserve fixed device storage only where its ABI requires it. |
| Whole-track membership | positional `TrackSeed::SurfaceMask` and binding position maps | A mask bit is a plan position, not a global `SurfaceId`. |
| Normalized measurement order | `SurfacePlanBinding` order and source | Measurement spans are indexed by ordered position and retain their source-qualified `ClusterRef`. |
| Index-table populated prefix | runtime active count and extent span | `IndexTableUtilsCore` retains its fixed maximum device capacity, but only the runtime prefix is configured and compared. |

`DetectorLayoutSet`/`DetectorLayoutView` remains the static descriptor and
sparse-topology authority. `SurfacePlanBinding` now retains the validated
ordered surface span it was built from, alongside the global-to-compact maps.
`SurfaceTrackingScratch::adoptPlan()` remains the mutable event-workspace
boundary. No `RuntimeTrackingPlan` wrapper, compatibility alias, detector
switch, or second runtime implementation was introduced.

`TrackingParameters::NLayers` remains only as an adapter-edge consistency
check against the adopted active count. It is not used as a common hot-loop
bound or as the size of a plan-sized host temporary after this slice.

## 2. Production changes

The shared implementation was converted in place:

- `TrackerTraits` material, measurement, reachable-surface, scattering-angle,
  first-cluster, and accepted-track working state now uses runtime-sized
  vectors/spans or scratch counts. The existing `TrackerTraits<NLayers>` type
  remains the explicitly deferred M7d boundary; its common topology/ROF
  dispatch has since been removed by M7c.
- `SurfaceTrackingScratch` derives initialization, cluster preparation, and
  vertexing extents from the adopted plan count. It now receives one
  non-owning runtime ROF view group and the sparse topology directly; fixed
  table builders remain at named application edges.
- `SurfacePlanBinding` exposes its validated ordered positions. The tracking
  interface, participant publication export, and combined workflow access that
  order rather than the layout's numeric identity.
- Index-table binding and comparison receive an explicit active count. The
  fixed device-portable index-table arrays retain maximum capacity and use a
  runtime populated prefix.
- Host-only refit measurement inputs and leg buffers use runtime spans/vectors.
  `LayerMeasurementSpans<NLayers>` was removed; `CellSeed` and `TrackSeed`
  remain their existing fixed-capacity device value types.
- The native `Propagator` calls, covariance handling, material policy,
  tracklet/cell/road algorithms, holes, ordering, output adapters, and
  workflow ownership are unchanged.

The only application-facing order change is representational: combined
workflow accessors now read the participant binding order. The values and
source-qualified exports are unchanged.

## 3. Remaining `NLayers` boundaries

The source guard in
`test/testM7bRuntimeCountAuthorityGuard.cxx` scans every common tracking
production header/source file containing `NLayers`. Each code-level use must
fall into exactly one of these categories:

| Guard category | Current boundary | Why it remains | Exit slice |
|---|---|---|---|
| Fixed device ABI/capacity | `Cell.h`'s fixed cell/metadata representation | A CA cell has three measurement slots and the value type must remain GPU-portable. | Retain; do not turn `CellSeed` into a dynamic object. |
| Temporary private operation implementation | `IndexTableConfiguration`, `TransitionPolicyOperations`, `NativeRefitDriver`, `NativeCylinderCylinderRefitDriver`, and `RefitLegAssembly` typed operation signatures | These are private typed operation/template boundaries. Their host traversal extents now come from runtime spans; the tag and operation templates are not public scheduling policy. | Reduce after M7d/M7e when one runtime core and adapter edge exist. |
| Adapter edge | `TrackingParameters::NLayers` validation, `DetectorTraits`, configuration/load policy, detector surface specs/catalogs, publication/refit output hooks, and the native-cylinder typed output loop | These objects still construct detector-specific typed output, ROF/configuration compatibility, or frozen adapter data. They are outside generic orchestration. | M7e and later adapter cleanup; no core routing is allowed here. |
| Explicitly deferred M7d compatibility boundary | `Tracker`, `TrackerTraits`, `CATracker`, `TransitionPolicyState`, and the remaining typed core bridge | The runtime-plan migration has removed common layer topology, ROF-table scratch ownership, and the `IndexTableUtils<N>` alias. The remaining boundary is the Tracker/TrackerTraits template/state-family composition. | M7d removes the common Tracker/TrackerTraits templates; adapter-only templates remain separately owned. |

The guard fails if a new production file or an unlisted line introduces an
unclassified `NLayers` use. It also rejects the plan-sized host arrays and
`.NLayers` count reads in `SurfaceTrackingScratch` and shared
`TrackerTraits` initialization. The guard does not scan frozen legacy ITS
code outside `common/tracking/{include,src}`; that exclusion is path-specific
and does not exempt any common-CA production file.

Fixed capacities are deliberately not conflated with active counts:

- `TrackSeed` keeps `MaxLayoutSurfaces` cluster slots and its positional
  `SurfaceMask` because it is a trivially-copyable GPU-portable whole-track
  value.
- `CellSeed` keeps its three-measurement CA-cell representation and the live
  `SeedMetadataBase<N>` base.
- `IndexTableUtilsCore` keeps its maximum device capacity and receives an
  active prefix/count. No dynamic device allocation is introduced.
- Host-only material, measurement, refit-slot, and scheduling temporaries are
  sized from the runtime plan instead of either detector's 7/10 width.

## 4. Sparse-plan proof

`testSurfacePlanBinding` builds a dense catalog with a sparse, non-identity
application order `{5, 2, 7, 1}` and a source-qualified binding. It verifies:

1. the binding preserves the exact ordered positions and maps each global
   surface to its compact position;
2. scratch surface, transition, and cell extents follow the binding/runtime
   topology counts;
3. normalized measurement views retain both position order and source; and
4. `TrackSeed::SurfaceMask` bits follow positions, so global surface 5 does
   not accidentally become mask bit 5.

The migrated index-table/refit fixtures continue to cover fixed-capacity
prefix behavior, forward/reverse legs, holes, and staging-before-commit
failure behavior.

## 5. Ranked deletion and simplification inventory

### Safe in the runtime-count slice

- Host `std::array<T, NLayers>` temporaries replaced by runtime vectors/spans.
- `LayerMeasurementSpans<NLayers>` removed after all common consumers moved to
  runtime span-of-spans inputs.
- Numeric-ID order reconstruction removed from the binding/participant load
  and publication path.
- Index-table configuration now accepts explicit active counts and validates
  the populated prefix.

### Blocked until M7d/M7e

- Dual `ROF*Table<N>` storage and templated scratch accessors: completed in
  M7c. Fixed-capacity table builders remain only at named adapter edges until
  their application template seam is reduced.
- `TrackingTopology<N>` layer graph and `IndexTableUtils<N>` compatibility
  alias: deleted from the common path in M7c; frozen legacy ITS compatibility
  remains outside the common-CA guard.
- `Tracker<NLayers>`, `TrackerTraits<NLayers>`, `CATracker` aliases, and
  their explicit core instantiations: blocked until M7b's count proof is
  followed by the single-body de-templating in M7d.
- `DetectorTraits<N>`, typed accepted-track staging, and remaining adapter
  templates: blocked on application-owned refit/output hooks and writer parity
  in M7e.

### Deferred physics or algorithm decisions

The following remain explicitly out of scope: changing Propagator equations or
covariance sanitization; changing material, hole, candidate-order, MinPt, or
MFT pre-cut behavior; changing transition-policy semantics; changing fixed
device widths/capacities; making currently rejected mixed cylinder/disk plans
legal; or deleting any legacy numerical path based only on similarity. Each
would require its own physics/algorithm decision and replay evidence.

## 6. Validation record

The durable build at
`/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement`
was reused with package `daily-20260717-0700-local1`. All 95 registered
ITS/MFT tests executed serially and passed; there were no `Not Run` tests.
The fixture manifest passed 43/43 both before and after replay.

The canonical standalone and combined replays produced the accepted anchors:

| Replay leg | Tracks | Content hash |
|---|---:|---|
| ITS standalone | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT standalone | 68 | `8106b08571ca593c6b76ff72b761a680` |
| ITS combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT combined | 68 | `8106b08571ca593c6b76ff72b761a680` |

The repository persisted-output comparator found the combined ITS leg
field-identical to the standalone leg. The MFT legs matched all initialized
products, including all 2,992 float-projected values with zero maximum
absolute/relative delta; the only excluded value was the known undefined
`MFTTrack.mInvQPtSeed` leaf. The same initialized-field comparison passed for
the retained accepted parent products used by the M6e3/M6g validation chain.
No other writer, sidecar, ROF, label, cluster-reference, or `CommonTrack`
content differed.

The source guard and `git diff --check` passed. The real device toolchain was
not available on this Darwin arm64 host (`nvcc`, `hipcc`, `nvidia-smi`, and
`rocminfo` were absent), so no GPU/device build or replay is claimed. The
acceptance anchors remain:

- ITS: 212 tracks, `46913a67a7e2fe7462e29df0db264fa8`;
- MFT: 68 tracks, `8106b08571ca593c6b76ff72b761a680`.

No GPU/device validation is claimed unless a real pinned `nvcc`/`hipcc`
toolchain and device are available. M7b does not alter frozen ITS workflows.

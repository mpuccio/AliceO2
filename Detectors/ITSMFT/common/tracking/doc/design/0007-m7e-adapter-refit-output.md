# M7e: typed refit and output at the adapter boundary

Status: complete, 2026-08-05
Base: M7d `87574740dd`
Implementation commits: `b54b264c3b`, `cc2a40b71d`

## Scope and architectural result

M7e removes detector-typed refit and publication ownership from the common
tracking core. `Tracker`, `TrackerTraits`, and `TrackingCandidate` now operate
on the runtime-plan data model: `TrackSeed`, `SurfaceKinematicState`,
`SurfaceMeasurement` spans, `CommonTrack`, source-qualified references, and
the existing plan/workspace state. The core publishes each accepted
`CommonTrack` into `TimeFrame` in deterministic acceptance order. It does not
construct or retain a typed accepted-track vector.

The single remaining operation seam is `TrackingOperationAdapter`. It is a
call-scoped adapter operation, not a coordinator or a second tracker:

| Operation | Common-core input/result | Adapter responsibility |
|---|---|---|
| Seed refit | `TrackSeed`, runtime measurement spans, `SurfaceKinematicState` result | Select the detector leaf operation and return a generic `TrackingCandidate` |
| Accepted completion | ordered generic candidates with their `CommonTrack` indices | Materialize only detector compatibility sidecars |
| Reset | no event data | Clear only adapter compatibility state; `TimeFrame`/scratch reset stays in the core/workflow contract |

The adapter boundary is intentionally operation-oriented. It carries no
`TrackITSExt`, MFT legacy track, detector pattern, writer, DPL, workflow, or
detector-ID type into `Tracker` or `TrackerTraits`. The MFT typed normalized
refit overload is isolated in `MFTAdapterRefit`; the generic forward leaf
operations remain in `MFTFwdTrackHelpers`. ITS native refit remains the shared
`Propagator` path, selected by the ITS adapter.

## Deleted compatibility ownership

The following completed bridges are deleted rather than renamed:

| Deleted item | Replacement/owner |
|---|---|
| `DetectorTraits<NLayers>` and `DetectorTraits.cxx` | `DetectorTrackingOperationAdapterSupport` and the explicit ITS/MFT adapter edge |
| `CATrackType<NLayers>` and typed accepted aliases | `TrackingCandidate` plus `TimeFrame::getCommonTracks()` |
| `LayerMeasurementSpans<NLayers>` | non-owning runtime measurement spans resolved by `TrackerTraits` |
| `AcceptedTrackShadowPublisher` | owner-thread generic `CommonTrack` publication followed by adapter sidecar completion |
| typed MFT refit in the generic forward helper | `MFTAdapterRefit.{h,cxx}` |

`MFTCATrack` remains only where the MFT adapter must preserve its typed
normalized-refit compatibility contract. It is not included by
`TrackerTraits`, `Tracker`, or the generic forward helper. `SeedMetadataBase<N>`
also remains because live `CellSeed` still uses it; it is not whole-track
compatibility state.

The ITS adapter now stages the existing shared-cluster selection from generic
candidate kinematics, the operation-local tracking parameters, and the
read-only scratch ROF view. The core retains the same acceptance/order
semantics but carries no shared-cluster flag in `TrackingCandidate` and sends
no typed or compatibility status through the operation callback. The adapter
commits the resulting sidecar only on the final accepted-result boundary.

## Adapter contracts

The ITS adapter converts the generic accepted sequence into
`ITSSharedClusterCompatibility` entries keyed by the committed CommonTrack
index. It preserves accepted order and shared-cluster status, and the existing
CommonTrack output adapter later reconstructs `TrackITS`/`TrackITSExt` fields
from CommonTrack, source-qualified measurements, and that sidecar.

The MFT adapter converts the same generic accepted sequence into
`MFTPublicationCompatibility` entries. The retained compatibility values are:

- `seedPattern`: the fixed-capacity `TrackSeed` hit mask;
- `invQPtSeed`: the generic seed q/pt value;
- legacy seed/output chi2 fields: their existing zero/default contract.

The old persisted `MFTTrack.mInvQPtSeed` value was an uninitialized legacy
writer byte. M7e does not make that undefined field a new physics or CommonTrack
contract. It remains the sole explicitly excluded byte artifact in parent
writer comparisons; all initialized MFT fields, references, ROFs, labels,
sidecars, and the approved float projection remain checked.

Both adapters use the same generic result completion call. There is no
detector switch, type-erased track, compatibility alias, or parallel tracking
implementation.

## Failure and lifecycle ownership

The core's existing all-or-nothing contracts are unchanged:

1. A failed normalized load leaves the prior event untouched.
2. A dropped timeframe or structural tracking failure resets the generic
   `TimeFrame`, scratch, and adapter compatibility state together.
3. Successful replacement first commits the new normalized event, then
   publishes new CommonTracks and sidecar entries.
4. Adapter completion is atomic at its sidecar boundary; malformed or
   non-monotonic CommonTrack indices fail without exposing a partial sidecar.

Raw ROFs, publication clocks, writer sidecars, and DPL event sequencing remain
at their existing workflow/application owners. M7e does not move those into
the generic engine.

## Guards and focused coverage

`testM7eAdapterBoundaryGuard` scans the generic Tracker/TrackerTraits,
CATracker, and operation-adapter declarations/definitions for typed track,
detector identity, DPL/workflow/writer, `DetectorTraits`, `CATrackType`, and
`LayerMeasurementSpans` dependencies. It also rejects production declarations
of the deleted typed forms outside the narrowly named adapter/frozen-legacy
edges. The generic MFT helper and the MFT adapter split are checked
independently.

`testDetectorAdapterCompatibility` covers both ITS and MFT generic accepted
result conversion, exact CommonTrack indexing, shared status, seed mask,
seed q/pt, and invalid native-refit rejection. Existing failure-contract,
normalized-refit, beam-position, M7b, M7c, and M7d tests remain registered and
were rebuilt against the new seam.

## Validation record

The durable build reused package `daily-20260717-0700-local1` and the existing
M7d build directory. The affected common library, ITS/MFT/combined workflow
libraries and executables, and tests built successfully. The complete serial
selector executed all 97 registered ITS/MFT tests: **97/97 passed**, with no
`Not Run` result.

The fixed fixture manifest passed **43/43 before and after** replay. Fresh
standalone and combined replays produced the required candidates:

| Replay leg | Tracks | Content hash |
|---|---:|---|
| ITS standalone | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS combined leg | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT standalone | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT combined leg | 68 | `8106b08571ca593c6b76ff72b761a680` |

The combined ITS and MFT products matched their standalone products
field-by-field, including ordered cluster references, ROFs, seed patterns,
labels, sidecars, and the MFT float-projected fields (2,992 values, zero
maximum absolute and relative delta). Comparison with the accepted M7d parent
matched every initialized ITS and MFT writer field/product. The only excluded
parent mismatch was `MFTTrack.mInvQPtSeed`, the known undefined byte artifact.

The pinned environment contained no `nvcc`, `hipcc`, `nvidia-smi`, or
`rocminfo`; no GPU/device validation is claimed.

## M7f follow-up and final cleanup

M7f completed the safe compile-time bridge deletions and guard closure. The
exact residual exception table, validation record, and remaining post-M7
decisions are in [design note 0008](0008-m7f-final-runtime-core-cleanup.md).
The following items remain intentionally outside that behavior-preserving
cleanup:

| Rank | Item | Gate |
|---:|---|---|
| 1 | Remove or reduce `TrackingOperationAdapter` once the generic accepted-result collection has a stable completion owner; delete the call-scoped seam only after the same sidecar/failure contract is expressed without a callback. | Deferred post-M7 ownership decision |
| 2 | Revisit `TrackingCandidate`'s adapter-only kinematics (`phi`, `eta`, `charge`) and any staging that remains after a stable generic accepted-result owner exists. | Safe structurally, but requires a separately reviewed ownership change |
| 3 | Delete `MFTAdapterRefit`'s typed `MFTCATrack` overload and remaining legacy state-export adapters once MFT typed refit/output is fully represented by the generic result path. | Blocked until a separately approved typed-output decision |
| 4 | Reduce `DetectorTrackingOperationAdapterSupport.h`, `MFTAdapterRefit.cxx`, and the adapter-side `NLayers` mixins after all typed leaf calls are removed. | Blocked until item 3 |
| 5 | Reduce `SurfacePlanTrackingParticipant<NLayers>` and `ITSMFTTrackingInterface<NLayers>` only after application-edge configuration/ROF parity is separately proven. | Deferred post-M7 adapter decision; M7f records the exact owner |
| 6 | Revisit remaining adapter-edge `NLayers` builders and fixed-capacity ROF compatibility objects. | Deferred compatibility decision |
| 7 | Revisit MFT leaf-helper and native cylinder test-oracle simplifications only if they change operation selection or numerical behavior. | Deferred physics/algorithm sign-off |

M7f did not remove the operation seam or de-template the remaining
workflow-facing adapter participants. It did not start a new Tracker/
TrackerTraits milestone or relocate workflow ownership.

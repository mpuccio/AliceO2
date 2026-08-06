# Design note 0002: M6 generic workspace — audit and deletion-oriented migration design

Status: M6g complete; design/audit and implementation record
Date: 2026-08-04
Scope: `Detectors/ITSMFT/common/tracking` (`o2::itsmft::tracking`)
Companion: [ADR 0007](../decisions/0007-generic-tracking-engine-boundary.md) decisions 2, 9, 10;
milestone M6 of [GenericTrackingEngineMigration.md](../GenericTrackingEngineMigration.md)
(this note supersedes that document's single-paragraph M6 sketch with the staged M6a–M6g
sequence below; see the edits to that file made alongside this one).
Companion design note: [design note 0001](0001-descriptor-driven-operation-boundary.md)
(the tracklet/cell/refit *algorithm* boundary; this note is the *container/ownership*
boundary one layer below it and does not revisit any of its conclusions).

The source inventory in §1 is the M6a audit snapshot. The implementation records
below document the subsequent behavior-preserving M6b–M6g cleanup. **No Propagator,
candidate physics, or fixture changes are authorized or made by this note.**
M5d's current candidate results (ITS 212, MFT 68 — [ADR 0008](../decisions/0008-native-refit-activation.md))
are read here only as context; they are not touched, re-derived, or re-approved by this
work, and remain not-yet-physics-approved.

## 1. Method

Every claim below is grounded in the current source, read in full, not asserted from
architecture docs alone:

- `include/ITSMFTTracking/LegacyTrackerScratch.h` / `src/LegacyTrackerScratch.cxx`
  (every public/protected data member and every method body, including
  `loadNormalizedSource()`/`resetScratch()`/`initialise()`/`setMemoryPool()`);
- `include/ITSMFTTracking/detail/DetectorTraversalBinding.h` (the full `build()` body and
  every accessor);
- `include/ITSMFTTracking/LegacyCATrackingParticipant.h` / `src/LegacyCATrackingParticipant.cxx`;
- `include/ITSMFTTracking/ITSMFTLegacyParticipantSet.h` / `src/ITSMFTLegacyParticipantSet.cxx`;
- `include/ITSMFTTracking/{TrackingParticipant,TrackingEngine,TimeFrame}.h` (the permanent
  generic boundary these four migrate toward);
- `include/ITSMFTTracking/MultiSourceTimeFrameLoader.h` / `.cxx` (already
  participant-count-generic atomic loading, M2b — unaffected by this note);
- `include/ITSMFTTracking/{ITSSharedClusterCompatibility,MFTPublicationCompatibility}.h`
  (the two detector-output sidecars);
- `include/ITSMFTTracking/Cell.h` (`SeedMetadataBase<N>`, `CellSeed`, `TrackSeedTpl<NLayers>`,
  `CATrackTypeHelper<NLayers>`) and `include/ITSMFTTracking/DetectorLayoutSet.h` (the
  already-dynamic, already-sparse plan type this note's workspace sizes itself from);
  `src/TrackerTraits.cxx` (every `mBinding->...`/`mScratch->...` call site, grepped and
  read in context, not assumed);
- `include/ITSMFTTracking/CATracker.h` (`Tracker<NLayers>::clustersToTracks()`'s own
  recoverable-failure reset path);
- [ADR 0007](../decisions/0007-generic-tracking-engine-boundary.md),
  [ADR 0008](../decisions/0008-native-refit-activation.md),
  [design note 0001](0001-descriptor-driven-operation-boundary.md), and
  [Architecture.md](../Architecture.md) §8–10 (sparse topology, indexing, policy boundary).

## 2. Verdict (summary; see §9 for the full reasoning)

**GO**, staged as six bounded implementation slices (M6b–M6g) behind this audit (M6a).
The generic workspace this note proposes is smaller than `LegacyTrackerScratch<NLayers>`,
not a rename of it: most of the class's content is *already* plan-sized at runtime (the
`DetectorLayoutSet`/`SparseTrackingTopology` it is built from is already dynamic — see
§4.1); the only real `NLayers`-as-compile-time-bound leaks are (a) the per-owned-surface
array width used throughout the legacy cluster/index-table cache and (b) one type,
`TrackSeedTpl<NLayers>`/`SeedMetadataBase<NLayers>` — and the codebase already documents
the second one as "a temporary legacy boundary, not a detector-selected state family"
(`Cell.h`, `TrackSeedTpl`'s own doc comment). `DetectorTraversalBinding`'s only two
detector-switched lines (an `ITS`/`MFT` allow-list and an
`ITS→CylinderCylinder`/`MFT→DiskDisk` tag lookup) are the only genuine ITS/MFT leak in
that class; everything else in it is already detector-neutral.

**Correction (this revision) — permanence of the ITS/MFT adapter *role* does not make
`LegacyCATrackingParticipant`/`ITSMFTLegacyParticipantSet`, as currently named and shaped,
permanent *classes*.** An earlier draft of this note conflated the two and declared both
classes "not deletion targets"; that was wrong, and is corrected throughout §3.3, §3.4,
§8–§10 below. The rule this revision applies instead: **the ITS/MFT application-adapter
responsibility survives M6 permanently (ADR 0007 decision 2).** M6f removed the old
workspace/binding participant name. M6g removed the remaining participant-set exception;
no `Legacy…`-named common-CA production class may survive. `LegacyCATrackingParticipant`
is deleted or renamed to a narrowly-scoped,
non-`Legacy` ITS/MFT participant type once its `LegacyTrackerScratch`/`DetectorTraversalBinding`
internals are gone (§3.3, §9 M6f). `ITSMFTLegacyParticipantSet` survives M6 **only if**,
at M6 completion, it is demonstrably an immutable application-plan/configuration builder —
no `loadEvent`/`process`/`execute`/`reset` behavior, no event-owned mutable state, no
coordinator role, no detector-tracking-algorithm ownership beyond constructing participant
adapters (§3.4). Audited against those four conditions as currently shaped, it **fails**
today (§3.4); it is retained under a new, non-`Legacy` name only if a further slice (M6g)
relocates its disqualifying state to the workflow adapter, and deleted/inlined into
`CombinedCATrackerSpec.cxx`'s composition otherwise. M6g took the latter path: the
publication/timing and per-event reset state moved to the combined DPL task, and no
coherent immutable holder remained after that relocation.

## 3. Ownership/classification table

Categories, as specified:
**(1)** generic per-participant CA working state; **(2)** shared event state already owned
by `TimeFrame`; **(3)** detector/output compatibility state, adapter-private; **(4)**
temporary migration bridge; **(5)** dead/redundant, delete rather than migrate.

Fields are grouped where an entire group shares owner/lifetime/reset/consumers —
individually listing all ~50 `LegacyTrackerScratch` data members would obscure, not
clarify, the actual boundary; per-field detail is in §4 and §7 where the classification
of an individual field diverges from its group.

### 3.1 `LegacyTrackerScratch<NLayers>`

| Group | Members | Owner / lifetime | Reset path | Load/failure path | Publication dependency | Consumers | Classification |
|---|---|---|---|---|---|---|---|
| **A. Legacy per-(owned-surface) cluster/index-table cache** | `mClusters`, `mUnsortedClusters`, `mTrackingFrameInfo`, `mClusterExternalIndices`, `mClusterSize`, `mROFramesClusters`, `mClusterLabels`, `mIndexTables`, `mIndexTableUtils`, `mUsedClusters`, `mNClustersPerROF`, `mMinR`/`mMaxR`, `mBogusClusters`, `mPositionResolution`, `mNTrackletsPerCluster`/`mNTrackletsPerClusterSum` | `LegacyTrackerScratch<NLayers>`, per-participant, event-scoped | `resetScratch()`/`resetROFrameData()` (`deepVectorClear`, `mClusterLabels` pointer-nulled not deep-cleared) | Committed only by `loadNormalizedSource()`'s legacy-backfill half of `MultiSourceTimeFrameLoader::loadEvent()`'s atomic transaction (§6); strong exception safety, allocator-identity gate before any commit | None directly | `TrackerTraits.cxx` hot loop: index-table binning (`prepareClusters()`), a `NormalizedMeasurementMismatch` cross-check against `TimeFrame::getNormalizedFrame()`, and `getClusterExternalIndex()`/`getClusterLabels()` MC-label indirection | **(1)** generic per-participant CA working state, but array-bound by `NLayers` where it should be bound by `ownedSurfaces().size()` (a runtime, plan-derived count — see §4.1). Not a `TimeFrame`-duplicate by policy (§3.1 note below) even though its *positional/covariance content* overlaps `TimeFrame`'s normalized frame. |
| **A′. Dead sub-path** | `loadROFrameData()`, `resetROFrameData()`, `prepareROFrameData()` | same | — | — | — | **Zero production call sites** in the new common tracker (grep-verified: no call in `ITSMFTLegacyParticipantSet.cxx`, `LegacyCATrackingParticipant.cxx`, `CATracker.cxx`, `TrackerTraits.cxx`, `TrackingInterface.cxx`). Referenced only by their own definitions and two tests (`testTimeFrameCovarianceLifecycle.cxx`, `testTimeFrameNormalizedSource.cxx`). The only committed loading path in production is `loadNormalizedSource()`. | **(5)** dead in the new common tracker — flagged as a deletion opportunity, **not** folded into the mandatory M6 migration (§12, Rank 1). |
| **B. Plan-sized CA construction/result transients** | `mTracklets`, `mTrackletsLookupTable`, `mTrackletLabels`, `mCells`, `mCellsLookupTable`, `mCellsNeighbours`, `mCellsNeighboursTopology`, `mCellsNeighboursLUT`, `mCellLabels`, `mTransitionPhiCuts`, `mTransitionMSAngles`, `mTrackerTopologies`/`mDefaultTrackingTopology`/`mVertexingTopology`/`mTrackingTopologyView` | `LegacyTrackerScratch<NLayers>`, resized every `initialise()` call (per iteration) from `mTrackingTopologyView.nCells`/`nTransitions` — **already a runtime count derived from the sparse `DetectorLayoutSet`/topology, not from `NLayers`** (§4.1) | `resetScratch()` (`deepVectorClear`) | Populated entirely inside the hot loop (`TrackerTraits::computeLayerTracklets`/`computeLayerCells`/etc.), not by loading | None directly | `TrackerTraits.cxx` end to end (tracklet/cell/road/refit orchestration, unchanged by this note — see design note 0001) | **(1)** generic per-participant CA working state. **This is the exact scope of the proposed `SurfaceTrackingScratch` (§4)** — already plan-indexed in every dimension except the surrounding class template parameter. |
| **C. Legacy per-detector result staging** | `mTracks` (`bounded_vector<CATrackType<NLayers>>`), `mTracksLabel` | same | `resetScratch()` | Populated by `TrackerTraits::acceptTracks()` | Read by the owning `LegacyCATrackingParticipant<NLayers>`'s adapter caller (`ITSMFTLegacyParticipantSet::getITSScratch()`/`getMFTScratch()`) for detector-typed output construction | The same `acceptTracks()` call site **already also** publishes a `CommonTrackShadowRecord` into `TimeFrame::getCommonTracks()`/`getTrackClusterIndices()` in the same transaction (`TrackerTraits.cxx:2028`, `mAcceptedTrackShadowPublisher.publish()`) — `CommonTrack` is **not** unpopulated in production today, contrary to `TimeFrame.h`'s own now-stale doc comment (§7.3) | **(4)** temporary migration bridge. `CATrackType<NLayers>` (via `CATrackTypeHelper<NLayers>`) is genuinely detector-output-typed (`TrackITSExt`-shaped for ITS, MFT-specific for MFT) — that destination stays **(3)** adapter-private — but the intermediate legacy staging container duplicates data `CommonTrack` already carries generically. M6e (§9) retires it once adapters build detector output from `CommonTrack` + `ParticipantPublicationExport` alone. |
| **D. Vertexer working scratch** | `mLines`, `mTrackletClusters`, `mNTrackletsPerROF`, `mTrackletsIndexROF`, `mLinesLabels`, `mTotalTracklets`/`mTotalLines`, `mNTotalLowPtVertices` | same | `resetScratch()` (all but `mNTrackletsPerROF`/`mTrackletsIndexROF`, which are only resized in `initialise()`) | Populated by the vertexer stage, upstream of tracklet/cell/road | None directly (result is `TimeFrame::getPrimaryVertices()`, already **(2)**) | Vertexer only | **(1)** generic per-participant CA working state. Never actually `NLayers`-sized (bound by ROF count and a fixed layer-pair index of 2, not `NLayers`) — it inherited the class template only by being declared inside the same struct. Moves into the same generic workspace automatically; **no separate design decision needed**, and this note does **not** propose touching vertexing algorithms (out of the tracklet/cell/road/refit scope this ADR/note governs). |
| **E. Memory/allocator/device plumbing** | `mExtMemoryPool`, `mMemoryPool`, `mExternalAllocator`, `mIsStaggered`, `mPValphaX` (declared, `initVector`-cleared in `setMemoryPool()`, no production write/read site found — flag for verification, not asserted dead) | same | `setMemoryPool()`/`setFrameworkAllocator()` | — | — | Every method | **(1)** generic, already detector-neutral in content and type. Unaffected by this migration beyond moving to the new class. |

**Note on Group A vs. `TimeFrame`'s normalized frame**: Group A's positional/covariance
content is largely redundant with `TimeFrame::getNormalizedFrame()`'s already-owned
per-`SurfaceId` `SurfaceMeasurement`s (confirmed by the `NormalizedMeasurementMismatch`
cross-check in `TrackerTraits.cxx` §2.6, which exists specifically because the two
representations must be kept in agreement). Eliminating that redundancy would mean
rewriting the tracklet/cell hot loop to index directly by `SurfaceId`/`SurfaceMeasurement`
instead of by legacy-layer-indexed `Cluster`/`TrackingFrameInfo` — an **algorithm input**
change, not a container-ownership change, and out of scope both for M6 (per the milestone's
own text: "delete `LegacyTrackerScratch<NLayers>` and the remaining `NLayers`-templated
legacy CA containers") and for this audit's explicit "do not redesign tracklet/cell/road
algorithms" constraint. Flagged here as a **future, separately-decided** opportunity
(tentatively "M7"), not part of M6's mandatory scope and not one of §12's ranked
candidates either (it needs an algorithm-input redesign §12's three items do not).

### 3.2 `DetectorTraversalBinding` responsibilities

| Responsibility | Detector-specific? | Classification | Disposition |
|---|---|---|---|
| `source.isValid()` check | No | **(1)** generic | Keep unchanged in the successor. |
| Ownership/surface-mask validation against the global layout (`ownedSurfaces.isSubsetOf(...)`, count check) | No | **(1)** generic | Keep unchanged. |
| `legacySurfaceOrder` → `mLegacyLayerBySurface` (global `SurfaceId` → compact position) mapping | No — this is exactly "position in `ownedSurfaces()`", a plan-derived index, not an ITS/MFT concept | **(1)** generic (named "legacy" only because its sole consumer today is Group A above) | Keep the mechanism; the successor type (§7) drops "legacy" from its own naming since it is the authoritative source of Group A's new plan-derived array bound too. |
| Global transition/cell topology validation (cross-boundary edge/skip/hit-surface subset checks) | No | **(1)** generic | Keep unchanged — already exactly the boundary a future ALICE 3 leg would need. |
| Compact scratch-slot assignment (`mScratchTransitionSlot`/`mScratchCellSlot`, `mGlobalTransitions`/`mGlobalCells`/`mGlobalRoadStartCells`/`mGlobalScheduledCells`) | No | **(1)** generic | Keep unchanged — this **is** the index `SurfaceTrackingScratch`'s Group-B containers are keyed by. |
| `detector != ITS && detector != MFT` → `UnsupportedDetector` gate | **Yes** | **(4)** temporary bridge | Delete. Ownership/topology validation above already fully determines correctness; this gate adds nothing an ALICE 3 leg would need and everything it would block. |
| `expectedKind`/`expectedPolicy` derived by switching on `detector` (`ITS→Cylinder/CylinderCylinder`, `MFT→Disk/DiskDisk`) | **Yes** | **(4)** temporary bridge (one more call site of the already-ADR-0007-classified-temporary `TransitionPolicyTag`) | Take `expectedKind`/`expectedPolicy` as a caller-supplied parameter, derived by the adapter from its own plan's surface kinds (mirrors M4b's existing `stateFamilyOf(SurfaceKind)`/decision 8's `transitionPolicyTagForSurfaceKind()`) — not a `build()`-internal detector switch. |

Net: **~90% of `DetectorTraversalBinding` is already detector-neutral.** The only leak is
two call sites inside `build()`. This is what makes M6b (§9) a small, bounded first slice.

### 3.3 `LegacyCATrackingParticipant<NLayers>` responsibilities

| Responsibility | Classification | Disposition |
|---|---|---|
| Owns `ScratchN mScratch` (`LegacyTrackerScratch<NLayers>`) | **(1)**/**(4)** (per §3.1) | Member *type* becomes `SurfaceTrackingScratch` (M6d/M6e); the participant class itself is unaffected. |
| Owns `std::unique_ptr<DetectorTraversalBinding> mBinding` | **(4)** (per §3.2) | Member *type* becomes the successor binding type (M6b/M6d). |
| Owns `TrackerTraits<NLayers> mTraits`, `Tracker<NLayers> mTracker` | **out of this audit's scope** | These remain `NLayers`-templated: they are the tracklet/cell/road/refit *algorithm orchestration* design note 0001 already classified and this note's own constraints forbid touching ("do not redesign tracklet/cell/road algorithms"). `LegacyCATrackingParticipant<NLayers>` therefore also **stays a class template** — its `NLayers` parameter selects which `Tracker`/`TrackerTraits` specialization it drives, independent of whether its scratch/binding members are generic. |
| Owns `ITSSharedClusterCompatibilityOwner<NLayers>`/`MFTPublicationCompatibilityOwner<NLayers>` sidecar mixins | **(3)** detector/output compatibility, adapter-private | Unaffected by M6 (§8) — these were never part of `LegacyTrackerScratch`/`DetectorTraversalBinding` and stay exactly where they are. |
| `loadTarget()` (`MultiSourceTimeFrameLoader::LoadTargetImpl<NLayers>`, bound to `mScratch`) | **(1)**/**(4)** | `LoadTargetImpl<NLayers>` moves from staging a `LegacyTrackerScratch<NLayers>` copy to staging the new scratch type's Group-A-equivalent containers; the atomic-loading *contract* (§6) is unaffected. |
| `track()`/`eventReset()`/`publicationExport()` (`TrackingParticipant` contract) | **(2)**-adjacent (touches `TimeFrame` only via the documented boundary) | Unaffected — already the permanent generic contract (ADR 0007 decision 5). |
| `getDropTFUponFailure()` | **(3)** adapter-private policy | Unaffected. |

**Corrected disposition**: the *concept* this class implements — a concrete, per-detector
`TrackingParticipant` driving `Tracker<NLayers>`/`TrackerTraits<NLayers>` — is permanent;
ITS and MFT (and any future ALICE 3 leg) always need exactly this kind of adapter. The
*class as currently named*, `LegacyCATrackingParticipant<NLayers>`, is not: its own name
records that it wraps `LegacyTrackerScratch<NLayers>` and `DetectorTraversalBinding`, both
of which M6d/M6e/M6f delete. Per the stated rule that no `Legacy…`-named class may survive
the completed M6 migration, once those two members are the generic `SurfaceTrackingScratch`/
`SurfacePlanBinding` types (end of M6e), the class itself is renamed — at M6f — to a
narrowly-scoped, non-`Legacy` ITS/MFT participant name (e.g. a shared
`CATrackingParticipant<NLayers>` template with `ITSCATrackingParticipant`/
`MFTCATrackingParticipant` aliases, mirroring today's `LegacyCATrackingParticipantITS`/`MFT`
aliases minus the "Legacy" — the exact name is an M6f implementation-slice decision, not
fixed by this design note). This is a rename, not a deletion of the adapter role — but the
*name* `LegacyCATrackingParticipant` must not exist in production code once M6 completes.

### 3.4 `ITSMFTLegacyParticipantSet` responsibilities

Evaluated against the four conditions this revision requires for a participant-set-shaped
class to survive M6 (no `loadEvent`/`process`/`execute`/`reset` behavior; no event-owned
mutable state; no coordinator role; no detector-tracking-algorithm ownership beyond
constructing participant adapters) — not asserted permanent by default, per the corrected
rule in §2.

| Responsibility | Classification | Passes the 4 conditions? | Disposition |
|---|---|---|---|
| Builds the one shared ITS+MFT `DetectorLayout` (`buildCombinedLayout()`), constructor-only | **(3)** adapter-private plan/config construction | **Yes** — pure construction-time factory work, nothing per-event | Consistent with "constructing participant adapters." |
| Fixed ITS=`ClusterSourceId{0}`/MFT=`ClusterSourceId{1}` contract (`validateSources()`) | **(3)** adapter-private, stateless | **Yes** — pure function of its arguments, no stored mutable state | Unaffected in kind. |
| `loadBindings()`: builds the per-event `AtomicLoadBinding` pair passed to `MultiSourceTimeFrameLoader::loadEvent()` | **(3)**/**(4)** borderline | **Marginal** — does not itself perform the load, but it is called once per event specifically to marshal that event's load transaction, not at construction | Defensible as "config/binding construction" *if* it stays a stateless, per-call marshaling function (as it is today — it returns a fresh array, stores nothing); flagged for scrutiny at M6g, not disqualifying by itself. |
| `adoptFrame()`/`setMemoryPool()`/`setBz()`/`setNThreads()`: bind-once setup forwarded to both participants | **(3)** adapter-private setup | **Yes**, if genuinely called once per run/instance (matching `LegacyCATrackingParticipant`'s own "bind-once" contract) rather than per event | To be confirmed by the M6g caller-site audit (§9) — not contradicted by anything read for this note, but not independently re-verified either. |
| `configureRofTables()`: reconfigures both participants' ROF-timing scratch, called once per TimeFrame | **(4)** temporary-bridge-adjacent | **No** — this is a per-TF (i.e. per-event) reconfiguration call, not one-time construction, even though the mutated state itself lives inside the participants' own scratch rather than in this class | Borderline "process" behavior under condition 1; relocate to the workflow adapter's own per-TF setup step at M6g, or keep only if redefined as pure pass-through with zero judgment/state here. |
| `clearPublicationSidecars()`: forwards a per-event reset call to both participants | **(4)** temporary bridge | **No** — explicit `reset`-shaped behavior, called every event | Disqualifying under condition 1 as currently shaped. Move the call itself into the workflow adapter (`CombinedCATrackerSpec.cxx`), which already owns the per-event loop and already calls the analogous `TrackingEngine::resetEvent()`. |
| `invalidatePublication()`/`markPublicationValid()`: mutate this class's own `mITSClock`/`mMFTClock`/`mPublicationValid` every event, sequencing publication validity across both participants as a function of the whole event's success/failure | **(4)** temporary bridge, genuinely coordinator-shaped | **No** — this is event-owned mutable state stored directly in `ITSMFTLegacyParticipantSet` itself (not delegated to a participant), and the very definition of "coordinator role": sequencing cross-participant state as a function of a combined event outcome. This is the clearest disqualifying finding in this table. | Relocate the publication/timing bridge (`ClockTimingPublicationView`, `mITSClock`/`mMFTClock`, `mPublicationValid`, and the three methods above) to the workflow adapter (`CombinedCATrackerSpec.cxx`), which already runs the per-event loop and already knows the combined outcome from `TrackingEngine::executeEvent()`'s own `EventResult`. |
| Owns the two `LegacyCATrackingParticipant<NLayers>` instances | **(3)**/pass-through | **Yes** | Unaffected in kind; the owned participants' *internal* member types and eventual name change per §3.3. |
| Caller readback (`getITSScratch()`/`getMFTScratch()`/`getITSSharedClusterCompatibility()`/etc.) | **(3)** adapter-private forwarding | **Yes** | Return types track whatever the renamed participant type's members become. |

**Verdict for this class, as currently shaped: fails the four-condition test.** The
publication/timing bridge (`mITSClock`/`mMFTClock`/`mPublicationValid` plus
`invalidatePublication()`/`markPublicationValid()`/`clearPublicationSidecars()`) is
concrete, cited, event-owned mutable state and reset/coordinator behavior living directly
inside this class today — not a hypothetical risk. §9's M6g slice is the explicit,
separately-scoped evaluation this correction requires: relocate that state/behavior to
`CombinedCATrackerSpec.cxx` (which already owns the event loop and already needs the
combined `EventResult` this bridge currently derives internally) and re-audit the
remaining class against the four conditions; if what remains is genuinely stateless
plan/config construction, retain it under a new, non-`Legacy` name; if any disqualifying
responsibility cannot be relocated without recreating a second event-owned/coordinator
layer elsewhere, delete it and inline its remaining construction logic directly into
`CombinedCATrackerSpec.cxx`'s own composition, per the explicit alternative this
correction specifies.

## 4. The generic workspace model

### 4.1 The one structural substitution

Every group in §3.1 that is genuinely `NLayers`-bound reduces to the **same single
substitution**, not two different mechanisms for "cluster cache" vs. "CA transients":

- Group B's outer containers are **already** sized from `mTrackingTopologyView.nCells`/
  `nTransitions` at runtime — themselves derived from `DetectorLayoutSet`/
  `SparseTrackingTopology` (`Architecture.md` §8: "contains only configured transitions
  and valid connected transition pairs"), which is **already dynamic and sparse-plan-based
  today**, for every detector, including the combined ITS+MFT case (`DetectorLayoutSet.h`
  §5-above). `NLayers` contributes nothing to their sizing; it only names the surrounding
  class template.
- Group A's per-layer `std::array<T, NLayers>` members carry exactly one slot per
  **owned surface**, i.e. `ownedSurfaces().size()` — a count the plan already exposes
  (`DetectorLayoutConfigurationKey::orderedSurfaces.size()`,
  `DetectorTraversalBinding::mLegacyLayerBySurface`'s own domain) and that a future
  ALICE 3 leg supplies identically, regardless of its own layer count.
- `TrackSeedTpl<NLayers>`/`SeedMetadataBase<NLayers>` (`Cell.h`) is the one true type-level
  leak: its `mClusters` is a `std::array<int, NLayers>`, one slot per owned surface,
  **already documented in its own source comment as "a temporary legacy boundary ...
  not a detector-selected state family."** `CellSeed` (`SeedMetadataBase<ClustersPerCell>`,
  a fixed constant) shows this was already solved once for cells; `TrackSeedTpl` is the
  one remaining instance of the same pattern.

So: **the smallest proposed generic workspace is a single non-templated
`SurfaceTrackingScratch`, sized once (at plan-adoption time, not per event) from two
runtime counts already available from the adopted plan/binding**:

```
nOwnedSurfaces = ownedSurfaces().size()             // replaces every std::array<T, NLayers>
{nTransitions, nCells} = topology.nTransitions/nCells  // already runtime today (Group B)
```

Every currently-`std::array<T, NLayers>`/`std::array<T, 2>`-inside-an-`NLayers`-template
member becomes a `std::vector<T>` (or `bounded_vector<T>`, matching today's allocator
discipline) sized once from `nOwnedSurfaces`; every currently-runtime-sized-but-templated-
container (Group B, Group D) keeps its existing sizing logic unchanged, just inside a
non-templated class.

### 4.2 The one exception: fixed-capacity, GPU-portable per-element types

`SeedMetadataBase<N>`/`CellSeed`/`TrackSeedTpl<NLayers>` are annotated `GPUhd()` throughout
`Cell.h` — they are used, or intended to be used, on device code paths, where heap
allocation (`std::vector`) is not available. **§4.1's "runtime-sized `std::vector`"
substitution does not apply to these per-element value types.** The generic replacement
for `TrackSeedTpl<NLayers>` (call it `TrackSeed`, mirroring `CellSeed`'s own existing
non-templated shape) must instead use a **single shared fixed-capacity constant**
(e.g. `kMaxSurfacesPerTrack`, chosen generously for every currently-known detector and
raised — as a plain constant, never a per-detector template parameter — if a future
ALICE 3 configuration needs more) plus a runtime "how many of these slots are actually
used" count, exactly the same pattern `CellSeed` already uses for its fixed 3-cluster
`ClustersPerCell` bound.

This is the one place "dynamic/sparse-plan-based rather than another `NLayers` template"
(the task's own framing) needs a caveat: dynamic at the **container** level (§4.1, sized
once from the plan, ordinary heap allocation, host-only), fixed-capacity at the
**per-element value-type** level (this section, GPU-portable, no allocation). Conflating
the two would either reintroduce a template (defeating the point) or break GPU
portability (a regression this note's own constraints forbid: "no ... GPU migration
decisions").

### 4.3 What stays out of `SurfaceTrackingScratch`

- `CATrackType<NLayers>` (`CATrackTypeHelper<NLayers>`) — genuinely detector-output-typed
  (`TrackITSExt`-shaped for ITS, MFT-specific for MFT). Category **(3)**, adapter-private,
  unaffected.
- `Tracker<NLayers>`/`TrackerTraits<NLayers>` themselves, and every internal `NLayers`-sized
  local (e.g. `LayerMeasurementSpans<NLayers>`) that lives inside `TrackerTraits<NLayers>`
  rather than inside the scratch — out of scope, per §3.3.
- The Group-A legacy cluster/index-table cache's **indexing scheme** (legacy-layer position
  vs. `SurfaceId`) is unchanged by M6 — only its array *bound* changes from `NLayers` to
  `ownedSurfaces().size()` (§3.1's note flags the future work; not part of M6).

## 5. Ownership between `TrackingEngine`, `TrackingParticipant`, `TimeFrame`, and adapters

Unchanged from ADR 0007 and unaffected by this migration — restated for completeness,
not revised:

- **`TrackingEngine`**: owns nothing per-detector; runs `track()` over an ordered
  `TrackingParticipant` schedule; owns the all-or-nothing failure/reset contract.
- **`TrackingParticipant`** (the interface): carries no loading operation, no ITS/MFT/
  `NLayers`/`SurfaceTrackingScratch` knowledge in its own public contract — unaffected,
  since `SurfaceTrackingScratch` lives inside *concrete* participants
  (`LegacyCATrackingParticipant<NLayers>` today), never in the interface.
- **`TimeFrame`**: unaffected by this migration. It already owns exactly the detector-
  neutral, permanent event state (normalized frame, `CommonTrack`/`TrackClusterReference`,
  vertices, beam/Bz). `SurfaceTrackingScratch` does not change what `TimeFrame` owns —
  it only changes what `LegacyTrackerScratch<NLayers>`'s successor owns.
- **Detector adapters** (`ITSMFTLegacyParticipantSet` today): own every genuinely
  detector-specific fact (§3.4) — unaffected in kind by this migration.

## 6. Atomic multi-source loading and reset semantics under the new workspace

`MultiSourceTimeFrameLoader::loadEvent()`/`LoadTarget`/`AtomicLoadBinding` (M2b, already
participant-count-generic — §1) are **unaffected in contract** by this migration.
`LoadTargetImpl<NLayers>` changes only its *staging target's type* (`SurfaceTrackingScratch`
instead of `LegacyTrackerScratch<NLayers>`), keeping the same stage-then-commit,
strong-exception-safety, same-allocator-swap discipline `loadNormalizedSource()` already
implements (§3.1 Group A). The whole-event transactionality
(`loadEvent()` stages every binding before committing any of them; a failure at any point
leaves `TimeFrame` and every live scratch untouched) is unchanged.

**Reset**: `TrackingEngine::resetEvent()`'s existing contract (reset every scheduled
participant, then wipe `TimeFrame` exactly once) is unaffected. One pre-existing wrinkle,
found during this audit and **not introduced by this migration**, is worth recording since
it bears on M6d/M6e's replay gates: `Tracker<NLayers>::clustersToTracks()`'s own
recoverable-failure path already calls `resetTimeFrameEvent(frame, scratch)`
(`LegacyTrackerScratch.h`) — which wipes the *shared* `TimeFrame` — **from inside a single
participant's own `track()`**, before `TrackingEngine::executeEvent()` ever sees the
non-Success outcome and calls its own `resetEvent()` (which wipes `TimeFrame` a second
time). `TimeFrame::wipe()` is idempotent on an already-empty owner, so this is not a
correctness bug — Gate 3's `DropTFUponFailure`/wipe/sentinel contract was accepted with
this exact behavior — but it is a documented departure from ADR 0007 decision 5's
layering ("a participant's `eventReset()` ... must never wipe ... `TimeFrame`'s storage"),
`LegacyTrackerScratch.h`'s own `resetTimeFrameEvent()` doc already flags it as "NOT the
future combined-owner policy." §12 (Rank 2) flags this as a **separate, optional simplification**
(stop wiping `TimeFrame` from inside `Tracker<NLayers>`'s own recoverable-drop path once
the engine's `resetEvent()` is the only caller that needs to), not mandatory M6 scope.

## 7. Replacing source-qualified bindings without leaking ITS/MFT

Per §3.2: the successor to `DetectorTraversalBinding` (call it `SurfacePlanBinding`) drops
the two detector-switched lines and instead accepts the caller's already-derived
`expectedKind`/`expectedPolicy` as parameters:

```
SurfacePlanBinding::build(globalLayout, source, ownedSurfaces, orderedSurfaces,
                          expectedKind, expectedPolicy)
```

The adapter (`ITSMFTLegacyParticipantSet`, or a future ALICE 3 equivalent) derives
`expectedKind`/`expectedPolicy` from its own plan's surface kinds — exactly the
information it already has when it calls `DetectorLayoutBuilder`/`buildDetectorLayoutSet()`
to build that plan in the first place, so this adds no new adapter-side state, only moves
an existing two-line decision from inside `build()` to just before the call. The fixed
`ClusterSourceId{0}`/`{1}` pairing was **already** isolated to
`ITSMFTLegacyParticipantSet::validateSources()`/`loadBindings()` before this note (M2c) —
not inside `DetectorTraversalBinding` — so no additional change is needed there.

## 8. Sidecars remain outside the generic workspace

`ITSSharedClusterCompatibility`/`MFTPublicationCompatibility` are unaffected by this
migration. They were never members of `LegacyTrackerScratch`/`DetectorTraversalBinding`
— they are owned directly by `LegacyCATrackingParticipant<NLayers>` via the
`*CompatibilityOwner<NLayers>` mixins (§3.3) and adopted into `TrackerTraits<NLayers>`,
exactly as today. **This note does not propose moving them into `SurfaceTrackingScratch`,
`TimeFrame`, or any new shared type** — they remain category **(3)**, adapter-private,
per the task's own explicit constraint ("do not move detector sidecars ... into
`TimeFrame`").

**Explicit clarification on the audit-goal's class list (corrected)**:
`LegacyCATrackingParticipant` and `ITSMFTLegacyParticipantSet` were named in the audit
scope alongside `LegacyTrackerScratch`/`DetectorTraversalBinding` for classification
purposes (§3.3, §3.4) — and, on that classification, *neither survives M6 under its
current name*. What is permanent is the **role** ADR 0007 decision 2 requires ("Detectors
are adapters" — every future detector, including ALICE 3, needs an analogous concrete
participant/adapter); it is not permanence of these two specific classes, as currently
named and shaped, and this note's earlier draft was wrong to conflate the two. The
corrected disposition (§3.3, §3.4, §9, §10):

- `LegacyCATrackingParticipant<NLayers>` is renamed at M6f, once its scratch/binding
  members are generic — the concrete-per-detector-participant *role* persists under a new,
  non-`Legacy` name.
- `ITSMFTLegacyParticipantSet` survives M6 only if M6g's relocation of its event-owned
  publication/timing/reset state leaves a class that demonstrably meets the four
  conditions in §3.4 (no `loadEvent`/`process`/`execute`/`reset` behavior, no event-owned
  mutable state, no coordinator role, no algorithm ownership beyond construction) — and is
  then renamed to drop `Legacy` too. If it does not meet those conditions after relocation,
  it is deleted and inlined into `CombinedCATrackerSpec.cxx`'s composition instead.

No `Legacy…`-named class may exist in production code once M6 (through M6g) completes —
see §11 for the explicit grep/deletion criteria this correction adds.

## 9. Staged implementation sequence

Mirrors the M5a→M5d granularity (design note, then bounded, replay-gated slices).

### M6a — this design/audit (current milestone)

No production code change. Deliverables: this note, plus the M6 section rewrite in
[GenericTrackingEngineMigration.md](../GenericTrackingEngineMigration.md).

### M6b — `SurfacePlanBinding` (first bounded implementation slice)

- **Scope**: add `SurfacePlanBinding` (§7) alongside `DetectorTraversalBinding`, unused by
  production. Copies `DetectorTraversalBinding::build()`'s body verbatim except: (a) drop
  the `detector != ITS && detector != MFT` gate; (b) accept `expectedKind`/`expectedPolicy`
  as parameters instead of switching on `detector`. Every other check (source validity,
  surface-mask/ownership, topology cross-boundary validation, scratch-slot assignment)
  copied unchanged. `DetectorTraversalBinding` itself is untouched and stays production-live.
- **Temporary bridge**: none additive — `SurfacePlanBinding` is not yet referenced by any
  production caller.
- **Acceptance/replay gate**: a focused test suite (ported from
  `testDetectorTraversalBinding*.cxx`) proves `SurfacePlanBinding::build()` produces
  byte-identical slot maps to `DetectorTraversalBinding::build()` for the current ITS/MFT
  combined layout, given the same `expectedKind`/`expectedPolicy` the old code derived
  internally. No production replay is affected (no production wiring yet).
- **Deletion/exit criterion**: none yet — `DetectorTraversalBinding` is deleted only at M6f.
- **Dependency**: none (this design note only).
- **Classification**: additive, behavior-preserving.

### M6c — `SurfaceTrackingScratch` and `TrackSeed`

- **Scope**: add the non-templated `SurfaceTrackingScratch` (§4.1) and non-templated
  `TrackSeed` (§4.2), unused by production. Sized from `SurfacePlanBinding`'s slot maps
  and `ownedSurfaces().size()`. Focused tests only (construction/sizing, reset, the
  allocator-swap discipline Group A's commit already requires) — mirrors M5a's "build the
  harness, don't wire it" pattern.
- **Implementation refinement (this milestone)**: `adoptPlan()` takes three plain runtime
  counts (`nOwnedSurfaces`, `nTransitions`, `nCells`) rather than a `SurfacePlanBinding`
  reference or including its (`detail/`-confined) header at all — the caller extracts
  `getOwnedSurfaces().count()`/`getGlobalTransitions().size()`/`getGlobalCells().size()`
  and passes them in. This is a stronger form of the decoupling `LegacyTrackerScratch<NLayers>`
  already has from `DetectorTraversalBinding` (neither type's header references the
  other's), not a deviation from it. Two further sizing inputs this section did not
  separately call out turned out to need their own scope decision, deferred rather than
  silently included: `IndexTableUtils<NLayers>` (Group A's index-table *binning
  configuration*, itself an `NLayers`-templated type with its own per-layer fixed arrays —
  genericizing it is out of this slice's authorized scope) and the `TrackingTopology<NLayers>`
  view objects (`mTrackerTopologies`/`mDefaultTrackingTopology`/`mVertexingTopology`/
  `mTrackingTopologyView`, also `NLayers`-templated) — Group B's `nTransitions`/`nCells`
  sizing is taken directly from the caller-supplied counts instead, which already
  supersedes needing an owned topology-view object for sizing purposes. Both are left for
  a later milestone to pick up explicitly, not silently dropped.
- **Temporary bridge**: none additive.
- **Acceptance/replay gate**: focused tests green; no production replay affected.
- **Deletion/exit criterion**: none yet.
- **Dependency**: M6b.
- **Classification**: additive, behavior-preserving.

### M6d — Wire MFT onto `SurfaceTrackingScratch`/`SurfacePlanBinding`

- **Scope**: `LegacyCATrackingParticipant<MFTNLayers>`'s `mScratch`/`mBinding` member types
  switch to the new types; `Tracker<MFTNLayers>`/`TrackerTraits<MFTNLayers>` (unchanged
  templates) adapt to the new scratch's accessor surface via the accessor seam the
  migration plan already anticipates. MFT first, per the existing migration-plan rationale
  (its refit is already fully normalized — [ADR 0008]). `LoadTargetImpl<MFTNLayers>`
  retargets its staging to the new scratch.
- **Temporary bridge**: the accessor seam both scratch types implement during this
  per-participant switch (both `LegacyTrackerScratch<MFTNLayers>` for ITS and
  `SurfaceTrackingScratch` for MFT coexist inside `ITSMFTLegacyParticipantSet` for the
  duration of M6d–M6e).
- **Acceptance/replay gate**: MFT standalone (`o2-mft-ca-tracker-workflow`) and combined
  (`o2-itsmft-combined-ca-tracker-workflow`) replays byte-identical to the M5d-era
  candidate baseline (MFT 68 tracks, hash `8106b08571ca593c6b76ff72b761a680` — [ADR 0008]
  addendum 2); ITS leg of the combined workflow untouched and still byte-identical to its
  own standalone replay (proving no cross-participant leak from the switch).
- **Deletion/exit criterion**: `LegacyCATrackingParticipant<MFTNLayers>` holds no
  `LegacyTrackerScratch<MFTNLayers>`/`DetectorTraversalBinding` member (grep-verified for
  the MFT instantiation only).
- **Dependency**: M6c.
- **Classification**: behavior-preserving cleanup (replay-gated).

**As-built correction (M6d implementation)**: the scope bullet's premise — that
`Tracker<MFTNLayers>`/`TrackerTraits<MFTNLayers>` stay "unchanged templates" and adapt via
some accessor-level seam alone — did not hold. Both class templates store `ScratchN`
(`= LegacyTrackerScratch<NLayers>`) and, for `TrackerTraits`, `const DetectorTraversalBinding*`
as their own internal type aliases/members, not as a caller-supplied, duck-typed parameter;
there was no accessor surface to retarget without also retargeting these types themselves.
The actual mechanism: `TrackerTraits`, `Tracker`, `LegacyCATrackingParticipant`,
`DetectorTraits`, and `refitTrackFwd` each gained two **defaulted** extra template
parameters, `ScratchT = LegacyTrackerScratch<NLayers>` and (where a binding is held)
`BindingT = DetectorTraversalBinding`. Every pre-existing instantiation (ITS; at the time of
this M6d correction, also the standalone MFT workflow's own
`ITSMFTTrackingInterface<MFTNLayers>`, untouched by M6d itself and still using
`TrackerTraits<10>`/`Tracker<10>` by their default args — **see M6e1 below, which migrates
this exact class next**) is therefore bit-for-bit unaffected by the new parameters'
existence, while the combined-workflow MFT participant explicitly instantiates
`TrackerTraits<10, SurfaceTrackingScratch, SurfacePlanBinding>` and
`Tracker<10, SurfaceTrackingScratch, SurfacePlanBinding>`. This preserves the milestone's own
intent (MFT switches, ITS/standalone-MFT stay on the legacy types, no per-call dispatch cost
in the hot loop) through a type-level rather than accessor-level seam. The temporary-bridge,
acceptance/replay-gate, and deletion/exit-criterion bullets above were all met as stated;
only the *mechanism* description needed correcting. Actual replay results: MFT standalone
68 tracks/hash `8106b08571ca593c6b76ff72b761a680`; ITS standalone 212 tracks/hash
`46913a67a7e2fe7462e29df0db264fa8`; combined workflow's ITS and MFT legs each bit-identical
to their own standalone replay — all four exactly matching [ADR 0008] addendum 2.

### M6e1 — Wire the standalone MFT interface onto `SurfaceTrackingScratch`/`SurfacePlanBinding`

Inserted after M6d and before the original M6e (renamed in no way at the time — M6e below
remains exactly what it always was, the ITS switch plus Group-C legacy-result-staging
retirement; M6e1 is a new, narrower slice this note did not originally anticipate). M6e itself
was later split in two (M6e2/M6e3, see M6e2's own note) once the ITS switch alone turned out to
need its own gated milestone; that split happened after this section was written and does not
change anything stated here. Reason for the insertion: M6d's
own scope was explicitly the *combined-workflow* MFT participant
(`LegacyCATrackingParticipant<MFTNLayers>`, driven by `ITSMFTLegacyParticipantSet`); M6d's own
audit never read `TrackingInterface.cxx` at all, since that file backs a structurally
different, independent owner — `ITSMFTTrackingInterface<MFTNLayers>`, which drives the
*standalone* MFT common-CA workflow (`o2-mft-ca-tracker-workflow`) and owns its own
`LegacyTrackerScratch<MFTNLayers>` directly, sharing only the `Tracker<10>`/`TrackerTraits<10>`
class templates with the combined participant (via the M6d seam's default template
arguments). After M6d, this second, independent live-production MFT path was still entirely
on the legacy types — a real gap, not a design omission this note had already scoped out.

- **Scope**: apply the identical M6d seam to `ITSMFTTrackingInterface<NLayers, ScratchT,
  BindingT>` (two new defaulted template parameters, same defaults). `ITSMFTTrackingInterfaceMFT`
  becomes the explicit `<MFTNLayers, SurfaceTrackingScratch, SurfacePlanBinding>` instantiation;
  `ITSMFTTrackingInterfaceITS` keeps the defaults, unaffected. `initialiseTracker()` builds and
  adopts a `SurfacePlanBinding` for MFT only (`ClusterSourceId{0}`, `SurfaceKind::Disk`/
  `TransitionPolicyTag::DiskDisk` as literal adapter constants — the identical
  adapter-derives-kind/policy pattern M6d already established, no detector switch added to
  `SurfacePlanBinding` itself) and calls `SurfaceTrackingScratch::adoptPlan()` (mirroring
  `LegacyCATrackingParticipant::adoptDetectorLayoutSet()`'s own identical step). `fillMFTOutputs()`
  (`CATrackerSpec.cxx`, MFT workflow) retargets from `LegacyTrackerScratchMFT` to
  `SurfaceTrackingScratch`.
- **Real gap in `SurfaceTrackingScratch` found and closed**: M6d's own claim that
  `initVertexingTopology()`/`initDefaultTrackingTopology()` have "zero production callers" was
  correct only for the files M6d's own audit actually read
  (`TrackerTraits.cxx`/`CATracker.cxx`/`LegacyCATrackingParticipant.cxx`) — never
  `TrackingInterface.cxx`, which calls `initDefaultTrackingTopology()` unconditionally once per
  event, for both ITS and MFT. `SurfaceTrackingScratch` gained this method (plus
  `initVertexingTopology()`, ported for structural parity even though it remains genuinely
  unreferenced) and the corresponding `mDefaultTrackingTopology`/`mVertexingTopology` members
  before this slice could work at all. General lesson for future slices: a prior milestone's
  "zero callers" claim is scoped to the files that milestone actually read, not to the whole
  codebase, unless it says so explicitly.
- **Temporary bridge**: none additive beyond what M6d already introduced — the same
  `ScratchT`/`BindingT` seam, now also instantiated by this second owner. The bare,
  default-argument `ITSMFTTrackingInterface<10>` stays explicitly instantiated in
  `TrackingInterface.cxx` solely because `testTrackingInterfaceLoadFailureContract.cxx`'s own
  pre-existing `NLayers`-generic load-failure-contract coverage links against it; this is not a
  second production abstraction layer, and no production caller reaches it any more.
- **Acceptance/replay gate**: standalone MFT (`o2-mft-ca-tracker-workflow`, MFT
  CommonTrack-output route enabled) replays exactly 68 tracks, hash
  `8106b08571ca593c6b76ff72b761a680`; standalone ITS (unaffected) exactly 212 tracks, hash
  `46913a67a7e2fe7462e29df0db264fa8`; combined workflow's MFT and ITS legs each bit-identical
  to their own standalone replay — all four exactly matching [ADR 0008] addendum 2 and M6d's
  own recorded numbers. `ctest -L itsmft`: 95/95 (up from M6d's 94). Fixture's 43 checksums
  unchanged before and after every replay.
- **Deletion/exit criterion**: `ITSMFTTrackingInterfaceMFT` holds no `LegacyTrackerScratch<MFTNLayers>`/
  `DetectorTraversalBinding` (grep-verified); combined MFT participant unaffected by this slice
  (already migrated in M6d).
- **Dependency**: M6d.
- **Classification**: behavior-preserving cleanup (replay-gated).

**State after M6e1**: combined MFT and standalone MFT both use `SurfaceTrackingScratch`/
`SurfacePlanBinding`; ITS (standalone and combined) still uses `LegacyTrackerScratch<7>`/
`DetectorTraversalBinding` entirely unchanged. Both storage models now coexist purely along
the ITS/MFT axis, not along any combined/standalone axis — the state M6e (next) inherits.
`mTracks`/`mTracksLabel` (Group C) were not retired by this slice; that remains M6e's own
scope, unstarted here.

### M6e2 — Wire both live ITS common-CA paths onto `SurfaceTrackingScratch`/`SurfacePlanBinding`

Inserted after M6e1, splitting the original M6e milestone in two — exactly the same kind of
insertion M6e1 itself performed on M6d/M6e (see that section's own note). Reason: M6e1's own
"State after M6e1" summary explicitly deferred the ITS switch and the Group C
(`mTracks`/`mTracksLabel`) retirement to "M6e (next)" as one combined slice; in practice the
ITS switch alone required enough new architectural work (below) to justify treating it as its
own gated milestone, leaving Group C retirement as a narrower follow-up (renamed M6e3, next
section) that starts from an already-fully-migrated ITS/MFT storage baseline instead of doing
both at once.

- **Scope**: apply the identical M6d/M6e1 seam to both remaining live ITS paths —
  `LegacyCATrackingParticipant<ITSNLayers>` (combined workflow) and
  `ITSMFTTrackingInterface<ITSNLayers, ScratchT, BindingT>` (standalone
  `o2-its-ca-tracker-workflow`). `LegacyCATrackingParticipantITS`/`ITSMFTTrackingInterfaceITS`
  become the explicit `<ITSNLayers, SurfaceTrackingScratch, SurfacePlanBinding>`
  instantiation; both bare, default-argument templates stay explicitly instantiated solely
  for `testTrackingInterfaceLoadFailureContract.cxx`'s own `NLayers`-generic coverage,
  unchanged from M6e1's own precedent. Each ITS participant/interface builds and adopts its
  own `SurfacePlanBinding` (`ClusterSourceId{0}`, `SurfaceKind::Cylinder`/
  `TransitionPolicyTag::CylinderCylinder` as literal adapter constants — no detector switch
  added to `SurfacePlanBinding` itself) and calls `SurfaceTrackingScratch::adoptPlan()` before
  any load attempt, mirroring M6d/M6e1's own identical step. `fillITSOutputs()`
  (`CATrackerSpec.cxx`, ITS workflow) retargets from `LegacyTrackerScratch<ITSNLayers>` to
  `SurfaceTrackingScratch`. `mTracks`/`mTracksLabel` (Group C) legacy result staging is
  **not** retired by this slice — both participants keep populating it in parallel with
  `AcceptedTrackShadowPublisher`/CommonTrack, exactly as M6d/M6e1 left MFT; that retirement is
  M6e3 (next section), unstarted here.
- **Three architectural gaps found and closed, all stemming from the same root cause**:
  `SurfaceTrackingScratch` had only ever actually served MFT in production (M6d/M6e1), so
  several of its members were still hardcoded to MFT's own shape in ways nothing had
  exercised for a second, differently-sized detector until this slice's ITS(7) instantiation
  actually compiled and ran:
  1. **`IndexTableUtils<nLayers>`** (§4.1's already-flagged deferred item) was a genuine class
     template — `IndexTableUtils<7>` and `IndexTableUtils<10>` were different C++ types, so a
     single `SurfaceTrackingScratch::IndexTableUtilsN` member could not serve both detectors
     at once, and neither could the (then-forward-declared) `IndexTableUtils<NLayers>` used
     by `detail/TransitionPolicyOperations.h`. Fixed by de-templating `IndexTableUtils<nLayers>`
     into a template *alias* for one shared, non-templated `IndexTableUtilsCore` class using
     `MaxLayoutSurfaces`-bounded (`SurfaceId.h`) fixed-capacity storage — the same established,
     reused GPU-portable-storage precedent `TrackSeed` (M6c, §4.2) already set, not a new one.
     `LegacyTrackerScratch<NLayers>::IndexTableUtilsN`'s own source line is unchanged and still
     compiles identically (now resolving to the same shared `IndexTableUtilsCore`, invisibly).
     `detail/TransitionPolicyOperations.h`'s stale forward declaration (`template <int NLayers>
     class IndexTableUtils;`, invalid once `IndexTableUtils<N>` became an alias template — alias
     templates cannot be forward-declared) was replaced with a real include.
  2. **Group C (`getTracks()`/`getNumberOfTracks()`)** is genuinely detector-output-typed
     (`CATrackType<7>` = `TrackITSExt`, `CATrackType<10>` = typed MFT output) and could not be
     alias-erased the way IndexTableUtils was. Fixed with the same dual-storage pattern this
     note already anticipated in principle (§3.1): `SurfaceTrackingScratch` now stores both
     `mTracksITS`/`mTracksMFT` (exactly one ever populated per instance) and exposes them via a
     compile-time-selecting *template* accessor `getTracks<NLayers>()` — never virtual/
     type-erased. Because `LegacyTrackerScratch<NLayers>::getTracks()` is, and must remain, a
     plain non-template method, the temporary `ScratchT` seam's own Group-C call sites inside
     the shared `Tracker<NLayers,ScratchT,BindingT>`/`TrackerTraits<NLayers,ScratchT,BindingT>`
     bodies (reachable with either scratch type) cannot use the same call-site text for both —
     resolved via small free dispatcher functions (`scratchTracks<NLayers>()`/
     `scratchNumberOfTracks<NLayers>()`) doing a narrow `if constexpr` on `ScratchT`'s identity.
     This selects on *scratch representation*, exactly what the seam already exists to select,
     never on ITS-vs-MFT identity, and collapses to a plain forwarding call once M6f removes
     `LegacyTrackerScratch<NLayers>` as anyone's `ScratchT`.
  3. **The ROF-overlap/ROF-vertex-lookup/ROF-mask/tracking-topology auxiliaries**
     (`ROFOverlapTable<N>`, `ROFVertexLookupTable<N>`, `ROFMaskTable<N>`,
     `TrackingTopology<N>` — §4.1's same deferred item) were not caught by M6d/M6e1's own
     "narrow, deliberate exception... this scratch is, for now, only ever used by MFT" framing
     being re-verified against ITS, because that framing's premise (MFT-only) was true at the
     time it was written and nothing had yet exercised the ITS(7) case. These types are
     genuinely `N`-sized (fixed-capacity internal storage — e.g.
     `TrackingTopology<N>::MaxTransitions = N*(N-1)/2` is 21 for ITS vs. 45 for MFT), so a
     single MFT(10)-shaped instance of any of them is wrong for ITS(7): not just a compile-time
     mismatch (which is what actually surfaced first, once the ITS(7) instantiation was
     compiled) but a real correctness bug had it compiled silently — wrong transition/cell
     counts, wrong ROF-layer array width. Fixed with the identical Group C dual-storage +
     template-accessor + free-dispatcher pattern, extended to these four groups (dual
     `mROFOverlapTableITS`/`MFT`, `mROFVertexLookupTableITS`/`MFT`,
     `mMultiplicityCutMaskITS`/`MFT` + `mUPCCutMaskITS`/`MFT` + their own active-pointer
     selectors, `mTrackerTopologiesITS`/`MFT` + `mDefaultTrackingTopologyITS`/`MFT` +
     `mVertexingTopologyITS`/`MFT`); `initialise()`/`prepareClusters()` (which read/write
     several of these) became templates on `NLayers` too, and `loadNormalizedSource()`'s own
     detector preflight (previously narrowed to MFT-only, matching the same now-outdated
     premise) was widened back to ITS-or-MFT, matching `LegacyTrackerScratch<NLayers>`'s own
     long-standing check — including restoring the ITS branch of its synthetic
     `TrackingFrameInfo` backfill (`m.frame.q/frameAngle/u/v`, ported byte-for-byte from
     `LegacyTrackerScratch::loadNormalizedSource()`), which the MFT-only version had never
     needed. Setters taking an `NLayers`-typed argument (`setROFOverlapTable()`, etc.) needed
     no dispatcher shim: ordinary template-argument deduction from the argument works
     identically whether the callee resolves to `SurfaceTrackingScratch`'s template method or
     `LegacyTrackerScratch<NLayers>`'s plain one, with no explicit `<>` and therefore no
     `.template` disambiguator at the call site.
- **Temporary bridge**: none additive beyond M6d/M6e1's own `ScratchT`/`BindingT` seam, now
  instantiated by ITS's two remaining owners too.
- **Acceptance/replay gate**: standalone ITS (`useDiamond=true`, diamond `(0,0,0)`,
  `pvRes=0.05`) replays exactly 212 tracks, hash `46913a67a7e2fe7462e29df0db264fa8`; standalone
  MFT (CommonTrack-output route) exactly 68 tracks, hash `8106b08571ca593c6b76ff72b761a680`;
  combined workflow's ITS and MFT legs each bit-identical to their own standalone replay — all
  four exactly matching [ADR 0008] addendum 2 and M6d/M6e1's own recorded numbers.
  `ctest -L itsmft`: 96/96 (up from M6e1's 95). Fixture checksums unchanged before and after
  every replay. `LegacyTrackerScratch.h`/`.cxx` and `detail/DetectorTraversalBinding.h`
  remain byte-for-byte zero-diff against the M6e1 base.
- **Deletion/exit criterion**: neither `LegacyCATrackingParticipantITS` nor
  `ITSMFTTrackingInterfaceITS` holds `LegacyTrackerScratch<ITSNLayers>`/
  `DetectorTraversalBinding` (grep/type-verified); no detector-specific switch reintroduced
  into `SurfacePlanBinding` (grep-verified, `testSurfacePlanBindingNoDetectorDependency.cxx`).
- **Dependency**: M6e1.
- **Classification**: behavior-preserving cleanup (replay-gated).

**State after M6e2**: standalone ITS, combined ITS, standalone MFT, and combined MFT all use
`SurfaceTrackingScratch`/`SurfacePlanBinding` — every live common-CA path is on the new storage
now, along both axes (detector and standalone/combined). `mTracks`/`mTracksLabel` (Group C)
were not retired by this slice; that remains M6e3's own scope, unstarted here. The temporary
`ScratchT`/`BindingT` seam itself also remains, unremoved until M6f.

### M6e3 — Retire legacy result staging (Group C)

Renamed from the original "M6e" (see M6e2's own note above for why): the ITS-switch half of
that milestone's original scope is now M6e2, already integrated; this section covers only
what remains.

- **Scope**: both participants' `mTracks`/`mTracksLabel` legacy result staging is
  retired — detector adapters (`ITSMFTLegacyParticipantSet`'s caller,
  `CombinedCATrackerSpec.cxx`, and the standalone `CATrackerSpec.cxx` workflows) build
  detector-typed output from `TimeFrame::getCommonTracks()`/`getTrackClusterIndices()` plus
  `ParticipantPublicationExport`/`CommonTrackPublicationExport` alone (already populated in
  parallel today via `AcceptedTrackShadowPublisher`, per §3.1 Group C — this slice removes
  the now-redundant legacy-typed copy, not adds new publication logic).
- **Temporary bridge**: none — this milestone does not touch the `ScratchT`/`BindingT` seam
  itself (that is M6f's own scope).
- **Acceptance/replay gate**: ITS and MFT standalone and combined replays byte-identical to
  the M6e2-era baseline (ITS 212 tracks, hash `46913a67a7e2fe7462e29df0db264fa8`; MFT 68
  tracks, hash `8106b08571ca593c6b76ff72b761a680` — [ADR 0008] addendum 2); **writer-level**
  output (not just internal track-count/hash) verified identical before/after the
  `mTracks`/`mTracksLabel` retirement, since this is the one part of this slice that changes
  an output-construction code path rather than only a container type.
- **Deletion/exit criterion**: no production reference to `CATrackType<NLayers>` inside
  `LegacyTrackerScratch`/its successor (the type itself survives, owned only by adapters
  building final output, per §4.3).
- **Dependency**: M6e2.
- **Classification**: behavior-preserving cleanup (replay-gated); output-construction path
  changes, so the writer-level check above is load-bearing, not optional.

**M6e3 validation record (2026-08-05)**: completed on
`codex/itsmft-m6e3-commontrack-output-retirement` with the pinned
`daily-20260717-0700-local1` environment.  A detached parent at
`1c277d2fda^` (`8b39ac8cfd`, already on the CommonTrack writer route) and
the M6e3 feature build both ran fresh standalone ITS and MFT replays.  The
field-level writer comparator reported exact ITS equality and exact MFT
non-projected equality; all 2,992 MFT forward-state/covariance/chi2 values
met the approved float-projection contract (observed max absolute and
relative deltas: zero).  The feature combined replay was field-equivalent to
each feature standalone leg.  ITS produced 212 tracks with
`46913a67a7e2fe7462e29df0db264fa8`; MFT produced 68 tracks with
`8106b08571ca593c6b76ff72b761a680`.  The read-only 43-file fixture checksum
manifest passed both before and after the campaign, and serial
`ctest -L itsmft --output-on-failure -j1` executed and passed all 96
registered tests.  The former output-selector parameter had no remaining
production, runtime-registration, or documentation reference, so its
removal is intentional rather than a retained dead toggle.  M6e3 deletes
only live `SurfaceTrackingScratch` Group C staging; `LegacyTrackerScratch`
and its Group C members remain intentionally untouched for M6f.

### M6f — Retire the temporary workspace/binding bridge

- **Scope**: delete `LegacyTrackerScratch<NLayers>` (header + source),
  `DetectorTraversalBinding` (`detail/DetectorTraversalBinding.h`), and the
  `TrackSeedTpl<NLayers>` whole-track representation. Retain `SeedMetadataBase<N>`
  because live `CellSeed` still requires it; do not delete that shared base merely
  because its former `TrackSeedTpl` user is gone. Rename
  `LegacyCATrackingParticipant<NLayers>` (and its `ITS`/`MFT` aliases and
  extern-template instantiations) to `SurfacePlanTrackingParticipant<NLayers>`.
  `ITSMFTLegacyParticipantSet` is **not** in this slice's scope — it is evaluated
  and disposed of separately at M6g.
- **Temporary bridge**: none — this milestone removes M6d/M6e's accessor-seam bridge.
- **Acceptance/replay gate**: every replay from M6d/M6e re-verified; `ctest -L itsmft`
  green; a dedicated guard (mirroring ADR 0008's `testNoLegacyFittingDependency.cxx`
  pattern; see §11) asserts zero occurrences of the retired workspace, binding,
  participant, whole-track seed, and compatibility-template spellings in
  `common/tracking/{include,src}`. The guard narrowly excludes frozen legacy ITS code
  outside that tree and permits only the unchanged `ITSMFTLegacyParticipantSet` M6g
  adapter name.
- **Deletion/exit criterion**: no production instantiation of
  `LegacyTrackerScratch<NLayers>` (grep-verified) — satisfies ADR 0007 decision 9 and
  the migration plan's own M6 exit criterion; `SurfaceTrackingScratch` has executed
  production traffic for both participants; zero old M6f bridge names remain in the
  common production include/src surface.
- **Dependency**: M6e3.
- **Classification**: behavior-preserving cleanup (replay-gated) for the deletions; the
  participant rename is mechanical (type identity only, no behavior change) but still
  replay-gated to prove nothing else moved. Any residual output delta requires separate
  approval under M5's decision (ADR 0007 decision 11/ADR 0008), same as the pre-existing
  migration-plan text for M6.

**M6f completion record (2026-08-05)**: the implementation used the integrated M6e3
head `ebcf6fc7608f2eb7562c248da008d799e3975515` on
`codex/itsmft-m6f-retire-legacy-workspace`. `SurfacePlanTrackingParticipant` is now
the narrow plan-driven participant name. The old workspace/binding headers, sources,
explicit instantiations, CMake entries, and binding-only tests are gone; useful
coverage was migrated to `SurfaceTrackingScratch`, `SurfacePlanBinding`, and
`TrackSeed`. The common-CA production guard covers only `include` and `src`, with
the frozen-ITS and M6g adapter exclusions stated above. `loadROFrameData`,
`resetROFrameData`, and `prepareROFrameData` had zero common production callers and
were removed with the old scratch; detector-specific raw-ROF workflow ownership and
compatibility state remain outside this deletion.

The durable build rebuilt the affected library, tests, and workflows. All 95
registered serial ITS/MFT tests passed. The fixture checksum manifest passed 43/43
before and after. Standalone and combined replays were mutually field-identical and
produced the accepted candidates: ITS 212 / `46913a67a7e2fe7462e29df0db264fa8`,
MFT 68 / `8106b08571ca593c6b76ff72b761a680`. Parent-output comparison matched all
initialized leaves and the approved MFT float-projected leaves; its sole mismatch
was the pre-existing uninitialized `TrackMFT::mInvQPtSeed` member (the parent and
M6f processes contained different run-dependent bytes). M6f did not rewrite that
compatibility field or alter the CommonTrack staging path, so the undefined legacy
value is recorded as a validation limitation rather than claimed byte equality.

### M6g — Evaluate, relocate, and dispose of `ITSMFTLegacyParticipantSet`

- **Scope/result**: the audit confirmed that the set's clocks, publication-validity
  flag, sidecar invalidation/clearing, ROF-table configuration, and combined-event
  publication bridge were event-owned behavior. They now live in
  `CombinedCATrackerDPL`/`CombinedCATrackerSpec.cxx`, which owns the DPL event lifecycle.
  The remaining combined layout, bindings, concrete ITS/MFT participants, and explicit
  `[ITS, MFT]` schedule were inlined into that task's construction path. No immutable
  plan holder remained that justified retaining or renaming the set.
- **Temporary bridge**: none — this was the final disposition, not an intermediate step.
- **Acceptance/replay gate**: the combined workflow was rebuilt against the post-M6f
  baseline; direct participant-composition tests cover load failure, dropped TF,
  structural tracking failure, successful replacement, publication invalidation, and
  source-qualified exports. The M6g guard asserts zero references to
  `ITSMFTLegacyParticipantSet` and the other retired bridge names in common production
  include/src and combined workflow include/src.
- **Deletion/exit criterion**: satisfied — the class/header/source/CMake entry are gone,
  and `CombinedCATrackerDPL` owns its former application composition directly. The
  generic tracking library retains no event-loop, timing/publication, or DPL source-order
  coordinator.
- **Dependency**: M6f.
- **Classification**: behavior-preserving cleanup (replay-gated).

## 10. "Not safe to delete yet" (M6-era refinement)

Supersedes the corresponding rows of
[GenericTrackingEngineMigration.md](../GenericTrackingEngineMigration.md)'s existing table
for this milestone; other rows (`TransitionPolicyTag` machinery, output sidecars) are
unchanged by this note and not repeated here.

| Artifact | Why it must stay | Removal gate |
|---|---|---|
| `LegacyTrackerScratch<NLayers>` | Retired after `SurfaceTrackingScratch` executed production traffic for **both** participants | **Deleted M6f** |
| `DetectorTraversalBinding` | Its plan-sizing/slot-assignment role is carried by `SurfacePlanBinding` for both participants | **Deleted M6f** |
| `TrackSeedTpl<NLayers>` | Retired whole-track representation; `TrackSeed` is now the sole common-CA whole-track seed | **Deleted M6f**; `SeedMetadataBase<N>` remains for `CellSeed` |
| `mTracks`/`mTracksLabel` (`CATrackType<NLayers>` legacy result staging) | Removed from live `SurfaceTrackingScratch` at M6e3 | **Already deleted M6e3** |
| `loadROFrameData()`/`resetROFrameData()`/`prepareROFrameData()` | Zero common production callers; detector-specific raw-ROF workflow APIs are a separate ownership path | **Deleted from common scratch M6f**; raw-ROF APIs remain outside common CA |
| `SurfacePlanTrackingParticipant<NLayers>` | Narrow plan-driven ITS/MFT participant adapter; no retired bridge ownership | **Remains**; it owns participant-local scratch/sidecars, not event-loop state |
| `ITSMFTLegacyParticipantSet` | Coordinator-shaped holder of combined application construction and event-owned publication/reset state | **Deleted M6g**; construction is inlined into the combined DPL task |

## 11. Grep/deletion criteria for `Legacy`-prefixed names

Every `Legacy`-containing production type identifier under
`Detectors/ITSMFT/common/tracking/{include,src}`, enumerated directly
(`grep -rhoE "\bLegacy[A-Za-z0-9_]*" include/ src/`), not assumed:

| Identifier(s) | Disposition | Excluded from the M6 guard? |
|---|---|---|
| `LegacyTrackerScratch`, `LegacyTrackerScratchITS`, `LegacyTrackerScratchMFT` | Deleted at M6f, replaced by `SurfaceTrackingScratch` | No — covered |
| `LegacyCATrackingParticipant`, `LegacyCATrackingParticipantITS`, `LegacyCATrackingParticipantMFT`, and the two `extern template class LegacyCATrackingParticipant<...>` instantiations | Renamed at M6f (§3.3, §9) | No — covered |
| `ITSMFTLegacyParticipantSet` | Deleted and inlined at M6g after the four-condition audit (§3.4, §9) | No — covered |
| `LegacyId`, `LegacyIndexMismatch` (`TrackerTraits.h`/`.cxx`) | Internal to `TrackerTraits<NLayers>`'s own hot-loop/topology-consistency machinery — a topology-ID type alias and a `TraversalFailureReason` enumerator, not a scratch/binding/participant-set artifact | **Yes** — `TrackerTraits<NLayers>`/`Tracker<NLayers>` stay `NLayers`-templated algorithm orchestration, explicitly out of M6's scope (§3.3, §4.3); renaming these belongs to whatever future milestone, if any, revisits `TrackerTraits`/`Tracker` naming, not M6 |
| `LegacyMaterialMismatch` (`mLayerMaterial`/`LegacyMaterialMismatch` et al.) | Already tracked in `GenericTrackingEngineMigration.md`'s own "Policy/legacy compatibility code" row, gated by its own M4–M6 slice | **Yes** — tracked there, not duplicated in this note's guard to avoid two sources of truth for the same item |

**Guard test**: M6g's acceptance gate adds a dedicated source/dependency guard (in the
style of ADR 0008's `testNoLegacyFittingDependency.cxx`) asserting zero occurrences of
`LegacyTrackerScratch`, `DetectorTraversalBinding`, `LegacyCATrackingParticipant`,
`TrackSeedTpl`, `ScratchT`, `BindingT`, or `ITSMFTLegacyParticipantSet` in every
production header/source under `Detectors/ITSMFT/common/tracking/{include,src}` and
the combined workflow include/src. Frozen legacy ITS code outside that common tree is
explicitly out of scope. Test files and historical documentation are not scanned. The
guard also checks that publication/timing state is declared in the workflow owner and
not in the generic tracking core.

## 12. Deletion/simplification opportunities flagged separately (not mandatory M6 scope)

Per the task's own instruction to flag these separately rather than silently broaden
scope. Ranked by how confidently actionable each is from this audit alone (Rank 1 =
grep-verified and mechanical; Rank 3 = needs more evidence before any action):

1. **Rank 1 — detector-specific raw-ROF APIs**: the common scratch methods
   `loadROFrameData()`/`resetROFrameData()`/`prepareROFrameData()` had zero common
   production callers and were deleted in M6f. The similarly named ITS/MFT raw-ROF
   workflow ownership APIs remain live and are intentionally outside this deletion.
2. **Rank 2 — `Tracker<NLayers>`'s own `TimeFrame::wipe()` call inside its
   recoverable-drop path** (§6): currently double-wipes `TimeFrame` on a recoverable
   single-participant failure (once inside `Tracker<NLayers>::clustersToTracks()`, once
   more inside the engine's `resetEvent()`). Idempotent today (not a correctness bug, not
   blocking M6), but a layering cleanup the codebase's own `resetTimeFrameEvent()` doc
   already anticipates ("NOT the future combined-owner policy"). Worth its own small,
   separately-decided slice once `TrackingEngine`/`ITSMFTLegacyParticipantSet` (or its
   M6g successor) are the only production callers of `Tracker<NLayers>::clustersToTracks()`
   (already true today, per §3.3) — not actioned here. Ranked below item 1 because it
   requires reasoning about call-site ordering, not just a dead-code grep.
3. **Rank 3 — `mNTotalLowPtVertices`** (former scratch protected member): declared
   and zero-initialized; no production read/write site found in this audit's
   grep pass. Flag for verification (not asserted dead — a wider grep across
   `ITStracking`/vertexer call sites than this audit performed would be needed before
   deleting it). Ranked last because, unlike items 1–2, this audit did not perform the
   wider grep that would confirm it.

None of these three block M6a's GO verdict or the M6f slice; none are included in
§9's acceptance gates.

## 13. Explicit non-actions (recap)

This note does not: modify the `Propagator` or its covariance sanitizer; change candidate
physics or run fixture generation; redefine any baseline; redesign any
tracklet/cell/road/refit algorithm or formula; move detector sidecars, raw ROFs, or
workflow output adaptation into `TimeFrame`; introduce a generic catalog copy, persisted
`SurfaceKindPair`, or public transition-policy API; begin M6g's
`ITSMFTLegacyParticipantSet` disposition; or touch vertexing
algorithms. `ctest -L itsmft` and the existing accepted/candidate replay baselines are
therefore unaffected by this milestone; no new replay was required or run.

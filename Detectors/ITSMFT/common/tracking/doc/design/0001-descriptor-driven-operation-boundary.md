# Design note 0001: descriptor-driven tracklet/cell/refit operation boundary

Status: M5a design, no production behavior change
Date: 2026-08-03
Scope: `Detectors/ITSMFT/common/tracking` (`o2::itsmft::tracking`)
Companion: [ADR 0007](../decisions/0007-generic-tracking-engine-boundary.md) decision 10
(one algorithmic flow for cylinders and disks); milestone M5 of
[GenericTrackingEngineMigration.md](../GenericTrackingEngineMigration.md).
Companion harness: `Detectors/ITSMFT/common/tracking/test/testNativeCylinderCylinderRefitDriver.cxx`
(existing, Gate 3 Slice B) and
`Detectors/ITSMFT/common/tracking/test/testITSNativeVsLegacyRefitCharacterizationHarness.cxx`
(new, M5a — see [Section 6](#6-ab-harness-summary)).

This note is deliberately narrower than ADR 0007: it does not revisit the
end-state boundary (permanent `TimeFrame`/`GenericTrack`, one `TrackingEngine`,
detectors as adapters) or the milestone sequencing. It only answers one
question ADR 0007 decision 10 leaves open: *exactly which functions, at
exactly what granularity, are the "narrow, explicit surface-operation
boundaries" family differences must be confined to*, and which of today's
`TransitionPolicyTag`-dispatched code is already at that boundary versus
still duplicated orchestration or legacy-only compatibility.

**M5a delivers this classification and the A/B refit harness only. It does
not activate native refit, does not change any cell/tracklet/refit formula,
does not reintroduce `TransitionPolicyTag` publicly, does not replace
`LegacyTrackerScratch`, and does not touch workflows/defaults/output
contracts.** Everything below is a design/characterization result, gated for
implementation by a separate, later M5 decision.

## 1. Method

Every claim below is grounded in the current source, not asserted from the
architecture docs alone:
- `Detectors/ITSMFT/common/tracking/include/ITSMFTTracking/detail/TransitionPolicyOperations.h`
  (the existing per-tag operation declarations and their documentation of
  which parts are shared vs. family-specific formulas);
- `Detectors/ITSMFT/common/tracking/include/ITSMFTTracking/{Barrel,Forward}SurfaceStateOperations.h`
  (the `barrel::`/`forward::` primitive signatures actually available);
- `Detectors/ITSMFT/common/tracking/src/TrackerTraits.cxx` (the orchestration
  bodies that call the above, read directly rather than assumed identical);
- `Detectors/ITSMFT/common/tracking/include/ITSMFTTracking/NativeCylinderCylinderRefitDriver.h`,
  `RefitLegAssembly.h` (the existing unwired native refit driver);
- `Detectors/ITSMFT/ITS/tracking/include/ITStracking/TrackHelpers.h` (the
  frozen legacy `fitTrack`/`refitTrack`/`refitTrackSeed` chain) and
  `Detectors/ITSMFT/common/tracking/src/DetectorTraits.cxx` (`refitSeedITS`,
  the production entry point into that chain).

## 2. Verdict

**The descriptor-driven operation boundary ADR 0007 decision 10 asks for
already exists in the tracklet/cell code, at the granularity this note
recommends keeping.** `TrackerTraits.cxx`'s per-policy orchestration
functions (`computeLayerTrackletsForPolicy<Tag>`, `computeLayerCellsForPolicy<Tag>`,
`findCellsNeighboursForPolicy<Tag>`, `findRoadsForPolicy<Tag>`) are each a
**single, already-shared template body** with no per-family orchestration
fork in their loop/candidate structure; `computeLayerCellsForPolicy`'s own
comment states this explicitly ("one unconditional call for both families --
no detector-ID/Tag branch here"). What remains keyed on `Tag` is:
1. a small, enumerable set of narrow leaf calls (`buildCellSeed<Tag>`,
   `attachHit<Tag>`, `cellsAreCompatible<Tag>`, `passesCellRoadPrecut<Tag>`,
   `projectSearchWindow<Tag>`, `refitHit<Tag>`) whose *bodies* already differ
   only in projection/propagation-target/rotation/material/representation
   (Section 3) -- exactly ADR decision 10's list, not a broader physics fork;
2. a handful of `if constexpr (Tag == ...)` sites inside the orchestration
   functions that construct a family-specific *input struct*
   (`TrackletProjectionState<Tag>`, `LayerScatteringInputs<Tag>`) before
   calling a leaf operation -- narrow by construction (a handful of scalar
   fields), not orchestration duplication;
3. the compile-time `Tag`/`NLayers` **template instantiation** mechanism
   itself, which produces two compiled specializations of the same source
   body. This is the one thing M5 actually removes: not duplicated logic,
   but a compile-time selector that could instead be a runtime branch on the
   transition's derived family/kind, letting one non-templated function
   (or one `NLayers`-generic function, se M6) serve both families.

The refit stage is materially different: the frozen legacy path
(`o2::its::track::fitTrack`/`refitTrack`/`refitTrackSeed`) and the native
path (`driveRefitLeg<Tag>`/`nativeRefitTrackCylinderCylinder`) are two
**independent implementations** of the same algorithm, not one shared body
with per-tag leaf calls -- see Section 4 and the harness (Section 6) for the
measured consequence of that independence.

## 3. The minimal descriptor-driven operation boundary

For each stage, the operations below are the leaf boundary: everything
upstream/downstream of them is already (or, per Section 5, should become)
one shared flow over `SurfaceDescriptor`/`SurfaceKinematicState`, dispatched
by inspecting the actual endpoint descriptor kind -- never by a persisted
tag, and never by widening these leaves' own scope.

### Tracklet projection

| Operation | Barrel | Forward | Irreducible difference |
|---|---|---|---|
| `TrackletProjectionState<Tag>` construction | `fromLayer,toLayer,meanDeltaR,targetMinR/MaxR,sourcePositionResolution,transitionMSAngle,transitionPhiCut` | `fromLayer,toLayer,fromZ,toZ,meanDeltaZ,sourceReferenceRadius,transitionMSAngle,transitionBendingAngle` | projection-target geometry (radius vs. z) |
| `projectSearchWindow<Tag, NLayers>` | LUT window in (row,col) from a cylindrical projection | LUT window in (x,y) from `MFTFwdTrackHelpers`' planar projection | projection formula |
| `TrackletSearchWindow<Tag>::acceptCandidate` | phi/tanLambda gate | x/y sigma gate | measurement-space gate shape |

Both are still legacy-layer-indexed (`fromLayer`/`toLayer` into
`TrackingParameters`/`IndexTableUtils<NLayers>`/`TimeFrame` per-layer
storage), not `SurfaceId`-indexed. This is the one piece of the tracklet
stage this note classifies as **temporary legacy branch**, not irreducible
family difference (Section 5) -- the projection *formula* is family-local,
but the *indexing scheme* it is expressed in is legacy, and a future
descriptor-driven rewrite would key it by `SurfaceId`/`SurfaceDescriptor`
directly.

### Triplet/cell construction and compatibility

| Operation | Signature (already descriptor/state-driven) | Irreducible difference |
|---|---|---|
| `buildCellSeed<Tag>` | `(SurfaceMeasurement×3, NominalSurfaceMaterial×3, bz, absCharge, pid, SurfaceKinematicState&, chi2&, Params, reason&)` | anchor formula (`barrel::buildSeed` vs `forward::buildSeed`), material-slot order, PID-aware energy-loss activation (disk only) |
| `attachHit<Tag>` | `(SurfaceKinematicState&, SurfaceMeasurement, NominalSurfaceMaterial, bz, chi2&, Params, reason&)` | rotate (barrel only -- forward has no `rotate`, alpha is always 0) + propagate target (X vs Z) + chi2 `< 0.f` rejection (barrel only) |
| `cellsAreCompatible<Tag>` | `(SurfaceKinematicState current, next, clusterIndices, bz, Params)` | `barrel::stateChi2` vs `forward::stateChi2`; disk additionally gates on raw cluster-index continuity |
| `passesCellRoadPrecut<Tag>` | `(GlobalPoint3F×3, layer indices, perLayerReferenceZ, Params)` | CylinderCylinder: compile-time no-op. DiskDisk: real geometric road cut | **legacy-only** (Section 5) |

These four already take `SurfaceMeasurement`/`SurfaceKinematicState`/
`NominalSurfaceMaterial` -- real descriptor-driven types, not legacy
`Cluster`/`TrackingFrameInfo` -- and are called from `TrackerTraits.cxx`
through one shared candidate loop per stage (Section 2). This is the
boundary this note recommends **keeping as the target shape**: a `Tag`
template parameter is a compile-time stand-in for "which of these four
leaves to call"; replacing it does not need to touch any of the four
signatures, only how the *caller* selects `Tag` (Section 5).

### Refit

| Operation | Native (`SurfaceKinematicState`) | Frozen legacy (`o2::track::TrackParCov`) |
|---|---|---|
| Per-hit step | `refitHit<Tag>` → `barrel::`/`forward::` rotate/propagate/correctForMaterial/predictedChi2/update | `fitTrack`'s inline rotate/propagateToX/correctForMaterial/getPredictedChi2Quiet/update, one shared loop body for both families (dispatch is via the caller picking `NLayers`/`layerxX0`, not a per-family branch inside `fitTrack` itself) |
| Leg orchestration | `driveRefitLeg<Tag>` (generic leg walker, Tag only selects `refitHit<Tag>`) | `refitTrack` (three-leg A/B/optional-C sequencing, `TrackFitContext<NLayers>`) |
| Whole-seed driver | `nativeRefitTrackCylinderCylinder<NLayers>` (Barrel only today) | `refitTrackSeed<NLayers>` (both families, via `DetectorTraits<NLayers>::refitSeed`'s `if constexpr (DetId == MFT)` branch to a *different* function, `refitTrackFwd`, not this chain) |

Material correction is the one place refit is **already** at the ADR-10
target: `barrel::correctForMaterial`/`forward::correctForMaterial` both
delegate to the same `material::calculateMaterialPhysics()`
(`MaterialPhysics.h`) -- the PID/absCharge-aware, direction-signed formula.
The frozen legacy `fitTrack` calls its own fixed-positive-sign,
non-PID-aware `correctForMaterial(*linRef, xOverX0, xOverX0*Radl*Rho, true)`
instead. This is not a leaf-operation gap to close by more decomposition; it
is a **named, already-documented intentional divergence** between the
native and legacy formulas (Section 5, and measured in Section 6).

## 4. Why refit is not "duplicated orchestration"

Unlike the tracklet/cell stage, `driveRefitLeg<Tag>`/
`nativeRefitTrackCylinderCylinder` and `fitTrack`/`refitTrack`/
`refitTrackSeed` are **two separately written, separately validated
implementations** of the same three-leg refit algorithm, not one shared body
with per-tag leaves. Both were verified (`testDriveRefitLeg.cxx`,
`testNativeCylinderCylinderRefitDriver.cxx`, and this note's own harness) to
reproduce the same *leg sequencing, acceptance-gate structure, and
transactionality contract* as the frozen legacy chain, but their per-hit
Kalman update arithmetic is independently implemented and does not agree
bit-for-bit even at zero material (Section 6). Collapsing these into one
body is exactly the deferred, separately-decided M5 activation step (ADR
0007 decision 11) -- this note's harness produces the evidence that decision
will need, without pre-empting it.

## 5. Classification

| Category | Items | Disposition |
|---|---|---|
| **Irreducible family-local leaf operations** (Section 3's tables; keep, key off descriptor/state family, never widen) | `barrel::`/`forward::` `rotate`/`propagate`/`predictedChi2`/`update`/`stateChi2`/`buildSeed`/`shiftReferenceToMeasurement`; `buildCellSeed<Tag>`/`attachHit<Tag>`/`cellsAreCompatible<Tag>`/`refitHit<Tag>` (the dispatch layer over the above); `TrackletProjectionState<Tag>`'s family-specific fields; `driveRefitLeg`'s per-leg direction/gate parameters | Keep as-is. These *are* ADR decision 10's narrow surface-operation boundaries already. |
| **Temporary legacy branches** (family-specific code that exists only to reproduce a frozen legacy quirk, not a genuine physics necessity) | `passesCellRoadPrecut<DiskDisk>` (hardcoded MFT road-cut constants); `LayerScatteringInputs<Tag>`/`layerMultipleScatteringAngle<Tag>` + `DiskDiskReferenceCoordinateView`/`bindLegacyMFTReferenceCoordinates()` (legacy MFT z-coordinate lookup); `clampTransitionCurvature<Tag>` (float-vs-double literal preserved bit-for-bit, not a physics choice); `TrackletProjectionState<Tag>`/`projectSearchWindow<Tag>`'s legacy-layer indexing (the *formula* is family-local, the *indexing scheme* is legacy); `SeedAnchor::Inner` (exists only for the frozen ITS reverse-anchor convention, see `SeedAnchor.h`); the entire `o2::its::track::fitTrack`/`refitTrack`/`refitTrackSeed` chain and `DetectorTraits::refitSeed`'s ITS branch | Stay exactly as they are until the specific replay-gated milestone that retires each (mostly M6, refit-formula unification is post-M5-decision). Not touched by M5a. |
| **Duplicated orchestration M5 removes** | `TrackerTraits<NLayers>::{computeLayerTrackletsForPolicy,computeLayerCellsForPolicy,findCellsNeighboursForPolicy,findRoadsForPolicy,processNeighbours}<Tag>`'s compile-time `Tag`/`NLayers` dual-instantiation (the bodies are already one shared template, per Section 2); `TransitionPolicyGrouping`/`dispatchTransitionPolicies`'s tag-keyed bucketing (pre-groups by tag instead of deriving the needed leaf at the point of use) | M5 implementation target: one non-Tag-templated (or `NLayers`-generic-only) function per stage, selecting the Section 3 leaf calls from the transition's own `SurfaceKind` at the call site. Not touched by M5a. |
| **Intentional behavior changes requiring separate approval** | Activating native refit (`nativeRefitTrackCylinderCylinder`) in place of the frozen legacy ITS refit (ADR 0007 decision 11); unifying the PID-aware vs. non-PID-aware material-correction formulas (Section 3); resolving the `clampTransitionCurvature` literal-precision quirk | Each requires its own recorded decision with A/B evidence before any production code changes. M5a's harness (Section 6) is the first such evidence for native-refit activation; it does not constitute the decision itself. |

## 6. A/B harness summary

See `testITSNativeVsLegacyRefitCharacterizationHarness.cxx` for the full
scenario matrix and `git log`/the test file's own header for provenance.
Structural invariants (success/failure classification, cluster pattern/hit
mask, timestamp, rollback-on-failure) match exactly in every scenario, as
expected: both paths consume the same seed/measurement inputs and apply the
same acceptance-gate *structure* (Section 3-4). Numeric state/covariance/chi2
are **not** asserted equal -- they are recorded as characterization, because
the two Kalman update implementations are independently written and do not
agree even at exactly zero material (first observed in
`testNativeCylinderCylinderRefitDriver.cxx`'s
`SingleHitZeroMaterialCovarianceDivergenceCharacterization`, ~1.6% single-hit
Y-Y covariance difference, ~0.2% chi2/parameter difference after a full
multi-hit leg). This harness extends that existing evidence with a wider,
JSON-reported scenario matrix (single-hit, multi-hit, holes, zero/nonzero
material, shift-reference on/off, repeat-refit-out on/off, chi2-gate
rejection, MinPt rejection) rather than superseding it.

## 7. Explicit non-actions (recap)

This note and its companion harness do not: wire native refit into any
production call site; change any cell/tracklet/refit formula; reintroduce
`TransitionPolicyTag` to any public/non-detail header; replace
`LegacyTrackerScratch<NLayers>`; alter any workflow, default, or output
contract. `ctest -L itsmft` and the existing accepted replay baselines are
therefore unaffected by this milestone; no new replay was required or run.

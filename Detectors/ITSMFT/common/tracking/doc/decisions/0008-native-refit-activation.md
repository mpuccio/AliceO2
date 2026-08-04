# ADR 0008: activate the native, descriptor-driven Propagator for cylinder and disk final refit

Status: Accepted
Date: 2026-08-03
Decision anchor: [ADR 0007](0007-generic-tracking-engine-boundary.md) decision 11 (this is
the "separate recorded decision, backed by A/B validation" that decision requires before
any descriptor-driven implementation slice may change production physics); milestone M5
(closed) and M5d of [GenericTrackingEngineMigration.md](../GenericTrackingEngineMigration.md).
Design context: [design note 0001](../design/0001-descriptor-driven-operation-boundary.md)
Section 6 (A/B characterization evidence) and Section 5 ("Intentional behavior changes
requiring separate approval").

## Context

M5a built and characterized (never activated) a native `SurfaceKinematicState` refit
driver for the Barrel family alongside the frozen legacy ITS refit
(`o2::its::track::fitTrack`/`refitTrack`/`refitTrackSeed`), and recorded that their
per-hit Kalman update arithmetic disagrees numerically even at exactly zero material
(design note 0001 Section 6). Until this decision, `DetectorTraits::refitSeed` still
called that frozen legacy chain for the barrel/ITS branch, and a second, independent
frozen legacy engine (`o2::mft::TrackFitter<TrackLTF>`, `MFTTracking/TrackFitter.h`) for
the forward/MFT branch -- two separately-implemented Kalman fitters, not the "one
tracklet/cell/road/refit flow" ADR 0007 decision 10 requires, and both still `TrackParCov`/
`TrackParCovFwd`-typed rather than `SurfaceKinematicState`-typed.

## Decision

1. Both `DetectorTraits::refitSeed` branches (ITS/barrel and MFT/forward) now call the
   same shared, descriptor-driven native refit orchestration
   (`fitTrackSeedLegs<NLayers>`, `NativeRefitDriver.h`, built on the new `Propagator`
   class, `Propagator.h`/`.cxx`) instead of either frozen legacy Kalman engine.
2. This is an intentional, approved numerical-departure commit range. Native and frozen
   output are not expected to agree, and this decision does not claim they do -- see
   design note 0001 Section 6 for the already-recorded barrel A/B evidence, and this
   milestone's own validation record (below) for the replay-level consequence.
3. `Propagator` is the new generic, `SurfaceKinematicState`-based counterpart to the
   pre-existing `o2::base::Propagator` (`DetectorsBase/Propagator.h`, left completely
   unmodified): it transports `SurfaceKinematicState` along the nominal field using
   this library's own nominal-material mechanism (`SurfaceDescriptor::material`, ADR
   0001) in place of a hit-based MatLUT lookup. It adds one new capability beyond what
   the M5a-era `barrel::`/`forward::` primitives already had: `Propagator::convertFamily`,
   a real coordinate/covariance transform between the Barrel and Forward
   representations, evaluated at the state's current reference surface (never a
   family-flag flip). `Propagator::propagateToMeasurement` composes compatibility
   checking, conversion-if-needed, rotate/propagate, material correction, the chi2
   gate, and the Kalman update into one transactional operation, selecting Barrel vs.
   Forward behavior from the target `SurfaceDescriptor::kind` alone -- never from
   `TransitionPolicyTag`, a persisted `SurfaceKindPair`, or a detector identity.
4. `Propagator::driveRefitLeg` is the shared (non-family-templated) counterpart of the
   pre-existing `detail::driveRefitLeg<Tag>`; `fitTrackSeedLegs<NLayers>` is the shared
   whole-seed three-leg (inward A / outward B / optional repeat C) driver built on it,
   serving both families with no branch of its own -- the "one tracklet/cell/road/refit
   flow" ADR 0007 decision 10 requires, now also true for the refit stage.
5. Removed from the new common tracker (`Detectors/ITSMFT/common/tracking`), proven dead
   at removal:
   - `DetectorTraits.cxx`'s call into `o2::its::track::refitTrackSeed`/`TrackFitContext`
     (the `#include "ITStracking/TrackHelpers.h"` that pulled it in);
   - `MFTFwdTrackHelpers.h`/`.cxx`'s `mftFwdPropagateToZ`/`mftFwdPredictedChi2`/
     `mftFwdStateChi2`/`mftFwdAttachCluster`/`mftFwdFitCellClusters` and the
     `o2::mft::TrackFitter<TrackLTF>`/`TrackLTFL` engine call they drove
     (`MFTTracking/TrackFitter.h`);
   - `TrackerTraits.cxx`'s per-iteration `o2::base::Propagator::Instance()`/raw
     `TrackingFrameInfo*`/`Cluster*` array staging that fed the two calls above (the
     `#include "DetectorsBase/Propagator.h"` and `#include "ITStracking/TrackHelpers.h"`
     they required).
   `DetectorTraits::refitSeed`'s public signature changed accordingly (dropped
   `tfInfos`/`unsortedClusters`/`propagator`, added `SurfaceCatalogView surfaceCatalog`)
   -- an internal common-tracking-library API change, not an output/writer/workflow
   contract change.
6. `legacy::exportBarrelTrackParCov`/`exportLegacyForwardTrackParCov`
   (`SurfaceKinematicStateLegacyAdapters.h`) are retained and still used at the very end
   of both branches: they are plain data-format adapters into `TrackITSExt`'s
   `TrackParCovF`-based storage and `TrackMFT`'s `TrackParCovFwd`-based storage
   respectively, not fitting algorithms, and every detector output type in this library
   is still exported this way. `mftLayerZ`/`mftLayerMSAngle`/`mftTrackletProject`/
   `mftTrackletSigmaXY` (`MFTFwdTrackHelpers.h`) are retained unchanged: they remain the
   live tracklet-*projection* primitives `TransitionPolicyOperations.cxx` and
   `TrackerTraits.cxx` call every iteration, a distinct CA-construction stage this
   milestone does not touch.
7. Frozen legacy ITS tracking/workflows (`Detectors/ITSMFT/ITS/tracking`,
   `o2::its::track::fitTrack`/`refitTrack`/`refitTrackSeed` themselves, and every
   `Detectors/ITSMFT/ITS/workflow*` entry point) remain completely untouched outside the
   new common tracker: they remain the external comparison reference for the later
   unified physics sign-off campaign, exactly as design note 0001 already established.
   `NativeCylinderCylinderRefitDriver.h`/`RefitLegAssembly.h` and the M5a A/B
   characterization harness (`testITSNativeVsLegacyRefitCharacterizationHarness.cxx`,
   `testNativeCylinderCylinderRefitDriver.cxx`) are likewise untouched: they remain a
   valid, still-passing characterization tool comparing the frozen legacy chain against
   the (Tag-templated, Barrel-only) M5a driver, independent of which driver production
   now calls.
8. One later, unified physics sign-off campaign -- comparing the now-fully-native common
   tracker against the frozen legacy ITS tracker and, where applicable, the retired MFT
   `TrackFitter` -- will approve or reject the resulting baseline. This commit range
   produces *candidate* results only; see the validation record below for what was
   actually measured, none of which is claimed as an accepted physics baseline.

## Non-goals

`M6` (native `SurfaceTrackingScratch`, deletion of `LegacyTrackerScratch<NLayers>`) is
untouched and unstarted; `DetectorTraversalBinding`, source-0/source-1 compatibility, and
the remaining legacy participant adapters are untouched. No workflow default, writer
contract, or output type changed. No GPU-readiness claim.

## Validation gates

- `Propagator`'s own focused tests (`testPropagator.cxx`): cylinder->cylinder and
  disk->disk propagate-to-measurement; compatible and incompatible state/measurement
  family handling; conversion followed by propagation, with full covariance/state
  validity; charge/PID/timestamp/reference preservation; MatLUT (nominal-material) and
  zero-material paths; holes; chi2-gate and minimum-pT failures; failure atomicity.
- A dedicated source-level guard test proves the new common tracker no longer includes
  `ITStracking/TrackHelpers.h`, `MFTTracking/TrackFitter.h`, or constructs
  `o2::its::track::TrackFitContext`/`o2::mft::TrackFitter` anywhere under
  `Detectors/ITSMFT/common/tracking`.
- `ctest --test-dir <build> -L itsmft --output-on-failure -j1`: **88/88 passed** (0
  failed), including the extended M5a harness (fitTrackSeedLegs reproduces
  nativeRefitTrackCylinderCylinder bit-for-bit for every scenario in its table) and the
  updated `testMFTNormalizedRefit.cxx`/`testDetectorTraitsTrackRepresentation.cxx`.
- `git diff --check <base> HEAD`: clean. `git clang-format --diff <base> HEAD`: clean.
- The durable fixture's 43 checksums (`sha256`) verified identical before and after
  every replay below -- this milestone touches no fixture, only tracking code.

### Candidate replay record (2026-08-03, `daily-20260717-0700-local1`, fixture
`pp-20ev-run303000-seed20260716-daily20260717`, `TIMESTAMP=1784207296000`,
`RUNNUMBER=303000`)

Not an accepted physics baseline -- candidate results only, per this ADR's own decision
2 and 8. Old-baseline column is the last M5-era accepted number (design note 0001 /
`AgentCoordination.md` Gate status), recorded for context, not as a target to match.

| Leg | Old (accepted, pre-M5d) | New (candidate, M5d) | Wall time | Peak RSS |
|---|---|---|---|---|
| ITS common-CA (`o2-its-ca-tracker-workflow`) | 203 tracks, hash `ee7f7c794d60f2362fd2564258b7887e` | 160 tracks, hash `4ee3e10ae920736e6b29ef35a5af4fe1`; MC eff. 78.9% (112/142), fake rate 8.75% (14/160), 0 clones | 19.81 s | 2.65 GB |
| MFT common-CA (`o2-mft-ca-tracker-workflow`) | 70 tracks, hash `24737e73b7146bf3bd35a90a2517c527` | 65 tracks, hash `f67de427ccf53a33354055d4fd8a1307`; MC eff. 50.5% (55/109), fake rate 1.5% (1/65), 0 clones | 7.90 s | 1.29 GB |
| Combined (`o2-itsmft-combined-ca-tracker-workflow`) | n/a (no combined-workflow M5-era pin) | ITS leg: 160 tracks, hash `4ee3e10ae920736e6b29ef35a5af4fe1` (bit-identical to the standalone ITS leg); MFT leg: 65 tracks, hash `f67de427ccf53a33354055d4fd8a1307` (bit-identical to the standalone MFT leg) | 6.69 s | 3.42 GB |

Structural invariants held in every leg: success/failure classification, cluster
pattern, timestamp, reference validity, publication ordering, and failure-path rollback
all passed their respective focused tests (`ctest -L itsmft`, above); the combined
workflow's per-detector output is bit-identical to each standalone single-detector
replay, proving no cross-detector topology was traversed -- the same structural proof
Gate 4's combined-coordinator acceptance evidence already established, now reproduced
with the native refit active. Track-count and efficiency deltas versus the old baseline
are the expected, disclosed consequence of activating two independently-implemented
native Kalman updates (barrel and forward) in place of two independent frozen legacy
engines; MFT's mean chi2 dropped from the legacy fitter's own scale to 0.20 (median
0.083), consistent with a well-converged fit, not a degenerate one. No numerical
tolerance was invented to declare native and frozen values equal, per this milestone's
own instruction.

## Addendum (2026-08-04): covariance-validity correction

A dedicated read-only investigation (never committed; disposable diagnostic worktrees
only) traced the ITS/MFT track-count drop recorded above to a proven root cause, not
"the Kalman implementations differ": `barrel::`/`forward::` `propagate()`/`rotate()`/
`update()` never sanitized their covariance output, unlike the retained
`o2::track::TrackParametrizationWithError<value_T>`, which calls `checkCovariance()`
(`abs()` every diagonal; clamp+rescale any over range) unconditionally after every one
of those three calls. A large single-step `propagate()` can leave the covariance matrix
non-positive-semi-definite (an off-diagonal correlation exceeding the Cauchy-Schwarz
bound) while every individual diagonal still looks valid; the very next `update()`'s
naive Kalman covariance subtraction then reveals this, mathematically and unavoidably
(reproduced with a Joseph-form, full-double-precision recomputation of the exact same
real captured inputs), as a negative diagonal -- which `FamilyMaterialOperations::
preflightValidate` correctly, but too late, rejects as `InvalidCovariance`. Confirmed
precision-independent (float32 and float64 hand-rederivations of the real captured
propagate step agree to 7 significant digits) and formula-correct (the hand-rederivation
matches real production output exactly): this is a missing operation, not a wrong one.

**Correction**: `sanitizeCovariance(state, maxDiagonal)` (`SurfaceKinematicState.h`) is
one detector-neutral, GPU-portable primitive implementing exactly `checkCovariance()`'s
policy (reusing, not duplicating, the abs/clamp/rescale logic `limitBarrelCovariance`
already had), invoked unconditionally after every successful `barrel::rotate`/
`propagate` (both overloads)/`update` and `forward::propagate<Model>` (both overloads,
all four models)/`update`. `FamilyMaterialOperations::preflightValidate` itself is
unchanged and unweakened: it still rejects a deliberately malformed *externally
supplied* state (`testCovarianceSanitization.cxx`), it just no longer sees one from the
Propagator's own normal operation. Conversion, propagation formulas/Jacobians, the
material model, chi2/min-pT gates, rollback, and acceptance/dedup are untouched.

### Validation gates (re-run after the correction)

- `ctest --test-dir <build> -L itsmft --output-on-failure -j1`: **89/89 passed** (0
  failed) -- the 88 gates above plus the 18 new focused
  `testCovarianceSanitization.cxx` cases (the sanitizer in isolation; both ITS and MFT
  real-data reproducers now sanitize instead of rejecting; the preceding large-step
  propagate preserves the invariant; every one of the eight call sites independently
  sanitizes; a deliberately malformed external state is still rejected; failure remains
  transactional). `testNoLegacyFittingDependency.cxx` still passes unmodified.
- `git diff --check <base> HEAD`: clean. `git clang-format --diff <base> HEAD`: clean.
- The durable fixture's 43 checksums verified identical before and after every replay
  below, including the disposable diagnostic-worktree replays used to confirm the
  reconciliation (this correction touches no fixture, only tracking code).

### Candidate replay record (2026-08-04, `daily-20260717-0700-local1`, same fixture)

Still not an accepted physics baseline, for the same reason as the original candidate
replay record above. Wall time was measured with a different methodology this session
(`/usr/bin/time -l` wrapping the whole reader|tracker pipe) than the original entries
above, so wall-time is not directly comparable across the two rows; peak RSS (dominated
by the fixed-size fixture/geometry/cluster load, not by this correction) lands within
noise of the original figures.

| Leg | M5d before this fix | M5d after this fix | Wall time | Peak RSS |
|---|---|---|---|---|
| ITS common-CA | 160 tracks, hash `4ee3e10ae920736e6b29ef35a5af4fe1`; eff. 78.9%, fake 8.75% | **208 tracks**, hash `fd155fde67409d6e330f69e179041a1e`; MC eff. 94.4% (134/142), fake rate 8.65% (18/208), 0 clones | 3.73 s | 2.65 GB |
| MFT common-CA | 65 tracks, hash `f67de427ccf53a33354055d4fd8a1307`; eff. 50.5%, fake 1.5% | **67 tracks**, hash `f3f988d25bdb6b423d69ffff01fb7dcc`; MC eff. 52.3% (57/109), fake rate 1.49% (1/67), 0 clones | 2.68 s | 1.29 GB |
| Combined | ITS leg 160 / MFT leg 65 (both bit-identical to standalone) | ITS leg 208 (hash `fd155fde67409d6e330f69e179041a1e`) / MFT leg 67 (hash `f3f988d25bdb6b423d69ffff01fb7dcc`), each bit-identical to its standalone single-detector replay | 4.00 s | 3.42 GB |

Legacy reference (unchanged, pre-M5d): ITS 203 tracks; MFT 70 tracks (same table,
top of this document).

**Candidate-level reconciliation** (disposable diagnostic instrumentation of the real
refit driver and the real, shared `TrackerTraits::acceptTracks`; never committed):
`InvalidCovariance` is *eliminated*, not merely reduced, in both detectors. ITS: 0/117
remaining refit failures are `MaterialFailure` (was the dominant failure category before
this fix, ~211/301); every remaining failure is an ordinary `PredictedChi2Failure` (116)
or `PropagationFailure` (1). MFT: 1/4 remaining refit failures is `MaterialFailure`, but
its `MaterialFailureReason` is `MomentumBelowMinimum`, not `InvalidCovariance` -- a
genuine, unrelated physics rejection; the other 3 are `PredictedChi2Failure` (1) and
`LegAcceptanceFailure` (2). MFT accounting verified exact: ground-truth
`TrackerTraits::acceptTracks` instrumentation reports 67 accepted, bit-identical to the
production track-writer's own output. Every remaining ITS/MFT delta versus the legacy
reference (208 vs. 203; 67 vs. 70) is an ordinary chi2/whole-leg-acceptance/material-
physics gate outcome or a dedup-competition effect between two independently-implemented
Kalman engines -- candidate physics, consistent with this ADR's own decision 2 and 8,
not evidence of a remaining defect.

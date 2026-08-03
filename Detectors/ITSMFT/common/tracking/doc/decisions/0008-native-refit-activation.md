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
  after every replay (this milestone touches no fixture, only tracking code).

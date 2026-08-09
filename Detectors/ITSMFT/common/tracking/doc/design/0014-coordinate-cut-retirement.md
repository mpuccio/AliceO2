# Coordinate-cut retirement

- Status: implemented and structurally validated; MFT physics sign-off pending
- Date: 2026-08-09
- Predecessor: [flat pair-list graph and one combined plan](0013-flat-pair-list-combined-plan.md)
- Validation: [coordinate-cut retirement](../validation/coordinate-cut-retirement.md)
- Pinned O2 package: `daily-20260717-0700-local1`

## Decision

Remove `TrackletMinAbsX` and `CellRoadRCut` from common tracking configuration,
kernel parameters, orchestration, and refit. Do not replace either with a
renamed disk parameter, a detector dispatch, or a fixed hidden cut.

The existing descriptor-selected cylinder/forward operation boundary remains.
Surface kind continues to choose coordinate mathematics inside tracklet, cell,
propagation, and refit leaves; it does not choose tracking policy.

## `TrackletMinAbsX`

The cut rejects a completed forward fit and every contributing measurement
when the absolute global X coordinate is below a configurable threshold. It
was also applied during tracklet formation in the first common-MFT port, but
only the final-refit application remains.

This is not a coordinate-neutral physical selection: rotating the global
transverse axes changes acceptance while leaving the detector, hit radii, and
track unchanged. It is not a singularity guard either. The forward seed and
propagation leaves guard their actual denominators and degenerate geometry
(layer-Z separation, transverse separation, slope, and radius) independently;
`x = 0` is valid when those invariants hold. The cut is therefore obsolete
MFT compatibility behavior and is deleted.

## `CellRoadRCut`

Legacy MFT `Tracker::findTracksCA()` used `ROADclsRCut` while constructing a
coarse road: two distant endpoint clusters defined a seed line and an
intermediate cluster was collected when its transverse distance to that line
was below the road radius, optionally scaled as a cone. That was an
implementation-specific search-window optimization, not a fitted-track
validity condition.

The common CA does not construct those roads. Its port applies three symmetric
line-extrapolation inequalities to every local three-hit cell after tracklet
search and cell-slope compatibility. This differs from both the call site and
the mathematics of the legacy optimization; retaining it as a generic cell
selection would preserve an accidental port constraint, not an algorithmic
invariant. The common cell leaves already perform coordinate-appropriate seed
construction, propagation, measurement updates, and the configured predicted
chi-square gate. The transplanted road pre-cut and its no-op cylinder partner
are therefore deleted.

Removing this cut can admit disk candidates that the current common path
rejects. Such changes are an expected characterization result, not an MFT
compatibility failure. ITS output remains the behavior-preservation gate.

## Preserved coordinate guards

The migration does not weaken genuine leaf validity checks. In particular,
disk tracklet projection retains its layer/vertex-Z denominator guards, and
forward seed construction retains its ordered-layer, nonzero-Z-separation,
and nonzero-transverse-separation guards. These stay inside the existing
forward operation leaves and are covered directly at those boundaries.

## Later-deletion inventory

The removal leaves `MFTFwdTrackHelpers.{h,cxx}` with no MFT selection policy;
it now performs normalized-measurement identity validation and then forwards
to the generic native refit driver. A later application-boundary slice should
move the identity check to normalized input validation, call the native refit
driver directly from the descriptor-selected refit leaf, and delete the
MFT-named wrapper. That slice may also narrow `DetectorRefitSupport.h` once
publication ownership moves to the detector workflows. Neither ownership
migration is part of this cut-retirement slice.

Legacy production MFT `MFTTrackingParam::ROADclsRCut` and its tracker use stay
unchanged because legacy workflow removal is explicitly outside this slice.

## Acceptance boundary

Focused tests must guard the absence of both retired names from common
production sources and configuration projection, reject source/detector tuning
dispatch, show cylinder and disk using the same stage orchestration, and
exercise the retained forward validity guards at leaf boundaries. The pinned
campaign must include the complete serial `itsmft` suite, strict preflight,
fixture integrity before and after replay, strict standalone/combined ITS
preservation, and standalone/combined MFT characterization with ordered
references, labels, timestamps where persisted, and efficiency/fake/clone
effects.

# Descriptor-driven mixed-family refit

## Decision

`detail::refitSurfaceSeed()` is the common tracker refit boundary. It validates
each attached seed position against its ordered surface, the event-owned
per-surface source mapping, normalized external index, and normalized
measurement before invoking `fitTrackSeedLegs()`.

The refit driver already advances through ordered measurements with its
private `detail::driveRefitLeg()` and `Propagator::propagateToMeasurement()`.
A
target descriptor selects the coordinate family at each measurement; a seed
is therefore not classified from its first hit and may contain Cylinder and
Disk measurements from distinct valid sources.

## Output boundary

A successful mixed refit is generic `TimeFrame` state. The existing ITS and
MFT output adapters remain detector-specific: a `CommonTrack` whose references
span their requested detector/source is rejected as `MixedDetector` before an
output request. This slice defines no mixed DPL product or routing rule.

## Conversion invariant

Family conversion preserves the fitted state's propagation anchor in the
paired `SurfaceLinearizationReference`. The linearization parameters can remain
nominal after a Kalman update, but their reference coordinate must agree with
the fitted state before the next propagation. This is required for an ordered
Cylinder-to-Disk-to-Cylinder leg.

## Scope

This removes first-hit detector/family dispatch from common refit and moves
identity checks to the generic boundary. It does not add a mixed typed output,
change workflow defaults, or alter disconnected ITS/MFT graph topology.

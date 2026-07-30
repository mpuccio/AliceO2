# ADR 0002: `CommonTrack` is the detector-neutral common-CA result type

Status: Accepted
Date: 2026-07-30

## Context

Architecture.md Sec 12 already specifies the target shape: "The TimeFrame
stores a generic internal result" (`InternalTrack` in that document's
sketch), and that `TrackITS`/`TrackMFT` are adapter outputs built from it,
never core storage types. Before this decision, no such type existed in the
common core: the only production track container is the legacy, per-`NLayers`
`CATrackType<NLayers>` (`TrackerTraits.cxx`, `DetectorTraits<NLayers>::TrackType`),
selected by which `TimeFrame<NLayers>` instantiation is running -- i.e. by
detector identity -- not by which surfaces/kinds a track actually touched.
Gate 4 (a single `TimeFrame` spanning ITS and MFT surfaces) cannot be built on
top of a result type whose very existence is parameterized by `NLayers`.

This decision establishes the result type and its `TimeFrame` storage only.
It does not populate it from CA seeds, does not touch `TrackerTraits`/
`CATracker` candidate/road/refit logic, and does not change any existing
output publication -- see "Non-goals" below.

## Decision

`CommonTrack` (`ITSMFTTracking/CommonTrack.h`) owns exactly:

```cpp
struct CommonTrack {
  SurfaceKinematicState innerState{};
  SurfaceKinematicState outerState{};
  float chi2{0.f};
  o2::its::TimeEstBC timestamp{};
  SurfaceMask hitSurfaces{};
  uint32_t firstClusterRef{0};
  uint32_t clusterRefEnd{0};
};
```

- `innerState`/`outerState` reuse `SurfaceKinematicState` directly (the same
  Barrel/Forward-tagged five-parameter/covariance state already used by
  `buildCellSeed`/`attachHit`/`refitHit`, `ITSMFTTracking/
  SurfaceKinematicState.h`). This is deliberately not a second state
  representation.
- `hitSurfaces` is the global 32-bit `SurfaceMask` (`ITSMFTTracking/
  SurfaceMask.h`), never the legacy 16-bit per-`NLayers` `LayerMask`.
- Cluster membership is not stored inline. `[firstClusterRef, clusterRefEnd)`
  is a half-open range of *positions* into a new, `TimeFrame`-owned flat
  array, `TimeFrame::getTrackClusterIndices()` (element type
  `SurfaceMeasurementIndex`, `ITSMFTTracking/SurfaceMeasurementIndex.h`).
  Each element of that array is, in turn, a canonical position into the
  single flattened `SurfaceMeasurement` array owned by the TimeFrame's
  normalized frame (`MultiSourceFrame`) -- never a detector-local layer index
  or a raw external cluster index. References are stored in traversal order,
  inner to outer.
- A range is valid iff `0 <= firstClusterRef <= clusterRefEnd <=
  trackClusterIndices.size()` (`isValidTrackRange()`, `CommonTrack.h`); the
  lower bound is automatic since both fields are unsigned.
- For a valid, completed track, `hitSurfaces` is the union of the `SurfaceId`
  of every measurement its range references. This is a consumer-side
  invariant (exercised by `testCommonTrack.cxx`), not something `CommonTrack`
  enforces by construction -- no code populates it from real seeds in this
  slice.
- `CommonTrack` carries no `DetID`, no `NLayers`, and no reference to
  `TrackITS`/`TrackITSExt`/`MFTCATrack`/`GeometryTGeo`/any workflow header.
  `TimeEstBC` (`DataFormatsITS/TimeEstBC.h`) is reused as-is: despite its
  header path, it is a common timing primitive already used throughout this
  library (`Cell.h`, `Tracklet.h`, `TransitionPolicyOperations.h`), not a
  detector output type.

`TimeFrame<NLayers>` gains two new members, `mCommonTracks` (`bounded_vector<
CommonTrack>`) and `mTrackClusterIndices` (`bounded_vector<
SurfaceMeasurementIndex>`), with `getCommonTracks()`/`getTrackClusterIndices()`
accessors. Neither element type depends on `NLayers`; this is a **temporary
bridge** holding a detector-neutral type inside the still-`NLayers`-templated
`TimeFrame`, exactly as the still-templated `TimeFrame` already holds other
non-`NLayers`-dependent members (e.g. `mPrimaryVertices`). Both containers are
event data: `TimeFrame::wipe()` clears them unconditionally, alongside every
other per-event CA artefact (`mTracks`, `mTracklets`, ...), and
`setMemoryPool()` initializes them like every other host-only bounded vector.
A `CommonTrack`'s range, and every `SurfaceMeasurementIndex` it reaches
through `trackClusterIndices`, are only meaningful together with the
`TimeFrame`'s normalized frame that was current when the track was built;
none of the three (`CommonTrack`, `trackClusterIndices`, the normalized
frame) is individually meaningful once any of the other two has been wiped or
reloaded.

`MultiSourceFrame`/`MultiSourceFrameView` gain a bounds-checked
`getMeasurement(SurfaceMeasurementIndex)` lookup -- the minimal support needed
to validate and dereference a stored `SurfaceMeasurementIndex` against the
flattened measurement store. This mirrors the existing per-surface
`getSurfaceMeasurements(SurfaceId)` accessor but indexes the flat array
directly rather than one surface's span within it.

## Adapter boundary

`CommonTrack` is core storage; `TrackITS`/`TrackITSExt`/`TrackMFT`/a future
ALICE-3 type are adapter outputs built from it by code outside the common
core (Architecture.md Sec 12), never the reverse. The adapter direction is:

```text
CommonTrack (+ TimeFrame::getTrackClusterIndices() + MultiSourceFrame)
        -> detector output adapter (not yet written)
        -> TrackITSExt / TrackMFT / ...
```

No such adapter exists yet: `TrackerTraits`/`CATracker` still populate the
legacy `CATrackType<NLayers>` exclusively, and every existing workflow output
is unchanged. Writing the adapter, and deciding exactly how/when a
`CommonTrack` is constructed from a completed CA road, is explicitly out of
scope for this decision (see "Non-goals").

## Non-goals

Populating `CommonTrack` from `TrackerTraits`/`CATracker` seeds or completed
roads; any change to `TrackerTraits`, `CATracker` candidate/road/refit logic,
`TrackITS`, `MFTCATrack`, legacy workflows, or current output publication;
real `ITSSurfaceSpec`/`MFTSurfaceSpec` constants; any change to topology
authority, hole enforcement, `TransitionPolicyTag`, or detector layout logic;
an ALICE-3 adapter. The existing detector-specific tracks remain the only
production outputs for now; this decision establishes only the common result
storage they will later be adapted from.

## Validation gates

- POD layout: `CommonTrack` is standard-layout and trivially-copyable
  (`static_assert` in `CommonTrack.h`).
- `isValidTrackRange()`: empty/default, single-, multi- and hole-containing
  ranges accepted; out-of-range and reversed ranges rejected
  (`testCommonTrack.cxx`).
- `hitSurfaces` equals the union of referenced measurements' `SurfaceId`,
  checked both for a single ITS source and for combined ITS+MFT sources
  sharing one `MultiSourceFrame` (`testCommonTrack.cxx`), including explicit
  cross-source/cross-detector identity checks on the resolved measurements.
- `TimeFrame::wipe()` clears `mCommonTracks`/`mTrackClusterIndices` together
  and both accept independent reload afterward (`testCommonTrack.cxx`).
- `CommonTrack.h`'s own include list contains no `DetID.h`, `TrackITS.h`/
  `TrackITSExt.h`, `MFTCATrack.h`, `GeometryTGeo.h`, or workflow header.

## Alternatives rejected

- **A second five-parameter/covariance state type for `innerState`/
  `outerState`**: rejected outright by the task defining this decision;
  `SurfaceKinematicState` is already the shared Barrel/Forward state and
  introducing a parallel representation would duplicate exactly the kind of
  machinery Architecture.md's Stage-A/Stage-B track-state strategy already
  consolidated.
- **Inline cluster storage on `CommonTrack` (a fixed-size array of
  references)**: rejected; a fixed-size in-struct array would reintroduce an
  `NLayers`-shaped bound on a type whose entire purpose is to have none. The
  flat `TimeFrame`-owned array plus a range is the same range-into-flat-array
  shape `SparseTrackingTopologyView`/`MultiSourceFrameView` already use
  elsewhere in this codebase for the same reason.
- **Storing raw `SurfaceMeasurement*`/cluster indices directly on
  `CommonTrack`**: rejected; `SurfaceMeasurementIndex` is the existing,
  already-defined canonical position type for the flattened measurement
  store (`ITSMFTTracking/SurfaceMeasurementIndex.h`), and reusing it keeps
  `CommonTrack`/`trackClusterIndices` trivially-copyable, GPU-representable
  PODs with no embedded pointers.

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

An initial revision of this decision briefly considered a single, flattened
global `SurfaceMeasurement` array (`MultiSourceFrame` owning one shared
array plus per-surface offset ranges, with `SurfaceMeasurementIndex`
reinterpreted as a position in that shared array). That reinterpretation was
rejected before integration: `SurfaceMeasurementIndex` remains what it always
was -- a strong index local to one surface's own measurement array -- and
`MultiSourceFrame` owns measurements per surface. The corrected model below
is the only one ever integrated.

## Decision

### `CommonTrack` and `TrackClusterReference`

`CommonTrack` (`ITSMFTTracking/CommonTrack.h`) owns exactly:

```cpp
struct CommonTrack {
  SurfaceKinematicState innerState{};
  SurfaceKinematicState outerState{};
  float chi2{0.f};
  CommonTrackTimestamp timestamp{};
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
- `timestamp` is `CommonTrackTimestamp` (`ITSMFTTracking/SurfaceTiming.h`): a
  new, common, `TFBC`-based half-open BC interval (`{TFBC begin; TFBC end;}`),
  standard-layout and trivially-copyable, with an `isValid()`/`isCompatible()`
  contract mirroring `ROFIntervalBC`'s own half-open-interval semantics
  exactly (adjacent/touching intervals do not intersect; an invalid interval
  is never compatible with anything). This replaces `o2::its::TimeEstBC`:
  `TimeEstBC`'s own `TimeStampWithError`/`TimeStamp` base hierarchy declares
  non-static data members at more than one level, which makes `TimeEstBC` --
  and therefore any aggregate containing it -- not standard-layout, despite
  being trivially-copyable. `CommonTrack` must be both: trivial-copyability
  alone only proves its bytes may be copied without invoking non-trivial
  special member functions, it says nothing about a single, well-defined,
  `offsetof`-usable member order consistent across host/device compilation.
  Both properties are asserted together on `CommonTrack`, `TrackClusterReference`
  and `CommonTrackTimestamp`, matching every other device-facing type in this
  library (`SurfaceMeasurement`, `StaticSurfaceDescriptor`,
  `DetectorLayoutView`, ...).
- Cluster membership is not stored inline. `[firstClusterRef, clusterRefEnd)`
  is a half-open range of *positions* into a new, `TimeFrame`-owned flat
  array, `TimeFrame::getTrackClusterIndices()` (element type
  `TrackClusterReference`):

  ```cpp
  struct TrackClusterReference {
    SurfaceId surface{};
    SurfaceMeasurementIndex index{};
  };
  ```

  Measurements are owned **per surface** (`MultiSourceFrame`, see below), so
  a bare `SurfaceMeasurementIndex` is never a complete identity on its own --
  it is only ever meaningful paired with the `SurfaceId` it is local to. A
  reference resolves through exactly one call:
  `normalizedFrame.getMeasurement(reference.surface, reference.index)`.
  References are stored in traversal order, inner to outer. Never a
  detector-local raw cluster index or external cluster ID: those remain
  `ClusterRef`'s job (source-qualified, `ITSMFTTracking/SurfaceMeasurement.h`),
  a different identity axis entirely.
- A range is valid iff `0 <= firstClusterRef <= clusterRefEnd <=
  trackClusterIndices.size()` (`isValidTrackRange()`, `CommonTrack.h`); the
  lower bound is automatic since both fields are unsigned.
- For a valid, completed track: every `TrackClusterReference` in
  `[firstClusterRef, clusterRefEnd)` resolves to an existing measurement
  whose own `SurfaceMeasurement::surface` equals that reference's `surface`
  field, and `hitSurfaces` is the union of those surfaces. This is a
  consumer-side invariant (exercised by `testCommonTrack.cxx`), not something
  `CommonTrack` enforces by construction -- no code populates it from real
  seeds in this slice.
- `CommonTrack`/`TrackClusterReference` carry no `DetID`, no `NLayers`, and no
  reference to typed detector output/`GeometryTGeo`/any
  workflow header.

`TimeFrame<NLayers>` gains two new members, `mCommonTracks` (`bounded_vector<
CommonTrack>`) and `mTrackClusterIndices` (`bounded_vector<
TrackClusterReference>`), with `getCommonTracks()`/`getTrackClusterIndices()`
accessors. Neither element type depends on `NLayers`; this is a **temporary
bridge** holding detector-neutral types inside the still-`NLayers`-templated
`TimeFrame`, exactly as the still-templated `TimeFrame` already holds other
non-`NLayers`-dependent members (e.g. `mPrimaryVertices`).

### `MultiSourceFrame`: per-surface measurement ownership

Measurements are owned **per `SurfaceId`** -- one `std::vector<
SurfaceMeasurement>` per surface (`MultiSourceFrame::mPerSurfaceMeasurements`)
-- never as a single flattened array with per-surface offset ranges.
`SurfaceMeasurementIndex` (`ITSMFTTracking/SurfaceMeasurementIndex.h`) is
unchanged and remains exactly what it always was: a strong `uint32_t` index
local to one surface's own measurement array, never a position in some
larger shared array.

`MultiSourceFrameView`'s device-facing surface accessor is a POD array of
per-surface pointer/count spans (`SurfaceMeasurementSpan{const
SurfaceMeasurement* data; uint32_t count;}`), indexed directly by `SurfaceId`
-- not a flattened array plus offset ranges. `getSurfaceMeasurements(SurfaceId)`/
`getSurfaceMeasurementCount(SurfaceId)` keep their existing per-surface
contract, unchanged, on both the owner and the view. A bounds-checked
`getMeasurement(SurfaceId, SurfaceMeasurementIndex)` (on both the owner and
the view) is the minimal support `TrackClusterReference` needs to validate
and dereference a stored reference; there is no global-index accessor.

**Cached-span safety.** `MultiSourceFrame::mSurfaceSpans` (the owner's cache
backing `getView()`) holds raw pointers into `mPerSurfaceMeasurements`. A
member-wise (default) copy would copy those pointers as-is, leaving the
copy's `mSurfaceSpans` pointing into the *original* object's storage instead
of its own newly-duplicated one -- silently wrong immediately after the copy,
not merely unsafe. `MultiSourceFrame`'s copy constructor and copy assignment
are therefore deleted rather than implemented incorrectly by default; no
caller needs copy (`TimeFrame` owns exactly one `MultiSourceFrame` per
normalized load and only ever moves it). Move is implemented via an explicit
`swap()` rather than the member-wise move `std::vector` would otherwise
generate implicitly: both are semantically equivalent here (moving a
vector-of-vectors never touches the innermost heap buffers `mSurfaceSpans`'
cached pointers reference), but `swap()` makes the `noexcept` guarantee
provable from the operator's own definition (`std::vector::swap` is
unconditionally noexcept: a pointer/size/capacity exchange only) rather than
from a conditional trait computed over every member's move operation, and it
leaves the moved-from side holding a fully self-consistent state (the
moved-into side's *old* content, spans included) rather than an
unspecified-but-valid moved-from state. This is what guarantees a successful
`TimeFrame::loadNormalizedSource()` commit (`mNormalizedFrame =
std::move(staged)`) can never leave either side's cached spans pointing into
the other side's storage.

### Reload/wipe lifecycle

`TimeFrame::loadNormalizedSource()`'s successful commit also clears
`mCommonTracks`/`mTrackClusterIndices`, in that same commit, since a
`CommonTrack`/`TrackClusterReference` set built against the just-replaced
normalized frame is meaningless once that frame is gone -- a stale reference
into a frame that no longer exists is worse than an empty one. A **failing**
call leaves the normalized frame, `mCommonTracks` and `mTrackClusterIndices`
all completely unchanged, matching that call's existing transactional
contract for every other `TimeFrame` member (preflight, decoding, and the
complete legacy backfill are staged before anything is committed).
`TimeFrame::wipe()` continues to clear all three (normalized frame,
`mCommonTracks`, `mTrackClusterIndices`) unconditionally, alongside every
other per-event CA artefact, and `setMemoryPool()` initializes the two new
containers like every other host-only bounded vector.

A `CommonTrack`'s range, and every `TrackClusterReference` it reaches through
`trackClusterIndices`, are only meaningful together with the `TimeFrame`'s
normalized frame that was current when the track was built; none of the
three is individually meaningful once any of the other two has been wiped or
reloaded.

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
typed detector output, legacy workflows, or current output publication;
real `ITSSurfaceSpec`/`MFTSurfaceSpec` constants; any change to topology
authority, hole enforcement, `TransitionPolicyTag`, or detector layout logic;
an ALICE-3 adapter. The existing detector-specific tracks remain the only
production outputs for now; this decision establishes only the common result
storage they will later be adapted from.

## Validation gates

- POD layout: `CommonTrack`, `TrackClusterReference` and
  `CommonTrackTimestamp` are all standard-layout and trivially-copyable
  (`static_assert`s in `CommonTrack.h`/`SurfaceTiming.h`).
- `MultiSourceFrame` is move-only, provably nothrow-movable
  (`static_assert(std::is_nothrow_move_constructible_v<MultiSourceFrame>)`/
  `..._assignable_v...`) and not copyable
  (`static_assert(!std::is_copy_constructible_v<MultiSourceFrame>)`/
  `..._assignable_v...`), `MultiSourceFrame.h`.
- `isValidTrackRange()`: empty/default, single-, multi- and hole-containing
  ranges accepted; out-of-range and reversed ranges rejected
  (`testCommonTrack.cxx`).
- Per-surface measurement storage and surface-local index bounds: an index
  valid for one surface but out of range for another (or for an empty
  surface) is rejected on both the owner and the view
  (`testCommonTrack.cxx`, `testMultiSourceLoading.cxx`).
- Cross-surface and cross-source `TrackClusterReference` resolution, and
  `hitSurfaces`-equals-referenced-surfaces (including per-reference
  `SurfaceMeasurement::surface` consistency), checked both for a single ITS
  source and for combined ITS+MFT sources sharing one `MultiSourceFrame`
  (`testCommonTrack.cxx`).
- `MultiSourceFrameView` resolves measurements on multiple surfaces
  correctly immediately after a real, successful load, cross-checked
  against the owner's own per-surface accessor
  (`testMultiSourceLoading.cxx`).
- A successful `TimeFrame::loadNormalizedSource()` reload clears
  `mCommonTracks`/`mTrackClusterIndices`; a failed one leaves both (and the
  normalized frame) completely unchanged; `TimeFrame::wipe()` clears all
  three together and each accepts independent reload afterward
  (`testCommonTrack.cxx`).
- `CommonTrack.h`'s own include list contains no `DetID.h`, `TrackITS.h`/
  `TrackITSExt.h`, a typed MFT output header, `GeometryTGeo.h`, or workflow
  header.

## Alternatives rejected

- **A single flattened global `SurfaceMeasurement` array with per-surface
  offset ranges, `SurfaceMeasurementIndex` reinterpreted as a position in
  that shared array**: rejected before integration (see "Context"). Measurements
  are owned per surface; `SurfaceMeasurementIndex` stays surface-local.
- **A second five-parameter/covariance state type for `innerState`/
  `outerState`**: rejected outright; `SurfaceKinematicState` is already the
  shared Barrel/Forward state and introducing a parallel representation
  would duplicate exactly the kind of machinery Architecture.md's
  Stage-A/Stage-B track-state strategy already consolidated.
- **Inline cluster storage on `CommonTrack` (a fixed-size array of
  references)**: rejected; a fixed-size in-struct array would reintroduce an
  `NLayers`-shaped bound on a type whose entire purpose is to have none. The
  flat `TimeFrame`-owned array plus a range is the same range-into-flat-array
  shape `SparseTrackingTopologyView` already uses elsewhere in this codebase
  for the same reason.
- **Storing a bare `SurfaceMeasurementIndex` (no `SurfaceId`) as track
  membership identity**: rejected; once measurements are owned per surface, a
  surface-local index alone is not a complete identity. `TrackClusterReference`
  pairs it with the `SurfaceId` it is local to.
- **`o2::its::TimeEstBC` for `CommonTrack::timestamp`**: rejected; not
  standard-layout (see "Decision"), and an ITS-namespaced type despite being
  reused elsewhere in this library. `CommonTrackTimestamp` is the common,
  `TFBC`-based replacement.
- **Relying on `MultiSourceFrame`'s implicitly-generated copy/move
  operations**: rejected; the implicit copy would silently misbind
  `mSurfaceSpans` to the wrong object's storage (see "Decision"). Copy is
  deleted; move is implemented via an explicit, provably-`noexcept` `swap()`.

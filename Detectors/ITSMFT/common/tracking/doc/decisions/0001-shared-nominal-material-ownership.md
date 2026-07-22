# ADR 0001: Nominal per-surface material is a field on SurfaceDescriptor

Status: Accepted
Date: 2026-07-22

## Context

Nominal per-surface material (a normal-incidence `xOverX0` and areal density
per surface) needs a runtime home so it can eventually be queried the same
way surface geometry already is.

The first accepted design (superseded by this one) gave material a fully
standalone ownership unit: a parallel `std::vector<NominalSurfaceMaterialBudget>`
held alongside `DetectorLayoutSet`'s surface catalog, populated from
identity-bearing `NominalSurfaceMaterialEntry` construction input that had to
be positionally validated against the catalog, exposed through
`DetectorLayoutView::nominalMaterial`/`getNominalMaterial()`, and kept current
through its own `MaterialCatalogEpoch` field on
`DetectorLayoutConfigurationKey`, independent of `DetectorGeometryEpoch`. That
design also added `DetectorLayoutViewStatus`/`isValid()` to `DetectorLayoutView`
so a legitimately empty-but-valid view could be distinguished from the
invalid/no-such-iteration sentinel -- needed only because a parallel array's
size-matching couldn't otherwise be trusted to imply validity.

Building it revealed the design was solving a problem it created for itself:
material is not a separate fact about the detector that happens to share an
index space with geometry -- it **is** part of the surface's description, no
different in kind from `referenceCoordinate` or `radialMin`/`radialMax`. Two
parallel arrays that must always agree in size and index, an identity-bearing
join to validate that agreement, a second currency dimension that must be
bumped in lockstep discipline with the first, and a validity flag needed only
because of the two arrays' potential disagreement -- none of that machinery
is needed once material lives on the same struct as the geometry it describes.

## Decision

`SurfaceDescriptor` (and `StaticSurfaceDescriptor`) directly contain their
nominal material as an appended field:

```cpp
struct NominalSurfaceMaterial {
  float xOverX0{0.f};
  float arealDensityGPerCm2{0.f};
};
```

Both fields describe normal-incidence nominal surface material and are
independently legal at zero (a surface with no material configured yet).
`NominalSurfaceMaterial` is defined exactly once (in `SurfaceDescriptor.h`)
and shared by both descriptor types -- there is no duplicate definition, no
identity join, no material pointer, and no separate epoch or currentness
mechanism for material.

Consequences of putting material directly on the descriptor:

- `SurfaceDescriptor` grows from 20 to 28 bytes (`StaticSurfaceDescriptor`
  from 32 to 36); `DetectorLayoutView` shrinks back down (no material
  pointer, no status byte) to its simpler pre-material shape. The 20-byte
  `SurfaceDescriptor` ABI was deliberately revised: preserving it would have
  required exactly the parallel-array machinery this decision removes, which
  produced more complexity than the 8 bytes it would have saved.
- `DetectorLayoutSet` owns **one** persistent surface catalog (unchanged from
  the prior slice) -- there is no second, parallel array to keep in sync,
  and no positional identity-bearing entry type to validate it against.
- `DetectorLayout` continues to own only iteration-specific topology; surface
  geometry-and-material and the full-catalogue cylinder/disk masks live once
  in `DetectorLayoutSet`, exactly as before.
- `DetectorLayoutView` accesses material through its `SurfaceDescriptor`
  entries (`layoutView.getSurface(id).material`), not through a parallel
  array or a dedicated accessor.
- Changing nominal material is changing the surface description, so it goes
  through the existing detector-layout invalidation/geometry-epoch mechanism
  (`TimeFrame::invalidateDetectorLayouts()`/`DetectorGeometryEpoch`) exactly
  like any other geometry change. There is no renamed epoch and no second
  invalidation method for material specifically.
- `DetectorLayoutViewStatus`/`isValid()` is removed along with it: it existed
  only to support the discarded parallel-material view contract and a
  theoretical valid-empty iteration that the production geometry contract
  (which rejects an empty detector, `DetectorSurfaceCatalogValidationError::
  EmptyDetector`) never actually produces. The default-constructed
  `DetectorLayoutView{}` (`surfaces == nullptr`) is the sentinel, matching
  every existing production/test call site's own null check.
- `SurfaceSpec` validation requires both material fields to be finite and
  non-negative, independently: zero is legal for either field on its own (a
  cylinder with `xOverX0 == 0.f` and a nonzero areal density, or the
  reverse, are both valid configurations); only negative, NaN or Inf values
  are rejected.

Runtime, endpoint-dependent `IntegratedMaterialBudget` (a query result
between two actual track-state endpoints, `MaterialPhysics.h`) remains a
query result and is never stored as the nominal surface property -- that
distinction is unchanged from the superseded design and is not revisited
here.

## Alternatives rejected

- **Parallel `NominalSurfaceMaterialBudget` array owned by `DetectorLayoutSet`,
  keyed by position against the surface catalog** (the prior accepted
  design): superseded by this decision. It required an identity-bearing
  construction/validation type, a second epoch, and a view-validity flag
  purely to manage the two arrays' agreement -- complexity in service of
  keeping the descriptor small, which this decision judges not worth it.
- **Standalone `NominalMaterialCatalog` owner with its own key and
  `isCurrentFor()`**: an earlier-still alternative to the parallel-array
  design above; rejected for the same reason and superseded along with it.
- **Material folded into `SurfaceCatalogView`**: breaks its documented
  narrow, geometry-only contract for every existing consumer.
- **Combined geometry+material provider abstraction now**: premature without
  a second, real (TGeo/LUT) material source to validate the abstraction
  against; still deferred.

## Non-goals

Concrete ITS/MFT material numbers; a compile-time `MaterialSpec`; a TGeo/LUT-
backed material provider or query; nominal-to-integrated conversion or
incidence-angle scaling; `TransitionPolicyBinding`/`TrackerTraits` material
consumption or wiring; propagation or candidate changes; any GPU-readiness
or CUDA/HIP-validation claim; any Stage-B completion claim.

## Validation gates

- POD layout: standard-layout, trivially-copyable, exact `sizeof`/`alignof`/
  field-offset assertions for `NominalSurfaceMaterial`, `SurfaceDescriptor`,
  `StaticSurfaceDescriptor` and `DetectorLayoutView`.
- Static-to-runtime projection (`toRuntimeSurfaceDescriptor()`) copies both
  material fields; independently zero fields are legal; negative/NaN/Inf are
  rejected by `SurfaceSpec` validation.
- `DetectorLayoutSet` iteration views all point to the same shared
  `SurfaceDescriptor` catalog, and material values are visible identically
  from every iteration's view; topology remains iteration-specific; no
  persistent per-iteration surface copy returns.
- `TimeFrame` layout construction and `wipe()`/currentness behavior are
  unchanged by this simplification.
- Full existing common ITSMFT tracking suite, plus the ITS
  `CommonTrackingParameters` and MFT publication-decision tests, remain
  green.

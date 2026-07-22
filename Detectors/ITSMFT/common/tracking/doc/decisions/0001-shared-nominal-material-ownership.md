# ADR 0001: Own nominal per-surface material once, alongside immutable surface geometry

Status: Accepted
Date: 2026-07-22

## Context

Nominal per-surface material (a normal-incidence `xOverX0` and areal density
per surface) needs a runtime home so it can eventually be queried the same
way surface geometry already is. The first candidate design gave it a fully
standalone owner with its own key, its own currentness check, and its own
identity validation, independent of `DetectorLayoutSet`. That duplicates
ownership machinery `DetectorLayoutSet` already has for geometry, and creates
two independently configured owners that could drift out of sync with each
other.

## Decision

Immutable surface geometry and nominal per-surface material form **one
shared detector-description ownership unit**, owned once by
`DetectorLayoutSet`:

- `std::vector<SurfaceDescriptor>` (existing) and
  `std::vector<NominalSurfaceMaterialBudget>` (new) are parallel, same-index
  arrays over the same dense `SurfaceId` space.
- The full-catalogue cylinder and disk masks are derived once from that
  shared array, not recomputed per iteration.
- `materialEpoch` is one more field on the existing
  `DetectorLayoutConfigurationKey`/currentness mechanism, not a second key.

Topology and transitions remain strictly **iteration-specific**, owned by
`DetectorLayout`. `DetectorLayout` instances do not own a surface, material,
or mask copy; they borrow a view of the shared arrays only while validating
at construction time.

`DetectorLayoutView` composes the shared geometry/material with the selected
iteration's topology at view-assembly time
(`DetectorLayoutSet::getLayoutView()`), and carries an explicit
`DetectorLayoutViewStatus` so a legitimately empty-but-valid view can never
be confused with the invalid/no-such-iteration sentinel.

`SurfaceCatalogView` remains geometry-only and is unchanged.

Runtime, endpoint-dependent `IntegratedMaterialBudget` (a query result
between two actual track-state endpoints) is never stored in
`DetectorLayout`/`DetectorLayoutSet`/`TimeFrame`. The shared nominal value is
its explicit, endpoint-independent fallback input — that consumption path is
not implemented in this slice.

## Alternatives rejected

- **Standalone `NominalMaterialCatalog` owner** with its own key and
  `isCurrentFor()`: duplicates `DetectorLayoutSet`'s role and risks drifting
  from it independently.
- **Material field on `SurfaceDescriptor`**: widens a locked, static-asserted
  ABI (20 bytes) for every geometry-only consumer.
- **Material folded into `SurfaceCatalogView`**: breaks its documented
  narrow, geometry-only contract for every existing consumer.
- **Combined geometry+material provider abstraction now**: premature without
  a second, real (TGeo/LUT) material source to validate the abstraction
  against; deferred.

## Consequences

- One allocation of the shared geometry+material+mask arrays per
  `DetectorLayoutSet` rebuild, not one per iteration.
- `DetectorLayoutView`'s ABI changes exactly once, fully audited and
  statically asserted, with an explicit validity discriminator.
- Geometry-epoch and material-epoch currency share one key and one
  invalidation path. Changing material content without bumping
  `materialEpoch` is a documented caller-contract violation the framework
  cannot independently detect (symmetric to the existing `geometryEpoch`
  contract).
- A wrap of either epoch drops the one shared stored owner, since both
  epochs now gate the same object.

## Non-goals

`TransitionPolicyBinding`/`TrackerTraits` material consumption;
`TrackingParameters::LayerxX0` replacement; a compile-time `MaterialSpec`;
a TGeo/LUT-backed material provider; nominal-to-integrated conversion or
incidence-angle scaling; any GPU-readiness or CUDA/HIP-validation claim.

## Validation gates

- POD layout: standard-layout, trivially-copyable, exact `sizeof`/`alignof`/
  field-offset assertions for `NominalSurfaceMaterialBudget` and the final
  `DetectorLayoutView`.
- Construction: max-surface-count, dense `SurfaceId`, positional
  identity-bearing material-entry validation, finite/non-negative budget
  fields, all existing `DetectorLayout` validation preserved.
- Currency: geometry-only and material-only epoch changes independently
  invalidate the stored owner; a failed rebuild leaves the previous owner
  stored but exposed as stale; `wipe()` preserves both epochs and the owner.
- Full existing common ITSMFT tracking suite remains green after this slice.

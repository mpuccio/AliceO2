# ADR 0005: static-plan ownership activation (Gate 4 B2 Slice 2)

Status: Accepted
Date: 2026-07-31

## Context

Gate 4 B2 Slice 1 added compile-time static surface catalogs
(`kITSStaticSurfaceCatalog`, `kMFTStaticSurfaceCatalog`,
`kITSMFTCombinedStaticSurfaceCatalog`) as inert, unwired projections. This
slice activates the per-detector static catalogs as the sole plan source
for production tracking: `ITSMFTTrackingInterface<NLayers>` now builds one
`DetectorLayoutSet` per instance, once, during `initialiseTracker()`, from
`kITSStaticSurfaceCatalog`/`kMFTStaticSurfaceCatalog`, and everything
downstream (`Tracker`, `TrackerTraits`, `TimeFrame::loadNormalizedSource`)
receives that plan explicitly instead of pulling it from a runtime
provider through `TimeFrame`.

The combined 17-surface catalog remains unactivated: `Tracker<7>`/
`Tracker<10>` per-detector workflows are still the only production
consumers, and mixed-detector tracking stays fail-closed, exactly as the
accepted B2 design requires.

## Decision

- `DetectorLayoutSet`/`DetectorLayoutBuilder` borrow a `SurfaceCatalogView`
  instead of owning/copying surfaces; `DetectorLayoutConfigurationKey`
  drops `geometryEpoch`/`catalogRequest` (nothing left to key a runtime
  fetch on).
- New pure function `buildDetectorLayoutSet(SurfaceCatalogView,
  orderedSurfaces, TransitionPolicyTag, trackingParameters)`
  (`DetectorLayoutSet.cxx`) builds a whole plan in one call, with no
  `TimeFrame` or provider involvement.
- `TimeFrame<NLayers>` becomes event-data only: it owns no catalog, no
  layout, no geometry epoch, and exposes no ensure/invalidate/current/get
  accessor for any of them. `loadNormalizedSource()` takes an explicit
  ordered-surfaces span and `SurfaceCatalogView` from its caller.
- `TrackerTraits<NLayers>::initialiseTimeFrame` takes the whole
  `DetectorLayoutSet` explicitly; `Tracker<NLayers>` gains
  `adoptDetectorLayoutSet()`, mirroring `adoptTimeFrame`'s bind-once
  pattern.
- `ITSMFTTrackingInterface<NLayers>` builds one immutable
  `std::optional<DetectorLayoutSet>` unconditionally whenever tracking
  parameters are configured, using a compile-time-selected static catalog
  and a trivial identity surface order (`orderedSurfaces[i] ==
  SurfaceId{i}`, true for both dense/local static catalogs). Build failure
  is fatal.
- `CATrackerSpec.cxx` (ITS and MFT workflows): the CCDB `GEOMTGEO` handler
  no longer invalidates or rebuilds any layout -- the static plan is
  immune to alignment updates by design. `GeometryTGeo::adopt`/
  `fillMatrixCache` are unchanged, since raw cluster decoding still needs
  them.

## Non-goals

Deleting `DetectorSurfaceCatalogProvider`/`GeometrySurfaceCatalogProvider`/
`CombinedSurfaceCatalogBuilder` or any of their tests (a later
provider/aggregation cleanup slice, per the accepted B2 sequencing: 1.
static per-detector + combined projection/tests; 2. owner/`TimeFrame` API
migration [this slice]; 3. provider/aggregation/test cleanup). Activating
the combined 17-surface catalog in any `Tracker` or workflow. Any change to
mixed-track rejection, `TransitionPolicyTag`, output adapters, tracking
physics, or workflow defaults.

## Validation gates

- `O2lib-ITSMFTTracking` and all direct ITS/MFT workflow consumers build
  clean.
- `ctest -L itsmft --output-on-failure`: 71/71 passed, both immediately
  after migration and again after a clang-format cleanup pass (which
  fixed several pre-existing, unrelated formatting violations found
  incidentally in touched files).
- `git diff --check`: clean; all touched `.cxx`/`.h` files pass
  `clang-format -style=file --dry-run -Werror`.
- Device-header check: `nvcc`/`hipcc` both absent from `PATH` in this
  environment -- recorded as the unavailable-toolchain limitation, not run.
- Detached validator (`o2-itsmft-nominal-geometry-validator`, both ITS/7
  and MFT/10, `--format json`): exit 0, byte-identical to the recorded
  Gate 4 B1 C1 lossless JSON. Expected and uninformative about the
  migration itself -- the validator touches no `DetectorLayoutSet`/
  `ITSMFTTrackingInterface` code path -- but re-verified anyway.
- Canonical ITS (203-track, `ee7f7c794d60f2362fd2564258b7887e`) and MFT
  (70-track, `24737e73b7146bf3bd35a90a2517c527`) common-CA replay,
  identical protocol to `gate3-slice3-its-ca-validation/manifest.json`:
  **MATCH** for both track count and content hash, output file sizes
  identical too (55007 B / 41683 B). Fixture `checksums.sha256`: 43/43 OK
  before and after both replays. Full record:
  `/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/gate4-b2-slice2-static-plan-activation/README.md`.

## Alternatives rejected

- **Keep `TimeFrame` as the plan owner and only swap its provider for a
  static one**: rejected -- would preserve the epoch/currency/invalidation
  machinery this slice's own task scope requires removing, and would leave
  the plan re-derivable per-TimeFrame instead of built once per tracker
  instance as designed.
- **Pass `TrackerTraits::initialiseTimeFrame` only the active iteration's
  view instead of the whole `DetectorLayoutSet`**: rejected -- the
  function also needs `getConfigurationKey().orderedSurfaces` (material
  staging) and `size()` (bounds check), both only available from the
  whole set.

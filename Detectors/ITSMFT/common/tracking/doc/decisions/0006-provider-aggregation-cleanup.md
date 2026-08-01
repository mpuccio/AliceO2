# ADR 0006: provider/aggregation/test cleanup (Gate 4 B2 Slice 3)

Status: Accepted
Date: 2026-08-01

## Context

Gate 4 B2 Slice 2 (ADR 0005) made `ITSMFTTrackingInterface` build and own
one immutable `DetectorLayoutSet` from the compile-time static catalogs
(`kITSStaticSurfaceCatalog`/`kMFTStaticSurfaceCatalog`, Slice 1)
unconditionally, for both production ITS/MFT workflows. From that point on,
the runtime geometry-backed catalog-provider system it replaced
(`DetectorSurfaceCatalogProvider` and its `ITS`/`MFT`/`Geometry` concrete
implementations, plus `CombinedSurfaceCatalogBuilder`, the runtime
multi-detector concatenation built on top of it) had no remaining caller:
Slice 2 deliberately kept them in the tree, inert, per its own scope
constraint ("do not delete providers ... yet"). This slice removes them now
that the replacement path has proven itself across two gates of build/test/
replay validation.

## Decision

### Deleted outright

`GeometrySurfaceCatalogProvider.{h,cxx}`, `ITSSurfaceCatalogProvider.{h,cxx}`,
`MFTSurfaceCatalogProvider.{h,cxx}`, `DetectorSurfaceCatalogProvider.h` (the
request/result/error API all three implemented),
`CombinedSurfaceCatalogBuilder.{h,cxx}`, `testCombinedSurfaceCatalogBuilder.cxx`,
`testGeometrySurfaceCatalogProviders.cxx`. `ITSMFTTrackingInterface`'s
constructor `catalogProvider` parameter and `mDetectorSurfaceCatalogProvider`
member (already unused since Slice 2) are removed with
`DetectorSurfaceCatalogProvider.h`; both ITS/MFT `CATrackerSpec.h` workflow
specs updated to the plain 3-arg constructor. `DetectorLayoutSetBuildResult`'s
`catalogError`/`catalogValidationError` fields and the
`DetectorSurfaceCatalogValidationError` enum, plus the
`MissingProvider`/`CatalogProviderFailure`/`InvalidCatalog` variants of
`DetectorLayoutSetBuildError`, are removed too: all were diagnostics only the
deleted runtime-provider path could ever populate, confirmed dead by grep
(no production or test code read any of them) before removal.

Per the Slice 1/2 deletion-ordering constraint, `CombinedSurfaceCatalogBuilder`
is deleted only now, in the same cleanup slice as the providers it built on
— never independently, and only after its replacement
(`testCombinedStaticSurfaceCatalogTopology.cxx`, landed in Slice 1) was
already in place and passing.

### Moved, not duplicated

`DetectorSurfaceCatalogAggregation.{h,cxx}` — the chip-active-area geometric
aggregation helper (`aggregateSurfaceGeometry()`) both deleted geometry-backed
providers used, and the detached validator's own
`GeometryTGeoSurfaceSource.cxx`/`NominalGeometryReport.cxx` also use — moves
from `Detectors/ITSMFT/common/tracking/src/` to
`Detectors/ITSMFT/common/tracking/validation/src/`, exactly as the Slice 1
validator `CMakeLists.txt` comment anticipated. Its focused test moves with
it into `validation/test/`. It compiles exactly once now, into
`ITSMFTNominalGeometryValidatorCore`; `O2::ITSMFTTracking` no longer lists it,
and the validator's `../src` include-directory workaround is removed.

### Also fixed: a pre-existing Slice 2 regression

`compileCommonTrackerOnboarding.cxx` (a link-only compile proof under
`Detectors/ITSMFT/ITS/tracking/test/`) called `TimeFrame::ensureDetectorLayouts()`/
`getDetectorLayouts()`/`getDetectorLayoutView()`, all removed in Slice 2, and
was independently confirmed broken (3 compile errors) before this slice
touched it — Slice 2's own build validation never covered this target, since
it falls outside "`O2lib-ITSMFTTracking` + direct ITS/MFT workflow
consumers". Rewritten onto the modern API: builds one `DetectorLayoutSet` via
`buildDetectorLayoutSet()` directly from `kITSStaticSurfaceCatalog` (no fake
provider needed any more, since the real static catalog is itself
synthetic/compile-time and requires no runtime geometry file), extracts
per-layer material the same way `TrackerTraits::initialiseTimeFrame()` does,
and drives `Tracker<7>` via `adoptDetectorLayoutSet()`.

## Non-goals

Activating the combined 17-surface catalog in any `Tracker` or workflow. Any
change to mixed-track rejection, `TransitionPolicyTag`, output adapters,
tracking physics, or workflow defaults.

## Validation gates

- `O2lib-ITSMFTTracking`, `O2lib-ITSMFTNominalGeometryValidatorCore`,
  `o2-itsmft-nominal-geometry-validator`, the onboarding-compile executable,
  and both direct ITS/MFT workflow consumers (`o2-its-ca-tracker-workflow`,
  `o2-mft-ca-tracker-workflow`) all build clean.
- Every `O2test-itsmft*` target (58) rebuilt explicitly and links clean.
- `ctest -L itsmft --output-on-failure`: 69/69 passed (71 in Slice 2, minus
  the 2 deleted test files).
- `git diff --check` clean; all touched files pass
  `clang-format -style=file --dry-run -Werror`.
- No buildable source/CMake/runtime-code reference to any deleted symbol
  remains (grep-verified across the whole tree).
- Device-header check: `nvcc`/`hipcc` both absent — recorded as the
  unavailable-toolchain limitation, not run.
- Detached validator (ITS/7, MFT/10, `--format json`): byte-identical to the
  recorded Gate 4 B1 C1 lossless JSON, as expected (the aggregation move is a
  verbatim relocation, no logic change).
- Canonical ITS (203-track, `ee7f7c794d60f2362fd2564258b7887e`) and MFT
  (70-track, `24737e73b7146bf3bd35a90a2517c527`) common-CA replay: **MATCH**
  for both track count and content hash, output file sizes identical to
  every prior gate (55007 B / 41683 B). Fixture `checksums.sha256`: 43/43 OK
  before and after. Full record:
  `/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/gate4-b2-slice3-provider-cleanup/README.md`.

## Alternatives rejected

- **Leave `catalogError`/`catalogValidationError` in place for wire/API
  stability**: rejected — unlike the `MultiSourceLoadError` enum values kept
  for genuine wire compatibility elsewhere in this codebase,
  `DetectorLayoutSetBuildResult` is purely an in-process return type with a
  single caller (`ITSMFTTrackingInterface::initialiseTracker()`); there is no
  serialized/cross-process representation to stay compatible with, so
  keeping permanently-`None` fields around would be pure clutter.
- **Delete `compileCommonTrackerOnboarding.cxx` instead of rewriting it**:
  rejected — it is not on this slice's explicit deletion list, and rewriting
  it onto the modern API both fixes the regression and keeps its original
  purpose (a link-only proof that the common-core onboarding path compiles
  end-to-end) intact.

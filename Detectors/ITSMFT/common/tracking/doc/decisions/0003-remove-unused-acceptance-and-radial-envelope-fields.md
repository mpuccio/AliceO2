# ADR 0003: remove `SurfaceAcceptance`/`radialMin`/`radialMax` from the permanent descriptor contract

Status: Accepted
Date: 2026-07-31

## Context

The Gate 4 B1 Slice 2 static `DetectorSpec` authoring proposal exposed that
`SurfaceAcceptance`/`CylinderZAcceptance`/`DiskRadialAcceptance`
(`StaticSurfaceDescriptor::nominalTrackingAcceptance`) and
`SurfaceDescriptor::radialMin`/`radialMax` have no tracking-runtime
consumer. Authoring real ITS/MFT static tables against the unmodified
shape would have forced either inventing an ITS z-acceptance with no
source and no consumer, or freezing MFT's chip-envelope radii into a field
presented as permanent tracker semantics, purely to satisfy an existing
compile-time concept -- not to serve any real tracking need. This decision
removes both instead of working around them.

A pre-deletion audit (grep across the whole repository) confirmed every
read/write of these items:

- `SurfaceAcceptance`/`CylinderZAcceptance`/`DiskRadialAcceptance`/
  `SurfaceAcceptanceKind`: only in `StaticSurfaceDescriptor.h` (definition)
  and `testSurfaceSpec.cxx` (its own test).
- `StaticSurfaceDescriptor::nominalTrackingAcceptance`: written/read only in
  `StaticSurfaceDescriptor.h` (`toRuntimeSurfaceDescriptor()`'s now-removed
  Disk-only copy branch), `SurfaceSpec.h`'s compile-time validator, and
  `testSurfaceSpec.cxx`.
- `SurfaceDescriptor::radialMin`/`radialMax`: written only by
  `GeometrySurfaceCatalogProvider.cxx`'s catalog construction and
  `StaticSurfaceDescriptor.h`'s projection; read only by ABI asserts and by
  tests of that same construction/projection code
  (`testGeometrySurfaceCatalogProviders.cxx`, `testSurfaceSpec.cxx`) and by
  the (now-deleted) manual `GeometrySurfaceCatalogValidation.h` tool.
  **Zero references** in `TimeFrame`, `TrackerTraits`,
  `DetectorLayoutBuilder`/`DetectorLayout.h` (`computeSurfaceKindMasks()`
  reads only `.kind`), `SparseTrackingTopology`,
  `TransitionPolicyOperations`, or any material/refit code.

`AggregatedSurfaceGeometry`'s own `radialMin`/`radialMax`
(`DetectorSurfaceCatalogAggregation.h`) are a different type and are
**not** part of this decision -- they remain the detached Gate 4 B1
nominal-geometry validator's diagnostic output, unchanged.

## Decision

- Delete `SurfaceAcceptance`, `SurfaceAcceptanceKind`,
  `CylinderZAcceptance`, `DiskRadialAcceptance`, and
  `StaticSurfaceDescriptor::nominalTrackingAcceptance`.
- Delete `SurfaceDescriptor::radialMin`/`radialMax`. New layout:
  `sizeof(SurfaceDescriptor) == 20` (was 28), `offsetof(..., material) ==
  12` (was 20).
- `StaticSurfaceDescriptor` shrinks to `{id, identity, kind,
  nominalReferenceCoordinate, material}`:
  `sizeof(StaticSurfaceDescriptor) == 24` (was 36),
  `offsetof(..., material) == 12` (was 24), `offsetof(...,
  material) == 12` (was 32).
- `toRuntimeSurfaceDescriptor()` simplifies to a flat, unconditional field
  copy -- no per-kind branch left, since there is nothing left to branch
  on.
- `SurfaceSpec.h`'s `validateSurfaceArray()` drops the acceptance-kind/
  bounds validation block entirely.
- `GeometrySurfaceCatalogProvider.cxx` stops forwarding
  `aggregation.surfaces[i].radialMin`/`radialMax` into `SurfaceDescriptor`
  construction; it still calls `aggregateSurfaceGeometry()` unchanged and
  still reads `.referenceCoordinate`.
- No replacement acceptance abstraction is introduced. This is a pure
  subtraction.
- The superseded old geometry-surface-catalog validator (an
  `O2::ITSMFTTracking`-linked executable + ROOT macro, predating and
  superseded by the Gate 4 B1 Slice 1 detached validator) is deleted
  entirely, not adapted to the new shape:
  `validateGeometrySurfaceCatalogs.cxx`, `validateGeometrySurfaceCatalogs.C`,
  `GeometrySurfaceCatalogValidation.h`, `SliceB2Validation.md`, and both
  CMake registrations (`o2_add_executable(geometry-surface-catalog-validation
  ...)`, `o2_add_test_root_macro(validateGeometrySurfaceCatalogs.C ...)`).

## Non-goals

Static `ITSSurfaceSpec`/`MFTSurfaceSpec` table authoring (blocked
separately on detector-owner sign-off); B2 static-plan activation or
deletion of the still-live runtime catalog providers/epoch mechanism
(`GeometrySurfaceCatalogProvider`, `ITSSurfaceCatalogProvider`,
`MFTSurfaceCatalogProvider`, `DetectorSurfaceCatalogProvider`,
`DetectorGeometryEpoch`) -- all untouched, deletable only once a static
plan replaces them; any change to `TimeFrame`, tracking algorithms,
workflows, CCDB, or alignment behavior; the detached Gate 4 B1 validator
(`validation/`), left fully unchanged.

## Validation gates

- Repository-wide grep: zero remaining references to any deleted symbol
  outside this decision's own edits; zero remaining references anywhere to
  the retired validator's names.
- `O2lib-ITSMFTTracking` and its direct workflow consumers
  (`o2-its-ca-tracker-workflow`, `o2-mft-ca-tracker-workflow`) build clean.
- Full ITSMFT ctest suite (69 tests: common tracking, the detached
  validator, ITS/MFT workflow, `gate0-baseline`/`gate3-slice3-its-ca-validation`/
  `cluster-attachment-validation` compiled macros) passes 100%.
- Detached validator: all its own tests still pass, and both real
  checksum-verified-fixture `--format json` invocations (ITS, MFT) still
  exit 0 and are byte-identical to the pre-decision recorded output --
  proving this decision has zero effect on that separate tool.
- GPU/device-header: `SurfaceDescriptor.h`/`StaticSurfaceDescriptor.h`/
  `IdTypes.h` are not included by any actual GPU-compiled target in this
  repository; no CUDA/HIP/OpenCL toolchain is available in this
  environment to exercise `GPUCommonDefAPI.h`'s device branches directly.
  The full host build already exercises the one branch reachable without a
  real GPU-compiler frontend (`!defined(GPUCA_GPUCODE)`); every surviving
  `GPUhdi()` annotation in both edited headers was manually re-verified
  intact.
- Canonical ITS (203-track, `ee7f7c794d60f2362fd2564258b7887e`) and MFT
  (70-track, `24737e73b7146bf3bd35a90a2517c527`) common-CA replay:
  **MATCH, confirmed**. The first attempt failed on a `TGrid::Connect`/
  AliEn-token environment error (`internal-dpl-ccdb-backend`; no AliEn
  certificate configured in that session) -- upstream of and unrelated to
  any file this decision touches. Once a valid AliEn token was available
  in the environment, the identical replay commands, from fresh output
  directories, against the same checksum-verified fixture (43/43 OK before
  and after), reproduced both canonical hashes exactly. Full record in
  `/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/gate4-acceptance-cleanup-c2/README.md`.

## Alternatives rejected

- **Keep the fields, populate them with placeholder/best-guess values**:
  rejected outright by the originating instruction -- inventing an ITS
  z-acceptance with no source, or presenting MFT's chip-envelope radii as
  permanent tracker semantics, was exactly the problem this decision
  exists to avoid.
- **Introduce a generic replacement acceptance abstraction**: rejected;
  explicitly out of scope, and nothing currently needs one.
- **Modify the old geometry-surface-catalog validator to drop its
  acceptance/radial checks instead of deleting it**: rejected; it predates
  and is fully superseded by the Gate 4 B1 Slice 1 detached validator, so
  keeping a second, `O2::ITSMFTTracking`-linked tool alive would preserve
  exactly the coupling the detached validator exists to avoid, for no
  remaining benefit.

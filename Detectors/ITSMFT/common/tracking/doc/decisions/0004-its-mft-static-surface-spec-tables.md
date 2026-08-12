# ADR 0004: candidate static `ITSSurfaceSpec`/`MFTSurfaceSpec` tables (Slice C3a)

Status: Accepted
Date: 2026-07-31

## Context

Gate 4 acceptance-cleanup C2 (ADR 0003) removed
`SurfaceAcceptance`/`radialMin`/`radialMax` from `StaticSurfaceDescriptor`/
`SurfaceDescriptor`, closing the two gaps the Gate 4 B1 Slice 2 authoring
report identified (an ITS z-acceptance with no source, and MFT
chip-envelope radii presented as permanent tracker semantics). With those
fields gone, authoring a candidate static table no longer requires
inventing or freezing either.

This slice authors those tables: `ITSSurfaceSpec` (7 surfaces) and
`MFTSurfaceSpec` (10 surfaces), header-only, using only the existing
`StaticSurfaceDescriptor`/`SurfaceSpec`/`toRuntimeSurfaceDescriptor()`
mechanism. Both are candidates awaiting detector-owner sign-off (per the
Slice 2 report's checklist), not yet activated: no runtime code
constructs a catalog from either, and no `TimeFrame`/`DetectorLayoutSet`
consumes them.

## Decision

### Source authority

- Geometry: the checksum-pinned unaligned
  `o2sim_geometry.root` (SHA-256
  `2a428746b3a0b57179d5ffe631afc9c4afb4ca41cc9baa948ff670099b9204e4`,
  fixture `pp-20ev-run303000-seed20260716-daily20260717`, package
  `daily-20260717-0700-local1`) the Gate 4 B1 Slice 1 detached validator
  already uses -- unchanged, re-verified in this slice.
- `nominalReferenceCoordinate` literals: transcribed verbatim (not
  re-rounded) from that validator's `--format json`
  `formatLosslessFloat()` output (Gate 4 acceptance-cleanup C1's lossless
  `std::to_chars` encoding), specifically
  `O2-validation-artifacts/itsmft/gate4-b1-slice1-nominal-geometry-validation/pp-20ev-run303000-seed20260716-daily20260717/acceptance-cleanup-c1-lossless-json/{its,mft}-report.json`.
  Not the earlier fixed-six-decimal output.
- `material`: `kNominalITSLayerX0`/`kNominalMFTLayerX0`
  (`ITSMFTDetectorDefinitions.h`), referenced by index, not
  re-literaled; `arealDensityGPerCm2` derived via the existing
  `xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho` formula
  (`ITStracking/Constants.h`), the same one
  `GeometrySurfaceCatalogProvider.cxx` already uses. No new detector
  constant is introduced anywhere in this decision.

### Shape

`include/ITSMFTTracking/ITSSurfaceSpec.h` / `MFTSurfaceSpec.h`, each a
single `struct { inline static constexpr std::array<StaticSurfaceDescriptor,
N> surfaces{...}; }` satisfying the `SurfaceSpec` concept, plus
`static_assert(SurfaceSpec<...>)`/`static_assert(SurfaceCount<...> == N)`.
No new wrapper type; existing mechanism only.

- `id`: dense, local `SurfaceId{0..N-1}` within each spec (ITS 0-6, MFT
  0-9) -- global rebasing across detectors is `ConcatenatedSurfaceSpec`'s
  job, exercised generically already in `testSurfaceSpec.cxx`, not
  performed here.
- `identity`: `{detectorId, detectorSurfaceIndex}` --
  `static_cast<uint8_t>(o2::detectors::DetID::ITS)` (0) /
  `o2::detectors::DetID::MFT` (8), local index equal to the surface's
  position.
- `kind`: `Cylinder` for every ITS
  surface, `Disk`/`CartesianXY` for every MFT surface.

## Non-goals

Static-plan activation (no `DetectorLayoutSet`/`TimeFrame`/
`ITSMFTTrackingInterface` change); deleting or otherwise touching the
still-live runtime catalog providers/epoch mechanism
(`GeometrySurfaceCatalogProvider`, `ITSSurfaceCatalogProvider`,
`MFTSurfaceCatalogProvider`, `DetectorSurfaceCatalogProvider`,
`DetectorGeometryEpoch`); any acceptance/radial-envelope field (removed in
ADR 0003, not reintroduced); a geometry-comparison tolerance or validator
comparison mode (the detached validator gains no new mode in this slice);
any `TimeFrame`, tracking-algorithm, workflow, CCDB, or alignment change;
detector-owner sign-off itself (still outstanding -- these are candidate
tables, not approved ones).

## Validation gates

- `O2lib-ITSMFTTracking` builds clean with both new headers included.
- New focused test, `testITSMFTSurfaceSpecProjection.cxx`: for all 17
  authored surfaces (a) the compiled `nominalReferenceCoordinate` literal
  reproduces the exact `float32` bit pattern of its source JSON token
  (`strtof` + `memcpy`'d `uint32_t` ==, zero tolerance -- not a geometry
  comparison, no geometry file touched); (b) `id`/`identity`/`kind`/
  `kind` matches the authored table exactly; (c) `material`
  matches `kNominalITSLayerX0`/`kNominalMFTLayerX0` and the `Radl*Rho`
  formula; (d) `toRuntimeSurfaceDescriptor()` preserves every field,
  including the reference-coordinate bit pattern, bit-exactly.
  `SurfaceSpecsCanBeConcatenated<ITSSurfaceSpec, MFTSurfaceSpec>` also
  holds (compile-time only, not exercised at runtime here).
- Full ITSMFT ctest suite: 70/70 passed (69 pre-existing + the new test).
- Detached validator re-run, both detectors, `--format json`: exit 0,
  byte-identical to the recorded C1 lossless output, against the
  checksum-reverified (43/43 OK before and after) fixture.
- Canonical ITS (203-track, `ee7f7c794d60f2362fd2564258b7887e`) and MFT
  (70-track, `24737e73b7146bf3bd35a90a2517c527`) common-CA replay,
  identical protocol to `gate3-slice3-its-ca-validation/manifest.json`:
  **MATCH**, confirming these new, unwired headers have zero effect on
  tracking output. Full record:
  `/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/gate4-b1-slice-c3a-static-spec-tables/README.md`.

## Alternatives rejected

- **Re-derive `nominalReferenceCoordinate` from the geometry file at test
  time (a mini validator-comparison mode)**: rejected; explicitly out of
  scope for this slice, and unnecessary -- the bit-exact
  literal-vs-JSON-token comparison already proves transcription fidelity
  without loading geometry a second time.
- **Re-literal the material constants instead of referencing
  `kNominalITSLayerX0`/`kNominalMFTLayerX0` by index**: rejected; would
  duplicate numbers that already exist and could silently drift from the
  source array.

# Header and graph cleanup: Campaign B validation

- Status: accepted
- Date: 2026-08-07
- Gate A base: `2bb3760a6b`
- Accepted head: `cb98190f6d`
- Package: `daily-20260717-0700-local1`
- Integration build: `O2-worktree-builds/itsmft-header-graph-cleanup`
- Frozen build: `O2-worktree-builds/itsmft-header-graph-cleanup-frozen-86113f2a14`

## Integrated slices

| Slice | Commit | Retired headers |
|---|---|---|
| descriptors and catalogs | `681cfceb3f` | `SurfaceCatalogView.h`, `StaticSurfaceDescriptor.h`, `ITSSurfaceSpec.h`, `MFTSurfaceSpec.h` |
| primitive identifiers and masks | `64f770d3b1` | `SurfaceMeasurementIndex.h`, `LayerMask.h` |
| state, refit, and seed-anchor gate | `cb98190f6d` | `SurfaceLinearizationReference.h`, `RefitLegAssembly.h`, `SeedAnchor.h` |

The repository and available downstream O2Physics, O2DPG, alidist, and
Strangeness-Tracker searches found no live `SeedAnchor` caller. The two local
anchor-taking implementations and their tests were the deletion set identified
by audit 0011, not production consumers. The anchor-less outer seed formulas
and independent physics fixtures remain unchanged.

No forwarding header, compatibility alias, retired production include, or
retired `SeedAnchor` failure value remains. The common tracking inventory after
this campaign is 52 public and four `detail/` headers.

Each slice passed `git diff --check`, pinned-environment
`git clang-format --diff`, changed-header compilation, a focused Ninja build,
and source-name-focused tests after its separate cherry-pick. The combined
Campaign B build covered common tracking, ITS, MFT, the combined workflow, CA
writer, Framework analysis/CCDB support, and every directory contributing an
`itsmft`-labelled test.

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 93
No Not Run entries.
```

The strict preflight passed with a valid AliEn token,
`ALICEO2_CCDB_NOTOKENCHECK` unset, the exact pinned package, and all reported
Geant4 data sets present. No CUDA or HIP compiler is available, so no real GPU
build is claimed.

## Frozen replay gate

Fresh outputs were created under
`/private/tmp/o2-itsmft-header-graph-cleanup-gate-b-20260807`. The fixture
manifest passed 43/43 checksums before and after replay.

The first feature and frozen ITS attempts exposed staged-library ABI
incompleteness: each build contained a local `ITStracking` library while the
reader loaded the pinned package's `ITSTrackingInterface`. Both attempts failed
before cluster publication or tracking. Building the matching interface from
each exact source tree resolved the symbol mismatch; retries were run in new
directories, and the failed directories were retained unchanged.

| Leg | Tracks | Content hash |
|---|---:|---|
| ITS feature standalone and combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS frozen standalone and combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT feature standalone and combined | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT frozen standalone and combined | 68 | `8106b08571ca593c6b76ff72b761a680` |

The complete extractor JSON was byte-identical across feature/frozen and
standalone/combined dimensions. Eight persisted-product comparisons reported
field-by-field equality for initialized writer output, ROFs, cluster indices,
labels, references, `CommonTrack` content, and sidecars. Every MFT comparison
covered 2,992 projected values with maximum absolute and relative delta zero;
only the established undefined `MFTTrack.mInvQPtSeed` byte remains excluded.

The named user stash remained the unchanged object
`6a90bbcd7e187673a7eeaedc2f8df07c471c09b4` throughout the gate.

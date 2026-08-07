# Header and graph cleanup: Campaign C validation

- Status: accepted
- Date: 2026-08-07
- Gate B base: `9403ee309e`
- Accepted head: `fa95ed125d`
- Package: `daily-20260717-0700-local1`
- Integration build: `O2-worktree-builds/itsmft-header-graph-cleanup`
- Frozen build: `O2-worktree-builds/itsmft-header-graph-cleanup-frozen-86113f2a14`

## Integrated slices

| Slice | Commit | Retired public headers |
|---|---|---|
| clock publication | `8d2df482f2` | `ClockTimingPublicationView.h` |
| ROF timing | `25212637d8` | `ROFTimingUniformity.h` |
| host decoding boundary | `30a22d9260` | `ClusterDecoder.h`, `DecodedCluster.h`, `SurfaceMeasurementAdapters.h` |
| multi-source Loader boundary | `7a906060b3` | `ClusterSource.h`, `MultiSourceLoading.h`, `TimeFrameLoadFailure.h` |
| implementation-header relocation | `fa95ed125d` | eight public paths moved to `detail/` |

The Loader remains the one host boundary that classifies source failures and
stages decoded measurements. `TimeFrame` still owns the tracking workspace,
and successful loads still replace that workspace atomically. The Loader does
not become a dependency of the workspace. Live legacy `IOUtils` entry points
remain available.

The relocation moved `CommonTrackShadow.h`, `SurfaceTrackingScratch.h`,
`DetectorPublicationAdapter.h`, `DetectorTrackingOperationAdapterSupport.h`,
`ITSSharedClusterCompatibility.h`, `MFTPublicationCompatibility.h`,
`MFTFwdTrackHelpers.h`, and `SurfaceKinematicStateLegacyAdapters.h` to
`ITSMFTTracking/detail/`. There are no forwarding headers at the old paths.
Production workflows do not include the relocated workspace or compatibility
headers directly. `SurfacePlanBinding` remains a detail implementation type.

The common tracking inventory is exactly 36 public and 12 `detail/` headers.
Each slice passed `git diff --check`, pinned-environment
`git clang-format --diff`, changed-header self-containment, a focused Ninja
build, and source-name-focused tests after its separate cherry-pick. The
combined Campaign C build covered common tracking, ITS, MFT, the combined
workflow, CA writer, tracking interface, Framework analysis/CCDB support, and
every directory contributing an `itsmft`-labelled test.

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 92
No Not Run entries.
```

The strict preflight passed with a valid AliEn token,
`ALICEO2_CCDB_NOTOKENCHECK` unset, the exact pinned package, and all reported
Geant4 data sets present. No CUDA or HIP compiler is available, so no real GPU
build is claimed.

## Frozen replay gate

Fresh outputs were created under
`/private/tmp/o2-itsmft-header-graph-cleanup-gate-c-20260807`. The fixture
manifest passed 43/43 checksums before and after replay.

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

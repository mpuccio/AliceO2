# M7f: final runtime-core cleanup

Status: complete, 2026-08-06
Base: accepted M7e `2026534db8`
Implementation commits: `76170fc7d0`, `700e8b48c8`

## Architectural result

M7f closes the redundant compile-time bridge left after M7e. Common tracking
now has one runtime-plan execution path composed from `DetectorLayoutView`,
`SurfacePlanBinding`, `SurfaceTrackingScratch`, `TrackSeed`, `TrackerTraits`,
`Tracker`, `TimeFrame`, and `CommonTrack`. The generic road-start loop derives
its extent from `SurfaceTrackingScratch::getNOwnedSurfaces()`; it no longer
uses `TrackingParameters::NLayers` through `CellsPerRoad()`.

This slice does not change propagation, covariance sanitization, material,
acceptance, holes, policy grouping, ordering, output format, raw-ROF
ownership, or workflow lifecycle.

## Production deletions and consolidation

| Deleted or reduced residue | Result |
|---|---|
| `TrackerTraits::mSurfaceToSlot`, its staged array, sentinel, and reset/commit logic | Transition and cell endpoints resolve through `SurfacePlanBinding::getOwnedSurfaceIndex()`; the binding is the sole sparse SurfaceId-to-position authority. |
| `IndexTableUtilsN` aliases | `IndexTableUtilsCore` is named directly; its fixed device capacity remains unchanged. |
| `assembleRefitLegSlots<NLayers>` and `fitTrackSeedLegs<NLayers>` | Host-only refit working state is runtime-span/vector based; the three-leg sequence is unchanged. |
| `refitITSSeed<NLayers>` and `refitMFTSeed<NLayers>` | Detector operation adapters call one non-templated helper; `refitDetectorSeed<DetId>` retains only the adapter detector choice. |
| `nativeRefitTrackCylinderCylinder<NLayers>` | The retained native cylinder oracle takes runtime measurement spans. |
| Dead `exportNativeRefitToTrackITSExt<NLayers>` | Typed output conversion remains in the adapter tests/publication boundary; no production consumer existed. |

The shared `TrackingParameters::NeighboursPerRoad()`, `CellsPerRoad()`, and
`TrackletsPerRoad()` accessors remain only because frozen `ITStracking` and GPU
consumers still call this configuration record. The common runtime tracker no
longer calls them; this is an explicitly named adapter/frozen-compatibility
exception, not a count authority.

## Exact residual exception table

The M7f guard scans only common tracking production `include/` and `src/`.
After comments are removed, every live `NLayers`/`ITSNLayers`/`MFTNLayers`
occurrence is classified by exact path as follows.

| Classification | Exact common path(s) | Role and deletion condition |
|---|---|---|
| Narrow adapter-owned compatibility | `include/ITSMFTTracking/Configuration.h`, `src/Configuration.cxx` | Configuration serialization, frozen configuration accessors, and DetID-to-layer construction at the application edge. Delete after all frozen ITS/GPU consumers take runtime plan counts. |
| Narrow adapter-owned compatibility | `include/ITSMFTTracking/TrackingConfigParam.h`, `ITSSurfaceSpec.h`, `MFTSurfaceSpec.h`, `ITSMFTDetectorDefinitions.h`, `ITSMFTDetectorDefinitions.h` | ITS/MFT application constants and static descriptor tables. They select descriptors, not generic traversal behavior. |
| Narrow adapter-owned compatibility | `src/DetectorLayoutSet.cxx`, `src/IndexTableConfiguration.cxx`, `src/TrackerTraits.cxx` | Edge validation that checks the supplied configuration count against the active plan count. The one `TrackerTraits` check is validation only and is not used as a loop or capacity bound. |
| Narrow adapter-owned compatibility | `include/ITSMFTTracking/DetectorPublicationAdapter.h` | Typed publication sidecars at the adapter edge. No typed accepted-track vector crosses into the core. |
| Narrow adapter-owned compatibility | `include/ITSMFTTracking/SurfacePlanTrackingParticipant.h`, `src/SurfacePlanTrackingParticipant.cxx` | ITS/MFT participant, ROF-table construction, typed refit/output edge, and publication sidecars. It owns plain non-templated `Tracker`/`TrackerTraits`; it is not a coordinator. |
| Narrow adapter-owned compatibility | `include/ITSMFTTracking/TrackingInterface.h`, `src/TrackingInterface.cxx` | Workflow-facing ITS/MFT interface and short-lived conversion from frozen detector ROF tables to runtime ROF views. |
| Fixed device ABI/capacity | `Cell.h`, `TrackSeed.h`, `IndexTableUtils.h` | `SeedMetadataBase<ClustersPerCell>`, three-measurement `CellSeed`, fixed `TrackSeed::MaxLayoutSurfaces`, and fixed index-table storage retain device-safe maximum capacity with runtime populated prefixes/masks. |
| Frozen legacy ITS exclusion outside common tracking | `Detectors/ITSMFT/ITS/tracking/{include,src,GPU}` | Not scanned by the common guard and not changed by M7f. The exclusion is this explicit path only; there is no blanket exclusion of common tracking. |
| Dead/redundant | Deleted in commits `76170fc7d0` and `700e8b48c8` | The SurfaceId bridge, aliases, layer-count refit façades, typed native export helper, and stale test call sites are gone. Any new occurrence fails `testM7fRuntimeCoreGuard`. |

No `TrackSeedTpl`, `DetectorTraversalBinding`, `LegacyTrackerScratch`,
`DetectorTraits<...>`, `CATrackType<...>`, `LayerMeasurementSpans<...>`,
`Tracker<...>`, `TrackerTraits<...>`, `CATrackerITS`, `CATrackerMFT`,
`TrackerITS`, or `TrackerMFT` remains in common production sources.

## Guard and coverage

`testM7fRuntimeCoreGuard` is the authoritative M7f source/dependency guard.
It strips comments, scans common production include/src, classifies every
residual layer-count token by exact path and role, rejects the deleted
vocabulary and legacy SurfaceId bridge names, and checks the generic core
orchestration files for detector IDs, ITS/MFT layer constants, ROF table
types, typed output, DPL, and writer dependencies. It also asserts the
runtime binding/scratch authorities and fixed device capacities are still
visible.

The existing sparse/non-identity plan, ROF lifecycle, runtime-count,
non-templated-core, adapter-boundary, failure/reset, CommonTrack, and workflow
contract tests remain registered. Refit tests were migrated to generic
runtime states, spans, masks, and timestamps; typed ITS objects remain only
in the frozen legacy-oracle comparison or adapter-specific compatibility
checks.

## Ranked deletion and simplification inventory

| Rank | Category | Item | Gate |
|---:|---|---|---|
| 1 | Safe during M7 | Remove the persistent SurfaceId-to-slot bridge, `IndexTableUtilsN`, and redundant `<NLayers>` host refit façades. | Completed M7f; behavior-preserving. |
| 2 | Safe during M7 | Keep runtime count authority visible at the road-start and refit-span boundaries and reject any newly introduced generic `.NLayers` loop/capacity use. | Completed by the M7f guard. |
| 3 | Blocked until a separately approved post-M7 owner decision | Reduce `TrackingOperationAdapter` only after accepted-result completion and sidecar ownership can be expressed with one stable seam and replay evidence. | Not changed by M7f. |
| 4 | Blocked until a separately approved post-M7 owner decision | Reduce `SurfacePlanTrackingParticipant<NLayers>` and `ITSMFTTrackingInterface<NLayers>` after all workflow-edge ROF/configuration consumers are independently runtime-view based. | Narrow adapter exception retained. |
| 5 | Deferred physics/algorithm choice | Simplify or remove the native cylinder oracle and MFT leaf compatibility helpers if doing so changes operation selection or numerical characterization. | Requires physics sign-off; not changed. |

M7f does not start another tracker de-templating slice, alter frozen ITS
workflows, or introduce a new wrapper, coordinator, policy taxonomy, or
detector dispatch abstraction.

## Validation record

The existing durable build directory and package were reused:

`/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement`
`daily-20260717-0700-local1`

The complete build passed all 200 scheduled build steps. The required serial
selector executed **98/98** registered ITS/MFT tests and passed every test;
there were no `Not Run` results.

The fixture manifest passed **43/43 before** and **43/43 after** replay.
Standalone and combined replay metrics were:

| Leg | Tracks | Content hash |
|---|---:|---|
| ITS standalone | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT standalone | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT combined | 68 | `8106b08571ca593c6b76ff72b761a680` |

The combined ITS and MFT writer products matched their corresponding
standalone products field-by-field. M7f output matched the accepted M7e
parent for initialized ITS and MFT writer fields, CommonTracks, references,
ROFs, labels, and sidecars. The only excluded value was the known undefined
`MFTTrack.mInvQPtSeed` byte artifact. The MFT projection comparison covered
2,992 values with zero maximum absolute and relative delta.

`nvcc`, `hipcc`, `nvidia-smi`, and `rocminfo` are absent from the pinned
environment; no GPU/device validation is claimed.

The user-owned `TripletTrackingRAndD.md` remains untouched in the named local
stash `stash@{0}` (`user WIP: TripletTrackingRAndD.md`) for later restoration.

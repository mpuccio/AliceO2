# ITSMFT fine-comb boundary cleanup validation

- Status: accepted
- Date: 2026-08-07
- Accepted parent: `e6b114fbfe93b6222343ed71dd12b13db686299d`
- Accepted implementation head: `f953d2b833`
- Package: `daily-20260717-0700-local1`
- Integration build: `O2-worktree-builds/itsmft-fine-comb-cleanup`
- Frozen build: `O2-worktree-builds/itsmft-header-graph-cleanup-frozen-86113f2a14`

## Outcomes

`ConfigKeyValuesPreflight` had one production caller and enforced ITS workflow
configuration policy. Its parser and exact diagnostics now live at that
workflow boundary; the common public header, source, and redundant standalone
test are deleted. Repository and O2Physics searches found no other caller.

Forward-state propagation now has one production boundary:
`Propagator::propagateForward`. It uses the established full Helix transport
for `abs(bz) > 0.01f` and the established Linear transport otherwise. The
Linear fallback is required because the Helix implementation contains
`1 / abs(B2C * bz)` and is singular at zero field. The public propagation-model
selector and the production-unused Quadratic and Optimized implementations are
deleted. Both fitted-state and paired-linearization-reference propagation are
transactional and retain their existing covariance sanitization and failure
classification.

`StateFamily` remains the necessary barrel/forward state-representation tag,
but it and `stateFamilyOf(SurfaceKind)` now belong to
`SurfaceKinematicState.h`; the standalone `StateFamily.h` is deleted.

`TrackingOperationAdapter` is deleted. `Tracker` stores a typed
`SeedRefitFunction`, selected once at the ITS or MFT workflow edge. The detail
boundary is named `DetectorRefitSupport.h` and exposes direct callback-shaped
ITS and MFT refit functions; the detector-ID refit dispatcher, virtual object,
heap allocations, and adapter-named beam dispatcher are gone. Distinct
`fitTrackSeedLegs` and `refitTrackFwd` behavior is unchanged.

The enforced inventory is 33 public and 12 detail headers. No graph,
surface-plan, hole, scheduling, ordering, cut, formula, default, writer, or
workflow-ownership behavior changed.

## Validation

The complete pinned build covered common tracking, ITS tracking/interface and
workflow, MFT workflow, combined workflow, CA writer, Framework analysis/CCDB
support, and every directory contributing an `itsmft`-labelled test.

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 93
No Not Run entries.
```

The whole series passed `git diff --check` and pinned-environment
`git clang-format --diff`. The strict preflight passed with a valid AliEn
token, `ALICEO2_CCDB_NOTOKENCHECK` unset, the pinned package resolved exactly,
and every reported Geant4 data set present. The configuration reports
`CMAKE_CUDA_COMPILER=NOTFOUND` and no HIP compiler, so no real accelerator
build is claimed.

Fresh feature and frozen-base outputs are under
`/private/tmp/o2-itsmft-fine-comb-gate-20260807`. The fixture checksum manifest
passed 43/43 before and after all six replays.

| Leg | Tracks | Content hash |
| --- | ---: | --- |
| ITS feature standalone and combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS frozen standalone and combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT feature standalone and combined | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT frozen standalone and combined | 68 | `8106b08571ca593c6b76ff72b761a680` |

All eight extractor JSON files are byte-identical across the relevant
feature/frozen and standalone/combined dimensions. Eight persisted-product
comparisons report field-by-field equality for initialized writer output,
ROFs, cluster indices, labels, references, `CommonTrack` content, and sidecars.
Each MFT comparison covers 2,992 projected values with maximum absolute and
relative delta zero. Only the established undefined
`MFTTrack.mInvQPtSeed` byte remains excluded.

The named user stash remains the unchanged object
`6a90bbcd7e187673a7eeaedc2f8df07c471c09b4`.

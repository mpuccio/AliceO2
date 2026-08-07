# L10 operation-tail validation

Status: implementation complete; final verdict: **INSUFFICIENT EVIDENCE**.
The combined L9+L10 replay gate remains open because strict CCDB
authentication was unavailable.

Branch: `codex/itsmft-postm7-l9-retire-interface`
L9+L10 base: `3d15feac7a`
Package: `daily-20260717-0700-local1`
Worktree: `/Users/mpuccio/alice/run3/O2-worktrees/m6e3-commontrack-output-retirement`
Build: `/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement`

Implementation commits: `8e8c89dc63` (production) and `d535cf04f1`
(tests/guards).

## Ownership result

The operation seam now contains only the call-scoped `refitSeed()` operation.
Publication completion and adapter reset are called directly by the standalone
ITS, standalone MFT, and combined workflow specifications. `Tracker::run()`
resets only the generic `TimeFrame` event state on failure. Per-iteration
accepted-result counts are returned as generic result metadata so workflow
publication can reproduce the previous cumulative sidecar staging order without
placing publication lifecycle back in the core seam.

The following compatibility paths were deleted after whole-repository caller
audits found no live production owner:

| Removed path | Replacement / retained owner |
|---|---|
| `TrackingFrameInfoAdapters.h` and `loadClusterTrackingFrameInfo()` | normalized `SurfaceMeasurement` decoding; live `TrackingFrameInfo` scratch storage remains only where the CA compatibility representation is still consumed |
| `NativeCylinderCylinderRefitDriver.h` | `NativeRefitDriver.h`, the live descriptor-driven native refit implementation |
| dedicated old-driver comparison test and characterization harness | durable P1 comparison artifacts and the live native refit/operation tests |
| common `TimeFrame` device-propagator no-op hook and inert member | no replacement; whole-repository audit found no common caller or override |

Detector publication adapters remain live and workflow-owned. The generic core
does not contain `completeAccepted()`, `resetAdapterState()`, typed output
vectors, detector timing tables, publication validity, or the deleted public
compatibility names.

## Guard and build evidence

The focused L10 guard is
`testL10OperationTailGuard.cxx`. It verifies the refit-only seam, absence of
the deleted headers/symbols in common include/src/test, and direct publication
calls at all three workflow edges. The affected durable-build targets compiled
successfully:

```text
O2lib-ITSMFTTracking
O2lib-ITSCAWorkflow
O2lib-MFTWorkflow
O2lib-ITSMFTCombinedCAWorkflow
O2test-itsmft-tracking-l10-operation-tail-guard
O2test-itsmft-tracking-surfacemeasurement-adapters
O2test-itsmft-tracking-timeframe-covariance-lifecycle
```

Focused tests passed after rebuilding the exact binaries:

```text
o2-test-itsmft-tracking-l10-operation-tail-guard
o2-test-itsmft-tracking-surfacemeasurement-adapters
o2-test-itsmft-tracking-timeframe-covariance-lifecycle
o2-test-itsmft-tracking-m7b-runtime-count-authority-guard
o2-test-itsmft-tracking-transitionpolicytag-containment
```

## Final L9+L10 campaign

The strict preflight was run before any replay, with no CCDB token-check
bypass:

```text
O2_PACKAGE=daily-20260717-0700-local1 \
  .agents/skills/alice-o2-environment/scripts/simulation-preflight.zsh \
  --package daily-20260717-0700-local1 --strict
```

```text
No certificate filename provided!
WARNING (status=nonzero-exit): alien-token-info exited 2.
Canonical CCDB/Grid reads (e.g. alignment objects via alice-ccdb.cern.ch)
require a currently-valid token; run alien-token-init if needed.
ALICEO2_CCDB_NOTOKENCHECK: unset (expected for canonical validation runs).
```

Because the required authentication was unavailable, no canonical replay was
started. This report therefore makes no current 43/43 checksum, track-count,
hash, combined-leg, or writer-parity claim. The historical L9 replay attempt
is recorded under
`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/l9-retire-interface-20260807/`;
it stopped at the same CCDB backend failure.

The complete rebuilt serial `itsmft` CTest gate did pass independently:

```text
ctest --test-dir /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
  -L itsmft --output-on-failure -j1
```

```text
100% tests passed, 0 tests failed out of 96
No Not Run entries.
```

An earlier run reported 18 failures because test executables were stale after
the `TimeFrame` layout change. Rebuilding the durable feature build and
rerunning the exact suite produced the 96/96 result above.

The pending canonical replay must run the standalone ITS, standalone MFT, and
combined replay protocol with no CCDB bypass.

Required acceptance values remain:

| Leg | Tracks | Hash |
|---|---:|---|
| ITS | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT | 68 | `8106b08571ca593c6b76ff72b761a680` |

The pending final result must also prove 43/43 fixture checksums before and after,
byte-identical combined legs, and initialized writer parity against the L8
expected output while excluding only the known undefined
`MFTTrack.mInvQPtSeed` byte artifact. If strict CCDB authentication is
unavailable, the exact preflight/replay failure must be recorded and the result
must remain `INSUFFICIENT EVIDENCE`.

## Post-L9 review checkpoint

The target composition is now expressible as:

```text
workflow owns compatibility
Loader acts on TimeFrame
Tracker acts on TimeFrame through TrackerTraits
workflow publishes generic results
```

No lifecycle owner, coordinator, interface, or publication callback framework
was introduced to achieve this composition.

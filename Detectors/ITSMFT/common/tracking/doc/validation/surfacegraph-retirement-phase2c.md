# SurfaceGraph retirement Phase 2c validation

Phase 2c migrates the `TimeFrame`/`Tracker` ownership chain and workflow
composition to immutable `SurfaceLayout` configuration plus the per-pass
`TraversalWorkspace`. No production execution path uses `SurfaceGraph`; the
legacy graph remains only for unconverted test/support code pending Phase 2d.

## Build and gate

Validation used the pinned package `daily-20260717-0700-local1` and build
directory
`/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`.
The existing build cache was verified to point at this scratch worktree, and a
full Ninja rebuild completed successfully after refreshing stale test
consumers of the lean `Edge`, `CellPathId`, and
`TrackerIterationConfiguration::layout` APIs. The normal Ninja output was
captured in `/private/tmp/triplet-tracking-rnd-scratch-ninja-rebuild-after-onboarding-layout-fix.log`.

The requested detached full serial CTest command was run exactly as follows:

```text
ctest --test-dir /Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch \
  -L itsmft --output-on-failure -j1
```

It produced a PID file but the process disappeared without creating either a
CTest log or exit-status file. This is a runner/process-supervision
termination, not a test result; it was not interpreted as pass or fail. The
labelled list obtained with `ctest -N -L itsmft` contained exactly 90 tests.
The complete gate was therefore established with 18 serial CTest index shards,
`-I 1,5,1` through `-I 86,90,1`, each with `-j1`. Every shard returned exit
status 0 and reported `100% tests passed, 0 tests failed out of 5`; aggregate
coverage was 90 passed, 0 failed, 0 Not Run, with the ranges forming the exact
disjoint partition 1--90. Final shard logs and status files are retained as
`/private/tmp/triplet-tracking-rnd-scratch-ctest-final-shard-{01..18}.{log,status}`.

The strict preflight passed with:

```text
.agents/skills/alice-o2-environment/scripts/simulation-preflight.zsh \
  --package daily-20260717-0700-local1 --strict
```

The read-only fixture is
`/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717`.
`shasum -a 256 -c checksums.sha256` passed **43/43** both before and after
the replay campaign; the closure logs are
`/private/tmp/triplet-tracking-rnd-scratch-checksums-{before,after}.log`.

The closure also ran `git diff --check` and the pinned
`git clang-format --diff HEAD`. The latter required permission to create its
temporary linked-worktree index, then passed without a diff.

The lean ABI cleanup left several test-only consumers with retired aggregate
initializers, `CellTopologyId` path lookups, or the old textual workspace
spelling. Those consumers were updated to the current lean API; no production
topology storage or test registration was weakened.

## Replay and parity

The four canonical products were replayed with run `303000`, condition
timestamp `1784207296000`, static diamond `(0,0,0)`, `pvRes=0.05`, and
`MFTCATrackerParam.nThreads=1`:

| Product | Tracks | Content hash |
| --- | ---: | --- |
| ITS standalone | 189 | `6343211326990c75370a76b06aad5840` |
| MFT standalone | 66 | `32555b198d9b094f3f3600ec619cd2e2` |
| ITS combined | 189 | `6343211326990c75370a76b06aad5840` |
| MFT combined | 94 | `96f4c632b7e0111501a63660774480ef` |

Fresh closure outputs, logs, metrics, and comparisons are retained under
`/private/tmp/itsmft-lean-topology-closure-20260815/`:

- `its-standalone-run2/` contains `o2trac_its_ca.root`, replay logs, and
  explicitly extracted `metrics.json`;
- `mft-standalone-run2/` contains `mfttracks.root`, replay logs, and
  explicitly extracted `metrics.json`;
- `combined-run3/` contains both output ROOT files, replay logs, and separate
  ITS/MFT metric JSON files;
- `metrics-extraction.log` and `field-comparisons.log` retain the explicit
  ROOT macro evidence.

The retained corresponding parent products are under
`/private/tmp/itsmft-target-z-interval-replay/` in `its/`, `mft/`, and
`combined/`. Each fresh leg was compared with its matching retained parent
using the checked-in
`test/gate4-b42b-validation/compare_common_ca_output.C` macro. ITS products
matched field-by-field. MFT comparisons covered **2,904** standalone and
**4,136** combined projected forward-state/covariance values; both reported
zero maximum absolute and relative deltas. The standalone/combined MFT
population difference is retained and is not treated as a parity failure.

The standalone replay wrappers' final validation only performs `ls` on the
output file. Therefore the checked-in metric macros
`test/gate3-slice3-its-ca-validation/extract_metrics_its_common_ca.C` and
`test/gate0-baseline/extract_metrics_common_ca.C`, together with the checked-in
field comparator above, were run explicitly; their logs and metrics are
retained beside the replay outputs.

No CUDA or HIP validation is claimed; `nvcc` and `hipcc` are unavailable in
the pinned environment.

# SurfaceGraph retirement Phase 2c validation

Phase 2c migrates the `TimeFrame`/`Tracker` ownership chain and workflow
composition to immutable `SurfaceLayout` configuration plus the per-pass
`TraversalWorkspace`. No production execution path uses `SurfaceGraph`; the
legacy graph remains only for unconverted test/support code pending Phase 2d.

## Build and gate

Validation used the pinned package `daily-20260717-0700-local1` and build
directory
`/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`.
The full serial gate passed **90/90**:

```text
ctest --test-dir /Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch \
  -L itsmft --output-on-failure -j1
```

The strict preflight passed with:

```text
.agents/skills/alice-o2-environment/scripts/simulation-preflight.zsh \
  --package daily-20260717-0700-local1 --strict
```

The read-only fixture is
`/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717`.
`shasum -a 256 -c checksums.sha256` passed **43/43** both before and after
the replay campaign.

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

Fresh outputs and logs are retained under
`/private/tmp/itsmft-workspace-plan-20260814/`:

- `its/o2trac_its_ca.root`, `its/its_ca_replay.log`, and
  `its/its_ca_replay.time.log`;
- `mft/mfttracks.root`, `mft/mft_ca_replay.log`, and
  `mft/mft_ca_replay.time.log`;
- `combined/o2trac_its_ca.root`, `combined/mfttracks.root`, and
  `combined/combined.log`.

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

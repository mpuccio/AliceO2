# P1 direct atomic TimeFrame-source loading

Status: historical P1 validation. The later TimeFrame ownership cleanup
superseded the transaction contract below: loading now resets and fills the
configured TimeFrame directly, and a failure leaves it empty.

At P1, `loadTimeFrameSources()` was the single public atomic loading operation
and published staged data through `TimeFrame::commitLoadedEvent()`.

The former `MultiSourceTimeFrameLoader` was a one-method façade around this
transaction. It was removed without changing source qualification, decoding,
ROF/timing ownership, rollback, or workflow behavior.

Validation used the reusable pinned build at
`/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`:

- tracking and ITS, MFT, and combined workflow consumers built successfully;
- the serial `ctest -L itsmft --output-on-failure -j1` gate passed 88/88;
- whitespace and clang-format checks passed.

## Replay acceptance

The canonical replay campaign used pinned package
`daily-20260717-0700-local1`, the reusable scratch build, and fixture
`pp-20ev-run303000-seed20260716-daily20260717`. The fixture checksum manifest
passed 43/43 before and after every replay. Artifacts are retained in
`/private/tmp/itsmft-p1-direct-atomic-loader-20260813/`.

| Product | Tracks | Content hash |
| --- | ---: | --- |
| ITS standalone | 189 | `6343211326990c75370a76b06aad5840` |
| MFT standalone | 66 | `32555b198d9b094f3f3600ec619cd2e2` |
| ITS combined | 189 | `6343211326990c75370a76b06aad5840` |
| MFT combined | 94 | `96f4c632b7e0111501a63660774480ef` |

ITS standalone and combined products compare field-for-field. The established
ROOT comparator also found each fresh product field-identical to the retained
pre-cleanup candidate product. MFT comparisons covered 2,904 standalone and
4,136 combined projected forward-state/covariance values, with zero absolute
and relative delta. The standalone/combined MFT population difference is the
existing common-policy composition and was not changed by this API cleanup.

# P1 direct atomic TimeFrame-source loading

`loadTimeFrameSources()` is the single public atomic loading operation. It
stages normalized measurements and the frame workspace, then commits both
through `TimeFrame::commitLoadedEvent()` only after all validation succeeds.
`loadSources()` remains the lower-level decode/staging primitive; it is not a
second public commit route.

The former `MultiSourceTimeFrameLoader` was a one-method façade around this
transaction. It was removed without changing source qualification, decoding,
ROF/timing ownership, rollback, or workflow behavior.

Validation used the reusable pinned build at
`/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`:

- tracking and ITS, MFT, and combined workflow consumers built successfully;
- the serial `ctest -L itsmft --output-on-failure -j1` gate passed 88/88;
- whitespace and clang-format checks passed.

No replay was required: the existing atomic transaction body and all workflow
inputs are unchanged; this slice only removes the forwarding class boundary.

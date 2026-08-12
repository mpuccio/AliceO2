# TimeFrame measurement ownership

`TimeFrame` owns per-surface event measurements, source ROF intervals, and
non-owning MC-label lookup pointers. `MeasurementView` is the POD kernel view.
`MultiSourceTimeFrameLoader` decodes into a non-movable staged `TimeFrame` and
swaps only its measurement payload after workspace backfill succeeds. Source
identity is the dense source index; source ROF ranges are addressed through
the offsets held by the measurement view. Detector identity remains input
validation provenance; tracking dispatch uses surface kind.

The reusable macro-off build completed after the follow-up correction. Its
serial ITS/MFT CTest gate passed 87/87 tests (split into 44 and 43 tests only
to stay below the command-runner window). Diff and formatting checks passed.
Replay was not repeated after this ownership-only correction.

# TimeFrame measurement ownership

`TimeFrame` owns per-surface event measurements, source ROF intervals, and
non-owning MC-label lookup pointers. `MeasurementView` is the POD kernel view.
`MultiSourceTimeFrameLoader` decodes into a staged `TimeFrame` and commits it
only after workspace backfill succeeds. Detector identity remains input
validation provenance; tracking dispatch uses surface kind.

The full build and replay gates are recorded with this slice after execution.

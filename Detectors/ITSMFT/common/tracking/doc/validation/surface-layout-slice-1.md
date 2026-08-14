# SurfaceLayout retirement slice 1

`SurfaceLayout` is an immutable, owning configuration contract for the
surface catalog, ordered surface IDs, component offsets, and static hole and
seeding policies. It deliberately has no execution topology. `SurfaceGraph`
and `SurfaceGraphBuilder` remain the unchanged temporary execution backend.

The combined ITS/MFT convention is expressed by component offset `7` for the
17-position catalog, making the former missing `(6,7)` connection explicit.

Validation for this slice builds the tracking library and the focused contract
test. No replay is required because no execution path is changed.

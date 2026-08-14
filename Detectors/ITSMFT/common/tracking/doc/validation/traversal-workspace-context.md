# Traversal workspace context migration

`TimeFrame` now owns an iteration-addressable `TraversalWorkspace` through
`SurfaceTrackingScratch`. Each entry owns the selected iteration's derived
kernel parameters, material and measurement views, reference coordinates,
accepted candidates, and validity marker. `Tracker` validates the selected
configuration, resets and prepares that entry transactionally, then builds a
call-scoped non-owning `TraversalWorkspaceView` for `TrackerTraits`.

`TrackerTraits` retains only its task arena and receives every traversal input
explicitly. It no longer adopts or caches a frame, scratch, binding,
parameters, magnetic field, allocator, derived spans, or accepted candidates.
Adapters consume each successful iteration's candidates from the corresponding
frame workspace entry. Any setup or traversal failure resets the frame, which
clears every entry before an adapter can publish output.

## Validation

- Full reusable pinned build completed, including ITS, MFT, and combined CA
  workflows.
- Serial `ctest -L itsmft --output-on-failure -j1`: 88/88 passed.
- Focused lifecycle tests cover iteration-specific workspace selection,
  real MFT-road refit exceptions, and clearing all workspace candidates after
  a later iteration fails.
- `git diff --check` and pinned `git clang-format --diff` passed.
- Strict `daily-20260717-0700-local1` preflight passed. The fixture checksum
  manifest passed 43/43 before and after replay.

| Product | Tracks | Content hash |
| --- | ---: | --- |
| ITS standalone | 189 | `6343211326990c75370a76b06aad5840` |
| MFT standalone | 66 | `32555b198d9b094f3f3600ec619cd2e2` |
| ITS combined | 189 | `6343211326990c75370a76b06aad5840` |
| MFT combined | 94 | `96f4c632b7e0111501a63660774480ef` |

Fresh outputs are field-identical to the retained candidate baseline. MFT
comparisons covered 2,904 standalone and 4,136 combined projected values;
both had zero absolute and relative delta. Artifacts are retained in
`/private/tmp/itsmft-traversal-workspace-20260814/`.

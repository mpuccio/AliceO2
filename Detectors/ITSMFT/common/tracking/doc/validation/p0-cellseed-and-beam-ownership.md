# P0 CellSeed and beam-policy ownership

`CellSeed` was the sole production specialization of `SeedMetadataBase<N>`.
The template was removed and its fixed three-cluster metadata is now owned
directly by `CellSeed`, preserving member order and its GPU value semantics.
`TrackSeed` remains a separate fixed-capacity device value.

The generic refit helper remains unchanged. The ITS and MFT beam helpers are
application policy, but moving them verbatim into the two workflow specs would
duplicate their common resolution/systematic extraction. Introducing a shared
workflow helper merely to retain that duplication would add an adapter boundary
without reducing it. This slice therefore deliberately does not relocate beam
policy; that ownership change needs a separately approved workflow-level design.

Validation uses the reusable pinned build at
`/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`.
Focused `testCellRepresentation` and `testTrackSeed` passed. The full serial
ITS/MFT suite, whitespace check, and pinned `git clang-format --diff` are
required before acceptance. No replay is required because this completed
portion is layout and ownership neutral.

# Phase 2d topology deletion validation

Phase 2d removes the obsolete owning topology model and migrates remaining
clients to `SurfaceLayout`, `TraversalTopologyView`, `TraversalWorkspace`,
endpoint-only `Edge`, and exact-two-ID `CellPath`.

## Source and build closure

- The required obsolete-name scan over ITSMFT production and test source is
  empty.
- `CMAKE_HOME_DIRECTORY` remains the scratch worktree and the existing build
  directory was reused.
- Full incremental Ninja output is in
  `/private/tmp/triplet-tracking-rnd-scratch-phase2d-build-pass5.log` and
  completed successfully.

## CTest supervision diagnosis

The durable full serial command was started with stdout/stderr redirected to
`/private/tmp/triplet-tracking-rnd-scratch-phase2d-ctest-full.log`, with PID
and exit-status files. The PID disappeared without an exit-status file. This
was classified as runner/process supervision termination, not a test result.

`ctest -N -L itsmft` registered 86 current tests after deletion of the four
obsolete graph-model test registrations. The complete labelled suite was then
run in 18 serial CTest index shards, ordinal ranges `1-5`, `6-10`, ..., `81-85`,
`86-86`. Logs are
`/private/tmp/triplet-tracking-rnd-scratch-phase2d-ctest-final-shard-01.log`
through `-18.log`; matching `.status` files record each range and exit code.
Every shard exited 0, reported 100% passed, zero failed, and zero Not Run.
The manifest has 86 entries and the ranges cover each ordinal exactly once.

## Replay closure

The four canonical replays use the pinned `daily-20260717-0700-local1`
environment and fresh external output directories. Metrics and retained-parent
field comparisons are recorded with the replay artifacts. Expected products
are unchanged: ITS standalone and combined each produce 189 tracks with hash
`6343211326990c75370a76b06aad5840`; MFT standalone produces 66 tracks with
hash `32555b198d9b094f3f3600ec619cd2e2`; MFT combined produces 94 tracks with
hash `96f4c632b7e0111501a63660774480ef`. Standalone and combined MFT are
compared only with their corresponding retained parents.

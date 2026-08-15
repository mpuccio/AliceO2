# Workspace-owned traversal topology

Phase 2b stores each pass's selected traversal topology in its
`TraversalWorkspace`. The kernel-facing `TraversalTopologyView` is
non-owning: its spans and pointers borrow the workspace's edge, path,
adjacency, schedule, and road-start storage for the lifetime of the pass.
Reset clears that storage and invalidates the workspace; failed selected-plan
preparation is therefore transactional and cannot expose partial topology.

`SurfaceGraph` remains temporarily at the configuration boundary for legacy
support. TrackerTraits and scratch execution consume the workspace view and
do not consume `SurfaceGraph` or `SurfaceGraphView` directly.

## Focused evidence

Using the fixed build directory
`/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`
and pinned package `daily-20260717-0700-local1`:

- The focused `stage/tests/o2-test-itsmft-tracking-traversal-workspace-plan`
  target built successfully.
- `ctest -R 'testTraversalWorkspacePlan' --output-on-failure -j1` passed
  **1/1**, covering all five workspace ownership, reset, transactional, and
  source-guard cases.
- Pinned `git clang-format --diff HEAD^ HEAD` and `git diff --check HEAD^ HEAD`
  passed.

Full serial CTest, replay/parity validation, and GPU validation are deferred
to the later Tracker/workflow migration and final retirement phases. CUDA/HIP
validation is not claimed.

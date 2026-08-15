# Pure traversal-topology derivation

Phase 2a introduces a pure, transactional `deriveTraversalTopology()` builder
from immutable `SurfaceLayout` configuration plus a pass-local disabled-surface
mask. The result owns the selected edges, exact two-ID `CellPath` records,
reverse path adjacency, static-order mappings, schedules, and road starts.
Disabled endpoints are omitted, admitted hole-policy bridges remain selectable,
and component boundaries are never crossed. Invalid input returns no topology.

The production slice is commit `186fd9a7e4`; its independent tests and header
guards are commit `f366f545aa`. The final formatting-only correction is
`54dfe3286a`.

## Focused evidence

Using the fixed build directory
`/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`
and pinned package `daily-20260717-0700-local1`:

```text
ctest --test-dir /Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch \
  -R 'testTraversalTopology|testHeaderGraphCleanupInventory' \
  --output-on-failure -j1
```

Result: **2/2 tests passed**. Coverage includes component isolation, full
A-B-C-D derivation, disabled-middle bridge selection, disabled endpoint
filtering, transactional failure, exact `CellPath` layout, and the direct
header inventory/immutability guards.

The staged production target
`stage/lib/libO2ITSMFTTracking.dylib` built successfully in the same build
directory. The pinned `git clang-format --diff 713e4085f1 HEAD` and
`git diff --check 713e4085f1 HEAD` checks both passed.

This record does not claim full serial CTest, replay/parity, or GPU validation;
those gates remain for the later migration and retirement phases.

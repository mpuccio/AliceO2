# Campaign D: graph pair-list prototype validation

This document records test-only evidence for the non-production pair-list
prototype. Validation used the pinned O2 package
`daily-20260717-0700-local1`.

## Exact parity cases

`testPairListGraphPrototype` compares vectors by element count and raw bytes
for every trivially-copyable record or ID. It compares transitions (including
their skipped-surface witnesses), transition IDs, cells, cell IDs, CSR offsets
and entries, scheduled cells, and road-start cells. `SurfaceTransition` has two
unspecified tail-padding bytes; its value fields are compared explicitly and
the records are zero-normalized before the raw comparison. All other cache and
ID arrays are compared raw. Per-component sequences are compared directly with
`SurfacePlanBinding`.

| Case | Components | E | C | Pair declarations | Result |
| --- | --- | ---: | ---: | ---: | --- |
| ITS-like, no holes | 7 cylinders | 6 | 5 | 6 × 4 = 24 B | pass |
| MFT-like, no holes | 10 disks | 9 | 8 | 9 × 4 = 36 B | pass |
| MFT-like, one hole | 10 disks | 10 | 10 | 9 × 4 = 36 B | pass |
| ITS + MFT, no holes | 7 cylinders + 10 disks | 15 | 13 | 15 × 4 = 60 B | pass |
| ITS-like, one hole | 7 cylinders | 7 | 7 | 6 × 4 = 24 B | pass |
| ITS-like, two allowed holes | 7 cylinders | 8 | 10 | 6 × 4 = 24 B | pass |
| non-monotonic IDs | 3, 1, 4, 0 | 3 | 2 | 3 × 4 = 12 B | pass |
| empty topology | no components | 0 | 0 | 0 B | pass |

The current topology byte count is, excluding vector capacity,

`E*sizeof(SurfaceTransition) + C*sizeof(SurfaceCellTopology) + (E+1)*sizeof(uint32_t) + C*sizeof(CellTopologyId)`.

With the current sizes (`SurfaceTransition` 12 B, `SurfaceCellTopology` 8 B,
`uint32_t` 4 B, `CellTopologyId` 2 B), the representative totals are:

| Case | Current topology bytes | Pair declaration bytes |
| --- | ---: | ---: |
| ITS-like, no holes (E=6, C=5) | 150 B | 24 B |
| MFT-like, no holes (E=9, C=8) | 228 B | 36 B |
| MFT-like, one hole (E=10, C=10) | 264 B | 36 B |
| ITS + MFT, no holes (E=15, C=13) | 374 B | 60 B |
| ITS-like, one hole (E=7, C=7) | 186 B | 24 B |
| ITS-like, two allowed holes (E=8, C=10) | 232 B | 24 B |

Pair declaration bytes use the actual prototype record
`PairListBasePair { uint16_t fromIndex; uint16_t toIndex; }` (4 B), verified by
`static_assert(sizeof(PairListBasePair) == 4)`. These figures exclude vector
capacity and all object/container overhead.

The independence test holds active surfaces/base pairs/hole policy fixed while
varying only the global seeding mask (only road starts change), holds active
surfaces/base pairs/seeding fixed while changing only hole policy (expanded
topology changes with current-builder parity), and removes one immediate base
edge to prove hole expansion does not cross missing adjacency. Active selection
is varied separately; global hole and seeding masks are validated against the
union of active components.

## Scope and limitations

The helper is host/test code under `tracking/test/`; no production header,
source, graph, binding, timeframe, tracker, or tracking behavior was changed.
The four authorities remain independent: active ordered surfaces, immediate
base pairs, hole budget/mask, and seeding mask. Derived runtime caches remain
required. This is evidence only and makes no production representation
recommendation.

Production authority, mixed cylinder/disk edges, and dynamic hole traversal are
explicitly deferred. The prototype validates only immediate forward adjacency
and expands a hole span when every immediate edge exists and every skipped
surface is allowed within the independent hole policy.

# Flat pair-list graph and one combined plan

- Status: implemented and validated
- Date: 2026-08-08
- Predecessor: [Campaign D pair-list prototype](../validation/campaign-d-pair-list-prototype.md)
- Validation: [production pair-list and combined-plan validation](../validation/production-pair-list-combined-plan.md)
- Pinned O2 package: `daily-20260717-0700-local1`

## Decision

`SurfaceGraphBuilder` now accepts one `SurfaceGraphDefinition`: an ordered
surface list, a flat list of immediate forward adjacency pairs, a separate
hole policy, and a separate seeding mask. Runtime transitions, skipped-surface
witnesses, cells, CSR, schedules, and road starts remain derived data.

An omitted immediate pair is a component boundary. The combined ITS+MFT plan
therefore uses one global order, `ITS[0,7) + MFT[7,17)`, and 15 immediate
pairs; the absent pair `(6,7)` keeps the detector families disconnected. No
subgraph API or source-selected plan remains.

The combined workflow owns one `Tracker`, `TrackerTraits`, `TimeFrame`, graph,
binding, workspace, schedule, and call to `Tracker::run`.
`TrackerTraits` may partition work by `SurfaceKind` internally to retain the
existing branch-free cylinder and disk kernels, but those partitions are not
plans and carry no source identity. An iteration may bind host parameters by
`SurfaceKind`; both entries describe the same global surface arrays and plan,
while preserving the established family-specific cuts and LUT configuration.

Input sources remain distinct only as data provenance. Source IDs continue to
qualify decoded measurements, errors, labels, refit checks, and publication;
they never select a graph, binding, tracker, workspace, or family parameters.

## Combined defaults

Every surface-indexed field is the ITS value sequence followed by the MFT value
sequence. Cylinder leaves retain ITS Sync scalar semantics and disk leaves
retain MFT Sync scalar semantics. `StartLayerMask` and the graph seeding mask
cover all 17 surfaces by default (`0x1ffff`), so both disconnected components
are eligible to produce roads and tracks. Detector-local hole masks are
projected into the global rank space; hole expansion never crosses the missing
ITS-to-MFT adjacency.

## Preserved boundaries

- Standalone ITS and MFT still construct a single chain through the same flat
  definition API.
- Cylinder-to-disk transitions remain invalid; mixed kinds are valid only as
  disconnected components in one graph.
- Propagation, cuts, hole semantics, candidate ordering, refit arithmetic, and
  output formats are unchanged. Timestamp uncertainty is clamped against only
  the surfaces hit by a track, so another disconnected component's ROF length
  cannot affect it.
- The one global workspace retains per-surface ROF views and detector-local
  layer metadata so normalized input provenance is not lost.

## Acceptance

Acceptance requires the production pair-list parity tests, explicit one-plan
API guards, nonzero combined ITS and MFT tracks, the complete serial `itsmft`
suite, strict authentication preflight, and unchanged fixture checksums.
Standalone replay and both combined components must retain the frozen baseline.
The all-surface start mask does not authorize one family's scalar cuts to leak
into the other disconnected component. A missing external prerequisite is
reported as a gate limitation rather than bypassed.

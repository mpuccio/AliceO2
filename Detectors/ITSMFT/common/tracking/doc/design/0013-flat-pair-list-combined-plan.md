# Flat pair-list graph and one combined plan

- Status: implemented; combined MFT output is a physics candidate
- Date: 2026-08-09
- Predecessor: [Campaign D pair-list prototype](../validation/campaign-d-pair-list-prototype.md)
- Validation: [production pair-list and combined-plan validation](../validation/production-pair-list-combined-plan.md)
- Pinned O2 package: `daily-20260717-0700-local1`

## Decision

`SurfaceGraphBuilder` accepts one `SurfaceGraphDefinition`: an ordered surface
list, a flat list of immediate forward adjacency pairs, a separate hole policy,
and a separate seeding mask. Runtime transitions, skipped-surface witnesses,
cells, CSR, schedules, and road starts remain derived data.

An omitted immediate pair is a component boundary. The combined ITS+MFT plan
uses the global order `ITS[0,7) + MFT[7,17)` and 15 immediate pairs. The absent
pair `(6,7)` keeps the components disconnected without a subgraph or
source-selected plan API.

The combined workflow owns one `Tracker`, `TrackerTraits`, `TimeFrame`, graph,
binding, workspace, schedule, and call to `Tracker::run`. It binds exactly one
`TrackingParameters` record per iteration. For the present ITS+MFT application,
that record takes every scalar tracking-selection value from the ITS
configuration and applies it to both disconnected components. Its
surface-indexed material, resolution, geometry, and timing arrays are the ITS
values followed by the MFT values in global surface order.

`SurfaceKind` selects the kinematic state and cylinder/disk leaf operations. It
is not a tuning-policy key and does not imply detector identity. Source IDs
remain data provenance for loading and publication; they never select tracking
parameters, graph topology, operations, or scheduling.

## Combined defaults

`StartLayerMask` and the graph seeding mask cover all 17 surfaces by default
(`0x1ffff`). Detector-local hole masks are projected into global rank space;
hole expansion cannot cross the missing ITS-to-MFT adjacency.

The single selection record deliberately changes the effective MFT selection
relative to standalone MFT. The resulting combined MFT population is therefore
a candidate for a later unified physics sign-off, not an accepted replacement
baseline. This decision does not introduce detector- or family-keyed tuning in
the generic core.

## Preserved boundaries

- Standalone ITS and MFT construct their chains through the same flat graph API.
- Cylinder-to-disk transitions remain invalid; the two kinds coexist only as
  disconnected components in the current combined graph.
- Cylinder and disk operation leaves are bound from the actual surface
  descriptors. Their propagation and refit arithmetic remain in the existing
  family-specific operation boundaries and `Propagator`.
- The shared workspace retains an owner ROF view and local layer for each
  surface. Timestamp clamping therefore uses the timing of the surfaces hit by
  a track, not the first source's timing view.
- Candidate ordering, output formats, and source publication ownership remain
  unchanged apart from the intentional effect of the unified scalar selection.

## Acceptance boundary

Structural acceptance requires the flat-graph tests, single-parameter API
guards, descriptor-selected operation tests, the 40-BC ITS / 80-BC MFT timing
regression, the complete serial `itsmft` suite, strict environment preflight,
and unchanged fixture checksums. Replay validation must characterize each
combined leg against its standalone output. Exact MFT population parity is not
required under the unified selection, but the changed population cannot become
an accepted physics baseline without a separate sign-off.

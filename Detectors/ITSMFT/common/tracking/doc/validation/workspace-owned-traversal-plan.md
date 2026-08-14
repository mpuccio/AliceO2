# Workspace-owned traversal plan

`SurfacePlanBinding` is deleted. `SurfaceGraph` remains immutable static
configuration; for each traversal pass `Tracker` derives selected active
surfaces, edges, cells, schedules, and compact slots directly into that
pass's `TraversalWorkspace`.

Disabled endpoints omit their edges. A graph-template bridge over a disabled
middle surface remains selectable when both endpoints are active; hole policy
therefore stays graph-owned rather than being reconstructed by the workspace.
Failed derivation is transactional: the workspace is reset and remains
invalid rather than retaining a partial schedule.

## Evidence

- Strict `daily-20260717-0700-local1` CCDB/Geant4 preflight passed.
- Rebuilt serial `ctest -L itsmft --output-on-failure -j1`: **88/88** passed.
- Fixture manifest: **43/43** before and after replay.
- Standalone ITS: 189 tracks, `6343211326990c75370a76b06aad5840`.
- Standalone MFT: 66 tracks, `32555b198d9b094f3f3600ec619cd2e2`.
- Combined: ITS 189; MFT 94, `96f4c632b7e0111501a63660774480ef`.

These are the current scratch candidate fences after the search-window and
phi--r migrations, not the older 184-track ITS fence. Artifacts are in
`/private/tmp/itsmft-workspace-plan-20260814/`.

The current run establishes population/hash parity. Field-level comparison
against the retained candidate products remains required before integration;
this record intentionally makes no unrun comparator claim.

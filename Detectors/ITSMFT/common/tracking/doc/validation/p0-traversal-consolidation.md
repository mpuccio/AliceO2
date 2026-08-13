# P0 traversal consolidation

`TrackerTraits` now reads its committed radii, material, kernel parameters,
and road-start schedule directly.  `LayerGeometryConfigView` and the two
road forwarding methods were deleted; they only repackaged state already
owned by `TrackerTraits`.

`TrackerTraversalPreparation` remains deliberately separate.  It contains
the frozen MFT reference-coordinate scattering policy and its independently
tested arithmetic.  Moving or generalizing that policy is outside this
structural cleanup.

Validation on the reusable pinned build:

- affected library and focused tests built;
- serial `ctest -L itsmft --output-on-failure -j1`: 88/88 passed;
- `git diff --check` and pinned `git clang-format --diff` passed.

No replay was required: this change removes forwarding and a span-only view
without changing transition preparation formulas, ordering, or selections.

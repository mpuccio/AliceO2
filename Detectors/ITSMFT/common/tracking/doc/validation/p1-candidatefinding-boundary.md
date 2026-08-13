# P1 CandidateFinding boundary

`CandidateFinding.h` now exposes only the generic candidate operations used by
`TrackerTraits`: tracklet-window projection and acceptance, cache binding, and
cell-seed construction. The cylinder and disk leaves have internal linkage in
`CandidateFinding.cxx`.

No candidate algorithm, descriptor selection, kernel/GPU boundary, graph,
timing, or output path changed. Numerical tests exercise the existing generic
operations with an explicit `SurfaceKind`; compact test-local fixture helpers
only keep their argument lists readable and export nothing.

Validation used the reusable pinned build:

- affected tracking library and ITS, MFT, and combined workflow consumers built;
- focused tracklet, cell, orchestration, and boundary tests passed;
- serial ITS/MFT CTest passed in two contiguous runs: the first 82 tests
  completed successfully before the external runner ended, then the remaining
  six ITS/MFT workflow/configuration tests passed;
- `git diff --check` and `git clang-format --diff` passed.

No replay was required: production execution remains the same generic calls;
this slice only removes external declarations and gives their source-local
coordinate leaves internal linkage.

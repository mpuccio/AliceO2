# Propagator state-operation boundary

## Scope

`Propagator` is now the only public state-evolution API for
`SurfaceKinematicState`. The former barrel and forward public operation
headers were removed. Coordinate-specific leaves remain implementation detail
in `detail/SurfaceStateOperations.h`.

Candidate finding retains three-hit seed construction. It is the sole
non-Propagator production client of the private family leaves, and uses only
their `buildSeed` operations. It delegates every subsequent measurement
attachment to `Propagator::attachMeasurement`.

The temporary dense `float[5][5]` working arrays were renamed
`DenseMatrix5`. They cannot be packed: Jacobians and the `J*C` intermediate
are not symmetric, although persisted state covariance remains packed.

`FamilyMaterialOperations.cxx` retains a private-leaf include because it
implements the two family-local covariance projections used by material
correction. It exposes no public state-evolution API and is called only by
`Propagator`.

## Validation

- Complete incremental build: 186 targets.
- Serial ITS/MFT CTest: 88/88 passed. The captured main run completed the
  first 85 tests before tool-output truncation; the final three workflow
  contract tests were rerun explicitly and passed.
- `git diff --check` and pinned `git clang-format --diff` passed.

No replay was run. This is a boundary-only migration: state-operation leaf
formulas, ordering, material inputs, graph construction, ROF handling, and
output publication are unchanged. Numerical coverage remains in the migrated
Propagator and family-leaf tests.

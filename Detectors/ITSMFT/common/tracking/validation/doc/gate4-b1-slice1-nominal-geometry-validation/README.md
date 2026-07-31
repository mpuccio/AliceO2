# Gate 4 B1 Slice 1: real-geometry validation run record

Closes the "not exercised against a real geometry file" gap noted at the end
of the Slice 1 implementation (commits `7e1730adec`, `9d0b603619`).
**Documentation only** — no source, test, CMake, table, tolerance, provider,
layout, or workflow file is changed by this commit. No `.root`
fixture/output is committed; every result lives outside the Git worktree
under the durable validation-artifact hierarchy referenced below.

## What ran

`o2-itsmft-nominal-geometry-validator` (the target added in the two commits
above), against the durable fixture's own unaligned geometry file:

- Input: `/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717/o2sim_geometry.root`
  — checksum-verified against that fixture's own `checksums.sha256` manifest
  before use; `o2sim_geometry-aligned.root` in the same directory was never
  read.
- Package: `O2_PACKAGE=daily-20260717-0700-local1` (pinned explicitly,
  matching the fixture's own `provenance.json`).
- Build: `O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/gate4-s0-sparse-hotloops-design-daily20260717`,
  staged `o2-itsmft-nominal-geometry-validator` from this worktree at commit
  `9d0b603619a166b7a49d8a6aeb4863e617dff771` (worktree clean at run time).

## Reproducible command

```bash
FIXDIR=/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717
BUILD=/Users/mpuccio/alice/run3/O2-worktree-builds/gate4-s0-sparse-hotloops-design-daily20260717

O2_PACKAGE=daily-20260717-0700-local1 O2_BUILD_DIR=$BUILD \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  o2-itsmft-nominal-geometry-validator --geometry "$FIXDIR" --detector ITS --surfaces 7 --format text

O2_PACKAGE=daily-20260717-0700-local1 O2_BUILD_DIR=$BUILD \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  o2-itsmft-nominal-geometry-validator --geometry "$FIXDIR" --detector MFT --surfaces 10 --format text
```

`--format json` run identically for both detectors alongside `--format text`.
All four invocations exited `0` (`ValidatorStatus::Ok`).

## Result

All 7 ITS and all 10 MFT derived per-surface entries (reference coordinate,
radial/z bounds), full stdout/stderr captures, a checksum/provenance record,
and explicit proof the unaligned (`preferAlignedFile=false`,
`applyMisalignment=false`) load route was taken are recorded under:

```
/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/gate4-b1-slice1-nominal-geometry-validation/pp-20ev-run303000-seed20260716-daily20260717/
```

(see that directory's own `README.md` for the full record; files there are
read-only).

## Scope note

This run reports geometry-derived values only. It establishes that the
detached validator target builds, links outside `O2::ITSMFTTracking`, and
runs end-to-end against real, checksum-verified, unaligned geometry for both
detectors. It does **not** author, compare against, or imply sign-off on any
`ITSSurfaceSpec`/`MFTSurfaceSpec` value or tolerance — that remains blocked
on detector-geometry-owner input, unchanged from the accepted Gate 4 B1
design.

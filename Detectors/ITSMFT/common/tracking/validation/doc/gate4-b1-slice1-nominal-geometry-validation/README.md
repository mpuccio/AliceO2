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

## Correction: `--format json` stdout contract (commits `2a698e9cd8`, `dff33a9341`)

The original run above used a post-hoc extraction (`grep '^{'`) to produce
clean JSON artifacts from `--format json` stdout, because
`GeometryManager::loadGeometry()`/`GeometryTGeo`'s own `FairLogger` console
sink wrote `[INFO]` lines to stdout interleaved around the tool's JSON
payload (and, for MFT, a further `[INFO] ~GeometryTGeo` line *after*
`main()` returns, during static teardown) — evidence of a bug, not an
acceptable CLI contract. Fixed by `JsonStdoutGuard`
(`validation/src/JsonStdoutGuard.{h,cxx}`, commit `2a698e9cd8`): redirects
stdout to stderr at construction, for the remainder of the process's life,
and writes the payload once through a duplicate of the true original stdout
fd. `--format text` is unchanged. Test coverage added in `dff33a9341`
(`testJsonStdoutContract.cxx`), capturing the complete stdout stream and
parsing it with RapidJSON for a synthetic success case and two typed
failures.

Re-ran the exact same command as above, with `--format json`, against the
identical checksum-verified `o2sim_geometry.root`:

```bash
FIXDIR=/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717
BUILD=/Users/mpuccio/alice/run3/O2-worktree-builds/gate4-s0-sparse-hotloops-design-daily20260717

O2_PACKAGE=daily-20260717-0700-local1 O2_BUILD_DIR=$BUILD \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  o2-itsmft-nominal-geometry-validator --geometry "$FIXDIR" --detector ITS --surfaces 7 --format json \
  1>its-report.json 2>its-report.stderr.log

O2_PACKAGE=daily-20260717-0700-local1 O2_BUILD_DIR=$BUILD \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  o2-itsmft-nominal-geometry-validator --geometry "$FIXDIR" --detector MFT --surfaces 10 --format json \
  1>mft-report.json 2>mft-report.stderr.log
```

Verified directly on the raw stdout capture (no extraction step): `wc -l`
is `1` for both; `jq -e .` and `python3 -m json.tool` both accept both files
as valid JSON; both processes exit `0`; content is byte-identical to the
original run's extracted values. Full raw captures recorded under, without
modifying anything already there:

```
/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/gate4-b1-slice1-nominal-geometry-validation/pp-20ev-run303000-seed20260716-daily20260717/json-stdout-fix-recheck/
```

## Gate 4 acceptance-cleanup C1: lossless geometry float encoding (commits `590888328b`, `f6991d434f`)

The fix above made `--format json` stdout exactly one JSON document, but
its float fields (`referenceCoordinate`, `radialMin`, `radialMax`) still
used fixed six-decimal formatting — not guaranteed to round-trip to the
exact source `float32` bit pattern (flagged explicitly in the Gate 4 B1
Slice 2 authoring report: MFT's largest magnitude, ~77.5, needs 8
significant digits at 6 decimal places, exceeding float32's ~7.22-digit
guarantee).

Replaced with `std::to_chars(float)`'s shortest round-tripping
representation, exposed as the new public `formatLosslessFloat()`
(`NominalGeometryReport.{h,cxx}`, commit `590888328b`). Confirmed available
at this build's exact deployment target (`-mmacosx-version-min=14.0`) via a
minimal compile probe using the production compile line's exact flags —
`std::to_chars(float)` works; `std::from_chars(float)` does not (Apple
libc++ marks it `availability(macos,introduced=26.0)`), which only affects
how the new test decodes tokens back (`std::strtof`), not production.
`--format text` is unchanged, still fixed six-decimal.

Test coverage added in `f6991d434f` (`testLosslessFloatRoundTrip.cxx`):
bitwise round-trip proof (original float32 bits → printed token →
`strtof`-parsed float32 bits, compared via `memcpy`'d `uint32_t`, zero
tolerance) for real geometry values from this fixture, `+0.0f`/`-0.0f`
(equal as floats, different bit patterns — RapidJSON's own `GetFloat()`
was tried first and found to silently drop `-0`'s sign bit, and to
classify some legitimately-shortest float tokens as integers rather than
floats; both real RapidJSON behaviors unrelated to what the encoder
promises, documented in the test file rather than worked around), adjacent-
to-1.0f max-precision values, extreme magnitudes/subnormals, and a
~20500-point deterministic sweep across the finite bit-pattern space.

Re-ran the same ITS/MFT commands, `--format json`, against the identical
checksum-verified `o2sim_geometry.root`. Both exited `0`; `wc -l` is `1`
for both; `jq -e .` and `python3 -m json.tool` both accept both files.
Content is numerically consistent (within 1e-3) with the pre-C1 recorded
values — confirming this changed only the float *representation*, not the
underlying geometry values. Raw captures, their own SHA-256, and the full
verification record are under, without modifying anything already there:

```
/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/gate4-b1-slice1-nominal-geometry-validation/pp-20ev-run303000-seed20260716-daily20260717/acceptance-cleanup-c1-lossless-json/
```

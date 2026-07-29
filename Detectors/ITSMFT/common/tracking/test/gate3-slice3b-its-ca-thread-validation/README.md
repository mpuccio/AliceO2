# Gate 3 Slice 3b: ITS common-CA thread/determinism validation (sibling to gate3-slice3-its-ca-validation)

Scripts and results characterizing `o2-its-ca-tracker-workflow`'s new
`ITSCommonCATrackerParam.nThreads` knob (added in commit `62e2c75bc8`,
after the accepted Gate 3 Slice 3 characterization ran) against the durable
Gate 0 fixture, per this validation task's scope: reproduce the accepted
Slice 3 canonical baseline with `nThreads=1` explicit, characterize
`nThreads=12` for semantic determinism and timing (without overclaiming a
scaling proof), obtain a labelled legacy-Sync reference (no diamond
override), rerun the MFT common-CA non-regression fence, and add a negative
`nThreads=0` check.

**This is characterization, not a bitwise-parity gate.** No physics
tolerance is invented here, and the legacy leg's real-vertexer output is
explicitly not compared against the common-CA legs' static-diamond output
as a pass/fail signal -- see `characterization_summary.md` §4.

**Production code is untouched.** This directory contains only validation
scripts, one small metrics extractor, and result documentation. No file
under `Detectors/ITSMFT/common/tracking/{src,include}`,
`Detectors/ITSMFT/ITS/workflow-ca`, `Detectors/ITSMFT/ITS/tracking`,
`gate0-baseline`, or `gate3-slice3-its-ca-validation` was modified. No
ROOT/binary fixture or replay output is committed; all replay/metrics
directories referenced below live under
`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/gate3-slice3b-its-ca-thread-validation/`,
outside the git worktree.

## Scripts (new, this directory)

- `replay_tracking_its_common_ca_nthreads.sh` -- sibling to
  `gate3-slice3-its-ca-validation/replay_tracking_its_common_ca.sh`, adding
  a required `ITS_CA_NTHREADS` env var (that script has no thread knob at
  all -- `ITSCommonCATrackerParam.nThreads` did not exist when it was
  written). Same explicit static-diamond constraint
  (`diamondPos=(0,0,0)`, `pvRes=0.05`). Produces `o2trac_its_ca.root`.
- `replay_tracking_its_legacy_sync.sh` -- runs the frozen
  `o2-its-reco-workflow` ITS-only in its **ordinary Sync/default-vertex
  mode**: only `ITSCATrackerParam.nThreads` is set, `useDiamond` stays at
  its struct default (`false`), so the real per-event Sync vertexer runs.
  This is a DIFFERENT leg from
  `gate3-slice3-its-ca-validation/replay_tracking_its_legacy_diamond.sh`
  (which forces `useDiamond=true`); neither file in that directory is
  touched. Produces `o2trac_its.root`. Explicitly closes stdin around the
  piped invocation as a defensive measure against a previously documented
  flaky hang for this same invocation shape (see that directory's README
  and this directory's `characterization_summary.md` §4) -- not observed
  to recur here, but guarded against regardless.
- `negative_check_its_common_ca_nthreads0.sh` -- runs the full piped
  `o2-its-cluster-reader-workflow | o2-its-ca-tracker-workflow` with an
  otherwise-valid static-diamond configuration plus
  `ITSCommonCATrackerParam.nThreads=0`, asserting the expected
  `LOGP(fatal)` fires and no track data is written. Unlike
  `gate3-slice3-its-ca-validation/negative_checks_its_common_ca.sh`'s two
  cases (which fire inside `ConfigPreflight`, before any device exists),
  this rejection lives inside `initialiseTracker()`, only reached on the
  tracker device's first real `run()` call -- so this check needs real
  piped data, and asserts zero TTree entries rather than file
  non-existence (the writer device creates the file regardless; see
  `characterization_summary.md` §6 for the full explanation). Wrapped in
  `timeout`/`gtimeout` as a safety bound.
- `extract_track_content_hash_its_legacy.C` -- new, minimal metrics
  extractor that adds a `trackContentHash` field (identical MD5 definition
  to `gate3-slice3-its-ca-validation/extract_metrics_its_common_ca.C`'s) for
  the legacy-Sync leg's `o2trac_its.root`, since
  `gate0-baseline/extract_metrics.C`'s `extractITS()` (reused unmodified for
  this leg's other metrics, see below) has no hash field of its own.

## Scripts reused UNCHANGED from sibling directories

- `gate3-slice3-its-ca-validation/extract_metrics_its_common_ca.C` -- for
  the common-CA legs' track-count/efficiency/fake-rate/hash metrics.
- `gate3-slice3-its-ca-validation/extract_metrics_its_legacy_diamond.sh` --
  for the legacy-Sync leg's track-count/efficiency/fake-rate metrics (this
  wrapper is generic over `FIXTURE_DIR`/`REPLAY_DIR`; it does not assume
  the replay used a diamond override, despite its name).
- `gate0-baseline/replay_tracking_common_ca.sh` +
  `gate0-baseline/extract_metrics_common_ca.C` -- for the MFT common-CA
  non-regression fence, exactly as in Gate 3 Slice 3.

None of the files in the two paragraphs above is read, written, or modified
by anything in this directory beyond an ordinary invocation (or, for the
`.C` macros, a ROOT interpreter load) with this session's own paths.

## Running

See `manifest.json` for the exact package/fixture/timestamp/output-location
values and full command lines used to produce the results in
`characterization_summary.md`. In short, all scripts are invoked through
`.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh` with
`O2_BUILD_DIR` pointed at a worktree build of this exact commit and
`O2_PACKAGE=daily-20260717-0700-local1` pinned explicitly (never an alias),
per that skill's reproducible-validation policy.

## Determinism / timing note

All 6 measured common-CA replays across both thread counts (3×
`nThreads=1`, 3× `nThreads=12`) produced byte-identical extracted-metrics
JSON: **203 tracks, `trackContentHash=ee7f7c794d60f2362fd2564258b7887e`**,
matching the accepted Gate 3 Slice 3 baseline exactly. No wall-time or
CPU-time benefit from `nThreads=12` was observed on this fixture; see
`characterization_summary.md` §3 for the overhead analysis and why this is
explicitly not treated as a scaling result.

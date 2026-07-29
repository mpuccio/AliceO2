# Gate 3 Slice 3b: ITS common-CA thread/determinism characterization + labelled legacy reference

Commit under test: `62e2c75bc8ff3b5c556910116a860bde2af4b424` ("ITSMFTTracking:
add explicit ITSCommonCATrackerParam.nThreads"), branch
`codex/itsmft-gate3-its-common-ca-thread-validation`. Package
`daily-20260717-0700-local1` (pinned explicitly, requested == resolved).
Fixture `pp-20ev-run303000-seed20260716-daily20260717`, checksums verified
before and after the full replay campaign. Full commands/paths in
`manifest.json`.

**This is validation-only characterization, not a bitwise-parity gate and
not a scaling proof.** No production tracker/workflow/config/CMake/legacy
code was touched by this task. No file under `gate0-baseline` or
`gate3-slice3-its-ca-validation` was modified -- both were only read or
invoked unchanged.

## 1. Build

`o2-its-ca-tracker-workflow`, `o2-its-cluster-reader-workflow`,
`o2-its-reco-workflow` (frozen legacy), `o2-mft-cluster-reader-workflow`,
`o2-mft-ca-tracker-workflow`, and the two dlopen'd DPL plugin libraries
(`libO2FrameworkAnalysisSupport`, `libO2FrameworkCCDBSupport`) all built
clean (1363/1363 ninja steps, exit 0) from this exact worktree/commit
against the pinned package, in a build directory dedicated to this task.

## 2. ITS common-CA canonical baseline (nThreads=1)

Explicit static diamond: `diamondPos=(0,0,0)`, `pvRes=0.05`,
`ITSCommonCATrackerParam.nThreads=1`. One warm-up + 3 measured fresh-directory
replays (2 required, 3 run for a more robust median given timing variance
observed -- see §3).

| Run | Output tracks | Track content hash | Wall time (real, s) | Peak RSS |
|---|---|---|---|---|
| run1 | 203 | `ee7f7c794d60f2362fd2564258b7887e` | 5.95 | 2.53 GiB (2649931776 B) |
| run2 | 203 | `ee7f7c794d60f2362fd2564258b7887e` | 5.72 | 2.53 GiB (2649800704 B) |
| run3 | 203 | `ee7f7c794d60f2362fd2564258b7887e` | 5.58 | 2.53 GiB (2649784320 B) |

**Median wall time 5.72 s, median peak RSS ≈2.53 GiB.** All three runs'
extracted metrics JSON are byte-identical to each other and to the accepted
Gate 3 Slice 3 baseline
(`gate3-slice3-its-ca-validation/characterization_summary.md` §3-4: **203
tracks, hash `ee7f7c794d60f2362fd2564258b7887e`**). The accepted semantic
baseline is reproduced exactly, one commit later (`nThreads` did not exist
in that baseline; explicitly pinning it to `1` here reproduces the same
result, confirming the added knob's default-equivalent value does not shift
tracking output). Raw `o2trac_its_ca.root` bytes differ run-to-run as
already established for this workflow (ROOT `TFile` UUID/timestamp
metadata, not content).

## 3. ITS common-CA thread characterization (nThreads=12)

Same static-diamond configuration, `ITSCommonCATrackerParam.nThreads=12`.
One warm-up + 3 measured fresh-directory replays.

| Run | Output tracks | Track content hash | Wall time (real, s) | User (s) | Peak RSS |
|---|---|---|---|---|---|
| run1 | 203 | `ee7f7c794d60f2362fd2564258b7887e` | 15.31 (outlier) | 6.02 | 2.53 GiB |
| run2 | 203 | `ee7f7c794d60f2362fd2564258b7887e` | 6.70 | 5.96 | 2.53 GiB |
| run3 | 203 | `ee7f7c794d60f2362fd2564258b7887e` | 6.03 | 6.06 | 2.53 GiB |

**Median wall time 6.70 s, median peak RSS ≈2.53 GiB.**

**Semantic determinism: all three nThreads=12 runs' metrics JSON are
byte-identical to each other AND to every nThreads=1 run above** (identical
203 tracks, identical `trackContentHash`). Changing the thread count did not
change tracking output at all in this fixture, at either thread count.

**This is explicitly NOT a scaling proof**, per task instruction #4:

- Median wall time did not improve at nThreads=12 (6.70 s) versus nThreads=1
  (5.72 s) -- if anything it was nominally higher, well within run-to-run
  noise given the small sample.
- `run1` at nThreads=12 took 15.31 s real vs. 6.02 s user -- a large
  real/user gap consistent with the process being descheduled or blocked
  (e.g. on CCDB I/O, disk, or another concurrent host process) rather than
  any tracker-internal effect; the other two nThreads=12 runs (6.70 s, 6.03
  s real) show real ≈ user, i.e. no such stall.
- Across ALL SIX measured common-CA replays (both thread counts), `user`
  CPU-seconds sits in a tight 5.96-6.06 s band regardless of `nThreads`.
  This is the clearest signal here: total CPU-seconds consumed does not
  grow with thread count the way it would if `computeLayerTracklets`-style
  work were meaningfully parallelized across 12 `tbb::task_arena` workers
  for this fixture. No separate tracker-only timing probe exists in this
  workflow to isolate reconstruction time from fixed per-run overhead
  (workflow startup, CCDB condition fetch, DPL device bring-up, geometry
  load) -- all of which this fixture's small size (20 events, 2304 ROFs)
  makes proportionally large. This validation does not claim the tracker
  fails to parallelize internally; it reports that, at THIS fixture's size,
  no wall-time or CPU-time benefit from nThreads=12 was observable end to
  end, and that fixed overhead plausibly dominates any parallel portion of
  the work.

## 4. Legacy ITS labelled reference (ordinary Sync, no diamond override)

Frozen `o2-its-reco-workflow`, ITS-only, in its **ordinary Sync
default-vertex mode** -- `ITSCATrackerParam.useDiamond` left at its struct
default (`false`); only `ITSCATrackerParam.nThreads=1` is set. This is a
**different, also non-equivalent** leg from
`gate3-slice3-its-ca-validation/replay_tracking_its_legacy_diamond.sh`
(which forces `useDiamond=true` to pair with a static-diamond common-CA
leg): this validation's legacy leg instead exercises the tracker's real
per-event vertexer, as it would run in production Sync reconstruction.

| Metric | Value |
|---|---|
| Input clusters | 4103 |
| Input ROFs | 2304 |
| Output tracks | 175 |
| Track content hash | `f69f50948c20fb7375cb5b7a6d126a3e` |
| Clusters/track (mean/median/min/max) | 6.22 / 7 / 4 / 7 |
| chi2 (mean/median/min/max) | 7.70 / 5.60 / 0.258 / 59.4 |
| MC reconstructable | 142 |
| Matched | 120 |
| Efficiency | 0.845 |
| Fake tracks / rate | 13 / 0.0743 |
| Clone tracks / rate | 0 / 0 |
| Wall time (real) | 6.68 s |
| Peak RSS | 639 MiB (669974528 B) |

Ran once, single attempt, no hang (see below). **This is a labelled
reference point only.** Per task instruction #5, no parity or pass/fail
comparison is drawn between this leg and the common-CA legs above: the
common-CA legs use a static synthetic diamond vertex constraint with no
per-event timing/position information; this leg uses the legacy tracker's
own real Sync vertexer (5 real vertices per the vertexer-stage log,
per-event timing-compatible tracklet finding). Track count (175 vs. 203),
efficiency (0.845 vs. 0.951), and hash necessarily differ for reasons that
have nothing to do with correctness of either leg -- they are different
reconstruction configurations by design, not a regression signal in either
direction. Note also that legacy `TrackerParamConfig<ITS>` (the
`ITSCATrackerParam` legacy namespace) has no `useDiamond` interaction with
the timestamp-envelope fix in commit `b47271fb13` at all -- that fix only
changed the common-CA `TrackerTraits::computeLayerTrackletsForPolicy`'s
diamond-timestamp derivation, which this leg never exercises (its
`useDiamond` stays `false`).

**No flaky hang observed on this leg.** `gate3-slice3-its-ca-validation`'s
README documents this same piped
`o2-its-cluster-reader-workflow | o2-its-reco-workflow` invocation shape
hanging on 2 of 4 attempts for its own (differently-configured)
legacy-diamond leg, root cause not identified there. This leg's single run
completed normally on its first attempt. `replay_tracking_its_legacy_sync.sh`
explicitly closes stdin (`< /dev/null`) around the whole piped invocation as
a defensive measure regardless -- this is NOT claimed to be a proven fix for
that hang (never reproduced with stdin open to compare), just a documented,
low-cost precaution layered on top of the existing `bash -o pipefail -c`
exit-code-masking guard.

## 5. MFT common-CA non-regression fence

`gate0-baseline/replay_tracking_common_ca.sh` +
`gate0-baseline/extract_metrics_common_ca.C`, unchanged, against the same
fixture:

| | Expected | Observed |
|---|---|---|
| Output tracks | 70 | 70 |
| Track content hash | `24737e73b7146bf3bd35a90a2517c527` | `24737e73b7146bf3bd35a90a2517c527` |

**Exact match. MFT non-regression fence held** -- this Gate 3 Slice 3b
validation work did not disturb the MFT common-CA path.

## 6. Negative check: `ITSCommonCATrackerParam.nThreads=0`

Not part of the two Slice-3-style `ConfigPreflight`-level negative checks
(those fire at `defineDataProcessing()`-time, before any device exists).
The `nThreads<=0` rejection added in commit `62e2c75bc8` lives inside
`ITSMFTTrackingInterface::initialiseTracker()`, only reached from
`CATrackerDPL::updateTimeDependentParams()` on the tracker device's FIRST
real `run()` call (confirmed by source inspection of
`Detectors/ITSMFT/ITS/workflow-ca/src/CATrackerSpec.cxx:184-189` before
writing this check) -- so this check needs the full piped invocation with
real data flowing, not a standalone no-input run, and devices (including
the writer) ARE constructed first.

Ran once, in a fresh `negative-nthreads0/` directory, with an otherwise
identical valid static-diamond configuration (`useDiamond=true,
diamondPos=(0,0,0), pvRes=0.05`) plus `ITSCommonCATrackerParam.nThreads=0`:

| | Result |
|---|---|
| Exit status | Non-zero (driver exit 1; `its-ca-tracker` child device itself exited 128/SIGABRT via `LOGP(fatal,...)`) |
| Expected fatal message | Present: `"ITS CA tracker requires ITSCommonCATrackerParam.nThreads > 0, got 0"` |
| `o2trac_its_ca.root` | **File exists** (created by the `its-ca-track-writer` device's own ordinary startup, independent of tracker success -- confirmed by direct inspection: this differs from the Slice 3 preflight-level checks, where no device is ever constructed and no file appears at all) |
| `o2trac_its_ca.root`'s `o2sim` tree entry count | **0** -- confirmed no track data was ever written |

**PASS**, with the caveat above about output-file existence recorded
explicitly rather than silently treated the same as the Slice 3 checks. No
baseline output was at risk: this ran in its own fresh directory, and the
fatal aborts the process during the very first `run()` call, upstream of
ever producing track data.

## 7. Fixture integrity

`shasum -a 256 -c checksums.sha256 --status` run against the fixture
directory both before the first replay of this session and again after the
full campaign (baseline + thread characterization + legacy reference + MFT
fence + negative check): all 40 listed files OK both times. Fixture never
regenerated, overwritten, or chmod'd; all replay/metrics output went to
fresh directories under
`O2-validation-artifacts/itsmft/gate3-slice3b-its-ca-thread-validation/`,
never inside the fixture directory.

## Known limits

1. **No isolated tracker-only timing.** As in §3: this workflow exposes no
   timing probe that separates tracking compute from workflow
   startup/CCDB/DPL/geometry-load overhead, and this fixture (20 events,
   2304 ROFs) is small enough that this overhead plausibly dominates total
   wall time at both thread counts tested. The observed lack of wall-time
   improvement at nThreads=12 is reported as-is and is NOT used to draw any
   conclusion about the tracker's internal parallel scalability on larger
   inputs.
2. **One timing outlier** (nThreads=12 run1, 15.31 s real vs. ~6 s user) is
   recorded, not discarded, with its likely explanation (external
   stall/contention, not tracker behavior) noted in §3; the median already
   in this report is robust to it (3 samples, not just 2).
3. **Legacy leg is a single run, not a repeated-determinism check** -- this
   validation's scope (per task) was to obtain one labelled reference point
   for the ordinary Sync/default-vertex configuration, not to re-establish
   determinism for the legacy tracker (already out of scope and not
   contradicted by anything observed here).
4. **No physics tolerance/threshold is defined** for either leg's own
   metrics (efficiency, fake rate, chi2 shape, etc.) -- all numbers above
   are recorded as labelled reference points for future comparison, not
   evaluated against a pass/fail bound, matching Gate 3 Slice 3's own
   stated scope.
5. All CCDB-network-dependency, no-local-snapshot, and installed
   sim/digi-stack-provenance caveats already documented for the shared
   fixture in `gate0-baseline/README.md`/`manifest.json` apply identically
   here, since the fixture is reused unmodified.

**Not marking this gate complete** -- the integration owner reviews this
evidence, per `AgentCoordination.md` §12 ("An agent may not declare a gate
passed based only on compilation").

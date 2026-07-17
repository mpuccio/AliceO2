# Gate 2 common-CA MFT characterization (sibling to the Gate 0 baseline)

Commit under test: `43011425b02780425f1929d1cdb7c7053cc6564e`
("ITSMFT: record explicit Gate 2 layout construction"), branch
`codex/itsmft-mft-ca-baseline`.

**This is not the Gate 0 baseline and does not modify it.** Gate 0's
accepted baseline (`baseline_summary.md`, `manifest.json`) exercises
`o2-mft-reco-workflow`, i.e. the legacy `O2::MFTTracking` package. Gate 2
modifies the separate common `o2::itsmft::tracking` core used by
`o2-mft-ca-tracker-workflow`. No accepted physics/performance baseline
previously covered that path. This file characterizes it, as a labelled,
reproducible reference point -- it does not claim Gate 2 completion, and it
does not require the common-CA numbers to match the legacy MFT numbers:
these are distinct algorithms (see `AgentCoordination.md` Gate 2 status).

Fixture: the same 20 unanchored pp events (`pythia8pp`/`TGeant4`),
PIPE+ITS+MFT only, run 303000, seed 20260716, reused unmodified from the
Gate 0 baseline (see `manifest.json` for full generation provenance). Only
`mftclusters.root`, `o2sim_Kine.root`, and
`o2simdigitizerworkflow_configuration.ini` are read; no ITS component runs
in this bounded replay. Full command lines, workflow-input-spec
verification, and output-contract verification are in
`manifest_common_ca.json`.

## Workflow under test

`o2-mft-cluster-reader-workflow | o2-mft-ca-tracker-workflow`, invoking
`o2::mft::CATrackerDPL` / `getCATrackerSpec`
(`Detectors/ITSMFT/MFT/workflow/src/CATrackerSpec.cxx`), which is the thin
MFT adapter around the Gate 2 common `o2::itsmft::tracking` core
(`Detectors/ITSMFT/common/tracking`). This is distinct from
`o2-mft-reco-workflow`, which uses the legacy prototype in
`Detectors/ITSMFT/MFT/tracking`.

Confirmed by direct source comparison (not assumed): `CATrackerSpec.cxx`'s
`InputSpec`s (`MFT/COMPCLUSTERS`, `MFT/PATTERNS`, `MFT/CLUSTERSROF`,
`MFT/CLUSTERSMCTR`) and `TrackWriterSpec.cxx`'s output contract
(`mfttracks.root`, tree `o2sim`, branches `MFTTrack`/`MFTTrackMCTruth`/etc.)
are byte-identical to the legacy `o2-mft-reco-workflow` path -- both share
the same `TrackWriterSpec.cxx` source, called with `useCA=true` here. This
is why `extract_metrics_common_ca.C`'s counting logic is an unmodified copy
of `extract_metrics.C`'s `extractMFT()`: the semantics genuinely match, per
task instruction #9, so no new extraction logic was needed, only a smaller
driver that does not also require an ITS replay output (the common-CA
replay directory has no ITS component to read).

`o2-mft-ca-tracker-workflow` has no `--clusters-from-upstream` or
`--cluster-rof-branch-only` flag: unlike `o2-mft-reco-workflow` it has no
embedded clusterizer or cluster-writer device, so there is no
upstream-vs-clusterize ambiguity to disambiguate and no `CLUSTERSMC2ROF`
requirement to route around.

## Metrics (canonical: single-threaded replay from fixed clusters)

| Metric | MFT common-CA | MFT legacy (Gate 0, for reference only) |
|---|---|---|
| Input clusters | 1679 | 1679 |
| Input ROFs | 2304 | 2304 |
| Output tracks | 91 | 149 |
| Clusters/track (mean / median / min / max) | 6.51 / 6 / 5 / 10 | 7.50 / 8 / 5 / 10 |
| chi2 (mean / median / min / max) | 0.2225 / 0.0835 / 0.000986 / 2.020 | 1.61 / 0.21 / 0.006 / 116.78 |
| MC reconstructable | 141 | 141 |
| Matched | 77 | 132 |
| Efficiency | 0.5461 | 0.9362 |
| Fake tracks / rate | 0 / 0 | 5 / 0.0336 |
| Clone tracks / rate | 0 / 0 | 0 / 0 |

Same fixture, same denominator convention (see below), same input cluster
count as the legacy row -- both consume `mftclusters.root` from the
identical fixture. The common-CA core currently reconstructs fewer tracks
with substantially lower efficiency and zero fakes; this is expected and
not a regression signal in either direction, since Gate 2 has not yet
reached "ITS and MFT algorithm regression report" / demonstrated parity
(`AgentCoordination.md` §6 Wave 2 exit criteria, still open). It is recorded
here purely as a labelled characterization baseline for future Gate 2
changes to be compared against, per task instruction #8.

MFT denominator (identical convention to the legacy baseline, applied to
the common-CA output): primary MC particle with >=4 correct MFT cluster
labels (no per-disk uniqueness check). Defined once in
`extract_metrics.C`'s `extractMFT()` and copied verbatim into
`extract_metrics_common_ca.C`'s `extractMFTCommonCA()` (only the JSON key
and the denominator-definition string differ, to keep common-CA and legacy
results unambiguous when both JSON files are read side by side).

Raw JSON for every run: `/private/tmp/o2-itsmft-gate0/metrics-common-ca/*.json`
(external, not committed -- regenerate with `extract_metrics_common_ca.C`).

## Repeatability (canonical: fixed input, fixed commit, 1 thread)

4 canonical replay runs (`canonical-1`..`canonical-4`,
`MFTCATrackerParam.nThreads=1`) plus 1 warm-up run, all from the identical
saved cluster fixture, in a separate replay directory tree
(`/private/tmp/o2-itsmft-gate0/replay-common-ca/`) from the legacy
baseline's `fixture-20ev/replay/`.

**Result: every integer and floating-point metric above was bit-identical
across all 5 single-threaded runs** (track count, clusters-per-track
summary, chi2 mean/median/min/max, efficiency, fake rate, clone rate). No
discrepancy observed.

`mfttracks.root`'s own file bytes differ run-to-run (verified with
`shasum -a 256` on all runs) even though the tracking output is
bit-identical: a ROOT `TFile` embeds a per-write UUID/timestamp in its
metadata, so a raw-file hash is not meaningful repeatability evidence here.
`extract_metrics_common_ca.C` therefore also computes `trackContentHash`:
an MD5 over the ordered per-track `(nPoints, chi2, x, y, z, phi, tanl,
invQPt)` tuples (`%.9g`-formatted), which is stable content-level evidence
independent of ROOT container metadata. **This hash was identical
(`826dc653cd936a472929c600c97c140b`) across all 6 runs**, including the
12-thread run below -- stronger evidence than the summary statistics alone,
since it is sensitive to e.g. track reordering that could otherwise leave
mean/median/min/max unchanged.

## Parallel characterization (`MFTCATrackerParam.nThreads=12`)

One replay run (`parallel-1`) at 12 threads (all local cores), compared
against the canonical single-threaded runs. `MFTCATrackerParam.nThreads` is
a genuine thread-count knob for the common-CA core (unlike the legacy
`MFTTrackingParam`, which has no `nThreads` field at all -- legacy MFT
tracking is inherently single-threaded, so this comparison has no legacy
MFT equivalent).

**Result: identical to canonical, no observed difference in any output
metric or the content hash.** Per task instruction #7 this is recorded as
characterization, not a scaling or determinism proof: no per-invocation
tracker-only compute-time breakdown is printed by
`o2-mft-ca-tracker-workflow` (unlike `o2-its-reco-workflow`'s "Tracker
summary: TOT=... s" line used by the Gate 0 baseline to show DPL/CCDB
startup dominates at this scale), so this run cannot itself demonstrate
that tracker time is a meaningful fraction of wall time at this 20-event,
single-TF scale; the identical wall time and RSS between `parallel-1` and
the canonical runs (see performance table below) suggest it likely does
not, consistent with the Gate 0 baseline's finding for the ITS tracker at
the same fixture size.

## Performance (median of 4 canonical runs, 1 discarded warm-up)

Hardware: Apple M3 Pro, 12 cores, 36 GiB RAM, macOS 26.3 (Darwin 25.3.0),
arm64 -- identical machine to the Gate 0 baseline. Build type:
RelWithDebInfo. Commit `43011425b0`. Thread count: 1
(`MFTCATrackerParam.nThreads=1`). Built in a dedicated build directory
configured from this worktree (`/private/tmp/o2-mft-ca-baseline-build`),
not the Gate 0 baseline's build directory.

| Run | Wall time | Peak RSS |
|---|---|---|
| MFT common-CA (cluster-reader \| ca-tracker-workflow), median of 4 canonical | 2.87 s | 1288.6 MB |
| MFT common-CA, range across 4 canonical | 2.66-3.30 s | 1288.2-1289.0 MB |
| MFT common-CA, 12-thread (`parallel-1`, characterization only) | 2.77 s | 1289.2 MB |
| MFT legacy (Gate 0, for reference only) | 2.54 s | 408.5 MB |

The common-CA core's peak RSS (~1.29 GB) is roughly 3x the legacy MFT
tracker's (~409 MB) on the identical fixture. This is a genuine
characterization finding, not a bug fixed by this work (no production code
was modified) -- worth flagging for whoever next tunes the common core's
memory footprint, since it may reflect the shared `o::itsmft::tracking`
core's index-table sizing (`LUTbinsU`/`LUTbinsV` = 64/128 for MFT) or other
common-core allocations not present in the legacy prototype. Wall time is
comparable between the two paths (~2.9 s vs ~2.5 s) but both are dominated
by DPL/CCDB process-startup overhead at this fixture size, not tracking
compute time, per the parallel-characterization note above.

Commands (see `manifest_common_ca.json` for the exact per-run invocation):
```
o2-mft-cluster-reader-workflow --with-mc --input-dir <fixture> \
    --mft-cluster-infile mftclusters.root |
o2-mft-ca-tracker-workflow -b --run --tracking-mode sync \
    --condition-not-after 1784207296000 \
    --configKeyValues 'MFTCATrackerParam.nThreads=1'
```
Wrapped with `/usr/bin/time -l`.

These numbers are characterization data / a trend baseline, not a portable
CI pass/fail threshold, and include DPL process-startup and CCDB-fetch
overhead, not just tracking compute time -- same caveat as the Gate 0
baseline.

## Known issues / limitations

1. **No tracker-only timing breakdown.** `o2-mft-ca-tracker-workflow` does
   not print a per-invocation tracker-compute-time summary (unlike
   `o2-its-reco-workflow`), so the 12-thread run's characterization-only
   status is inferred from identical wall time/RSS rather than demonstrated
   directly the way the Gate 0 baseline did for ITS.
2. **Substantially lower efficiency and fewer tracks than the legacy MFT
   path on the identical fixture/denominator.** Expected at this point in
   Gate 2 (algorithm parity is an open Wave 2 exit criterion per
   `AgentCoordination.md`), not evidence of a defect in this
   characterization work; recorded as a labelled comparison only, per task
   instruction #8.
3. **~3x higher peak RSS than the legacy MFT tracker.** Recorded as an
   observation for future tuning; out of scope to investigate or fix here
   (no production code was modified by this baseline work).
4. All limitations already documented in `baseline_summary.md` for the
   shared fixture (CCDB ETag network dependency, no local CCDB snapshot,
   installed sim/digi stack provenance gap, parallel characterization being
   inconclusive at this fixture size) apply identically here, since the
   fixture and its generation are entirely reused and unmodified.

## Files produced by this characterization

- `replay_tracking_common_ca.sh` -- replay script (sibling to
  `replay_tracking.sh`).
- `extract_metrics_common_ca.C` -- metric extractor (sibling to
  `extract_metrics.C`; reuses its `extractMFT()` logic verbatim).
- `manifest_common_ca.json` -- full provenance (sibling to `manifest.json`).
- `baseline_summary_common_ca.md` -- this file (sibling to
  `baseline_summary.md`).

No existing Gate 0 file (`README.md`, `manifest.json`,
`baseline_summary.md`, `generate_fixture.sh`, `replay_tracking.sh`,
`extract_metrics.C`) was modified. No production tracking code, parameters,
workflows, kernels, or outputs were modified. No binary fixture or
generated output is committed; all replay/metrics directories referenced
above live under `/private/tmp/o2-itsmft-gate0/` outside the git worktree.

**Not marking Gate 2 complete** -- the integration owner reviews this
evidence, per `AgentCoordination.md` §12 ("An agent may not declare a gate
passed based only on compilation").

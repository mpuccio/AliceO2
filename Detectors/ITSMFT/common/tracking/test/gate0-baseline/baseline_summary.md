# Gate 0 physics/performance baseline summary

Commit under test: `e4d35cb1a70c0254130e87bc3d0caab3cd567e75`
("ITSMFT: accept normalized measurement and TimeFrame contracts" -- production
ITS/MFT tracking is still the pre-refactor baseline implementation; see
`manifest.json`).

Fixture: 20 unanchored pp events (`pythia8pp`/`TGeant4`), PIPE+ITS+MFT only,
run 303000, seed 20260716. Full provenance (commit, build paths, CCDB
objects/ETags, exact commands) in `manifest.json`. Twenty pp events are a
regression smoke sample, not a statistically authoritative physics
validation.

Gate 0 does **not** cover: architecture decisions, tracking parameter
changes, or any production code review -- those are separately owned per
`AgentCoordination.md`. This is characterization data only.

## Metrics (canonical: single-threaded replay from fixed clusters)

Denominators are defined explicitly per-detector in `extract_metrics.C` and
repeated in each metric JSON's `denominatorDefinition` field.

| Metric | ITS | MFT |
|---|---|---|
| Input clusters | 5378 | 1679 |
| Input ROFs | 2304 | 2304 |
| Output tracks | 252 | 149 |
| Clusters/track (mean / median / min / max) | 6.10 / 7 / 4 / 7 | 7.50 / 8 / 5 / 10 |
| chi2 (mean / median / min / max) | 7.61 / 5.43 / 0.05 / 105.86 | 1.61 / 0.21 / 0.006 / 116.78 |
| MC reconstructable | 201 | 141 |
| Matched | 169 | 132 |
| Efficiency | 0.8408 | 0.9362 |
| Fake tracks / rate | 17 / 0.0675 | 5 / 0.0336 |
| Clone tracks / rate | 0 / 0 | 0 / 0 |

ITS denominator: primary MC particle with a correct cluster label on all 7
ITS layers (bitmask `0x7f`), matching the existing convention in
`Detectors/ITSMFT/ITS/macros/test/CheckTracksCA.C`.

MFT denominator: primary MC particle with >=4 correct MFT cluster labels
(no per-disk uniqueness check -- simpler than the ITS convention since no
MFT `GeometryTGeo` layer lookup was implemented here; see
`extract_metrics.C`).

Raw JSON for every run: `/private/tmp/o2-itsmft-gate0/metrics/*.json`
(external, not committed -- regenerate with `extract_metrics.C`).

## Repeatability (canonical: fixed input, fixed commit, 1 thread)

4 canonical replay runs (`canonical-1`..`canonical-4`, `ITSCATrackerParam.nThreads=1`)
plus 1 warm-up run, all from the identical saved cluster fixture.

**Result: every integer and floating-point metric above was bit-identical
across all 5 single-threaded runs** (track counts, clusters-per-track
summary, chi2 mean/median/min/max, efficiency, fake rate, clone rate). No
discrepancy observed.

## Parallel characterization (`ITSCATrackerParam.nThreads=12`)

One replay run (`parallel-1`) at 12 threads (all local cores), compared
against the canonical single-threaded runs.

**Result: identical to canonical, no observed difference in any output
metric.** This is expected to be inconclusive rather than a strong
determinism claim: per-invocation logs show the ITS tracker's own compute
time is ~0.04s of the ~3.6s measured wall time (`Tracker summary: Processed
1 TFs ... TOT=0.04 s`) -- DPL/CCDB startup overhead dominates at this
20-event, single-TF scale, leaving too little actual parallel tracking work
to exercise thread-count-dependent floating-point summation-order effects.
A larger fixture would be needed to characterize this properly; per policy
this is recorded as an observation, not treated as a determinism gate.

## Performance (median of 4 canonical runs, 1 discarded warm-up)

Hardware: Apple M3 Pro, 12 cores, 36 GiB RAM, macOS 26.3 (Darwin 25.3.0),
arm64. Build type: RelWithDebInfo (FairMQ-reported; O2 build type not
independently re-verified). Commit `e4d35cb1a7`. Thread count: 1
(`ITSCATrackerParam.nThreads=1`).

| Detector | Median wall time | Range | Median peak RSS | Range |
|---|---|---|---|---|
| ITS (cluster-reader \| its-reco-workflow) | 3.64 s | 3.51-3.93 s | 670.3 MB | 670.0-671.0 MB |
| MFT (cluster-reader \| mft-reco-workflow) | 2.54 s | 2.44-2.58 s | 408.5 MB | 408.0-408.7 MB |

Commands (see `replay_tracking.sh`):
```
o2-its-cluster-reader-workflow --with-mc --input-dir <fixture> \
    --its-cluster-infile o2clus_its.root |
o2-its-reco-workflow -b --run --clusters-from-upstream --tracking-mode sync \
    --cluster-rof-branch-only --condition-not-after 1784207296000 \
    --configKeyValues 'ITSCATrackerParam.nThreads=1' \
    --its-track-writer '--outfile o2trac_its.root'
```
(and the MFT equivalent). Wrapped with `/usr/bin/time -l`.

These numbers are characterization data / a trend baseline, not a portable
CI pass/fail threshold, and include DPL process-startup and CCDB-fetch
overhead, not just tracking compute time (see parallel-characterization
note above for the actual tracker-only compute time).

## Known issues / limitations

1. **`--hbfutils-config <absolute-ini-path>` hangs.** Passing FIXTURE_DIR's
   `o2simdigitizerworkflow_configuration.ini` to `o2-its-reco-workflow` via
   an explicit absolute path reproducibly hung the driver indefinitely
   before it spawned any device, when reading clusters through
   `o2-its-cluster-reader-workflow` piped upstream. Root cause not
   identified further given the time budget. Workaround (used by
   `replay_tracking.sh`): copy the ini file into `REPLAY_DIR` under its
   *default* relative filename and omit `--hbfutils-config` entirely.
2. **`o2-sim`/`o2-sim-digitizer-workflow` come from a separately-built,
   possibly-modified installed stack** (`o2-sim` reports version `1.2.0
   (cf19214f6e (+local changes))` -- an older commit than the one under
   test, with unspecified local modifications), per task instructions
   (simulation/digitization is not part of the ITS/MFT tracking refactor).
   This is a provenance gap for the simulation/digitization stage
   specifically, not for the tracking code being characterized.
3. **No local CCDB snapshot.** Reproducibility depends on the recorded
   ETags in `manifest.json` remaining available on `alice-ccdb.cern.ch`.
   Verified stable across fixture generation and 5 replay runs in this
   session; not verified stable over longer time spans.
4. **Local build omitted two dlopen'd DPL plugin libraries** by default
   (`libO2FrameworkAnalysisSupport.dylib`, `libO2FrameworkCCDBSupport.dylib`)
   since `ninja` was only asked to build the named executables; the CCDB
   one is not just cosmetic (its absence silently breaks condition
   fetching, i.e. tracking data-completeness) -- see README "Required
   environment" section. Both were built explicitly for this baseline.
5. **Machine-specific Geant4 data env vars.** `GEANT4`'s `init.sh` did not
   export `G4*DATA` variables despite the data files being present on
   disk; `o2-sim` aborts without them. Worked around locally by exporting
   them explicitly; not fixed upstream (out of scope for this role).
6. **Parallel characterization is inconclusive** at this fixture scale
   (see above) -- do not read the identical parallel/canonical numbers as
   a validated determinism guarantee for larger, tracker-compute-dominated
   workloads.

## Remaining blockers for future runs

- A CI-facing invocation of this pipeline still needs a decision on where
  fixtures/CCDB objects are cached (this session ran everything under
  `/private/tmp`, ephemeral).
- Thread-scaling characterization needs a larger fixture (more events
  and/or higher interaction rate) so that tracker compute time is not
  dwarfed by DPL/CCDB startup overhead.
- The MFT reconstructable-track denominator does not check per-disk
  uniqueness (unlike the ITS 7-layer bitmask); adding an MFT
  `GeometryTGeo`-based layer lookup would align the two conventions if a
  future gate needs that.

Not marking Gate 0 complete -- the integration owner reviews this evidence.

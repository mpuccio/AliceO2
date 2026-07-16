# Gate 0 physics/performance baseline

Scripts to generate a small ITS+MFT tracker-input fixture and characterize
the baseline (pre-refactor) reconstruction, per
`Detectors/ITSMFT/common/tracking/doc/AgentCoordination.md` §12/§15 and
`Architecture.md` Gate 0.

Not wired into CTest: this is an expensive, network-dependent workflow
(CCDB access, `o2-sim`), not a unit test. It is meant to be run by hand or
from CI infrastructure that explicitly opts in.

## What this does, and does not, characterize

The reconstruction binaries under test (`o2-its-reco-workflow`,
`o2-mft-reco-workflow`, `o2-mft-assessment-workflow`) must come from the
commit under review — see `manifest.json` for the exact commit and build
paths used to produce `baseline_summary.md`. `o2-sim` and
`o2-sim-digitizer-workflow` are simulation/digitization infrastructure, not
part of the ITS/MFT tracking refactor; they may legitimately come from a
separately-built/installed O2 stack (see manifest for that provenance too).

The **tracker input boundary** — compact clusters, patterns, cluster ROF
records, MC cluster labels — is the reusable fixture. Tracks are a baseline
*output*, not an input: comparisons across code changes must replay tracking
from the fixed cluster-level fixture (`replay_tracking.sh`), not regenerate
simulation and clustering each time.

## Scripts

- `generate_fixture.sh` — sim (`o2-sim`, PIPE+ITS+MFT only) → digitization
  (ITS,MFT only) → one initial `o2-its-reco-workflow` /
  `o2-mft-reco-workflow` pass. This initial reco pass is what actually
  produces the fixture's cluster files (`o2clus_its.root`,
  `mftclusters.root`) — ITS/MFT clusterization is fused into the same
  binary as tracking, there is no separate clusterizer executable.
- `replay_tracking.sh` — reruns *only* tracking from a fixture's saved
  clusters (`o2-its-cluster-reader-workflow | o2-its-reco-workflow
  --clusters-from-upstream`, and the MFT equivalent). This is what
  canonical-repeatability, parallel-characterization, and performance runs
  all call.

Both scripts assume the caller has already put the right `o2-*` binaries on
`PATH`/`LD_LIBRARY_PATH`/`DYLD_LIBRARY_PATH` (see manifest for the exact
`source .../init.sh` + build-dir prefixing used for the recorded baseline).
Neither script hardcodes a machine path.

## Required environment for this specific machine's build

Two DPL plugin libraries (`libO2FrameworkAnalysisSupport.dylib`,
`libO2FrameworkCCDBSupport.dylib`) are loaded via `dlopen` at runtime and are
*not* linked dependencies of `o2-its-reco-workflow` / `o2-mft-reco-workflow`
/ `o2-mft-assessment-workflow`, so a `ninja` build that only targets those
three executables silently omits them (the workflow keeps running with the
CCDB-fetch device dead, which is a correctness risk, not just a log-noise
issue). Build them explicitly alongside the three targets:

```
ninja stage/lib/libO2FrameworkAnalysisSupport.dylib \
      stage/lib/libO2FrameworkCCDBSupport.dylib
```

Also, on this machine, `GEANT4`'s `etc/profile.d/init.sh` does not export the
`G4*DATA` environment variables even though the data files are present
under `$GEANT4_ROOT/share/Geant4/data/`; `o2-sim` aborts
(`G4ENSDFSTATEDATA environment variable must be set`) without them. Export
them from that directory before calling `generate_fixture.sh` (see
`manifest.json` for the exact paths recorded for this baseline run).

`replay_tracking.sh`'s two DPL pipelines (`cluster-reader-workflow |
reco-workflow`) each run inside their own `bash -o pipefail -c "..."`.
This is required, not decorative: the outer script's own `set -o
pipefail` only applies to pipelines written directly in that script, not
to a nested `bash -c` subshell it spawns, which starts with pipefail off
by default. Without it, a crash in the upstream reader would be masked by
the downstream reco-workflow's exit code, since a plain (non-pipefail)
pipe reports only the last command's status.

Finally, `replay_tracking.sh` copies `FIXTURE_DIR`'s
`o2simdigitizerworkflow_configuration.ini` into `REPLAY_DIR` and relies on
`--hbfutils-config`'s *default* value (the same filename, resolved
relative to cwd) rather than passing an explicit path. Passing that ini
file's absolute path via `--hbfutils-config <path>` reproducibly hung
`o2-its-reco-workflow` indefinitely before it spawned any device when
reading clusters through the piped `o2-its-cluster-reader-workflow`; root
cause not identified further given the time budget. Do not "simplify" this
back to an explicit path without re-testing.

## Fresh-directory enforcement

Both scripts refuse to run if their output directory already exists and is
non-empty (`FIXTURE_DIR` for `generate_fixture.sh`, `REPLAY_DIR` for
`replay_tracking.sh`): they exit 1 with a message rather than silently
mixing a previous run's files with a new one. Remove the directory (or
point at a fresh path) to regenerate/re-replay. Both scripts also validate
every required output file as non-empty before exiting 0 -- a partially
failed run cannot look like a success just because some `ls` printed file
names.

## Variables

`generate_fixture.sh`: `FIXTURE_DIR`, `RUNNUMBER`, `SEED`, `TIMESTAMP`
(CCDB `condition-not-after`, ms), `NEVENTS`; optional `SIM_WORKERS`,
`SHM_SEGMENT_SIZE`, `INTERACTION_RATE_HZ`.

`replay_tracking.sh`: `FIXTURE_DIR`, `REPLAY_DIR`, `TIMESTAMP`, `RUNNUMBER`
(sanity-checked against `FIXTURE_DIR`'s HBFUtils ini); optional
`ITS_NTHREADS`, `SHM_SEGMENT_SIZE`, `TIME_CMD`.

## Determinism model

- **Canonical baseline**: fixed input fixture, fixed commit, one thread
  (`ITS_NTHREADS=1`). Integer outputs (cluster/track/ROF counts, label
  classifications) are expected exact between repeated canonical runs.
- **Parallel characterization**: same fixture/commit, a fixed practical
  thread count. Differences from canonical are observations, not a gate.
- Floating-point outputs (chi2, pt, eta, covariance) are compared as
  distributions/summaries with a documented tolerance, never bitwise.
- Runtime/memory are trend data (median + observed range across ≥3 runs
  after one warm-up), not a pass/fail CI threshold.

See `baseline_summary.md` (generated per run) for the actual recorded
numbers and `manifest.json` for full provenance.

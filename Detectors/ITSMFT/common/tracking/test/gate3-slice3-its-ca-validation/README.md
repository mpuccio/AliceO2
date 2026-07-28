# Gate 3 Slice 3: ITS common-CA validation (sibling to gate0-baseline)

Scripts and results characterizing `o2-its-ca-tracker-workflow` (the opt-in
ITS common-CA workflow built in Gate 3 workflow-onboarding Slice 2, commit
`4a2bf6c364` at the time this validation ran) against the durable Gate 0
fixture, per this validation task's scope: build/verify/characterize the ITS
common-CA leg, pair it with a legacy ITS leg under an equivalent explicit
diamond constraint, check determinism, check negative (fail-closed)
behavior, and run the existing MFT common-CA non-regression fence unchanged.

**This is characterization, not a bitwise-parity gate.** No physics
tolerance is invented here, and no claim is made that the common-CA leg's
numbers match the legacy leg's, or that the explicit static-diamond
constraint used on both legs is equivalent to the legacy real-vertexer or
CCDB MeanVertex-override Sync defaults.

**Production code is untouched.** This directory contains only validation
scripts, a metrics extractor, and result documentation. No file under
`Detectors/ITSMFT/common/tracking/{src,include}`,
`Detectors/ITSMFT/ITS/workflow-ca`, `Detectors/ITSMFT/ITS/tracking`, or
`gate0-baseline` was modified. No ROOT/binary fixture or replay output is
committed; all replay/metrics directories referenced below live under
`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/gate3-slice3-its-ca-validation/`,
outside the git worktree.

## Scripts

- `replay_tracking_its_common_ca.sh` -- replays `o2-its-cluster-reader-workflow
  | o2-its-ca-tracker-workflow` from the fixture's saved ITS clusters, with
  an explicit `ITSCommonCATrackerParam.useDiamond=true` +
  `diamondPos[0..2]` + `pvRes` constraint (values recorded, not left at
  struct defaults). Produces `o2trac_its_ca.root`.
- `replay_tracking_its_legacy_diamond.sh` -- pairs the above: replays
  `o2-its-cluster-reader-workflow | o2-its-reco-workflow` (ITS-only, no
  MFT -- this validation's build scope is legacy `o2-its-reco-workflow`,
  not legacy MFT tracking) with the identical explicit diamond values via
  the legacy `ITSCATrackerParam` namespace. Produces `o2trac_its.root`.
  Sibling to `gate0-baseline/replay_tracking.sh`, not a replacement --
  no file in `gate0-baseline` is read, written, or modified.
- `negative_checks_its_common_ca.sh` -- runs `o2-its-ca-tracker-workflow`
  standalone under two configurations that must both fatal before any
  device is constructed and before any output file is created: (1) default
  configuration (`useDiamond` stays false), (2) a legacy
  `ITSCATrackerParam.*` `--configKeyValues` override. Asserts exit code,
  expected fatal-message substring, and output-file absence for both.
- `extract_metrics_its_common_ca.C` -- metrics extractor for
  `o2trac_its_ca.root`: track count, clusters/track and chi2 summaries, an
  ordered semantic track-content digest (see field list below), and
  MC-label availability/alignment. Sibling to
  `gate0-baseline/extract_metrics_common_ca.C` (the MFT equivalent); same
  denominator convention as `gate0-baseline/extract_metrics.C`'s
  `extractITS()`.
- `extract_metrics_its_legacy_diamond.sh` -- thin wrapper that loads
  `gate0-baseline/extract_metrics.C` via ROOT's `.L` (unmodified, not
  copied) and calls its `extractITS()` directly, since that file's own
  `extract_metrics()` entry point unconditionally also requires an MFT
  replay output this validation's legacy leg does not produce.

The existing `gate0-baseline/replay_tracking_common_ca.sh` +
`gate0-baseline/extract_metrics_common_ca.C` (MFT common-CA) are reused
UNCHANGED for the MFT non-regression fence (see `manifest.json` and the
characterization summary below) -- neither file lives in this directory.

## Track-content digest field list

`extract_metrics_its_common_ca.C`'s `trackContentHash` is an MD5 over the
ordered per-track tuple `(nClusters, chi2, x, alpha, y, z, snp, tgl, q2pt)`,
`%.9g`-formatted. `nClusters`/`chi2` are `o2::its::TrackITS`'s own summary
accessors (`getNumberOfClusters()`/`getChi2()`); the other seven are exactly
the fields `o2::track::TrackParametrization::hash()` itself hashes over --
i.e. this is the framework's own notion of a barrel track's parametrized
state, not an ad hoc field choice, and it matches the *documented* ITS
track-representation fields the task asked this digest to cover. It does
NOT cover every branch the writer emits (`ITSTrackClusIdx`/`ITSTracksROF`
are excluded) and it is NOT a hash of `o2trac_its_ca.root`'s raw bytes,
which differ run-to-run purely from ROOT `TFile` UUID/timestamp metadata
even when the semantic content is identical (verified explicitly below).

## Running

See `manifest.json` for the exact package/fixture/timestamp/output-location
values and full command lines used to produce the results in
`characterization_summary.md`. In short, all scripts are invoked through
`.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh` with
`O2_BUILD_DIR` pointed at a worktree build of this exact commit and
`O2_PACKAGE=daily-20260717-0700-local1` pinned explicitly (never an alias),
per that skill's reproducible-validation policy.

## Known flaky-hang observation (operational, not a script bug)

During this validation, `replay_tracking_its_legacy_diamond.sh`'s piped
`o2-its-cluster-reader-workflow | o2-its-reco-workflow` invocation hung
indefinitely (no CPU progress, driver stuck in its event loop after only
printing DPL-initialisation log lines) on 2 of 4 attempts, all under the
same `REPLAY_DIR` tree
(`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/gate3-slice3-its-ca-validation/its-legacy-diamond-1`).
The identical command (same `--configKeyValues`, same fixture) run
individually in scratch directories, and a subsequent bare retry of the
same script/`REPLAY_DIR`, completed normally in well under a minute each
time. No orphaned FairMQ shared-memory segments (`ipcs -m`) or stuck
processes were found before any retry. Root cause not identified further
within this task's scope (no production code was touched to investigate
it); recorded here as an operational flakiness observation for whoever next
runs this script -- retry once on an apparent hang before treating it as a
real failure.

# Slice C traversal activation validation

This note records the remaining acceptance evidence for Gate 2 Slice C. It
characterizes the deterministic sparse cell schedule at `9c3dc9ba6e` and the
`findCellsNeighbours` activation at `caa1b6e299`. No production code was
changed during this validation.

All O2 builds, tests, ROOT commands, workflows, and measurements were run
through `.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh` with
`ALIBUILD_WORK_DIR=/Users/mpuccio/alice/sw`, architecture `osx_arm64`, package
`latest-run3-o2`, and ROOT 6.36.10. Both compared builds are
`RelWithDebInfo`. The retained 20-event fixture, run 303000, condition cutoff
1784207296000, and `MFTCATrackerParam.nThreads=1` were used throughout.
The host is Darwin 25.3.0 arm64 on the same Apple M3 Pro, 12-core, 36 GiB
machine recorded by the canonical baseline manifest.

The compared staged builds were:

- canonical common-CA baseline commit `43011425b0`, retained build
  `/private/tmp/o2-mft-ca-baseline-build`;
- Slice C commit `caa1b6e299`, build
  `/private/tmp/o2-topology-cell-neighbours-build`.

The canonical build's CMake cache identifies the original source checkout and
commit recorded by `manifest_common_ca.json`; that checkout is no longer
present, so this comparison uses the retained, previously validated staged
binaries rather than rebuilding the canonical revision.

The Slice C workflow and runtime plugin targets used here were built with:

```sh
O2_BUILD_DIR=/private/tmp/o2-topology-cell-neighbours-build \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  ninja -C /private/tmp/o2-topology-cell-neighbours-build \
    O2exe-mft-ca-tracker-workflow O2exe-mft-cluster-reader-workflow \
    O2lib-FrameworkAnalysisSupport O2lib-FrameworkCCDBSupport
```

## Fail-closed negative path

The negative check used the fresh directory
`/private/tmp/o2-itsmft-slicec-validation/negative-invalid-nsigma`. It injected
an infinite MFT policy multiplier, which survives configuration mapping as an
invalid bound policy parameter:

```sh
O2_BUILD_DIR=/private/tmp/o2-topology-cell-neighbours-build \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  bash -lc 'cd /private/tmp/o2-itsmft-slicec-validation/negative-invalid-nsigma; \
    set -o pipefail; \
    o2-mft-cluster-reader-workflow --with-mc \
      --input-dir /private/tmp/o2-itsmft-gate0/fixture-20ev \
      --mft-cluster-infile mftclusters.root --shm-segment-size 4000000000 | \
    o2-mft-ca-tracker-workflow -b --run --tracking-mode sync \
      --shm-segment-size 4000000000 \
      --condition-not-after 1784207296000 \
      --configKeyValues \
      "MFTCATrackerParam.nThreads=1;MFTCATrackerParam.nSigmaCut=inf"'
```

The pipeline exited 128. `workflow.log` identifies the failure as
`CA traversal initialization failed at iteration 0 (reason=7)`; reason 7 is
`TraversalFailureReason::InvalidPolicyParameters`. The tracker did not send a
timeframe to the writer.

The independently initialized DPL track-writer process creates its ROOT
container before the producer runs, so a 7,793-byte `mfttracks.root` crash
artifact exists. A ROOT check found a readable `o2sim` tree with exactly zero
entries. Thus there is no apparently successful timeframe output, although
the filename is not physically absent. The nonzero workflow status, explicit
typed traversal diagnostic, and zero published tree entries are the practical
fail-closed evidence for this multi-process workflow. Transactional temporary
file/rename behavior in the generic writer is outside Slice C.

The inspection command was:

```sh
O2_BUILD_DIR=/private/tmp/o2-topology-cell-neighbours-build \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  root -l -b -q -e \
  'auto f=TFile::Open("/private/tmp/o2-itsmft-slicec-validation/negative-invalid-nsigma/mfttracks.root"); \
   auto t=dynamic_cast<TTree*>(f->Get("o2sim")); \
   std::cout << "entries=" << (t ? t->GetEntries() : -1) << std::endl;'
```

## Focused performance characterization

One warm-up was run for each build, followed by 20 measured pairs. Odd pairs
ran baseline then Slice C; even pairs reversed the order. Every run used the
same fixture, environment, one-thread setting, replay script, and command;
only `O2_BUILD_DIR` and the fresh output directory changed. All 42 workflows
reached the replay script's output validation, produced nonempty normal-run
ROOT output, and had no workflow or traversal failure marker.

The per-run command was:

```sh
FIXTURE_DIR=/private/tmp/o2-itsmft-gate0/fixture-20ev \
REPLAY_DIR=/private/tmp/o2-itsmft-slicec-validation/perf-paired20-a/<run> \
TIMESTAMP=1784207296000 RUNNUMBER=303000 MFT_CA_NTHREADS=1 \
TIME_CMD='/usr/bin/time -l' O2_BUILD_DIR=<compared-build> \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  bash Detectors/ITSMFT/common/tracking/test/gate0-baseline/replay_tracking_common_ca.sh
```

Times cover the entire cluster-reader, CA tracker, track-writer, DPL startup,
and condition-fetch pipeline. RSS is the command-level maximum resident set
reported by macOS `/usr/bin/time -l`, converted using decimal MB. Spread is
reported as the linearly interpolated 25th--75th percentile interval (IQR),
with the full observed range retained.

| Variant | Warm-up wall / RSS | Median wall | Wall IQR; range | Median peak RSS | RSS IQR; range |
|---|---:|---:|---:|---:|---:|
| Canonical baseline | 10.86 s / 1288.487 MB | 8.635 s | 8.448--9.853 s; 8.30--73.46 s | 1288.618 MB | 1288.450--1288.852 MB; 1288.307--1289.454 MB |
| Slice C | 8.38 s / 1288.569 MB | 8.465 s | 8.408--8.728 s; 8.30--68.70 s | 1288.634 MB | 1288.499--1288.782 MB; 1288.258--1289.929 MB |

Slice C minus baseline has a paired median wall delta of -0.050 s, IQR
-0.858 to +0.238 s, and full range -65.08 to +57.55 s. The paired median RSS
delta is +0.098 MB, IQR -0.209 to +0.225 MB, and full range -1.180 to
+1.245 MB. Comparing the two unpaired medians gives -0.170 s (-2.0%) wall
time and +0.016 MB peak RSS. These data show no focused performance or memory
regression at this resolution.

The full ranges contain deliberately retained external/startup outliers:
68.70 s for Slice C pair 14 and 73.46 s for baseline pair 20, plus 25.31 s
for baseline pair 8. CPU work and RSS remained near the other runs. Because
the workflow has no tracker-only timing and condition/DPL startup dominates
this small fixture, this is characterization of the complete command, not a
scaling proof or a portable performance threshold.

## Once-per-iteration cache evidence

The focused cache test was rerun directly:

```sh
O2_BUILD_DIR=/private/tmp/o2-topology-cell-neighbours-build \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  /private/tmp/o2-topology-cell-neighbours-build/stage/tests/\
o2-test-itsmft-tracking-timeframe-detector-layouts \
  --run_test=traversal_cache_groups_and_binds_once_across_repeated_neighbour_calls \
  --log_level=test_suite
```

It passed. The test observes one grouping construction and one disk-policy
binding after initialization, calls `findCellsNeighbours(0)` twice, and
confirms both counters remain one. This directly covers repeated per-vertex
neighbor invocation within one iteration; the production workflow does not
export these internal counters.

All output directories and raw logs referenced above are external validation
artifacts under `/private/tmp` and are not committed.

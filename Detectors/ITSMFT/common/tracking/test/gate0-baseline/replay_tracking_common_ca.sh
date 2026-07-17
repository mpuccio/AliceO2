#!/usr/bin/env bash
# Gate 2 common-CA MFT baseline: replay o2-mft-ca-tracker-workflow (the
# o2::itsmft::tracking common CA core, NOT the legacy O2::MFTTracking package
# exercised by replay_tracking.sh) from the same fixed cluster-level fixture
# used by the accepted Gate 0 baseline.
#
# This is a SIBLING to replay_tracking.sh, not a replacement: it exists
# because Gate 0's accepted baseline exercises o2-mft-reco-workflow (legacy
# O2::MFTTracking), while Gate 2 changes the separate common
# o2::itsmft::tracking core used by o2-mft-ca-tracker-workflow. No accepted
# baseline previously covered that path. Legacy metrics/manifest/summary
# files are never modified by this script or its outputs.
#
# Unlike o2-mft-reco-workflow, o2-mft-ca-tracker-workflow has no embedded
# clusterizer or cluster-writer device -- its only inputs are
# MFT/COMPCLUSTERS, MFT/PATTERNS, MFT/CLUSTERSROF (and MFT/CLUSTERSMCTR with
# MC), identical InputSpecs to the legacy TrackerSpec, so it pipes directly
# after o2-mft-cluster-reader-workflow with no --clusters-from-upstream or
# --cluster-rof-branch-only flag (neither exists on this workflow; there is
# no competing clusterizer or cluster-writer device to disambiguate).
#
# Required variables:
#   FIXTURE_DIR   directory produced by generate_fixture.sh (read-only here;
#                 shared with the legacy replay, never written to).
#   REPLAY_DIR    output directory for this replay. Must not already exist
#                 as a non-empty directory -- refuses to run into a dirty
#                 REPLAY_DIR rather than silently mixing runs. Remove it, or
#                 point at a fresh path, to re-replay.
#   TIMESTAMP     same fixed CCDB condition-not-after value used at fixture
#                 generation.
#   RUNNUMBER     same fixed run number used at fixture generation
#                 (sanity-checked against FIXTURE_DIR's HBFUtils ini).
#
# HBFUtils ini handling: identical workaround to replay_tracking.sh --
# --hbfutils-config pointed at an absolute path reproducibly hung the ITS
# reco driver before it spawned any device when reading clusters through a
# piped cluster-reader-workflow (root cause not identified further given the
# time budget there); this script inherits the same defensive pattern of
# copying the ini into REPLAY_DIR under its default relative name and
# omitting --hbfutils-config, even though it has only been directly observed
# to hang for the ITS binary, not this one.
#
# Optional:
#   MFT_CA_NTHREADS   MFTCATrackerParam.nThreads (default 1). This is a
#                     genuine common-CA thread-count knob; the legacy
#                     MFTTrackingParam has no such field (legacy MFT tracking
#                     is inherently single-threaded), so replay_tracking.sh
#                     has no MFT thread-count equivalent to compare against.
#   SHM_SEGMENT_SIZE  bytes (default 4000000000)
#   TIME_CMD          wrapper prefixed to the reco invocation, e.g.
#                      "/usr/bin/time -l" to capture wall time + peak RSS.
#                      Left empty by default (no wrapping).

set -euo pipefail

: "${FIXTURE_DIR:?set FIXTURE_DIR}"
: "${REPLAY_DIR:?set REPLAY_DIR}"
: "${TIMESTAMP:?set TIMESTAMP}"
: "${RUNNUMBER:?set RUNNUMBER}"
MFT_CA_NTHREADS="${MFT_CA_NTHREADS:-1}"
SHM_SEGMENT_SIZE="${SHM_SEGMENT_SIZE:-4000000000}"
TIME_CMD="${TIME_CMD:-}"

for bin in o2-mft-cluster-reader-workflow o2-mft-ca-tracker-workflow; do
  command -v "$bin" >/dev/null 2>&1 || { echo "missing $bin on PATH" >&2; exit 1; }
done
[[ -f "$FIXTURE_DIR/mftclusters.root" ]] || { echo "missing $FIXTURE_DIR/mftclusters.root" >&2; exit 1; }
HBFUTILS_INI_SRC="$FIXTURE_DIR/o2simdigitizerworkflow_configuration.ini"
[[ -f "$HBFUTILS_INI_SRC" ]] || { echo "missing $HBFUTILS_INI_SRC" >&2; exit 1; }
grep -q "runNumber=$RUNNUMBER" "$HBFUTILS_INI_SRC" || { echo "RUNNUMBER=$RUNNUMBER does not match $HBFUTILS_INI_SRC" >&2; exit 1; }

if [[ -d "$REPLAY_DIR" ]] && [[ -n "$(ls -A "$REPLAY_DIR" 2>/dev/null)" ]]; then
  echo "REPLAY_DIR '$REPLAY_DIR' already exists and is not empty; refusing to reuse it. Remove it or pick a fresh path." >&2
  exit 1
fi
mkdir -p "$REPLAY_DIR"
cp "$HBFUTILS_INI_SRC" "$REPLAY_DIR/o2simdigitizerworkflow_configuration.ini"
cd "$REPLAY_DIR"

echo "[replay_tracking_common_ca] MFT common-CA replay: nthreads=$MFT_CA_NTHREADS"
# Only the last command in a DPL pipe gets -b/--run; the upstream reader
# must stay a bare invocation so it dumps its workflow spec (via stdout) for
# merging instead of executing as its own standalone driver. The inner
# "bash -c" starts a fresh shell that does NOT inherit this script's
# "set -o pipefail" (each bash process gets its own default shell options),
# so it is set explicitly here too -- otherwise a failure in the upstream
# reader would be masked by the downstream tracker's exit code, which is the
# only one a plain pipe reports.
$TIME_CMD bash -o pipefail -c "
  o2-mft-cluster-reader-workflow --with-mc \
    --input-dir '$FIXTURE_DIR' --mft-cluster-infile mftclusters.root \
    --shm-segment-size $SHM_SEGMENT_SIZE |
  o2-mft-ca-tracker-workflow -b --run --tracking-mode sync \
    --shm-segment-size $SHM_SEGMENT_SIZE \
    --condition-not-after $TIMESTAMP \
    --configKeyValues 'MFTCATrackerParam.nThreads=$MFT_CA_NTHREADS'
" > mft_ca_replay.log 2> mft_ca_replay.time.log
[[ -s mfttracks.root ]] || { echo "mft common-CA replay produced no mfttracks.root, see mft_ca_replay.log" >&2; exit 1; }

echo "[replay_tracking_common_ca] validating replay outputs"
[[ -s mfttracks.root ]] || { echo "missing or empty required replay output: mfttracks.root" >&2; exit 1; }
ls -la mfttracks.root

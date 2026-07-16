#!/usr/bin/env bash
# Gate 0 baseline: replay ITS+MFT tracking from a fixed cluster-level fixture.
#
# Reruns only the tracker (clusterizer skipped via --clusters-from-upstream),
# reading the compact clusters/patterns/ROFs/labels saved by
# generate_fixture.sh. This is what canonical repeatability, parallel
# characterization, and performance runs all call, so that "rerun" always
# means "rerun tracking from fixed clusters", not "regenerate simulation".
#
# Required variables:
#   FIXTURE_DIR   directory produced by generate_fixture.sh (read-only here)
#   REPLAY_DIR    output directory for this replay (created, should be empty)
#   TIMESTAMP     same fixed CCDB condition-not-after value used at generation
#   RUNNUMBER     same fixed run number used at generation (sanity-checked
#                 against FIXTURE_DIR's HBFUtils ini).
#
# HBFUtils ini handling: the run-number/orbit-derived CCDB query timestamp
# (needed so e.g. GLO/Config/GRPECS resolves against the run's actual
# historical time instead of wall-clock "now", which 404s) only resolves
# correctly when --hbfutils-config is left at its *default* value and the
# ini file with that exact default name exists as a relative path in cwd.
# Passing the ini file via an explicit --hbfutils-config <absolute-path>
# reproducibly hung o2-its-reco-workflow indefinitely before it spawned any
# device (root cause not identified further, given time budget) -- so this
# script copies FIXTURE_DIR's ini into REPLAY_DIR under its default name
# instead of pointing at it directly.
# Optional:
#   ITS_NTHREADS      ITSCATrackerParam.nThreads (default 1)
#   SHM_SEGMENT_SIZE  bytes (default 4000000000)
#   TIME_CMD          wrapper prefixed to each reco invocation, e.g.
#                      "/usr/bin/time -l" to capture wall time + peak RSS.
#                      Left empty by default (no wrapping).

set -euo pipefail

: "${FIXTURE_DIR:?set FIXTURE_DIR}"
: "${REPLAY_DIR:?set REPLAY_DIR}"
: "${TIMESTAMP:?set TIMESTAMP}"
: "${RUNNUMBER:?set RUNNUMBER}"
ITS_NTHREADS="${ITS_NTHREADS:-1}"
SHM_SEGMENT_SIZE="${SHM_SEGMENT_SIZE:-4000000000}"
TIME_CMD="${TIME_CMD:-}"

for bin in o2-its-cluster-reader-workflow o2-its-reco-workflow o2-mft-cluster-reader-workflow o2-mft-reco-workflow; do
  command -v "$bin" >/dev/null 2>&1 || { echo "missing $bin on PATH" >&2; exit 1; }
done
[[ -f "$FIXTURE_DIR/o2clus_its.root" ]] || { echo "missing $FIXTURE_DIR/o2clus_its.root" >&2; exit 1; }
[[ -f "$FIXTURE_DIR/mftclusters.root" ]] || { echo "missing $FIXTURE_DIR/mftclusters.root" >&2; exit 1; }
HBFUTILS_INI_SRC="$FIXTURE_DIR/o2simdigitizerworkflow_configuration.ini"
[[ -f "$HBFUTILS_INI_SRC" ]] || { echo "missing $HBFUTILS_INI_SRC" >&2; exit 1; }
grep -q "runNumber=$RUNNUMBER" "$HBFUTILS_INI_SRC" || { echo "RUNNUMBER=$RUNNUMBER does not match $HBFUTILS_INI_SRC" >&2; exit 1; }

mkdir -p "$REPLAY_DIR"
cp "$HBFUTILS_INI_SRC" "$REPLAY_DIR/o2simdigitizerworkflow_configuration.ini"
cd "$REPLAY_DIR"

echo "[replay_tracking] ITS replay: nthreads=$ITS_NTHREADS"
# Note: only the last command in a DPL pipe gets -b/--run; the upstream
# reader must stay a bare invocation so it dumps its workflow spec (via
# stdout) for merging instead of executing as its own standalone driver.
$TIME_CMD bash -c "
  o2-its-cluster-reader-workflow --with-mc \
    --input-dir '$FIXTURE_DIR' --its-cluster-infile o2clus_its.root \
    --shm-segment-size $SHM_SEGMENT_SIZE |
  o2-its-reco-workflow -b --run --clusters-from-upstream --tracking-mode sync \
    --cluster-rof-branch-only \
    --shm-segment-size $SHM_SEGMENT_SIZE \
    --condition-not-after $TIMESTAMP \
    --configKeyValues 'ITSCATrackerParam.nThreads=$ITS_NTHREADS' \
    --its-track-writer '--outfile o2trac_its.root'
" > its_replay.log 2> its_replay.time.log
grep -q "Processed [1-9][0-9]* TFs" its_replay.log || { echo "its-tracker processed 0 TFs, see its_replay.log" >&2; exit 1; }

echo "[replay_tracking] MFT replay"
$TIME_CMD bash -c "
  o2-mft-cluster-reader-workflow --with-mc \
    --input-dir '$FIXTURE_DIR' --mft-cluster-infile mftclusters.root \
    --shm-segment-size $SHM_SEGMENT_SIZE |
  o2-mft-reco-workflow -b --run --clusters-from-upstream \
    --cluster-rof-branch-only \
    --shm-segment-size $SHM_SEGMENT_SIZE \
    --condition-not-after $TIMESTAMP
" > mft_replay.log 2> mft_replay.time.log
[[ -s mfttracks.root ]] || { echo "mft replay produced no mfttracks.root, see mft_replay.log" >&2; exit 1; }

echo "[replay_tracking] outputs:"
ls -la o2trac_its.root mfttracks.root 2>&1 || true

#!/usr/bin/env bash
# Gate 3 Slice 3b: legacy ITS labelled reference. Replays the frozen legacy
# o2-its-reco-workflow (O2::ITStracking) ITS-only, in its ORDINARY Sync
# default-vertex mode -- i.e. NO static-diamond override
# (ITSCATrackerParam.useDiamond stays at its struct default, false), so the
# legacy tracker runs its own real per-event vertexer exactly as it would in
# production Sync reconstruction.
#
# This is a DIFFERENT leg from
# gate3-slice3-its-ca-validation/replay_tracking_its_legacy_diamond.sh
# (which forces ITSCATrackerParam.useDiamond=true to pair with the common-CA
# static-diamond leg under an "equivalent explicit constraint" framing). No
# file under gate3-slice3-its-ca-validation is read, written, or modified by
# this script.
#
# Per task scope: this is a LABELLED REFERENCE ONLY. Its real per-event
# vertexer is not comparable to the common-CA leg's static-diamond
# constraint (see replay_tracking_its_common_ca_nthreads.sh /
# gate3-slice3-its-ca-validation's own static-diamond legs). Do not treat
# any hash/track-count relationship between this leg and the common-CA legs
# as a pass/fail signal.
#
# Required variables:
#   FIXTURE_DIR   directory produced by gate0-baseline/generate_fixture.sh
#                 (read-only here).
#   REPLAY_DIR    output directory for this replay. Must not already exist
#                 as a non-empty directory.
#   TIMESTAMP     same fixed CCDB condition-not-after value used at fixture
#                 generation.
#   RUNNUMBER     same fixed run number used at fixture generation
#                 (sanity-checked against FIXTURE_DIR's HBFUtils ini).
#
# HBFUtils ini handling: identical workaround to
# gate0-baseline/replay_tracking.sh -- --hbfutils-config pointed at an
# absolute path reproducibly hung the ITS reco driver before it spawned any
# device; this script copies the ini into REPLAY_DIR under its default
# relative name and omits --hbfutils-config.
#
# Closed-stdin / pipefail robustness: gate3-slice3-its-ca-validation's
# README records that this exact FIXTURE_DIR/piped-ITS-reco invocation
# pattern hung indefinitely (no CPU progress, stuck after DPL init log
# lines) on 2 of 4 attempts under an unrelated REPLAY_DIR tree, root cause
# not identified. As a defensive measure (not a proven fix -- the hang was
# never reproduced with stdin already closed, so this is not claimed to
# eliminate it) this script explicitly redirects the whole piped
# invocation's stdin from /dev/null, on top of the existing
# "bash -o pipefail -c" wrapper that already ensures a failure in the
# upstream reader is not masked by the downstream tracker's exit code.
#
# Optional:
#   ITS_NTHREADS      ITSCATrackerParam.nThreads (default 1). No vertex
#                     constraint (diamond) config is set at all here -- only
#                     nThreads, so the tracker's real Sync vertexer default
#                     path is exercised unmodified.
#   SHM_SEGMENT_SIZE  bytes (default 4000000000)
#   TIME_CMD          wrapper prefixed to the reco invocation, e.g.
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

for bin in o2-its-cluster-reader-workflow o2-its-reco-workflow; do
  command -v "$bin" >/dev/null 2>&1 || { echo "missing $bin on PATH" >&2; exit 1; }
done
[[ -s "$FIXTURE_DIR/o2clus_its.root" ]] || { echo "missing or empty $FIXTURE_DIR/o2clus_its.root" >&2; exit 1; }
HBFUTILS_INI_SRC="$FIXTURE_DIR/o2simdigitizerworkflow_configuration.ini"
[[ -s "$HBFUTILS_INI_SRC" ]] || { echo "missing or empty $HBFUTILS_INI_SRC" >&2; exit 1; }
grep -q "runNumber=$RUNNUMBER" "$HBFUTILS_INI_SRC" || { echo "RUNNUMBER=$RUNNUMBER does not match $HBFUTILS_INI_SRC" >&2; exit 1; }

if [[ -d "$REPLAY_DIR" ]] && [[ -n "$(ls -A "$REPLAY_DIR" 2>/dev/null)" ]]; then
  echo "REPLAY_DIR '$REPLAY_DIR' already exists and is not empty; refusing to reuse it. Remove it or pick a fresh path." >&2
  exit 1
fi
mkdir -p "$REPLAY_DIR"
cp "$HBFUTILS_INI_SRC" "$REPLAY_DIR/o2simdigitizerworkflow_configuration.ini"
cd "$REPLAY_DIR"

echo "[replay_tracking_its_legacy_sync] ordinary Sync/default-vertex mode (no diamond override); nthreads=$ITS_NTHREADS"

# Only the last command in a DPL pipe gets -b/--run; the upstream reader
# must stay a bare invocation so it dumps its workflow spec (via stdout) for
# merging instead of executing as its own standalone driver. The inner
# "bash -o pipefail -c" starts a fresh shell that does NOT inherit this
# script's own "set -o pipefail" (each bash process gets its own default
# shell options), so it is set explicitly here too -- otherwise a failure in
# the upstream reader would be masked by the downstream tracker's exit code,
# which is the only one a plain pipe reports. stdin is explicitly closed
# (</dev/null) as a defensive measure against the known flaky-hang
# observation documented above.
$TIME_CMD bash -o pipefail -c "
  o2-its-cluster-reader-workflow --with-mc \
    --input-dir '$FIXTURE_DIR' --its-cluster-infile o2clus_its.root \
    --shm-segment-size $SHM_SEGMENT_SIZE |
  o2-its-reco-workflow -b --run --clusters-from-upstream --tracking-mode sync \
    --cluster-rof-branch-only \
    --shm-segment-size $SHM_SEGMENT_SIZE \
    --condition-not-after $TIMESTAMP \
    --configKeyValues 'ITSCATrackerParam.nThreads=$ITS_NTHREADS' \
    --its-track-writer '--outfile o2trac_its.root'
" < /dev/null > its_legacy_sync_replay.log 2> its_legacy_sync_replay.time.log
[[ -s o2trac_its.root ]] || { echo "its legacy-sync replay produced no o2trac_its.root, see its_legacy_sync_replay.log" >&2; exit 1; }
grep -q "Processed [1-9][0-9]* TFs" its_legacy_sync_replay.log || { echo "its-tracker processed 0 TFs, see its_legacy_sync_replay.log" >&2; exit 1; }

echo "[replay_tracking_its_legacy_sync] validating replay outputs"
ls -la o2trac_its.root

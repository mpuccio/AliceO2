#!/usr/bin/env bash
# Gate 3 Slice 3b: replay o2-its-ca-tracker-workflow with an explicit,
# recorded ITSCommonCATrackerParam.nThreads value, for thread/determinism
# characterization against the accepted Slice 3 canonical baseline.
#
# This is a SIBLING to
# gate3-slice3-its-ca-validation/replay_tracking_its_common_ca.sh, not a
# replacement or a modification of it -- that file has no nThreads knob
# (ITSCommonCATrackerParam.nThreads did not exist when Slice 3 ran; it was
# added afterwards in commit 62e2c75bc8 "ITSMFTTracking: add explicit
# ITSCommonCATrackerParam.nThreads"). No file under gate3-slice3-its-ca-validation
# is read, written, or modified by this script. It reads the same durable
# Gate 0 fixture (read-only) that gate0-baseline/generate_fixture.sh
# produced.
#
# Static-diamond compatibility note: since commit b47271fb13
# ("ITSMFTTracking: derive static-diamond timestamp from the real per-ROF
# interval envelope"), the diamond vertex's per-ROF timestamp compatibility
# is derived automatically inside TrackerTraits from each ROF's own real
# interval envelope -- there is no separate "diamond timestamp" CLI/config
# knob to set here. TIMESTAMP below is only the ordinary
# --condition-not-after CCDB gate, identical in meaning to every other
# replay script in this validation family.
#
# Required variables:
#   FIXTURE_DIR     directory produced by gate0-baseline/generate_fixture.sh
#                   (read-only here).
#   REPLAY_DIR      output directory for this replay. Must not already exist
#                   as a non-empty directory -- refuses to run into a dirty
#                   REPLAY_DIR rather than silently mixing runs.
#   TIMESTAMP       same fixed CCDB condition-not-after value used at
#                   fixture generation.
#   RUNNUMBER       same fixed run number used at fixture generation
#                   (sanity-checked against FIXTURE_DIR's HBFUtils ini).
#   ITS_CA_NTHREADS explicit ITSCommonCATrackerParam.nThreads value for this
#                   replay (no default -- must be passed explicitly so every
#                   invocation in this characterization records exactly
#                   which value it used).
#
# HBFUtils ini handling: identical workaround to the Slice 3 script --
# --hbfutils-config pointed at an absolute path reproducibly hung the ITS
# reco driver before it spawned any device when reading clusters through a
# piped cluster-reader-workflow; this script copies the ini into REPLAY_DIR
# under its default relative name and omits --hbfutils-config.
#
# Optional:
#   DIAMOND_POS_X/Y/Z   ITSCommonCATrackerParam.diamondPos[0..2] (cm),
#                       default 0/0/0.
#   PV_RES              ITSCommonCATrackerParam.pvRes (cm), default 0.05.
#   SHM_SEGMENT_SIZE    bytes (default 4000000000)
#   TIME_CMD            wrapper prefixed to the reco invocation, e.g.
#                        "/usr/bin/time -l" to capture wall time + peak RSS.
#                        Left empty by default (no wrapping).

set -euo pipefail

: "${FIXTURE_DIR:?set FIXTURE_DIR}"
: "${REPLAY_DIR:?set REPLAY_DIR}"
: "${TIMESTAMP:?set TIMESTAMP}"
: "${RUNNUMBER:?set RUNNUMBER}"
: "${ITS_CA_NTHREADS:?set ITS_CA_NTHREADS explicitly}"
DIAMOND_POS_X="${DIAMOND_POS_X:-0}"
DIAMOND_POS_Y="${DIAMOND_POS_Y:-0}"
DIAMOND_POS_Z="${DIAMOND_POS_Z:-0}"
PV_RES="${PV_RES:-0.05}"
SHM_SEGMENT_SIZE="${SHM_SEGMENT_SIZE:-4000000000}"
TIME_CMD="${TIME_CMD:-}"

for bin in o2-its-cluster-reader-workflow o2-its-ca-tracker-workflow; do
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

CKV="ITSCommonCATrackerParam.useDiamond=true;ITSCommonCATrackerParam.diamondPos[0]=$DIAMOND_POS_X;ITSCommonCATrackerParam.diamondPos[1]=$DIAMOND_POS_Y;ITSCommonCATrackerParam.diamondPos[2]=$DIAMOND_POS_Z;ITSCommonCATrackerParam.pvRes=$PV_RES;ITSCommonCATrackerParam.nThreads=$ITS_CA_NTHREADS"
echo "[replay_tracking_its_common_ca_nthreads] configKeyValues: $CKV"

# Only the last command in a DPL pipe gets -b/--run; the upstream reader
# must stay a bare invocation so it dumps its workflow spec (via stdout) for
# merging instead of executing as its own standalone driver. The inner
# "bash -c" starts a fresh shell that does NOT inherit this script's own
# "set -o pipefail" (each bash process gets its own default shell options),
# so it is set explicitly here too.
$TIME_CMD bash -o pipefail -c "
  o2-its-cluster-reader-workflow --with-mc \
    --input-dir '$FIXTURE_DIR' --its-cluster-infile o2clus_its.root \
    --shm-segment-size $SHM_SEGMENT_SIZE |
  o2-its-ca-tracker-workflow -b --run --tracking-mode sync \
    --shm-segment-size $SHM_SEGMENT_SIZE \
    --condition-not-after $TIMESTAMP \
    --configKeyValues '$CKV'
" > its_ca_replay.log 2> its_ca_replay.time.log
[[ -s o2trac_its_ca.root ]] || { echo "its common-CA replay produced no o2trac_its_ca.root, see its_ca_replay.log" >&2; exit 1; }

echo "[replay_tracking_its_common_ca_nthreads] validating replay outputs"
ls -la o2trac_its_ca.root

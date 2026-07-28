#!/usr/bin/env bash
# Gate 3 Slice 3: replay o2-its-ca-tracker-workflow (the opt-in ITS
# common-CA workflow built in Gate 3 Slice 2, o2::its::ca /
# o2::itsmft::tracking::ITSMFTTrackingInterface<7>) from the same fixed
# cluster-level fixture used by the accepted Gate 0 baseline.
#
# This is a SIBLING to gate0-baseline/replay_tracking.sh and
# gate0-baseline/replay_tracking_common_ca.sh (the MFT common-CA
# characterization script), not a replacement for either: no file in
# gate0-baseline is read, written, or modified by this script. It reads the
# same durable fixture (read-only) that gate0-baseline/generate_fixture.sh
# produced.
#
# o2-its-ca-tracker-workflow has no embedded clusterizer or cluster-writer
# device (identical InputSpecs to the legacy ITS TrackerSpec: ITS/
# COMPCLUSTERS, ITS/PATTERNS, ITS/CLUSTERSROF, ITS/CLUSTERSMCTR with MC), so
# it pipes directly after o2-its-cluster-reader-workflow with no
# --clusters-from-upstream/--cluster-rof-branch-only flag -- neither exists
# on this workflow.
#
# Vertex/beam constraint: this workflow's ConfigPreflight
# (ITSCAWorkflow/ConfigPreflight.h) fatals unless
# ITSCommonCATrackerParam.useDiamond=true, since the common core has no real
# per-event vertexing capability for ITS. DIAMOND_POS_X/Y/Z and PV_RES below
# are passed explicitly (not left at their struct defaults) so the exact
# static-diamond constraint used is recorded here rather than implied. This
# is characterization against a static synthetic vertex constraint, NOT a
# claim of equivalence to the legacy real-vertexer Sync default or a CCDB
# MeanVertex beam-position override -- see replay_tracking_its_legacy_diamond.sh
# for the paired legacy leg using the identical explicit values via
# ITSCATrackerParam (legacy namespace, same field names/semantics).
#
# Required variables:
#   FIXTURE_DIR   directory produced by gate0-baseline/generate_fixture.sh
#                 (read-only here).
#   REPLAY_DIR    output directory for this replay. Must not already exist
#                 as a non-empty directory -- refuses to run into a dirty
#                 REPLAY_DIR rather than silently mixing runs. Remove it, or
#                 point at a fresh path, to re-replay.
#   TIMESTAMP     same fixed CCDB condition-not-after value used at fixture
#                 generation.
#   RUNNUMBER     same fixed run number used at fixture generation
#                 (sanity-checked against FIXTURE_DIR's HBFUtils ini).
#
# HBFUtils ini handling: identical workaround to gate0-baseline's replay
# scripts -- --hbfutils-config pointed at an absolute path reproducibly hung
# the ITS reco driver before it spawned any device when reading clusters
# through a piped cluster-reader-workflow; this script copies the ini into
# REPLAY_DIR under its default relative name and omits --hbfutils-config.
#
# Optional:
#   DIAMOND_POS_X/Y/Z   ITSCommonCATrackerParam.diamondPos[0..2] (cm),
#                       default 0/0/0.
#   PV_RES              ITSCommonCATrackerParam.pvRes (cm), default 0.05.
#   EXTRA_CONFIGKEYVALUES  additional ';'-separated ITSCommonCATrackerParam.*
#                       (or other namespace) key=value pairs appended after
#                       the mandatory diamond-constraint keys. Left empty by
#                       default.
#   SHM_SEGMENT_SIZE    bytes (default 4000000000)
#   TIME_CMD            wrapper prefixed to the reco invocation, e.g.
#                        "/usr/bin/time -l" to capture wall time + peak RSS.
#                        Left empty by default (no wrapping).

set -euo pipefail

: "${FIXTURE_DIR:?set FIXTURE_DIR}"
: "${REPLAY_DIR:?set REPLAY_DIR}"
: "${TIMESTAMP:?set TIMESTAMP}"
: "${RUNNUMBER:?set RUNNUMBER}"
DIAMOND_POS_X="${DIAMOND_POS_X:-0}"
DIAMOND_POS_Y="${DIAMOND_POS_Y:-0}"
DIAMOND_POS_Z="${DIAMOND_POS_Z:-0}"
PV_RES="${PV_RES:-0.05}"
EXTRA_CONFIGKEYVALUES="${EXTRA_CONFIGKEYVALUES:-}"
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

CKV="ITSCommonCATrackerParam.useDiamond=true;ITSCommonCATrackerParam.diamondPos[0]=$DIAMOND_POS_X;ITSCommonCATrackerParam.diamondPos[1]=$DIAMOND_POS_Y;ITSCommonCATrackerParam.diamondPos[2]=$DIAMOND_POS_Z;ITSCommonCATrackerParam.pvRes=$PV_RES"
if [[ -n "$EXTRA_CONFIGKEYVALUES" ]]; then
  CKV="$CKV;$EXTRA_CONFIGKEYVALUES"
fi
echo "[replay_tracking_its_common_ca] configKeyValues: $CKV"

# Only the last command in a DPL pipe gets -b/--run; the upstream reader
# must stay a bare invocation so it dumps its workflow spec (via stdout) for
# merging instead of executing as its own standalone driver. The inner
# "bash -c" starts a fresh shell that does NOT inherit this script's
# "set -o pipefail" (each bash process gets its own default shell options),
# so it is set explicitly here too -- otherwise a failure in the upstream
# reader would be masked by the downstream tracker's exit code, which is the
# only one a plain pipe reports.
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

echo "[replay_tracking_its_common_ca] validating replay outputs"
ls -la o2trac_its_ca.root

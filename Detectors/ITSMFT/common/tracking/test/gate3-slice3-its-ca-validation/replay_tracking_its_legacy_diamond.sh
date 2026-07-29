#!/usr/bin/env bash
# Gate 3 Slice 3: replay the frozen legacy o2-its-reco-workflow (O2::ITStracking)
# ITS-only, from the same fixed cluster-level fixture, with an explicit
# static-diamond vertex/beam constraint (ITSCATrackerParam.useDiamond=true)
# instead of its real per-event Sync vertexer default.
#
# This is a SIBLING to gate0-baseline/replay_tracking.sh, not a replacement:
# no file in gate0-baseline is read, written, or modified by this script.
# Unlike gate0-baseline/replay_tracking.sh, this script runs ONLY the ITS
# leg (no MFT) and does not require o2-mft-reco-workflow to be built --
# Gate 3 Slice 3's build scope is legacy o2-its-reco-workflow, not legacy
# MFT tracking. It exists purely to pair diagnostically with
# replay_tracking_its_common_ca.sh, passing the same diamondPos/pvRes CLI
# values (ITSCATrackerParam is the legacy namespace equivalent of
# ITSCommonCATrackerParam in field names only) -- NOT an equivalent
# vertex/beam-constraint configuration: this leg's diamond vertex carries no
# timestamp and yields zero tracks (see characterization_summary.md Sec.3),
# while the common-CA leg derives a separate, common-only, ROF-local
# timestamped diamond. Never to claim the common-CA leg reproduces the
# legacy real-vertexer/CCDB MeanVertex-override Sync default, nor that the
# two legs' constraints are equivalent to each other.
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
# HBFUtils ini handling: identical workaround to gate0-baseline/replay_tracking.sh.
#
# Optional:
#   DIAMOND_POS_X/Y/Z   ITSCATrackerParam.diamondPos[0..2] (cm), default 0/0/0.
#                       Must match replay_tracking_its_common_ca.sh's
#                       DIAMOND_POS_X/Y/Z for the two legs to be comparable.
#   PV_RES              ITSCATrackerParam.pvRes (cm), default 0.05. Must
#                       match replay_tracking_its_common_ca.sh's PV_RES.
#   ITS_NTHREADS        ITSCATrackerParam.nThreads (default 1).
#   SHM_SEGMENT_SIZE    bytes (default 4000000000)
#   TIME_CMD            wrapper prefixed to the reco invocation, e.g.
#                        "/usr/bin/time -l". Left empty by default.

set -euo pipefail

: "${FIXTURE_DIR:?set FIXTURE_DIR}"
: "${REPLAY_DIR:?set REPLAY_DIR}"
: "${TIMESTAMP:?set TIMESTAMP}"
: "${RUNNUMBER:?set RUNNUMBER}"
DIAMOND_POS_X="${DIAMOND_POS_X:-0}"
DIAMOND_POS_Y="${DIAMOND_POS_Y:-0}"
DIAMOND_POS_Z="${DIAMOND_POS_Z:-0}"
PV_RES="${PV_RES:-0.05}"
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

CKV="ITSCATrackerParam.useDiamond=true;ITSCATrackerParam.diamondPos[0]=$DIAMOND_POS_X;ITSCATrackerParam.diamondPos[1]=$DIAMOND_POS_Y;ITSCATrackerParam.diamondPos[2]=$DIAMOND_POS_Z;ITSCATrackerParam.pvRes=$PV_RES;ITSCATrackerParam.nThreads=$ITS_NTHREADS"
echo "[replay_tracking_its_legacy_diamond] configKeyValues: $CKV"

# See gate0-baseline/replay_tracking.sh for why the inner "bash -o pipefail -c"
# is required (a nested bash -c does not inherit this script's own
# "set -o pipefail").
$TIME_CMD bash -o pipefail -c "
  o2-its-cluster-reader-workflow --with-mc \
    --input-dir '$FIXTURE_DIR' --its-cluster-infile o2clus_its.root \
    --shm-segment-size $SHM_SEGMENT_SIZE |
  o2-its-reco-workflow -b --run --clusters-from-upstream --tracking-mode sync \
    --cluster-rof-branch-only \
    --shm-segment-size $SHM_SEGMENT_SIZE \
    --condition-not-after $TIMESTAMP \
    --configKeyValues '$CKV' \
    --its-track-writer '--outfile o2trac_its.root'
" > its_legacy_diamond_replay.log 2> its_legacy_diamond_replay.time.log
[[ -s o2trac_its.root ]] || { echo "its legacy-diamond replay produced no o2trac_its.root, see its_legacy_diamond_replay.log" >&2; exit 1; }
grep -q "Processed [1-9][0-9]* TFs" its_legacy_diamond_replay.log || { echo "its-tracker processed 0 TFs, see its_legacy_diamond_replay.log" >&2; exit 1; }

echo "[replay_tracking_its_legacy_diamond] validating replay outputs"
ls -la o2trac_its.root

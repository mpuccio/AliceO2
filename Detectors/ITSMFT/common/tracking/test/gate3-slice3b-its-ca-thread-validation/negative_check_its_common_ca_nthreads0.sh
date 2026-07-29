#!/usr/bin/env bash
# Gate 3 Slice 3b: negative check for ITSCommonCATrackerParam.nThreads=0.
#
# Unlike gate3-slice3-its-ca-validation/negative_checks_its_common_ca.sh's
# two cases (both rejected inside ConfigPreflight, at
# defineDataProcessing()-time, before any device is constructed), the
# nThreads<=0 check added in commit 62e2c75bc8 ("ITSMFTTracking: add
# explicit ITSCommonCATrackerParam.nThreads") lives in
# ITSMFTTrackingInterface::initialiseTracker(), which is only reached from
# CATrackerDPL::updateTimeDependentParams() on the FIRST real
# ProcessingContext the device receives (see
# Detectors/ITSMFT/ITS/workflow-ca/src/CATrackerSpec.cxx:184-189) -- i.e.
# this fatal is NOT observable from a standalone invocation with no piped
# input (confirmed by source inspection before writing this script: no
# amount of waiting fires it without a TF actually arriving). This script
# therefore runs the FULL piped
# o2-its-cluster-reader-workflow | o2-its-ca-tracker-workflow, with an
# otherwise-valid explicit static-diamond configuration (useDiamond=true,
# diamondPos=(0,0,0), pvRes=0.05 -- the same values used everywhere else in
# this validation) plus ITSCommonCATrackerParam.nThreads=0, and asserts a
# fatal fires on/after the first TF, before any o2trac_its_ca.root is
# created.
#
# Reads the same durable, read-only Gate 0 fixture as every other script in
# this validation; writes only into OUT_DIR (a fresh directory, never the
# fixture and never a baseline REPLAY_DIR). No baseline output is at risk:
# OUT_DIR is required fresh/empty like every other REPLAY_DIR in this
# validation.
#
# Output-file caveat (found empirically while writing this script, unlike
# gate3-slice3-its-ca-validation's two ConfigPreflight-level negative checks
# where no device is ever constructed): the its-ca-track-writer DPL device
# opens/creates o2trac_its_ca.root as part of its own ordinary startup,
# independently of whether the tracker device ever successfully produces a
# TF for it to write -- so a FILE existing here is expected, not a failure
# signal. What this script actually asserts is that the "o2sim" TTree inside
# it has ZERO entries (confirmed by direct inspection before writing this
# assertion: the tracker's LOGP(fatal) aborts the process during its FIRST
# run() call, upstream of ever calling the writer's output ports), i.e. no
# track data was ever produced or written -- not that no file appears on
# disk at all.
#
# A `timeout` wrapper is used (not just LOGP(fatal)'s own SIGABRT) as a
# defensive bound against gate3-slice3-its-ca-validation/README.md's
# documented flaky-hang observation for this same
# cluster-reader-workflow-piped-into-ITS-tracker invocation shape (root
# cause not identified there either) -- if the expected fatal does not fire
# promptly this script reports a timeout failure rather than hanging
# indefinitely.
#
# Required variables:
#   FIXTURE_DIR   directory produced by gate0-baseline/generate_fixture.sh
#                 (read-only here).
#   OUT_DIR       output directory for this check. Must not already exist
#                 as a non-empty directory.
#   TIMESTAMP     same fixed CCDB condition-not-after value used elsewhere
#                 in this validation.
#   RUNNUMBER     same fixed run number used at fixture generation
#                 (sanity-checked against FIXTURE_DIR's HBFUtils ini).
#
# Optional:
#   SHM_SEGMENT_SIZE  bytes (default 4000000000)
#   TIMEOUT_SECONDS   seconds before giving up and reporting failure
#                     (default 180). Requires GNU `timeout` (or macOS
#                     Homebrew coreutils `gtimeout`) on PATH.
#
# Exit status: 0 if the fatal fired with the expected message before any
# output file appeared; non-zero otherwise (message on stderr identifying
# which assertion failed).

set -uo pipefail

: "${FIXTURE_DIR:?set FIXTURE_DIR}"
: "${OUT_DIR:?set OUT_DIR}"
: "${TIMESTAMP:?set TIMESTAMP}"
: "${RUNNUMBER:?set RUNNUMBER}"
SHM_SEGMENT_SIZE="${SHM_SEGMENT_SIZE:-4000000000}"
TIMEOUT_SECONDS="${TIMEOUT_SECONDS:-180}"

TIMEOUT_BIN=""
for cand in timeout gtimeout; do
  if command -v "$cand" >/dev/null 2>&1; then
    TIMEOUT_BIN="$cand"
    break
  fi
done
[[ -n "$TIMEOUT_BIN" ]] || { echo "missing 'timeout'/'gtimeout' on PATH" >&2; exit 1; }

for bin in o2-its-cluster-reader-workflow o2-its-ca-tracker-workflow root; do
  command -v "$bin" >/dev/null 2>&1 || { echo "missing $bin on PATH" >&2; exit 1; }
done
[[ -s "$FIXTURE_DIR/o2clus_its.root" ]] || { echo "missing or empty $FIXTURE_DIR/o2clus_its.root" >&2; exit 1; }
HBFUTILS_INI_SRC="$FIXTURE_DIR/o2simdigitizerworkflow_configuration.ini"
[[ -s "$HBFUTILS_INI_SRC" ]] || { echo "missing or empty $HBFUTILS_INI_SRC" >&2; exit 1; }
grep -q "runNumber=$RUNNUMBER" "$HBFUTILS_INI_SRC" || { echo "RUNNUMBER=$RUNNUMBER does not match $HBFUTILS_INI_SRC" >&2; exit 1; }

if [[ -d "$OUT_DIR" ]] && [[ -n "$(ls -A "$OUT_DIR" 2>/dev/null)" ]]; then
  echo "OUT_DIR '$OUT_DIR' already exists and is not empty; refusing to reuse it." >&2
  exit 1
fi
mkdir -p "$OUT_DIR"
cp "$HBFUTILS_INI_SRC" "$OUT_DIR/o2simdigitizerworkflow_configuration.ini"
cd "$OUT_DIR"

EXPECTED_SUBSTRING="ITS CA tracker requires ITSCommonCATrackerParam.nThreads > 0, got 0"
CKV="ITSCommonCATrackerParam.useDiamond=true;ITSCommonCATrackerParam.diamondPos[0]=0;ITSCommonCATrackerParam.diamondPos[1]=0;ITSCommonCATrackerParam.diamondPos[2]=0;ITSCommonCATrackerParam.pvRes=0.05;ITSCommonCATrackerParam.nThreads=0"
echo "[negative_check_its_common_ca_nthreads0] configKeyValues: $CKV"

"$TIMEOUT_BIN" "${TIMEOUT_SECONDS}s" bash -o pipefail -c "
  o2-its-cluster-reader-workflow --with-mc \
    --input-dir '$FIXTURE_DIR' --its-cluster-infile o2clus_its.root \
    --shm-segment-size $SHM_SEGMENT_SIZE |
  o2-its-ca-tracker-workflow -b --run --tracking-mode sync \
    --shm-segment-size $SHM_SEGMENT_SIZE \
    --condition-not-after $TIMESTAMP \
    --configKeyValues '$CKV'
" < /dev/null > run.log 2>&1
status=$?

echo "[negative_check_its_common_ca_nthreads0] exit=$status"
overall_status=0

if [[ "$status" -eq 124 ]]; then
  echo "  FAIL: timed out after ${TIMEOUT_SECONDS}s waiting for the expected fatal (see run.log)" >&2
  overall_status=1
elif [[ "$status" -eq 0 ]]; then
  echo "  FAIL: expected a non-zero (fatal) exit, got 0" >&2
  overall_status=1
fi

if ! grep -qF "$EXPECTED_SUBSTRING" run.log; then
  echo "  FAIL: expected log to contain: $EXPECTED_SUBSTRING" >&2
  overall_status=1
fi

if [[ -e "$OUT_DIR/o2trac_its_ca.root" ]]; then
  # The writer device creates this file as part of its own ordinary
  # startup regardless of whether the tracker ever succeeds (see header
  # comment) -- what matters is that its "o2sim" tree has zero entries,
  # i.e. no track data was ever written to it.
  nentries="$(ROOT_HIST=0 root -l -b -q -e "TFile f(\"$OUT_DIR/o2trac_its_ca.root\"); TTree* t=(TTree*)f.Get(\"o2sim\"); std::cout << (t ? t->GetEntries() : -1) << std::endl;" 2>/dev/null | tail -1)"
  if [[ "$nentries" != "0" ]]; then
    echo "  FAIL: o2trac_its_ca.root's o2sim tree has $nentries entries (expected 0 -- no track data should have been written)" >&2
    overall_status=1
  fi
fi

if [[ "$overall_status" -eq 0 ]]; then
  echo "  OK: non-zero exit, expected fatal message present, no track data written (output file, if any, has zero entries)"
fi

exit "$overall_status"

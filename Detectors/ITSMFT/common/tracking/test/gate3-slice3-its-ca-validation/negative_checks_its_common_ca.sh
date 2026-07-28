#!/usr/bin/env bash
# Gate 3 Slice 3: negative preflight checks for o2-its-ca-tracker-workflow.
#
# Both checks below invoke o2-its-ca-tracker-workflow standalone (no piped
# cluster-reader input): ConfigPreflight.h's three checks run inside
# defineDataProcessing(), i.e. at workflow-spec construction time, strictly
# before any DataProcessorSpec/device is built or any input is read -- so a
# rejected configuration is observable with no input pipe at all. This was
# confirmed interactively before writing this script: piping input changes
# nothing about where/when the fatal fires.
#
# Check 1 (no-diamond): default configuration (ITSCommonCATrackerParam.useDiamond
# stays false) must fatal via requireDiamondVertexConstraintOrFatal(),
# before any device runs and before o2trac_its_ca.root is created.
#
# Check 2 (legacy-param-override): a --configKeyValues string carrying a
# legacy 'ITSCATrackerParam.*' key must fatal via
# applyConfigKeyValuesOrFatal()'s RejectedLegacyITSNamespace path, BEFORE
# o2::conf::ConfigurableParam::updateFromString() is ever called on that
# string (see ConfigPreflight.cxx) -- i.e. the legacy override is never
# even applied, let alone allowed to reach tracking/output.
#
# Required variables:
#   TIMESTAMP   same fixed CCDB condition-not-after value used elsewhere in
#               this validation slice (only used to keep the invocation
#               realistic; the preflight fatals before any CCDB access).
#   OUT_DIR     directory to run both checks in (created fresh; must not
#               already exist as non-empty). Two subdirectories
#               (no-diamond/, legacy-param-override/) are created inside it.
#
# Exit status: 0 if both checks observed the expected fatal + no-output
# behavior, non-zero (with a message identifying which check failed) if
# either check's exit code, log content, or output-file absence did not
# match what is asserted above.

set -euo pipefail

: "${TIMESTAMP:?set TIMESTAMP}"
: "${OUT_DIR:?set OUT_DIR}"

command -v o2-its-ca-tracker-workflow >/dev/null 2>&1 || { echo "missing o2-its-ca-tracker-workflow on PATH" >&2; exit 1; }

if [[ -d "$OUT_DIR" ]] && [[ -n "$(ls -A "$OUT_DIR" 2>/dev/null)" ]]; then
  echo "OUT_DIR '$OUT_DIR' already exists and is not empty; refusing to reuse it." >&2
  exit 1
fi
mkdir -p "$OUT_DIR"

overall_status=0

run_negative_case() {
  local name="$1" expected_substring="$2"; shift 2
  local dir="$OUT_DIR/$name"
  mkdir -p "$dir"
  local log="$dir/run.log"
  local status=0
  ( cd "$dir" && o2-its-ca-tracker-workflow -b --run --tracking-mode sync \
      --condition-not-after "$TIMESTAMP" "$@" > run.log 2>&1 < /dev/null ) || status=$?

  echo "[negative_checks_its_common_ca] case=$name exit=$status"
  if [[ "$status" -eq 0 ]]; then
    echo "  FAIL: expected a non-zero (fatal) exit, got 0" >&2
    overall_status=1
  fi
  if ! grep -qF "$expected_substring" "$log"; then
    echo "  FAIL: expected log to contain: $expected_substring" >&2
    overall_status=1
  fi
  if [[ -e "$dir/o2trac_its_ca.root" ]]; then
    echo "  FAIL: o2trac_its_ca.root was created; negative check must produce no usable output" >&2
    overall_status=1
  fi
  if [[ "$status" -ne 0 ]] && grep -qF "$expected_substring" "$log" && [[ ! -e "$dir/o2trac_its_ca.root" ]]; then
    echo "  OK: non-zero exit, expected fatal message present, no output file"
  fi
}

run_negative_case "no-diamond" \
  "requires ITSCommonCATrackerParam.useDiamond=true"

run_negative_case "legacy-param-override" \
  "rejects legacy 'ITSCATrackerParam' --configKeyValues override" \
  --configKeyValues 'ITSCATrackerParam.useDiamond=true'

exit "$overall_status"

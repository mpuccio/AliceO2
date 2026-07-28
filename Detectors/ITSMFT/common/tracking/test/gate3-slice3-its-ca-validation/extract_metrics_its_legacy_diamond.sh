#!/usr/bin/env bash
# Gate 3 Slice 3: extract metrics for the legacy ITS-only diamond-constraint
# replay (replay_tracking_its_legacy_diamond.sh's o2trac_its.root), reusing
# gate0-baseline/extract_metrics.C's extractITS() function UNMODIFIED and
# UNCOPIED (loaded via ROOT's ".L", not duplicated into this directory) --
# same counting/denominator convention as the accepted Gate 0 baseline, no
# new extraction logic.
#
# This exists only because extract_metrics.C's top-level extract_metrics()
# driver unconditionally also calls extractMFT() (which needs mfttracks.root
# in replayDir); the paired legacy leg in this Gate 3 Slice 3 validation is
# ITS-only (see replay_tracking_its_legacy_diamond.sh's own header for why),
# so this wrapper calls extractITS() directly instead of going through
# extract_metrics()'s ITS+MFT-combined entry point. gate0-baseline's own
# files (extract_metrics.C, baseline_summary.md, manifest.json) are neither
# read-written-to nor modified beyond a read-only ".L" load.
#
# Usage:
#   FIXTURE_DIR=<fixture> REPLAY_DIR=<legacy diamond replay dir> \
#     OUT_JSON=<out.json> ./extract_metrics_its_legacy_diamond.sh

set -euo pipefail

: "${FIXTURE_DIR:?set FIXTURE_DIR}"
: "${REPLAY_DIR:?set REPLAY_DIR}"
: "${OUT_JSON:?set OUT_JSON}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXTRACT_METRICS_C="$SCRIPT_DIR/../gate0-baseline/extract_metrics.C"
[[ -s "$EXTRACT_METRICS_C" ]] || { echo "missing $EXTRACT_METRICS_C" >&2; exit 1; }
[[ -s "$REPLAY_DIR/o2trac_its.root" ]] || { echo "missing or empty $REPLAY_DIR/o2trac_its.root" >&2; exit 1; }

command -v root >/dev/null 2>&1 || { echo "missing root on PATH" >&2; exit 1; }

root -l -b -q \
  -e ".L $EXTRACT_METRICS_C" \
  -e "o2::base::GeometryManager::loadGeometry(std::string(\"$FIXTURE_DIR\") + \"/o2sim\");" \
  -e "std::ofstream out(\"$OUT_JSON\"); out << \"{\\n\"; extractITS(\"$FIXTURE_DIR\", \"$REPLAY_DIR\", out); out << \"\\n}\\n\"; out.close(); std::cout << \"[extract_metrics_its_legacy_diamond] wrote $OUT_JSON\" << std::endl;"

[[ -s "$OUT_JSON" ]] || { echo "extraction produced no output: $OUT_JSON" >&2; exit 1; }

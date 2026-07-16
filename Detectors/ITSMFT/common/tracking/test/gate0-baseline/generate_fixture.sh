#!/usr/bin/env bash
# Gate 0 baseline: generate an ITS+MFT tracker-input fixture (sim -> digi -> initial reco).
#
# Not wired into CTest. Produces files outside the repository; nothing here
# writes into the git worktree. See README.md for the full protocol.
#
# Required environment (not hardcoded, since this is machine-specific):
#   PATH/LD_LIBRARY_PATH/DYLD_LIBRARY_PATH must already resolve o2-sim,
#   o2-sim-digitizer-workflow, o2-its-reco-workflow, o2-mft-reco-workflow,
#   o2-grp-simgrp-tool built from the commit under test.
#
# Required variables:
#   FIXTURE_DIR   output directory. Must not already exist as a non-empty
#                 directory -- this script refuses to run into a dirty
#                 FIXTURE_DIR rather than silently mixing outputs from a
#                 previous run with this one. Remove it first, or point at
#                 a fresh path, to regenerate.
#   RUNNUMBER     fixed ALICE run number
#   SEED          fixed o2-sim seed
#   TIMESTAMP     fixed CCDB condition-not-after value, ms since epoch
#   NEVENTS       number of pp events to generate
# Optional:
#   SIM_WORKERS   o2-sim -j (default 1)
#   SHM_SEGMENT_SIZE  bytes (default 4000000000)
#   INTERACTION_RATE_HZ  (default 400000)

set -euo pipefail

: "${FIXTURE_DIR:?set FIXTURE_DIR}"
: "${RUNNUMBER:?set RUNNUMBER}"
: "${SEED:?set SEED}"
: "${TIMESTAMP:?set TIMESTAMP}"
: "${NEVENTS:?set NEVENTS}"
SIM_WORKERS="${SIM_WORKERS:-1}"
SHM_SEGMENT_SIZE="${SHM_SEGMENT_SIZE:-4000000000}"
INTERACTION_RATE_HZ="${INTERACTION_RATE_HZ:-400000}"

for bin in o2-sim o2-sim-digitizer-workflow o2-its-reco-workflow o2-mft-reco-workflow; do
  command -v "$bin" >/dev/null 2>&1 || { echo "missing $bin on PATH" >&2; exit 1; }
done

if [[ -d "$FIXTURE_DIR" ]] && [[ -n "$(ls -A "$FIXTURE_DIR" 2>/dev/null)" ]]; then
  echo "FIXTURE_DIR '$FIXTURE_DIR' already exists and is not empty; refusing to reuse it. Remove it or pick a fresh path." >&2
  exit 1
fi
mkdir -p "$FIXTURE_DIR"
cd "$FIXTURE_DIR"

echo "[generate_fixture] sim: n=$NEVENTS run=$RUNNUMBER seed=$SEED j=$SIM_WORKERS"
o2-sim -j "$SIM_WORKERS" -n "$NEVENTS" -g pythia8pp -e TGeant4 -m PIPE ITS MFT \
  --run "$RUNNUMBER" --seed "$SEED" --field -5 -o o2sim \
  > sim.log 2>&1

echo "[generate_fixture] digitization (ITS,MFT only)"
o2-sim-digitizer-workflow -b --run --onlyDet ITS,MFT \
  --interactionRate "$INTERACTION_RATE_HZ" \
  --configKeyValues "HBFUtils.runNumber=${RUNNUMBER}" \
  --shm-segment-size "$SHM_SEGMENT_SIZE" \
  --condition-not-after "$TIMESTAMP" \
  > digi.log 2>&1

echo "[generate_fixture] initial ITS reco (clusterizer+tracker+vertexer, produces fixture clusters)"
/usr/bin/time -l o2-its-reco-workflow -b --run --tracking-mode sync \
  --shm-segment-size "$SHM_SEGMENT_SIZE" \
  --condition-not-after "$TIMESTAMP" \
  --configKeyValues "ITSCATrackerParam.nThreads=1" \
  --resources-monitoring 1 --resources-monitoring-file its_reco_initial.resources.json \
  > its_reco_initial.log 2> its_reco_initial.time.log

echo "[generate_fixture] initial MFT reco (clusterizer+CA tracker, produces fixture clusters)"
/usr/bin/time -l o2-mft-reco-workflow -b --run \
  --shm-segment-size "$SHM_SEGMENT_SIZE" \
  --condition-not-after "$TIMESTAMP" \
  --resources-monitoring 1 --resources-monitoring-file mft_reco_initial.resources.json \
  > mft_reco_initial.log 2> mft_reco_initial.time.log

echo "[generate_fixture] validating fixture outputs"
for f in o2clus_its.root mftclusters.root o2trac_its.root mfttracks.root \
         o2sim_Kine.root o2sim_geometry.root o2sim_geometry-aligned.root \
         o2simdigitizerworkflow_configuration.ini; do
  [[ -s "$f" ]] || { echo "missing or empty required fixture output: $f" >&2; exit 1; }
done
ls -la o2clus_its.root mftclusters.root o2trac_its.root mfttracks.root o2sim_Kine.root

#!/bin/zsh
#
# Preflight checks for batch ALICE simulation (o2-sim / TGeant4): AliEn
# token status, coherent aliBuild provenance for o2-sim/ROOT/Geant4/
# Geant4VMC/VMC, and every G4*DATA dataset variable Geant4 needs -- derived
# MECHANICALLY from `geant4-config --datasets`, never hand-transcribed. See
# SKILL.md for usage.
#
# This script only inspects and reports; it does not run a simulation. Use
# --print-exports to get eval-able G4*DATA assignments for a real run.

set -e

SCRIPT_DIR="${0:A:h}"
A2E_PROG="simulation-preflight.zsh"
source "${SCRIPT_DIR}/lib/env-common.zsh"

fail() { print -u2 -- "${A2E_PROG}: $*"; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  simulation-preflight.zsh [--package NAME] [--print-exports] [--strict] [--help]

Options:
  --package NAME     O2_PACKAGE for the runner (alias or exact resolved
                      name). Default: inherits O2_PACKAGE from the
                      environment, else "latest".
  --print-exports     Also print `export VAR=value` lines for every G4*DATA
                       variable, suitable for `eval "$(... --print-exports)"`
                       or for passing individually to run-in-o2-env.zsh.
  --strict             Exit non-zero if the AliEn token is missing/expired,
                        or if ALICEO2_CCDB_NOTOKENCHECK is set to a truthy
                        value (canonical/validation fixture generation must
                        never rely on that bypass). Without --strict this
                        script only reports, it never fails on token state.
  --help                Show this message and exit 0.

What this reports:
  - AliEn token visibility and expiry (`alien-token-info`).
  - aliBuild provenance (version/revision/alidist-hash/root) for O2, ROOT,
    GEANT4, GEANT4_VMC, VMC -- read directly from the env vars aliBuild's
    generated profile exports for each dependency, never guessed from a
    directory listing or hand-transcribed.
  - o2-sim's own self-reported version string, best-effort, bounded by
    `timeout` so a driver process can never hang the preflight (skipped
    with a note if `timeout` is not on PATH).
  - Every G4*DATA variable `geant4-config --datasets` says Geant4 needs,
    and whether its directory actually exists. This is the ONLY source for
    these variable names in this skill; do not hand-maintain a separate
    list (see SKILL.md's G4LEDATA/G4EMLOWDATA incident note for why).
  - Whether ALICEO2_CCDB_NOTOKENCHECK is set, and to what.

This script never closes stdin itself and never backgrounds anything; it is
safe to run interactively. Batch/backgrounded simulation and reconstruction
commands need their OWN stdin closed explicitly (`< /dev/null`) -- see
SKILL.md's "standalone backgrounded DPL workflow" section. Do not add a
blanket stdin-closing mode here that would also apply to interactive use.
EOF
}

PRINT_EXPORTS=0
STRICT=0

while (( $# > 0 )); do
  case "$1" in
    --package) O2_PACKAGE="$2"; shift 2 ;;
    --print-exports) PRINT_EXPORTS=1; shift ;;
    --strict) STRICT=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) print -u2 -- "${A2E_PROG}: unrecognized argument: $1"; usage 1>&2; exit 2 ;;
  esac
done

if ! a2e_discover_environment; then
  a2e_err "$A2E_ERROR"
  if [[ -d "${O2_INSTALL_ROOT:-}" ]]; then
    a2e_err "available O2 packages under ${O2_INSTALL_ROOT}:"
    a2e_list_available_packages "${O2_INSTALL_ROOT}" 1>&2
  fi
  exit 2
fi
a2e_source_profile

STATUS=0

print -- "package requested: ${O2_PACKAGE}"
print -- "package resolved:  ${O2_PACKAGE_RESOLVED}"
print --

# --- AliEn token ------------------------------------------------------
print -- "-- AliEn token --"
TOKEN_OK=0
if command -v alien-token-info >/dev/null 2>&1; then
  TOKEN_OUTPUT="$(alien-token-info 2>&1 < /dev/null || true)"
  if [[ -n "$TOKEN_OUTPUT" ]]; then
    print -r -- "$TOKEN_OUTPUT"
    if print -r -- "$TOKEN_OUTPUT" | grep -q "^EXPIRE"; then
      TOKEN_OK=1
    fi
  else
    print -- "alien-token-info produced no output (no token, or not logged in)."
  fi
else
  print -- "alien-token-info not found on PATH after sourcing ${O2_PACKAGE_RESOLVED}."
fi
if (( ! TOKEN_OK )); then
  print -- "WARNING: no valid AliEn token detected. Canonical CCDB/Grid reads (e.g. alignment objects via alice-ccdb.cern.ch) require one; run alien-token-init first."
  (( STRICT )) && STATUS=1
fi
print --

# --- ALICEO2_CCDB_NOTOKENCHECK -----------------------------------------
print -- "-- ALICEO2_CCDB_NOTOKENCHECK --"
if [[ -n "${ALICEO2_CCDB_NOTOKENCHECK:-}" ]]; then
  print -- "SET to '${ALICEO2_CCDB_NOTOKENCHECK}' -- this bypasses the alien-token fatal check in CcdbApi::initTGrid()."
  print -- "Do not use this for canonical/validation fixture generation; get a real token instead (see alien-token-init)."
  if (( STRICT )); then
    case "${ALICEO2_CCDB_NOTOKENCHECK}" in
      0|""|false|False|FALSE) ;;
      *) STATUS=1 ;;
    esac
  fi
else
  print -- "unset (expected for canonical validation runs)."
fi
print --

# --- aliBuild package provenance ---------------------------------------
print -- "-- aliBuild package provenance --"
for prefix in O2 ROOT GEANT4 GEANT4_VMC VMC; do
  a2e_package_provenance "$prefix" || true
done

if command -v timeout >/dev/null 2>&1; then
  print -- "o2-sim self-report (bounded, best-effort):"
  timeout --kill-after=2 3 o2-sim --version < /dev/null 2>&1 | head -2 || true
else
  print -- "o2-sim self-report: skipped (\`timeout\` not on PATH; relying on aliBuild O2 provenance above only)."
fi
print --

# --- G4*DATA, derived mechanically --------------------------------------
print -- "-- G4*DATA (from geant4-config --datasets) --"
if ! command -v geant4-config >/dev/null 2>&1; then
  print -u2 -- "${A2E_PROG}: geant4-config not found on PATH after sourcing ${O2_PACKAGE_RESOLVED}."
  exit 2
fi

DATASETS="$(geant4-config --datasets)"
[[ -n "$DATASETS" ]] || fail "geant4-config --datasets produced no output."

ALL_DATASETS_OK=1
EXPORT_LINES=()
while IFS=' ' read -r _name var path; do
  [[ -n "$var" ]] || continue
  if [[ -d "$path" ]]; then
    print -- "OK   ${var}=${path}"
  else
    print -- "MISSING ${var}=${path} (directory does not exist)"
    ALL_DATASETS_OK=0
  fi
  EXPORT_LINES+=("export ${var}=${path}")
done <<< "$DATASETS"

if (( ! ALL_DATASETS_OK )); then
  print -u2 -- "${A2E_PROG}: one or more G4*DATA directories are missing; a batch simulation would fail or silently use wrong defaults."
  STATUS=1
fi

if (( PRINT_EXPORTS )); then
  print --
  print -- "-- exports --"
  for line in "${EXPORT_LINES[@]}"; do
    print -r -- "$line"
  done
fi

exit $STATUS

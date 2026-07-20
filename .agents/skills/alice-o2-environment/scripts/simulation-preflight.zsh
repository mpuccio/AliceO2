#!/bin/zsh
#
# Preflight checks for batch ALICE simulation (o2-sim / TGeant4): AliEn
# token status, aliBuild provenance for O2/ROOT/Geant4/Geant4VMC/VMC, and
# every G4*DATA dataset variable Geant4 needs -- derived MECHANICALLY from
# `geant4-config --datasets`, never hand-transcribed. See SKILL.md for
# usage.
#
# This script only inspects and reports; it NEVER executes o2-sim (there is
# no safe metadata-only "just print the version" mode -- see the provenance
# section below) and creates no files of its own. Use --print-exports to
# get eval-able G4*DATA assignments for a real run.

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
                      environment, else "latest" (the current moving
                      package alias selected by this installation -- not a
                      claim that it is "newest" or "working").
  --print-exports     Also print `export VAR=value` lines for every G4*DATA
                       variable, suitable for `eval "$(... --print-exports)"`
                       or for passing individually to run-in-o2-env.zsh.
  --strict             Exit non-zero if the AliEn token is missing,
                        unparseable, or expired (see "AliEn token" below),
                        or if ALICEO2_CCDB_NOTOKENCHECK is set to a truthy
                        value (canonical/validation fixture generation must
                        never rely on that bypass). Without --strict this
                        script only reports, it never fails on token state.
  --help                Show this message and exit 0.

What this reports:
  - AliEn token status. Runs `alien-token-info`, honors its exit status,
    and ALSO independently checks the EXPIRE timestamp it prints against
    the current time -- alien-token-info's own exit status is 0 for any
    parseable certificate regardless of whether it is currently within its
    validity period, so exit status alone is not suffient to call a token
    "valid".
  - aliBuild provenance (version/revision/alidist-hash/root) for O2, ROOT,
    GEANT4, GEANT4_VMC, VMC -- read directly from the env vars aliBuild's
    generated profile exports for each dependency, never guessed from a
    directory listing or hand-transcribed. This is the MANDATORY and ONLY
    source of o2-sim/ROOT/Geant4/Geant4VMC version/provenance this script
    uses -- it never executes o2-sim itself (no observed invocation of this
    o2-sim build has a safe "print version and exit" mode: `o2-sim
    --version` starts real simulator initialization -- binding IPC
    sockets, contacting CCDB, and in one observed case leaving behind a
    core dump -- rather than just printing a version and exiting).
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
    --package)
      (( $# >= 2 )) || fail "missing value for $1"
      [[ "$2" != -* ]] || fail "missing value for $1 (got option-like '$2')"
      O2_PACKAGE="$2"; shift 2 ;;
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
if a2e_check_alien_token; then
  print -- "OK: ${A2E_TOKEN_DETAIL}"
else
  print -- "WARNING (status=${A2E_TOKEN_STATUS}): ${A2E_TOKEN_DETAIL}"
  print -- "Canonical CCDB/Grid reads (e.g. alignment objects via alice-ccdb.cern.ch) require a currently-valid token; run alien-token-init if needed."
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
# This is the ONLY source of o2-sim/ROOT/Geant4/Geant4VMC provenance this
# script uses. It deliberately never executes o2-sim: this installation's
# o2-sim has no observed safe metadata-only invocation -- `o2-sim
# --version` starts real simulator initialization (IPC sockets, CCDB
# contact, and in one observed case a leaked core dump) rather than
# printing a version and exiting. If a genuinely safe, documented
# metadata-only mechanism is ever confirmed for a given o2-sim build, it
# can be added here explicitly and bounded by `timeout`; until then, do not
# invoke o2-sim from this preflight under any flag.
print -- "-- aliBuild package provenance (mandatory and only version source; o2-sim itself is never executed) --"
for prefix in O2 ROOT GEANT4 GEANT4_VMC VMC; do
  a2e_package_provenance "$prefix" || true
done
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

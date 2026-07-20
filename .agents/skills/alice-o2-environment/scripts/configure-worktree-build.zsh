#!/bin/zsh
#
# Deterministically configure (or safely reuse) a CMake/Ninja build for an
# O2 source worktree, through the aliBuild environment runner. See
# SKILL.md for usage. Concrete reference: this is the exact configuration
# shape validated for the Gate 3 index-table paired A/B build (base
# 3696b2cc57 / feature commit, both against O2_PACKAGE=daily-20260717-0700-local1).

set -e

SCRIPT_DIR="${0:A:h}"
A2E_PROG="configure-worktree-build.zsh"
source "${SCRIPT_DIR}/lib/env-common.zsh"

fail() { print -u2 -- "${A2E_PROG}: $*"; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  configure-worktree-build.zsh --build DIR [--source DIR] [options]

Required:
  --build DIR         Build directory. Must not be inside --source. Created
                       if it does not exist. May also be given via
                       O2_BUILD_DIR; --build takes precedence.

Optional:
  --source DIR         Source worktree root. Defaults to
                        `git rev-parse --show-toplevel` from the current
                        directory when that succeeds ("safe" default);
                        otherwise --source is required.
  --build-type TYPE     CMAKE_BUILD_TYPE. Default: RelWithDebInfo.
  --package NAME        O2_PACKAGE for the runner (alias or exact resolved
                         name). Default: inherits O2_PACKAGE from the
                         environment, else "latest". Pin an exact resolved
                         name for reproducible/paired validation builds.
  --check-only           Only report whether --build is fresh, compatible,
                          or mismatched; never invokes CMake and never
                          creates --build if it does not exist.
  --help                  Show this message and exit 0.

Behavior:
  - Refuses (exit 2) if --build already exists with a CMakeCache.txt whose
    CMAKE_HOME_DIRECTORY resolves to a DIFFERENT source worktree. It never
    deletes, resets, or otherwise modifies that build directory -- pick a
    different --build path, or remove the stale one yourself if that is
    really what you want.
  - If --build has a CMakeCache.txt matching --source, reuses it: CMake's
    own configure step is safe and idempotent to re-run against an existing
    matching cache (this is the "safe reuse" path), so it configures again
    to pick up any CMakeLists.txt changes rather than silently trusting a
    possibly-stale cache.
  - If --build does not exist or is empty, configures fresh.
  - Always configures CMAKE_INSTALL_PREFIX=<build>/stage and -G Ninja, and
    always verifies CMAKE_HOME_DIRECTORY in the resulting cache identifies
    --source before reporting success.
  - Prints the exact source, build, build type, stage path, requested
    package, and resolved package -- before running CMake, so the values
    are visible even if CMake then fails.
EOF
}

SOURCE_DIR=""
BUILD_DIR="${O2_BUILD_DIR:-}"
BUILD_TYPE="RelWithDebInfo"
CHECK_ONLY=0

while (( $# > 0 )); do
  case "$1" in
    --source) SOURCE_DIR="$2"; shift 2 ;;
    --build) BUILD_DIR="$2"; shift 2 ;;
    --build-type) BUILD_TYPE="$2"; shift 2 ;;
    --package) O2_PACKAGE="$2"; shift 2 ;;
    --check-only) CHECK_ONLY=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) print -u2 -- "${A2E_PROG}: unrecognized argument: $1"; usage 1>&2; exit 2 ;;
  esac
done

if [[ -z "$SOURCE_DIR" ]]; then
  if SOURCE_DIR="$(git rev-parse --show-toplevel 2>/dev/null)"; then
    : # safe default: cwd is inside a git worktree
  else
    fail "no --source given and the current directory is not inside a git worktree (git rev-parse --show-toplevel failed); pass --source explicitly."
  fi
fi
[[ -d "$SOURCE_DIR" ]] || fail "--source does not exist: $SOURCE_DIR"
SOURCE_DIR="${SOURCE_DIR:A}"
[[ -f "${SOURCE_DIR}/CMakeLists.txt" ]] || fail "--source does not look like a CMake project root (no CMakeLists.txt): $SOURCE_DIR"

[[ -n "$BUILD_DIR" ]] || fail "no --build given and O2_BUILD_DIR is unset; a build directory must be explicit."
# Resolve as far as possible without requiring the directory to exist yet.
if [[ -d "$BUILD_DIR" ]]; then
  BUILD_DIR="${BUILD_DIR:A}"
else
  BUILD_DIR="${BUILD_DIR:h:A}/${BUILD_DIR:t}"
fi

case "${BUILD_DIR}/" in
  "${SOURCE_DIR}/"*) fail "--build must be outside --source; got build=${BUILD_DIR} nested inside source=${SOURCE_DIR}." ;;
esac

CACHE_FILE="${BUILD_DIR}/CMakeCache.txt"
MODE="fresh"
if [[ -f "$CACHE_FILE" ]]; then
  EXISTING_HOME="$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' "$CACHE_FILE" | head -1)"
  if [[ -n "$EXISTING_HOME" ]]; then
    EXISTING_HOME_RESOLVED="${EXISTING_HOME:A}"
    if [[ "$EXISTING_HOME_RESOLVED" == "$SOURCE_DIR" ]]; then
      MODE="reuse"
    else
      print -u2 -- "${A2E_PROG}: refusing to touch ${BUILD_DIR}: it is configured from a DIFFERENT source."
      print -u2 -- "${A2E_PROG}:   existing CMAKE_HOME_DIRECTORY: ${EXISTING_HOME_RESOLVED}"
      print -u2 -- "${A2E_PROG}:   requested --source:            ${SOURCE_DIR}"
      print -u2 -- "${A2E_PROG}: pick a different --build path, or remove the stale directory yourself. Nothing was deleted or reset."
      exit 2
    fi
  fi
fi

if ! a2e_discover_environment; then
  a2e_err "$A2E_ERROR"
  if [[ -d "${O2_INSTALL_ROOT:-}" ]]; then
    a2e_err "available O2 packages under ${O2_INSTALL_ROOT}:"
    a2e_list_available_packages "${O2_INSTALL_ROOT}" 1>&2
  fi
  exit 2
fi

STAGE_DIR="${BUILD_DIR}/stage"
print -- "source:            ${SOURCE_DIR}"
print -- "build:             ${BUILD_DIR}"
print -- "build type:        ${BUILD_TYPE}"
print -- "stage:             ${STAGE_DIR}"
print -- "mode:              ${MODE}"
print -- "package requested: ${O2_PACKAGE}"
print -- "package resolved:  ${O2_PACKAGE_RESOLVED}"

if (( CHECK_ONLY )); then
  print -- "check-only: no CMake invocation performed."
  exit 0
fi

mkdir -p "$BUILD_DIR"

RUNNER="${SCRIPT_DIR}/run-in-o2-env.zsh"
O2_PACKAGE="$O2_PACKAGE" O2_BUILD_DIR="$BUILD_DIR" "$RUNNER" -- \
  cmake -S "$SOURCE_DIR" -B "$BUILD_DIR" -G Ninja \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DCMAKE_INSTALL_PREFIX="$STAGE_DIR"

FINAL_HOME="$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' "$CACHE_FILE" | head -1)"
FINAL_HOME_RESOLVED="${FINAL_HOME:A}"
if [[ "$FINAL_HOME_RESOLVED" != "$SOURCE_DIR" ]]; then
  fail "post-configure sanity check failed: CMAKE_HOME_DIRECTORY is '${FINAL_HOME_RESOLVED}', expected '${SOURCE_DIR}'."
fi

print -- "configured OK: CMAKE_HOME_DIRECTORY verified as ${SOURCE_DIR}"

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
  --build DIR         Build directory. Must not be inside --source. May
                       also be given via O2_BUILD_DIR; --build takes
                       precedence.

Optional:
  --source DIR         Source worktree root. Defaults to
                        `git rev-parse --show-toplevel` from the current
                        directory when that succeeds ("safe" default);
                        otherwise --source is required.
  --build-type TYPE     CMAKE_BUILD_TYPE. Default: RelWithDebInfo.
  --package NAME        O2_PACKAGE for the runner (alias or exact resolved
                         name). Default: inherits O2_PACKAGE from the
                         environment, else "latest" (the current moving
                         package alias selected by this installation).
                         Pin an exact resolved name for reproducible/paired
                         validation builds.
  --check-only           Only report whether --build is fresh, compatible,
                          or mismatched; never invokes CMake and never
                          creates or writes into --build.
  --help                  Show this message and exit 0.

--build directory states (checked in this order, nothing is ever deleted
or reset by this tool):
  - does not exist                                -> fresh, permitted.
  - exists and is empty                            -> fresh, permitted.
  - exists, non-empty, no CMakeCache.txt            -> REFUSED (exit 2):
    not a build directory this tool created; pick a different path or
    clear it yourself.
  - CMakeCache.txt present but missing
    CMAKE_HOME_DIRECTORY                            -> REFUSED (exit 2):
    malformed/incomplete cache.
  - CMakeCache.txt's CMAKE_HOME_DIRECTORY resolves
    to a DIFFERENT source                           -> REFUSED (exit 2).
  - CMakeCache.txt's CMAKE_HOME_DIRECTORY matches
    --source                                        -> reused: CMake's own
    configure step is safe and idempotent to re-run against a matching
    cache, so it reconfigures to pick up any CMakeLists.txt changes rather
    than silently trusting a possibly-stale cache.

Always configures -G Ninja, -DCMAKE_INSTALL_PREFIX=<build>/stage, and
-DBUILD_TESTING=ON (O2's CMakeLists.txt already calls CMake's own
`include(CTest)`, which defaults BUILD_TESTING to ON; this is set
explicitly anyway so behavior does not depend on an undocumented CMake
default). Always verifies CMAKE_HOME_DIRECTORY in the resulting cache
identifies --source before reporting success. Prints the exact source,
build, build type, stage path, requested package, and resolved package --
before running CMake, so the values are visible even if CMake then fails.

Note on CTest test names: O2's `o2_add_test()` registers each CTest test
under its SOURCE FILE's relative path (e.g.
"Detectors/ITSMFT/common/tracking/test/testLayerMask.cxx"), not the Ninja
target or staged executable name. Use `ctest -N` with no filter, or filter
by `-L <label>` (labels match `o2_add_test(... LABELS ...)`, typically the
detector/component name), rather than guessing a `-R` pattern from a target
name.
EOF
}

SOURCE_DIR=""
BUILD_DIR="${O2_BUILD_DIR:-}"
BUILD_TYPE="RelWithDebInfo"
CHECK_ONLY=0

while (( $# > 0 )); do
  case "$1" in
    --source)
      (( $# >= 2 )) || fail "missing value for $1"
      [[ "$2" != -* ]] || fail "missing value for $1 (got option-like '$2')"
      SOURCE_DIR="$2"; shift 2 ;;
    --build)
      (( $# >= 2 )) || fail "missing value for $1"
      [[ "$2" != -* ]] || fail "missing value for $1 (got option-like '$2')"
      BUILD_DIR="$2"; shift 2 ;;
    --build-type)
      (( $# >= 2 )) || fail "missing value for $1"
      [[ "$2" != -* ]] || fail "missing value for $1 (got option-like '$2')"
      BUILD_TYPE="$2"; shift 2 ;;
    --package)
      (( $# >= 2 )) || fail "missing value for $1"
      [[ "$2" != -* ]] || fail "missing value for $1 (got option-like '$2')"
      O2_PACKAGE="$2"; shift 2 ;;
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

BUILD_STATE="nonexistent"
if [[ -d "$BUILD_DIR" ]]; then
  if [[ -z "$(ls -A "$BUILD_DIR" 2>/dev/null)" ]]; then
    BUILD_STATE="empty"
  else
    BUILD_STATE="nonempty"
  fi
fi

CACHE_FILE="${BUILD_DIR}/CMakeCache.txt"
MODE="fresh"
if [[ "$BUILD_STATE" == "nonempty" ]]; then
  if [[ ! -f "$CACHE_FILE" ]]; then
    print -u2 -- "${A2E_PROG}: refusing to use ${BUILD_DIR}: it already exists, is non-empty, and has no CMakeCache.txt (not a build directory this tool created)."
    print -u2 -- "${A2E_PROG}: pick a different --build path, or clear it yourself if that is really what you want. Nothing was modified."
    exit 2
  fi

  EXISTING_HOME="$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' "$CACHE_FILE" | head -1)"
  if [[ -z "$EXISTING_HOME" ]]; then
    print -u2 -- "${A2E_PROG}: refusing to use ${BUILD_DIR}: it has a CMakeCache.txt but no CMAKE_HOME_DIRECTORY entry (malformed or incomplete cache)."
    print -u2 -- "${A2E_PROG}: pick a different --build path, or remove the stale directory yourself. Nothing was deleted or modified."
    exit 2
  fi

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
print -- "build directory:   ${BUILD_STATE} -> mode=${MODE}"
print -- "package requested: ${O2_PACKAGE}"
print -- "package resolved:  ${O2_PACKAGE_RESOLVED}"

if (( CHECK_ONLY )); then
  print -- "check-only: no CMake invocation performed, nothing created."
  exit 0
fi

mkdir -p "$BUILD_DIR"

RUNNER="${SCRIPT_DIR}/run-in-o2-env.zsh"
O2_PACKAGE="$O2_PACKAGE" O2_BUILD_DIR="$BUILD_DIR" "$RUNNER" -- \
  cmake -S "$SOURCE_DIR" -B "$BUILD_DIR" -G Ninja \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DCMAKE_INSTALL_PREFIX="$STAGE_DIR" \
    -DBUILD_TESTING=ON

FINAL_HOME="$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' "$CACHE_FILE" | head -1)"
FINAL_HOME_RESOLVED="${FINAL_HOME:A}"
if [[ "$FINAL_HOME_RESOLVED" != "$SOURCE_DIR" ]]; then
  fail "post-configure sanity check failed: CMAKE_HOME_DIRECTORY is '${FINAL_HOME_RESOLVED}', expected '${SOURCE_DIR}'."
fi

print -- "configured OK: CMAKE_HOME_DIRECTORY verified as ${SOURCE_DIR}"

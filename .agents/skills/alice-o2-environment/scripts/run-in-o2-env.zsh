#!/bin/zsh

set -e

fail()
{
  print -u2 -- "alice-o2-environment: $*"
  exit 2
}

if [[ -z "${ALIBUILD_WORK_DIR:-}" ]]; then
  fail "ALIBUILD_WORK_DIR is not exported. Start the agent from a configured shell or export it before invoking this runner."
fi

if [[ -z "${ALIBUILD_ARCH_PREFIX:-}" ]]; then
  if ! command -v aliBuild >/dev/null 2>&1; then
    fail "ALIBUILD_ARCH_PREFIX is unset and aliBuild is not on PATH, so the architecture cannot be derived."
  fi
  ALIBUILD_ARCH_PREFIX="$(aliBuild architecture)"
  export ALIBUILD_ARCH_PREFIX
fi

O2_PACKAGE="${O2_PACKAGE:-latest-run3-o2}"
O2_PROFILE="${ALIBUILD_WORK_DIR}/${ALIBUILD_ARCH_PREFIX}/O2/${O2_PACKAGE}/etc/profile.d/init.sh"

if [[ ! -r "${O2_PROFILE}" ]]; then
  print -u2 -- "alice-o2-environment: O2 profile not found: ${O2_PROFILE}"
  install_root="${ALIBUILD_WORK_DIR}/${ALIBUILD_ARCH_PREFIX}/O2"
  if [[ -d "${install_root}" ]]; then
    print -u2 -- "alice-o2-environment: available O2 packages:"
    for package in "${install_root}"/*(N:t); do
      print -u2 -- "  ${package}"
    done
  fi
  exit 2
fi

# Generated aliBuild profiles still use WORK_DIR internally. Keep
# ALIBUILD_WORK_DIR as the portable external contract and bridge it here.
WORK_DIR="${ALIBUILD_WORK_DIR}"
export WORK_DIR

source "${O2_PROFILE}"

if [[ -n "${O2_BUILD_DIR:-}" ]]; then
  [[ -d "${O2_BUILD_DIR}" ]] || fail "O2_BUILD_DIR does not exist: ${O2_BUILD_DIR}"
  if [[ -d "${O2_BUILD_DIR}/stage/bin" ]]; then
    export PATH="${O2_BUILD_DIR}/stage/bin:${PATH}"
  fi
  if [[ -d "${O2_BUILD_DIR}/stage/lib" ]]; then
    if [[ "$(uname -s)" == "Darwin" ]]; then
      export DYLD_LIBRARY_PATH="${O2_BUILD_DIR}/stage/lib:${DYLD_LIBRARY_PATH:-}"
    else
      export LD_LIBRARY_PATH="${O2_BUILD_DIR}/stage/lib:${LD_LIBRARY_PATH:-}"
    fi
  fi
fi

if [[ "${1:-}" == "--check" ]]; then
  print -- "ALIBUILD_WORK_DIR=${ALIBUILD_WORK_DIR}"
  print -- "ALIBUILD_ARCH_PREFIX=${ALIBUILD_ARCH_PREFIX}"
  print -- "O2_PACKAGE=${O2_PACKAGE}"
  print -- "O2_PROFILE=${O2_PROFILE}"
  print -- "O2_BUILD_DIR=${O2_BUILD_DIR:-}"
  print -- "ninja=$(command -v ninja 2>/dev/null || true)"
  print -- "root=$(command -v root 2>/dev/null || true)"
  print -- "o2-sim=$(command -v o2-sim 2>/dev/null || true)"
  exit 0
fi

if [[ "${1:-}" == "--" ]]; then
  shift
fi

(( $# > 0 )) || fail "no command supplied; use --check or -- <command> [args...]"

exec "$@"

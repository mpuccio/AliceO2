#!/bin/zsh
#
# Run a command inside the aliBuild-sourced O2 environment. See SKILL.md for
# usage. `--check` and `--help` never source a target command; every other
# invocation sources the resolved O2 profile in this same process and then
# execs the given command, so environment state never has to survive across
# separate tool calls.

set -e

SCRIPT_DIR="${0:A:h}"
A2E_PROG="run-in-o2-env.zsh"
source "${SCRIPT_DIR}/lib/env-common.zsh"

fail() { print -u2 -- "alice-o2-environment: $*"; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  run-in-o2-env.zsh --check
  run-in-o2-env.zsh -- <command> [args...]
  run-in-o2-env.zsh --help

Sources the aliBuild O2 profile for O2_PACKAGE (default: latest) and either
reports the resolved environment (--check) or execs <command> inside it.

Environment variables:
  ALIBUILD_WORK_DIR    required. aliBuild sw/ prefix (e.g. .../alice/sw).
  ALIBUILD_ARCH_PREFIX optional. Derived via `aliBuild architecture` if unset.
  O2_PACKAGE           optional. An alias (e.g. "latest") or an exact
                        resolved package name (e.g. "daily-20260717-0700-local1").
                        Default: latest. Pin an exact name for reproducible
                        validation; never rely on an alias for that.
  O2_BUILD_DIR          optional. A configured worktree build directory;
                        when set, its stage/bin and stage/lib are prepended
                        to PATH/DYLD_LIBRARY_PATH (or LD_LIBRARY_PATH) so
                        branch-built executables take precedence over
                        installed ones.

--check prints the requested package alias AND its fully resolved target
(never just one or the other), plus the tool paths actually on PATH after
sourcing. It never scans the install root or chooses a package by
modification time, and never silently substitutes a different package when
the requested one is missing -- it fails with a diagnostic listing instead.
EOF
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage
  exit 0
fi

if ! a2e_discover_environment; then
  a2e_err "$A2E_ERROR"
  if [[ -d "${O2_INSTALL_ROOT:-}" ]]; then
    a2e_err "available O2 packages under ${O2_INSTALL_ROOT}:"
    a2e_list_available_packages "${O2_INSTALL_ROOT}" 1>&2
  fi
  exit 2
fi

a2e_source_profile

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
  print -- "O2_PACKAGE_RESOLVED=${O2_PACKAGE_RESOLVED}"
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

(( $# > 0 )) || fail "no command supplied; use --check, --help, or -- <command> [args...]"

exec "$@"

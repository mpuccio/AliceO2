# alice-o2-environment shared library: environment discovery, O2 package
# resolution, and profile sourcing shared by every script in this skill.
#
# Sourced, never executed directly. Provider-neutral plain zsh -- no
# dependency on any particular agent host or UI.
#
# Contract for callers:
#   - Set A2E_PROG to the caller's own program name before sourcing (used
#     only in diagnostic prefixes); defaults to "alice-o2-environment".
#   - Call a2e_discover_environment first. On failure it returns 1 and sets
#     A2E_ERROR; the caller decides how to report and exit.
#   - Call a2e_source_profile only after a2e_discover_environment succeeds.
#   - Every exported variable this library sets is a plain, documented name
#     (ALIBUILD_WORK_DIR, ALIBUILD_ARCH_PREFIX, O2_PACKAGE,
#     O2_PACKAGE_RESOLVED, O2_PROFILE, O2_INSTALL_ROOT) -- no hidden state.

: "${A2E_PROG:=alice-o2-environment}"

a2e_err() { print -u2 -- "${A2E_PROG}: $*"; }

# a2e_resolve_package_target <alias> <install-root>
# Prints the real installed directory name behind an O2-style package alias
# (e.g. "latest" -> "daily-20260717-0700-local1"), following any symlink
# chain physically via `cd -P`. This is the ONLY package-resolution
# mechanism in this skill: it never scans the install root for candidates
# and never chooses by modification time or any other heuristic -- an alias
# resolves to whatever it is actually linked to, or resolution fails.
a2e_resolve_package_target() {
  local alias="$1" root="$2" real
  [[ -d "${root}/${alias}" ]] || return 1
  real="$(cd -P -- "${root}/${alias}" 2>/dev/null && pwd -P)" || return 1
  print -r -- "${real:t}"
}

# a2e_list_available_packages <install-root>
# Diagnostic listing only, printed on failure to help a human/agent pick a
# valid --package value. Never used to select a package automatically.
a2e_list_available_packages() {
  local root="$1" package
  [[ -d "$root" ]] || return 0
  for package in "${root}"/*(N:t); do
    print -r -- "  ${package}"
  done
}

# a2e_discover_environment
# Populates (and exports) ALIBUILD_WORK_DIR, ALIBUILD_ARCH_PREFIX,
# O2_PACKAGE, O2_PACKAGE_RESOLVED, O2_PROFILE, O2_INSTALL_ROOT. Does not
# source the profile. On any failure, sets A2E_ERROR to an actionable
# message and returns 1; O2_INSTALL_ROOT is still set in that case so the
# caller can list available packages.
#
# Package default: O2_PACKAGE defaults to "latest" -- the current moving
# package alias selected by this installation. This script has no way to
# know, and does not claim, that whatever "latest" happens to point to at
# any given time is newest or working; it only reports what it actually
# resolves to (O2_PACKAGE_RESOLVED). Reproducible validation work must pass
# an explicit --package/O2_PACKAGE naming the exact resolved package, never
# an alias.
a2e_discover_environment() {
  A2E_ERROR=""

  if [[ -z "${ALIBUILD_WORK_DIR:-}" ]]; then
    A2E_ERROR="ALIBUILD_WORK_DIR is not exported. Start the agent from a configured shell or export it before invoking this tool."
    return 1
  fi

  if [[ -z "${ALIBUILD_ARCH_PREFIX:-}" ]]; then
    if ! command -v aliBuild >/dev/null 2>&1; then
      A2E_ERROR="ALIBUILD_ARCH_PREFIX is unset and aliBuild is not on PATH, so the architecture cannot be derived."
      return 1
    fi
    ALIBUILD_ARCH_PREFIX="$(aliBuild architecture)"
    export ALIBUILD_ARCH_PREFIX
  fi

  O2_PACKAGE="${O2_PACKAGE:-latest}"
  O2_INSTALL_ROOT="${ALIBUILD_WORK_DIR}/${ALIBUILD_ARCH_PREFIX}/O2"
  O2_PROFILE="${O2_INSTALL_ROOT}/${O2_PACKAGE}/etc/profile.d/init.sh"
  export O2_PACKAGE O2_PROFILE O2_INSTALL_ROOT

  if [[ ! -r "${O2_PROFILE}" ]]; then
    A2E_ERROR="O2 profile not found for package '${O2_PACKAGE}': ${O2_PROFILE}"
    return 1
  fi

  if O2_PACKAGE_RESOLVED="$(a2e_resolve_package_target "${O2_PACKAGE}" "${O2_INSTALL_ROOT}")"; then
    export O2_PACKAGE_RESOLVED
  else
    # The profile is readable (checked above) but physical resolution still
    # failed (race between the two checks, a non-symlink path `cd -P`
    # couldn't handle, or similar). This must be a hard failure, not a
    # fallback to reporting the alias as if it were its own resolution --
    # this library's whole contract is that O2_PACKAGE_RESOLVED is always a
    # real, physically-verified target, never a guess.
    A2E_ERROR="could not physically resolve package '${O2_PACKAGE}' under ${O2_INSTALL_ROOT} even though its profile is readable; refusing to report an unresolved alias as if it were resolved."
    return 1
  fi

  return 0
}

# a2e_source_profile
# Sources the resolved O2 profile into the current shell. Must be called
# after a2e_discover_environment succeeds. Bridges ALIBUILD_WORK_DIR to the
# legacy WORK_DIR name generated O2 profiles still expect internally, so
# callers only ever need to set the portable ALIBUILD_WORK_DIR name.
a2e_source_profile() {
  WORK_DIR="${ALIBUILD_WORK_DIR}"
  export WORK_DIR
  source "${O2_PROFILE}"
}

# a2e_package_provenance <pkg-env-prefix>
# Prints one line of aliBuild-exported provenance for an already-sourced
# package (e.g. "ROOT", "GEANT4", "GEANT4_VMC", "VMC", "O2"): version,
# revision, alidist recipe hash, and install root. These come directly from
# the env vars aliBuild's generated profiles export for every dependency
# (<PKG>_VERSION/_REVISION/_HASH/_ROOT) -- never hand-transcribed, never
# guessed from a directory listing. Prints nothing (not an error) if the
# package's env vars are not set, e.g. because the profile does not depend
# on it.
a2e_package_provenance() {
  local prefix="$1"
  local version revision hash root
  version="$(eval "print -r -- \${${prefix}_VERSION:-}")"
  revision="$(eval "print -r -- \${${prefix}_REVISION:-}")"
  hash="$(eval "print -r -- \${${prefix}_HASH:-}")"
  root="$(eval "print -r -- \${${prefix}_ROOT:-}")"
  [[ -n "$root" ]] || return 0
  print -r -- "${prefix}: version=${version} revision=${revision} hash=${hash} root=${root}"
}

# a2e_check_alien_token
# Runs `alien-token-info`, HONORS its exit status (the previous version of
# this check discarded it with `|| true`), and independently verifies the
# EXPIRE timestamp it prints against the current time.
#
# This second step is necessary, not decorative: `alien-token-info` invokes
# xjalienfs's CertInfo(), which parses and prints certificate fields
# unconditionally and returns 0 for ANY parseable certificate -- it never
# checks whether "now" falls inside [BEGIN, EXPIRE]. Only the unexposed
# `token-verify` subcommand (no standalone `alien-token-verify` binary
# ships in this install) does real chain/expiry verification. So an exit
# status of 0 alone does NOT mean the token is currently valid; comparing
# the EXPIRE field it does print, portably, against `date`, is the
# strongest check actually available from what this tool exposes.
#
# Sets A2E_TOKEN_STATUS to one of: ok, missing-command, nonzero-exit,
# malformed-output, expired. Sets A2E_TOKEN_DETAIL to a human-readable
# detail line. Prints alien-token-info's own stdout+stderr (if any) so the
# caller can relay it. Returns 0 only for "ok".
a2e_check_alien_token() {
  A2E_TOKEN_STATUS="ok"
  A2E_TOKEN_DETAIL=""

  if ! command -v alien-token-info >/dev/null 2>&1; then
    A2E_TOKEN_STATUS="missing-command"
    A2E_TOKEN_DETAIL="alien-token-info not found on PATH."
    return 1
  fi

  local output rc
  output="$(alien-token-info 2>&1 < /dev/null)"
  rc=$?
  [[ -n "$output" ]] && print -r -- "$output"

  if (( rc != 0 )); then
    A2E_TOKEN_STATUS="nonzero-exit"
    A2E_TOKEN_DETAIL="alien-token-info exited ${rc}."
    return 1
  fi

  local dn_line begin_line expire_line
  dn_line="$(print -r -- "$output" | grep -m1 '^DN')"
  begin_line="$(print -r -- "$output" | grep -m1 '^BEGIN')"
  expire_line="$(print -r -- "$output" | grep -m1 '^EXPIRE')"
  if [[ -z "$dn_line" || -z "$begin_line" || -z "$expire_line" ]]; then
    A2E_TOKEN_STATUS="malformed-output"
    A2E_TOKEN_DETAIL="alien-token-info exited 0 but its output is missing DN/BEGIN/EXPIRE fields."
    return 1
  fi

  local expire_str expire_epoch now_epoch
  expire_str="${expire_line#EXPIRE >>> }"
  if [[ "$(uname -s)" == "Darwin" ]]; then
    expire_epoch="$(date -j -u -f '%Y-%m-%d %H:%M:%S' "$expire_str" +%s 2>/dev/null)"
  else
    expire_epoch="$(date -u -d "$expire_str" +%s 2>/dev/null)"
  fi
  if [[ -z "$expire_epoch" ]]; then
    A2E_TOKEN_STATUS="malformed-output"
    A2E_TOKEN_DETAIL="could not parse EXPIRE timestamp '${expire_str}' as UTC."
    return 1
  fi
  now_epoch="$(date -u +%s)"
  if (( expire_epoch <= now_epoch )); then
    A2E_TOKEN_STATUS="expired"
    A2E_TOKEN_DETAIL="token EXPIRE ${expire_str} UTC has already passed."
    return 1
  fi

  A2E_TOKEN_DETAIL="token valid until ${expire_str} UTC."
  return 0
}

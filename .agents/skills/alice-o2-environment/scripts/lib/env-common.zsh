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
# Package default: O2_PACKAGE defaults to "latest" (the moving alias),
# matching the policy that ordinary unpinned work should track the newest
# working daily build. Reproducible validation work must pass an explicit
# --package/O2_PACKAGE naming the exact resolved package, never an alias.
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
    # The alias resolved to a readable profile above but somehow cannot be
    # stat'd as a directory now (race, or a non-symlink path) -- fall back
    # to reporting the alias itself as its own resolution rather than
    # silently guessing at another package.
    O2_PACKAGE_RESOLVED="${O2_PACKAGE}"
    export O2_PACKAGE_RESOLVED
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

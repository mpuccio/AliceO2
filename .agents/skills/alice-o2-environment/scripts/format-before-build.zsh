#!/bin/zsh
#
# Apply the repository's git-clang-format result to the working tree without
# staging changes. New, untracked C/C++ files are formatted directly because
# git-clang-format cannot see them.

set -eu

PROG="format-before-build.zsh"

fail() { print -u2 -- "${PROG}: $*"; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  format-before-build.zsh [--source /path/to/worktree]
  format-before-build.zsh --help

Formats every tracked C/C++ change relative to HEAD using git clang-format,
then formats untracked C/C++ files directly. Changes are applied only to the
working tree: the script never stages, commits, resets, or restores files.

Run this immediately before every source build and test invocation.
EOF
}

source_root=""
while (( $# > 0 )); do
  case "$1" in
    --source)
      (( $# >= 2 )) || fail "--source requires a path"
      source_root="$2"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      fail "unknown argument: $1"
      ;;
  esac
done

if [[ -z "${source_root}" ]]; then
  source_root="$(git rev-parse --show-toplevel 2>/dev/null)" ||
    fail "current directory is not inside a Git worktree; pass --source"
fi
source_root="$(cd "${source_root}" && pwd -P)" || fail "source directory does not exist: ${source_root}"
[[ "$(git -C "${source_root}" rev-parse --show-toplevel 2>/dev/null)" == "${source_root}" ]] ||
  fail "--source is not the root of a Git worktree: ${source_root}"
git -C "${source_root}" rev-parse --verify HEAD >/dev/null 2>&1 || fail "worktree has no HEAD commit"
command -v clang-format >/dev/null 2>&1 || fail "clang-format is not on PATH"

changed=0
patch_file="$(mktemp "${TMPDIR:-/tmp}/o2-git-clang-format.XXXXXX")" || fail "cannot create temporary patch"
verify_file="$(mktemp "${TMPDIR:-/tmp}/o2-git-clang-format-verify.XXXXXX")" || fail "cannot create verification patch"
cleanup() { rm -f "${patch_file}" "${verify_file}"; }
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

cd "${source_root}"

# git-clang-format deliberately ignores untracked files, so format those with
# the exact clang-format resolved in the same environment first.
while IFS= read -r -d '' file; do
  [[ -f "${file}" ]] || continue
  if ! clang-format --dry-run --Werror -- "${file}" >/dev/null 2>&1; then
    clang-format -i -- "${file}"
    changed=1
  fi
done < <(git ls-files --others --exclude-standard -z -- \
  '*.c' '*.cc' '*.cpp' '*.cxx' '*.h' '*.hh' '*.hpp' '*.hxx' '*.cu' '*.cuh')

# --diff works with mixed staged and unstaged edits. It returns 1 when it emits
# a patch, so accept both 0 (clean) and 1 (formatting needed). Applying its
# patch changes only the working tree and preserves the user's index exactly.
if git clang-format --diff HEAD >"${patch_file}"; then
  format_status=0
else
  format_status=$?
  (( format_status == 1 )) || fail "git clang-format failed with status ${format_status}"
fi
if grep -q '^diff --git ' "${patch_file}"; then
  git apply --whitespace=nowarn "${patch_file}"
  changed=1
fi

if git clang-format --diff HEAD >"${verify_file}"; then
  verify_status=0
else
  verify_status=$?
  (( verify_status == 1 )) || fail "git clang-format verification failed with status ${verify_status}"
fi
if grep -q '^diff --git ' "${verify_file}"; then
  fail "git-clang-format still reports changes after formatting"
fi

if (( changed )); then
  print -- "${PROG}: formatting applied; rebuild before running tests"
else
  print -- "${PROG}: clean"
fi

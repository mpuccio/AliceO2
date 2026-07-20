---
name: alice-o2-environment
description: Run ALICE O2 builds, tests, reconstruction, simulation, and command-line tools inside the repository's aliBuild environment. Use whenever a command needs O2, ROOT, Ninja, detector workflows, staged branch binaries, or libraries that are unavailable in the agent's default shell.
---

# ALICE O2 environment

Run environment-dependent commands through `scripts/run-in-o2-env.zsh`. The
runner discovers the installation from `ALIBUILD_WORK_DIR`, sources O2, and
executes the command in that same process. Do not source the profile in one
tool call and run the command in another; environment changes do not
persist across separate tool calls.

This skill is provider-neutral: every script is ordinary POSIX-ish zsh with
no dependency on a particular agent UI, discovered relative to this
repository (`.agents/skills/alice-o2-environment/scripts/...`). It works
identically from Claude Code and Codex. `agents/openai.yaml` supplies only
Codex's display metadata; it does not change behavior.

## Scripts in this skill

- `scripts/run-in-o2-env.zsh` -- source the O2 profile and exec a command,
  or `--check`/`--help`.
- `scripts/configure-worktree-build.zsh` -- deterministically configure (or
  safely reuse) a CMake/Ninja build for a source worktree.
- `scripts/simulation-preflight.zsh` -- verify AliEn token, aliBuild
  provenance, and every G4\*DATA dataset before a batch simulation.
- `scripts/lib/env-common.zsh` -- shared environment-discovery and
  package-resolution functions sourced by all three scripts above; not
  meant to be invoked directly. If you need version/package-resolution
  logic in a new script, source this rather than re-deriving it.

Every script supports `--help`.

## Environment / package check

```bash
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh --check
```

Prints `ALIBUILD_WORK_DIR`, `ALIBUILD_ARCH_PREFIX`, the **requested**
`O2_PACKAGE` and its **fully resolved target** (`O2_PACKAGE_RESOLVED`), the
profile path, and the resolved `ninja`/`root`/`o2-sim` paths. If this fails,
report its diagnostics. Do not hardcode a developer home directory; require
`ALIBUILD_WORK_DIR` to be exported by the process that started the agent.

The runner derives `ALIBUILD_ARCH_PREFIX` with `aliBuild architecture` when
unset, and bridges `ALIBUILD_WORK_DIR` to the legacy `WORK_DIR` name used
inside generated O2 profiles; callers should not need to set both.

### Package selection policy

- **Ordinary unpinned work** (day-to-day builds, tests, exploring the repo)
  should just use the default: `O2_PACKAGE` unset resolves to `latest` --
  the current moving package alias selected by this installation. This
  skill has no way to know, and does not claim, that whatever `latest`
  happens to point to at any given time is newest or working; it only
  reports what it actually resolves to (`O2_PACKAGE_RESOLVED`).
- **Reproducible/canonical validation work** (fixture generation, paired
  A/B replay, anything whose result you will compare against a recorded
  baseline) must set `O2_PACKAGE` to the **exact resolved package name**
  (e.g. `daily-20260717-0700-local1`), never an alias. Record that exact
  name in whatever you are validating.
- Package resolution never scans the install directory and never chooses by
  modification time -- it only ever follows the alias you asked for
  (`cd -P`, i.e. a real symlink resolution) or fails loudly. A missing
  package is never silently substituted with another one; the runner exits
  2 and lists what is actually installed.

```bash
# ordinary work
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- ninja --version

# pinned, reproducible
O2_PACKAGE=daily-20260717-0700-local1 \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh --check
```

## Fresh worktree configuration

```bash
.agents/skills/alice-o2-environment/scripts/configure-worktree-build.zsh \
  --build /path/to/durable/build \
  [--source /path/to/worktree] [--package NAME] [--build-type RelWithDebInfo]
```

`--source` defaults to `git rev-parse --show-toplevel` from the current
directory when that succeeds; otherwise pass it explicitly. `--build` is
always required (directly or via `O2_BUILD_DIR`) and must be outside
`--source`. It configures `-G Ninja -DCMAKE_INSTALL_PREFIX=<build>/stage`,
verifies `CMAKE_HOME_DIRECTORY` in the resulting cache actually matches
`--source`, and prints the exact source/build/build-type/stage/requested-
package/resolved-package before touching CMake. If `--build` already has a
cache configured from a **different** source it refuses (exit 2) without
deleting or resetting anything -- pick a fresh `--build` path, or remove the
stale directory yourself if that's really what you want. If the cache
already matches `--source`, it safely reconfigures in place (CMake's own
configure step is idempotent). Add `--check-only` to see fresh/reuse/
mismatch without invoking CMake at all.

## Narrow Ninja build

```bash
O2_BUILD_DIR=/path/to/build \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  ninja -C /path/to/build O2lib-ITSMFTTracking
```

## Focused ctest

```bash
O2_BUILD_DIR=/path/to/build \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  ctest --test-dir /path/to/build --output-on-failure -R 'testLayerMask'
```

Prefer a build configured from the active worktree via
`configure-worktree-build.zsh`, which explicitly passes
`-DBUILD_TESTING=ON` (O2's own `include(CTest)` already defaults it to ON,
but the helper does not rely on that default going unchanged). Verify the
source directory in `CMakeCache.txt` before using an unfamiliar build
directory (the configure script does this for you).

**`-R`/test-name pitfall:** O2's `o2_add_test()` registers each CTest test
under its **source file's relative path**
(`Detectors/ITSMFT/common/tracking/test/testLayerMask.cxx`), not the Ninja
target name (`O2test-itsmft-tracking-layermask`) or the staged executable
name (`o2-test-itsmft-tracking-layermask`). A `-R` pattern built from
either of the latter two matches nothing and `ctest` reports "No tests were
found!!!" -- which looks exactly like a build/registration problem but
isn't one. Filter by the actual test name substring (e.g. `testLayerMask`)
or by `-L <label>` (labels match `o2_add_test(... LABELS ...)`, typically
the detector/component name, e.g. `-L itsmft`); run `ctest -N` with no
filter first if unsure what's registered.

## Staged vs. installed execution

- Set `O2_BUILD_DIR` to a configured worktree build when you need **this
  branch's** reconstruction/tracking code (the runner prepends
  `<build>/stage/bin` to `PATH` and `<build>/stage/lib` to
  `DYLD_LIBRARY_PATH`/`LD_LIBRARY_PATH`, so branch binaries take precedence
  over installed ones).
- Omit `O2_BUILD_DIR` to use the **installed** aliBuild executables
  directly -- this is what fixture/simulation generation normally wants,
  since `o2-sim`/`o2-sim-digitizer-workflow` are simulation infrastructure,
  not the code under test.
- For a **paired A/B comparison** (e.g. base commit vs. feature commit),
  configure two separate builds against the **identical, explicitly pinned**
  `O2_PACKAGE`, and run each replay with its own `O2_BUILD_DIR`. Never let
  one leg default to an alias while the other pins a name -- both must
  resolve to the same package.

```bash
# installed simulation binary (no O2_BUILD_DIR)
O2_PACKAGE=daily-20260717-0700-local1 \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  o2-sim -j 1 -n 20 -g pythia8pp -e TGeant4 -m PIPE ITS MFT --run 303000 --seed 20260716 --field -5 -o o2sim

# staged branch reconstruction binary
O2_BUILD_DIR=/path/to/build O2_PACKAGE=daily-20260717-0700-local1 \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  o2-mft-ca-tracker-workflow -b --run --tracking-mode sync
```

## G4 dataset setup and token preflight

```bash
.agents/skills/alice-o2-environment/scripts/simulation-preflight.zsh \
  --package daily-20260717-0700-local1
```

Reports AliEn token status, `ALICEO2_CCDB_NOTOKENCHECK` state, aliBuild
provenance (version/revision/alidist-hash/root) for O2/ROOT/GEANT4/
GEANT4\_VMC/VMC, and every `G4*DATA` variable Geant4 needs. The dataset
variable names and paths are derived **mechanically** from `geant4-config
--datasets` -- never hand-transcribed. Add `--print-exports` to get
`export VAR=value` lines for a real simulation run, and `--strict` to fail
if the token is missing/expired or `ALICEO2_CCDB_NOTOKENCHECK` is set
truthy (canonical validation must never rely on that bypass).

```bash
eval "$(.agents/skills/alice-o2-environment/scripts/simulation-preflight.zsh \
  --package daily-20260717-0700-local1 --print-exports | grep '^export ')"
```

## Running a standalone backgrounded DPL workflow

A DPL driver invoked standalone (not as the downstream end of a `|` pipe)
still probes stdin for a piped workflow spec. Under some backgrounded/
non-interactive execution contexts that probe blocks indefinitely instead
of detecting "no pipe" and proceeding, because stdin is left open rather
than closed. Close it explicitly for any standalone batch invocation:

```bash
O2_PACKAGE=daily-20260717-0700-local1 \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  o2-sim-digitizer-workflow -b --run --onlyDet ITS,MFT \
  < /dev/null > digi.log 2>&1
```

This is a **batch-only** concern. Do not close stdin for interactive or
piped commands (`reader-workflow | tracker-workflow -b --run`) -- only the
last command in such a pipe runs with `-b --run` at all, and the upstream
end must stay a bare, non-batch invocation so it dumps its spec for
merging. `simulation-preflight.zsh` deliberately never closes stdin itself
for this reason; add `< /dev/null` at the call site of the specific
standalone command that needs it.

## Building the DPL plugin libraries

`libO2FrameworkAnalysisSupport` and `libO2FrameworkCCDBSupport` are loaded
via `dlopen` at runtime and are **not** linked dependencies of the
reconstruction/simulation executables, so a Ninja build that only targets
those executables silently omits them -- which breaks CCDB condition
fetching, not just logging. Build them explicitly whenever replay
readiness matters. The shared-library suffix is platform-dependent --
`.dylib` on macOS/Darwin (this environment), `.so` on Linux (most CI/cluster
environments this skill also has to work in); pick the one matching
`uname -s`, or use the plain CMake/Ninja **target name** (no suffix), which
resolves to the correct artifact on either platform:

```bash
# Darwin
O2_BUILD_DIR=/path/to/build \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  ninja -C /path/to/build \
    stage/lib/libO2FrameworkAnalysisSupport.dylib \
    stage/lib/libO2FrameworkCCDBSupport.dylib

# Linux
O2_BUILD_DIR=/path/to/build \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  ninja -C /path/to/build \
    stage/lib/libO2FrameworkAnalysisSupport.so \
    stage/lib/libO2FrameworkCCDBSupport.so

# portable: CMake target names, no platform-specific suffix
O2_BUILD_DIR=/path/to/build \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  ninja -C /path/to/build \
    O2lib-FrameworkAnalysisSupport \
    O2lib-FrameworkCCDBSupport
```

## Durable artifacts

Never store ROOT fixtures, generated workflow outputs, or reusable builds
inside this skill or inside a Git worktree. Keep them as separate external
roots, chosen via arguments or environment variables at call time -- not
hardcoded into any script here, and not tied to one developer's home
directory layout:

- **Reusable worktree builds** -- pass an explicit `--build`/`O2_BUILD_DIR`
  outside the source tree (e.g. a sibling `O2-worktree-builds/<name>`
  directory next to your `O2-worktrees/` checkouts).
- **Durable validation fixtures** -- a package-qualified directory name,
  since a fixture generated under one aliBuild O2 package is not
  interchangeable with one generated under another even with identical
  run/seed/timestamp/modules. For example
  `pp-20ev-run303000-seed20260716-daily20260717` names the run, seed, event
  count, and generating package explicitly; it must never be confused with
  or overwrite a differently-provenanced fixture (including an older,
  differently-packaged one covering the same nominal run/seed). Promote a
  fixture into its durable location only after every required file is
  verified non-empty and checksummed, and make it read-only afterward.
- **Per-run validation/replay artifacts** -- a fresh directory per run under
  its own external root, never reused across runs and never written inside
  the fixture directory that fed it.

Run simulations and replays from an explicit artifact directory outside the
source tree. Never let generated ROOT, topology, CCDB-cache, or log files
appear in the Git worktree.

## Known incidents (read before debugging something that looks like these)

- **AliEn token required for canonical CCDB/Grid access.** A plain CCDB
  HTTP read (e.g. `GLO/Config/GRPMagField`) can succeed with no token at
  all, which is misleading: `CcdbApi::initTGrid()` still requires a valid
  AliEn token for reads/writes against `alice-ccdb.cern.ch` more broadly
  (e.g. `ITS/Calib/Align`), and fails fatally without one unless
  `ALICEO2_CCDB_NOTOKENCHECK=1` is set. That variable is a debug bypass,
  not a fix -- it silently degrades alignment objects to empty rather than
  fetching them. Run `simulation-preflight.zsh` (or `alien-token-info`
  directly) before assuming a CCDB-adjacent failure is a code or package
  problem.
- **One specific package's Geant4 physics-list construction is broken,
  independently of the token issue above.** The `run3-local1` O2 package
  observed in this environment crashes (SIGSEGV) inside Geant4's
  `FTFP_BERT_EMV` physics-list construction during `o2-sim`, reproduced
  deterministically across multiple runs, with and without a valid token,
  with default and enlarged stack size. This is unresolved and specific to
  that one package -- `daily-20260717-0700-local1` runs the identical
  simulation successfully. Do not assume switching packages "fixes" a
  Geant4 crash in general, and do not assume a working package elsewhere
  means every package works.
- **Backgrounded standalone digitization blocks on stdin**, a third,
  unrelated issue: a DPL driver invoked standalone under some backgrounded
  execution contexts hangs indefinitely reading stdin rather than detecting
  no pipe is connected. Fixed by closing stdin explicitly (see above), and
  by that alone -- it does not explain, and did not contribute to, the
  Geant4 crash above; they were isolated independently by reproducing each
  with the other's fix already in place.

## Guardrails

- Avoid interactive `alienv enter`; use the runner.
- Keep the package (requested and resolved), architecture, build directory,
  and output directory explicit in validation reports.
- Do not install or download a new O2 stack when discovery fails unless the
  user authorizes it.
- Treat missing optional CMake dependencies separately from an environment
  initialization failure.
- Never delete or reset a build directory automatically to work around a
  mismatched cache; `configure-worktree-build.zsh` refuses instead.

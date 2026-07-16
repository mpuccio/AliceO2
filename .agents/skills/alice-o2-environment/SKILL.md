---
name: alice-o2-environment
description: Run ALICE O2 builds, tests, reconstruction, simulation, and command-line tools inside the repository's aliBuild environment. Use whenever a command needs O2, ROOT, Ninja, detector workflows, staged branch binaries, or libraries that are unavailable in the agent's default shell.
---

# ALICE O2 environment

Run environment-dependent commands through `scripts/run-in-o2-env.zsh`. The
runner discovers the installation from `ALIBUILD_WORK_DIR`, sources O2, and
executes the command in that same process. Do not source the profile in one
tool call and run the command in another; environment changes do not persist.

## Check the environment

```bash
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh --check
```

If this fails, report its diagnostics. Do not hardcode a developer home
directory. Require `ALIBUILD_WORK_DIR` to be exported by the process that
started the agent.

The runner derives `ALIBUILD_ARCH_PREFIX` with `aliBuild architecture` when it
is unset. It uses `O2_PACKAGE=latest-run3-o2` by default. Override that variable
when the task requires another aliBuild O2 package. The runner also bridges
`ALIBUILD_WORK_DIR` to the legacy `WORK_DIR` name used inside generated O2
profiles; callers should not need to set both.

## Build and test a worktree

Set `O2_BUILD_DIR` to the configured build directory. The runner prepends its
`stage/bin` and host library directory so branch-built executables take
precedence over installed ones.

```bash
O2_BUILD_DIR=/path/to/build \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  ninja -C /path/to/build O2lib-ITSMFTTracking
```

```bash
O2_BUILD_DIR=/path/to/build \
  .agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
  ctest --test-dir /path/to/build --output-on-failure -R 'test-name'
```

Prefer a build configured from the active worktree. Verify the source directory
in `CMakeCache.txt` before using an unfamiliar build directory.

## Run installed and staged tools

Omit `O2_BUILD_DIR` when the installed aliBuild executable is intentionally
required, for example fixture simulation. Set it when reconstruction must use
the active branch's staged binary.

Run simulations from an explicit artifact directory outside the source tree,
normally under `/private/tmp` or `/tmp`. Never let generated ROOT, topology,
CCDB-cache, or log files appear in the Git worktree.

## Guardrails

- Avoid interactive `alienv enter`; use the runner.
- Keep the package, architecture, build directory, and output directory
  explicit in validation reports.
- Do not install or download a new O2 stack when discovery fails unless the
  user authorizes it.
- Treat missing optional CMake dependencies separately from an environment
  initialization failure.

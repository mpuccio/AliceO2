# ITS/MFT Tracking Refactoring: Agent Coordination Manifest

Status: Active planning document  
Integration owner: Codex architecture/integration agent  
Integration branch: `codex/itsmft-integration`  
Architecture: [Architecture.md](Architecture.md)

## 1. Objective

Coordinate multiple coding agents working on the shared ITS/MFT CA tracking refactoring without duplicating work, changing incompatible interfaces, or creating unreviewable cross-cutting commits.

This document defines roles, file ownership, work waves, integration gates, communication rules, and the required handoff format. It is operational rather than architectural; architectural choices belong in `Architecture.md` and the decision log below.

## 2. Coordination principles

1. Contract first: public interfaces are agreed before dependent agents implement against them.
2. One active owner per production file in each wave.
3. One worktree and branch per agent; agents do not share a mutable checkout.
4. The integration owner alone updates the integration branch.
5. Commits are small, buildable, and ordered by dependency.
6. Compatibility is implemented with aliases/adapters, not copied code.
7. Physics retuning is excluded from structural commits.
8. Existing user changes are preserved and never overwritten by an agent.
9. A gate closes only when its acceptance checks are recorded as passing.
10. Architecture changes require a decision-log entry before implementation spreads to other branches.

## 3. Agent roster

Names are stable role identifiers. The concrete tool or human assigned to a role may change.

| Role | Suggested branch | Responsibilities | Starts |
|---|---|---|---|
| Architecture/integration | `codex/itsmft-integration` | Contracts, dependency boundaries, reviews, integration, combined validation | Wave 0 |
| Baseline/validation | `codex/itsmft-baseline-tests` | Characterization tests, unit tests, benchmarks, regression reports | Wave 0 |
| Core layout/data model | `codex/itsmft-core-layout` | Surface IDs/descriptors, masks, sparse topology, shared primitive ownership | Wave 0/1 |
| TimeFrame/input | `codex/itsmft-timeframe-input` | Multi-source loading, cluster references, timing, normalized measurements | Wave 1 |
| CA policies | `codex/itsmft-ca-policies` | Policy boundary, barrel/disk operations, detector-free orchestration | Wave 2 |
| ITS migration | `codex/itsmft-its-migration` | ITS production adapter, CPU parity, duplicate removal, GPU impact report | Wave 2/3 |
| MFT adapter | `codex/itsmft-mft-adapter` | MFT output/refit adapter, workflow migration, `MFTCATrack` removal | Wave 3 |

The architecture/integration role should remain lightly loaded with direct implementation so it can review and unblock other roles.

## 4. Ownership map

Ownership applies only while a role's wave is active. Exact new filenames are assigned in the kickoff issue/message.

| Area | Primary owner | Other agents may |
|---|---|---|
| `LayerMask`, surface IDs and descriptors | Core layout | Review and request interface changes |
| `TrackingTopology` and layout builders | Core layout | Add tests through validation owner |
| Common cluster/tracklet/cell primitives | Core layout | Add compatibility aliases after approval |
| `TimeFrame`, normalized measurements | TimeFrame/input | Consume established core interfaces |
| Cluster loading, dictionaries, labels | TimeFrame/input | Add detector adapters in owned adapter files |
| ROF timing and overlap integration | TimeFrame/input | Add tests through validation owner |
| `TrackerTraits` orchestration | CA policies | Add instrumentation through agreed hooks |
| Projection/index/fitting policies | CA policies | ITS/MFT agents may implement adapter-specific backends |
| ITS production workflow and compatibility | ITS migration | Modify only ITS-owned adapter/workflow files |
| MFT workflow and output conversion | MFT adapter | Modify only MFT-owned adapter/workflow files |
| Tests, fixtures, regression scripts | Baseline/validation | Production owners may add local unit tests after coordination |
| Cross-target CMake dependency changes | Architecture/integration | Propose a patch or dependency requirement |
| Architecture and decision log | Architecture/integration | Submit proposed amendments |

When a task requires edits in another role's active files, the agent must stop and send one of:

- An interface request.
- A minimal proposed diff for the owner to apply.
- A request to transfer ownership for explicitly listed files.

Silent cross-ownership edits are not accepted for integration.

## 5. Worktree and branch policy

Each role receives a dedicated Git worktree based on the integration branch. Example layout outside the primary checkout:

```text
O2-worktrees/
  baseline-tests/
  core-layout/
  timeframe-input/
  ca-policies/
  its-migration/
  mft-adapter/
```

Branch names use the `codex/` prefix regardless of whether Codex or Claude Code performs the work, unless repository policy requires otherwise. The branch identifies the workstream, not the tool.

Agents must:

- Check `git status` before starting and before handing off.
- Never reset, restore, or overwrite changes they did not create.
- Rebase or merge only when instructed by the integration owner.
- Avoid generated formatting changes outside owned files.
- Avoid amending commits after their hashes have been handed off unless requested.

## 6. Wave plan

### Wave 0: contracts and baseline

Active roles:

- Architecture/integration.
- Baseline/validation.
- Core layout feasibility.

Deliverables:

- Architecture RFC accepted for implementation.
- Decision log initialized.
- Current ITS and MFT build/run commands recorded.
- Initial characterization tests or a documented fixture plan.
- Proposed public types for surface IDs, descriptors, masks, and topology views.

Exit criteria:

- No unresolved decision blocks Wave 1 interfaces.
- Baseline failures are classified as pre-existing or introduced.

### Wave 1: foundations and TimeFrame

Active roles:

- Core layout/data model.
- TimeFrame/input.
- Baseline/validation.
- Architecture/integration.

Deliverables:

- Surface/layout model.
- Mask supporting the chosen maximum surface count.
- Explicit sparse topology.
- Normalized cluster and qualified reference proposal/implementation.
- Single-detector loading retained through adapters.
- Unit tests for the new foundations.

Exit criteria:

- ITS-only and MFT-only layouts reproduce current reachability.
- A 17-surface disconnected layout constructs successfully.
- Current single-detector TimeFrames still load and expose equivalent clusters.

### Wave 2: common CA orchestration

Active roles:

- CA policies.
- ITS migration preparation.
- Baseline/validation.
- Architecture/integration.

Deliverables:

- Transition-policy interface.
- Cylinder-cylinder and disk-disk policies.
- Detector-wide branches removed from main CA loops.
- Existing refitters available behind adapters.
- ITS and MFT algorithm regression report.

Exit criteria:

- Both detector layouts run through the same traversal implementation.
- Differences from baseline are understood and accepted.

### Wave 3: production workflow migration

Active roles:

- ITS migration.
- MFT adapter.
- Baseline/validation.
- Architecture/integration.

Deliverables:

- ITS CPU workflow uses the common core.
- MFT CA workflow uses the common core.
- External outputs and configurable parameters remain compatible.
- Duplicate CPU implementation removal plan executed.
- CMake dependency direction corrected.

Exit criteria:

- Production workflows build and pass agreed regression tests.
- No shared core public header depends on detector tracking implementations.

### Wave 4: combined disconnected TimeFrame

Deliverables:

- ITS and MFT sources loaded into one TimeFrame.
- Independent timing handled in one common time domain.
- Two disconnected topology components processed in one invocation.
- Results routed to the correct detector output adapters.

Exit criteria:

- Combined results agree with running the two layouts separately within tolerance.
- Memory and runtime changes are recorded.

### Wave 5: mixed-surface state research

This wave starts only after a separate track-state RFC is approved. It is not implicitly authorized by earlier waves.

## 7. Commit protocol

Each commit should do one conceptual job and remain reviewable. Preferred examples:

```text
ITSMFT: add strongly typed surface identifiers
ITSMFT: widen tracking surface mask to 32 bits
ITSMFT: add sparse transition topology
ITSMFT: qualify cluster references by input source
ITSMFT: add disk-disk transition policy
ITS: adapt CPU tracker to common TimeFrame
MFT: convert common tracks to TrackMFT output
```

Do not combine:

- Structural refactoring and cut retuning.
- Data-model changes and workflow output changes.
- Broad renaming and behavioral changes.
- Multiple ownership areas without prior agreement.

Every production commit must state tests run in its commit message or handoff report.

## 8. Agent kickoff message template

The integration owner sends each agent a bounded task using this template:

```text
Role: <stable role name>
Branch/worktree: <branch and absolute worktree path>
Base commit: <integration commit>

Objective:
<one concrete outcome>

Read first:
- Detectors/ITSMFT/common/tracking/doc/Architecture.md
- Detectors/ITSMFT/common/tracking/doc/AgentCoordination.md
- <role-specific files>

Owned files:
- <explicit paths or path patterns>

Do not modify:
- <other active ownership areas>

Required deliverables:
- <code/tests/design note>

Acceptance checks:
- <commands or observable conditions>

Handoff:
Return commit hashes, changed interfaces, tests/results, assumptions,
known limitations, and follow-up requests. Do not merge or push unless asked.
```

Objectives must be achievable in one bounded work cycle. “Refactor TimeFrame” is too broad; “introduce detector-qualified `ClusterRef` and adapt accessors without changing loading” is appropriate.

## 9. Handoff report template

Every agent returns:

```text
Role:
Branch:
Base commit:
Commits:

Outcome:

Interfaces added or changed:

Files modified:

Tests and results:

Behavioral or physics differences:

Assumptions:

Known limitations:

Requested decisions or follow-up work:

Worktree status:
```

If no commit was created, the report must say why and list any working-tree modifications.

## 10. Integration procedure

For every handoff, the integration owner:

1. Confirms the branch is based on the expected integration commit.
2. Reviews ownership compliance and unrelated changes.
3. Reviews public interfaces against `Architecture.md`.
4. Runs focused tests from the handoff.
5. Cherry-picks or merges the smallest dependency-complete commit set.
6. Runs integration-level build/tests.
7. Updates the gate status and decision log.
8. Notifies dependent agents of the new integration commit.

If substantial rework is needed, the original role receives a concrete follow-up task rather than the integration owner silently rewriting the subsystem.

## 11. Cross-agent communication

Agents communicate through short written messages tied to stable role names. Messages should contain one of:

- `DECISION NEEDED`: options, recommendation, and blocking impact.
- `INTERFACE REQUEST`: current interface, requested change, and consumer.
- `OWNERSHIP REQUEST`: exact files and duration.
- `BASE UPDATE`: new integration commit and required rebase/merge instruction.
- `HANDOFF`: the report from Section 9.

Agents should continue independent work when a question is non-blocking. They must stop before making an assumption that changes a public contract or another role's scope.

## 12. Testing and acceptance ownership

Production owners write focused unit tests for behavior introduced in their component. The baseline/validation role owns:

- Shared fixtures.
- End-to-end regression comparison.
- Physics-output tolerance definitions.
- Performance and memory reports.
- Classification of pre-existing failures.

An agent may not declare a gate passed based only on compilation. Gate status is assigned by the integration owner after validation results are available.

## 13. Decision log

The integration owner records accepted decisions here or links a dedicated follow-up RFC.

| ID | Status | Decision | Rationale | Affected waves |
|---|---|---|---|---|
| D001 | Accepted | Use explicit detector-qualified `SurfaceId` metadata; never infer detector from surface count | Required for arbitrary layouts | 1+ |
| D002 | Accepted | Use a 32-bit `SurfaceMask` with at most 32 global surfaces | Covers the 17-surface target with 15 spare bits while preserving a compact GPU-friendly POD | 1 |
| D003 | Accepted | Use an explicit sparse directed topology with 16-bit strong IDs and CSR lookup views | Avoids combinatorial storage and numeric-order assumptions while remaining device-friendly | 1+ |
| D004 | Accepted | First combined layout uses disconnected ITS and MFT subgraphs | Decouples container integration from mixed-state propagation | 4 |
| D005 | Accepted | Use a 72-byte trivially-copyable `SurfaceMeasurement` with detector-qualified sensor/surface identity and `{ClusterSourceId, external index}` cluster identity | Gives CPU/device code one normalized measurement while keeping labels and detector decoding at the adapters; disk covariance must be decoded in explicit x/y surface axes rather than copied from synthetic MFT `TrackingFrameInfo` semantics | 1 |
| D006 | Accepted | Compose a common tracking TimeFrame, separate `ITSVertexingState`, optional non-owning `VertexConstraintView`, and framework/device allocation state | Keeps ITS-only products out of MFT/common ownership and permits incremental GPU migration through POD views and temporary compatibility facades | 1/3 |
| D007 | Accepted | Carry one validated policy tag per transition; derive the state family and cell policy, dispatch only outside hot loops, and launch one specialized GPU kernel per active `(stage, family)` | Supports data-driven disconnected layouts without virtual/device-indirect calls or per-candidate branching; parameter and typed-state view ABIs remain separate implementation contracts | 2/3 |
| D008 | Deferred | Select mixed cylinder-disk track state | Requires a dedicated RFC after common production migration | 5 |

Status values are `Open`, `Proposed`, `Accepted`, `Superseded`, or `Deferred`. Only the integration owner marks a decision accepted after maintainer agreement.

## 14. Gate status

| Gate | Status | Evidence |
|---|---|---|
| Gate 0: baseline | Complete | Legacy characterization tests pass; the committed 20-event ITS/MFT fixture protocol records bit-identical single-thread replay metrics, physics definitions, CCDB provenance, wall time, and peak RSS |
| Gate 1: foundations | Complete | Layout, normalized measurements, geometry-backed ITS/MFT decode adapters, Stage-A transition policy validation, transactional multi-source loading with TF-relative timing, and the normalized-loading compatibility boundary in the common single-detector TimeFrame are integrated; ITS and MFT parity tests cover legacy backfill and the accepted real-geometry validation covers 7,057 fixture clusters |
| Gate 2: common CA traversal | In progress | D007, device-compatible cylinder/disk policy parameter and state-family traits, host outer-family dispatch, and the first host CPU policy operation (`cellsAreCompatible`) are integrated. The legacy duplicate formula and the detector branch in that hot loop are removed; associating the sparse policy-tagged topology with the production TimeFrame, completing topology-driven TrackerTraits/CATracker orchestration, and demonstrating ITS/MFT algorithm parity remain. |
| Gate 3: production migration | Blocked by Gate 2 | |
| Gate 4: combined disconnected tracking | Blocked by Gate 3 | |
| Gate 5: mixed-surface tracking | Deferred | Requires track-state RFC |

Gate 0 acceptance is a regression and performance baseline, not a claim of
offline or scaling coverage. Reproduction currently depends on the recorded
CCDB objects remaining network-accessible; the 12-thread replay is
characterization rather than a meaningful scaling test at this fixture size;
and the documented MFT reconstructability denominator is intentionally simpler
than the geometry-aware ITS definition.

## 15. Kickoff record and next task

The first three bounded assignments should be:

1. **Baseline/validation**: inventory existing relevant tests and build targets; propose and implement the smallest unit tests that characterize current 7- and 10-surface topology/mask behavior. Do not change production code.
2. **Core layout feasibility**: propose concrete identifier widths, mask width, sparse topology storage, and device view layouts. Implement only after D002 is accepted.
3. **Architecture/integration**: review the current uncommitted `DetectorTraits.cxx` seed-conversion fix, establish the integration branch/worktrees, and resolve D002, D005, and D006 with maintainers.

The first follow-on implementation was the normalized measurement primitives and ITS/MFT adapter parity fixtures described by D005. It did not change TimeFrame ownership, workflows, kernels, or physics behavior. Conversion helpers receive explicit source, sensor, surface, ROF, shape, and covariance-axis information; the MFT path does not blindly project the legacy synthetic `TrackingFrameInfo` into disk coordinates.

Gate 1 closed after integration of the common TimeFrame normalized-loading
compatibility boundary. Existing ITS-only and MFT-only inputs now expose
equivalent clusters through the normalized owner/view while backfilling the
legacy common TimeFrame structures without changing production workflows, CA
orchestration, GPU kernels, vertexing ownership, or physics behavior. The
single-detector boundary tests and the separately accepted real-geometry
validation provide the evidence needed to start Gate 2.

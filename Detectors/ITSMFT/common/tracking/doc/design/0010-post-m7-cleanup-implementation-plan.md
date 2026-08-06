# Post-M7 cleanup implementation plan

- Status: L1 complete; L2 remains the next authorized slice
- Date: 2026-08-06
- Architecture source: [post-M7 intentionality audit](0009-post-m7-intentionality-cleanup-audit.md)
- Implementation owner: Luna agents under maintainer review
- Physics status: behavior-preserving structural cleanup only

## 1. Purpose

This document turns design note 0009 into a dependency-ordered execution
contract. It does not authorize implementation by itself and contains no
production change.

The target architecture has one reusable generic tracking entity:

```text
DPL/application event owner
  owns raw ROFs, timing construction, publication validity, typed output
  owns the lifetime of one configured TimeFrame
       |
       +-- Loader::load(TimeFrame&, normalized inputs/runtime ROF views)
       |
       +-- Tracker::run(TimeFrame&, TrackerTraits&)
                                      |
                                      +-- CPU or GPU kernel implementation
```

`TimeFrame` owns the data on which the components act:

- one immutable `std::vector<SurfaceGraph>` configured once;
- per-iteration parameters and source-qualified graph partitions;
- allocator/pool identity and allocator-backed CA workspaces;
- normalized current-event measurements and runtime timing/mask views; and
- generic `CommonTrack`, cluster-reference, label, and support results.

`Loader` and `Tracker` do not own those data. `TrackerTraits` contains the
target-architecture implementations of the kernels invoked by `Tracker`.
There is no public `TrackingPlan` class: vector order is iteration order, and
private TimeFrame records may keep each graph paired with its parameters and
partitions without creating another public abstraction.

## 2. Non-negotiable boundaries

Every implementation slice must preserve all of the following:

- no change to Propagator/refit arithmetic, covariance sanitization, material,
  MinPt, chi2 or acceptance cuts, holes, road/candidate ordering, or policy
  grouping;
- no change to raw-ROF ownership, timing values, workflow defaults,
  publication clocks/validity, writer content, sidecars, or output formats;
- no detector identity, source-position convention, DPL/workflow/writer type,
  typed output track, or public transition taxonomy in generic core APIs;
- `StateFamily` remains representation metadata only;
- fixed device capacities remain fixed where their ABI requires them;
- frozen `Detectors/ITSMFT/ITS/tracking/` code remains untouched;
- no runtime-plan wrapper, coordinator, manager, service locator, callback
  framework, or central detector dispatcher is introduced; and
- no triplet R&D, mixed-cylinder/disk enablement, ALICE3 integration, or
  physics-sign-off change is bundled with cleanup.

The named user stash containing `TripletTrackingRAndD.md` must never be
restored, modified, dropped, or committed by these slices.

## 3. Luna execution protocol

No repository-specific Luna protocol was found, so this section is the
campaign protocol.

1. Run slices sequentially. Each branch starts from the current integration
   head containing every previously accepted slice, never from another
   unaccepted Luna branch.
2. Use the worktree and durable build assigned by the maintainer. Confirm a
   clean worktree and that `CMAKE_HOME_DIRECTORY` names that worktree before
   editing. Reconfigure the assigned build only when required; do not create a
   substitute build to avoid a failure.
3. Treat each subsection in Section 6 as a maximum scope. A Luna agent may
   make the slice smaller, but may not absorb work from the next slice.
4. Preserve user changes and stashes. Do not clean, reset, restore, or rewrite
   unrelated work.
5. Prefer deletion and direct ownership. Any temporary bridge must be the one
   explicitly listed for that slice and must have a named deletion slice.
6. Use separate reviewable commits where a slice changes production, tests,
   and documentation:
   - production migration/deletion;
   - tests and ownership/source guards; and
   - validation record/status update.
7. Stop on a physics delta, an ownership ambiguity that changes this plan, an
   unavailable required replay dependency, or a failure that cannot be fixed
   without widening scope. Diagnose and record the blocker; do not improvise a
   new abstraction.
8. A slice is not complete until its deletion criterion is true. Compiling a
   forwarding bridge is not completion.

### 3.1 Required handoff record

Every Luna agent must leave a concise handoff in the slice's validation note
or commit message containing:

| Field | Required content |
|---|---|
| Provenance | base and final commits, branch, package, source/build paths |
| Scope | files changed and files intentionally left for the next slice |
| Ownership result | owner before and after; temporary bridge, if any |
| Deletions | removed types/files/APIs and repository-search evidence |
| Validation | exact build, test, formatting, checksum, and replay commands/results |
| Output parity | standalone hashes, combined equality, writer comparison |
| Device status | real CUDA/HIP result, or explicit toolchain absence |
| Exceptions | only previously approved exclusions; no newly invented one |
| User state | confirmation that unrelated changes and named stash are untouched |

## 4. Common validation contract

Call the following production-slice contract **R**:

1. build all affected common libraries/tests and ITS, MFT, and combined
   workflow targets in the assigned durable build;
2. run the complete selector serially:

   ```sh
   ctest --test-dir <durable-build> -L itsmft --output-on-failure -j1
   ```

   Every registered test must execute and pass; `Not Run` is a failure;
3. run the focused tests and source/dependency guards named by the slice;
4. run `git diff --check <slice-base> HEAD` and
   `git clang-format --diff <slice-base> HEAD`;
5. verify all 43 fixture checksums before and after replay;
6. require standalone native results:
   - ITS: 212 tracks / `46913a67a7e2fe7462e29df0db264fa8`;
   - MFT: 68 tracks / `8106b08571ca593c6b76ff72b761a680`;
7. require each combined leg to be byte-identical to its standalone product;
8. compare initialized writer content, `CommonTrack`s, cluster references,
   ROFs, labels, and sidecars with the accepted parent. Exclude only the known
   undefined `MFTTrack.mInvQPtSeed` byte artifact; and
9. perform a real device build only when the pinned CUDA/HIP toolchain exists.
   Never infer or fabricate a GPU result.

Documentation-only status updates run Markdown link/heading/table validation
and `git diff --check`; they do not rerun R.

## 5. Dependency order

```text
L0 dead typed path
 |
L1 canonical Tracker API
 |
L2 dead forwarding/reset residue
 |
L3 SurfaceGraph vocabulary and single graph representation
 |
L4 TimeFrame configuration ownership and per-iteration partitions
 |
L5 TimeFrame workspace/reset ownership
 |
L6 direct Loader -> TimeFrame transaction
 |
L7 workflow timing/publication ownership extraction
 |
L8 Tracker orchestration; Engine/Participant retirement
 |
L9 TrackingInterface/SurfacePlanTrackingParticipant retirement
 |
L10 operation seam and public compatibility tail
```

L0-L2 are low-risk deletion/API cleanup. L3-L9 are sequential because they
mutate overlapping ownership and lifecycle paths. They must not be assigned
concurrently. L10 is intentionally last: narrowing refit/publication seams
before the entity/component lifecycle is stable would create another bridge.

## 6. Implementation slices

### L0 — delete the dead typed-MFT refit/export path

**Suggested branch:** `luna/itsmft-postm7-l0-dead-mft-refit`

**Outcome:** delete the unused common typed-MFT refit/export files; migrate
the meaningful fixture to the already-live generic native-refit result path.

**Primary scope:**

- the three unused common typed-MFT refit/export files
- common tracking CMake source lists
- `testMFTNormalizedRefit.cxx`
- source guards and stale comments that whitelist these types

**Temporary bridge:** none.

**Deletion criterion:** repository production/tests contain no common typed
MFT refit/export symbol, typed refit overload, CMake entry, or replacement
wrapper. MFT publication remains unchanged. **Completed in L0:** the focused
guard and repository search prove this absence.

**Gate:** R plus focused native-refit failure/result tests. This is the first
Luna implementation task and must not include any lifecycle cleanup.

### L1 — make `Tracker.{h,cxx}` the canonical implementation files

**Suggested branch:** `luna/itsmft-postm7-l1-canonical-tracker`

**Outcome:** move the existing non-templated `Tracker` declaration and
definition from stale `CATracker` file names to `Tracker.{h,cxx}` and delete
the reversed forwarding relationship.

**Primary scope:** `CATracker.{h,cxx}`, `Tracker.h`, CMake lists, direct
includes, public-header guards, and comments/tests that name the old file.

**Temporary bridge:** an include-forwarder is permitted only within the
production commit while callers are being changed; it must be deleted before
the slice's final commit.

**Deletion criterion:** one canonical header/source pair; no `CATracker`
filename, include, alias, declaration, or CMake entry remains in common
tracking.

**Gate:** R plus a public-header compile/dependency guard.

**Status: complete (2026-08-06).** `Tracker.h` and `Tracker.cxx` are now the
single implementation pair. Direct ITS/MFT adapter includes, CMake sources,
dependency guards, and the failure-contract fixture use the canonical names;
the old files and forwarding relationship are deleted. The focused guard
checks the old filenames, includes, standalone identifier, and CMake entries.
The validation record is [L1 canonical Tracker validation](../validation/l1-canonical-tracker.md).

### L2 — remove dead scratch and fixed-source forwarding residue

**Suggested branch:** `luna/itsmft-postm7-l2-dead-forwarders`

**Outcome:** delete only already-proven dead or forwarding-only state before
the ownership migration obscures it.

**Primary actions:**

- delete `SurfaceTrackingScratch::mPValphaX` and its allocator/reset/swap
  bookkeeping;
- retain one scratch reset spelling and delete `resetScratch()` forwarding;
- delete test-only fixed ITS/MFT loader/reset wrappers
  `loadITSAndMFT()`/`resetITSAndMFTEvent()` after migrating tests to the
  generic load/reset behavior; and
- correct obsolete comments/test names without changing assertions.

**Temporary bridge:** none. Separate commits are preferred for the dead
member, reset alias, and fixed-source helpers.

**Deletion criterion:** repository searches find none of the removed member
or helper names, and no source-0/1 helper survives in common loading.

**Gate:** R plus allocator/swap/reset, atomic-load retry, dropped-TF, and
combined-source-isolation tests.

### L3 — consolidate the public graph representation

**Suggested branch:** `luna/itsmft-postm7-l3-surface-graph`

**Outcome:** one public `SurfaceGraph`/`SurfaceGraphView` concept replaces the
overlap among `DetectorLayout`, `DetectorLayoutView`, and public sparse
topology ownership. `DetectorLayoutConfigurationKey` is deleted. There is no
`TrackingPlan` class.

**Primary scope:** layout/topology headers and builders, graph fixtures,
standalone/combined graph construction, Tracker/traits graph consumers, and
the final runtime-core guard.

**Temporary bridge:** existing owners may hold `std::vector<SurfaceGraph>`
until L4 moves that vector into `TimeFrame`. No alias preserving the old type
names and no second graph representation may survive L3.

**Deletion criterion:**

- one owning graph and one device-facing graph view;
- sparse adjacency is at most private graph storage, not a competing public
  model;
- no stale configuration-key/equality cache or always-true `rebuilt` result;
- combined construction owns one graph, not copied ITS/MFT graphs; and
- vector order alone defines iteration order.

**Gate:** R plus non-contiguous `SurfaceId`, sparse adjacency, holes/seeding,
graph validation, combined partition, and device-view layout tests.

### L4 — move static tracking configuration into `TimeFrame`

**Suggested branch:** `luna/itsmft-postm7-l4-timeframe-configuration`

**Outcome:** the reusable `TimeFrame` owns its graph vector, per-iteration
parameters, mandatory source-qualified graph partitions/bindings, allocator
identity, and workspace capacity requirements through one fallible atomic
configuration boundary.

**Initialization contract:** configuration occurs once before event loading.
The frame cannot be loaded or tracked while partially configured. Graph count,
parameter count, and partition validity are checked together. Every iteration
has its own partition data; iteration-0 binding reuse, numeric traversal, and
identity/source-0 fallback mapping are deleted.

**Temporary bridge:** participant/interface wrappers may still own the
physical scratch containers, but they must borrow all sizing/configuration
from `TimeFrame`. L5 deletes this split. Tracker-owned parameter/pool handles
must be gone by L4 exit.

**Deletion criterion:** one configured frame is the only static configuration
authority; no Tracker, participant, interface, or layout-set object owns a
second graph/parameter/partition/pool authority.

**Gate:** R plus construction failure/retry, allocator identity, graph-count
mismatch, per-iteration differing-hole/start-mask, non-identity partition,
and partial-configuration rejection tests.

### L5 — internalize workspace and generic reset in `TimeFrame`

**Suggested branch:** `luna/itsmft-postm7-l5-timeframe-workspace`

**Outcome:** mutable CA workspace becomes TimeFrame-owned entity state.
`SurfaceTrackingScratch` is either merged into `TimeFrame` or retained only as
a private TimeFrame partition-workspace implementation detail; it is not a
separately configured public entity.

One TimeFrame event reset clears normalized event data, timing/mask views,
workspace contents, generic results, labels, and event-derived support state.
It preserves graphs, parameters, partitions, allocator identity, and reserved
capacity. Loader staging failure remains non-mutating.

**Temporary bridge:** wrapper methods may forward to the TimeFrame reset only
until L8/L9 deletes the wrappers. They may not clear independent generic
state.

**Deletion criterion:** no participant/interface-owned generic scratch, no
independent scratch configuration/reset authority, and exactly one generic
event-state reset primitive.

**Gate:** R plus allocator reuse, successful replacement, recoverable drop,
structural failure, exception, stale timing/mask/result, and reset-count tests.

### L6 — make Loader act directly on `TimeFrame`

**Suggested branch:** `luna/itsmft-postm7-l6-direct-loader`

**Outcome:** expose one atomic `Loader::load(TimeFrame&, inputs)` composition
(using the existing loader implementation/name unless a direct rename is
separately justified). Delete the one-implementation virtual load-target
hierarchy, participant load-target members, public scratch targeting, and
friendship-only forwarding.

The workflow/adapter constructs raw-ROF-derived runtime timing/mask views;
Loader consumes those views and normalized sources but never owns raw ROFs or
publication state.

**Temporary bridge:** none after the final commit. Do not add a generic input
manager or callback-based loader.

**Deletion criterion:** every production load reaches TimeFrame directly and
atomically; partial staging cannot mutate the live frame; no `LoadTarget`
object survives.

**Gate:** R plus multi-source atomicity, allocator failure, partial decode,
retry/replacement, source qualification, and timing-view lifetime tests.

### L7 — place timing and publication compatibility at workflow edges

**Suggested branch:** `luna/itsmft-postm7-l7-application-ownership`

**Outcome:** fixed ITS/MFT ROF tables, timing construction/backing storage,
publication adapters, clocks/validity, and typed sidecars are owned by the
ITS/MFT/combined application or DPL task, not participant/interface/core
objects. TimeFrame borrows only event-lifetime runtime views and owns only
generic results.

**Temporary bridge:** existing participant/interface methods may accept the
workflow-owned objects until L8/L9 removes those wrappers. They may not own or
reset them.

**Deletion criterion:** generic TimeFrame/Loader/Tracker/TrackerTraits and
their public headers contain no fixed detector table, detector publication
sidecar, typed output, clock/validity owner, or workflow reset decision.

**Gate:** R plus empty/first/last ROF, diamond timing, mask parity, load
failure, dropped TF, structural failure, replacement, sidecar invalidation,
and writer-publication tests.

### L8 — make Tracker the sole non-owning orchestrator

**Suggested branch:** `luna/itsmft-postm7-l8-retire-engine`

**Outcome:** `Tracker::run(TimeFrame&, TrackerTraits&)` reads the frame-owned
iteration/partition order and sequences initialization, tracklets, cells,
neighbours, roads/refit, generic result commit, and failure handling.
`Tracker` stores no frame configuration or event state. `TrackerTraits`
remains a distinct borrowed CPU/GPU kernel implementation.

Delete `TrackingEngine`, `TrackingParticipant`, `ParticipantId`, and the
combined `SurfacePlanTrackingParticipant` composition after the combined DPL
task invokes Loader + TimeFrame + Tracker directly.

**Temporary bridge:** standalone `TrackingInterface` may call the new Tracker
API until L9. It may not retain a second orchestration/reset path.

**Deletion criterion:** no engine/participant class or renamed schedule
executor remains; one Tracker invocation handles explicit frame partition
order; TimeFrame supplies the only generic reset mechanics.

**Gate:** R plus exact schedule order, non-identity partition order, combined
source isolation, CPU traits substitution, success/recoverable/structural/
exception outcomes, and single-reset tests.

### L9 — retire `TrackingInterface` and complete component composition

**Suggested branch:** `luna/itsmft-postm7-l9-retire-interface`

**Outcome:** standalone ITS and MFT workflows use the same explicit
workflow-owned composition as combined tracking:

```text
workflow owns TimeFrame lifetime + application compatibility
Loader acts on TimeFrame
Tracker acts on TimeFrame through TrackerTraits
workflow publishes successful generic results
```

Delete `TrackingInterface`, remaining `SurfacePlanTrackingParticipant`
residue, aliases/mixins, and the float failure sentinel. Preserve existing
typed workflow API and failure classification at the workflow edge.

**Temporary bridge:** none. Do not replace the interface with a facade,
coordinator, context, service, or manager.

**Deletion criterion:** standalone and combined use one entity/component
composition; no interface/participant lifecycle owner remains; raw input and
publication lifecycle are visibly workflow-owned.

**Gate:** R plus standalone configuration/lifecycle, dropped-TF behavior,
combined parity, sidecar reset, and exact writer tests.

### L10 — narrow the operation seam and close public compatibility residue

**Suggested branch:** `luna/itsmft-postm7-l10-operation-tail`

**Outcome:** move `completeAccepted()` and `resetAdapterState()` out of
`TrackingOperationAdapter` into application publication handling. Retain only
the exact call-scoped refit operation if it still has more than one real
implementation; otherwise delete the seam. Do not create callbacks or central
dispatch.

Separately resolve the remaining proven-dead public compatibility tail from
design note 0009: no-op common device hook, common
`TrackingFrameInfoAdapters` path, old comparison driver, detector-spec header
placement, and stale tests/comments. Each deletion requires its stated
whole-repository or device gate and may be split into smaller Luna branches.

**Temporary bridge:** none beyond a retained, documented, call-scoped refit
operation with concrete CPU/GPU ownership.

**Deletion criterion:** the core operation seam performs one algorithmic
operation only; no publication/reset lifecycle crosses it; every remaining
public compatibility header has a live named owner and caller.

**Gate:** R plus native-refit ordering/holes/rejection tests, ITS shared-status
and MFT compatibility checks, failure-stale-state tests, public-header guards,
and a real device build when available.

## 7. Review checkpoints

Maintainer review is mandatory after L3, L5, L7, and L9:

| Checkpoint | Review question |
|---|---|
| After L3 | Is there exactly one graph representation, with no plan wrapper or detector-shaped topology? |
| After L5 | Is TimeFrame genuinely the sole entity, or did state merely move behind forwarding methods? |
| After L7 | Are timing/publication/raw-ROF responsibilities visibly application-owned without duplicating event state? |
| After L9 | Can every workflow be described only as owner + Loader + TimeFrame + Tracker + traits + publication adapter? |

A checkpoint may reorder or stop later slices, but it must not retroactively
declare a failed deletion criterion acceptable.

## 8. First Luna assignment

The first assignment is L0 only. Its copy-ready intent is:

> Delete the unused common MFT typed refit/export path and migrate its focused
> test to the live generic native-refit result. Do not alter lifecycle,
> ownership, physics, publication, or output. Prove the symbols/files/CMake
> entries are absent, run contract R, commit production, tests/guards, and the
> validation record separately, and leave the named user stash untouched.

L1 must not begin until L0 is integrated and its replay evidence is accepted.

# L3 SurfaceGraph consolidation

Status: implementation and validation record for Gate 4 M7/L3
Branch: `codex/itsmft-postm7-l3-surface-graph`
Package: `daily-20260717-0700-local1`

## Result

L3 consolidates the common public graph representation into:

- `SurfaceGraph`, the host-owned graph containing copied surface descriptors,
  ordered surface IDs, sparse transitions/cells, CSR cell lookup, masks, and
  validation state;
- `SurfaceGraphView`, a standard-layout, trivially-copyable POD whose pointers
  borrow the immutable graph storage for host/device traversal.

`SurfaceGraphBuilder` is the construction helper. It is not a second graph
owner: it validates declarations, constructs one graph, and returns that
graph. The old `DetectorLayout*` and public sparse-topology headers,
configuration-key residue, and forwarding accessors are deleted. Vector order
of `std::vector<SurfaceGraph>` remains the iteration order until L4.

The combined ITS+MFT composition now builds one graph and derives both
source-qualified bindings from it. The graph remains disconnected at the
detector boundary unless an explicit declaration supplies an edge.

## Ownership boundary with L4

`TimeFrame` is passive reusable storage. L3 does not add graph-building policy
or graph mutation to it. Until L4, existing application/interface owners hold
the graph vector and pass it to the tracker components.

The L4 direction is deliberately narrow:

```text
Tracker::initialize(TimeFrame&, static declarations/configuration)
    -> build and validate locally
    -> one fallible atomic TimeFrame configuration commit
```

The workflow supplies declarations, conditions, and I/O; it does not build the
operational graph. `Loader` stages and atomically installs event-derived input
only. L4 may add passive TimeFrame storage/commit/reset primitives, but must
not add a TimeFrame method that constructs, selects, or validates topology.

`SurfacePlanBinding` remains the source-qualified partition and compact-slot
owner in L3. Moving that configuration and the physical workspace into
TimeFrame is intentionally deferred to L4/L5.

## Validation coverage

The focused L3 tests cover:

- non-contiguous global `SurfaceId` lookup through the graph-owned compact
  index, non-identity ordered traversal, sparse disconnected adjacency,
  holes, and seeding masks;
- reproduction of the owning graph through the device-safe POD view;
- fail-closed graph validation for invalid membership and topology;
- one combined graph shared by ITS/MFT bindings; and
- a source/dependency guard proving that deleted layout/topology ownership
  names are absent from common production include/src.

No Propagator, refit, covariance, topology-connectivity, publication,
workflow, or physics behavior is changed by L3.

## Exact residual scope

The following are intentionally not L3 work:

| Area | L3 disposition | Next boundary |
| --- | --- | --- |
| `TimeFrame` graph/configuration storage | remains passive and unchanged | L4 atomic configuration commit |
| `SurfacePlanBinding` ownership | remains application/interface-owned | L4 moves static partition authority |
| `SurfaceTrackingScratch` ownership | unchanged | L5 workspace ownership |
| raw ROFs and publication lifecycle | workflow-owned | out of L3 |
| detector-specific static surface declarations | adapter/application data | future adapter cleanup |

## Validation record

Package `daily-20260717-0700-local1` and the existing durable build were used.
The public graph-view layout change required rebuilding the downstream ITS,
MFT, and combined workflow targets before replay. A first attempt with stale
workflow executables failed closed with “cluster dictionary is not available”;
the retry below is the validated run after that rebuild.

```sh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
ninja -C /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
  O2lib-ITSMFTTracking O2lib-ITSCAWorkflow O2exe-its-ca-tracker-workflow \
  O2lib-MFTTracking O2lib-MFTWorkflow O2exe-mft-ca-tracker-workflow \
  O2lib-ITSMFTCombinedCAWorkflow O2exe-itsmft-combined-ca-tracker-workflow
```

The final serial suite was:

```sh
O2_BUILD_DIR=/Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
O2_PACKAGE=daily-20260717-0700-local1 \
.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh -- \
ctest --test-dir /Users/mpuccio/alice/run3/O2-worktree-builds/m6e3-commontrack-output-retirement \
  -L itsmft --output-on-failure -j1
```

Result: **100% tests passed, 0 tests failed out of 99**; every registered test
executed and there were no `Not Run` entries. Formatting checks used
`git diff --check c738b2a3c7 HEAD` and
`git clang-format --diff c738b2a3c7 HEAD`.

The fixture
`/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717`
passed `shasum -a 256 -c checksums.sha256` for all 43 entries before and after
replay. Replay artifacts are under:
`/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/m7f-l3-surface-graph-20260806/`.
The successful retry products are in `its-standalone-retry`,
`mft-standalone-retry`, and `combined-retry`, using run `303000`,
condition-not-after `1784207296000`, ITS diamond `(0,0,0)`, `pvRes=0.05`, and
one MFT thread.

| Leg | Tracks | Content hash |
| --- | ---: | --- |
| ITS standalone | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| ITS combined | 212 | `46913a67a7e2fe7462e29df0db264fa8` |
| MFT standalone | 68 | `8106b08571ca593c6b76ff72b761a680` |
| MFT combined | 68 | `8106b08571ca593c6b76ff72b761a680` |

The established `test/gate4-b42b-validation/compare_common_ca_output.C`
reported field-by-field equality for each L3 standalone product against the
accepted L2 parent product and for each standalone/combined pair. The MFT
comparison covered 2,992 forward-state/covariance values with maximum absolute
and relative delta zero. This covers initialized writer, CommonTrack,
cluster-reference, ROF, label, and sidecar content; only the known undefined
`MFTTrack.mInvQPtSeed` byte artifact is excluded.

No CUDA or HIP compiler was available in the pinned environment (`nvcc` and
`hipcc` were absent), so no device-build result is claimed. The named user
stash `stash@{0}` (`user WIP: TripletTrackingRAndD.md`) remained untouched.

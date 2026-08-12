# Header intentionality and surface-graph simplification audit

- Status: audit/design complete; no implementation authorized
- Date: 2026-08-07
- Audited revision: `4a3c628ee5c2a5173593720c4d1bf9bdbb674f3f`
- Branch: `codex/itsmft-header-graph-audit`
- Scope: every header under `include/ITSMFTTracking/`, including `detail/`
- Architecture prerequisites: [post-M7 audit](0009-post-m7-intentionality-cleanup-audit.md) and [cleanup plan](0010-post-m7-cleanup-implementation-plan.md)
- Execution status: [header and graph cleanup campaign](0012-header-and-graph-cleanup-execution.md)

This note is an audit and design only. It changes no C++, CMake, tests,
workflows, physics, configuration, defaults, comments, or public API. Every
implementation item below remains separately reviewable and must pass its
stated gate.

## 1. Fixed architectural boundary and non-goals

The audit preserves these responsibilities:

- `TimeFrame` is passive reusable storage. It owns configuration, surface
  graphs, workspace/capacity, normalized event data, and generic results. It
  neither initializes itself nor applies policy to itself.
- `Tracker` is the non-owning orchestrator. It validates supplied declarations,
  atomically installs configuration into a `TimeFrame`, and runs tracking.
- `TrackerTraits` remains the CPU/GPU kernel implementation seam. It is not a
  merge candidate for `Tracker`.
- `MultiSourceTimeFrameLoader` is non-owning and atomically loads normalized
  event input into a configured `TimeFrame`.
- Workflows/adapters own raw ROFs, timing construction and backing storage,
  publication validity, detector compatibility, typed output, and writers.
- A `SurfaceGraph` is detector-neutral and may contain multiple detector
  identities. ITS and MFT are disconnected today because their two submitted
  subgraphs contain no cross-detector edges.
- Cylinder and disk paths stay structurally parallel where physics permits.
- Device-facing value types and POD views remain free of host-only ownership,
  exception, STL-allocation, and workflow dependencies.

Explicit non-goals are changing propagation, material, covariance, cuts,
track/candidate order, hole behavior, seeding behavior, accepted results,
output bytes, detector defaults, mixed-family physics, or adding a utility
umbrella, manager, dispatcher, plan wrapper, or new detector-specific core.

## 2. Method and measured inventory

The inventory was produced from the filesystem, not from the approximate
count. Direct includes, declarations, host/device guards, line counts, symbol
uses, production construction sites, workflow consumers, and focused tests
were searched at the audited revision.

| Location | Headers | Lines |
|---|---:|---:|
| `include/ITSMFTTracking/*.h` | 61 | 9,516 |
| `include/ITSMFTTracking/detail/*.h` | 6 | 1,944 |
| **Total** | **67** | **11,460** |

The proposed end state is 36 root-level headers and 12 `detail/` headers, 48
headers total. Sixteen files disappear through cohesive merges, `SeedAnchor`
is evidence-gated for deletion, and the five temporary `TransitionPolicy`
headers are replaced by three stage-coherent detail headers. Eight
implementation-only files move to `detail/`. This is a target map, not
authorization to land all reductions together.

Classification keys used below are exactly the requested categories:

1. **Keep** as an independent concept (for an existing `detail/` header, keep
   as an independent detail concept rather than promote it);
2. **Detail**: move to `detail/`;
3. **Merge** into the specifically named header;
4. **Delete**/inline/replace;
5. **Evidence**: do not decide until the named evidence exists.

## 3. Complete header inventory and disposition

| Header | Class | Intentional responsibility and disposition |
|---|---:|---|
| `BarrelSurfaceStateOperations.h` | 1 | Cylinder state-operation leaf; keep separate from disk formulas and from orchestration. |
| `Cell.h` | 1 | Device-portable CA cell and whole-track seed values. Fixed capacity is intentional. |
| `ClockTimingPublicationView.h` | 3 | Merge into `CommonTrackOutputAdapter.h`; it exists only to translate output timestamps. |
| `ClusterDecoder.h` | 3 | Merge into `ClusterDecoding.h`; decoder interface and geometry adapter are one host decoding boundary. |
| `ClusterDecoding.h` | 1 | Destination for the complete host decoding boundary; keep separate from device measurement values. |
| `ClusterSource.h` | 3 | Merge into `IOUtils.h`; the input declaration has no lifecycle outside a loader call. |
| `CommonTrack.h` | 1 | Device-portable generic result and membership reference. |
| `CommonTrackOutputAdapter.h` | 1 | Host application output staging; remains outside the generic tracking core. |
| `CommonTrackShadow.h` | 2 | TrackerTraits-only transactional result-commit implementation. |
| `ConfigKeyValuesPreflight.h` | 1 | DPL-free configuration-string validation boundary. |
| `Configuration.h` | 1 | Algorithm iteration/vertexing parameters shared with frozen consumers. |
| `DecodedCluster.h` | 3 | Merge into `ClusterDecoding.h`; transient result of that one host pipeline. |
| `DetectorPublicationAdapter.h` | 2 | Detector compatibility completion used only by workflow composition and tests. |
| `DetectorTrackingOperationAdapterSupport.h` | 2 | Detector application refit/beam support, not a generic public concept. |
| `ForwardSurfaceStateOperations.h` | 1 | Disk state-operation leaf; keep structurally parallel to the cylinder leaf. |
| `IOUtils.h` | 5 | Mixed legacy conversion and live bounded decoding. Decide only after downstream API evidence; do not grow it. |
| `ITSSharedClusterCompatibility.h` | 2 | Workflow-owned ITS publication sidecar and transaction. |
| `ITSSurfaceSpec.h` | 3 | Merge into `ITSMFTDetectorDefinitions.h`; it is catalog data, not an independent facility. |
| `IndexTableConfiguration.h` | 1 | Host binding/validation kept apart from the device-portable LUT type; remove its temporary policy-tag template in the retirement campaign. |
| `IndexTableUtils.h` | 1 | Device-portable LUT geometry and search operations. |
| `LayerMask.h` | 3 | Merge into `SurfaceMask.h`; the 16-bit positional compatibility mask and 32-bit surface mask form one mask contract. |
| `MCLabelAccumulator.h` | 1 | Coherent host-side MC-label reduction algorithm; not merely an alias or flag. |
| `MFTFwdTrackHelpers.h` | 2 | MFT compatibility projection/refit helpers used by implementation code only. |
| `MFTPublicationCompatibility.h` | 2 | Workflow-owned MFT publication sidecar and transaction. |
| `MFTSurfaceSpec.h` | 3 | Merge into `ITSMFTDetectorDefinitions.h`; same reason as the ITS spec. |
| `MaterialPhysics.h` | 1 | Scalar material kernel and diagnostics; separate physics contract. |
| `MultiSourceFrame.h` | 1 | Owning normalized event frame plus device-facing read-only view. |
| `MultiSourceLoading.h` | 3 | Merge into `IOUtils.h`; errors/results/free staging belong to loader lifecycle. |
| `IOUtils.h` | 1 | Architectural non-owning atomic Loader; destination for its declarations and failures. |
| `NativeRefitDriver.h` | 1 | Family-neutral whole-seed refit sequencing over `Propagator`. |
| `ITSMFTDetectorDefinitions.h` | 1 | Shared detector-default data used by configuration and static catalogs; merging would widen a high-fan-in include. |
| `Propagator.h` | 1 | Descriptor-driven propagation and refit-leg operation boundary. |
| `ROFTimingUniformity.h` | 3 | Merge into `SurfaceTiming.h`; one validation/result used only to construct timing configuration. |
| `ROFViews.h` | 1 | Device-portable, non-owning runtime timing/mask views. |
| `RefitLegAssembly.h` | 3 | Merge into `NativeRefitDriver.h`; its only production consumer is that driver. |
| `SeedAnchor.h` | 4 | Delete with the unused anchor-taking overloads after downstream evidence; `Inner` is test/historical only. |
| `SurfaceKind.h` | 1 | Small but real device-shared representation tag; isolation prevents pulling the 259-line state into kernel-operation headers. It is not a replacement topology/dispatch policy. |
| `ITSMFTDetectorDefinitions.h` | 1 | Compile-time ITS/MFT catalog data and runtime projection; absorbs both detector spec files. |
| `StaticSurfaceDescriptor.h` | 3 | Merge into `SurfaceSpec.h`; the type exists solely to define/validate compile-time specs. |
| `SurfaceCatalogView.h` | 3 | Merge into `SurfaceDescriptor.h`; descriptor plus narrow borrowed catalog view is one device concept. |
| `SurfaceDescriptor.h` | 1 | Runtime surface identity, kind, reference coordinate, material, and catalog view. |
| `SurfaceGraph.h` | 1 | Owning immutable graph plus device POD view; representation may simplify but remains a concept. |
| `SurfaceGraphBuilder.h` | 1 | Host declaration-to-graph validation/build boundary; must not enter device headers. |
| `IdTypes.h` | 1 | Strong device-portable identifiers; absorbs `SurfaceMeasurementIndex`. |
| `SurfaceKinematicState.h` | 1 | Device-portable fitted state and typed views; absorbs its covariance-free companion. |
| `SurfaceKinematicStateLegacyAdapters.h` | 2 | Host-only legacy import/export used at adapter/test boundaries. |
| `SurfaceLinearizationReference.h` | 3 | Merge into `SurfaceKinematicState.h`; it is explicitly meaningful only beside that state. |
| `SurfaceMask.h` | 1 | Device-portable generic mask and positional conversion; absorbs `LayerMask`. |
| `SurfaceMeasurement.h` | 1 | Device-portable normalized measurement and its identity/geometry components. |
| `SurfaceMeasurementAdapters.h` | 3 | Merge into `ClusterDecoding.h`; the two projections consume only `DecodedCluster`. |
| `SurfaceMeasurementIndex.h` | 3 | Merge into `IdTypes.h`; duplicate strong-identifier mechanics with a 32-bit value. |
| `SurfaceSpec.h` | 1 | Compile-time spec concept/validation/concatenation; absorbs static descriptor types. |
| `SurfaceStateOperationResult.h` | 1 | Shared failure vocabulary avoids coupling barrel, disk, material, and propagator headers. |
| `SurfaceTiming.h` | 1 | Device-safe intervals/timestamps plus guarded host builders; absorbs uniformity validation. |
| `SurfaceTrackingScratch.h` | 2 | Private TimeFrame-owned workspace implementation, not independently configured public state. |
| `TimeFrame.h` | 1 | Passive reusable entity and sole generic state owner. |
| `TimeFrameLoadFailure.h` | 3 | Merge into `IOUtils.h`; workflow classification is part of loader failure contract. |
| `Tracker.h` | 1 | Non-owning initialization/execution orchestrator. |
| `TrackerTraits.h` | 1 | CPU/GPU kernel seam; must remain distinct from `Tracker`. |
| `TrackingConfigParam.h` | 1 | ROOT-visible workflow configuration registration/data. |
| `TrackingOperationAdapter.h` | 1 | Narrow call-scoped application refit operation. |
| `detail/SurfacePlanBinding.h` | 5 | Keep private for now, but obtain evidence for replacing public TimeFrame exposure with a narrower frame-owned partition view. |
| `detail/TransitionPolicy.h` | 4 | Delete. The tag is a one-to-one restatement of `SurfaceKind`; use `SurfaceKind` for operation selection and `SurfaceKind` only for state representation. |
| `detail/TransitionPolicyBinding.h` | 4 | Inline live host binding into `TrackerTraits.cxx`, relocate hit-attachment configuration to `detail/CellFinding.h`, and delete the unwired material-mode preflight after confirming no downstream caller. |
| `detail/TransitionPolicyDispatch.h` | 4 | Delete. Derive deterministic transition/cell/road/neighbour schedules directly in `SurfacePlanBinding`; no tag grouping or dispatcher remains. |
| `detail/TransitionPolicyOperations.h` | 4 | Replace with stage-coherent `detail/TrackletFinding.h` and `detail/CellFinding.h`; delete superseded policy refit helpers after downstream evidence. Do not split architecture by detector. |
| `detail/TransitionPolicyState.h` | 4 | Replace the useful fixed-layout values with `detail/TrackingKernelParameters.h`; delete `TransitionPolicyTraits` and policy-tag mapping. Preserve device ABI validation. |

## 4. Cohesive consolidation map

Paths in consumer lists are relative to `Detectors/ITSMFT/` unless stated
otherwise. Tests are exact direct header consumers at the audited revision.

### 4.1 `RefitLegAssembly.h` into `NativeRefitDriver.h` (independent mechanical slice)

- **Responsibility/lifecycle:** assembling the three traversal-ordered refit
  legs immediately before the native driver executes them.
- **Why accidental:** `assembleRefitLegSlots()` has one production include,
  `NativeRefitDriver.h`; its file comment still says “unwired” although all
  three driver legs call it.
- **API/dependencies:** retain the inline function and tests; replace includes,
  then delete the old header. No added dependency because the driver already
  includes `Cell.h` and measurement types transitively.
- **Host/device:** both declarations are under `!GPUCA_GPUCODE`; unchanged.
- **Exact consumers:** `common/tracking/include/ITSMFTTracking/NativeRefitDriver.h`
  and `common/tracking/test/testRefitLegAssembly.cxx`.
- **Smallest slice:** move the function byte-for-byte, update the test include,
  delete the file and its stale “unwired/Gate 3” prose.
- **Deletion/gate:** no include or file remains; run header self-containment,
  `testRefitLegAssembly`, `testDriveRefitLeg`, `testPropagator`, the affected
  common target, serial `itsmft`, formatting, and replay parity.

### 4.2 Identifier and mask consolidation

`SurfaceMeasurementIndex.h` belongs in `IdTypes.h`; `LayerMask.h` belongs in
`SurfaceMask.h`.

- **Responsibility:** strongly typed device identifiers; fixed-width generic
  and positional masks.
- **Why accidental:** the measurement index repeats the strong-wrapper shape;
  `SurfaceMask.h` already includes `LayerMask.h` solely for positional
  conversion. Neither tiny file owns an independent lifecycle.
- **API/dependencies:** names, widths, layout assertions, and semantics remain.
  `CommonTrack.h`, `MultiSourceFrame.h`, and `testSurfaceSpec.cxx` change the
  index include. Mask direct consumers are `Cell.h`, `Configuration.h`,
  `SurfaceMask.h`, `SurfaceTrackingScratch.cxx`, `TrackerTraits.cxx`,
  `testLayerMask.cxx`, and `testSurfaceMask.cxx`.
- **Host/device:** preserve every `GPUhd/GPUhdi` annotation and ABI assertion;
  no STL or host facility enters either destination.
- **Smallest slices:** land identifiers and masks as separate commits.
- **Deletion/gate:** old includes/files absent; compile device-facing header
  fixtures and run identifier, cell, layer-mask, surface-mask, common-track,
  and normalized-frame tests plus a real device build when available.

### 4.3 State representation consolidation

Merge `SurfaceLinearizationReference.h` into `SurfaceKinematicState.h`.

- **Responsibility:** fitted state and its covariance-free, paired
  linearization point.
- **Why accidental:** the reference explicitly has no meaning alone and
  already includes the full state header.
- **API/dependencies:** direct consumers are both family operation headers,
  `Propagator.h`, `detail/TransitionPolicyOperations.h`,
  `SurfaceLinearizationReference.cxx`, and tests
  `testBarrelLinearizationOperations`, `testCovarianceSanitization`,
  `testForwardLinearizationOperations`, and
  `testSurfaceLinearizationReference`. No consumer gains a new dependency.
- **Host/device:** retain the 32-byte POD assertions and the host-only factory
  guard; do not mix legacy adapters into the destination.
- **Slice/criterion/gate:** move declarations unchanged, update the `.cxx` and
  includes, delete the file; run state, linearization, covariance, propagator,
  refit, header, and device compilation gates.

### 4.4 Descriptor/spec/catalog consolidation

Merge `SurfaceCatalogView.h` into `SurfaceDescriptor.h`,
`StaticSurfaceDescriptor.h` into `SurfaceSpec.h`, and both
`ITSSurfaceSpec.h`/`MFTSurfaceSpec.h` into `ITSMFTDetectorDefinitions.h`.

- **Responsibility:** runtime descriptor plus borrowed catalog view;
  compile-time descriptor plus validation; concrete application catalog data.
- **Why accidental:** every destination already includes its fragment, and
  each removed header is a type/data fragment used to build only the
  destination concept.
- **API/dependencies:** catalog-view direct consumers are
  `compileCommonTrackerOnboarding.cxx`, `MFTFwdTrackHelpers.h`,
  `MultiSourceLoading.h`, `NativeRefitDriver.h`, `Propagator.h`,
  `SurfaceGraph.h`, `SurfaceGraphBuilder.h`, `SurfaceTrackingScratch.h`,
  `TrackingOperationAdapter.h`, `detail/TransitionPolicyOperations.h`,
  `testCombinedStaticSurfaceCatalogTopology.cxx`, and
  `testMFTNormalizedRefit.cxx`. Static descriptor direct consumers are the
  two spec headers, `ITSMFTDetectorDefinitions.h`, and `SurfaceSpec.h`. Each
  detector spec is included only by `ITSMFTDetectorDefinitions.h` and
  `testITSMFTSurfaceSpecProjection.cxx`.
- **Host/device:** runtime descriptor/view and static descriptor remain POD;
  concrete catalog arrays remain host compile-time application data. Do not
  pull detector catalog data into `SurfaceDescriptor.h`.
- **Smallest slices:** descriptor/view, static descriptor/spec, then concrete
  detector catalogs—three commits.
- **Deletion/gate:** self-contained destinations, no old paths, exact static
  layout/catalog projection tests, graph/loading tests, downstream header
  compile, and device compile for runtime PODs.

### 4.5 Timing publication and timing validation

Merge `ClockTimingPublicationView.h` into `CommonTrackOutputAdapter.h` and
`ROFTimingUniformity.h` into `SurfaceTiming.h`.

- **Responsibility:** output timestamp/ROF conversion; source timing
  construction and validation.
- **Why accidental:** the clock view is embedded only in output contexts;
  uniformity produces one `ROFTimingConfig` and declares itself loading detail.
- **API/dependencies:** clock direct consumers are ITS/MFT/combined tracker
  spec headers, `CommonTrackOutputAdapter.h`, `testL7AdapterOwnership.cxx`,
  and `testM6dMFTMigration.cxx`. Uniformity direct consumers are ITS, MFT,
  and combined `CATrackerSpec.cxx` plus `testROFTimingUniformity.cxx`.
- **Host/device:** the clock view stays host-only and must not enter
  `SurfaceTiming` or `ROFViews`; uniformity stays inside the existing host
  guard in `SurfaceTiming.h`.
- **Slices/criterion/gate:** separate commits; remove old paths and run output
  adapter/publication/timing tests and exact writer/replay parity.

### 4.6 Host decoding boundary

Make `ClusterDecoding.h` the cohesive host API by absorbing
`ClusterDecoder.h`, `DecodedCluster.h`, and
`SurfaceMeasurementAdapters.h`. Move `SurfaceMeasurementDecodeResult` out of
`IOUtils.h` in the same dependency-untangling slice; leave legacy/general
`IOUtils` functions pending evidence.

- **Responsibility:** bounded pattern parsing, typed decode failure, decoded
  facts, normalized projections, decoder interface, and geometry-backed
  implementation.
- **Why accidental:** the current chain is inverted: `ClusterDecoder.h`
  includes the 276-line `IOUtils.h`, while `IOUtils.h` includes only the tiny
  cursor/error header. Decoded facts and their two projections are never used
  outside decoding or decoder fixtures.
- **API/dependencies:** decoder direct consumers are ITS/MFT/combined workflow
  spec headers and tests `testCombinedTrackingComposition`, `testCommonTrack`,
  `testComputeLayerCellsOrchestration`,
  `testComputeLayerTrackletsOrchestration`, `testGate4B31OwnershipContract`,
  `testIndexTableActivation`, `testM6e2ITSWorkspaceMigration`,
  `testMultiSourceLoading`, `testMultiSourceTimeFrameLoader`,
  `testTimeFrameCovarianceLifecycle`, `testTimeFrameLifecycle`,
  `testTimeFrameNormalizedSource`, and `testTrackerFailureContract`.
  `DecodedCluster` is also directly consumed by `IOUtils.cxx`,
  `testSurfaceMeasurement`, and `testSurfaceMeasurementAdapters`; the latter
  projection header has the same list except the cell-orchestration test and
  adds no production consumer beyond `IOUtils.cxx`.
- **Host/device:** the merged header is explicitly host-only. It may produce
  `SurfaceMeasurement`, but no decoder, geometry, dictionary, exception, or
  STL owner enters a device view.
- **Smallest slice:** first move result/decode declarations to break the
  `ClusterDecoder` -> `IOUtils` dependency, then move types/projections and
  update includes; do not delete `convertCompactClusters` in this slice.
- **Deletion/gate:** three old headers absent, no include cycle, public-header
  self-containment, all listed decoder/loading/covariance tests, ITS/MFT
  geometry-backed loading, downstream build, and replay parity.

### 4.7 Atomic loader boundary

Merge `ClusterSource.h`, `MultiSourceLoading.h`, and
`TimeFrameLoadFailure.h` into `IOUtils.h`.

- **Responsibility:** one call-scoped source declaration, staging result/error,
  atomic load component, and workflow failure classification.
- **Why accidental:** all three fragments exist only before/during/after one
  Loader call; the 32-line class currently relies on accidental includes for
  `TimeFrame`, `LoadSourcesResult`, and `SurfaceCatalogView`.
- **API/dependencies:** source direct consumers are ITS/MFT/combined spec
  headers, combined spec source, `CombinedTrackingTestSupport.h`, and tests
  `testCombinedTrackingComposition`, `testCommonTrack`, `testM6dMFTMigration`,
  `testM6e2ITSWorkspaceMigration`, `testMultiSourceLoading`, and
  `testMultiSourceTimeFrameLoader`. Loading-result direct consumers are
  `SurfaceTrackingScratch.h/.cxx`, `TimeFrameLoadFailure.h/.cxx`, and tests
  `testCommonTrack`, both compute-orchestration tests,
  `testMultiSourceLoading`, `testTimeFrameLifecycle`,
  `testTimeFrameNormalizedSource`, and `testTrackerFailureContract`.
  Failure direct consumers are ITS/MFT spec header/source, combined spec
  source, `TimeFrameLoadFailure.cxx`,
  `testCombinedTrackingComposition`, and `testTimeFrameLoadFailure`.
- **Host/device:** whole boundary stays host-only; forward-declare `TimeFrame`
  and hide `SurfaceTrackingScratch` so workflows stop paying for its 28
  includes.
- **Smallest slice:** assemble declarations without semantic changes, move
  implementations only as naming requires, update direct includes, then
  delete three files.
- **Deletion/gate:** one Loader header, no transitive scratch dependency,
  atomic success/failure/retry/replacement, malformed-pattern, timing
  classification, source isolation, allocator failure, workflow drop/fatal
  tests, downstream build, and replay parity.

### 4.8 `SeedAnchor` deletion

- **Responsibility/lifecycle:** the unused anchor-selecting overload of the
  family-local seed builders.
- **Why accidental:** production calls only the anchor-less outer convention;
  `Inner` exists solely as a frozen-ITS comparison fixture and the header's
  35-line history explains why it survived.
- **API/dependencies:** direct includes are only the barrel and forward
  operation headers. Symbol consumers additionally are their `.cxx` files,
  `SurfaceStateOperationResult.h`,
  `testBarrelLinearizationOperations.cxx`, and
  `testForwardLinearizationOperations.cxx`; design notes retain historical
  evidence.
- **Host/device:** the enum is laid out for a hypothetical future device call,
  but no live device or production call exists. Do not delete on that
  inference alone—perform the downstream/API search.
- **Smallest slice:** prove no external caller, archive the two comparison
  assertions in validation Markdown if still useful, delete the overloads,
  `InvalidSeedAnchor`, includes, tests that exist only for the overload, and
  the header. Do not alter the outer seed formula.
- **Deletion/gate:** whole-project/downstream symbol search is empty; outer
  seed byte/physics fixtures, cell building, refit, device compilation, full
  suite, and replay parity pass.

### 4.9 Retire the `TransitionPolicy` header cluster (recommended first campaign)

The five headers are a coupled temporary migration layer, not five enduring
concepts. `TransitionPolicyTag` has exactly two enabled values and each is
derived one-to-one from `SurfaceKind`; `TransitionPolicyGrouping` groups by
that derived value; and `TrackerTraits` retains
`mActiveTag`, policy-specific parameter optionals, and `*ForPolicy<Tag>`
methods. This is the tag-keyed duplicated orchestration that the
[descriptor-driven operation-boundary classification](0001-descriptor-driven-operation-boundary.md#5-classification)
says M5 removes.

- **Responsibility/lifecycle:** retain three real responsibilities only:
  deterministic source-binding schedules built once from immutable topology;
  compact validated kernel parameters bound once per iteration; and narrow
  tracklet/cell numerical leaves selected from endpoint `SurfaceKind` at the
  `TrackerTraits` seam. A transition-policy identity has no lifecycle of its
  own.
- **Why accidental:** the stored tag already disappeared from
  `SurfaceTransition`; every remaining tag is reconstructed from a surface
  kind. The grouping, traits, two parameter optionals, wrappers, tests, and
  comments mutually justify one another but encode no third choice that is
  absent from `SurfaceKind`/`SurfaceKind`.
- **Public API/dependencies:** remove `expectedPolicy` and
  `IncompatibleExpectedPolicyKind` from `SurfacePlanBinding`; remove
  `TransitionPolicyTag`, `TransitionPolicyTraits`,
  `TransitionPolicyGrouping`, `dispatchTransitionPolicies`, `mActiveTag`,
  and all `*ForPolicy` names. `IndexTableConfiguration` becomes an explicit
  `SurfaceKind`/family operation rather than a policy-tag template. No
  workflow, graph declaration, result, or adapter gains a replacement policy
  type.
- **Cohesive destinations:** `SurfacePlanBinding.h` owns validation plus its
  filtered transition, cell, road-start, and neighbour schedule; new
  `detail/TrackingKernelParameters.h` owns only the compact device-facing
  barrel/forward parameter records and ABI checks; new
  `detail/TrackletFinding.h` owns projection/search-window and transition-cut
  preparation leaves; new `detail/CellFinding.h` owns cell seed, hit
  attachment, compatibility, and road-precut leaves. Host conversion from
  `TrackingParameters` is private `TrackerTraits.cxx` initialization code.
  Existing `BarrelSurfaceStateOperations.h`,
  `ForwardSurfaceStateOperations.h`, `Propagator.h`, and
  `NativeRefitDriver.h` remain the state/refit leaves; do not duplicate them.
- **Host/device boundary:** the compact replacement parameter records retain
  standard-layout/trivially-copyable, size, alignment, offset, finite-value,
  and `GPUhdi` contracts. Graph scheduling and `TrackingParameters` binding
  remain host-only. Do not move STL, exceptions, configuration vectors, or
  scheduling into either stage-operation header's device-visible portion.
- **Exact production consumers:** `TrackerTraits.h/.cxx`,
  `IndexTableConfiguration.h/.cxx`, `SurfacePlanBinding.h`,
  `TransitionPolicyOperations.cxx`, and the one stale comment in
  `ForwardSurfaceStateOperations.h`. The five headers also directly support
  focused policy-tag/state/dispatch/operation tests and policy-named uses in
  surface-layout, binding, index-table, refit-hit/leg, tracklet/cell
  orchestration, material-preflight, TimeFrame, and migration guards. Rename
  tests only when they still protect a replacement behavior; delete tests
  whose sole assertion is policy-tag containment, mapping, or dispatch count.
- **Smallest safe slices:** P1 removes `TransitionPolicyGrouping`,
  `dispatchTransitionPolicies`, and redundant `expectedPolicy`, moving the
  byte-identical deterministic schedules into `SurfacePlanBinding`; P2
  replaces the compact parameter records and host binders; P3 converts
  tracklet operations, then cell operations, from tag templates to explicit
  structurally parallel leaves; P4 deletes unwired/superseded refit and
  material-preflight policy helpers, then deletes the last policy header,
  source name, test name, guard exception, and comment.
- **Deletion criterion and validation gate:** repository and downstream
  searches find no `TransitionPolicy` identifier or filename. For every
  current ITS, MFT, and disconnected combined declaration, transition IDs,
  cell IDs, compact slots, scheduled cells, road starts, CSR, accepted
  candidates, result order, and output bytes remain identical. Gate binding
  error classification, index-table configuration, parameter rejection,
  tracklet/cell/neighbour/road orchestration, state/refit/material tests,
  header self-containment, real device compilation, complete common/ITS/MFT
  builds, and replay parity. Mixed cylinder/disk edges remain rejected until
  separately specified and physics-validated.

## 5. Public-to-detail relocation map

These moves reduce public signaling without merging unrelated facilities.

| Header(s) | Exact production consumers | Consequence and gate |
|---|---|---|
| `CommonTrackShadow.h` | `TrackerTraits.cxx` only | Move unchanged; update `testCommonTrack`. Gate transactional reserve/append/rollback and result ordering. |
| `SurfaceTrackingScratch.h` | internal common headers/sources (`MFTFwdTrackHelpers`, Loader, Tracker, Traits, operation adapter, `SurfaceTrackingScratch.cxx`, `TimeFrame.cxx`, `TrackerTraits.cxx`) | TimeFrame remains owner; workflows must not include it. Move after Loader untangling. Gate all workspace/allocator/load/tracker tests and device header compile. |
| `DetectorPublicationAdapter.h` | ITS/MFT/combined tracker spec headers | Private workflow-composition compatibility. Tests: `CombinedTrackingTestSupport`, `testDetectorAdapterCompatibility`. Gate publication ordering/sidecars/writers. |
| `DetectorTrackingOperationAdapterSupport.h` | ITS/MFT/combined tracker spec sources | Private application refit/beam helpers. Tests: `CombinedTrackingTestSupport`, `testBeamPositionOwnership`, `testDetectorAdapterCompatibility`, `testMFTNormalizedRefit`. |
| `ITSSharedClusterCompatibility.h`, `MFTPublicationCompatibility.h` | `CommonTrackOutputAdapter.h`, `DetectorPublicationAdapter.h`, combined spec header | Private workflow-owned sidecars. Exact tests: `CombinedTrackingTestSupport`, `testCombinedTrackingComposition`, `testCommonTrack`, `testL7AdapterOwnership`, plus `testTrackerFailureContract` for ITS. Gate failure invalidation and writer parity. |
| `MFTFwdTrackHelpers.h` | `DetectorTrackingOperationAdapterSupport.h`, `MFTFwdTrackHelpers.cxx`, `TrackerTraits.cxx`, `TransitionPolicyOperations.cxx` | Private MFT compatibility arithmetic. Tests: combined composition, compute-tracklets orchestration, transition-policy operations. No formula move or family-specific architecture. |
| `SurfaceKinematicStateLegacyAdapters.h` | `CommonTrackOutputAdapter.h`, `DetectorTrackingOperationAdapterSupport.h` | Host-only legacy edge. Exact tests: barrel/forward state operations, barrel linearization, detector adapter, family material, and state. Preserve conversion bytes. |

`detail/SurfacePlanBinding.h` needs separate evidence. `TimeFrame.h` currently
exposes `BindingSet` and `getBinding()` using a detail type, while Tracker and
Traits consume its source/ordered-surface/compact-ID mappings. A future slice
should test whether a nested, read-only TimeFrame partition view can expose
only those spans without copying, policy tags, or another public header. Do
not merge the 288-line host builder into `TimeFrame.h` or introduce a public
`TrackingPlan` merely to remove the include.

## 6. Dependency and compile-cost findings

Direct reverse include counts identify the main public fan-in nodes:
`Configuration.h` 36, `TimeFrame.h` 32, `SurfaceTrackingScratch.h` 26,
`SurfaceGraphBuilder.h` 26, `TrackingConfigParam.h` 26, and
`TrackerTraits.h` 15. The most expensive fan-out headers are
`SurfaceTrackingScratch.h` (28 includes, 381 lines), `TrackerTraits.h` (19,
451), and `CommonTrackOutputAdapter.h` (large typed-output dependencies).

Concrete hotspots and corrections:

- `Cell.h` includes all of `Configuration.h` only for constants/types. Do not
  fix this with a utility header; first measure whether a forward/narrow
  parameter include is possible after masks are consolidated.
- `ClusterDecoder.h` pulls `IOUtils.h` and therefore configuration, geometry,
  patterns, measurements, and vectors into every decoder mock. Section 4.6
  removes this inverted chain.
- `IOUtils.h` pulls the full scratch implementation to get
  types it does not own. Section 4.7 removes that chain.
- `Tracker.h` includes `TrackerTraits.h`, scratch, graph builder, binding, and
  TimeFrame. Forward declarations can narrow this only after inline/API uses
  are audited; do not use a new facade.
- `TimeFrame.h` publicly includes the detail binding and, through it, the
  temporary policy chain. Retire that chain first; resolve the remaining
  binding exposure evidence before attempting a PImpl or nested view.
- `detail/TransitionPolicyOperations.h` is 912 lines and mixes tracklet, cell,
  transition-preparation, and superseded refit responsibilities. Replace it
  by stage, not by detector or family; keep cylinder/disk leaf structure
  parallel and preserve the compact device-parameter boundary.

## 7. File-by-file comment audit

License headers are mandatory boilerplate and excluded from this classification.
Line ranges refer to the audited revision. `K`, `C`, `M`, and `D` mean Keep,
Condense, Move to Markdown, and Delete. Replacement text follows `=>`.

| File | Substantial blocks and action |
|---|---|
| `BarrelSurfaceStateOperations.h` | K 22-38 operation invariants; C 42-69 material/state-chi2 => “Transactional; diagnostics identify preflight, scalar-kernel, or projection failure.” M 75-261 legacy transcription/history; retain per-symbol anchor, pairing, failure, and transaction rules in at most 2-4 lines. |
| `Cell.h` | D 42-45 Stage-B history. K 68-72 q/pT invariant. M 139-166 M6 narrative. C 179-217 mask/bounds => “Seed indices are compact plan positions; `SurfaceMask` is authoritative.” K/C 278-291 device-copyability reason => retain the `TimeEstBC` caveat in 4 lines. |
| `ClockTimingPublicationView.h` | C 19-21 => “Host-only immutable output view; delegates ROF semantics to `LayerTiming`.” |
| `ClusterDecoder.h` | K/C 24-34 host-only/lifetime contract. D 48-52 obvious production-adapter narration. K 69-71 fail-before-geometry rule. |
| `ClusterDecoding.h` | K 21-23 typed-boundary distinction and 37-41 bounded-cursor invariant; condense each to two lines. |
| `ClusterSource.h` | K/C 27-33 ownership/ID contract; K 43-51 exact ROF partition failure contract; C 38-41 and 60-62 field narration. |
| `CommonTrack.h` | C 25-35 => “Reference index is local to `surface`; resolve through `MultiSourceFrame::getMeasurement`.” M 48-99 architecture/lifetime essay; retain range, traversal order, and same-frame lifetime in 6 lines. C 110-119 ABI rationale. K 123-127 range invariant. |
| `CommonTrackOutputAdapter.h` | K 7-9 host/DPL boundary, 58-59 source-local ownership, 168-170 deterministic ordering, 194-195 out-of-span publication rule, and 358-360 persisted-chi2 rule; condense 168-170 wording to remove obsolete scratch reference. |
| `CommonTrackShadow.h` | K 25-26, 74-78 transactional owner-thread/rollback contract, 102-103 reserve-before-append invariant. D 41-43 ITS-specific historical rationale. |
| `ConfigKeyValuesPreflight.h` | C 36-57 to a 5-line parse/fail-closed contract; move option history/diagnostics examples to Markdown. |
| `Configuration.h` | C 89-92 compatibility accessors. M 120-122 and 147-148 Gate histories; replace with “Retained for frozen adapter consumers; generic runtime does not execute this feature.” |
| `DecodedCluster.h` | K/C 18-26 axis/projection facts, folded beside fields in `ClusterDecoding.h`. |
| `DetectorPublicationAdapter.h` | C 26-28 ownership summary. M/D 150-154 previous typed-staging history; keep only the undefined-byte compatibility invariant if still live. |
| `DetectorTrackingOperationAdapterSupport.h` | K/C 40-43 => “Detector edge produces only generic `TrackingCandidate`; typed output remains outside core.” |
| `ForwardSurfaceStateOperations.h` | K 24-25 model values, 33-36 host/float rule, 49-55 measurement/update semantics, 66-100 physics invariants. M 106-267 Stage-B/legacy narrative; retain concise seed anchor, pairing, propagation-model, failure, and transaction contracts per symbol. |
| `IOUtils.h` | M 51-81 configuration ownership history; retain ITS-no-op/MFT-source rule in 3 lines. K/C 105-144 decode-once and result validity. D 156-157 obvious function name. C 206-209 bounded-pattern invariant. |
| `ITSSharedClusterCompatibility.h` | D/M 21-25 “temporary/pending” history. K/C 49-51, 96-101, 140-141 accepted-sequence, seal timing, and transaction rollback contracts. |
| `ITSSurfaceSpec.h` | D the 25-line banner/history before declarations; catalog values and names are self-describing. |
| `IndexTableConfiguration.h` | K 13-21 host/device boundary; C 48-63 bind/equality contracts; K 86-89 checked-overflow contract. |
| `IndexTableUtils.h` | K/C 48-67 axis/device-capacity contract, 143-154 configuration equality/capacity facts, 198-199 and 280-281 family coordinate rules. |
| `LayerMask.h` | No substantial maintenance comment beyond file boilerplate; add none during merge. |
| `MCLabelAccumulator.h` | K/C 27-29 call/ordinal contract and 40-41 repeated-label behavior. |
| `MFTFwdTrackHelpers.h` | C 36-37 half-disk indexing invariant. M 129-139 removed-function history. K/C 146-154 normalized-measurement/source identity contract. |
| `MFTPublicationCompatibility.h` | D/M 23-27 temporary/history wording. K/C 42-44 and 84-86 accepted-sequence and transaction contracts. |
| `MFTSurfaceSpec.h` | D the 25-line banner/history before declarations; catalog values and names are self-describing. |
| `MaterialPhysics.h` | K 17-19 host-only boundary, 24-43 semantic value contracts, 65-73 fail-closed diagnostics. C 88-91 ABI caveat. C/M 113-172 formula narrative: retain units, flags, failure/transaction semantics; move derivation/history to Markdown. |
| `MultiSourceFrame.h` | K/C 32-55 device view/ownership, 63-94 bounds/null rules, 141-180 cached-pointer move/swap lifetime, 196-207 lookup/labels, 221-241 staged construction/cache rebuild. M 118-135 architecture history; retain owner/lifetime rule only. |
| `MultiSourceLoading.h` | K/C 34-43 decoder/kind failure distinctions, 74-79 timing detail, 85-94 atomic failure contract. M/D 46-54 Gate/pre-slice history; describe current not-configured and reserved-value semantics only. |
| `IOUtils.h` | K 4-6 component ownership and 24-25 atomicity; combine into one class comment. |
| `NativeRefitDriver.h` | M 27-34 M5d history. K/C 38-53 covariance physics rationale, 77-80 q/pT formula, 93-106 leg ordering/transaction. D 153-154 section banner after names make order clear. |
| `ITSMFTDetectorDefinitions.h` | K/C 15-37 source/units/default provenance; remove milestone wording. |
| `Propagator.h` | M 26-46 milestone/ADR narrative. K/C 53-134 conversion, propagation, hole, direction, and transactional contracts; remove decorative section rulers. |
| `ROFTimingUniformity.h` | D 15-17 “detail/not API” once merged. K/C 33-42 exact uniformity fields and empty-input result. |
| `ROFViews.h` | K 34-36 borrowed backing lifetime and 366-368 assembled event-view ownership; ordinary print/accessor code needs no narration. |
| `RefitLegAssembly.h` | M 24-36 history. K/C 37-63 traversal order, hole encoding, mapping precondition => 6-line contract beside merged function. D “unwired”. |
| `SeedAnchor.h` | M all file-level M4/M5 history and 42-56 frozen convention; archive with deletion evidence. |
| `SurfaceKind.h` | M 17-26 migration history. K/C 27-32 and 42-45 => “Representation tag only; never topology or dispatch policy. Derive from `SurfaceKind`.” |
| `ITSMFTDetectorDefinitions.h` | K/C 63-68 combined global-ID concatenation invariant. |
| `StaticSurfaceDescriptor.h` | C 21-22 open detector identity. D 34-35 and 74-75 cross-file narration. K 52-55 validated-spec projection precondition. |
| `SurfaceCatalogView.h` | C 21-36 => “Borrowed descriptor catalog with optional global-ID-to-compact-index map; no topology or ownership.” K 40-41 sparse lookup semantics. |
| `SurfaceDescriptor.h` | K/C 25-29 nominal material ownership/zero validity and 42-44 device-shared geometry scope. |
| `SurfaceGraph.h` | K/C 55-57 borrowed immutable device view and 144-146 immutable owner. Remove “only” claims if a pair view is prototyped. |
| `SurfaceGraphBuilder.h` | No substantial comments; declarations need one future sentence separating active selection, holes, seeding, and adjacency. |
| `IdTypes.h` | K 53-54 dense source-ID/invalid-value rule. |
| `SurfaceKinematicState.h` | K 23-25 parameter conventions, 36-37 shallow family check, 161-162/194-195 null-view rule, 205-209 signed q/pT invariant. M 59-120 validation history/empirical narrative; retain mathematical invariant and explicit non-PSD limitation in 10 lines, move reproducer/history to Markdown. |
| `SurfaceKinematicStateLegacyAdapters.h` | K/C 24-27 host-only edge restriction and 113-116 float-to-legacy output rule; remove migration-stage wording. |
| `SurfaceLinearizationReference.h` | C/M 21-54 alternatives/history; retain paired-state lifetime, omitted particle hypothesis, parameter conventions, and POD status in 8 lines. K/C 75-84 factory failure/transaction contract. |
| `SurfaceMask.h` | K 91-96 positional mapping and precondition; condense to 3 lines. |
| `SurfaceMeasurement.h` | No substantial blocks; concise field names/types are sufficient. |
| `SurfaceMeasurementAdapters.h` | K/C 34-36 disk axis/covariance mapping; add an equivalent one-line cylinder axis rule when merged. |
| `SurfaceMeasurementIndex.h` | K one-line “position local to one surface” beside merged type. |
| `SurfaceSpec.h` | K/C 58-59 static lifetime, 93-94 per-detector dense index, 123-129 definition-vs-valid concept distinction. |
| `SurfaceStateOperationResult.h` | C/M 29-93: keep one short discriminator sentence for each non-obvious failure; move legacy call-chain and milestone history to Markdown. Delete `InvalidSeedAnchor` with its API. |
| `SurfaceTiming.h` | K/C 25-100 signed/half-open/overflow/intersection/config-sign invariants, 122-144 result self-consistency, 168-172 explicit origin, 196-199 widening, 221-223 cross-source intersection. These protect arithmetic; shorten without moving. |
| `SurfaceTrackingScratch.h` | K 90-103 allocator destruction order (remove duplicated sentence), 117-124 capacity precondition, 131-144 reset/reseat contract, 160-182 allocator-safe stage/swap, 232-234 borrowed ROF view. M file-level milestone/group narrative and 185-186, 318-340 group labels/history; private member groups may keep terse labels. |
| `TimeFrame.h` | M file-level Gate history; replace with current passive-owner contract. K/C 91-104 view invalidation, 108-127 atomic commit/reset, 150-162 result same-frame lifetime, 170-180 allocator destruction order. D/M 193-206 “unpopulated by B3.1” and old `loadNormalizedSource` history; describe current members only. |
| `TimeFrameLoadFailure.h` | K 15-16 host boundary, 29-34 recoverable-vs-structural rule, 47-76 typed failure payload, 84-89 test-backed exhaustiveness limitation. Condense repeated wording. |
| `Tracker.h` | D/M 38-55 stale standalone-interface/history. K/C 62-70 complete result/boundaries, 116-128 adapter lifetime and all-iteration failure/reset contract. |
| `TrackerTraits.h` | K 44-48 forward-declaration dependency rule. C/M 60-123 failure reasons: keep discriminating failure conditions, move milestone narratives. D/replace 145-166 stale ownership (“TimeFrame owns workspace/config; Traits borrows per call”). C 178-180, 203-228 cache lifetime, 236-324 compact-ID/schedule/operation-binding invariants after policy retirement. M 355-439 M4-M6/Gate histories; retain only source-of-truth, units, lifetime, and non-owning rules. |
| `TrackingConfigParam.h` | D/M 35-36 “header only for now”; verify current use. C 79-121 namespace/thread semantics. M 127-137 Gate-4 activation history; retain default/ownership semantics only. |
| `TrackingOperationAdapter.h` | K/C 25-28 generic candidate and 42-44 call-scoped refit-only seam. |
| `detail/SurfacePlanBinding.h` | M file-level 11-47 migration history. K/C 96-103 order mapping, 126-130 stored order, 155-178 boundary validation/filtering, 230-235 scheduled-order distinction. |
| `detail/TransitionPolicy.h` | M the file-level migration history to this document. D tag definition, enabled-tag test, compatibility mapping, inverse kind mapping, and their assertions when the file is deleted; `SurfaceKind` is the replacement source of truth. |
| `detail/TransitionPolicyBinding.h` | K/relocate the host/device boundary and material/config validity contracts beside the surviving cell/config initialization code. M Stage-B and activation narrative. D policy/tag/preflight prose and the unwired correction-mode preflight when its code is removed. Replacement wording: “Host-only validation binds one iteration's configuration before traversal state is committed.” |
| `detail/TransitionPolicyDispatch.h` | K/relocate seeding-endpoint, deterministic schedule, cycle, and stable-order contracts to `SurfacePlanBinding`. D policy groups, per-tag spans, dispatch-count rules, and historical decision IDs with the deleted abstraction. |
| `detail/TransitionPolicyOperations.h` | M milestone/history portions of 50-53, 128-181, 246-267, 309-393, 442-506, 573-610, 650-690, 775-854. K/relocate formula, units, ordering, host/device, failure, hole, and unchecked-precondition facts in 195-231, 277-298, 435-439, 619-640, 699-756, 788-799, 864-903 beside the replacement stage leaves, at most 3-8 lines per operation. D policy-boundary and tag-dispatch narration. |
| `detail/TransitionPolicyState.h` | K/relocate the device finite check, parameter units/validity, and ABI locks to `TrackingKernelParameters.h`. M compatibility/Stage-B history. D `TransitionPolicyTraits`, expected-kind/tag assertions, and policy terminology. |

### 7.1 Safe comment-cleanup slices

1. Remove stale milestone/history comments in headers already moving or
   merging; move unique rationale into this note or the existing ADRs.
2. Clean passive-entity/component ownership wording in `TimeFrame`, Loader,
   `Tracker`, `TrackerTraits`, scratch, and binding. This is text-only but must
   be reviewed together because cross-file lifetime statements must agree.
3. Condense value-type/device comments (`Cell`, identifiers/masks,
   descriptors, state/reference, measurements, graph/view). Re-run device
   header compilation because comments sit beside ABI assertions.
4. Retire policy comments with the policy code, relocating only real schedule,
   device, formula, unit, failure, and physics invariants. Then condense
   operation comments one structurally parallel barrel/forward pair at a
   time. Physics reviewers must confirm no invariant was lost.
5. Clean adapter/loading/output comments and verify failure/writer contracts.

## 8. SurfaceGraph simplification investigation

### 8.1 What is primary today

Each immutable `SurfaceGraph` persists:

- descriptors, ordered surfaces, kind masks, seeding mask, and a fixed
  global-ID-to-compact-index table;
- `SurfaceTransition {from, to, skippedSurfaces, flags}`;
- `SurfaceCellTopology {firstTransition, secondTransition, hitSurfaces}`; and
- cell-successor CSR (`cellsByFirstTransitionOffsets` and
  `cellsByFirstTransition`).

`SurfaceGraphBuilder` starts from ordered active subgraphs, `maxHoles`, a
hole-surface mask, and seeding mask. It creates every forward `(from,to)` pair
whose intervening ordered surfaces fit the hole budget/mask, then creates
every composable pair of transitions whose combined skipped count fits the
budget. `finalize()` builds CSR. Thus transitions, cells, CSR, and hole
expansion are all primary persistent topology today.

### 8.2 Current consumers and what they truly need

| Consumer | True requirement | Derivable from ordered pairs? | When |
|---|---|---|---|
| `Tracker::initialize` / `TimeFrame::commitConfiguration` | validated counts and immutable per-iteration storage | Yes | Once before atomic configuration commit. |
| `SurfacePlanBinding::build` | source-owned surfaces/order, filtered transitions/cells, global-to-compact slots, road-start and scheduled cell IDs | Yes | Once per iteration/source after graph derivation. |
| `SurfacePlanBinding` schedule derivation | cycle/rank order, source-filtered cells, scheduled cells, road starts | Yes | Once at initialization; no policy grouping is required. |
| `SurfaceTrackingScratch::adoptPlan/initialise` | compact transition/cell counts and endpoint lookup | Yes | Capacity at configuration; views at iteration bind. |
| `TrackerTraits::computeLayerTracklets` | deterministic transition sequence and endpoints | The pair list is sufficient | Cache binding-filtered pair indices once. |
| `TrackerTraits::computeLayerCells` | every composable two-transition path and its three surfaces | Yes | Join pairs on `first.to == second.from`; cache before events. |
| `TrackerTraits::findCellsNeighbours` | successors whose first transition equals the current cell's second transition | Yes, but CSR is hot-loop useful | Derive/cache CSR once; do not scan all cells in candidate loops. |
| `TrackerTraits::findRoads` | cells whose second-transition endpoint is seed-eligible, in deterministic order | Yes | Derive from cells plus the seeding declaration once. |
| Tests/fixtures | direct mutation/error injection of transitions/cells/views | Not a runtime need | Rewrite fixtures against declarations/derived caches; keep targeted invalid-view tests privately. |

No production consumer requires cells or neighbour CSR to be authored by the
caller or serialized as independent topology. They require fast immutable
execution views. That supports making a pair list primary and cells/CSR
derived caches, while retaining those caches for hot-loop performance.

### 8.3 Proposed representation boundary

The simplest representation worth prototyping is:

```text
SurfaceGraph declaration (persistent authority)
  descriptors + stable SurfaceIds
  ordered active surface indices/subgraphs
  pair list: (fromCompactSurfaceIndex, toCompactSurfaceIndex)
  separate seeding declaration
  separate hole-policy declaration

Tracker initialization (local, fallible)
  validate/canonicalize pairs
  derive transition execution records
  derive cells and deterministic cell order
  derive cell-successor CSR
  derive binding/road-start/schedule caches directly from surface kinds
  atomically install all immutable data in TimeFrame

device/runtime view
  compact pair/transition array + only the derived arrays hot loops need
```

This keeps `TimeFrame` passive: `Tracker` performs derivation locally and
commits only a complete configuration. “Derived” does not mean “recomputed per
event.” Cells, schedules, and CSR should be cached per immutable graph (or per
source binding where ownership filtering changes compact IDs).

### 8.4 Memory, initialization, and device implications

Today each transition costs 12 bytes and each cell 8 bytes, plus roughly
`4*(E+1) + 2*C` bytes of CSR/index data, excluding vector capacity. A compact
pair is 4 bytes with two `uint16_t` indices. `skippedSurfaces`, cell hit masks,
and both transition IDs in a cell are derivable for ordered-chain policies;
dropping them from persistent authority reduces memory. If all derived caches
remain stored for runtime, peak configured memory changes less than the file
format suggests; the principal win is one source of truth and cheaper host-to-
device declaration transfer. A second phase may pack cache records only after
measurement.

Initialization gains an O(E + C) derivation and CSR build. This is configure-
once work. The temporary `TransitionPolicyGrouping` already scans transitions,
builds a topological rank, groups IDs, and sorts scheduled cells. Its rank,
road-start, and stable-order mechanics survive as binding/graph derivation;
the tag buckets and dispatcher do not. Pair-based adjacency can reduce the
current repeated edge scans, but initialization speed is secondary to
clarity. Candidate hot loops must retain O(1) indexed endpoint and CSR
successor access; deriving cells or scanning pairs inside them is not
acceptable.

The device view remains standard-layout/trivially-copyable pointer/count
data. Host `std::vector`, exceptions, graph/schedule builders, and output
facilities must not enter it. Compact indices are preferable on device, while
the catalog retains stable `SurfaceId` identity for measurements/results.

### 8.5 Stable identity and deterministic traversal

Keep `SurfaceId` stable in declarations and normalized data. Map it once to a
compact surface index. Define pair/transition identity as canonical pair-list
position after validation. Preserve current deterministic behavior by:

1. retaining declaration subgraph order and ordered-surface rank;
2. emitting pairs in `(subgraph, fromRank, toRank)` order;
3. rejecting duplicate/self/invalid pairs;
4. emitting cells in `(firstTransitionId, secondTransitionId)` order;
5. building CSR stably from that cell order; and
6. sorting scheduled cells by target rank with `CellTopologyId` as the stable
   tie-break, matching the current grouping contract.

Global `TransitionId`/`CellTopologyId` are currently internal configuration
identities, but tests depend on their order. A prototype must compare every
current ID sequence, binding slot, road-start list, scheduled list, and CSR
entry before considering the mapping safely replaceable.

### 8.6 Disconnected and future configurations

The current combined graph contains one ITS cylinder subgraph and one MFT disk
subgraph. The builder iterates each subgraph independently, so no pair crosses
between them. A single pair list naturally represents the same disconnected
components without detector branching.

Multiple detector identities within a connected component also need no graph
special case: endpoint descriptors retain detector identity while topology
uses compact indices. Current execution, however, rejects mixed-kind
transitions (`SurfaceGraph::validate`) and mixed-kind subgraphs (builder
preflight), because only same-kind cylinder and disk CA operation leaves
exist. A future cylinder-to-disk pair is representable but cannot be enabled
until a mixed-family operation/cut/material contract and physics validation
exist. Do not reintroduce a policy tag to model that future operation, and do
not introduce separate graph architectures for cylinders and disks.

## 9. Hole policy is not adjacency

Four decisions must remain explicit and separate:

| Decision | Current input | Meaning |
|---|---|---|
| Active-surface selection | `NLayers` prefix / submitted ordered subgraph | Which surfaces participate in this iteration. |
| Permitted skipped surfaces | `MaxHoles` + `HoleLayerMask` | Which missing active surfaces and how many a candidate may tolerate. |
| Seeding/start eligibility | `StartLayerMask` -> `seedingSurfaces` | Which cell traversal endpoints may start roads. |
| Adjacency | generated transitions/cells | Which surface-to-surface CA work is structurally available. |

A hole is missing accepted measurement support on an otherwise active
surface. It is not an inactive surface and not merely an absent edge. Current
code both expands edges from hole policy and later checks the final
`TrackSeed` hit mask with `isAllowed(MaxHoles, HoleLayerMask)`; any new model
must demonstrate why both checks remain equivalent or which one changes.

### 9.1 Viable models

1. **Preferred structural model for prototype: base adjacency plus derived
   hole expansion.** Persist ordered active surfaces and immediate/base pairs;
   keep `{maxSkipped, skippableSurfaceMask}` separate. At Tracker
   initialization, derive execution pairs (including allowed skips), a
   deterministic skip witness/mask, cells, and CSR. This cleanly separates
   authority while preserving current hot-loop shape. It is straightforward
   for the current ordered disconnected chains; an arbitrary DAG needs an
   unambiguous path/witness rule.
2. **Explicit execution pairs plus separate hole validation.** The declaration
   supplies every allowed `(from,to)` work pair, including long pairs; hole
   policy separately validates/annotates the skipped ordered span. This is
   simple and deterministic but leaves hole expansion encoded in the authored
   pair set, so the conceptual separation is weaker.
3. **Dynamic missing-hit traversal.** Persist only immediate pairs and let road
   construction cross missing measurements at runtime under a hole-state
   budget. This minimizes topology caches but changes candidate construction,
   ordering, and likely performance. It is an algorithm/physics project, not
   a header cleanup, and is not the preferred first model.

No behavioral model is recommended for production yet. Required physics
validation includes every current iteration's exact transition/cell/CSR and
road-start sequence; zero/one/multiple allowed holes; holes at inner, middle,
outer, and seeding surfaces; minimum track length; shared-cluster and refit
acceptance; ITS/MFT standalone and combined candidate/result ordering and
hashes; writer parity; and future mixed-component fixtures. Model 1 is only
the simplest design prototype to compare against those gates.

## 10. Recommended sequence, risks, and stop conditions

1. Run Section 4.9 Slice P1: remove tag grouping/dispatch and redundant
   `expectedPolicy` while preserving every binding and schedule sequence.
2. Complete policy retirement through P2-P4 in independently replay-gated
   stage slices; no policy identifier should survive as a compatibility alias.
3. Land the mechanical `RefitLegAssembly` merge from Section 4.1, then
   consolidate identifiers/masks and descriptor/spec fragments in small,
   device-verified commits.
4. Consolidate timing and host decoding, then the Loader boundary. This
   unlocks moving scratch and compatibility implementation headers to detail.
5. Perform comment cleanup independently of representation changes, deleting
   policy narration with the retired machinery.
6. Build a non-production pair-list derivation fixture that consumes current
   graph declarations and compares every derived ID/cache byte-for-byte.
7. Only after exact structural parity, decide whether `SurfaceGraph` stores
   pair authority plus caches or keeps the current representation.
8. Treat any hole-model or mixed-kind execution difference as a separate
   physics-reviewed campaign.

Primary risks are downstream public includes, accidentally preserving the
policy tag under a family-dispatch alias, accidental host dependencies in
device headers, changing stable ID/order while deriving schedules/caches,
hiding a second topology in `SurfacePlanBinding`, widening
`CommonTrackOutputAdapter` into a utility umbrella, and removing comments that
encode actual physics or transactional guarantees. Stop if a slice needs a
behavior/default change, cannot state one owner and deletion criterion, or
lacks downstream/device evidence.

## 11. Documentation validation contract

For this audit, validate all relative Markdown links, unique headings, table
delimiter/row column counts, inventory counts, and `git diff --check`. No full
build or replay is appropriate because the commits contain Markdown only.

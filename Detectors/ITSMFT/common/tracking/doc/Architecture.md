# Shared ITS/MFT Cellular Automaton Tracking Architecture

Status: Draft RFC; current topology terminology follows immutable `DetectorLayout` and pass-local `TraversalTopology`
Owners: ITSMFT tracking maintainers  
Scope: CPU tracking first, with GPU-compatible data views retained  

Post-Gate-4 target: [ADR 0007](decisions/0007-generic-tracking-engine-boundary.md), plan in [GenericTrackingEngineMigration.md](GenericTrackingEngineMigration.md)

## 1. Purpose

The target is one Cellular Automaton (CA) tracking implementation that can operate on a TimeFrame containing cylindrical layers, forward disks, or a mixture of both. ITS and MFT remain detector-specific at their input and output boundaries, while cluster bookkeeping, topology traversal, CA artefact construction, and track selection are shared.

The migration must satisfy two goals simultaneously:

1. Remove the duplicated ITS-derived implementation currently present in `ITStracking` and `ITSMFTTracking`.
2. Make detector geometry an explicit property of surfaces and transitions rather than inferring it from the number of layers.

The first combined milestone is one TimeFrame containing ITS and MFT surfaces with two disconnected tracking subgraphs. Tracks crossing between cylindrical and disk surfaces are a later milestone because they require a common state representation or explicit state conversions.

## 2. Current constraints

The current common implementation is a functional MFT prototype, but several assumptions prevent it from being the final shared design:

- Detector identity is inferred from `NLayers` (7 means ITS, 10 means MFT).
- `LayerMask` has 16 bits, while ITS plus MFT requires at least 17 surfaces.
- Topology identifiers are 8 bits and topology storage enumerates all possible combinations.
- A TimeFrame has one detector ID, dictionary, label source, index-coordinate mode, and output track type.
- The main CA loops branch on the detector instead of dispatching on a surface transition.
- The loader assumes a single ROF stream and identical ROF counts across all surfaces.
- The common library depends on the production ITS tracking library, leaving two algorithm implementations.

These are architectural constraints rather than isolated implementation defects. Compatibility shims may preserve current APIs during migration, but new code must not extend these assumptions.

## 3. Architectural invariants

The following are mandatory design rules.

1. Detector identity must never be inferred from a layer or surface count.
2. The common CA core must not depend on ITS or MFT workflow code.
3. Geometry decoding must terminate at an input adapter; the core consumes normalized measurements.
4. Hits must be addressed by a global `SurfaceId` and a detector-qualified `ClusterRef`.
5. Tracking topology must be an explicit, sparse, directed pass-local structure derived from an immutable layout.
6. Indexing and timing are surface properties.
7. Projection, compatibility, attachment, and fitting behavior are transition or tracking-policy properties.
8. `TrackITS` and `TrackMFT` are adapter outputs, not core storage types.
9. Production ITS behavior must be characterized before replacing its implementation.
10. CPU layouts should retain explicit device-friendly views so GPU migration does not require another data-model rewrite.

## 4. Dependency model

The desired dependency direction is:

```text
                       ITSMFT tracking core
                      /                    \
              ITS reconstruction       MFT reconstruction
              and output adapters       and output adapters
                      \                    /
                       detector workflows
```

The core may depend on common ALICE data formats, math, propagation interfaces, allocators, and GPU utility definitions. It must not link `O2::ITStracking` or `O2::MFTTracking` for its fundamental data structures.

During migration, temporary adapters may call the existing ITS or MFT refitter. Such dependencies must sit behind an adapter/policy target and must not leak into the core public headers.

## 5. Terminology

- **Surface**: a logical tracking measurement surface. It may be cylindrical, disk-like, or another future geometry.
- **Layout**: the immutable ordered collection of surfaces, component boundaries, and static hole/seed policy.
- **Edge**: one pass-local allowed search step from one surface to another.
- **Cell path**: two connected pass-local edges used to construct cells.
- **Cluster source**: one detector-qualified input stream containing clusters, patterns, ROFs, dictionary, and optional labels.
- **Tracking policy**: operations required by CA construction and fitting for one state/geometry combination.
- **External track**: a detector data-format object such as `TrackITS` or `TrackMFT`.

Inside the common core, new APIs should use `surface` rather than `layer`. Existing detector-facing APIs may retain layer or disk terminology.

## 6. Core data model

### 6.1 Identifiers

Identifiers must be explicit and strongly typed where practical:

```cpp
using SurfaceId = uint16_t;
using EdgeId = uint16_t;
using CellPathId = uint16_t;
using ClusterSourceId = uint16_t;

struct ClusterRef {
  ClusterSourceId source;
  uint32_t index;
};
```

The exact widths are implementation choices, but the combined 17-surface layout and a sparse topology must be representable without special cases. Invalid sentinels must be defined centrally.

### 6.2 Surface descriptor

Each surface describes geometry, material, measurement, timing, and indexing independently:

```cpp
enum class SurfaceKind : uint8_t {
  Undefined,
  Cylinder,
  Disk
};

struct SurfaceDescriptor {
  SurfaceId id;
  o2::detectors::DetID::ID detector;
  SurfaceKind kind;
  uint16_t detectorSurfaceIndex;
  float referenceCoordinate; // nominal radius for cylinders, nominal z for disks
  float radialMin;
  float radialMax;
  MaterialDescriptor material;
  ResolutionDescriptor resolution;
  TimingDescriptor timing;
  IndexingDescriptor indexing;
};
```

`SurfaceKind` is defined once in `IdTypes.h`. `Undefined` is permitted only
for a default-constructed state or reference; validated surface graphs contain
only `Cylinder` or `Disk` descriptors.

Geometry constants should be populated by ITS and MFT layout builders. Tracking parameters may override cuts and uncertainties, but should not duplicate the detector geometry definition.

### 6.3 Detector layout

`DetectorLayout` owns the dense ordered tracking layers, component boundaries,
and static hole-layer policy. `LayerId` is exactly a layer's position in this
container. The layout is immutable while a `TimeFrame` is processed; it owns
no pass-local edges, paths, adjacency, or schedules.

```cpp
class DetectorLayout {
 public:
  gsl::span<const SurfaceDescriptor> getLayers() const;
  gsl::span<const uint16_t> getComponentOffsets() const;
  LayerMask getHoleLayers() const;
};
```

The layout must support:

- ITS-only cylindrical layouts.
- MFT-only disk layouts.
- A combined layout with disconnected ITS and MFT subgraphs.
- Future layouts whose pass-local topology contains cylinder-disk edges.

`Tracker` combines the immutable layout with one iteration's configuration to
derive its `TraversalTopology`: active layers, edges, `CellPath`s, adjacency,
and road schedules. Algorithms use that topology for reachability rather than
inferring it from adjacent numeric layer IDs.

### 6.4 Surface mask

Replace the 16-bit layer mask with a surface mask capable of representing the supported maximum surface count. A 32-bit POD is sufficient for the immediate combined layout; a 64-bit mask provides more extension space at modest cost.

The selected width and binary layout must be recorded because it affects GPU views and serialized/debug representations. Mask operations must distinguish:

- Number of hits: population count.
- Traversed span or path length.
- Permitted holes along a topology path.

These quantities must not be inferred interchangeably.

### 6.5 Normalized cluster measurement

The core cluster representation is a 72-byte, standard-layout, trivially-copyable
`SurfaceMeasurement`. It contains only detector-independent information:

- Global position.
- Measurement coordinates and covariance in the surface frame.
- Sensor identity qualified by detector/source.
- Global `SurfaceId`.
- External `ClusterRef`.
- Cluster size or shape metadata needed by output adapters.
- Timing/ROF association.

ITS `TrackingFrameInfo` may initially be adapted into this representation, but it must not remain the semantic definition for all future surface types.

External identity is the pair `{ClusterSourceId, external index}`. A source identifies
an input stream, not a detector, so two streams from the same detector and equal
external indices remain distinct. Source IDs are dense and TimeFrame-local, remain
stable for the TimeFrame lifetime, and reserve the all-ones value as invalid. Sensor
identity independently contains both the detector ID and detector-local sensor ID.

Surface-frame coordinates have an explicit convention selected by the surface
descriptor:

- Cylinder: normal coordinate `q=xTF`, measured coordinates `u=yTF`, `v=zTF`, and
  `frameAngle=alpha`.
- Disk: normal coordinate `q=z`, measured coordinates `u=x`, `v=y`, and
  `frameAngle=0`.

The 2D covariance always corresponds to `(u,v)`. In particular, the MFT adapter must
decode this information from detector geometry and cluster errors. It must not infer
disk-frame semantics by copying the current synthetic MFT `TrackingFrameInfo`, whose
stored position axes and covariance axes do not form this contract. MC labels are not
embedded in each measurement; host-side label lookup is keyed by `ClusterRef`.

## 7. TimeFrame responsibilities

The common TimeFrame owns per-surface input and CA artefacts:

- Unsorted and sorted normalized clusters.
- External cluster references.
- Per-surface index tables.
- Used-cluster state.
- Tracklets, cells, neighbours, roads, and internal tracks.
- Per-surface timing and overlap views.
- The active layout and tracking iteration configuration.

It must not intrinsically own:

- A single detector ID.
- A single detector dictionary or geometry singleton.
- ITS-only vertexing artefacts that are unused by tracking.
- Detector output track containers.

ITS vertexing state is composed with the tracking TimeFrame. It owns primary vertices,
vertex labels and lookup tables, and vertexing-only scratch, while reading a common
measurement/timing view. Tracking receives an optional non-owning
`VertexConstraintView`; an empty view supports MFT and vertex-free operation.
Compatibility accessors and the inherited GPU TimeFrame may remain temporary facades
while CPU and device consumers migrate to POD views.

Framework allocation state owns memory resources, device mirrors, streams, and pinning
registrations, but not the semantic host data. Lifetimes must guarantee that common
tracking and vertexing owners outlive all host/device views and asynchronous work that
uses them.

### 7.1 Loading multiple sources

Loading accepts one or more detector-qualified sources:

```cpp
struct ClusterSourceView {
  ClusterSourceId id;
  o2::detectors::DetID::ID detector;
  gsl::span<const CompClusterExt> clusters;
  gsl::span<const unsigned char> patterns;
  gsl::span<const ROFRecord> rofs;
  const TopologyDictionary* dictionary;
  const MCTruthContainer<MCCompLabel>* labels;
  const ClusterDecoder& decoder;
};
```

The decoder maps a detector-local sensor/layer to a global surface and normalized
measurement. Decoder polymorphism is confined to this host loading boundary. The
source retains ownership of detector-specific inputs; the common TimeFrame retains no
compact clusters, dictionaries, geometry singletons, or detector output types. Pattern
cursors and ROF ordinals are maintained independently per source.

### 7.2 Timing

Per-source ROF timing is validated at the loading boundary as checked, signed 64-bit
TimeFrame-relative bunch-crossing intervals. The intervals are not stored in
`TimeFrame`; raw ROFs and timing tables remain workflow-owned. Loading must not
require equal ROF counts, lengths, or delays across detectors. A measurement's
`sourceROF` is interpreted together with the source in its `ClusterRef`.

Cluster and track timestamps must remain convertible to detector-specific output ROF records. Tests must cover triggered and continuous readout configurations.

## 8. Sparse tracking topology

The topology contains only configured transitions and valid connected transition pairs. It must provide device-friendly views for:

- Transition lookup.
- Successor transitions or cells.
- Cells beginning with a transition.
- Path/hole validation.

The topology builder validates:

- Surface IDs exist.
- Transitions are directed and non-duplicated.
- Connected cell transitions share the pivot surface.
- Unsupported surface-kind/state combinations are rejected before tracking.
- Configured minimum track lengths are reachable in the graph.

For the first combined milestone, the ITS and MFT transition sets remain disconnected.

## 9. Indexing

Index tables are per surface, not per TimeFrame. An indexing descriptor defines:

- Coordinate mapping, such as phi-z or x-y.
- Row and column ranges.
- Bin counts.
- Periodicity and boundary handling.

Projection into the target index space belongs to the transition policy. The generic table stores and retrieves bins without knowing whether coordinates represent phi-z, x-y, or a future mapping.

## 10. CA algorithm and policy boundary

The shared `Tracker` and `TrackerTraits` own orchestration only:

1. Iterate enabled transitions and time overlaps.
2. Build tracklets from policy-provided projections and windows.
3. Combine compatible tracklets into cells.
4. Build cell-neighbour relations.
5. Traverse roads through the explicit topology.
6. Invoke policy fitting/attachment operations.
7. Rank and select internal tracks.

The policy boundary supplies operations equivalent to:

```cpp
SearchWindow projectSearchWindow(...);
bool buildCellSeed(..., SeedState& out, float& chi2);
bool cellsAreCompatible(...);
bool attachHit(..., SeedState& inOut, float& chi2);
bool finalRefit(..., InternalTrack& out);
```

Initial policy implementations are:

- Cylinder-cylinder using the ITS barrel state and propagation behavior.
- Disk-disk using the MFT forward state and fitting behavior.

The core loops must not branch on `DetID`. Dispatch may be compile-time or through compact policy tags in a transition descriptor, provided device execution remains practical.

### 10.1 Policy dispatch

Each configured transition carries a validated `TransitionPolicyTag`. The initial tags
are cylinder-cylinder and disk-disk; each maps to exactly one Stage-A state family.
The state family is derived from the tag rather than stored independently. A cell
inherits its policy from its two transitions, which must have equal tags, so the cell
does not store a second copy of the tag.

Host owners pre-group transition, cell-topology, road, and seed work by policy family.
CPU orchestration selects a template-specialized policy at the outer family or topology
loop. It must not dispatch inside candidate, cluster, neighbour, road, or refit loops.
GPU execution uses a separate specialized kernel launch per active `(stage, family)`;
single-family layouts therefore keep one launch per stage, while a disconnected
cylinder-plus-disk layout may issue two. Device kernels contain no policy branch.

Dispatch is an enum/tag switch to statically compiled policy functions. Device virtual
calls, function pointers, `std::variant` visitation, and generated dispatch tables are
not part of this contract. Device views use explicit pointer/count PODs, not host
`gsl::span` objects.

The tag is the single source of policy identity. Parameter storage and typed
family-state owner/view layouts are separate contracts: they must not duplicate tag or
family fields, and an untyped flat-float parameter block is not accepted without a
field-level ABI and bounds validation. Layout construction rejects invalid tags,
surface-kind/tag mismatches, mixed-policy cells, and all cylinder-disk transitions
until D008 is resolved. Configuration-dependent reachability and parameter-range
checks occur when configuration is bound to a valid layout, not in the topology owner
alone.

Post Gate 4, `TransitionPolicyTag` is classified as a temporary legacy
implementation detail to be contained behind private code — see
[ADR 0007](decisions/0007-generic-tracking-engine-boundary.md); the dispatch
mechanics above remain valid while it exists.

## 11. Track-state strategy

The migration has two track-state stages.

### Stage A: shared traversal, separate state families

Cylinder-cylinder paths use the barrel state and disk-disk paths use the forward state. The layout validator rejects transitions requiring a state-family change. This is sufficient for:

- One common codebase.
- One combined TimeFrame.
- Independent ITS and MFT tracking in disconnected subgraphs.

### Stage B: mixed-surface tracks

Before enabling cylinder-disk transitions, adopt one of:

1. A common global/curvilinear state usable on both surface kinds.
2. Explicit, tested barrel-to-forward and forward-to-barrel conversions at transition boundaries.

This choice requires a dedicated RFC covering covariance transport, material correction, magnetic-field behavior, numerical stability, GPU impact, and the final reference surface. It must not be hidden inside a detector trait.

## 12. Internal and external tracks

The TimeFrame stores a generic internal result:

```cpp
struct InternalTrack {
  TrackState state;
  float chi2;
  TimeEstBC timestamp;
  LayerMask hitMask;
  bounded_vector<SurfaceHitRef> hits;
  TrackFlags flags;
};
```

Detector output adapters perform:

- Internal barrel track to `TrackITSExt`.
- Internal forward track to `TrackMFT`, cluster-index vector, and seed pattern.
- MC-label creation or transfer.
- Detector-specific cluster ordering and ROF grouping.

The native forward refit now returns the generic `SurfaceKinematicState`/
`TrackingCandidate` result directly; no typed common-MFT refit/export record is
retained. MFT publication compatibility remains at its application edge.

## 13. Configuration

Configuration is divided into three layers:

1. **Layout configuration**: surfaces, geometry, timing, material defaults, indexing.
2. **Algorithm configuration**: cuts, holes, iteration steps, memory and threading.
3. **Detector presets**: ITS and MFT values for sync, async, cosmics, and future modes.

Per-surface vectors must be validated against the layout. Road and level requirements must derive from configured topology paths rather than `NLayers - N` formulas.

Existing configurable parameter names should remain supported through adapters during migration. Renaming user-facing keys is out of scope unless separately approved.

## 14. GPU considerations

CPU migration precedes GPU migration, but every common owner/view pair must define:

- Standard-layout and trivially-copyable guarantees where needed.
- Explicit identifier and mask widths.
- Host and device views without owning STL containers.
- Stable topology and surface descriptor arrays.
- Policy dispatch suitable for device compilation.

No CPU abstraction requiring virtual calls in inner tracking loops should be introduced without a device-side representation. Host-side adapter polymorphism is acceptable during loading and output conversion.

## 15. Migration sequence and gates

### Gate 0: baseline

- ITS and MFT characterization tests exist.
- Current performance and memory measurements are recorded.
- Known physics differences are documented.

### Gate 1: foundations

- Surface/layout model is available.
- Masks and topology support at least 17 surfaces.
- ITS-only and MFT-only layouts reproduce current topology.
- Common primitives no longer require duplicated definitions.
- Single-detector normalized loading, timing, catalog ownership, and legacy
  compatibility backfill are transactional and validated for ITS and MFT.

### Gate 2: common CA traversal

- Cylinder-cylinder and disk-disk policy/state boundaries use one host
  orchestration and one-shot outer dispatch.
- The first production CA operation is driven by validated sparse topology
  without a detector branch in its hot loop.
- Legacy graph/index parity, failure behavior, physics output, and performance
  are characterized and accepted before subsequent operations migrate.

### Gate 3: production migration

- Explicit opt-in ITS and MFT CPU workflows exercise the complete common core
  while the legacy workflows and production defaults remain unchanged.
- Remaining CA operations migrate to the common topology/policy boundary.
- Output, failure, physics, determinism, and resource contracts are validated
  through independent legacy/common jobs with identical inputs and conditions.
- Both implementations remain buildable on the same development branch for an
  extended A/B validation period. Default switching and legacy retirement
  require a separate decision based on that evidence.
- New adapters and workflow targets follow the dependency direction in Section
  4. Existing compatibility dependencies are documented and time-boxed; they
  are not replaced by new duplicated primitives merely to avoid touching the
  frozen legacy implementation.

### Gate 4: combined disconnected tracking

- One TimeFrame contains all ITS and MFT surfaces.
- Both disconnected topology components run in one invocation.
- Outputs are split correctly by detector and time.

### Gate 5: mixed-surface tracking

- A separate track-state RFC is approved.
- Cross-kind transitions are implemented and validated.
- Physics and performance acceptance criteria are met.

## 16. Compatibility and deprecation

Compatibility should be provided through aliases and thin adapters, not copied implementations. Every compatibility layer must identify:

- The old API it preserves.
- The new API it forwards to.
- The gate at which it can be removed.

Removal of public headers, configuration keys, workflow outputs, or ROOT dictionary types requires an explicit compatibility review.

## 17. Testing requirements

Each architectural component requires focused tests:

- Mask capacity, hit count, path span, and holes.
- Sparse topology construction and invalid layouts.
- Phi-z and x-y indexing boundaries and periodicity.
- Multiple source loading, sensor qualification, labels, and cluster sizes.
- Different per-surface ROF timing and overlap queries.
- Cylinder and disk policy unit tests.
- Determinism across serial and TBB execution.
- ITS and MFT end-to-end regression.
- Combined disconnected layout.
- Host/device view layout tests before GPU migration.

Physics comparisons must define tolerances rather than requiring bitwise equality where refitting order or floating-point evaluation legitimately changes.

## 18. Out of scope for the initial migration

- Changing detector data-format outputs.
- Retuning tracking cuts without a separately documented reason.
- Replacing the legacy MFT fitter before the disk adapter is stable.
- Enabling cylinder-disk tracks before the track-state RFC.
- Rewriting ITS vertexing unless required to separate it cleanly from tracking state.

## 19. Open decisions

The following must be resolved and recorded in the decision log:

1. Accepted: a 32-bit `LayerMask` supports at most 32 global surfaces in the initial shared layout.
2. Compile-time maximum layout versus fully runtime-sized host storage.
3. Exact normalized measurement representation.
4. Ownership boundary between tracking TimeFrame and ITS vertexing state.
5. Policy dispatch representation for GPU execution.
6. Internal track state container during Stage A.
7. Long-term mixed-surface track-state representation.
8. Migration and removal schedule for existing public ITS tracking types.

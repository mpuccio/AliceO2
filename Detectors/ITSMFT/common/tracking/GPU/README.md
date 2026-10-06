# GPU common-CA tracking

`gpu::TrackerTraitsGPU` derives from `TrackerTraits`, following the legacy ITS
GPU traits pattern. It overrides the tracklet, cell, neighbour and road stages
(the road stage includes seed selection and the final refit). The ITS CA
workflow selects the traits with `--tracklet-backend cpu|cuda|hip` (default:
`cpu`). Requesting an unavailable backend fails explicitly.

CUDA sources are built when O2 is configured with `ENABLE_CUDA=ON`. With
`ENABLE_HIP=ON`, `o2_add_hipified_library` translates the same `.cu` sources,
following the legacy ITS tracker build strategy. There is no separately
maintained HIP implementation. Both libraries can be loaded in one process:
their classes live in the inline namespaces `gpu::cuda` and `gpu::hip`
(selected by `ITSMFT_TRACKING_HIP`, which the HIP library exports).

## Structure

The layout follows the legacy ITS GPU tracker, in one header and one source:

- `TimeFrameGPU` (`ITSMFTTrackingGPU/TimeFrameGPU.h`) owns everything kept on
  the device: the timeframe data (sorted clusters, ROF cluster boundaries,
  index tables built on the device from the cluster bins, surface
  measurements, used-cluster flags, ROF masks, vertices), exposed to kernels as
  the same `FrameView` the CPU stages read, and every stage's results
  (`LinkedResults` of tracklets, cells and neighbours, road results, track
  seeds, refitted tracks). The timeframe is uploaded once, on one stream per
  layer, before its first tracklet stage; ROF masks are loaded per iteration
  and used flags again after every batch of accepted tracks
  (`usedClustersChanged`), like the legacy `load*Device` calls.
- The handlers (same header, `cuda/TrackingKernels.cu`) run one stage step on
  a `TimeFrameGPU`, as the legacy `TrackingKernels` handlers. Each leases one of
  a few streams with its scratch buffers, handed out in order and reused.
- Device memory comes from `cudaMalloc` by default: buffers only grow and are
  kept across timeframes. Inside the GPU reconstruction framework
  (`TrackerTraits::setFrameworkAllocator` with the `ExternalAllocator` that
  `GPUChainITS::GetITSMFTFrameworkAllocator` provides) it comes from the
  framework's pool instead, with the legacy tracker's scopes: the timeframe
  data are plain allocations, which the framework clears between timeframes,
  and stage results and scratch are allocated on the memory stack that the
  traits push before the tracklet stage of every pass and pop after its road
  stage (tags `ITSITER<iteration>`). The pool frees nothing individually, so a
  buffer that grows leaves its previous memory behind until that release.
- `TrackerTraitsGPU` orchestrates the stages. Host scheduling, shared-cluster
  arbitration (`acceptTracks`) and output publication stay on the CPU, so
  cluster-use decisions stay synchronized between passes.

Every CA step has one host/device implementation in
`ITSMFTTracking/TrackingKernels.h`, written for one source element:
`forEachTracklet` (a source cluster), `forEachCell` (a first-edge tracklet),
`forEachNeighbour` (a source cell), `roadRange`/`makeSeedInput`/`makeRoadTarget`
(a road seed), `forEachRoad` (a road job) and `prepareRefitJob` (a track seed).
The CPU stages run them in `tbb::parallel_for` into slab sinks; the handlers run
them in one generic grid-stride kernel into atomic appenders. The stages differ
only in that loop and in how the results are ordered. `TripletFitting.h`,
`Propagator.h`, `MaterialPhysics.h` and `RefitDriver.h` hold the underlying
physics. `GPUhdi()`/`GPUdi()` compile them as inline C++ on the CPU and as
device functions under CUDA/HIP.

## Stages

Nothing is downloaded between stages except counts; the host reads tracklets,
cells or the neighbour graph only to validate them or to build MC artefact
labels. Output buffers are sized by the capacity estimator and grown on
overflow (`runOnSlab`), as in the legacy tracker.

- Tracklets: one thread per source cluster walks the target index tables and
  appends tracklets per edge; the device sorts them by (first, second) cluster,
  removes duplicates and builds the lookup table. Edges run concurrently on
  independent streams.
- Cells: one thread per first-edge tracklet appends the cells (triplet factors
  and Jacobians included), which are sorted by (first, second) tracklet, the
  host's order, with their lookup table. Paths run concurrently.
- Neighbours: targets are visited in order of their outer layer. Each target's
  neighbours among its source paths are appended while the target levels are
  raised (`atomicMax`), then sorted by (target cell, source path, source cell)
  with their lookup table. Targets sharing an outer layer run concurrently.
  Cell measurements are resolved from the resident clusters.
- Roads: for every start path, the jobs and targets of each extension are built
  on the device from the cells and neighbour graph, exactly as the host builds
  them; extensions chain in place through alternating result buffers, and the
  seeds passing the seed selector are appended in order. Each start level's
  seeds are refitted on the device (measurement slots built from the resident
  clusters), the accepted tracks sorted (more clusters first, then lower chi2,
  ties in seed order, as the stable host sort) and only they are downloaded, in
  `TrackingCandidate` layout, for acceptance.

Count/scan/fill kernels remain where the output order must match the host's
emission order directly (road extensions, seed selection); elsewhere sorting
restores the host order after atomic appends. GPU allocation and execution
failures throw rather than silently dropping candidates.

## Validation

With validation enabled, the traits run each CPU stage too and compare. They
always pass the device results on to later stages; the CPU results never
replace them.

`--validate-tracklets` runs the CPU stage and checks tracklet identity/order,
timestamps and lookup tables exactly and `tanLambda`/`phi` within
`2e-6 * max(1, abs(cpu), abs(gpu))`. This allows single-precision math-library
rounding differences.
`--validate-cells` runs the CPU stage and checks cell identity, order, metadata
and lookup tables exactly, plus every factor/Jacobian entry within
`2e-5 * max(1, abs(cpu), abs(gpu))`. `--validate-neighbours` checks the cell
levels and the neighbour graph (identities, order, lookup tables) exactly.
`--validate-roads` checks extension identities, ordering, clusters and timestamps
exactly. Seed parameters/covariances use an absolute tolerance of `1e-4` plus a
relative tolerance of `2e-3`; chi2 uses `1e-2` plus `5e-3` relative. Non-finite
values fail validation. GPU seeds are passed to subsequent extensions and final
GPU refitting without replacement by CPU validation results.
`--validate-refit` compares the CPU and device refits of every seed, checks
acceptance decisions and metadata exactly, and uses the final-track numerical
tolerances below for both states and chi2.
Use the validation options with `--tracklet-backend cuda` or `hip`.

All shared physics calculations, including triplet automatic derivatives, use
single precision. The small-curvature geometry computes `1-index` directly,
using a series for `1-angle*cot(angle)`, instead of subtracting nearly equal
floats. The adjacent-fit residual uses float FMA. These changes retain the
existing exact-helix and rotation tests without double intermediates.
Seed chord lengths use the same float square-root arithmetic on both paths;
this avoids amplification of CPU/device `hypot` rounding during seed fitting.
Propagation uses float transcendental functions and the standard float
Bethe–Bloch helper. Sine and cosine share float range reduction and FMA
polynomials over the tracking-angle range; device intrinsics handle exceptionally
large arguments. This avoids the library FP64 slow path without the accuracy
loss of using direct intrinsics for ordinary tracking angles. The compiled CUDA
kernels were checked with `cuobjdump`: no FP64 arithmetic or conversion
instructions. HIP compilation/runtime verification still requires an AMD host.

## Reproduce the event comparison

In an O2 environment containing the newly built workflow and PyROOT:

```sh
python3 Detectors/ITSMFT/common/tracking/test/scripts/validateHybridTracking.py \
  --input "$HOME/alice/run3/tracking" \
  --ccdb "$HOME/alice/run3/events/ccdb" \
  --output /tmp/itsmft-gpu-validation \
  --tracker /path/to/build/stage/bin/o2-its-ca-tracker-workflow
```

The script clusters `batch_000` through `batch_019`, then runs full CPU tracking
and GPU tracking with tracklet, cell, neighbour, road and refit validation. Inputs are read-only; products,
logs and `report.json` go under `--output`. Existing local calibration snapshots
supply the dictionary, Alpide parameters and material LUT; each batch's own GRP
and aligned geometry override those entries. Both runs use identical conditions,
198-BC ROFs, a 64-BC bias and four CPU threads. These timing values match the
supplied sample. Use a different output directory for different configurations.

Comparison requires nonempty outputs and matches tracks by their exact ordered
cluster identities within each ROF, allowing chi2-sort permutations. ROF
metadata, track metadata, timestamps and MC labels must match exactly, including
multiplicities. Fitted parameters, both covariances and chi2 use finite-value
absolute-plus-relative tolerances:

| Quantity | Absolute tolerance | Relative tolerance |
| --- | ---: | ---: |
| Reference coordinate | 1e-5 | 1e-4 |
| Frame angle | 1e-6 | 1e-4 |
| Track parameters | 1e-5 | 2e-3 |
| Covariance elements | 1e-8 | 2e-2 |
| Chi2 | 1e-2 | 5e-3 |

`test/scripts/testTrackingComparison.py <cpu-output.root>` exercises permitted
numerical differences and within-ROF permutations, and verifies rejection of
large errors, NaNs, changed metadata, clusters and labels.
`--reuse-clusters` reuses clustering products. `--mode async --vertex-source truth`
exercises multiple passes with truth vertices; `--first` and `--last` limit the
batch range. `--keep-going` checks the remaining batches after a comparison
failure, records failures in the report, and still exits unsuccessfully.
`--backend hip` selects the translated backend.

## Validation on 2026-09-22

CUDA 13.0, NVIDIA RTX 3080, GCC 14.2:

| Run | Batches | Matched final tracks | Validated GPU tracklets | Validated GPU cells | Validated GPU neighbours | Validated GPU road extensions | Validated GPU refits |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Sync, diamond vertices | 000–019 | 550,279 of 550,287 | 7,941,854 | 2,453,289 | 1,395,859 | 1,366,330 | 792,984 |
| Async, truth vertices (three passes) | 000 | 20,195 | 971,950 | 350,556 | 199,532 | 71,987 | 40,645 |

Seventeen of twenty sync batches pass the complete comparison. All five GPU
stage validations pass in every batch. The final comparison retains these
failures without widening tolerances:

- Batch 007: two entries choose different outer clusters (26,759 tracks in both
  outputs), as in the preceding hybrid implementation.
- Batch 013: three entries choose different outer clusters (28,576 tracks in both
  outputs); competing four-hit tracks have chi2 near 0.2158.
- Batch 011: three entries of one track in overlapping ROFs have outer covariance
  element 9 equal to about 5.7984e-6 on CPU and 6.00417e-6 on GPU, exceeding the
  2% relative tolerance. Their parameters and chi2 pass. Both outputs contain
  23,181 tracks.

All other final tracks pass the identity, numerical, metadata and label checks.
The async batch passes the complete comparison. Floating-point bytes and
chi2-sort order are not required to match. The shared float seed arithmetic can
also change a few CPU selections relative to the previous `hypot` implementation.
Reports and logs are in
`/tmp/itsmft-gpu-validation` and `/tmp/itsmft-gpu-validation-async` on the test host.

All fourteen targeted CTest suites passed. They cover CPU cell cuts, triplet
numerics, propagation/material physics, CPU cell orchestration, CUDA parity,
boundaries and allocation recovery,
tracker failure contracts, workflow configuration and DPL publication. GPU orchestration covers cylinder and disk surfaces, hole edges, signed
slopes, timestamps and lookup tables. The GPU cell tests also cover disjoint
time intervals, degenerate geometry, empty inputs, invalid ranges, and the
real-event small-curvature rounding regression. Neighbour tests cover ordering
across multiple blocks, timing and shared-hit cuts, sparse paths, chi2 cuts,
invalid fits/ranges, empty inputs and allocation-failure recovery. Road tests
cover branching across multiple blocks, cylinder/disk and mixed propagation,
zero and signed fields, nonzero material, timing and level rejection, lookup
failures and allocation recovery. The batch_015 material-rounding fixture
checks inverse momentum within the seed tolerance. Refit tests cover cylinder, disk and mixed-surface
seeds, positive/negative/zero fields, material, holes, optional third
legs, MinPt rejection, malformed inputs and allocation-failure recovery.
The shared refit rejects non-finite and negative measurement variances before
propagation. Workspace tests additionally check multi-block device scans,
concurrent calls, allocation reuse, cell-cache invalidation at the same host
addresses, and resident road-state handoff. CUDA Compute Sanitizer reports zero memory errors for the new
seed/refit tests. HIP build/runtime validation remains pending:
this host has no ROCm/HIP toolchain or AMD device.

## Timing CPU and GPU tracking

Use the prepared event fixtures from the validation script, but disable all
validation for timing:

```sh
python3 Detectors/ITSMFT/common/tracking/test/scripts/benchmarkTracking.py \
  --fixtures /tmp/itsmft-gpu-validation \
  --output /tmp/itsmft-tracking-benchmark-sync \
  --tracker /path/to/build/stage/bin/o2-its-ca-tracker-workflow
```

The script defaults to three repetitions of all twenty batches with four CPU
threads, alternates backend order, and runs workflows sequentially. Initial
CPU/GPU runs prime filesystem caches and are excluded. Each measured run starts
a fresh workflow and GPU context: these are first-timeframe measurements, not
steady-state kernel throughput. For async timing add
`--mode async --vertex-source truth`; `--first`/`--last` limit the sample.

`report.json` records commands, counts, individual timings and per-batch medians.
The tracker timer covers preparation, tracking, GPU allocation/transfers and
synchronization, and MC labeling. Separate workflow wall times also include
process startup, conditions loading, cluster input and ROOT output. GPU validation
would repeat CPU calculations and must not be included in a speed comparison.

### Baseline measurements before workspace reuse (2026-09-22)

Ryzen 9 5950X (four tracking threads), RTX 3080, CUDA 13.0, driver
580.105.08. Validation disabled; ROOT output and MC labeling enabled.

| Workload | CPU tracking | CUDA tracking | CUDA / CPU |
| --- | ---: | ---: | ---: |
| Sync, first TF, all 20 batches: sum of per-batch medians | 7.502 s | 8.791 s | 1.172 |
| Sync, first TF: median across the 20 batch medians | 374.77 ms | 438.97 ms | 1.171 |
| Sync, repeated batch 000, excluding first TF | 416.16 ms | 338.72 ms | 0.814 |
| Async, first TF, batch 000 with truth vertices | 523.38 ms | 612.01 ms | 1.169 |

Each first-TF result uses three repetitions with alternating backend order.
The warmed result uses six identical copies of batch 000 in each process,
three processes per backend, excluding the first TF: 15 measurements per
backend. It measures this repeated input, not a diverse stream of events.
Thus CUDA is 17.2% slower for the sync first-TF comparison, but 18.6% lower
latency (1.23x throughput) for the warmed batch-000 sample. Mean workflow wall
time, using per-batch medians, is 3.757 s for CPU and 4.005 s for CUDA in sync.
Every repetition produced the expected nonzero track counts.

An independent Nsight Systems profile of CUDA batch 000 (first TF) measured:

| Operation | GPU execution time |
| --- | ---: |
| Tracklet kernels | 0.429 ms |
| Cell kernels | 0.687 ms |
| Neighbour kernels | 0.480 ms |
| Road kernels, including seed construction | 3.265 ms |
| Final refit kernels | 0.728 ms |
| All kernels, 84 launches | 5.590 ms |
| Host-to-device transfers, 709.668 MB | 39.712 ms |
| Device-to-host transfers, 63.284 MB | 3.198 ms |

The same profile recorded 241 `cudaMalloc` and 241 `cudaFree` calls, consuming
35.858 ms of CPU API time together, plus 49 stream creation/destruction pairs.
The 231 `cudaMemcpyAsync` calls consumed 56.005 ms of CPU API time. CPU API and
GPU execution times overlap and must not be added. The profiled tracker took
519.92 ms; use the unprofiled measurements above for the speed comparison.
These baseline measurements motivated the workspace reuse and transfer
reductions described above.

Raw logs and reports are under `/tmp/itsmft-tracking-benchmark-sync`,
`/tmp/itsmft-tracking-benchmark-stream` and `/tmp/itsmft-tracking-benchmark-async`.
The sync directory also contains `timings.csv` and `environment.json`.
The profiler capture is `/tmp/itsmft-tracking-profile.nsys-rep`, with exported
summaries in `/tmp/itsmft-profile-stats.csv`. The warmed fixture is a six-entry
ROOT tree made by repeating the batch-000 cluster tree, with unchanged
geometry and conditions; the source events are untouched.

### Workspace optimization validation

The workspace implementation preserves the preceding GPU results exactly:
all final track, cluster-index, ROF and MC-label branches are bit-for-bit
identical, including ordering, across batches 000–019 and async batch 000.
The existing CPU/GPU comparison flags in batches 007, 011 and 013 remain
unchanged; no numerical tolerance was widened. All fourteen targeted suites
pass. Compute Sanitizer reports zero errors and zero leaked bytes for all five
GPU suites, and CUDA disassembly still contains no FP64 instructions.

### Timing after workspace reuse (2026-09-22)

Same machine, four tracking threads, validation disabled, three repetitions,
and the same event/condition fixtures as the baseline. CUDA initialization
remains inside the first tracking call; the gain does not rely on moving setup
outside the timer.

| GPU workload | Before | After | Time reduction |
| --- | ---: | ---: | ---: |
| Sync, first TF, sum of 20 per-batch medians | 8.791 s | 8.046 s | 8.5% |
| Sync, warmed repeated batch 000 | 338.72 ms | 289.68 ms | 14.5% |
| Async, first TF, batch 000 | 612.01 ms | 549.80 ms | 10.2% |

The contemporaneous CPU measurements were 7.418 s for the twenty sync batches,
434.11 ms for warmed batch 000, and 523.21 ms for async batch 000. The warmed
sample therefore gives about 1.50x CPU throughput; first-TF GPU latency remains
about 8.5% above CPU for sync and 5.1% above CPU for async. The warmed input is
one repeated batch, so this is not a claim about arbitrary event streams.
All 120 sync measurements and the stream/async samples preserve the expected
track counts. Their GPU output branches are bit-for-bit identical to baseline,
including all six timeframes in the repeated-input test.

The independent batch-000 profile shows the data-movement reduction:

| Per first timeframe | Before | After |
| --- | ---: | ---: |
| Host-to-device bytes | 709.668 MB | 211.890 MB |
| Device-to-host bytes | 63.284 MB | 51.429 MB |
| Copy calls | 231 | 182 |
| Device allocation calls | 241 | 48 |
| Stream creations | 49 | 1 |

In the six-timeframe profile, all 48 allocations occur before the first
frame's last kernel; subsequent frames perform no device allocations.
The remaining host scheduling and intermediate downloads still prevent full
event residency on the device. Prefix scans add GPU kernels but eliminate
host offset arrays and their uploads.

Updated logs/reports are in `/tmp/itsmft-workspace-benchmark-sync`,
`/tmp/itsmft-workspace-benchmark-stream` and `/tmp/itsmft-workspace-benchmark-async`;
the sync directory includes per-batch `timings.csv` and hardware metadata.
Profiles are `/tmp/itsmft-workspace-profile.nsys-rep` and
`/tmp/itsmft-workspace-stream-profile.nsys-rep`. Exact output comparisons are
recorded in `/tmp/itsmft-workspace-final-exact-comparison.log`.

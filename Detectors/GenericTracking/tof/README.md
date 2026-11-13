# ALICE3 TOF matching

The optional `--tof-matching` stage of `o2-gtrk-trk-reco-workflow` consumes
`GTRK/TRKTRACKS` and `TF3/CLUSTERS`, and publishes `GTRK/TRKTOFMATCH`.
It runs after tracking on the host, including when tracking itself uses a GPU.

Matching has three passes, with exclusive ownership of each cluster:

1. Match tracks to the outer TOF using spatial chi2 and compatibility with the
   track's coarse readout interval under at least one mass hypothesis.
2. For tracks with an outer match, match inner clusters using spatial chi2 and
   the measured outer-minus-inner time compared with the expected flight time.
3. Match tracks without an outer match to the remaining inner clusters.

Each pass sorts accepted candidates by score, input track row and original
cluster index, and assigns them greedily. The second pass scores spatial plus
minimum timing chi2 across the compatible hypotheses. The first and third
passes score spatial chi2. A track whose outer match has no compatible inner
cluster retains its outer timing; it does not enter the third pass.

## Times and output

Times use an integer absolute BC and a float offset in nanoseconds from the
start of that BC. BC differences are formed before conversion to double;
absolute timestamps are never flattened into float. This preserves better than
1 ps representation precision across a span of seconds.

The time reference is the track's **innermost fitted point**, the reference of
`TrackRecord<TRKFT3>::state`. It is not a vertex time. The inner flight can be
negative when that point is outside the inner TOF. Electron, muon, pion, kaon
and proton hypotheses are kept separately, in that order (`kPDG` in `Match.h`).

There is one `tof::Match` per input track, in input order, including unmatched
tracks and failed refits. `inner` and `outer` address the original timeframe's
flat TF3 cluster array; -1 means unmatched. Each hypothesis has a validity flag,
BC/fine reference time, sigma, flight times and, for two matches, timing chi2.
The spatial chi2 and accepted-candidate count are retained for each layer.

Flight times integrate `ds/(beta*c)` in bounded helix steps, using the existing
material energy-loss model. Spatial transport carries covariance and scattering.
Outer spatial matching starts from the outward fit; inner matching starts from
the inward fit. Flight-time derivatives with respect to the inward state give
the shared track error. This correlation is included in the inner/outer timing
gate and in the combined time, rather than counting it twice as independent
noise. Active-area acceptance is checked by the spatial prediction. The
independent trajectory used to integrate flight time may cross the same plane
outside the chip; it does not apply a second active-area veto. No TOF spatial
measurement is added to the tracker fit.

## Inputs and running

The temporary cluster bridge is isolated in commit `b2c6a695c6`, for replacement
by Mario's BC/fine-time implementation. It changes the cluster ROOT version to
3. Existing clusters with integer-truncated times must be regenerated from
digits. The ROOT reader rejects older cluster versions. `src/Clusters.cxx` is
the matching adapter to the cluster API.

The first implementation supports segmented inner and outer barrel TOF. It
uses each chip's placed plane and active bounds, and the cluster pattern's
pixel centre of gravity. Clusters without a usable pattern are counted as
skipped; retained hits keep their original indices. Pixel variances are
pitch-squared/12; the single-cluster timing sigma is configurable.

Keep the tracking, TOF, HBF timing, field and geometry from the same production.
The matching device uses the tracking conditions, including `TRK/Param/MatLUT`.
The default tracker LUT ends at r=83 cm and |z|=225 cm, before the full outer
TOF acceptance. Build an extended LUT from the same geometry, for example:

```sh
o2-gtrk-build-matlut -d TRKFT3 -g o2sim_geometry-aligned.root \
  --trk-rmax-cm 100 --trk-zmax-cm 360 -o gtrk_tof_matlut.root
```

Provide that LUT as `TRK/Param/MatLUT` through the production's existing CCDB or
snapshot/remapping setup. The new extent options preserve the default tracker
region unless explicitly supplied. The matcher rejects TOF sensor positions
outside the LUT coverage.

Add these options to the existing tracking command, keeping its input and
condition arguments:

```sh
o2-gtrk-trk-reco-workflow --tof-matching \
  --tof-cluster-infile tf3clusters.root --tof-time-resolution-ps 20 -b
```

Alternatively, `--tof-clusters-from-upstream` consumes `TF3/CLUSTERS` from a
preceding workflow instead of adding the ROOT reader. The default output is
`o2_gtrk_tof.root`, tree `o2sim`, branch `TRKTOFMatch`; `--disable-root-output`
disables this writer. Without `--tof-matching`, tracking adds no TOF devices.

The initial cuts are `--tof-max-spatial-chi2 25`, `--tof-max-time-chi2 25`, and
`--tof-time-margin-ns 0.1`. Flight integration uses `--tof-max-step-cm 2` and
`--tof-max-path-cm 1500`. These are starting values, not tuned physics cuts.

## Validation and current scope

Build targets include `O2lib-GTrkTOF`, `O2test-generic-tracking-tof-matching`,
`O2test-iotof-cluster-timing`, `O2exe-gtrk-trk-reco-workflow` and
`O2exe-gtrk-build-matlut`. Run the registered tests with:

```sh
ctest --test-dir BUILD --output-on-failure \
  -R 'test(TOFMatching|ClusterTiming).cxx'
```

Set `GTRK_TOF_GEOMETRY` to a matching ALICE3 geometry to additionally compare the
cluster adapter with TGeo sensor transformations for both barrel layers.
The tests cover BC boundary averaging, 1 ps differences at large BC counters,
ROOT serialization, pass priority, cluster exclusivity, competing tracks,
timing rejection, correlated uncertainties, signed/zero fields, analytic
helix flight times, material, finite acceptance, path failures and recovery of
the reference time through the complete predictor/matcher combination.

This is an exhaustive CPU reference matcher. It has no spatial index or GPU
matching kernels, and uses deterministic greedy assignment rather than a
maximum-cardinality assignment. Transport uses uniform nominal Bz and the
tracker's material approximation. Timing errors include the initial track
covariance and cluster timing error; additional stochastic flight-time
variance generated along the material path is not yet included. Mass
hypotheses share the input fitted momentum rather than refitting tracking for
each species. First-crossing transport and the path limit constrain loopers.

Pile-up MC studies of matching efficiency, purity and time pulls, cut tuning,
material-step convergence on real tracks, and indexing performance remain
necessary before production use. Synthetic and geometry tests do not establish
those results.

## Initial efficiency studies

`study/tofCoverage.py` joins a single timeframe's tracking QA with a source's
TF3 transport hits by source/event/track identity. It reports hit availability
among good, non-clone, valid refits, with generated pT/eta cuts and a separate
physical-primary selection. It rejects multi-timeframe QA files and duplicate
identities. This is a denominator/reach check, not matching efficiency.
Supply the chip partition from the production's TOF geometry; do not infer it
from a different layout. It writes a per-track CSV, a JSON summary and pT plots.
The Python requirements are numpy, uproot and matplotlib (no awkward needed).

`o2-gtrk-prepare-tof-clusters` runs the native IOTOF digitizer and clusterer
on an existing production's source-0 hits and collision context. It loads
`IOTOFBase`, `ChipSpecificsParam` and `TF3DigitizerParam` from the production INI,
preserves MC labels, and writes `tf3digits.root`, `tf3clusters.root`, a topology
map, effective settings and a preparation summary in a new output directory.
Response maps are disabled, matching the native workflow default; the flat
response uses the production settings. The driver requires continuous readout
and unique source-0 events. Use the full production workflow for other sources,
reused events or calibrated response maps.

```sh
o2-gtrk-prepare-tof-clusters --production-config sgn_configuration.ini \
  --geometry o2sim_geometry-aligned.root --context collisioncontext.root \
  --tof-hits sgn_HitsTF3.root --output-dir /path/to/new/tof-output --seed 42
```

`o2-gtrk-tof-hit-study` reruns GTRK on a frozen tracking-input dump and samples
unique good physical-primary tracks using a fixed seed. With `--tof-clusters`,
it reads native TF3 clusters and their MC labels, converts them through the
same adapter as the matching workflow, and runs the actual three-pass matcher.
It uses the stored BC/fine times without added smearing. Truth is used only for
evaluation after matching; a correct match carries the particle's MC identity.
Clusters with multiple distinct identities and assignments to them are counted
separately. The tracker retains its original LUT and configuration; matching
uses TGeo material queries. All cluster clutter is present; reconstructed-track
competition is limited to the sampled tracks.

```sh
o2-gtrk-tof-hit-study --input TRK_tf1.input \
  --geometry o2sim_geometry-aligned.root --tracking-material tracker_matlut.root \
  --tof-clusters /path/to/tof-output/tf3clusters.root \
  --max-tracks 100 --min-pt 0.5 --max-abs-eta 1 --output study/clusters
```

The optional transport-hit proxy mode replaces `--tof-clusters` with
`--tof-hits sgn_HitsTF3.root --context collisioncontext.root`. It converts all
source-0 transport hits into single-pixel measurements and adds 20 ps Gaussian
timing. It does not model detection efficiency, multi-pixel response or cluster
merging, and therefore does not measure digitized-cluster efficiency.

Keep `.ini` and `.truth` beside the frozen input. Set `ALICE3_MAGFIELD_MACRO`
to the production's field, and `--first-orbit` to the timeframe anchor.
Source events reused in several collision occurrences are explicitly rejected.
The default 100 ns hit-time preselection is an additional study selection;
correct hits excluded by it are counted and cause a nonzero study status.
The default kinematic selection is intended for an initial barrel probe.

The denominator for each layer is sampled reconstructed tracks with at least
one usable measurement of their own MC particle in that layer. Correct matches,
wrong matches and purity are separate. The guided-inner denominator is tracks
with a correct outer association and an available inner measurement. The CSV
keeps each truth identity and the per-track decisions. Correct-hit propagation,
mean position inside active bounds, covariance transport, spatial chi2 and
coarse-time compatibility are also tested after matching, independently of
cluster ownership, to distinguish transport/gate losses from assignment losses.
These diagnostic gates use the forward crossing for the selected primary sample.

Before quoting production efficiency, repeat with all reconstructed
competitors, representative pile-up and held-out samples, and validate time
pulls at the explicitly defined track reference.

An initial local pp-TF check used 100 fixed-seed sampled primary refits out of
11,138 eligible reconstructed particles, generated pT >= 0.5 GeV/c and |eta| < 1,
with 243,605 usable synthetic measurements from 264,401 transported hits.
On identical inputs, removing the second chip-bounds veto from the independent
flight-time projection changed correct outer associations from 63/89 to 78/89:
15 recovered identities, none lost. Correct inner associations stayed 90/97,
with one wrong inner association and no wrong outer association observed.
The guided-inner result was 73/77; the remaining-track pass found 17/20.
Outer losses after this fix were six spatial means outside active bounds,
four spatial-chi2 rejections and one flight/coarse-time rejection. These are
small-sample diagnostic results under the approximations above, not production
efficiencies or tuned defaults.

The same 100 tracks were subsequently matched to native digitized TF3 clusters
from the same hits, context and production settings (20 ps resolution, 10 ps
TDC bins, 98% flat efficiency, native dead-time handling). Native clustering
produced 238,845 usable clusters; 91 carried multiple distinct particle labels.
There was no extra timing smearing or truth-based candidate pruning.

| Layer | Tracks with own cluster | Correct matches | Wrong matches, all tracks | Purity |
| --- | ---: | ---: | ---: | ---: |
| Outer | 88 | 77 (87.5%) | 1 | 77/78 (98.7%) |
| Inner | 94 | 88 (93.6%) | 1 | 88/89 (98.9%) |

The wrong outer match belonged to a track without its own outer cluster.
There were no assignments to ambiguous clusters and no true clusters lost by
the study's time preselection. The guided-inner result was 71/75. These are
conditional efficiencies among tracks with an available cluster, in a sample
of 100 reconstructed primary tracks. Only those 100 tracks compete for the
complete timeframe's clusters; this remains a preliminary efficiency probe.

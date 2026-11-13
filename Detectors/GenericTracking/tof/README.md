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
noise. No TOF spatial measurement is added to the tracker fit.

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

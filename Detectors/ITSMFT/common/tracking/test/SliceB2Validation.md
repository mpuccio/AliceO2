# Slice B2 geometry and observational replay validation

This records the reproducible validation of the geometry-backed ITS/MFT
surface catalogs and the observational MFT production configuration. It does
not mark Gate 2 complete.

The checks used a fresh `RelWithDebInfo` build configured from this worktree at
`/private/tmp/o2-layout-providers-build`. Every O2, ROOT, and workflow command
was run through
`.agents/skills/alice-o2-environment/scripts/run-in-o2-env.zsh`.

## Real geometry

The committed macro was run against the retained geometry fixture without
modifying it:

```sh
root -l -b -q \
  'Detectors/ITSMFT/common/tracking/test/validateGeometrySurfaceCatalogs.C("/private/tmp/o2-itsmft-gate0/fixture-20ev/o2sim")'
```

It loaded both detector geometry singletons, populated their L2G caches, and
verified dense catalogs with 7 ITS layers and 10 MFT disk planes, finite
ordered bounds, zero flags, populated chip buckets, and documented
plausibility windows. The reported values, in centimetres, were:

| Detector surface | reference coordinate | radial minimum | radial maximum |
|---|---:|---:|---:|
| ITS 0 | 2.32597 | 2.23456 | 2.60234 |
| ITS 1 | 3.13535 | 3.00804 | 3.39994 |
| ITS 2 | 3.91624 | 3.78343 | 4.16012 |
| ITS 3 | 19.5882 | 19.3533 | 19.8943 |
| ITS 4 | 24.5272 | 24.3032 | 24.8082 |
| ITS 5 | 34.3546 | 34.1432 | 34.6070 |
| ITS 6 | 39.3106 | 39.1032 | 39.5540 |
| MFT 0 | -45.2889 | 2.30291 | 12.2258 |
| MFT 1 | -46.7111 | 2.30291 | 12.2258 |
| MFT 2 | -48.5889 | 2.30291 | 12.2258 |
| MFT 3 | -50.0111 | 2.30291 | 12.2258 |
| MFT 4 | -52.3889 | 2.30291 | 13.6984 |
| MFT 5 | -53.8111 | 2.30291 | 13.6984 |
| MFT 6 | -67.6889 | 3.40291 | 16.7088 |
| MFT 7 | -69.1111 | 3.40291 | 16.7088 |
| MFT 8 | -76.0889 | 3.80291 | 17.3105 |
| MFT 9 | -77.5111 | 3.80291 | 17.3105 |

The fixture contained 24,120 ITS chips and 936 MFT chips. MFT geometry
construction and L2G extraction use the common `SegmentationAlpide` sensor
definition, confirming that the provider's standard asymmetric ALPIDE active
footprint assumption applies to MFT.

## Common-CA observational replay

The existing fixture checks confirmed the retained inputs were present and
non-empty and that the HBF configuration contains run 303000. No simulation
or digitization was regenerated. The committed
`replay_tracking_common_ca.sh` was run single-threaded into the fresh external
directory
`/private/tmp/o2-itsmft-gate0/replay-common-ca-slice-b2-final`, using
timestamp 1784207296000 and this worktree's staged binaries.

`extract_metrics_common_ca.C` produced exactly the accepted documented
metrics:

- input clusters 1679; input ROFs 2304; output tracks 91;
- clusters/track: n 91, mean 6.50549, median 6, minimum 5, maximum 10;
- chi2: n 91, mean 0.222481, median 0.0835225, minimum 0.000986104,
  maximum 2.01957;
- MC reconstructable 141; matched 77; efficiency 0.546099;
- fake tracks/rate 0/0; clone tracks/rate 0/0;
- `trackContentHash` `826dc653cd936a472929c600c97c140b`.

The hash covers the extractor's documented ordered per-track tuple only. This
validation does not claim byte-identical ROOT files or equality of fields
outside the documented metrics and hash. The successful geometry-layout
diagnostic appeared exactly once in `mft_ca_replay.log`; no layout-failure
diagnostic appeared. The extracted values therefore confirm that constructing
the observational layout leaves the existing legacy CA result unchanged.

## Geometry-update limitation

Invalidation after a successful detector-specific `GeometryTGeo` adoption and
cache fill prevents a previously built catalog from being presented as
current. This is not hot-reload support: detector-specific `GeometryTGeo`
adoption currently fatals on replacement, and GRP-managed aligned geometry
has no singleton/cache rebuild path. Genuine hot geometry reload remains
outside Slice B2; no `GEOMALIGN` invalidation hook is added.

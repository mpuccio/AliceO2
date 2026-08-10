# Local triplet fitting R2 validation

> Historical characterization of the fixed-point material experiment. The
> implementation is preserved at `ae27032797`; it is not the active fast
> estimator.

- Date: 2026-08-10
- Package: `daily-20260717-0700-local1`
- Build: `/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`
- Fixture: `/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717`
- Artifacts: `/Users/mpuccio/alice/run3/O2-validation-artifacts/itsmft/triplet-fitting-r2-20260810-sqUdOL`
- Design: [momentum-dependent triplet material](../design/0019-momentum-dependent-local-triplet-material.md)

## Scope

R2 adds the family-leaf surface normal, incidence-scaled path material, fitted
momentum, shared material-kernel evaluation, and safeguarded fixed-point
iteration to the standalone R1 primitive. Temporary environment-gated
instrumentation evaluated the fit beside the existing cell gates and wrote
candidate TSV files. It never affected a decision and was removed before the
final build. `TrackerTraits` has no `TripletFitting` include or call site.

## Candidate cut-flow

The diagnostic traversed every time-compatible pair of tracklets considered
for a cell. The longitudinal and transverse columns below are independent
decisions; `both` is their intersection. Existing `CellSeed` construction is
unchanged.

| Stage | ITS | standalone MFT |
|---|---:|---:|
| Candidate tracklet pairs | 1,699 | 8,084 |
| Existing transverse gate passes | 1,267 | 1,024 |
| Existing longitudinal gate passes | 1,502 | 1,063 |
| Both existing direction gates pass | 1,223 | 582 |
| Existing accepted cells | 1,211 | 559 |
| Material-aware fit succeeds, all pairs | 1,550 | 816 |
| Below scalar-kernel momentum range | 1 | 1,520 |
| Other scalar material failure | 146 | 5,733 |
| No fixed-point convergence | 2 | 15 |

Among existing cells, the fit succeeds for `1,209/1,211` ITS cells and
`512/559` MFT cells. The remaining two ITS and one MFT non-convergences are
unlabelled combinations. The 47 other MFT failures comprise one candidate
below the scalar kernel's momentum range and 45 candidates that fail its
material traversal plus one non-convergence. Five material-failing MFT cells
have a common valid MC label. Four have MC `pT <= 0.00497 GeV/c`; the fifth
has MC `pT=0.507438 GeV/c` but the local three-hit estimate is
`0.000777 GeV/c`. This is direct evidence that a material-domain failure
cannot yet be promoted blindly to a production rejection.

A hypothetical common `chi2 < NSigmaCut^2 = 25` applied only after successful
fitting would retain 1,201 of the 1,209 successfully fitted ITS cells and all
512 successfully fitted MFT cells. This is characterization, not a proposed
selection: fit-domain failures and candidate competition still require a
production contract.

## Momentum and pull characterization

Cluster labels were associated to the fixture `MCTrack` by source, event, and
track ID. For successfully fitted, currently accepted cells:

| Quantity | ITS | standalone MFT |
|---|---:|---:|
| Associated cells | 1,088 | 448 |
| median fitted-pT relative residual | -0.0064 | -0.6183 |
| p10 / p90 fitted-pT relative residual | -0.3518 / +0.5989 | -0.9478 / +1.6250 |
| median absolute-curvature pull | +0.0470 | +0.0295 |
| p10 / p90 absolute-curvature pull | -1.3250 / +0.8389 | -0.0259 / +0.1029 |

The pull uses the fitted curvature variance and the production-vertex MC
momentum, so it is a characterization rather than a calibrated detector pull.
MFT's small pulls coexist with a poor point estimate because its short radial
lever arm and material term produce a very broad curvature uncertainty. The
three-hit estimate is consequently not a sound common `TrackletMinPt`
replacement for disks in this state.

## Cost

Per-candidate wall-clock timing covered observation-ready material iteration
but excluded TSV locking and formatting. Successful-fit distributions were:

| Host cost | ITS | standalone MFT |
|---|---:|---:|
| median | 0.833 us | 1.250 us |
| p90 | 1.417 us | 3.750 us |
| p99 | 4.417 us | 11.458 us |
| maximum | 24.625 us | 22.042 us |

The successful iteration-count medians are 3 and 4; p99 values are 21 and 52.
The tail is why this primitive should not replace the cheap preselection
without a separate throughput and placement decision.

## Non-invasive replay evidence

The diagnostic replay retained the current production products exactly:

| Product | Tracks | Content hash | Efficiency | Fake | Clone |
|---|---:|---|---:|---:|---:|
| ITS standalone | 184 | `932e38144a0575f495b29b532012f468` | 0.830986 | 15 / 0.0815217 | 1 / 0.00704225 |
| MFT standalone | 69 | `f6dee3f7a5f7def6b55900dbac734ef0` | 0.522936 | 2 / 0.0289855 | 0 / 0 |

Fixture checksums passed 43/43 before and after. After removal of diagnostics,
the focused fit/header tests passed 2/2 and the complete pinned serial label
passed 94/94 in 94.95 seconds. The final source guard proves that no
`TripletFitting` include or call remains in `TrackerTraits`.

## Verdict

**Standalone R2 primitive accepted; production integration deferred.** The
neutral material contract and solver are suitable for continued triplet R&D,
but current evidence is insufficient to replace the two cheap gates or apply
the fitted momentum as a common disk/cylinder physics selection. The next
step should extend the lever arm through multi-triplet/track fitting and then
repeat candidate-level pull, cut-flow, and timing validation.

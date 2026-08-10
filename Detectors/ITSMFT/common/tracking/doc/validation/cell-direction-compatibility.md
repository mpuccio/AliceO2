# Shared cell-direction compatibility validation

- Status: safe implementation candidate; unified physics sign-off required
- Date: 2026-08-10
- Production commit: `ec9deeb178`
- Tests and guards: `58d87eb8e2`
- Design: [shared cell-direction compatibility](../design/0016-cell-direction-compatibility.md)
- Complete ITS population inventory: [cell-direction population delta](cell-direction-population-delta.md)
- Pinned O2 package: `daily-20260717-0700-local1`
- Build: `/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`
- Fresh artifacts: `/private/tmp/o2-cell-direction-harmonization-20260810-E2xZk8`

## Implementation evidence

The common core now constructs three neutral `(r,z)` observations and tests

```text
K = (z0-z1)(r1-r2) - (z1-z2)(r0-r1)
chi2 = K^2 / Var(K)
chi2 < NSigmaCut^2
```

The analytic variance includes `Var(r)`, `Cov(r,z)`, and `Var(z)` for every
observation. Disk leaves project the decoded global `x/y` covariance onto
`r` and use exact descriptor `z`. Cylinder leaves use exact descriptor `r`
and tracking-frame measured `z` with `Cvv`. Invalid, non-finite, non-PSD, or
zero-final-variance inputs reject transactionally.

Focused tests cover the exact analytic derivative/covariance contraction
with non-zero cross terms, global-Z rotation invariance, disk radial
projection, cylinder fixed-radius projection, zero residual for both
families, the common `NSigmaCut` boundary, and closed transactional failure.
Source guards prove that:

- `CellFinding.h` and `TrackletFinding.h` are absent;
- `computeDiskCellDirectionCompatibilityChi2` and
  `cellDirectionsAreCompatible(SurfaceKind,...)` are absent;
- the shared compatibility body contains no `SurfaceKind`, family name,
  switch, or `if constexpr`;
- `TrackerTraits` constructs three observations and makes one compatibility
  call without reading `Tracklet::tanLambda`;
- cylinder and disk cell seeds continue through the same generic candidate
  loop.

## Test and environment campaign

After rebuilding all affected test binaries against the shortened compact
parameter layout, the clean complete serial campaign passed:

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 93
Total Test time (real) = 105.07 sec
```

The final focused rerun passed 7/7. `git diff --check` passed and
`git clang-format --diff` reported no changes for both production and test
commits.

Strict preflight resolved the exact pinned package, found the AliEn token
valid from `2026-08-07 05:11:10 UTC` through
`2026-09-05 09:33:43 UTC`, confirmed that
`ALICEO2_CCDB_NOTOKENCHECK` was unset, and found every dataset reported by
`geant4-config --datasets`. Replays used condition cap `1784207296000`.
The fixture checksum manifest passed 43/43 before replay and 43/43 after all
primary, repeat, and diagnostic runs.

## Stage-level cause

Temporary environment-gated instrumentation was used only for diagnosis,
then removed before the final build. The first ITS change is cell creation;
tracklet formation is exactly unchanged.

| ITS stage, standalone and combined component | Prior fixed pull | Shared collinearity | Change |
|---|---:|---:|---:|
| Accepted tracklets | 1,886 | 1,886 | 0 |
| Accepted cells | 1,455 | 136 | -1,319 |
| Neighbour links | 1,049 | 29 | -1,020 |
| Accepted/published tracks | 212 | 17 | -195 |

A temporary diagnostic selector using the former cylinder expression
`abs(tanLambda01-tanLambda12)/0.007 < 5` reproduced all prior stage counts
and the complete prior 212-track TSV byte-for-byte. This proves that graph,
timing, tracklet search, neighbour linking, road finding, refit, and
publication do not create the initial divergence. The changed cell gate does;
later reductions and substitutions are downstream consequences.

The disk leaf is unchanged from the preserved uncommitted starting
correction, and MFT products remain field-identical. Current component counts
are:

| MFT composition | Tracklets | Cells | Neighbour links | Tracks |
|---|---:|---:|---:|---:|
| Standalone | 2,723 | 601 | 314 | 69 |
| Combined MFT component | 4,705 | 631 | 315 | 100 |

The combined aggregate is 6,591 tracklets, 767 cells, 344 neighbour links,
and 117 tracks; subtracting the measured ITS component gives the MFT row
above. No neighbour, road, graph, hole, or timing code changed in this slice.

## Numerical representative cases

The acceptance threshold is `chi2 < 25`. `oldPull` below is characterization
of the removed cylinder decision and is not present in production.

| Candidate refs, inner to outer | Old tanLambda pair / pull | Descriptor `(r,z)` observations | `K` | `Var(K)` | chi2 | New decision |
|---|---|---|---:|---:|---:|---|
| `639,812,1010` | `(-0.3741004,-0.3744212)` / `0.0458317` | `(19.58824,-7.349940)`, `(24.52716,-9.196802)`, `(34.35460,-12.877872)` | `0.0305864` | `7.53209e-4` | `1.24205` | accept |
| `812,1010,1294` | `(-0.3744212,-0.3750362)` / `0.0878487` | `(24.52716,-9.196802)`, `(34.35460,-12.877872)`, `(39.31064,-14.737780)` | `0.0345615` | `1.26357e-3` | `0.945334` | accept |
| `470,639,812` | `(-0.3716724,-0.3741004)` / `0.346856` | `(3.916242,-1.423988)`, `(19.58824,-7.349940)`, `(24.52716,-9.196802)` | `-0.323775` | `9.70130e-4` | `108.058` | reject |
| `273,470,639` | `(-0.3718146,-0.3716724)` / `0.0203082` | `(3.135354,-1.166676)`, `(3.916242,-1.423988)`, `(19.58824,-7.349940)` | `0.594916` | `2.58556e-4` | `1368.85` | reject |
| `24,273,470` | `(-0.3714343,-0.3718146)` / `0.0543254` | `(2.325965,-0.930319)`, `(3.135354,-1.166676)`, `(3.916242,-1.423988)` | `0.0236971` | `1.71507e-6` | `327.421` | reject |
| `117,327,515` | `(-2.060187,-2.059882)` / `0.0435625` | `(2.325965,-4.604550)`, `(3.135354,-6.598718)`, `(3.916242,-7.959709)` | `-0.455653` | `1.51942e-6` | `136644.5` | reject |

The first two accepted rows form current track #0
`1294,1010,812,639`. It is the four-hit prefix of old track #23
`1294,1010,812,639,470,273,24`; the third row is its first rejected
extension cell. Old track #0 `693,515,327,117` loses both cells; the last row
is one of them.

The old stored tracklet slopes use per-cluster global radius. The new
cylinder observations use the descriptor's nominal fixed layer radius, as
the requested measurement model specifies. On real stave geometry those
radii are not identical. The last row illustrates the effect: the old
slopes agree near `-2.06`, while slopes reconstructed from the displayed
fixed-radius observations are about `-2.464` and `-1.743`. No coordinate
conversion or denominator error is hidden here; the fixed-surface model and
small measured-Z uncertainty make the normalized pull large.

## Population characterization

Standalone and combined ITS products are field-identical both before and
after this change. Therefore one complete 212-to-17 association inventory
covers both compositions.

| ITS product | Tracks / hash | Matched / 142 | Efficiency | Fake / rate | Clone / rate | Hit distribution |
|---|---:|---:|---:|---:|---:|---|
| Prior | 212 / `46913a67a7e2fe7462e29df0db264fa8` | 135 | `0.950704` | 16 / `0.0754717` | 0 / 0 | 46x4, 16x5, 29x6, 121x7 |
| Current | 17 / `b194927a8306cfac57ada54fc55c0c6f` | 10 | `0.0704225` | 3 / `0.176471` | 0 / 0 | 13x4, 3x5, 1x7 |

Ordered ROF plus ordered cluster references retain 2 tracks exactly, remove
210, and add 15. Every addition is a shorter or changed candidate sharing a
contiguous reference run with a prior track; the changed gate alters road
availability and candidate competition rather than creating unrelated hits.
All 227 association rows, including hit count, ROF/BC/orbit, raw label
validity/fake bit, chi2, ordered references, and best
prefix/suffix/internal relation, are in the linked population appendix.

MFT is unchanged from the preserved disk-collinearity starting point:

| MFT product | Tracks / hash | Matched / 109 | Efficiency | Fake / rate | Clone / rate |
|---|---:|---:|---:|---:|---:|
| Standalone | 69 / `f6dee3f7a5f7def6b55900dbac734ef0` | 57 | `0.522936` | 2 / `0.0289855` | 0 / 0 |
| Combined | 100 / `98f9730d6fca9e738c3f20afc66d296d` | 78 | `0.715596` | 3 / `0.03` | 6 / `0.0550459` |

The current MFT TSVs are byte-identical to their immediate predecessors.
Standalone and combined still share 68 of 69 standalone sequences; combined
uses a five-hit suffix where standalone publishes the corresponding six-hit
extension, and the wider combined population remains explained principally
by the intentional common `MinTrackLength=4` policy. Equality is not an
acceptance requirement.

Fresh clean repeats are field-identical to every primary product for ITS
standalone, MFT standalone, combined ITS, and combined MFT. This comparison
includes ordered references, ROFs, labels, fake bits, chi2, and hit counts.

## Verdict

**Safe implementation candidate; unified physics sign-off required.** The
code implements the requested common mathematical and covariance contract,
removes both fixed family direction scales, keeps family knowledge in
descriptor-selected observation leaves, and introduces no detector/source
dispatch. MFT is stable. The severe deterministic ITS efficiency loss is
fully localized to the requested fixed-radius, measurement-uncertainty-only
cell statistic and is recorded rather than tuned away. Whether nominal
cylinder surfaces provide a sufficiently faithful physical observation for
this selection is the explicit question for later unified physics review.

# Adjacent-triplet neighbour-fit validation

- Date: 2026-08-10
- Package: `daily-20260717-0700-local1`
- Branch: `codex/itsmft-triplet-rnd-scratch`
- Parent candidate: global measurements and stored triplet factors
- Verdict: **safe candidate physics change**

## Result

Neighbour finding now combines the two stored Equation-19 factors over one
common curvature. The old cylinder/disk state propagation and state-chi2
functions are gone. `TrackerTraits` traverses one neutral scheduled-cell span,
resolves four `GlobalMeasurement`s through the cells' cluster references, and
makes one `fitAdjacentTripletFactors` call. There is no `SurfaceKind`, detector,
or source dispatch in the compatibility operation.

For row order `(theta0,phi0,theta1,phi1)`, the fit constructs

```text
C = D_MS + H V_hit H^T
D_MS = diag(a0, a0/sin(theta0)^2, a1, a1/sin(theta1)^2)
a_i = deltaTheta_i^2 at TrackletMinPt
kappa = -(rho^T C^-1 Psi)/(rho^T C^-1 rho)
chi2 = Psi^T C^-1 Psi - (rho^T C^-1 Psi)^2/(rho^T C^-1 rho).
```

The four-observation `V_hit` retains every `xx,xy,xz,yy,yz,zz` term. The two
shared hits occur in both factor Jacobians, so their contractions generate the
cross-triplet covariance blocks. MS is a physical space-angle variance;
azimuth receives the geometric `1/sin(theta)^2` conversion. There is no
energy-loss term, iterative material correction, uncertainty floor, new magic
threshold, or slope division. The existing common
`maxChi2ClusterAttachment=60` acceptance contract is retained.

## Covariance validity correction

The first replay collapsed ITS from 184 to 30 tracks. Instrumentation showed
that 822 of 898 time-compatible cell pairs failed before chi2: cylinder global
covariances are physically rank deficient, stored as floats, and their zero
minors can become slightly negative when promoted to double. Only 76 pairs
reached the fit, of which 73 passed.

The PSD test now uses a term-scaled tolerance based on float epsilon (16 eps
for 2x2 principal minors, 32 eps for the determinant), matching the actual
`GlobalCovariance3F` storage precision. It does not add covariance or alter the
fit. A focused boundary test accepts a real float-roundoff violation and still
rejects a 1% non-PSD perturbation. With this correction every current replay
observation and adjacent covariance is valid.

## Current neighbour cut-flow

Tracklet formation and cell construction are unchanged by activation of the
adjacent fit. Temporary, subsequently removed instrumentation measured the
first changing stage:

| Leg | time-compatible pairs | valid observations | valid 4x4 fits | chi2 <= 60 | rejected |
|---|---:|---:|---:|---:|---:|
| ITS standalone | 898 | 898 | 898 | 889 | 9 |
| MFT standalone | 372 | 372 | 372 | 353 | 19 |
| combined, total | 1271 | 1271 | 1271 | 1243 | 28 |
| combined ITS component | 898 | 898 | 898 | 889 | 9 |
| combined MFT component | 373 | 373 | 373 | 354 | 19 |

The combined MFT component has one additional input pair, and that pair is
accepted. This is consistent with the pre-existing combined scalar-policy
composition; the compatibility algorithm and threshold are identical.

Representative rejected ITS values span `chi2=62.5184` to `145.467`.
Representative MFT cases are:

| Ordered references | chi2 | decision |
|---|---:|---|
| `152,193,212,237` | 60.0479 | reject |
| `459,489,524,549` | 60.7218 | reject |
| `320,358,383,416` | 66.2705 | reject |
| `306,334,367,396` | 106.594 | reject |
| `424,481,513,544` | 173.649 | reject |

The combined values differ from standalone only at the per-mille level due to
the established combined scalar inputs; all classifications above remain the
same. No fit failed through a singular geometry, invalid factor, non-PSD
covariance, or Cholesky failure.

## Published populations

The correct parent ITS reference is the accepted 184-track
`932e38144a0575f495b29b532012f468` product. A stale intermediate 185-track
artifact exists in the global-measurement scratch directory and is not used as
the comparison reference.

| Product | Parent tracks/hash | Current tracks/hash | Matched / reconstructable | Efficiency | Fake | Clone |
|---|---|---|---:|---:|---:|---:|
| ITS standalone | 184 / `932e38144a0575f495b29b532012f468` | 184 / `e6da9d94faed581d5bff044993698e30` | 118 / 142 | 0.830986 | 15 (0.0815217) | 1 (0.00704225) |
| ITS combined | same as standalone | 184 / `e6da9d94faed581d5bff044993698e30` | 118 / 142 | 0.830986 | 15 (0.0815217) | 1 (0.00704225) |
| MFT standalone | 69 / `f6dee3f7a5f7def6b55900dbac734ef0` | 66 / `32555b198d9b094f3f3600ec619cd2e2` | 54 / 109 | 0.495413 | 2 (0.030303) | 0 |
| MFT combined | 100 / `98f9730d6fca9e738c3f20afc66d296d` | 94 / `96f4c632b7e0111501a63660774480ef` | 75 / 109 | 0.688073 | 3 (0.0319149) | 3 (0.0275229) |

ITS retains 183 ordered cluster-reference sequences exactly. The only
substitution is fake-for-fake in ROF 0:

| Parent | Current | Relation |
|---|---|---|
| 5 hits, label `0x8000000300000aee`, refs `1505,1233,956,780,602` | 4 hits, label `0x800000030000003f`, refs `780,602,434,236` | two shared references; both labels set/fake |

MFT retains 65 sequences exactly. One valid same-label track is shortened by
one leading hit, and three valid-label five-hit tracks disappear:

| Change | ROF | Label (raw decimal) | Seed | Parent references | Current relation |
|---|---:|---:|---:|---|---|
| substitution | 0 | 10737418505 | 510 -> 254 | `545,512,480,424,416,383,358,319` | exact seven-hit suffix `512,480,424,416,383,358,319` |
| loss | 0 | 12884902109 | 992 | `565,544,511,479,423` | none published |
| loss | 0 | 12884902077 | 992 | `247,229,203,186,142` | none published |
| loss | 2 | 30064771095 | 992 | `899,877,860,845,834` | none published |

The lost tracks are not rejected at refit/publication: the only production
change is neighbour-link acceptance, and the output differences therefore
originate at the neighbour graph/level and subsequent road competition. Their
own final four-hit subsequences do not appear verbatim in the rejected-link
list; the loss is an indirect topological consequence of removing other
incompatible links, not a covariance-construction failure hidden later in the
pipeline.

## Standalone versus combined MFT

All 66 standalone tracks have a current combined relation: 65 exact sequences
and one valid same-label six-hit to five-hit suffix substitution. Combined has
28 additional unmatched sequences: 27 four-hit tracks and one six-hit track;
27 have valid non-fake labels and one is fake. Their ROF distribution is
`13,3,6,1,5` over ROFs 0 through 4. This is the established consequence of the
ITS-derived combined application scalar policy, not a graph/start-mask/hole,
timing, source dispatch, or shared-workspace discrepancy. The combined
neighbour cut-flow corroborates it: only one additional MFT adjacent pair is
present, while the same 19 compatibility rejections remain.

The suffix substitution is label `38654707716`, seed `63 -> 31`, from
`1154,1144,1123,1112,1095,1082` to
`1144,1123,1112,1095,1082`. The complete combined-only population is:

| ROF | Hits | Label raw | Fake | Seed | References |
|---:|---:|---:|:---:|---:|---|
| 0 | 4 | 17179869291 | no | 960 | `249,233,205,188` |
| 0 | 4 | 17179873460 | no | 60 | `173,162,119,88` |
| 0 | 4 | 17179869229 | no | 15 | `397,368,335,306` |
| 0 | 4 | 10737418558 | no | 15 | `398,370,337,308` |
| 0 | 4 | 17179869373 | no | 15 | `394,364,333,302` |
| 0 | 4 | 12884901917 | no | 240 | `517,482,431,414` |
| 0 | 4 | 17179869495 | no | 30 | `154,95,64,32` |
| 0 | 4 | 17179869289 | no | 15 | `98,68,34,4` |
| 0 | 4 | 17179869225 | no | 15 | `396,367,334,305` |
| 0 | 4 | 10737418506 | no | 960 | `244,223,177,181` |
| 0 | 6 | 10737418522 | no | 63 | `169,158,118,87,59,25` |
| 0 | 4 | 9223372051887161452 | yes | 15 | `375,342,312,283` |
| 0 | 4 | 15032385605 | no | 960 | `562,539,506,475` |
| 1 | 4 | 19327352842 | no | 960 | `631,628,624,620` |
| 1 | 4 | 19327352842 | no | 30 | `613,608,605,600` |
| 1 | 4 | 21474836581 | no | 960 | `630,627,623,618` |
| 2 | 4 | 23622320206 | no | 15 | `710,701,693,685` |
| 2 | 4 | 27917287551 | no | 30 | `724,704,695,687` |
| 2 | 4 | 27917287543 | no | 30 | `716,706,698,690` |
| 2 | 4 | 30064771095 | no | 15 | `807,798,790,780` |
| 2 | 4 | 27917287485 | no | 960 | `889,881,868,851` |
| 2 | 4 | 23622320605 | no | 15 | `812,802,796,784` |
| 3 | 4 | 32212254951 | no | 960 | `960,958,953,947` |
| 4 | 4 | 38654705803 | no | 15 | `1005,996,984,975` |
| 4 | 4 | 38654705729 | no | 120 | `1175,1155,1141,1124` |
| 4 | 4 | 38654705760 | no | 15 | `1006,997,985,976` |
| 4 | 4 | 38654705840 | no | 960 | `1201,1197,1174,1166` |
| 4 | 4 | 38654705839 | no | 15 | `1127,1114,1096,1083` |

## Tests, determinism, and cost

Focused tests cover the closed-form Equation-19 solution, full shared-hit
cross covariance, space-angle MS projection, exact-helix zero quality,
rotation invariance, float PSD boundary behavior, transactional failures, and
the cylinder/disk shared orchestration path. Source guards remove the legacy
family functions and prohibit surface/detector/source dispatch in the new
method.

The instrumentation-free pinned serial suite passed 94/94 tests in 92.05 s. Standalone
ITS, standalone MFT, and combined ITS+MFT were each repeated from the fixed
fixture; their metric JSON and content hashes are identical between repeats.
The persisted-product comparator also reports field-by-field equality,
including references, ROFs, seed patterns, and labels (2,904 standalone and
4,136 combined MFT float values, with zero maximum delta).
The fixture checksum manifest passed all 43 entries after replay (and passed
before this campaign). Strict CCDB/Geant4 preflight had already confirmed the
package provenance, valid token, unset bypass, and all datasets.

The 20,000-call host microbenchmarks report approximately 504 ns for local
factor construction and 213 ns for the adjacent two-factor solve. These
exclude traversal and cache effects; they are not workflow-throughput claims.
`TripletFitFactor` remains 88 bytes and `CellSeed` remains 220 bytes, so
activation adds no representation cost beyond the already-reviewed payload.

## Verdict and outlook

This is a **safe candidate**, not behavior preservation. The covariance is
mathematically complete for the available independent global-hit errors and
the TrackletMinPt MS hypothesis; validity failures are now zero; the rejected
tails are modest (1.0% ITS, 5.1% standalone MFT); and output changes are
deterministic and attributable to the intended replacement of propagated
family state chi2 by the common triplet-factor chi2.

The MFT efficiency decrease (57 to 54 reconstructed references standalone,
78 to 75 combined) and the combined clone reduction (6 to 3) remain candidate
physics for unified sign-off. The next useful study is neighbour-threshold and
pull calibration against truth, followed by extending the fitted common
curvature/direction state into longer road fits if needed. Energy loss belongs
there or in refit, not in this fast local compatibility. Hole logic, pair-list
topology, road harmonization, and legacy workflow removal remain outside this
slice.

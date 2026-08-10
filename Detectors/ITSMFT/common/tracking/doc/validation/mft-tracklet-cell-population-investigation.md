# MFT population investigation after tracklet/cell harmonization

- Date: 2026-08-10
- Verdict: **requires targeted correction**
- Base: `be70d57e29`
- Candidate: `cca168a093` (`bbf9b27579`, `0bd80c179b`, `cca168a093`)
- Pinned package: `daily-20260717-0700-local1`
- Fixture: `/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/pp-20ev-run303000-seed20260716-daily20260717`
- Current build: `/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`
- Diagnostic artifacts: `/private/tmp/o2-itsmft-mft-population-audit-20260810-zgBqip`

## Executive finding

The disk conversion itself is correct: accepted disk tracklets now store

```text
tanLambda = (source z - target z) / (source r - target r) = dz/dr
```

and old and current runs form exactly the same ordered tracklet pairs in each composition. The new same-radius guard is also inactive on this fixture.

The population loss begins in the unchanged common cell comparison:

```text
abs(first.tanLambda - second.tanLambda) / 0.007 < 5
```

Before the harmonization, the disk value was `(measured delta z) / (nominal disk delta z)`, approximately a signed constant near `-1`; `0.007` therefore acted on deviations of that nominal-Z ratio. It now acts as an absolute uncertainty for physical `dz/dr`. Those are not the same statistical quantity. No covariance or geometry scaling was introduced when the meaning changed.

All 55 historical standalone sequences absent from the current product first lose a required triple at this slope comparison. Fifty-three have valid non-fake MC labels and two have fake labels. The matched-track efficiency falls from `57/109 = 0.522936` to `19/109 = 0.174312`. This is concrete evidence that retaining the numerical `0.007` scale is not yet a justified coordinate-neutral physical selection, even though replacing the old disk pseudo-slope with `dz/dr` is the intended semantic correction.

The bounded correction should keep `tanLambda = dz/dr` and the leaf guard. It should establish an uncertainty-normalized cell compatibility contract, derived from measurements/geometry through the existing operation boundary, without reintroducing detector tuning or a generic disk branch. No correction is implemented in this diagnostic.

## Replay and association method

A source snapshot was made with `git archive be70d57e29` and built separately. Current and base standalone and combined workflows were replayed from the pinned environment against the same durable fixture. Temporary environment-gated diagnostics recorded accepted tracklet pairs, accepted cells, neighbour-link totals, road/refit/publication counts, and disk candidate arithmetic. The instrumentation was removed, the current library was rebuilt cleanly, and focused tests were rerun.

The first base combined attempt stopped before tracking because the base workflow loaded a staged current `libO2ITSTrackingInterface`. Matching base cluster-reader/workflow dependencies were then built; the successful base combined replay is `base-combined-2`. The failed attempt contributes no population evidence.

Exact track identity means equal ROF and equal ordered cluster-reference sequence. Partial associations require the same ROF and are ranked by shared ordered references, with label equality used as supporting evidence. ROF metadata were compared directly from the persisted products.

## Stage-level cut flow

### Standalone MFT

| Stage | Base | Current | Identity/result |
|---|---:|---:|---|
| Accepted disk tracklets | 2,723 | 2,723 | exact pair set unchanged |
| Accepted cells | 1,012 | 211 | 211 retained, 801 removed, 0 added |
| Neighbour links | 315 | 83 | downstream of cell loss |
| Road candidates | 77 | 23 | downstream of cell loss |
| Filtered road seeds | 77 | 23 | no additional current filtering |
| Refit attempts | 77 | 23 | all current roads attempted |
| Successful refits | 72 | 23 | all current attempts succeed |
| Published tracks | 69 | 23 | 14 exact retained, 55 removed, 9 added/substituted |

### Combined MFT component

| Stage | Base | Current | Identity/result |
|---|---:|---:|---|
| Accepted disk tracklets | 4,705 | 4,705 | exact pair set unchanged |
| Accepted disk cells | 1,427 | 214 | 214 retained, 1,213 removed, 0 added |
| Disk neighbour links | 320 | 83 | downstream of cell loss |
| Disk road candidates | 113 | 50 | downstream of cell loss |
| Filtered disk road seeds | 112 | 50 | one base seed filtered; none current |
| Disk refit attempts | 112 | 50 | all current roads attempted |
| Successful disk refits | 103 | 50 | all current attempts succeed |
| Published MFT tracks | 100 | 50 | 24 exact retained, 76 removed, 26 added |

The first changed aggregate is cell construction in both compositions. The graph, timing, neighbour, road, refit, and publication stages do not create the initial divergence.

## Numerical cell cases

The search-window threshold is `transChi2 < NSigmaCut^2 = 25`. In each example both tracklets are accepted in old and current runs; the first missing object is the indicated cell.

| Case | Ordered old track | First triple | Old inputs and decision | Current inputs and decision |
|---|---|---|---|---|
| Typical truncation | `149,133,110,79,48,18` | `18,48,79` | measured `dZ=(1.42220306,1.87779999)` cm; nominal `dZ=(-1.40000153,-1.89999771)` cm; pseudo-tgl `(-1.01585817,-0.988316953)`; pull `3.93446`, pass; search chi2 `(3.82338e-5,8.10570e-5)` | `dR=(-0.126790047,-0.168716908)` cm; tgl `(-11.2169933,-11.1298866)`; pull `12.4438`, reject. Current publishes the five-hit prefix `149,133,110,79,48`. |
| Marginal rejection | `429,412,387,355,324` | `355,387,412` | measured `dZ=(1.42219925,2.37779999)` cm; nominal `dZ=(-1.40000153,-2.40000153)` cm; pseudo-tgl `(-1.01585555,-0.990749359)`; pull `3.58660`, pass; search chi2 `(9.19129e-5,2.41124e-5)` | `dR=(-0.0801587105,-0.133744955)` cm; tgl `(-17.7422924,-17.7786140)`; pull `5.18880`, reject. |
| Small but non-singular radial steps | `152,136,113,82,51,21` | `21,51,82` | measured `dZ=(1.42220306,1.87779999)` cm; nominal `dZ=(-1.40000153,-1.89999771)` cm; pseudo-tgl `(-1.01585817,-0.988316953)`; pull `3.93446`, pass; search chi2 `(0.129545,0.228699)` | `dR=(-0.0219125748,-0.0326793194)` cm; tgl `(-64.9035110,-57.4614182)`; pull `1063.16`, reject. Both radial steps exceed the guard by more than four orders of magnitude. |

Across the 55 first-failing triples, the old normalized difference ranges from `2.49279` to `3.93446` and therefore always passes. The current difference ranges from `5.18880` to `1063.16` (median `10.6994`): 26 are in `[5,10)`, 27 in `[10,100)`, and 2 are at least 100. The strict decision is `< 5`.

These values exclude an accidental sign or reciprocal conversion. The current code and diagnostic values use measured `dz/dr`; the old code uses measured `dz/nominalDeltaZ`. Tracklet identities and transverse search chi2 values are unchanged.

## Same-radius guard audit

The disk leaf computes

```text
deltaR = hypot(source.x, source.y) - hypot(target.x, target.y)
radiusDefined = abs(deltaR) > 1e-6 cm
accept = transChi2 < nSigmaCut^2 && radiusDefined
```

The rejection precondition is therefore finite accepted-window geometry with `abs(deltaR) <= 1e-6 cm`. At that boundary `dz/dr` is undefined or numerically singular; the test is local to disk candidate acceptance.

| Workflow | Candidates inside transverse cut | Guard rejects | Minimum `abs(deltaR)` | Maximum `abs(deltaR)` |
|---|---:|---:|---:|---:|
| Current standalone | 2,723 | 0 | `2.28404999e-4` cm | `6.36911964` cm |
| Current combined | 4,705 | 0 | `2.28404999e-4` cm | `6.36911964` cm |

The closest observed candidate is 228.405 times the tolerance. Since old/current accepted tracklet pair sets are identical and the guard rejects no candidate, none of the removed tracks hit the guard. No distinct-surface, physically valid fixture trajectory is rejected as same-radius.

## Standalone versus combined

The current comparison is exact: all 23 standalone sequences occur in combined, with no standalone removal, and combined has 27 additional sequences. ROF metadata have 0 differences across all 2,304 records.

| Explanation class | Combined-only sequences | Evidence |
|---|---:|---|
| Intentional ITS-derived common scalar policy | 27 | every sequence has exactly four hits; combined `MinTrackLength=4`, standalone MFT `MinTrackLength=5` |
| Graph/start-mask/hole-policy difference | 0 | all 27 have all 3/3 constituent pairs and 2/2 constituent cells in standalone; both use `MaxHoles=0` |
| Timing/ROF composition | 0 | exact ROF/BC/orbit association; 0/2,304 ROF metadata differences |
| Candidate competition in shared workspace | 0 | every standalone sequence survives exactly; no competing combined-only sequence displaces it |
| Unexplained | 0 | all combined-only sequences are accounted for by minimum length |

The combined replay uses the common ITS-derived `PVres=0.05` cm while standalone MFT uses `0.01` cm, producing 4,705 versus 2,723 disk tracklets. The standalone tracklet set is a subset of combined. At cell level all 211 standalone cells are shared and combined has three extra cells (`75,108,132`, `426,479,507`, and `1085,1098,1114`); none is needed by the 27 combined-only tracks. Thus the persisted 23-versus-50 composition is explained specifically by minimum track length, not the wider search, graph, holes, timing, or competition.

## Product metrics

| Product | Base tracks / hash | Current tracks / hash | Matched / 109 | Efficiency | Fake | Clone |
|---|---:|---:|---:|---:|---:|---:|
| Standalone MFT | 69 / `f6dee3f7a5f7def6b55900dbac734ef0` | 23 / `c65075052aa08b7be8a2283cdedc79ea` | 57 -> 19 | 0.522936 -> 0.174312 | 2 -> 0 | 0 -> 0 |
| Combined MFT | 100 / `4dacc3058740cce83656a3bcb71def95` | 50 / `45e46e0475a5d6479054319e8daba176` | 78 -> 41 | 0.715596 -> 0.376147 | 3 -> 0 | 6 -> 1 |

The fresh hashes reproduce the reviewed campaign exactly. All output labels are set and valid. The loss is not merely fake/clone removal: 53 of the 55 absent standalone sequences have non-fake labels.

## Complete standalone loss table

This table covers every base standalone sequence absent from current. Nine rows have a same-label related current prefix, suffix, or internal subsequence; the other 46 have no current candidate sharing a cluster in the same ROF. Every row first fails cell creation at the displayed slope decision.

| Old # | Hits | ROF / BC / orbit | Label | Seed | Ordered references | First missing stage and numerical decision | Best current relation |
|---:|---:|---|---|---|---|---|---|
| 1 | 6 | 0 / 0 / 0 | `0x28000001f` valid=1 fake=0 | `0x3f` | `149,133,110,79,48,18` | cell `18,48,79`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-11.217,-11.1299) pull=12.44; dR=(-0.12679,-0.168717) cm; dZ=(1.4222,1.8778) cm | #1 current prefix; drops 1 trailing; label `0x28000001f`; refs `149,133,110,79,48` |
| 2 | 5 | 0 / 0 / 0 | `0x3000007ea` valid=1 fake=0 | `0x1f` | `127,104,73,43,10` | cell `10,43,73`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-16.1193,-16.0253) pull=13.44; dR=(-0.0882297,-0.117177) cm; dZ=(1.4222,1.8778) cm | none |
| 4 | 5 | 0 / 0 / 0 | `0x4000000c1` valid=1 fake=0 | `0x3e` | `429,412,387,355,324` | cell `355,387,412`; old tgl=(-1.01586,-0.990749) pull=3.587; new=(-17.7423,-17.7786) pull=5.189; dR=(-0.0801587,-0.133745) cm; dZ=(1.4222,2.3778) cm | none |
| 5 | 6 | 0 / 0 / 0 | `0x280000afe` valid=1 fake=0 | `0x3f` | `152,136,113,82,51,21` | cell `21,51,82`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-64.9035,-57.4614) pull=1063; dR=(-0.0219126,-0.0326793) cm; dZ=(1.4222,1.8778) cm | none |
| 6 | 6 | 0 / 0 / 0 | `0x400000027` valid=1 fake=0 | `0xfc` | `199,182,163,155,96,67` | cell `67,96,155`; old tgl=(-1.01586,-0.990749) pull=3.587; new=(-11.5431,-11.3921) pull=21.57; dR=(-0.123208,-0.208724) cm; dZ=(1.4222,2.3778) cm | none |
| 7 | 6 | 0 / 0 / 0 | `0x300000052` valid=1 fake=0 | `0x3f` | `141,121,92,65,31,1` | cell `1,31,65`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-7.89798,-7.9465) pull=6.932; dR=(-0.180072,-0.236305) cm; dZ=(1.4222,1.8778) cm | none |
| 8 | 6 | 0 / 0 / 0 | `0x28000001d` valid=1 fake=0 | `0x3f` | `143,126,101,71,38,9` | cell `38,71,101`; old tgl=(-0.988317,-1.01586) pull=3.934; new=(-12.7458,-12.4888) pull=36.72; dR=(-0.147327,-0.113878) cm; dZ=(1.8778,1.4222) cm | none |
| 9 | 9 | 0 / 0 / 0 | `0x300000033` valid=1 fake=0 | `0x1ff` | `546,513,481,425,417,384,359,320,294` | cell `294,320,359`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-9.16053,-9.11626) pull=6.324; dR=(-0.155253,-0.205984) cm; dZ=(1.4222,1.8778) cm | #5 current prefix; drops 3 trailing; label `0x300000033`; refs `546,513,481,425,417,384` |
| 11 | 5 | 0 / 0 / 0 | `0x400000069` valid=1 fake=0 | `0x3e0` | `270,227,202,184,164` | cell `184,202,227`; old tgl=(-1.01585,-0.996828) pull=2.718; new=(-9.18564,-9.22246) pull=5.259; dR=(-0.154828,-0.756609) cm; dZ=(1.4222,6.9778) cm | none |
| 12 | 7 | 0 / 0 / 0 | `0x40000006d` valid=1 fake=0 | `0x7f` | `192,147,138,106,83,45,22` | cell `22,45,83`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-6.2204,-6.176) pull=6.343; dR=(-0.228635,-0.304048) cm; dZ=(1.4222,1.8778) cm | #2 current prefix; drops 2 trailing; label `0x40000006d`; refs `192,147,138,106,83` |
| 13 | 6 | 0 / 0 / 0 | `0x40000002a` valid=1 fake=0 | `0x3f` | `151,135,112,81,50,20` | cell `20,50,81`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-16.2276,-16.1647) pull=8.982; dR=(-0.0876412,-0.116167) cm; dZ=(1.4222,1.8778) cm | none |
| 14 | 5 | 0 / 0 / 0 | `0x4000000c0` valid=1 fake=0 | `0x7c` | `485,454,435,392,362` | cell `362,392,435`; old tgl=(-1.01586,-0.990749) pull=3.587; new=(-14.4171,-14.3321) pull=12.14; dR=(-0.0986466,-0.165907) cm; dZ=(1.4222,2.3778) cm | none |
| 15 | 10 | 0 / 0 / 0 | `0x38000006e` valid=1 fake=0 | `0x3ff` | `251,235,209,191,144,132,103,78,42,17` | cell `42,78,103`; old tgl=(-0.988317,-1.01586) pull=3.934; new=(-10.4031,-10.2846) pull=16.92; dR=(-0.180504,-0.138284) cm; dZ=(1.8778,1.4222) cm | none |
| 16 | 10 | 0 / 0 / 0 | `0x30000003b` valid=1 fake=0 | `0x3ff` | `556,536,502,472,421,402,378,344,315,285` | cell `285,315,344`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-12.5691,-12.4069) pull=23.17; dR=(-0.11315,-0.151351) cm; dZ=(1.4222,1.8778) cm | none |
| 18 | 5 | 0 / 0 / 0 | `0x380000047` valid=1 fake=0 | `0x3e` | `446,433,376,341,313` | cell `313,341,376`; old tgl=(-0.988317,-1.01586) pull=3.934; new=(-13.0505,-12.9447) pull=15.11; dR=(-0.143888,-0.109867) cm; dZ=(1.8778,1.4222) cm | none |
| 19 | 9 | 0 / 0 / 0 | `0x380000044` valid=1 fake=0 | `0x1ff` | `534,503,471,451,405,377,347,314,287` | cell `314,347,377`; old tgl=(-0.988317,-1.01586) pull=3.934; new=(-7.66809,-7.72528) pull=8.17; dR=(-0.244885,-0.184097) cm; dZ=(1.8778,1.4222) cm | none |
| 20 | 6 | 0 / 0 / 0 | `0x20000005e` valid=1 fake=0 | `0x3f` | `148,131,109,77,47,16` | cell `16,47,77`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-10.059,-10.1189) pull=8.564; dR=(-0.141387,-0.185573) cm; dZ=(1.4222,1.8778) cm | none |
| 21 | 8 | 0 / 0 / 0 | `0x280000109` valid=1 fake=0 | `0x1fe` | `545,512,480,424,416,383,358,319` | cell `358,383,416`; old tgl=(-1.01586,-0.990749) pull=3.587; new=(-8.7512,-8.79224) pull=5.863; dR=(-0.162515,-0.270443) cm; dZ=(1.4222,2.3778) cm | none |
| 23 | 9 | 0 / 0 / 0 | `0x4000000c2` valid=1 fake=0 | `0x1ff` | `548,516,484,430,415,388,357,325,296` | cell `296,325,357`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-7.5208,-7.56448) pull=6.241; dR=(-0.189103,-0.248239) cm; dZ=(1.4222,1.8778) cm | none |
| 24 | 5 | 0 / 0 / 0 | `0x3000000dd` valid=1 fake=0 | `0x3e0` | `565,544,511,479,423` | cell `423,479,511`; old tgl=(-0.998403,-1.01585) pull=2.493; new=(-8.69499,-8.75516) pull=8.596; dR=(-1.59607,-0.162441) cm; dZ=(13.8778,1.4222) cm | none |
| 25 | 7 | 0 / 0 / 0 | `0x280000126` valid=1 fake=0 | `0x7f` | `174,153,139,120,91,60,28` | cell `28,60,91`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-5.52796,-5.57662) pull=6.952; dR=(-0.257275,-0.336727) cm; dZ=(1.4222,1.8778) cm | none |
| 26 | 7 | 0 / 0 / 0 | `0x28000001a` valid=1 fake=0 | `0x7f` | `486,453,434,391,361,329,298` | cell `298,329,361`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-11.4903,-11.6161) pull=17.97; dR=(-0.123774,-0.161655) cm; dZ=(1.4222,1.8778) cm | none |
| 27 | 5 | 0 / 0 / 0 | `0x3000000bd` valid=1 fake=0 | `0x3e0` | `247,229,203,186,142` | cell `142,186,203`; old tgl=(-0.998403,-1.01585) pull=2.493; new=(-11.841,-11.6362) pull=29.26; dR=(-1.17202,-0.122222) cm; dZ=(13.8778,1.4222) cm | none |
| 28 | 5 | 0 / 0 / 0 | `0x28000019a` valid=1 fake=0 | `0xf8` | `210,190,146,130,108` | cell `130,146,190`; old tgl=(-1.01586,-0.998403) pull=2.494; new=(-5.84485,-5.92159) pull=10.96; dR=(-0.243325,-2.34359) cm; dZ=(1.4222,13.8778) cm | none |
| 29 | 9 | 0 / 0 / 0 | `0x38000001b` valid=1 fake=0 | `0x1ff` | `237,212,193,150,134,111,80,49,19` | cell `19,49,80`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-14.2303,-14.3199) pull=12.8; dR=(-0.099942,-0.131133) cm; dZ=(1.4222,1.8778) cm | none |
| 30 | 7 | 0 / 0 / 0 | `0x280000091` valid=1 fake=0 | `0x7f` | `474,422,404,380,346,317,286` | cell `286,317,346`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-6.76134,-6.80586) pull=6.359; dR=(-0.210343,-0.275909) cm; dZ=(1.4222,1.8778) cm | #4 current internal subsequence; drops 2; label `0x280000091`; refs `422,404,380,346,317` |
| 31 | 5 | 0 / 0 / 0 | `0x80000004000010b3` valid=1 fake=1 | `0x1f` | `137,114,84,53,26` | cell `26,53,84`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-13.4786,-12.2769) pull=171.7; dR=(-0.105516,-0.152954) cm; dZ=(1.4222,1.8778) cm | none |
| 33 | 6 | 0 / 0 / 0 | `0x3000000f6` valid=1 fake=0 | `0x1f8` | `582,521,487,455,440,393` | cell `393,440,455`; old tgl=(-0.990749,-1.01586) pull=3.587; new=(-8.13837,-8.09785) pull=5.789; dR=(-0.292171,-0.175627) cm; dZ=(2.3778,1.4222) cm | none |
| 34 | 10 | 0 / 0 / 0 | `0x280000898` valid=1 fake=0 | `0x3ff` | `566,543,514,476,427,410,386,351,323,292` | cell `351,386,410`; old tgl=(-1.01586,-0.990749) pull=3.587; new=(-7.42901,-7.4976) pull=9.799; dR=(-0.191439,-0.317142) cm; dZ=(1.4222,2.3778) cm | none |
| 35 | 6 | 1 / 198 / 0 | `0x480000029` valid=1 fake=0 | `0x3f` | `651,647,644,641,638,635` | cell `641,644,647`; old tgl=(-1.01586,-0.990749) pull=3.587; new=(-7.68329,-7.62863) pull=7.808; dR=(-0.185103,-0.311694) cm; dZ=(1.4222,2.3778) cm | none |
| 36 | 5 | 1 / 198 / 0 | `0x480000008` valid=1 fake=0 | `0x1f` | `609,602,601,595,596` | cell `595,601,602`; old tgl=(-0.988317,-1.01586) pull=3.934; new=(-5.74235,-5.78143) pull=5.583; dR=(-0.327009,-0.245995) cm; dZ=(1.8778,1.4222) cm | none |
| 37 | 9 | 1 / 198 / 0 | `0x480000007` valid=1 fake=0 | `0x1ff` | `668,660,656,650,646,643,640,637,634` | cell `634,637,640`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-12.5457,-12.4284) pull=16.75; dR=(-0.113362,-0.151089) cm; dZ=(1.4222,1.8778) cm | none |
| 38 | 7 | 1 / 198 / 0 | `0x480000054` valid=1 fake=0 | `0x7f` | `658,653,652,645,642,639,636` | cell `645,652,653`; old tgl=(-0.990749,-1.01586) pull=3.587; new=(-8.87111,-8.94742) pull=10.9; dR=(-0.268039,-0.158951) cm; dZ=(2.3778,1.4222) cm | #10 current suffix; drops 2 leading; label `0x480000054`; refs `652,645,642,639,636` |
| 39 | 5 | 2 / 396 / 0 | `0x680000042` valid=1 fake=0 | `0x1f` | `817,809,799,791,781` | cell `781,791,799`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-8.03701,-8.11191) pull=10.7; dR=(-0.176957,-0.231487) cm; dZ=(1.4222,1.8778) cm | none |
| 40 | 5 | 2 / 396 / 0 | `0x700000017` valid=1 fake=0 | `0x3e0` | `899,877,860,845,834` | cell `834,845,860`; old tgl=(-0.998403,-1.01585) pull=2.493; new=(-13.8389,-13.9259) pull=12.43; dR=(-1.00281,-0.102126) cm; dZ=(13.8778,1.4222) cm | none |
| 42 | 7 | 2 / 396 / 0 | `0x680000075` valid=1 fake=0 | `0x3f8` | `772,763,752,742,723,718,708` | cell `752,763,772`; old tgl=(-0.996828,-1.01586) pull=2.719; new=(-8.43791,-8.38306) pull=7.835; dR=(-0.826959,-0.169652) cm; dZ=(6.9778,1.4222) cm | #12 current suffix; drops 1 leading; label `0x680000075`; refs `763,752,742,723,718,708` |
| 43 | 7 | 2 / 396 / 0 | `0x680000065` valid=1 fake=0 | `0x7f` | `853,824,819,810,801,792,783` | cell `783,792,801`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-7.24813,-7.30761) pull=8.497; dR=(-0.196217,-0.256965) cm; dZ=(1.4222,1.8778) cm | none |
| 45 | 7 | 2 / 396 / 0 | `0x5800001de` valid=1 fake=0 | `0x3f8` | `898,878,861,846,835,816,808` | cell `808,816,835`; old tgl=(-0.990749,-1.01586) pull=3.587; new=(-14.8768,-14.8308) pull=6.575; dR=(-0.159832,-0.0958948) cm; dZ=(2.3778,1.4222) cm | none |
| 46 | 7 | 2 / 396 / 0 | `0x6800000b5` valid=1 fake=0 | `0x3f8` | `766,757,735,736,719,713,702` | cell `702,713,719`; old tgl=(-0.990749,-1.01586) pull=3.587; new=(-7.51979,-7.3862) pull=19.09; dR=(-0.316206,-0.192548) cm; dZ=(2.3778,1.4222) cm | #14 current prefix; drops 2 trailing; label `0x6800000b5`; refs `766,757,735,736,719` |
| 47 | 6 | 2 / 396 / 0 | `0x5800001db` valid=1 fake=0 | `0x3f0` | `885,876,858,844,823,829` | cell `823,844,858`; old tgl=(-0.998403,-1.01585) pull=2.493; new=(-5.94926,-5.9904) pull=5.876; dR=(-2.33269,-0.237412) cm; dZ=(13.8778,1.4222) cm | none |
| 48 | 7 | 2 / 396 / 0 | `0x6800000b0` valid=1 fake=0 | `0x1fc` | `764,753,745,731,728,709,700` | cell `709,728,731`; old tgl=(-0.990749,-1.01586) pull=3.587; new=(-7.77271,-7.82276) pull=7.149; dR=(-0.305916,-0.181803) cm; dZ=(2.3778,1.4222) cm | none |
| 49 | 9 | 2 / 396 / 0 | `0x68000007a` valid=1 fake=0 | `0x1ff` | `760,749,739,729,715,705,697,688,682` | cell `682,688,697`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-15.0283,-15.2191) pull=27.25; dR=(-0.094635,-0.123385) cm; dZ=(1.4222,1.8778) cm | none |
| 50 | 10 | 2 / 396 / 0 | `0x6800000b1` valid=1 fake=0 | `0x3ff` | `771,762,751,741,722,717,707,699,691,683` | cell `691,699,707`; old tgl=(-0.988317,-1.01586) pull=3.934; new=(-8.06107,-7.99539) pull=9.382; dR=(-0.232947,-0.177877) cm; dZ=(1.8778,1.4222) cm | none |
| 51 | 7 | 2 / 396 / 0 | `0x58000003c` valid=1 fake=0 | `0x3f8` | `890,883,869,854,825,830,811` | cell `811,830,825`; old tgl=(-0.990749,-1.01586) pull=3.587; new=(-11.9138,-11.9958) pull=11.72; dR=(-0.199584,-0.118558) cm; dZ=(2.3778,1.4222) cm | none |
| 52 | 6 | 2 / 396 / 0 | `0x680000385` valid=1 fake=0 | `0x3f0` | `874,875,843,840,821,815` | cell `815,821,840`; old tgl=(-1.01586,-0.998403) pull=2.494; new=(-6.5429,-6.85386) pull=44.42; dR=(-0.217365,-2.02481) cm; dZ=(1.4222,13.8778) cm | none |
| 54 | 5 | 3 / 594 / 0 | `0x780000043` valid=1 fake=0 | `0x1f` | `938,929,926,921,915` | cell `915,921,926`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-13.5434,-13.7517) pull=29.77; dR=(-0.105011,-0.13655) cm; dZ=(1.4222,1.8778) cm | none |
| 55 | 10 | 3 / 594 / 0 | `0x780000042` valid=1 fake=0 | `0x3ff` | `959,957,950,946,936,933,930,928,923,917` | cell `917,923,928`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-15.4292,-15.6599) pull=32.95; dR=(-0.0921762,-0.119912) cm; dZ=(1.4222,1.8778) cm | #16 current internal subsequence; drops 5; label `0x780000042`; refs `957,950,946,936,933` |
| 57 | 5 | 4 / 792 / 0 | `0x88000003b` valid=1 fake=0 | `0x3e` | `1152,1140,1122,1109,1092` | cell `1092,1109,1122`; old tgl=(-0.988317,-1.01586) pull=3.934; new=(-6.96801,-6.91421) pull=7.686; dR=(-0.269489,-0.205692) cm; dZ=(1.8778,1.4222) cm | none |
| 60 | 5 | 4 / 792 / 0 | `0x88000009d` valid=1 fake=0 | `0x3e0` | `1072,1067,1048,1039,1026` | cell `1026,1039,1048`; old tgl=(-0.998403,-1.01585) pull=2.493; new=(-9.44556,-9.3378) pull=15.39; dR=(-1.46924,-0.152305) cm; dZ=(13.8778,1.4222) cm | none |
| 61 | 6 | 4 / 792 / 0 | `0x900000804` valid=1 fake=0 | `0x3f` | `1154,1144,1123,1112,1095,1082` | cell `1082,1095,1112`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-1.23808,-1.19821) pull=5.696; dR=(-1.14872,-1.56718) cm; dZ=(1.4222,1.8778) cm | #22 current prefix; drops 1 trailing; label `0x900000804`; refs `1154,1144,1123,1112,1095` |
| 63 | 6 | 4 / 792 / 0 | `0x900000063` valid=1 fake=0 | `0x3f0` | `1216,1206,1186,1178,1163,1145` | cell `1145,1163,1178`; old tgl=(-1.01586,-0.998403) pull=2.494; new=(-10.1985,-10.2723) pull=10.54; dR=(-0.139452,-1.351) cm; dZ=(1.4222,13.8778) cm | none |
| 64 | 10 | 4 / 792 / 0 | `0x900000016` valid=1 fake=0 | `0x3ff` | `1057,1053,1045,1035,1014,1008,1000,990,980,970` | cell `980,990,1000`; old tgl=(-0.988317,-1.01586) pull=3.934; new=(-10.4302,-10.481) pull=7.263; dR=(-0.180036,-0.135693) cm; dZ=(1.8778,1.4222) cm | none |
| 66 | 9 | 4 / 792 / 0 | `0x88000003c` valid=1 fake=0 | `0x1ff` | `1209,1189,1181,1157,1149,1132,1118,1102,1087` | cell `1087,1102,1118`; old tgl=(-1.01586,-0.988317) pull=3.934; new=(-8.53858,-8.78866) pull=35.73; dR=(-0.166562,-0.213662) cm; dZ=(1.4222,1.8778) cm | none |
| 67 | 5 | 4 / 792 / 0 | `0x980000098` valid=1 fake=0 | `0x1f` | `1009,1002,991,981,971` | cell `981,991,1002`; old tgl=(-0.988317,-1.01586) pull=3.934; new=(-4.52297,-4.41345) pull=15.65; dR=(-0.41517,-0.322242) cm; dZ=(1.8778,1.4222) cm | none |
| 68 | 9 | 4 / 792 / 0 | `0x8000000980000059` valid=1 fake=1 | `0x3fe` | `1220,1211,1192,1182,1158,1150,1133,1121,1103` | cell `1121,1133,1150`; old tgl=(-1.01586,-0.990749) pull=3.587; new=(-6.83552,-6.72773) pull=15.4; dR=(-0.20806,-0.353433) cm; dZ=(1.4222,2.3778) cm | none |

## Complete standalone additions and substitutions

All nine current additions are shorter related candidates rather than unrelated new tracks. Together with the 55 removals, these produce the net loss of 46 tracks.

| Current # | Hits | ROF / BC / orbit | Label | Seed | Ordered references | Best old relation |
|---:|---:|---|---|---|---|---|
| 1 | 5 | 0 / 0 / 0 | `0x28000001f` valid=1 fake=0 | `0x3e` | `149,133,110,79,48` | old #1 current prefix; drops 1 trailing; label `0x28000001f`; refs `149,133,110,79,48,18` |
| 2 | 5 | 0 / 0 / 0 | `0x40000006d` valid=1 fake=0 | `0x7c` | `192,147,138,106,83` | old #12 current prefix; drops 2 trailing; label `0x40000006d`; refs `192,147,138,106,83,45,22` |
| 4 | 5 | 0 / 0 / 0 | `0x280000091` valid=1 fake=0 | `0x3e` | `422,404,380,346,317` | old #30 current internal subsequence; drops 2; label `0x280000091`; refs `474,422,404,380,346,317,286` |
| 5 | 6 | 0 / 0 / 0 | `0x300000033` valid=1 fake=0 | `0x1f8` | `546,513,481,425,417,384` | old #9 current prefix; drops 3 trailing; label `0x300000033`; refs `546,513,481,425,417,384,359,320,294` |
| 10 | 5 | 1 / 198 / 0 | `0x480000054` valid=1 fake=0 | `0x1f` | `652,645,642,639,636` | old #38 current suffix; drops 2 leading; label `0x480000054`; refs `658,653,652,645,642,639,636` |
| 12 | 6 | 2 / 396 / 0 | `0x680000075` valid=1 fake=0 | `0x1f8` | `763,752,742,723,718,708` | old #42 current suffix; drops 1 leading; label `0x680000075`; refs `772,763,752,742,723,718,708` |
| 14 | 5 | 2 / 396 / 0 | `0x6800000b5` valid=1 fake=0 | `0x3e0` | `766,757,735,736,719` | old #46 current prefix; drops 2 trailing; label `0x6800000b5`; refs `766,757,735,736,719,713,702` |
| 16 | 5 | 3 / 594 / 0 | `0x780000042` valid=1 fake=0 | `0x1f0` | `957,950,946,936,933` | old #55 current internal subsequence; drops 5; label `0x780000042`; refs `959,957,950,946,936,933,930,928,923,917` |
| 22 | 5 | 4 / 792 / 0 | `0x900000804` valid=1 fake=0 | `0x3e` | `1154,1144,1123,1112,1095` | old #61 current prefix; drops 1 trailing; label `0x900000804`; refs `1154,1144,1123,1112,1095,1082` |

## Complete current combined-only inventory

These are all 27 current combined sequences not present in current standalone.

| Combined # | Hits | ROF / BC / orbit | Label | Seed | Ordered references | Components present in standalone | Explanation |
|---:|---:|---|---|---|---|---|---|
| 0 | 4 | 0 / 0 / 0 | `0x4000000c2` valid=1 fake=0 | `0x78` | `484,430,415,388` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 1 | 4 | 0 / 0 / 0 | `0x38000001b` valid=1 fake=0 | `0xf0` | `212,193,150,134` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 2 | 4 | 0 / 0 / 0 | `0x3000000dd` valid=1 fake=0 | `0x3c0` | `565,544,511,479` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 3 | 4 | 0 / 0 / 0 | `0x280000898` valid=1 fake=0 | `0xf0` | `514,476,427,410` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 4 | 4 | 0 / 0 / 0 | `0x380000044` valid=1 fake=0 | `0x1e0` | `534,503,471,451` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 5 | 4 | 0 / 0 / 0 | `0x300000052` valid=1 fake=0 | `0x3c` | `141,121,92,65` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 8 | 4 | 0 / 0 / 0 | `0x4000000bd` valid=1 fake=0 | `0xf` | `394,364,333,302` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 9 | 4 | 0 / 0 / 0 | `0x280000898` valid=1 fake=0 | `0xf` | `386,351,323,292` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 10 | 4 | 0 / 0 / 0 | `0x280000126` valid=1 fake=0 | `0x1e` | `139,120,91,60` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 11 | 4 | 0 / 0 / 0 | `0x30000001d` valid=1 fake=0 | `0xf0` | `517,482,431,414` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 15 | 4 | 0 / 0 / 0 | `0x38000006e` valid=1 fake=0 | `0x1e0` | `235,209,191,144` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 16 | 4 | 0 / 0 / 0 | `0x20000005e` valid=1 fake=0 | `0x1e` | `131,109,77,47` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 19 | 4 | 0 / 0 / 0 | `0x28000010a` valid=1 fake=0 | `0x3c0` | `244,223,177,181` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 23 | 4 | 1 / 198 / 0 | `0x480000029` valid=1 fake=0 | `0xf` | `644,641,638,635` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 25 | 4 | 2 / 396 / 0 | `0x58000004e` valid=1 fake=0 | `0xf` | `710,701,693,685` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 26 | 4 | 2 / 396 / 0 | `0x68000007f` valid=1 fake=0 | `0x1e` | `724,704,695,687` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 27 | 4 | 2 / 396 / 0 | `0x6800000b0` valid=1 fake=0 | `0xf0` | `753,745,731,728` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 28 | 4 | 2 / 396 / 0 | `0x680000042` valid=1 fake=0 | `0x1e` | `817,809,799,791` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 29 | 4 | 2 / 396 / 0 | `0x6800000b1` valid=1 fake=0 | `0x1e0` | `762,751,741,722` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 30 | 4 | 2 / 396 / 0 | `0x680000077` valid=1 fake=0 | `0x1e` | `716,706,698,690` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 31 | 4 | 2 / 396 / 0 | `0x68000003d` valid=1 fake=0 | `0x3c0` | `889,881,868,851` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 34 | 4 | 2 / 396 / 0 | `0x5800001db` valid=1 fake=0 | `0x3c0` | `885,876,858,844` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 38 | 4 | 3 / 594 / 0 | `0x7800000e7` valid=1 fake=0 | `0x3c0` | `960,958,953,947` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 40 | 4 | 4 / 792 / 0 | `0x900000016` valid=1 fake=0 | `0x1e0` | `1053,1045,1035,1014` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 42 | 4 | 4 / 792 / 0 | `0x900000041` valid=1 fake=0 | `0x78` | `1175,1155,1141,1124` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 43 | 4 | 4 / 792 / 0 | `0x900000063` valid=1 fake=0 | `0x1e0` | `1206,1186,1178,1163` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |
| 46 | 4 | 4 / 792 / 0 | `0x9000000b0` valid=1 fake=0 | `0x3c0` | `1201,1197,1174,1166` | pairs 3/3; cells 2/2 | common MinTrackLength=4 versus standalone MFT 5 |

## Validation and repository hygiene

- Pinned base/current standalone and combined replays completed; counts and hashes match the reviewed products.
- Persisted association checks report standalone `retained=14 removed=55 added=9`, combined `retained=24 removed=76 added=26`, and current composition `retained=23 removed=0 added=27`.
- All three association checks report `ROF_METADATA differences=0` over 2,304 records.
- The fixture checksum manifest passed 43/43 before diagnostics and 43/43 afterward.
- The instrumentation-free current library and five focused tests were rebuilt. `testTrackletFinding`, `testComputeLayerTrackletsOrchestration`, `testComputeLayerCellsOrchestration`, `testTrackingKernelParameters`, and `testTraversalOperationBindingContainment` pass 5/5.
- Temporary production instrumentation was removed. No production source, CMake, workflow default, graph/hole policy, or physics parameter was changed.
- ITS replay was not rerun and no new ITS-preservation claim is made.

## Verdict and boundary for follow-up

**Verdict: requires targeted correction.**

The semantic correction (`tanLambda = dz/dr`) and same-radius guard are sound. The defect is the unproven reuse of the old raw `0.007` cell scale after the disk quantity changed meaning. It is responsible for the first rejection of every absent historical standalone sequence and for the large loss of non-fake matched population. This is neither a combined-policy effect nor a neighbour/road effect.

Do not restore the nominal-Z pseudo-slope and do not add an MFT-specific compatibility value. Before neighbour/road harmonization, define and validate a coordinate-neutral, uncertainty-aware compatibility measure for consecutive physical slopes through the existing measurement/state-operation boundary, then repeat strict ITS preservation and MFT characterization.

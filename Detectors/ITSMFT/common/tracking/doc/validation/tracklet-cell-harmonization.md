# Tracklet and cell harmonization validation

- Status: implementation candidate; unified MFT physics sign-off pending
- Date: 2026-08-09
- Design: [tracklet and cell harmonization](../design/0015-tracklet-cell-harmonization.md)
- Predecessor: [coordinate-cut retirement](coordinate-cut-retirement.md)
- Input commit: `be70d57e29`
- Production commit: `bbf9b27579`
- Test commit: `0bd80c179b`
- Pinned O2 package: `daily-20260717-0700-local1`
- Build: `/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`
- Fresh artifacts: `/private/tmp/o2-itsmft-tracklet-cell-harmonization-20260809-VV6rxF`

## Implemented contract

`Tracklet::tanLambda` now means `delta z / delta r` for both surface
families. Cylinder arithmetic is unchanged. The disk leaf calculates global
transverse radii from the accepted measurements and rejects
`abs(delta r) <= 1e-6 cm`, where the quotient is undefined.

Tracklet and cell construction each traverse the binding's single ordered
list once. Generic `TrackerTraits` orchestration no longer partitions those
lists or selects typed cylinder/disk implementations. Projection, candidate
acceptance, row wrapping, and cell-seed construction select coordinate math
only inside their leaf operations. There is one scalar
`TrackingKernelParameters` record rather than one record per surface kind.

The flat pair graph, hole policy, scalar values, `(6,7)` disconnection, ROF
ownership/timing fix, and neighbour/road traversal are unchanged.

## Build and tests

The directly affected operation, kernel-layout, orchestration, and source
containment tests passed 5/5. They exercise the physical disk `delta z /
delta r` value and its zero-`delta r` leaf guard, preserve the cylinder
arithmetic, prove one global tracklet loop and one global cell loop, and
confine surface selection to the leaf source.

The complete required serial suite passed:

```text
ctest -L itsmft --output-on-failure -j1
100% tests passed, 0 tests failed out of 93
Total Test time (real) = 119.39 sec
```

The ITS/MFT test targets, replay executables, and required DPL support
libraries were built. `CMAKE_CUDA_COMPILER` is `NOTFOUND`; host and
device-guard compilation pass, but this campaign makes no real CUDA/HIP
compiler claim.

## Preflight and fixture integrity

Strict preflight resolved the exact pinned package, found a valid AliEn token
from `2026-08-07 05:11:10 UTC` through `2026-09-05 09:33:43 UTC`, confirmed
that `ALICEO2_CCDB_NOTOKENCHECK` was unset, and found every required Geant4
dataset. Replays used CCDB condition cap `1784207296000` and no token bypass.

The durable fixture's 43-entry checksum manifest passed before replay and
again after every primary, repeat, and comparison run:

```text
/Users/mpuccio/alice/run3/O2-validation-data/itsmft/fixtures/
  pp-20ev-run303000-seed20260716-daily20260717/checksums.sha256
43/43 checksums passed before replay
43/43 checksums passed after all replays
```

## Replay result

The preceding coordinate-cut candidate is the characterization reference,
not an MFT acceptance baseline.

| Product | Preceding tracks / hash | Current tracks / hash | Matched / reconstructable | Efficiency | Fake tracks / rate | Clone tracks / rate |
|---|---:|---:|---:|---:|---:|---:|
| ITS standalone | 212 / `46913a67a7e2fe7462e29df0db264fa8` | 212 / `46913a67a7e2fe7462e29df0db264fa8` | 135 / 142 | 0.950704 | 16 / 0.0754717 | 0 / 0 |
| ITS combined | 212 / `46913a67a7e2fe7462e29df0db264fa8` | 212 / `46913a67a7e2fe7462e29df0db264fa8` | 135 / 142 | 0.950704 | 16 / 0.0754717 | 0 / 0 |
| MFT standalone | 69 / `f6dee3f7a5f7def6b55900dbac734ef0` | 23 / `c65075052aa08b7be8a2283cdedc79ea` | 19 / 109 | 0.174312 | 0 / 0 | 0 / 0 |
| MFT combined | 100 / `4dacc3058740cce83656a3bcb71def95` | 50 / `45e46e0475a5d6479054319e8daba176` | 41 / 109 | 0.376147 | 0 / 0 | 1 / 0.00917431 |

ITS standalone and combined products match their preceding products and each
other field-for-field, including ordered references, labels, ROFs, states,
and covariances. This is strict ITS preservation.

Fresh repeat MFT standalone and combined replays match their respective
track counts, hashes, and every compared persisted field. The MFT comparison
excludes the known unwritten `mInvQPtSeed` value, as in earlier campaigns.

## MFT association characterization

Exact identity here means equal ROF and equal ordered cluster-reference list.
The artifact logs record every added and removed track with output index, ROF,
seed pattern, raw MC label and fake bit, chi2, and ordered references.

| Comparison | Retained | Removed | Added | ROF metadata differences |
|---|---:|---:|---:|---:|
| preceding standalone -> current standalone | 14 | 55 | 9 | 0 / 2304 |
| preceding combined -> current combined | 24 | 76 | 26 | 0 / 2304 |
| current standalone -> current combined | 23 | 0 | 27 | 0 / 2304 |

The standalone additions contain seven five-hit and two six-hit tracks. Its
removals contain 16 five-hit, 13 six-hit, 11 seven-hit, one eight-hit, eight
nine-hit, and six ten-hit tracks. Nine removed tracks have a current track
with the same raw label and ROF sharing references; the log identifies each
best suffix or subset explicitly.

The combined additions contain 17 four-hit, seven five-hit, and two six-hit
tracks. Its removals contain 19 four-hit, 17 five-hit, 14 six-hit, 11
seven-hit, one eight-hit, eight nine-hit, and six ten-hit tracks. Twenty-five
removed tracks have a same-label, same-ROF current overlap. Current combined
contains all 23 standalone associations exactly and adds 27 four-hit tracks;
this is structural evidence, not a requirement that the two workflows agree.

All changed associations are in ROFs 0 through 4. Across all 2304 ROFs, BC,
orbit, RO-frame, and flag metadata are unchanged, so no timing or ownership
change is implicated. The complete inventories are
`mft-standalone-associations.log`, `mft-combined-associations.log`, and
`mft-current-structure.log` in the fresh artifact directory.

## Attribution and disposition

The removed disk expression divided the measured `delta z` by the nominal
disk `delta z`. For adjacent disks it was therefore close to a signed
constant, causing differences between consecutive tracklet slopes to be
artificially near zero. The existing common cell selection

```text
abs(first.tanLambda - second.tanLambda) / 0.007 < 5
```

now compares actual `delta z / delta r` values. Non-collinear MFT triples
that passed the legacy near-constant ratio consequently stop forming cells;
some surviving sequences are shorter subsets of preceding tracks. This
directly explains the reference substitutions and population reduction.
The new same-radius rejection is only the mathematical validity guard for
the quotient and is tested at the disk leaf.

No MFT-specific tuning or compatibility branch was introduced to restore
the preceding counts. The deterministic MFT result remains a candidate for
later unified physics sign-off.

## Follow-up inventory

- Move disk projection inputs out of the legacy MFT reference-Z table, then
  remove that table when descriptors own the required state.
- Harmonize neighbour and road traversal after state-family attachment and
  compatibility are wholly leaf-owned; then remove their cylinder/disk
  wrappers and partitions.
- Remove `MFTFwdTrackHelpers.{h,cxx}` with the remaining detector-specific
  MFT common-tracker facilities.

Those items deliberately remain outside this slice.

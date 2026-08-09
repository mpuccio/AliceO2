# Tracklet and cell harmonization

- Status: accepted for implementation; physics validation pending
- Date: 2026-08-09
- Predecessor: [coordinate-cut retirement](0014-coordinate-cut-retirement.md)
- Pinned O2 package: `daily-20260717-0700-local1`

## Decision

`Tracklet::tanLambda` has one meaning for every surface family:

```text
tanLambda = delta z / delta r
```

Here `r` is the global transverse radius of the measurement. Cylinder
tracklets retain their existing radius arithmetic so ITS output remains
bitwise stable. Disk tracklets calculate the same quantity from the two
accepted global measurements instead of dividing the Z difference by the
nominal disk separation.

The disk calculation is undefined when the two measurements have the same
transverse radius. The `abs(delta r) > 1e-6 cm` check is therefore a
coordinate-validity guard. It belongs to disk candidate acceptance, behind
the existing tracklet-operation boundary, and is not a tracking parameter or
a generic `SurfaceKind::Disk` branch.

## Orchestration boundary

Tracklet enumeration and cell enumeration are common CA operations. They
must traverse the binding's ordered transition and cell lists once, without
pre-partitioning those lists merely to select cylinder or disk code. A narrow
tracklet leaf owns projection, candidate acceptance, and index-row wrapping;
a narrow cell leaf owns seed construction. Surface kind may select the
coordinate mathematics inside those leaves only.

The shared cell slope comparison remains

```text
abs(first.tanLambda - second.tanLambda) / sigma < nSigma
```

With the unified `tanLambda` contract it compares one physical quantity for
both families. The former disk value was approximately a signed constant
because both numerator and denominator represented the same nominal Z
separation. Its resulting MFT acceptance was compatibility behavior, not a
valid common slope selection. MFT changes from correcting that value are to
be characterized, not tuned away.

## Scope boundary

This slice does not change the flat pair-list graph, hole policy, scalar
configuration policy, the `(6,7)` disconnection, ROF timing ownership, or
road topology. Neighbour walking and refit remain separate follow-up work
except where a cell leaf signature can be narrowed without changing their
orchestration.

## Later-deletion inventory

- The legacy MFT reference-Z table remains necessary for the current disk
  projection and uncertainty calculation, but no longer defines tracklet
  slope. Remove it when those projection inputs become descriptor-owned.
- The cylinder/disk neighbour and road wrappers remain after this slice.
  Once state-family attachment and compatibility are selected wholly inside
  their leaf operations, remove those orchestration partitions and wrappers.
- `MFTFwdTrackHelpers.{h,cxx}` remains in the inventory established by design
  note 0014.

## Validation contract

Focused tests must exercise the physical disk slope and its zero-`delta r`
guard directly, preserve cylinder arithmetic, prove that both surface kinds
use one tracklet loop and one cell loop, and confine their surface selection
to leaf translation units. ITS replay is a strict preservation gate. MFT
standalone and combined replay are characterization evidence, including
ordered track/reference differences and efficiency, fake, and clone effects
where available.

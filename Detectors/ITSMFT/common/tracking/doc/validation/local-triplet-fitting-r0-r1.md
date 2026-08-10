# Local triplet fitting R0/R1 validation

- Date: 2026-08-10
- Package: `daily-20260717-0700-local1`
- Worktree build: `/Users/mpuccio/alice/run3/O2-worktree-builds/triplet-tracking-rnd-scratch`
- Design: [detector-neutral local triplet fitting](../design/0018-local-triplet-fitting-r0-r1.md)

## Scope and output contract

This validation covers the standalone R1 local fitter. There is deliberately
no call from `TrackerTraits`, no configuration projection, and no workflow
switch. Consequently, the cell population and published ITS/MFT tracks are
unchanged by construction. Physics replays and fixture checksum campaigns
were not repeated: they cannot exercise an uncalled primitive and would add
no evidence beyond the source guard and complete test label.

The user stash remained untouched.

## Numerical coverage

The focused test covers ten cases:

1. cylinder `(u,v)` covariance lifted to the full rank-two global covariance;
2. disk global `(x,y)` covariance with exact descriptor `z`;
3. an exact non-axis-aligned helix, including curvature and transverse
   momentum;
4. simultaneous global rotation around Z, translation, and covariance
   rotation;
5. non-zero `xy`, `xz`, and `yz` covariance contributions;
6. the removable straight/zero-bending limit;
7. a non-helical kink and the expected quality reduction after adding MS;
8. non-PSD, non-finite, repeated-point, and invalid-process inputs failing
   transactionally;
9. source guards proving no family/detector/source dispatch in the common fit
   and no `TrackerTraits` call site;
10. deterministic host-cost characterization.

All cases passed.

## Performance characterization

The pinned RelWithDebInfo binary ran 20,000 identical fits after construction
of the observations. The measured cost on this host was:

```text
245.4375 ns/fit
```

This is an R1 microbenchmark, not a production throughput claim. It excludes
surface-observation construction, material iteration, candidate traversal,
and cache/occupancy effects. R2 must measure those effects on actual cell
candidates before choosing whether the fit replaces the two cheap gates or
runs only for their survivors.

## Test campaign

```text
focused TripletFitting + header inventory: 2/2 passed
complete serial ctest -L itsmft -j1:         94/94 passed
```

The first full label run correctly failed the exact public-header inventory
at `34 != 33`. The inventory was updated to name `TripletFitting.h`; its
focused rerun and the complete rerun then passed.

## R2 outlook

The standalone primitive is ready for candidate-level characterization but
not yet for a physics cut. The next slice must bind middle-surface material,
derive MS at a fitted/reference momentum through the existing common material
kernel, iterate the local fit explicitly, and collect candidate cut-flow,
pulls, and timing without changing publication. Only that evidence can decide
the production placement and common `NSigmaCut` contract.

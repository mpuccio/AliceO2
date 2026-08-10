# Detector-neutral local triplet fitting (R0/R1)

- Status: prototype design
- Date: 2026-08-10
- Predecessor: [transverse tracklet-direction compatibility](0017-transverse-tracklet-direction-compatibility.md)
- Reference: A. Schoening, [A General Track Fit based on Triplets](https://arxiv.org/abs/2406.05240)

## Scope

This slice implements the R0 data-contract audit and the R1 standalone local
triplet primitive described by the deferred triplet-tracking R&D note. It does
not wire the fit into cell filtering, replace `CellSeed`, alter the neighbour
graph, or change workflow output. Production integration requires a separate
R2 characterization of candidate rejection, cost, and momentum pulls.

## Neutral observation contract

The common fit receives three ordered observations. Each contains a global
position `(x,y,z)` and the full symmetric global covariance

```text
C = ( Cxx Cxy Cxz )
    ( Cxy Cyy Cyz )
    ( Cxz Cyz Czz ).
```

Zero eigenvalues are valid: a coordinate fixed by a surface is expressed by
zero variance. Inputs must be finite and positive semidefinite. The common
fit contains no `SurfaceKind`, detector ID, source ID, or geometry lookup.

The existing normalized measurement is sufficient:

- Cylinder: `(q,u,v)` uses measured `(u,v)` with `q` exact. With
  `e_u=(-sin(alpha),cos(alpha),0)` and `e_v=(0,0,1)`, the leaf constructs
  `C_global = [e_u e_v] C_uv [e_u e_v]^T` and obtains the point from the same
  frame transform.
- Disk: decoder covariance is global `(x,y)`. The leaf copies its `2x2`
  block, uses descriptor `z`, and sets every covariance component involving
  `z` to zero.

No decoder, static plan, or `SurfaceMeasurement` extension is required.

## Local-fit contract

The primitive implements the uniform-solenoid triplet parameters of Section
6.1 and the local fit of Section 4 of the reference. From the three points it
constructs the circle-reference parameters

```text
(PhiTilde, ThetaTilde, rhoPhi, rhoTheta)
```

and the hit gradients

```text
hPhi_i   = grad_i(PhiTilde)   + kappaRef grad_i(rhoPhi)
hTheta_i = grad_i(ThetaTilde) + kappaRef grad_i(rhoTheta),
kappaRef = -PhiTilde / rhoPhi.
```

The gradients are evaluated by forward automatic differentiation over all
nine hit coordinates. This is analytic propagation through the implemented
geometry, not a covariance-scaled finite difference and not a configurable
resolution approximation.

Contracting every gradient with the complete per-hit covariance gives
`GammaThetaTheta`, `GammaThetaPhi`, and `GammaPhiPhi`. The neutral process
input is the variance `sigmaMS^2` of the physical space angle at the middle
observation. Its spherical-coordinate projections are

```text
sigmaTheta^2 = sigmaMS^2
sigmaPhi^2   = sigmaMS^2 / sin(thetaHat)^2,
```

where `thetaHat` is the endpoint-chord polar angle. The local curvature,
variance, and one-degree-of-freedom quality then follow Equations 46--48.
The fit is conditional on the supplied process variance; momentum-dependent
material iteration belongs to R2 and is not hidden inside this primitive.

The result reports geometric 3D curvature in `1/cm`. A helper converts it to
transverse momentum with the fitted curvature and reference polar direction;
magnetic field is not an input to the geometric fit itself.

## Numerical validity

Repeated points, collinear transverse points with an undefined circle
reference, recurling/asin-domain geometries, non-finite values, non-PSD
covariances, vanishing polar projection, and non-positive fit covariance fail
transactionally. Near zero bending, algebraic series for the otherwise
removable `0/0` terms are selected by a machine-precision-scale numerical
branch. This is not a physics cut.

## R2 entry criteria

Before production filtering, the prototype must demonstrate:

- correct synthetic curvature, covariance, rotation/translation invariance,
  and failure behavior;
- stable pulls with realistic ITS and MFT measurement covariances;
- an explicit momentum-dependent material iteration using the fitted
  momentum, PID, charge, and middle-surface material;
- candidate-level cut-flow and CPU cost against the two current cheap gates;
- a decision on whether the fit replaces both gates or runs only on their
  survivors.


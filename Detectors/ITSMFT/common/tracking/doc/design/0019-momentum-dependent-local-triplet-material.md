# Momentum-dependent material for the local triplet fit (R2)

- Status: archived experiment
- Date: 2026-08-10
- Predecessor: [detector-neutral local triplet fitting](0018-local-triplet-fitting-r0-r1.md)
- Preserved implementation: `ae27032797`
- Successor: [fast triplet estimator](0020-fast-triplet-estimator.md)

## Experimental decision

R2 kept the geometric triplet fit detector-neutral and added one explicit
fixed-point layer around it. It is not called from `TrackerTraits` and does
not yet replace either cheap cell-direction gate.

The family leaf converts the middle surface to a neutral material observation:

```text
(unit surface normal, normal-incidence x/X0, normal-incidence areal density).
```

For a cylinder the normal is the local radial frame direction
`(cos(alpha), sin(alpha), 0)`. For a disk it is `(0,0,1)`. No detector or
source identity enters either leaf.

## Incidence and momentum contract

The common code derives two 3D segment directions from the same signed circle
curvature and arc geometry used by the local fit. At the middle point, the
incoming chord direction is rotated by its positive half-bend and the
outgoing chord direction by its negative half-bend. Each receives its exact
segment `tanLambda = dz / transverseArcLength`; their normalized bisector is
the reference direction `tHat`.

Normal-incidence material is converted to path material only through

```text
incidenceCosine = abs(dot(tHat, nHat))
path x/X0       = nominal x/X0 / incidenceCosine
path density    = nominal density / incidenceCosine.
```

A zero or non-finite incidence is invalid. There is no incidence floor and no
family-specific threshold.

The local fit returns 3D curvature in `1/cm`. For non-zero field and charge,

```text
p  = absCharge * abs(B2C * Bz / curvature)
pT = p * referenceSinTheta.
```

The shared scalar material kernel then evaluates the physical Highland
space-angle variance using this momentum, PID, absolute charge, and
path-integrated middle-surface material. The local fit consumes that variance
through its existing neutral process-noise input. Energy-loss output is not
applied to the fitted momentum: the kernel's pre-material momentum defines
the scattering kick, while its normal validity checks still reject a particle
that cannot traverse the stated material.

## Fixed-point solution

The required equation is

```text
sigmaMS2 = Highland(p(sigmaMS2)).
```

The first eight updates use direct substitution. Later updates use equal
under-relaxation between the current and newly evaluated variance. This
retains the same physical fixed point while stabilizing oscillating real
candidates. Convergence is relative at `4.76837158203125e-7`; the solve is
bounded at 128 material evaluations and otherwise returns `NoConvergence`.
These are numerical solver controls, not configurable physics selections.

Zero material reproduces the measurement-only fit. Exact zero bending has
infinite momentum and zero scattering without dividing by curvature. All
other invalid geometry, PID/charge/field, material-domain, and convergence
failures are explicit and transactional. Material failures distinguish the
kernel's lower momentum range, excessive scattering, and stopping in
material from other evaluation failures.

## Integration decision

The real-candidate R2 characterization demonstrates a usable ITS momentum
estimate but a broad MFT three-hit transverse-momentum residual. Its iterative
material wrapper is therefore not part of the active prototype. The exact
implementation and tests remain recoverable at the commit named above.

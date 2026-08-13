# Generic tracklet search window

`TrackletSearchWindow` is the single data-only contract between
surface-coordinate projection and common tracklet traversal:

```cpp
struct TrackletSearchWindow {
  int4 bins;
  float prediction[2];
  float variance[3]; // uu, uv, vv
};
```

It contains neither a surface-family tag nor a callback. The traversal reads
the bins directly and passes `TrackingKernelParameters::nSigmaCut` separately
to the common covariance-normalized acceptance. Coordinate conversion is
selected at the descriptor boundary; detector/source identity has no role.

For disks, the projected coordinates are `(x,y)`. The disk leaf projects to
the event-derived target mean z `(targetMinZ + targetMaxZ) / 2`. Treating the
target z interval as uniform, it adds the secant-propagated covariance
`(p(targetMaxZ)-p(targetMinZ))(...)^T / 12`, including the XY cross term.
For cylinders the coordinates are `(z,phi)`. The cylinder leaf
projects `z` to the target surface mean radius
`(targetMinR + targetMaxR) / 2`; its radial interval is modeled as uniform,
so it adds `tanLambda^2 * (targetMaxR-targetMinR)^2 / 12` to the longitudinal
prediction variance. The azimuth residual is wrapped modulo `2 pi`, and its
variance is `(transitionPhiCut / NSigmaCut)^2`, so the established transition
cut becomes a one-axis `NSigmaCut` covariance scale.

This deliberately changes barrel candidate selection from independent strict
`z`/`phi` cuts to one strict two-dimensional ellipse. It is a physics change,
not a legacy-equivalence refactor, and requires replay characterization before
acceptance.

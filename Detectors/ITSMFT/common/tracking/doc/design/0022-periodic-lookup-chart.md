# Periodic lookup chart

Every tracking surface has a finite `SurfaceChartRange` for its LUT's
non-periodic coordinate.  Cylinders use `[-zMax,zMax]`; disks use
`[rMin,rMax]`.  The static ITS ranges are the common tracking configuration's
existing `LayerZ` defaults.  The MFT ranges reproduce the existing immutable
`RMin/RMax` geometry constants; a test compares both tables.

`detail/LookupCoordinates.h` transforms a global measurement and its full XY
covariance into `(phi, transverse)`, retaining the off-diagonal term.  Phi is
always periodic.  Disk conversion rejects `r == 0` rather than introducing a
numeric floor.  Its window helper bounds the correlated ellipse
conservatively in the chart, clamps only the non-periodic coordinate to the
descriptor range, and records phi wrapping.

This is an inert diagnostic contract.  The production disk LUT and its XY
acceptance remain unchanged.  A later slice must compare candidate populations
and choose the new acceptance explicitly; this contract does not claim that
the two ellipses are equivalent.

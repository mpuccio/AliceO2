# Descriptor-driven neighbour and road traversal

## Decision

Transition preparation, neighbour attachment, and road-start scheduling use
the graph's ordered surfaces directly. A transition may join cylinder and disk
surfaces. Coordinate-specific math is selected from the target surface
descriptor at the operation leaf.

The curvature bound is transition-local: it starts from the configured
`TrackletMinPt` bound and is clamped using the minimum and maximum endpoint
radii. It is not carried between transitions or graph components.

## Current limit

The retained detector-specific refit callback still selects one legacy refit
from the first seed hit and validates one source for that refit. Consequently
mixed-family roads can be constructed, but mixed-family accepted-track refit
is not yet supported. This is an explicit later migration, not a reason for
the graph/road layer to reject a mixed edge.

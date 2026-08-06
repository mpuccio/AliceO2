// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACESTATEOPERATIONRESULT_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACESTATEOPERATIONRESULT_H_

#include <cstdint>

namespace o2::itsmft::tracking
{

enum class OperationFailureReason : uint8_t {
  SourceFamilyMismatch = 0,
  NonFiniteInput = 1,
  NonFiniteOutput = 2,
  InvalidCovariance = 3,
  UnreachableTarget = 4,
  PropagationFailure = 5,
  MaterialFailure = 6,
  PredictedChi2Failure = 7,
  UpdateFailure = 8,
  RotationFailure = 9,
  AlphaMismatch = 10,
  ReferenceCoordinateMismatch = 11,
  // forward::buildSeed's own established strict-boundary rejections
  // (the retained detail::mftFwdFitCellClusters / buildCellSeed<DiskDisk>
  // initializer treats an insufficient inner/outer z-ordering margin, or any
  // inner-middle/inner-outer separation below its 1e-6f minimum, as a hard
  // rejection before any direction estimate is computed -- not a NaN/Inf
  // artifact of the arithmetic itself). Distinct from NonFiniteInput (raw
  // inputs are finite; the *geometry* is degenerate) and from
  // NonFiniteOutput (no output was even attempted yet). barrel::buildSeed
  // never raises this: its retained oracle (o2::its::track::buildTrackSeed)
  // has no analogous early rejection, only NonFiniteOutput-detected numeric
  // fallbacks.
  SeedGeometryDegenerate = 12,
  // barrel::/forward::buildSeed(SeedAnchor, ...): the anchor argument is
  // not one of the enum's locked values (SeedAnchor::Inner/Outer).
  // Distinct from NonFiniteInput -- the raw measurement/bz/absCharge
  // inputs may be perfectly well-formed; it is the anchor selector itself
  // that is invalid -- so an out-of-range SeedAnchor is never
  // misclassified as a numeric-input problem.
  InvalidSeedAnchor = 13,
  // driveRefitLeg<Tag>: a present (non-hole) slot's SurfaceId/catalog
  // association could not be validated -- covers an invalid
  // measurement.surface, a catalog with a null surfaces pointer and a
  // nonzero nSurfaces, an out-of-range measurement.surface.value(), or a
  // resolved SurfaceDescriptor whose own id does not match
  // measurement.surface. Distinct from every reason above: no propagation,
  // material, or chi2 arithmetic has been attempted yet for this slot.
  InvalidSurfaceCatalogAssociation = 14,
  // Gate 3 Slice B (native CylinderCylinder refit driver, unwired): one leg's
  // final acceptance check failed after driveRefitLeg<Tag> itself already
  // succeeded for that leg -- reproducing the frozen ITS fitTrack's own
  // trailing `|Q2Pt| < maxQoverPt && chi2 < maxChi2NDF*(nCl*2-5)` return
  // condition (ITSMFTTracking/TrackHelpers.h), evaluated once per leg with
  // that leg's own maxQoverPt (VeryBig for the two inward-index legs, 50.f
  // for the outward-index leg) and driveRefitLeg's own per-leg
  // acceptedHitCount/chi2 outputs. Distinct from PredictedChi2Failure (a
  // per-hit rejection raised inside refitHit/driveRefitLeg itself, before
  // this whole-leg check is even reached).
  LegAcceptanceFailure = 15,
  // Gate 3 Slice B: the frozen ITS refitTrackSeed's trailing
  // `if (minPt > 0.f && track.getPt() < minPt) return false;` check
  // (TrackHelpers.h), evaluated once after the outward-index leg using
  // the active traversal count minus seed.getHitLayerMask().count(). Distinct from
  // LegAcceptanceFailure: this is a refitTrackSeed-level check keyed on the
  // seed's own attached-cluster count, not a per-leg fitTrack acceptance
  // condition.
  MinPtFailure = 16,
  // Gate 3 Slice B hardening: nativeRefitTrackCylinderCylinder does not
  // reproduce the frozen seedTrackForRefit's conditional mid-track geometric
  // reseed (ncl < reseedIfShorter && ncl > 2, re-deriving the initial
  // parametrization via buildTrackSeed/selectReseedMidLayer from raw
  // Cluster/TrackingFrameInfo, ITSMFTTracking/TrackHelpers.h). Any nonzero
  // reseedIfShorter is rejected unconditionally -- before leg A or any output
  // is touched -- rather than silently running the non-reseeded algorithm
  // for a configuration where legacy could have taken a different starting
  // point. Distinct from every reason above: no rotate/propagate/material/
  // chi2 arithmetic has been attempted for any leg.
  ReseedNotSupported = 17,
  // M5d Propagator: state-representation-family conversion
  // (Propagator::convertFamily, Propagator.h) could not produce a valid
  // target-family state -- an unrecognized target family, a source
  // Snp/direction too close to the +-1 boundary for the Barrel->Forward
  // asin() step, or a non-finite converted parameter/covariance entry.
  // Distinct from SourceFamilyMismatch
  // (that reason means "no conversion was attempted and the families simply
  // differ"; this one means "a conversion was attempted and failed").
  FamilyConversionFailure = 18
};

static_assert(sizeof(OperationFailureReason) == sizeof(uint8_t));

} // namespace o2::itsmft::tracking

#endif // ALICEO2_ITSMFT_TRACKING_SURFACESTATEOPERATIONRESULT_H_

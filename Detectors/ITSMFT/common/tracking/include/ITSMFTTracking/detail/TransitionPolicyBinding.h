// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYBINDING_H_
#define ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYBINDING_H_

#include "ITSMFTTracking/detail/TrackingKernelParameters.h"
#include "ITSMFTTracking/detail/TransitionPolicy.h"

// Host-only: TrackingParameters (Configuration.h) owns std::vector members
// and is not itself device-compatible; only the bound kernel record crosses
// into device-facing policy operations. This binding step has no device
// counterpart and must never be compiled for device code.
#ifndef GPUCA_GPUCODE

#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include <gsl/span>

namespace o2::itsmft::tracking
{

/// Single shared definition of "is this a MatCorrType value production
/// configuration currently recognizes at all" (USEMatCorrNONE/TGeo/LUT).
/// AttachHitPolicyConfigView::isValid() and the Stage-B material-mode
/// preflight below both call this rather than each encoding their own copy,
/// so the two contracts cannot silently diverge (e.g. if a new MatCorrType
/// value is ever added and only one call site is updated).
inline bool isRecognizedMatCorrType(o2::base::PropagatorF::MatCorrType corrType) noexcept
{
  return corrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE ||
         corrType == o2::base::PropagatorF::MatCorrType::USEMatCorrTGeo ||
         corrType == o2::base::PropagatorF::MatCorrType::USEMatCorrLUT;
}

/// Host view of the authoritative per-surface nominal material (resolved by
/// TrackerTraits::initialiseTimeFrame() into TrackerTraits::mLayerMaterial,
/// see TrackerTraits.h) and the existing propagation-correction configuration
/// consumed by attachHit<Tag>. The view borrows one iteration's already-
/// resolved material array plus the TrackingParameters correction type, and
/// is bound/validated with the typed family Params before traversal starts.
struct AttachHitPolicyConfigView {
  gsl::span<const NominalSurfaceMaterial> layerMaterial;
  o2::base::PropagatorF::MatCorrType corrType{o2::base::PropagatorF::MatCorrType::USEMatCorrNONE};

  bool isValid(size_t expectedLayers) const noexcept
  {
    if (layerMaterial.size() < expectedLayers || !isRecognizedMatCorrType(corrType)) {
      return false;
    }
    for (size_t layer = 0; layer < expectedLayers; ++layer) {
      const auto& material = layerMaterial[layer];
      if (!isFiniteParam(material.xOverX0) || material.xOverX0 < 0.f ||
          !isFiniteParam(material.arealDensityGPerCm2) || material.arealDensityGPerCm2 < 0.f) {
        return false;
      }
    }
    return true;
  }
};

inline AttachHitPolicyConfigView bindAttachHitPolicyConfig(gsl::span<const NominalSurfaceMaterial> layerMaterial,
                                                           const TrackingParameters& params) noexcept
{
  return {layerMaterial, params.CorrType};
}

/// Stage-B Slice C (design report Sec 8/11; Architecture.md Sec 11): result
/// of the pure material-correction-mode preflight below. `Supported` and
/// `Unsupported` both mean the CorrType value is recognized
/// (isRecognizedMatCorrType()); `InvalidMode` means it is not, and callers
/// must defer to the existing configuration-validation failure
/// classification (AttachHitPolicyConfigView::isValid()) rather than
/// treating an invalid value as this preflight's own failure reason.
enum class MaterialCorrectionModeSupport : uint8_t {
  Supported,
  Unsupported,
  InvalidMode
};

/// Stage-B Slice C: pure, host-only preflight for whether one active
/// transition policy's Stage-A native SurfaceKinematicState path
/// (TrackletFinding.h) supports one iteration's configured
/// MatCorrType. Additive and unwired: no production caller exists yet --
/// TrackerTraits::initialiseTimeFrame() does not call this in this slice,
/// and production behavior is unchanged. When a later activation slice
/// wires this in, an `Unsupported` result becomes
/// TraversalFailureReason::UnsupportedMaterialCorrectionMode
/// (TrackerTraits.h): a structural/configuration failure that is always
/// wiped and rethrown by TraversalException, regardless of
/// DropTFUponFailure -- never the dropped-TF sentinel.
///
/// CylinderCylinder native tracking currently supports only
/// USEMatCorrNONE; USEMatCorrLUT/USEMatCorrTGeo remain supported only by the
/// untouched legacy ITS tracker, not by the new native path, so they report
/// `Unsupported` here (recognized, but not yet activated). DiskDisk's native
/// path always uses descriptor-based nominal material regardless of
/// CorrType -- matching the existing common-MFT semantics -- so it is never
/// constrained by this preflight and always reports `Supported` for every
/// recognized CorrType.
///
/// A CorrType value isRecognizedMatCorrType() does not recognize reports
/// `InvalidMode` for every Tag, including DiskDisk: an invalid mode is never
/// silently treated as DiskDisk's always-supported case.
///
/// noexcept, host-only (like AttachHitPolicyConfigView above), and takes/
/// returns only enums by value: no owning container, no TimeFrame or
/// candidate-state dependency, no material spans, no detector identity, no
/// NLayers. Deterministic and side-effect-free -- safe to call once per
/// active policy tag during iteration initialization
/// (the removed tag-grouping/dispatch layer),
/// and idempotent under repeated calls with identical arguments.
///
/// The primary template is `= delete`d: instantiating this for
/// TransitionPolicyTag::Invalid (or any future tag without its own
/// specialization) is a hard compile-time error -- "attempt to use a deleted
/// function" -- rather than a silent fallback to `Supported` or
/// `Unsupported`. This is a stronger, SFINAE-observable guarantee than the
/// declared-but-undefined (link-error) convention used by the
/// TransitionPolicyTag-dispatched operations in TrackletFinding.h
/// (cellsAreCompatible, attachHit, buildCellSeed), chosen here because this
/// preflight's entire purpose is to reject unsupported/invalid configuration
/// before any device-facing code is reached.
template <TransitionPolicyTag Tag>
MaterialCorrectionModeSupport checkMaterialCorrectionModeSupport(o2::base::PropagatorF::MatCorrType corrType) noexcept = delete;

template <>
inline MaterialCorrectionModeSupport checkMaterialCorrectionModeSupport<TransitionPolicyTag::CylinderCylinder>(o2::base::PropagatorF::MatCorrType corrType) noexcept
{
  if (!isRecognizedMatCorrType(corrType)) {
    return MaterialCorrectionModeSupport::InvalidMode;
  }
  return corrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE
           ? MaterialCorrectionModeSupport::Supported
           : MaterialCorrectionModeSupport::Unsupported;
}

template <>
inline MaterialCorrectionModeSupport checkMaterialCorrectionModeSupport<TransitionPolicyTag::DiskDisk>(o2::base::PropagatorF::MatCorrType corrType) noexcept
{
  if (!isRecognizedMatCorrType(corrType)) {
    return MaterialCorrectionModeSupport::InvalidMode;
  }
  return MaterialCorrectionModeSupport::Supported;
}

/// Host view of TrackingParameters::LayerRadii for the transition-preparation
/// slice (layerMultipleScatteringAngle<Tag> / prepareTransitionScatteringAndBending).
/// `layerMaterial` is deliberately not rebound here: it is borrowed directly
/// from the caller's already-bound-and-validated AttachHitPolicyConfigView
/// rather than re-validated independently, so this struct cannot define a
/// second, divergent numeric contract for the same data. isValid() therefore
/// checks span bounds only -- it intentionally adds no numeric constraint on
/// LayerRadii (legacy TimeFrame::initialise() has none: degenerate/zero radii
/// flow through to the existing floating-point behavior unrejected, and this
/// slice must not silently start rejecting them).
struct LayerGeometryConfigView {
  gsl::span<const float> layerRadii;
  gsl::span<const NominalSurfaceMaterial> layerMaterial;

  bool isValid(size_t expectedLayers) const noexcept
  {
    return layerRadii.size() >= expectedLayers && layerMaterial.size() >= expectedLayers;
  }
};

inline LayerGeometryConfigView bindLayerGeometryConfig(const TrackingParameters& params,
                                                       const AttachHitPolicyConfigView& attachHitConfig) noexcept
{
  return {gsl::span<const float>{params.LayerRadii.data(), params.LayerRadii.size()}, attachHitConfig.layerMaterial};
}

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYBINDING_H_ */

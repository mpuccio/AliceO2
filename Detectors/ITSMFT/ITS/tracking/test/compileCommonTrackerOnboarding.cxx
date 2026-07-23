// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include <array>
#include <memory>
#include <vector>

#include <oneapi/tbb/task_arena.h>

#include "ITSCommonTracking/CommonTrackingParameters.h"
#include "ITSMFTTracking/CATracker.h"
#include "ITSMFTTracking/DetectorSurfaceCatalogProvider.h"
#include "ITSMFTTracking/ITSSurfaceCatalogProvider.h"
#include "ITSMFTTracking/TransitionPolicyBinding.h"

namespace
{
using namespace o2::itsmft::tracking;

class CompileOnlyITSCatalogProvider final : public DetectorSurfaceCatalogProvider
{
 public:
  DetectorSurfaceCatalogResult buildCatalog(const DetectorSurfaceCatalogRequest& request) const final
  {
    if (request.detector != o2::detectors::DetID::ITS || request.firstSurface != SurfaceId{0} || request.detectorSurfaceCount != 7) {
      return {{}, DetectorSurfaceCatalogError::InvalidRequest};
    }
    std::vector<SurfaceDescriptor> catalog;
    catalog.reserve(7);
    for (uint16_t layer = 0; layer < 7; ++layer) {
      catalog.push_back({SurfaceId{layer}, layer, static_cast<uint8_t>(o2::detectors::DetID::ITS),
                         SurfaceKind::Cylinder, 0, static_cast<float>(layer + 1), 0.f, 100.f});
    }
    return {std::move(catalog), DetectorSurfaceCatalogError::None};
  }
};

int initializeCommonITSTracker(DetectorSurfaceCatalogProvider& provider)
{
  auto parameters = o2::its::commontracking::makeTrackingParameters(o2::its::TrackingMode::Sync);
  if (parameters.empty()) {
    return 1;
  }

  TimeFrame<7> frame;
  const std::vector<SurfaceId> ordered{SurfaceId{0}, SurfaceId{1}, SurfaceId{2}, SurfaceId{3}, SurfaceId{4}, SurfaceId{5}, SurfaceId{6}};
  const DetectorSurfaceCatalogRequest request{o2::detectors::DetID::ITS, SurfaceId{0}, 7};
  if (!frame.ensureDetectorLayouts(&provider, request, ordered, TransitionPolicyTag::CylinderCylinder, parameters).ok()) {
    return 2;
  }

  // Resolve the authoritative per-surface nominal material from the
  // TimeFrame-owned surface catalogue, mirroring
  // TrackerTraits<NLayers>::initialiseTimeFrame()'s own accessor chain
  // (TrackerTraits.cxx): TimeFrame::getDetectorLayouts() for the ordered
  // SurfaceId mapping, TimeFrame::getDetectorLayoutView() for the resolved
  // per-iteration view, and SurfaceDescriptor::material off each ordered
  // surface -- never a hardcoded or duplicated material value. `layerMaterial`
  // must outlive the borrowed span handed to bindAttachHitPolicyConfig()
  // below, so it stays alive in this same scope through that call.
  const auto* layouts = frame.getDetectorLayouts();
  if (layouts == nullptr) {
    return 4;
  }
  const auto layout = frame.getDetectorLayoutView(0);
  const auto& orderedSurfaces = layouts->getConfigurationKey().orderedSurfaces;
  if (orderedSurfaces.size() < 7) {
    return 4;
  }
  std::array<NominalSurfaceMaterial, 7> layerMaterial{};
  for (int layer = 0; layer < 7; ++layer) {
    const auto surfaceId = orderedSurfaces[layer];
    if (!surfaceId.isValid() || surfaceId.value() >= layout.nSurfaces) {
      return 4;
    }
    layerMaterial[layer] = layout.getSurface(surfaceId).material;
  }

  frame.initDefaultTrackingTopology(parameters.front(), 7);
  frame.initTrackerTopologies(parameters);

  auto pool = std::make_shared<BoundedMemoryResource>();
  frame.setMemoryPool(pool);
  TrackerTraits<7> traits;
  traits.setMemoryPool(pool);
  traits.adoptTimeFrame(&frame);
  traits.updateTrackingParameters(parameters);
  std::shared_ptr<tbb::task_arena> arena;
  traits.setNThreads(1, arena);

  const auto cylinderParameters = bindTransitionPolicyParams<TransitionPolicyTag::CylinderCylinder>(parameters.front());
  const auto materialParameters = bindAttachHitPolicyConfig(
    gsl::span<const NominalSurfaceMaterial>(layerMaterial.data(), layerMaterial.size()), parameters.front());
  if (!cylinderParameters.isValid() || !materialParameters.isValid(7)) {
    return 3;
  }

  Tracker<7> tracker{&traits};
  tracker.adoptTimeFrame(frame);
  tracker.setMemoryPool(pool);
  tracker.setParameters(parameters);
  return 0;
}
} // namespace

int main()
{
  // The geometry-backed production provider is part of the same compile/link
  // proof. The deterministic provider performs the no-event initialization
  // without requiring a runtime geometry file.
  [[maybe_unused]] ITSSurfaceCatalogProvider productionProvider;
  CompileOnlyITSCatalogProvider provider;
  return initializeCommonITSTracker(provider);
}

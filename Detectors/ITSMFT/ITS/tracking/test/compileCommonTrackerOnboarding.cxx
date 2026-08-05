// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include <array>
#include <memory>
#include <vector>

#include <gsl/span>

#include <oneapi/tbb/task_arena.h>

#include "ITSCommonTracking/CommonTrackingParameters.h"
#include "ITSMFTTracking/CATracker.h"
#include "ITSMFTTracking/DetectorLayoutSet.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/SurfaceCatalogView.h"
#include "ITSMFTTracking/detail/TransitionPolicyBinding.h"

namespace
{
using namespace o2::itsmft::tracking;

// kITSStaticSurfaceCatalog is dense and local (StaticDetectorCatalogs.h):
// surface i's id is always SurfaceId{i}. Mirrors identitySurfaceOrder() in
// TrackingInterface.cxx, the production plan builder this link-only proof
// exercises the same way.
constexpr std::array<SurfaceId, ITSNLayers> identityOrder()
{
  std::array<SurfaceId, ITSNLayers> order{};
  for (int i = 0; i < ITSNLayers; ++i) {
    order[i] = SurfaceId{static_cast<uint16_t>(i)};
  }
  return order;
}

int initializeCommonITSTracker()
{
  auto parameters = o2::its::commontracking::makeTrackingParameters(o2::its::TrackingMode::Sync);
  if (parameters.empty()) {
    return 1;
  }

  static constexpr auto kOrderedSurfaces = identityOrder();
  auto planResult = buildDetectorLayoutSet(
    SurfaceCatalogView{kITSStaticSurfaceCatalog.data(), static_cast<uint32_t>(kITSStaticSurfaceCatalog.size())},
    gsl::span<const SurfaceId>{kOrderedSurfaces}, parameters);
  if (!planResult.ok()) {
    return 2;
  }
  const auto& plan = *planResult.layout;

  // Resolve the authoritative per-surface nominal material from the built
  // plan, mirroring TrackerTraits<NLayers>::initialiseTimeFrame()'s own
  // accessor chain (TrackerTraits.cxx): DetectorLayoutSet::getConfigurationKey()
  // for the ordered SurfaceId mapping, DetectorLayoutSet::getLayoutView() for
  // the resolved per-iteration view, and SurfaceDescriptor::material off each
  // ordered surface -- never a hardcoded or duplicated material value.
  // `layerMaterial` must outlive the borrowed span handed to
  // bindAttachHitPolicyConfig() below, so it stays alive in this same scope
  // through that call.
  const auto layout = plan.getLayoutView(0);
  const auto& orderedSurfaces = plan.getConfigurationKey().orderedSurfaces;
  if (orderedSurfaces.size() < ITSNLayers) {
    return 4;
  }
  std::array<NominalSurfaceMaterial, ITSNLayers> layerMaterial{};
  for (int layer = 0; layer < ITSNLayers; ++layer) {
    const auto surfaceId = orderedSurfaces[layer];
    if (!surfaceId.isValid() || surfaceId.value() >= layout.nSurfaces) {
      return 4;
    }
    layerMaterial[layer] = layout.getSurface(surfaceId).material;
  }

  // Gate 4 B3.1: the permanent, non-templated TimeFrame (event data:
  // vertices, beam state, Bz, CommonTrack/TrackClusterReference storage,
  // normalized measurements) and the temporary, per-detector
  // SurfaceTrackingScratch (the common CA scratch/topology containers) is
  // constructed after the permanent TimeFrame so C++'s reverse declaration
  // order tears the scratch down first. Neither owner stores a reference to
  // the other; this function is what binds both.
  TimeFrame frame;
  SurfaceTrackingScratch scratch;

  auto pool = std::make_shared<BoundedMemoryResource>();
  frame.setMemoryPool(pool);
  scratch.setMemoryPool(pool);
  scratch.adoptPlan(orderedSurfaces.size(), layout.topology.nTransitions, layout.topology.nCells);
  scratch.initDefaultTrackingTopology<ITSNLayers>(parameters.front(), ITSNLayers);
  scratch.initTrackerTopologies<ITSNLayers>(parameters);
  TrackerTraits<ITSNLayers> traits;
  traits.setMemoryPool(pool);
  traits.adoptScratch(&scratch);
  traits.adoptFrame(&frame);
  traits.updateTrackingParameters(parameters);
  std::shared_ptr<tbb::task_arena> arena;
  traits.setNThreads(1, arena);

  const auto cylinderParameters = bindTransitionPolicyParams<TransitionPolicyTag::CylinderCylinder>(parameters.front());
  const auto materialParameters = bindAttachHitPolicyConfig(
    gsl::span<const NominalSurfaceMaterial>(layerMaterial.data(), layerMaterial.size()), parameters.front());
  if (!cylinderParameters.isValid() || !materialParameters.isValid(ITSNLayers)) {
    return 3;
  }

  Tracker<ITSNLayers> tracker{&traits};
  tracker.adoptScratch(scratch);
  tracker.adoptFrame(frame);
  tracker.adoptDetectorLayoutSet(plan);
  tracker.setMemoryPool(pool);
  tracker.setParameters(parameters);
  return 0;
}
} // namespace

int main()
{
  return initializeCommonITSTracker();
}

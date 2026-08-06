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
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
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
  const auto& orderedSurfaces = kOrderedSurfaces;
  std::array<NominalSurfaceMaterial, ITSNLayers> layerMaterial{};
  for (int layer = 0; layer < ITSNLayers; ++layer) {
    const auto surfaceId = orderedSurfaces[layer];
    if (!surfaceId.isValid()) {
      return 4;
    }
    layerMaterial[layer] = kITSStaticSurfaceCatalog[layer].material;
  }

  // TimeFrame is the reusable owner of event data and configured workspace;
  // Tracker initialization installs that workspace before the kernel seam is
  // used.
  TimeFrame frame;

  auto pool = std::make_shared<BoundedMemoryResource>();
  TrackerTraits traits;
  traits.setMemoryPool(pool);
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

  Tracker tracker{&traits};
  tracker.setSource(ClusterSourceId{0});
  TrackerInitialization configuration;
  configuration.catalog = SurfaceCatalogView{kITSStaticSurfaceCatalog.data(), static_cast<uint32_t>(kITSStaticSurfaceCatalog.size())};
  configuration.memoryPool = pool;
  TrackerIterationConfiguration iteration;
  SurfaceGraphSubgraph subgraph;
  subgraph.orderedSurfaces.assign(orderedSurfaces.begin(), orderedSurfaces.end());
  subgraph.maxHoles = parameters.front().MaxHoles;
  subgraph.holeSurfaces = positionalSurfaceMask(parameters.front().HoleLayerMask, orderedSurfaces, ITSNLayers);
  subgraph.seedingSurfaces = positionalSurfaceMask(parameters.front().StartLayerMask, orderedSurfaces, ITSNLayers);
  iteration.graphSubgraphs.push_back(std::move(subgraph));
  iteration.parameters.push_back(parameters.front());
  SurfaceMask owned;
  for (const auto surface : orderedSurfaces) {
    owned.set(surface);
  }
  iteration.bindings.push_back(SurfacePlanBinding::Declaration{ClusterSourceId{0}, owned,
                                                               std::vector<SurfaceId>{orderedSurfaces.begin(), orderedSurfaces.end()},
                                                               SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder});
  configuration.iterations.push_back(std::move(iteration));
  const auto result = tracker.initialize(frame, configuration);
  if (!result.ok()) {
    return 5;
  }
  traits.setMemoryPool(frame.getMemoryPool());
  const auto* capacity = frame.getWorkspaceCapacity(0, ClusterSourceId{0});
  if (capacity == nullptr) {
    return 6;
  }
  return 0;
}
} // namespace

int main()
{
  return initializeCommonITSTracker();
}

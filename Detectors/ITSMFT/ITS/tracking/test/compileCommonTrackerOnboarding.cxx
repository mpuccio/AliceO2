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
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/detail/CandidateFinding.h"

namespace
{
using namespace o2::itsmft::tracking;

int initializeCommonITSTracker()
{
  auto parameters = o2::its::commontracking::makeTrackingParameters(o2::its::TrackingMode::Sync);
  if (parameters.empty()) {
    return 1;
  }

  std::array<NominalSurfaceMaterial, ITSNLayers> layerMaterial{};
  for (int layer = 0; layer < ITSNLayers; ++layer) {
    layerMaterial[layer] = kITSStaticSurfaceCatalog[layer].material;
  }

  // TimeFrame is the reusable owner of event data and configured workspace;
  // Tracker initialization installs that workspace before the kernel seam is
  // used.
  TimeFrame frame;

  auto pool = std::make_shared<BoundedMemoryResource>();
  TrackerTraits traits;
  std::shared_ptr<tbb::task_arena> arena;
  traits.setNThreads(1, arena);

  TrackingKernelParameters cylinderParameters;
  const auto materialParameters = bindAttachHitConfig(
    gsl::span<const NominalSurfaceMaterial>(layerMaterial.data(), layerMaterial.size()), parameters.front());
  if (!cylinderParameters.isValid() || !materialParameters.isValid(ITSNLayers)) {
    return 3;
  }

  Tracker tracker;
  TrackerInitialization configuration;
  configuration.catalog = SurfaceCatalogView{kITSStaticSurfaceCatalog.data(), static_cast<uint32_t>(kITSStaticSurfaceCatalog.size())};
  configuration.memoryPool = pool;
  configuration.layout = makeDetectorLayout();
  configuration.parameters = parameters;
  const auto result = tracker.initialize(frame, configuration);
  if (!result.ok()) {
    return 5;
  }
  if (frame.getLayout().size() != ITSNLayers) {
    return 6;
  }
  return 0;
}
} // namespace

int main()
{
  return initializeCommonITSTracker();
}

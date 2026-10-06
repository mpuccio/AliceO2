// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#ifndef ALICEO2_ITSMFT_TRACKERTRAITSGPU_H_
#define ALICEO2_ITSMFT_TRACKERTRAITSGPU_H_

#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTrackingGPU/TimeFrameGPU.h"

namespace o2::itsmft::tracking
{
namespace gpu
{
inline namespace ITSMFT_TRACKING_GPU_BACKEND
{
// Every CA stage on the device, on the data TimeFrameGPU keeps there. With
// validation a stage also runs on the CPU and the results are compared.
class TrackerTraitsGPU final : public TrackerTraits
{
 public:
  const char* getName() const noexcept final;
  bool isGPU() const noexcept final { return true; }
  void beginTimeframe() final { mTimeFrameGPU.initialise(0); } // loaded again by the first tracklet stage
  // The device builds its own tables; the host ones serve only the CPU
  // tracklet stage that validation runs.
  bool usesHostIndexTables() const noexcept final { return validateTracklets(); }
  void setFrameworkAllocator(ExternalAllocator* allocator) final { mTimeFrameGPU.setFrameworkAllocator(allocator); }

 protected:
  void computeLayerTracklets(IterationContext& context, int iteration, int iVertex) final;
  void computeLayerCells(IterationContext& context, int iteration) final;
  void findCellsNeighbours(IterationContext& context, int iteration) final;
  void findRoads(IterationContext& context, int iteration) final;
  // As the legacy tracker: what the stages allocate on the device is scoped
  // to one pass.
  void beginPass(IterationContext& context) final { mTimeFrameGPU.pushMemoryStack(context.iteration); }
  void endPass(IterationContext& context) final { mTimeFrameGPU.popMemoryStack(context.iteration); }
  void usedClustersChanged(IterationContext& context) final { mTimeFrameGPU.loadUsedClusters(context.frame); }
  void adoptHostTracklets(IterationContext& context) final;

 private:
  // Copies the cells, with their levels, and the neighbour graph to the host.
  void downloadRoadGraph(IterationContext& context);

  TimeFrameGPU mTimeFrameGPU;
  bounded_vector<TrackingCandidate> mTracks; // reused for every start level
};
} // namespace ITSMFT_TRACKING_GPU_BACKEND
} // namespace gpu

std::unique_ptr<TrackerTraits> createTrackerTraitsCUDA();
std::unique_ptr<TrackerTraits> createTrackerTraitsHIP();
} // namespace o2::itsmft::tracking
#endif

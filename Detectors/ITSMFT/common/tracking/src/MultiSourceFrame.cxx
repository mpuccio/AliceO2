// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/MultiSourceFrame.h"

#include "SimulationDataFormat/MCTruthContainer.h"

namespace o2::itsmft::tracking
{

void MultiSourceFrame::clear() noexcept
{
  mMeasurements.clear();
  mSurfaceRanges.clear();
  mSources.clear();
  mROFIntervals.clear();
  mSourceROFOffsets.clear();
  mLabelSources.clear();
}

MultiSourceFrameView MultiSourceFrame::getView() const noexcept
{
  MultiSourceFrameView view;
  view.measurements = mMeasurements.data();
  view.nMeasurements = static_cast<uint32_t>(mMeasurements.size());
  view.surfaceRanges = mSurfaceRanges.data();
  view.nSurfaces = static_cast<uint32_t>(mSurfaceRanges.size());
  view.rofIntervals = mROFIntervals.data();
  view.sourceROFOffsets = mSourceROFOffsets.data();
  view.nSources = static_cast<uint32_t>(mSources.size());
  return view;
}

gsl::span<const SurfaceMeasurement> MultiSourceFrame::getSurfaceMeasurements(SurfaceId surface) const
{
  if (!surface.isValid() || surface.value() >= mSurfaceRanges.size()) {
    return {};
  }
  const auto& range = mSurfaceRanges[surface.value()];
  if (range.getEntries() == 0) {
    return {};
  }
  return {mMeasurements.data() + range.getFirstEntry(), range.getEntries()};
}

gsl::span<const ROFIntervalBC> MultiSourceFrame::getSourceIntervals(ClusterSourceId source) const
{
  if (!source.isValid() || source.value() + 1 >= mSourceROFOffsets.size()) {
    return {};
  }
  const auto first = mSourceROFOffsets[source.value()];
  const auto last = mSourceROFOffsets[source.value() + 1];
  if (last == first) {
    return {};
  }
  return {mROFIntervals.data() + first, last - first};
}

gsl::span<const o2::MCCompLabel> MultiSourceFrame::getLabels(ClusterRef cluster) const
{
  if (!cluster.isValid() || cluster.source.value() >= mLabelSources.size()) {
    return {};
  }
  const auto* container = mLabelSources[cluster.source.value()];
  if (container == nullptr) {
    return {};
  }
  return container->getLabels(cluster.index);
}

void MultiSourceFrame::assignLoadedData(std::vector<SurfaceMeasurement>&& measurements,
                                        std::vector<SurfaceMeasurementRange>&& surfaceRanges,
                                        std::vector<SourceMetadata>&& sources,
                                        std::vector<ROFIntervalBC>&& rofIntervals,
                                        std::vector<uint32_t>&& sourceROFOffsets,
                                        std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>&& labelSources)
{
  mMeasurements = std::move(measurements);
  mSurfaceRanges = std::move(surfaceRanges);
  mSources = std::move(sources);
  mROFIntervals = std::move(rofIntervals);
  mSourceROFOffsets = std::move(sourceROFOffsets);
  mLabelSources = std::move(labelSources);
}

} // namespace o2::itsmft::tracking

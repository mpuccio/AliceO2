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
  mPerSurfaceMeasurements.clear();
  mPerSurfaceGlobalMeasurements.clear();
  mSurfaceSpans.clear();
  mSources.clear();
  mROFIntervals.clear();
  mSourceROFOffsets.clear();
  mLabelSources.clear();
}

void MultiSourceFrame::rebuildSurfaceSpans()
{
  mSurfaceSpans.resize(mPerSurfaceMeasurements.size());
  for (size_t s = 0; s < mPerSurfaceMeasurements.size(); ++s) {
    const auto& measurements = mPerSurfaceMeasurements[s];
    const auto& globals = mPerSurfaceGlobalMeasurements[s];
    mSurfaceSpans[s] = SurfaceMeasurementSpan{globals.data(), measurements.data(), static_cast<uint32_t>(measurements.size())};
  }
}

MultiSourceFrameView MultiSourceFrame::getView() const noexcept
{
  MultiSourceFrameView view;
  view.surfaces = mSurfaceSpans.data();
  view.nSurfaces = static_cast<uint32_t>(mSurfaceSpans.size());
  view.rofIntervals = mROFIntervals.data();
  view.sourceROFOffsets = mSourceROFOffsets.data();
  view.nSources = static_cast<uint32_t>(mSources.size());
  return view;
}

gsl::span<const SurfaceMeasurement> MultiSourceFrame::getSurfaceMeasurements(SurfaceId surface) const
{
  if (!surface.isValid() || surface.value() >= mPerSurfaceMeasurements.size()) {
    return {};
  }
  const auto& measurements = mPerSurfaceMeasurements[surface.value()];
  return {measurements.data(), measurements.size()};
}

gsl::span<const GlobalMeasurement> MultiSourceFrame::getGlobalMeasurements(SurfaceId surface) const
{
  if (!surface.isValid() || surface.value() >= mPerSurfaceGlobalMeasurements.size()) {
    return {};
  }
  const auto& measurements = mPerSurfaceGlobalMeasurements[surface.value()];
  return {measurements.data(), measurements.size()};
}

const GlobalMeasurement* MultiSourceFrame::getGlobalMeasurement(SurfaceId surface, SurfaceMeasurementIndex index) const noexcept
{
  if (!surface.isValid() || surface.value() >= mPerSurfaceGlobalMeasurements.size() || !index.isValid()) {
    return nullptr;
  }
  const auto& measurements = mPerSurfaceGlobalMeasurements[surface.value()];
  return index.value() < measurements.size() ? &measurements[index.value()] : nullptr;
}

const SurfaceMeasurement* MultiSourceFrame::getSurfaceMeasurement(SurfaceId surface, SurfaceMeasurementIndex index) const noexcept
{
  if (!surface.isValid() || surface.value() >= mPerSurfaceMeasurements.size() || !index.isValid()) {
    return nullptr;
  }
  const auto& measurements = mPerSurfaceMeasurements[surface.value()];
  return index.value() < measurements.size() ? &measurements[index.value()] : nullptr;
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

void MultiSourceFrame::assignLoadedData(std::vector<std::vector<GlobalMeasurement>>&& perSurfaceGlobalMeasurements,
                                        std::vector<std::vector<SurfaceMeasurement>>&& perSurfaceMeasurements,
                                        std::vector<SourceMetadata>&& sources,
                                        std::vector<ROFIntervalBC>&& rofIntervals,
                                        std::vector<uint32_t>&& sourceROFOffsets,
                                        std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>&& labelSources)
{
  mPerSurfaceGlobalMeasurements = std::move(perSurfaceGlobalMeasurements);
  mPerSurfaceMeasurements = std::move(perSurfaceMeasurements);
  mSources = std::move(sources);
  mROFIntervals = std::move(rofIntervals);
  mSourceROFOffsets = std::move(sourceROFOffsets);
  mLabelSources = std::move(labelSources);
  rebuildSurfaceSpans();
}

} // namespace o2::itsmft::tracking

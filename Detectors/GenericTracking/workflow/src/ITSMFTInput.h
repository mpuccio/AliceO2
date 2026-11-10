// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef GTRK_WORKFLOW_ITSMFTINPUT_H_
#define GTRK_WORKFLOW_ITSMFTINPUT_H_

#include "GenericTracking/Alice/Clusters.h"
#include "GenericTracking/Base/FrameTiming.h"
#include "DataFormatsITSMFT/ClusterID.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "Framework/InputSpec.h"
#include "Framework/InputRecord.h"
#include "Framework/ProcessingContext.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include <gsl/span>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace o2::gtrk::workflow
{

inline std::string clusterInputBinding(const char* name, int stream, int streams)
{
  return streams == 1 ? name : std::format("{}_{}", name, stream);
}

/// The same subscription and ordering for the converter and the output adaptor.
inline void addITSMFTClusterInputs(std::vector<o2::framework::InputSpec>& inputs, o2::header::DataOrigin origin, bool useMC, int streams)
{
  using o2::framework::Lifetime;
  for (int stream = 0; stream < streams; ++stream) {
    inputs.emplace_back(clusterInputBinding("compclusters", stream, streams), origin, "COMPCLUSTERS", stream, Lifetime::Timeframe);
    inputs.emplace_back(clusterInputBinding("patterns", stream, streams), origin, "PATTERNS", stream, Lifetime::Timeframe);
    inputs.emplace_back(clusterInputBinding("rofs", stream, streams), origin, "CLUSTERSROF", stream, Lifetime::Timeframe);
    if (useMC) {
      inputs.emplace_back(clusterInputBinding("labels", stream, streams), origin, "CLUSTERSMCTR", stream, Lifetime::Timeframe);
    }
  }
}

/// GTRK indices address the streams concatenated in layer order. O2 references
/// retain the input stream and its local index, including when a stream is empty.
class ITSMFTInput
{
 public:
  using Labels = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
  struct Stream {
    std::span<const o2::itsmft::CompClusterExt> compact;
    std::span<const unsigned char> patterns;
    std::span<const o2::itsmft::ROFRecord> rofs;
    std::unique_ptr<const Labels> labels;
  };
  struct Source {
    std::size_t stream;
    std::uint32_t index;
  };

  explicit ITSMFTInput(std::vector<Stream> streams) : mStreams(std::move(streams))
  {
    if (mStreams.empty() || mStreams.size() > o2::itsmft::MaxClusLayers) {
      throw std::invalid_argument("ITS/MFT cluster input: invalid stream count");
    }
    for (const auto& stream : mStreams) {
      if (stream.compact.size() > static_cast<std::size_t>(o2::itsmft::ClusIndexMask) + 1U ||
          stream.compact.size() > std::numeric_limits<std::uint32_t>::max() - mOffsets.back()) {
        throw std::invalid_argument("ITS/MFT cluster input: too many clusters");
      }
      if (stream.labels && stream.labels->getIndexedSize() != stream.compact.size() && !stream.compact.empty()) {
        throw std::invalid_argument("ITS/MFT cluster input: labels do not cover the clusters");
      }
      mOffsets.push_back(mOffsets.back() + static_cast<std::uint32_t>(stream.compact.size()));
    }
  }

  [[nodiscard]] std::span<const Stream> streams() const { return mStreams; }
  [[nodiscard]] std::uint32_t size() const { return mOffsets.back(); }
  [[nodiscard]] std::uint32_t offset(std::size_t stream) const { return mOffsets.at(stream); }

  [[nodiscard]] Source source(std::uint32_t index) const
  {
    if (index >= size()) {
      throw std::out_of_range("ITS/MFT cluster input: cluster index out of range");
    }
    const auto stream = static_cast<std::size_t>(std::upper_bound(mOffsets.begin(), mOffsets.end(), index) - mOffsets.begin() - 1);
    return {stream, index - mOffsets[stream]};
  }

  [[nodiscard]] int externalIndex(std::uint32_t index) const
  {
    const auto from = source(index);
    return o2::itsmft::composeClusID(static_cast<int>(from.stream), static_cast<int>(from.index));
  }

  [[nodiscard]] const o2::itsmft::CompClusterExt& cluster(std::uint32_t index) const
  {
    const auto from = source(index);
    return mStreams[from.stream].compact[from.index];
  }

  [[nodiscard]] auto labels(std::uint32_t index) const
  {
    const auto from = source(index);
    return mStreams[from.stream].labels->getLabels(from.index);
  }

 private:
  std::vector<Stream> mStreams;
  std::vector<std::uint32_t> mOffsets{0U};
};

/// Decode streams separately, then concatenate clusters and their labels in
/// the same order. The run's frame lengths also cover single-ROF inputs.
inline std::vector<Cluster> convertITSMFTInput(const ITSMFTInput& input, const o2::itsmft::TopologyDictionary& dictionary,
                                              std::span<const FrameTiming> timings, std::uint32_t orbit)
{
  if (timings.size() < input.streams().size()) {
    throw std::invalid_argument("ITS/MFT input: missing per-layer timing");
  }
  std::vector<Cluster> result;
  result.reserve(input.size());
  for (std::size_t stream = 0; stream < input.streams().size(); ++stream) {
    const auto& source = input.streams()[stream];
    if (source.rofs.size() > 1U && alice::readoutFrames(source.rofs).lengthBC != timings[stream].lengthBC) {
      throw std::runtime_error(std::format("ITS/MFT stream {}: ROF spacing differs from the run's frame length {} BC", stream, timings[stream].lengthBC));
    }
    const auto local = alice::localClusters(source.compact, source.patterns, dictionary);
    const auto clusters = alice::toClusters(local, source.rofs, {.lengthBC = timings[stream].lengthBC}, orbit);
    result.insert(result.end(), clusters.begin(), clusters.end());
  }
  return result;
}

inline ITSMFTInput::Labels concatenateLabels(const ITSMFTInput& input)
{
  ITSMFTInput::Labels labels;
  for (const auto& stream : input.streams()) {
    if (!stream.compact.empty()) {
      labels.mergeAtBack(*stream.labels);
    }
  }
  return labels;
}

inline ITSMFTInput readITSMFTInput(o2::framework::ProcessingContext& pc, bool useMC, int streams)
{
  std::vector<ITSMFTInput::Stream> result;
  result.reserve(streams);
  for (int stream = 0; stream < streams; ++stream) {
    const auto compact = pc.inputs().get<gsl::span<o2::itsmft::CompClusterExt>>(clusterInputBinding("compclusters", stream, streams).c_str());
    const auto patterns = pc.inputs().get<gsl::span<unsigned char>>(clusterInputBinding("patterns", stream, streams).c_str());
    const auto rofs = pc.inputs().get<gsl::span<o2::itsmft::ROFRecord>>(clusterInputBinding("rofs", stream, streams).c_str());
    result.push_back({.compact = {compact.data(), compact.size()}, .patterns = {patterns.data(), patterns.size()}, .rofs = {rofs.data(), rofs.size()},
                      .labels = useMC ? pc.inputs().get<const ITSMFTInput::Labels*>(clusterInputBinding("labels", stream, streams).c_str()) : nullptr});
  }
  return ITSMFTInput{std::move(result)};
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_ITSMFTINPUT_H_

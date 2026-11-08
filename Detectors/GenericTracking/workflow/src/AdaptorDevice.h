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

/// \file AdaptorDevice.h
/// \brief The O2-format adaptor of the ITS and the MFT: the tracking device's
///        products as what O2's own tracker of the detector publishes.
///
/// The ITS's: ITS/TRACKS, TRACKCLSID, ITSTrackROF, VERTICES, VERTICESROF,
/// IRFRAMES, with MC TRACKSMCTR, VERTICESMCTR, VERTICESMCPUR; the MFT's:
/// MFT/TRACKS, TRACKCLSID, MFTTrackROF, with MC TRACKSMCTR; both META/<DET>TRACKER
/// once. Records whose refit failed are left out. The parameters are the
/// states' own (O2Parameters.h), the labels the study's rule (TrackLabel.h),
/// a track's time and readout frame ITS's rule on its clusters' frames.

#ifndef GTRK_WORKFLOW_ADAPTORDEVICE_H_
#define GTRK_WORKFLOW_ADAPTORDEVICE_H_

#include "GenericTracking/Alice/Clusters.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Workflow/AdaptorSpec.h"
#include "GenericTracking/Workflow/Conditions.h"
#include "GenericTracking/Workflow/DataSpecs.h"
#include "GenericTracking/Workflow/GeometryProvider.h"
#include "GenericTracking/Workflow/GridConfig.h"
#include "GenericTracking/Workflow/O2Parameters.h"
#include "GenericTracking/Workflow/RecoParam.h"
#include "GenericTracking/Workflow/TrackLabel.h"
#include "GenericTracking/Workflow/TrackerSpec.h"

#include "CommonConstants/LHCConstants.h"
#include "CommonDataFormat/IRFrame.h"
#include "CommonUtils/ConfigurableParam.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITS/Vertex.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/DeviceSpec.h"
#include "Framework/Logger.h"
#include "Framework/Task.h"
#include "Framework/TimingInfo.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITSMFTTracking/ROFViews.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include <TMap.h>
#include <TObjString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <span>
#include <string>
#include <type_traits>
#include <vector>

namespace o2::gtrk::workflow
{

namespace adaptor
{

/// A particle's label, the fake flag unless \a good; O2's noise label for no
/// particle.
[[nodiscard]] inline o2::MCCompLabel labelOf(std::uint64_t particle, bool good)
{
  return particle == kNoTrackParticle ? o2::MCCompLabel(true) : labelOfParticle(particle, !good);
}

/// A timeframe's clusters as the adaptor needs them: each one's readout frame
/// (its record in CLUSTERSROF) and that frame's start in BC from the
/// timeframe's start, without the delay and bias.
struct ClusterFrames {
  std::vector<std::uint32_t> record;  ///< per cluster
  std::vector<std::int64_t> recordBC; ///< per record
};

[[nodiscard]] inline ClusterFrames clusterFrames(std::span<const o2::itsmft::ROFRecord> rofs, std::size_t clusters, std::uint32_t tfFirstOrbit)
{
  ClusterFrames result;
  result.record.assign(clusters, std::numeric_limits<std::uint32_t>::max());
  result.recordBC.reserve(rofs.size());
  for (std::size_t r = 0; r < rofs.size(); ++r) {
    const auto& ir = rofs[r].getBCData();
    result.recordBC.push_back((static_cast<std::int64_t>(ir.orbit - tfFirstOrbit) * o2::constants::lhc::LHCMaxBunches) + ir.bc);
    for (int i = 0; i < rofs[r].getNEntries(); ++i) {
      result.record[static_cast<std::size_t>(rofs[r].getFirstEntry() + i)] = static_cast<std::uint32_t>(r);
    }
  }
  return result;
}

/// A track's time as ITS takes it (`TrackerTraits`): the intersection of its
/// clusters' frames `[start + delay + bias, + length)` in BC from the
/// timeframe's start; where they do not intersect -- a track across frames,
/// which only a margin lets the tracker make -- the record's box, the frames
/// widened by the margin. [lower, upper) in BC.
struct BCRange {
  std::int64_t lower = 0;
  std::int64_t upper = 0;
};

template <typename Record>
[[nodiscard]] BCRange trackTime(const Record& record, std::span<const std::uint32_t> clusterIndex, const ClusterFrames& frames, const FrameTiming& timing)
{
  BCRange nominal{.lower = std::numeric_limits<std::int64_t>::min(), .upper = std::numeric_limits<std::int64_t>::max()};
  for (std::uint8_t i = 0; i < record.count; ++i) {
    const auto start = frames.recordBC[frames.record[clusterIndex[record.clusters.index[i]]]] + timing.firstStartBC;
    nominal.lower = std::max(nominal.lower, start);
    nominal.upper = std::min(nominal.upper, start + timing.lengthBC);
  }
  if (nominal.lower < nominal.upper) {
    return nominal;
  }
  return {.lower = bcOfPS(record.timePS), .upper = bcOfPS(static_cast<std::int64_t>(record.timePS) + record.widthPS)};
}

/// The earliest frame (its CLUSTERSROF record) among a track's clusters.
template <typename Record>
[[nodiscard]] std::uint32_t earliestRecord(const Record& record, std::span<const std::uint32_t> clusterIndex, const ClusterFrames& frames)
{
  std::uint32_t earliest = std::numeric_limits<std::uint32_t>::max();
  for (std::uint8_t i = 0; i < record.count; ++i) {
    earliest = std::min(earliest, frames.record[clusterIndex[record.clusters.index[i]]]);
  }
  return earliest;
}

/// Each cluster's pixel count, the converter's decoding (`localClusters`).
[[nodiscard]] inline std::vector<std::uint32_t> pixelCounts(std::span<const o2::itsmft::CompClusterExt> compact, std::span<const unsigned char> patterns,
                                                            const o2::itsmft::TopologyDictionary& dictionary)
{
  const auto local = alice::localClusters(compact, patterns, dictionary);
  std::vector<std::uint32_t> result(local.size());
  std::ranges::transform(local, result.begin(), [](const auto& cluster) { return cluster.nPixels; });
  return result;
}

/// The track's particle by the shared rule (TrackLabel.h), its rows'
/// clusters' labels, innermost first; \a carried gets, per row, whether it
/// carries that particle.
template <typename Record>
[[nodiscard]] o2::MCCompLabel trackLabel(const Record& record, std::span<const std::uint32_t> clusterIndex,
                                         const o2::dataformats::MCTruthContainer<o2::MCCompLabel>& labels, std::span<bool> carried)
{
  constexpr std::size_t kRows = std::extent_v<std::remove_cvref_t<decltype(record.clusters.index)>>;
  std::array<RowParticles, kRows> rows{};
  std::array<std::uint32_t, kRows> track{};
  for (std::uint8_t i = 0; i < record.count; ++i) {
    const auto all = labels.getLabels(clusterIndex[record.clusters.index[i]]);
    rows[i] = rowParticles(std::span<const o2::MCCompLabel>{all.data(), all.size()});
    track[i] = i;
  }
  const std::span<const std::uint32_t> order{track.data(), record.count};
  const auto majority = majorityParticle(rows, order);
  const bool good = majority.particle != kNoTrackParticle && majority.carried == record.count;
  for (std::uint8_t i = 0; i < record.count; ++i) {
    carried[i] = majority.particle != kNoTrackParticle && rows[i].carries(majority.particle);
  }
  return labelOf(majority.particle, good);
}

/// The description of O2's tracker metadata of \a Det, META/<DET>TRACKER.
template <DetectorLayout Det>
[[nodiscard]] o2::header::DataDescription metaDescription()
{
  o2::header::DataDescription result;
  result.runtimeInit((detectorOrigin<Det>().template as<std::string>() + "TRACKER").c_str());
  return result;
}

/// META/<DET>TRACKER once, on the first timeslice: the configuration that
/// made the tracks, `GTRKReco<DET>` as JSON, as O2's trackers publish theirs.
template <DetectorLayout Det>
void publishConfiguration(o2::framework::ProcessingContext& pc, bool& published)
{
  if (published) {
    return;
  }
  published = true;
  if (pc.services().get<const o2::framework::DeviceSpec>().inputTimesliceId != 0) {
    return;
  }
  const auto name = RecoParamOf<Det>::type::Instance().getName();
  TMap md;
  md.SetOwnerKeyValue();
  md.Add(new TObjString(name.c_str()), new TObjString(o2::conf::ConfigurableParam::asJSON(name).c_str()));
  pc.outputs().snapshot(o2::framework::Output{"META", metaDescription<Det>(), 0}, md);
}

/// The device's own inputs and conditions, both detectors'.
template <DetectorLayout Det>
class AdaptorTask : public o2::framework::Task
{
 public:
  using Record = TrackRecord<Det>;
  using Labels = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
  static constexpr bool kITS = std::same_as<Det, detectors::ITS>;

  AdaptorTask(bool useMC, std::shared_ptr<Conditions<Det>> conditions) : mUseMC(useMC), mConditions(std::move(conditions)) {}

  void init(o2::framework::InitContext& ic) final { mConditions->init(ic); }

  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) final
  {
    if (matcher == o2::framework::ConcreteDataMatcher(detectorOrigin<Det>(), "CLUSDICT", 0)) {
      mDictionary = static_cast<const o2::itsmft::TopologyDictionary*>(obj);
      return;
    }
    mConditions->finaliseCCDB(matcher, obj);
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    mConditions->checkUpdates(pc);
    pc.inputs().get<o2::itsmft::TopologyDictionary*>("dictionary");
    adaptor::publishConfiguration<Det>(pc, mConfigurationPublished);
    if (!mGeometry) {
      // The run's geometry, as the tracker's (TrackingInterface::configure).
      mGeometry = std::make_unique<TrackingGeometry<Det>>(*makeGeometry<Det>(mConditions->measurementDepthsCM()));
    }
    const auto orbit = pc.services().get<o2::framework::TimingInfo>().firstTForbit;
    const auto records = pc.inputs().get<gsl::span<Record>>("tracks");
    const auto rows = pc.inputs().get<gsl::span<std::uint32_t>>("clusterIndex");
    const auto rofs = pc.inputs().get<gsl::span<o2::itsmft::ROFRecord>>("rofs");
    const auto compact = pc.inputs().get<gsl::span<o2::itsmft::CompClusterExt>>("compclusters");
    const auto patterns = pc.inputs().get<gsl::span<unsigned char>>("patterns");
    std::unique_ptr<const Labels> labels;
    if (mUseMC) {
      labels = pc.inputs().get<const Labels*>("labels");
    }
    const std::span<const std::uint32_t> clusterIndex{rows.data(), rows.size()};
    const std::span<const o2::itsmft::ROFRecord> frameRecords{rofs.data(), rofs.size()};
    const auto frames = adaptor::clusterFrames(frameRecords, compact.size(), orbit);
    const auto pixels = compact.empty() ? std::vector<std::uint32_t>{}
                                        : adaptor::pixelCounts({compact.data(), compact.size()}, {patterns.data(), patterns.size()}, *mDictionary);
    const auto timing = mConditions->frameTiming().front();
    const double bz = nominalBzKiloGauss();

    // The fitted tracks, in the tracker's order.
    std::vector<std::uint32_t> fitted;
    fitted.reserve(records.size());
    for (std::uint32_t track = 0; track < records.size(); ++track) {
      if (records[track].status == PropagationStatus::kOk) {
        fitted.push_back(track);
      }
    }
    const auto layerOf = [&](std::uint32_t row) { return mGeometry->surface(compact[clusterIndex[row]].getSensorID()).layer; };
    if constexpr (kITS) {
      runITS(pc, records, fitted, clusterIndex, frameRecords, frames, pixels, labels.get(), timing, bz, orbit, layerOf);
    } else {
      runMFT(pc, records, fitted, clusterIndex, frameRecords, frames, pixels, labels.get(), bz, orbit, layerOf);
    }
  }

 private:
  template <typename LayerOf>
  void runITS(o2::framework::ProcessingContext& pc, gsl::span<const Record> records, std::span<const std::uint32_t> fitted, std::span<const std::uint32_t> clusterIndex,
              std::span<const o2::itsmft::ROFRecord> rofs, const adaptor::ClusterFrames& frames, std::span<const std::uint32_t> pixels, const Labels* labels,
              const FrameTiming& timing, double bz, std::uint32_t orbit, const LayerOf& layerOf)
  {
    using o2::framework::Output;
    using AlpideParam = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>;
    const auto& origin = detectorOrigin<Det>();
    auto& irFrames = pc.outputs().make<std::vector<o2::dataformats::IRFrame>>(Output{origin, "IRFRAMES", 0});
    auto& clusterIds = pc.outputs().make<std::vector<int>>(Output{origin, "TRACKCLSID", 0});
    auto& tracks = pc.outputs().make<std::vector<o2::its::TrackITS>>(Output{origin, "TRACKS", 0});
    auto& trackROFs = pc.outputs().make<std::vector<o2::itsmft::ROFRecord>>(Output{origin, "ITSTrackROF", 0});
    auto& vertices = pc.outputs().make<std::vector<o2::its::Vertex>>(Output{origin, "VERTICES", 0});
    auto& vertexROFs = pc.outputs().make<std::vector<o2::itsmft::ROFRecord>>(Output{origin, "VERTICESROF", 0});

    const auto length = static_cast<std::int64_t>(timing.lengthBC);
    const auto delay = static_cast<std::int64_t>(AlpideParam::Instance().getROFDelayInBC(0));
    using BC = o2::itsmft::tracking::ROFTimingLayer::BCType;
    // ITS's frame of a time: its lower edge's, 0 before the first.
    const o2::itsmft::tracking::ROFTimingLayer rofTiming{.mROFLength = static_cast<BC>(length),
                                                         .mROFDelay = static_cast<BC>(delay),
                                                         .mROFBias = static_cast<BC>(AlpideParam::Instance().getROFBiasInBC(0))};
    // A vertex's frame as a track's: the lower edge of its time with the half
    // width at most half a frame. Its box's own lower edge lies a margin
    // before its frame starts, in the frame before.
    const auto vertexFrame = [&](const o2::its::Vertex& vertex) {
      const auto time = vertex.getTimeStamp().makeSymmetrical();
      return static_cast<std::int64_t>(
        rofTiming.getROF(o2::its::TimeStamp(time.getTimeStamp(), std::min(time.getTimeStampError(), static_cast<float>(length) / 2.F))));
    };

    // A cluster in more than one track's innermost place: ITS's shared flag.
    std::vector<std::uint8_t> innermost(clusterIndex.empty() ? 0U : *std::ranges::max_element(clusterIndex) + 1U, 0U);
    for (const auto track : fitted) {
      auto& uses = innermost[clusterIndex[records[track].clusters.index[0]]];
      uses = static_cast<std::uint8_t>(std::min(uses + 1, 2));
    }

    struct Converted {
      o2::its::TrackITS track;
      o2::MCCompLabel label;
      std::uint32_t record = 0;
      std::int64_t frame = 0;
    };
    std::vector<Converted> converted;
    converted.reserve(fitted.size());
    std::size_t turned = 0U;
    std::array<bool, TrackRefit<Det>::kMaxClusters> carried{};
    for (const auto track : fitted) {
      const auto& record = records[track];
      const auto geometry = o2::its::GeometryTGeo::Instance();
      const double sensorAlpha = geometry->getSensorRefAlpha(static_cast<int>(record.state.surface));
      const double alpha = frameAlpha(record.state, sensorAlpha);
      turned += std::cos(static_cast<double>(record.state.phi()) - sensorAlpha) < 0. ? 1U : 0U;
      auto& out = converted.emplace_back();
      out.record = track;
      out.track = o2::its::TrackITS{toTrackParCov(record.state, mGeometry->planeOf(record.state), alpha, bz)};
      const double outerAlpha = frameAlpha(record.outwardState, geometry->getSensorRefAlpha(static_cast<int>(record.outwardState.surface)));
      out.track.getParamOut() = toTrackParCov(record.outwardState, mGeometry->planeOf(record.outwardState), outerAlpha, bz);
      out.track.setChi2(record.chi2);
      out.track.setNumberOfClusters(record.count);
      if (labels != nullptr) {
        out.label = adaptor::trackLabel(record, clusterIndex, *labels, std::span<bool>{carried.data(), record.count});
      }
      std::uint32_t pattern = 0U;
      for (std::uint8_t i = 0; i < record.count; ++i) {
        const auto layer = layerOf(record.clusters.index[i]);
        pattern |= 1U << layer;
        pattern |= labels != nullptr && carried[i] ? 1U << (16U + layer) : 0U;
        out.track.setClusterSize(layer, static_cast<int>(pixels[clusterIndex[record.clusters.index[i]]]));
      }
      out.track.setPattern(pattern);
      out.track.setSharedClusters(innermost[clusterIndex[record.clusters.index[0]]] > 1U);
      const auto time = adaptor::trackTime(record, clusterIndex, frames, timing);
      const float half = std::min(static_cast<float>(time.upper - time.lower) / 2.F, static_cast<float>(length) / 2.F);
      const float centre = static_cast<float>(time.lower + time.upper) / 2.F;
      out.track.getTimeStamp() = o2::its::TimeStamp(centre, half);
      out.frame = rofTiming.getROF(out.track.getTimeStamp());
    }
    // ITS's order: by the time's lower edge, the tracker's order within.
    std::ranges::stable_sort(converted, {}, [](const Converted& c) { return c.track.getTimeStamp().getTimeStamp() - c.track.getTimeStamp().getTimeStampError(); });

    // The own vertexer's vertices, at the beam line's x and y.
    const auto shipped = pc.inputs().get<MeasurementGeometry<Det>>("measurementGeometry");
    const auto seeds = pc.inputs().get<gsl::span<VertexSeed>>("vertices");
    const auto contributors = pc.inputs().get<gsl::span<std::uint32_t>>("contributors");
    for (std::size_t v = 0; v < seeds.size(); ++v) {
      const auto& seed = seeds[v];
      const std::array<float, o2::its::Vertex::kNCov> cov{seed.varTransverse, 0.F, seed.varTransverse, 0.F, 0.F, seed.varZ};
      auto& vertex = vertices.emplace_back(o2::math_utils::Point3D<float>(static_cast<float>(shipped.beamLine.x), static_cast<float>(shipped.beamLine.y), seed.z), cov,
                                           static_cast<ushort>(std::min<std::uint32_t>(contributors[v], std::numeric_limits<ushort>::max())), 0.F);
      const auto lower = bcOfPS(seed.timePS);
      vertex.setTimeStamp(o2::its::TimeEstBC(static_cast<std::uint32_t>(lower), static_cast<std::uint32_t>(bcOfPS(static_cast<std::int64_t>(seed.timePS) + seed.widthPS) - lower)));
    }
    // By frame, as the frame records count them, then by time and contributors.
    std::ranges::stable_sort(vertices, [&](const o2::its::Vertex& a, const o2::its::Vertex& b) {
      const auto frameA = vertexFrame(a);
      const auto frameB = vertexFrame(b);
      if (frameA != frameB) {
        return frameA < frameB;
      }
      return a.getTimeStamp().lower() != b.getTimeStamp().lower() ? a.getTimeStamp().lower() < b.getTimeStamp().lower() : a.getNContributors() > b.getNContributors();
    });

    // The clock's frames of the timeframe, as ITS's (`setBCData`), and more
    // if a track or vertex lies beyond them.
    std::int64_t frameCount = (o2::constants::lhc::LHCMaxBunches / length) * mConditions->orbitsPerTimeframe();
    for (const auto& c : converted) {
      frameCount = std::max(frameCount, c.frame + 1);
    }
    for (const auto& vertex : vertices) {
      frameCount = std::max(frameCount, vertexFrame(vertex) + 1);
    }
    for (auto* frameRecords : {&trackROFs, &vertexROFs}) {
      frameRecords->resize(static_cast<std::size_t>(frameCount));
      for (std::int64_t i = 0; i < frameCount; ++i) {
        auto& rof = (*frameRecords)[static_cast<std::size_t>(i)];
        rof.setBCData(o2::InteractionRecord(static_cast<std::uint16_t>(((i * length) % o2::constants::lhc::LHCMaxBunches) + delay),
                                            orbit + static_cast<std::uint32_t>(i * length / o2::constants::lhc::LHCMaxBunches)));
        rof.setROFrame(static_cast<std::uint32_t>(i));
        rof.setNEntries(0);
        rof.setFirstEntry(-1);
      }
    }
    // As ITS: the standard vertexing mode on the input's frames, with the vertexer.
    if (recoParam<Det>().vertexer) {
      for (std::size_t i = 0; i < std::min(rofs.size(), vertexROFs.size()); ++i) {
        vertexROFs[i].setFlag(o2::itsmft::ROFRecord::VtxStdMode);
      }
    }

    std::vector<int> entries(static_cast<std::size_t>(frameCount) + 1U, 0);
    tracks.reserve(converted.size());
    for (auto& c : converted) {
      const auto& record = records[c.record];
      c.track.setFirstClusterEntry(static_cast<int>(clusterIds.size()));
      for (std::uint8_t i = record.count; i-- > 0;) { // outermost first
        clusterIds.push_back(static_cast<int>(clusterIndex[record.clusters.index[i]]));
      }
      tracks.push_back(c.track);
      ++entries[static_cast<std::size_t>(c.frame)];
    }
    std::exclusive_scan(entries.begin(), entries.end(), entries.begin(), 0);
    for (std::size_t i = 0; i < trackROFs.size(); ++i) {
      trackROFs[i].setFirstEntry(entries[i]);
      trackROFs[i].setNEntries(entries[i + 1] - entries[i]);
      trackROFs[i].setFlags(vertexROFs[i].getFlags());
      // Every frame is processed: no multiplicity selection here.
      auto& irFrame = irFrames.emplace_back(trackROFs[i].getBCData(), trackROFs[i].getBCData() + (length - 1));
      irFrame.info = static_cast<std::uint64_t>(trackROFs[i].getNEntries());
    }
    std::fill(entries.begin(), entries.end(), 0);
    for (const auto& vertex : vertices) {
      ++entries[static_cast<std::size_t>(vertexFrame(vertex))];
    }
    std::exclusive_scan(entries.begin(), entries.end(), entries.begin(), 0);
    for (std::size_t i = 0; i < vertexROFs.size(); ++i) {
      vertexROFs[i].setFirstEntry(entries[i]);
      vertexROFs[i].setNEntries(entries[i + 1] - entries[i]);
    }

    if (mUseMC) {
      auto& trackLabels = pc.outputs().make<std::vector<o2::MCCompLabel>>(Output{origin, "TRACKSMCTR", 0});
      trackLabels.reserve(converted.size());
      for (const auto& c : converted) {
        trackLabels.push_back(c.label);
      }
      // No vertex truth: the products do not say which clusters contributed.
      pc.outputs().make<std::vector<o2::MCCompLabel>>(Output{origin, "VERTICESMCTR", 0}, vertices.size());
      pc.outputs().make<std::vector<float>>(Output{origin, "VERTICESMCPUR", 0}, vertices.size(), 0.F);
    }
    LOGP(info, "{}: the timeframe at orbit {}: {} tracks ({} failed refits left out, {} in a frame turned by pi), {} vertices, {} readout frames", Det::name, orbit,
         tracks.size(), records.size() - fitted.size(), turned, vertices.size(), trackROFs.size());
  }

  template <typename LayerOf>
  void runMFT(o2::framework::ProcessingContext& pc, gsl::span<const Record> records, std::span<const std::uint32_t> fitted, std::span<const std::uint32_t> clusterIndex,
              std::span<const o2::itsmft::ROFRecord> rofs, const adaptor::ClusterFrames& frames, std::span<const std::uint32_t> pixels, const Labels* labels, double bz,
              std::uint32_t orbit, const LayerOf& layerOf)
  {
    using o2::framework::Output;
    const auto& origin = detectorOrigin<Det>();
    auto& trackROFs = pc.outputs().make<std::vector<o2::itsmft::ROFRecord>>(Output{origin, "MFTTrackROF", 0}, rofs.begin(), rofs.end());
    auto& clusterIds = pc.outputs().make<std::vector<int>>(Output{origin, "TRACKCLSID", 0});
    auto& tracks = pc.outputs().make<std::vector<o2::mft::TrackMFT>>(Output{origin, "TRACKS", 0});

    struct Converted {
      o2::mft::TrackMFT track;
      o2::MCCompLabel label;
      std::uint32_t record = 0;
      std::uint32_t frame = 0;
    };
    std::vector<Converted> converted;
    converted.reserve(fitted.size());
    std::size_t crossing = 0U;
    std::array<bool, TrackRefit<Det>::kMaxClusters> carried{};
    for (const auto track : fitted) {
      const auto& record = records[track];
      auto& out = converted.emplace_back();
      out.record = track;
      const auto inner = toTrackParCovFwd(record.state, mGeometry->planeOf(record.state), bz, record.chi2);
      static_cast<o2::track::TrackParCovFwd&>(out.track) = inner;
      out.track.setOutParam(toTrackParCovFwd(record.outwardState, mGeometry->planeOf(record.outwardState), bz, record.outwardChi2));
      out.track.setCA(true);
      out.track.setNumberOfPoints(record.count);
      if (labels != nullptr) {
        out.label = adaptor::trackLabel(record, clusterIndex, *labels, std::span<bool>{carried.data(), record.count});
      }
      for (std::uint8_t i = 0; i < record.count; ++i) {
        out.track.setClusterSize(layerOf(record.clusters.index[i]), static_cast<int>(pixels[clusterIndex[record.clusters.index[i]]]));
      }
      // TrackMFT has no time: a track across frames is timed by the earlier.
      out.frame = adaptor::earliestRecord(record, clusterIndex, frames);
      for (std::uint8_t i = 0; i < record.count; ++i) {
        if (frames.record[clusterIndex[record.clusters.index[i]]] != out.frame) {
          ++crossing;
          break;
        }
      }
    }
    std::ranges::stable_sort(converted, {}, &Converted::frame);

    std::vector<int> entries(trackROFs.size() + 1U, 0);
    tracks.reserve(converted.size());
    for (auto& c : converted) {
      const auto& record = records[c.record];
      c.track.setExternalClusterIndexOffset(static_cast<int>(clusterIds.size()));
      for (std::uint8_t i = 0; i < record.count; ++i) { // innermost first
        clusterIds.push_back(static_cast<int>(clusterIndex[record.clusters.index[i]]));
      }
      tracks.push_back(c.track);
      ++entries[c.frame];
    }
    std::exclusive_scan(entries.begin(), entries.end(), entries.begin(), 0);
    for (std::size_t i = 0; i < trackROFs.size(); ++i) {
      trackROFs[i].setFirstEntry(entries[i]);
      trackROFs[i].setNEntries(entries[i + 1] - entries[i]);
    }
    if (mUseMC) {
      auto& trackLabels = pc.outputs().make<std::vector<o2::MCCompLabel>>(Output{origin, "TRACKSMCTR", 0});
      trackLabels.reserve(converted.size());
      for (const auto& c : converted) {
        trackLabels.push_back(c.label);
      }
    }
    LOGP(info, "{}: the timeframe at orbit {}: {} tracks ({} failed refits left out, {} across readout frames, timed by the earlier), {} readout frames", Det::name,
         orbit, tracks.size(), records.size() - fitted.size(), crossing, trackROFs.size());
  }

  bool mUseMC = false;
  bool mConfigurationPublished = false;
  std::shared_ptr<Conditions<Det>> mConditions;                ///< the run's geometry, field, frames and orbits per timeframe
  std::unique_ptr<const TrackingGeometry<Det>> mGeometry;      ///< the planes the tracker's states are bound to
  const o2::itsmft::TopologyDictionary* mDictionary = nullptr; ///< the CCDB cache's: the clusters' pixel counts
};

} // namespace adaptor

template <DetectorLayout Det>
o2::framework::DataProcessorSpec getAdaptorSpec(bool useMC)
{
  using namespace o2::framework;
  constexpr bool kITS = std::same_as<Det, detectors::ITS>;
  const auto origin = detectorOrigin<Det>();
  std::vector<InputSpec> inputs{
    InputSpec{"tracks", gDataOrigin, description<Det>(kTracks), 0, Lifetime::Timeframe},
    InputSpec{"clusterIndex", gDataOrigin, description<Det>(kMeasurements), static_cast<std::uint32_t>(MeasurementColumn::clusterIndex), Lifetime::Timeframe},
    InputSpec{"rofs", origin, "CLUSTERSROF", 0, Lifetime::Timeframe},
    InputSpec{"compclusters", origin, "COMPCLUSTERS", 0, Lifetime::Timeframe},
    InputSpec{"patterns", origin, "PATTERNS", 0, Lifetime::Timeframe},
    InputSpec{"dictionary", origin, "CLUSDICT", 0, Lifetime::Condition, ccdbParamSpec(origin.template as<std::string>() + "/Calib/ClusterDictionary")}};
  std::vector<OutputSpec> outputs{{origin, "TRACKS", 0, Lifetime::Timeframe},
                                  {origin, "TRACKCLSID", 0, Lifetime::Timeframe},
                                  {"META", adaptor::metaDescription<Det>(), 0, Lifetime::Sporadic}};
  if constexpr (kITS) {
    inputs.emplace_back("measurementGeometry", gDataOrigin, description<Det>(kMeasurementGeometry), 0, Lifetime::Timeframe);
    inputs.emplace_back("vertices", gDataOrigin, description<Det>(kVertexSeeds), 0, Lifetime::Timeframe);
    inputs.emplace_back("contributors", gDataOrigin, description<Det>(kVertexContributors), 0, Lifetime::Timeframe);
    outputs.emplace_back(origin, "ITSTrackROF", 0, Lifetime::Timeframe);
    outputs.emplace_back(origin, "VERTICES", 0, Lifetime::Timeframe);
    outputs.emplace_back(origin, "VERTICESROF", 0, Lifetime::Timeframe);
    outputs.emplace_back(origin, "IRFRAMES", 0, Lifetime::Timeframe);
  } else {
    outputs.emplace_back(origin, "MFTTrackROF", 0, Lifetime::Timeframe);
  }
  if (useMC) {
    inputs.emplace_back("labels", origin, "CLUSTERSMCTR", 0, Lifetime::Timeframe);
    outputs.emplace_back(origin, "TRACKSMCTR", 0, Lifetime::Timeframe);
    if constexpr (kITS) {
      outputs.emplace_back(origin, "VERTICESMCTR", 0, Lifetime::Timeframe);
      outputs.emplace_back(origin, "VERTICESMCPUR", 0, Lifetime::Timeframe);
    }
  }
  auto conditions = std::make_shared<Conditions<Det>>(inputs);
  return DataProcessorSpec{.name = deviceName<Det>("adaptor"),
                           .inputs = inputs,
                           .outputs = outputs,
                           .algorithm = AlgorithmSpec{adaptFromTask<adaptor::AdaptorTask<Det>>(useMC, conditions)}};
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_ADAPTORDEVICE_H_

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

// The device kernels the device tests run, one entry per backend.

#ifndef GTRK_TEST_DEVICE_DEVICETESTKERNELS_H_
#define GTRK_TEST_DEVICE_DEVICETESTKERNELS_H_

#include "ToyBarrel8.h"
#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/Neighbour.h"
#include "GenericTracking/Tracking/LayerRowMap.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/Triplet.h"
#include "GenericTracking/Tracking/TripletCut.h"
#include "GenericTracking/Tracking/Trig.h"
#include "GenericTracking/Tracking/GPU/DeviceTracker.h"

#include "SyntheticStage.h"

#include "FlatObject.h"
#include "UniformMaterial.h"

#include <cstddef>
#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

namespace o2::gtrk
{

/// Every function of Trig.h at one argument pair, in both precisions: in
/// the deterministic build the host's bits (PortableMath.h).
struct TrigRow {
  double sinD = 0.;
  double cosD = 0.;
  double sincosSinD = 0.;
  double sincosCosD = 0.;
  double asinD = 0.;
  double atanD = 0.;
  double atan2D = 0.;
  float sinF = 0.F;
  float cosF = 0.F;
  float sincosSinF = 0.F;
  float sincosCosF = 0.F;
  float asinF = 0.F;
  float atanF = 0.F;
  float atan2F = 0.F;
  float logF = 0.F;
};

/// TrigRow of (y, x), on the host or, in a kernel, on the device.
GPUhdi() TrigRow trigRow(double y, double x)
{
  TrigRow row;
  row.sinD = sinOf(y);
  row.cosD = cosOf(y);
  sincosOf(y, row.sincosSinD, row.sincosCosD);
  row.asinD = asinOf(y);
  row.atanD = atanOf(y);
  row.atan2D = atan2Of(y, x);
  const auto yF = static_cast<float>(y);
  const auto xF = static_cast<float>(x);
  row.sinF = sinOf(yF);
  row.cosF = cosOf(yF);
  sincosOf(yF, row.sincosSinF, row.sincosCosF);
  row.asinF = asinOf(yF);
  row.atanF = atanOf(yF);
  row.atan2F = atan2Of(yF, xF);
  row.logF = logOf(xF);
  return row;
}

/// What a flat synthetic run reports:
/// the run's report and the elements it listed.
struct FlatRunInfo {
  FlatRunReport report;
  std::uint32_t listed = 0U;
};

/// What a synthetic GpuStage run reports besides its output.
struct SyntheticRunInfo {
  bool predictionMissed = false;
  std::size_t capacityAfter = 0U; ///< the staging capacity the stage keeps
  std::size_t scratchAfter = 0U;  ///< the arena's scratch end once run returned
  std::uint32_t listed = 0U;      ///< the elements a listed run generated
};

/// A small O2 FlatObject for the tests of DeviceFlatObject: \a count values
/// in its flat buffer, reached through a pointer the relocation must fix. ITS
/// host methods exist outside GPU translation units only, as FlatObject's do.
class ToyFlatObject : public o2::gpu::FlatObject
{
 public:
#ifndef GPUCA_GPUCODE
  void construct(std::span<const std::int32_t> values)
  {
    startConstruction();
    mCount = static_cast<std::int32_t>(values.size());
    finishConstruction(static_cast<std::int32_t>(values.size_bytes()));
    mValues = reinterpret_cast<std::int32_t*>(mFlatBufferPtr); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast): the flat buffer holds them
    std::ranges::copy(values, mValues);
  }
  void cloneFromObject(const ToyFlatObject& object, char* newFlatBufferPtr)
  {
    o2::gpu::FlatObject::cloneFromObject(object, newFlatBufferPtr);
    mCount = object.mCount;
    mValues = reinterpret_cast<std::int32_t*>(mFlatBufferPtr); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
  }
  void setFutureBufferAddress(char* futureFlatBufferPtr)
  {
    mValues = reinterpret_cast<std::int32_t*>(futureFlatBufferPtr); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    o2::gpu::FlatObject::setFutureBufferAddress(futureFlatBufferPtr);
  }
  void setActualBufferAddress(char* actualFlatBufferPtr)
  {
    mValues = reinterpret_cast<std::int32_t*>(actualFlatBufferPtr); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    o2::gpu::FlatObject::setActualBufferAddress(actualFlatBufferPtr);
  }
#endif
  [[nodiscard]] GPUhdi() std::int32_t count() const { return mCount; }
  [[nodiscard]] GPUhdi() std::int32_t value(std::int32_t index) const { return mValues[index]; }

 private:
  std::int32_t* mValues = nullptr;
  std::int32_t mCount = 0;
};

/// Defined in DeviceTestHost.cxx, a host translation unit.
template <>
struct FlatObjectPort<ToyFlatObject> {
  static std::size_t bufferBytes(const ToyFlatObject& source);
  static void image(const ToyFlatObject& source, char* futureBuffer, FlatObjectImage& image);
};

using ToyBarrel8Det = detectors::ToyBarrel8;
inline constexpr std::size_t kToyBarrel8Layers = static_cast<std::size_t>(Detector<ToyBarrel8Det>::nLayers);

/// What the road glue left, read back.
struct RoadGlueResult {
  std::vector<std::uint32_t> successorAbove;
  std::vector<std::uint32_t> roads;
  std::vector<std::uint32_t> starts;
  std::vector<std::uint32_t> offsets;
  std::uint32_t chained = 0U;
};

inline constexpr std::size_t kToyBarrel8Clusters = maxChainLayers<ToyBarrel8Det>;
using ToyBarrel8Records = BoundedColumns<SelectionRecordOf<kToyBarrel8Clusters>::RecordT>;

/// What the selection glue left, read back.
struct SelectionGlueResult {
  std::vector<std::uint32_t> ranked;
  std::vector<std::uint32_t> selected;
  std::vector<std::uint8_t> mask;
  std::uint32_t duplicates = 0U;
  std::uint32_t mismatches = 0U;
  std::uint32_t rounds = 0U;
  std::uint32_t readbacks = 0U;
  std::size_t scratchAfter = 0U; ///< the arena's scratch end once the selection returned
};

} // namespace o2::gtrk
namespace o2::base
{
class MatLayerCylSet;
}
namespace o2::gtrk
{

struct PathBoundsResult {
  std::vector<std::uint8_t> measurementPrefix;
  std::vector<std::uint8_t> measurementSuffix;
  std::vector<std::uint8_t> tripletPrefix;
  std::vector<std::uint8_t> tripletSuffix;
  std::vector<std::uint32_t> permutation;
};

struct PropagationProbe {
  Propagator prop;
  BoundTrackState before;
  TrackingSurface from;
  TrackingSurface to;
  PropagationDirection direction = PropagationDirection::kForward;
};
struct PropagationAnswer {
  BoundTrackState after;
  PropagationResult result;
};

/// The lane groups' script: per item, a group of
/// lanes takes `trips` votes and broadcasts, and on odd items nested ones
/// inside each, so that the groups of one warp take different paths and
/// some none at all. Each lane folds what the collectives give it into two
/// words; the host knows what every lane must get.
namespace lane_group_script
{
struct Payload {
  std::uint32_t a = 0U;
  std::uint32_t b = 0U;
  std::uint32_t c = 0U;
};
GPUhdi() std::uint32_t trips(std::uint32_t item) { return item % 5U; }
GPUhdi() bool vote(std::uint32_t lane, std::uint32_t trip, std::uint32_t item) { return ((lane * 7U) + (trip * 13U) + item) % 3U == 0U; }
GPUhdi() std::uint32_t source(std::uint32_t trip, std::uint32_t item, std::uint32_t lanes) { return (trip + item) % lanes; }
GPUhdi() Payload payload(std::uint32_t lane, std::uint32_t trip, std::uint32_t item) { return {.a = (lane * 1000U) + trip, .b = item, .c = (trip * 3U) + lane}; }
GPUhdi() std::uint32_t innerTrips(std::uint32_t trip, std::uint32_t item) { return item % 2U == 1U ? (item + trip) % 3U : 0U; }
GPUhdi() bool innerVote(std::uint32_t lane, std::uint32_t inner) { return lane % (inner + 2U) == 0U; }
GPUhdi() std::uint32_t innerValue(std::uint32_t lane, std::uint32_t inner) { return (lane * 31U) + inner; }
GPUhdi() std::uint32_t innerSource(std::uint32_t inner, std::uint32_t trip, std::uint32_t lanes) { return (inner + trip) % lanes; }

/// The script on \a group's lane for \a item: the folded ballots and the
/// folded broadcasts.
template <typename Group>
GPUhdi() std::array<std::uint64_t, 2> run(const Group& group, std::uint32_t item)
{
  const std::uint32_t lane = group.lane();
  std::array<std::uint64_t, 2> folded{};
  for (std::uint32_t trip = 0U; trip < trips(item); ++trip) {
    folded[0] ^= group.ballot(vote(lane, trip, item)) * (trip + 1U);
    const auto got = group.broadcast(payload(lane, trip, item), source(trip, item, Group::kLanes));
    folded[1] += got.a + (7U * got.b) + (13U * got.c);
    for (std::uint32_t inner = 0U; inner < innerTrips(trip, item); ++inner) {
      folded[0] ^= group.ballot(innerVote(lane, inner)) << (inner + 1U);
      folded[1] += group.broadcast(innerValue(lane, inner), innerSource(inner, trip, Group::kLanes));
    }
  }
  return folded;
}
} // namespace lane_group_script

struct DeviceTestEntry {
  /// trigRow of every (ys[i], xs[i]) on the device.
  void (*trig)(DeviceRuntime&, DeviceStream&, const double* ys, const double* xs, TrigRow* rows, std::uint32_t count) = nullptr;
  /// GpuStage<TrackletT> over the synthetic generator, its output prepared in
  /// the caller's scope of \a arena and, if complete, downloaded to \a result.
  StageStatus (*runSynthetic)(DeviceRuntime&, DeviceStream&, DeviceArena&, synthetic::Kind, std::uint32_t elements, std::size_t stagingCapacity, StageOutput<TrackletT>& result, SyntheticRunInfo& info) = nullptr;
  /// The same over a list: every element that emits and every third of the
  /// rest (`GpuStage::run` with a list).
  StageStatus (*runSyntheticListed)(DeviceRuntime&, DeviceStream&, DeviceArena&, synthetic::Kind, std::uint32_t elements, std::size_t stagingCapacity, StageOutput<TrackletT>& result, SyntheticRunInfo& info) = nullptr;
  /// The flat form (`flat::runFlat`) over the
  /// synthetic candidates of \a kind with `synthetic::Accept{divisor}`: all
  /// the elements, or the listed ones (as `runSyntheticListed` lists them);
  /// chunks of at most \a maxChunkPairs pairs (0: the arena decides).
  StageStatus (*runFlatSynthetic)(DeviceRuntime&, DeviceStream&, DeviceArena&, synthetic::Kind, std::uint32_t divisor, std::uint32_t elements, bool listedRun, std::size_t maxChunkPairs, StageOutput<TrackletT>& result, FlatRunInfo& info) = nullptr;
  /// The device tracker for the ToyBarrel8 fixtures.
  std::unique_ptr<DeviceTracker<detectors::ToyBarrel8, UniformMaterial>> (*makeToyBarrel8Tracker)(DeviceRuntime&, DeviceBlock arena) = nullptr;
  /// The slice's glue passes on host inputs copied to the device, their
  /// results copied back.
  StageStatus (*permutation)(DeviceRuntime&, DeviceStream&, DeviceArena&, const IterationParameters<ToyBarrel8Det>&, const LayerRowMap<ToyBarrel8Det>&, const StageOutput<TrackletT>& tracklets, std::span<const std::uint8_t> prefix, std::span<const std::uint8_t> suffix, std::vector<std::uint32_t>& result) = nullptr;
  StageStatus (*pathBounds)(DeviceRuntime&, DeviceStream&, DeviceArena&, const IterationParameters<ToyBarrel8Det>&, const LayerRowMap<ToyBarrel8Det>&, const StageOutput<TrackletT>&, const StageOutput<TripletT>&, std::span<const std::uint32_t>, PathBoundsResult&) = nullptr;
  StageStatus (*endIndex)(DeviceRuntime&, DeviceStream&, DeviceArena&, const StageOutput<TripletT>& triplets, std::uint32_t trackletCount, std::vector<std::uint32_t>& offsets, std::vector<std::uint32_t>& indices) = nullptr;
  StageStatus (*chained)(DeviceRuntime&, DeviceStream&, DeviceArena&, std::span<const std::uint32_t> neighbourOffsets, std::vector<std::uint32_t>& result) = nullptr;
  /// The shared device algorithms and storage on host inputs copied to the
  /// device.
  StageStatus (*selectIndices)(DeviceRuntime&, DeviceStream&, DeviceArena&, std::span<const std::uint8_t> mask, std::vector<std::uint32_t>& result) = nullptr;
  StageStatus (*scan32)(DeviceRuntime&, DeviceStream&, DeviceArena&, std::span<const std::uint32_t> input, std::vector<std::uint32_t>& result) = nullptr;
  StageStatus (*scan64)(DeviceRuntime&, DeviceStream&, DeviceArena&, std::span<const std::uint64_t> input, std::vector<std::uint64_t>& result) = nullptr;
  /// A chain of \a length values, the first 1: each pass raises value i to
  /// value i-1 plus one, raising the change flag; a value above \a threshold
  /// raises the second, sticky flag.
  StageStatus (*sweepChain)(DeviceRuntime&, DeviceStream&, DeviceArena&, std::uint32_t length, std::uint32_t threshold, std::uint32_t maxChangingSweeps, std::vector<std::uint32_t>& values, std::array<std::uint32_t, 2>& flags, std::uint32_t& changingSweeps) = nullptr;
  /// \a rows tracklets of \a source from \a sourceOffset copied into a device
  /// destination of \a destinationRows rows, all bits set, at
  /// \a destinationOffset; the destination read back.
  StageStatus (*copyColumns)(DeviceRuntime&, DeviceStream&, DeviceArena&, const StageOutput<TrackletT>& source, std::size_t sourceOffset, std::size_t destinationRows, std::size_t destinationOffset, std::size_t rows, std::vector<Tracklet>& destination) = nullptr;
  /// \a source cut into chunks of \a chunkRows rows, appended on the scratch
  /// end in one scope, gathered on the base end and read back; false if the
  /// arena refused a chunk or the gather.
  bool (*gatherChunks)(DeviceRuntime&, DeviceStream&, DeviceArena&, const StageOutput<TrackletT>& source, std::span<const std::size_t> chunkRows, std::vector<Tracklet>& gathered) = nullptr;
  /// \a object ported to the device and its values read there; false if the
  /// arena refused it.
  /// \a object uploaded, read back; uploaded again elsewhere from the kept image, read back into \a moved.
  bool (*flatObject)(DeviceRuntime&, DeviceStream&, DeviceArena&, const ToyFlatObject& object, std::vector<std::int32_t>& readBack, std::vector<std::int32_t>& moved) = nullptr;
  /// O2's material LUT ported with DeviceFlatObject, and getMatBudget of each
  /// segment (x0, y0, z0, x1, y1, z1) evaluated on the device: its mean
  /// density, mean X/X0 and length.
  StageStatus (*matBudgets)(DeviceRuntime&, DeviceStream&, DeviceArena&, const o2::base::MatLayerCylSet& lut, std::span<const std::array<float, 6>> segments, std::vector<std::array<float, 3>>& result) = nullptr;
  /// The selection glue over \a records into a copy of \a mask, in an arena
  /// scope closed on return.
  StageStatus (*selectionGlue)(DeviceRuntime&, DeviceStream&, DeviceArena&, const ToyBarrel8Records& records, std::span<const std::uint8_t> mask, std::uint32_t sharedInnermost, SelectionGlueResult& result) = nullptr;
  /// Exact helix propagation on the device for independent closure validation.
  StageStatus (*propagation)(DeviceRuntime&, DeviceStream&, DeviceArena&, std::span<const PropagationProbe>, std::vector<PropagationAnswer>&) = nullptr;
  /// `lane_group_script::run` for \a items items in groups of \a lanes lanes
  /// (`forEachGroup`, `WarpGroup`): two words per lane, item by item.
  StageStatus (*laneGroups)(DeviceRuntime&, DeviceStream&, DeviceArena&, std::uint32_t lanes, std::uint32_t items, std::vector<std::uint64_t>&) = nullptr;
  /// A slice's road counts, independently scheduled ordered counts, then its
  /// starts and offsets, on host inputs copied to the device: the triplets
  /// with their final levels, the neighbours and the core rows per layer.
  StageStatus (*orderedRoadGlue)(DeviceRuntime&, DeviceStream&, DeviceArena&, const IterationParameters<ToyBarrel8Det>&, const StageOutput<TrackletT>& tracklets, const StageOutput<TripletT>& triplets, const StageOutput<NeighbourT>& neighbours, const std::array<MeasurementRange, kToyBarrel8Layers>& core, RoadGlueResult& result) = nullptr;
};

const DeviceTestEntry* testKernelsCUDA();
const DeviceTestEntry* testKernelsHIP();

} // namespace o2::gtrk

#endif // GTRK_TEST_DEVICE_DEVICETESTKERNELS_H_

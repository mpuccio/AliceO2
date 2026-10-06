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

#include <cuda_runtime.h>
#if defined(__HIPCC__)
#include <hipcub/hipcub.hpp>
namespace tracking_cub = hipcub;
#else
#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_scan.cuh>
#include <cub/device/device_select.cuh>
namespace tracking_cub = cub;
#endif
#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstring>
#include <exception>
#include <limits>
#include <new>
#include <string>
#include <thread>
#include <utility>
#include "ITSMFTTrackingGPU/TimeFrameGPU.h"

namespace o2::itsmft::tracking::gpu
{
inline namespace ITSMFT_TRACKING_GPU_BACKEND
{
namespace
{
void check(cudaError_t error)
{
  if (error != cudaSuccess) {
    throw std::runtime_error{std::string{"GPU tracking: "} + cudaGetErrorString(error)};
  }
}
cudaStream_t toStream(void* stream) { return static_cast<cudaStream_t>(stream); }
cudaEvent_t toEvent(void* event) { return static_cast<cudaEvent_t>(event); }

template <typename F>
__global__ void forEachIndex(size_t n, F f)
{
  for (size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x; i < n; i += size_t(gridDim.x) * blockDim.x) {
    f(i);
  }
}
template <typename F>
void launch(cudaStream_t stream, size_t n, F f, unsigned threads = 128)
{
  if (n) {
    const auto blocks = static_cast<unsigned>(std::min<size_t>((n + threads - 1) / threads, 65535));
    forEachIndex<<<blocks, threads, 0, stream>>>(n, f);
    check(cudaGetLastError());
  }
}

// Stores every item emitted up to capacity and counts them all, as the
// legacy kernels do; the order is fixed afterwards by sorting.
template <typename T>
struct Appender {
  T* items;
  int capacity;
  int* counter;
  GPUdi() void operator()(const T& item) const
  {
    const int slot = atomicAdd(counter, 1);
    if (slot < capacity) {
      items[slot] = item;
    }
  }
};

struct RoadRange {
  int path, begin, end;
};
// A road job built on the device: its seed is initial[seed] (initialize),
// the previous extension's result of job, or seeds[seed].
struct RoadDeviceJob {
  size_t job, firstTarget, endTarget, seed;
  bool initialize;
};
template <typename Emit>
GPUdi() void runRoadJob(const RoadDeviceJob& job, const SeedInput* initial, const TrackSeed* seeds, const RoadResult* previous,
                        const RoadTarget* targets, float bz, float maxChi2, Emit emit)
{
  const TrackSeed seed = job.initialize ? TrackSeed{} : previous ? previous[job.job].emission.seed : seeds[job.seed];
  forEachRoad(seed, job.initialize ? initial + job.seed : nullptr, targets, job.firstTarget, job.endTarget, bz, maxChi2, emit);
}

void requireLoaded(const TimeFrameGPU& frame)
{
  if (!frame.isLoaded()) {
    throw std::logic_error{"GPU tracking: the TimeFrameGPU is not loaded"};
  }
}
} // namespace

struct Workspace {
  static constexpr unsigned MaxSlots = 6;
  // A stream and the scratch buffers it hands out in order, reused by every
  // handler call on it.
  struct Slot {
    cudaStream_t stream{};
    std::vector<DeviceBuffer> buffers;
    size_t next{};
  };
  Workspace(int deviceId, const DeviceMemory* deviceMemory) : device{deviceId}, memory{deviceMemory} {}
  ~Workspace()
  {
    int previous = device;
    cudaGetDevice(&previous);
    cudaSetDevice(device);
    for (auto& slot : slots) {
      cudaStreamSynchronize(slot.stream);
      slot.buffers.clear();
      cudaStreamDestroy(slot.stream);
    }
    cudaSetDevice(previous);
  }
  int device;
  const DeviceMemory* memory;
  std::once_flag created; // on first use: constructing traits does no GPU work
  std::vector<Slot> slots;
  std::vector<bool> free;
  std::mutex mutex;
  std::condition_variable released;
  std::atomic<size_t> allocations{0}, reservedBytes{0}, uploadedBytes{0};
};

namespace
{
// One slot of the workspace for one handler call. Handlers run on distinct
// slots concurrently; extra callers wait for a free one.
struct Lease {
  explicit Lease(TimeFrameGPU& frame) : w{frame.workspace()}
  {
    check(cudaSetDevice(w.device)); // TBB workers have independent current devices
    std::call_once(w.created, [this] {
      w.slots.resize(std::clamp(std::thread::hardware_concurrency(), 1u, Workspace::MaxSlots));
      w.free.assign(w.slots.size(), true);
      for (auto& slot : w.slots) {
        check(cudaStreamCreateWithFlags(&slot.stream, cudaStreamNonBlocking));
      }
    });
    std::unique_lock lock{w.mutex};
    w.released.wait(lock, [this] { return std::find(w.free.begin(), w.free.end(), true) != w.free.end(); });
    index = std::find(w.free.begin(), w.free.end(), true) - w.free.begin();
    w.free[index] = false;
    stream = w.slots[index].stream;
    w.slots[index].next = 0;
  }
  ~Lease()
  {
    if (std::uncaught_exceptions() > exceptions) {
      cudaStreamSynchronize(stream); // before host staging unwinds
    }
    {
      std::lock_guard lock{w.mutex};
      w.free[index] = true;
    }
    w.released.notify_one();
  }
  Lease(const Lease&) = delete;
  Lease& operator=(const Lease&) = delete;

  template <typename T>
  T* take(size_t n)
  {
    auto& slot = w.slots[index];
    if (slot.next == slot.buffers.size()) {
      slot.buffers.emplace_back(w.memory, DeviceMemory::Scope::Pass);
    }
    auto& buffer = slot.buffers[slot.next++];
    const size_t before = buffer.capacity;
    if (buffer.reserve(std::max<size_t>(n, 1) * sizeof(T))) {
      ++w.allocations;
      w.reservedBytes += buffer.capacity - before;
    }
    return static_cast<T*>(buffer.data);
  }
  template <typename T>
  T* upload(const T* from, size_t n)
  {
    auto* to = take<T>(n);
    if (n) {
      check(cudaMemcpyAsync(to, from, n * sizeof(T), cudaMemcpyHostToDevice, stream));
      w.uploadedBytes += n * sizeof(T);
    }
    return to;
  }
  template <typename T>
  T get(const T* from)
  {
    T value;
    check(cudaMemcpyAsync(&value, from, sizeof(T), cudaMemcpyDeviceToHost, stream));
    sync();
    return value;
  }
  void zero(void* data, size_t bytes) { check(cudaMemsetAsync(data, 0, bytes, stream)); }
  void sync() { check(cudaStreamSynchronize(stream)); }
  template <typename F>
  void launch(size_t n, F f, unsigned threads = 128)
  {
    gpu::launch(stream, n, f, threads);
  }
  // Runs a cub algorithm with scratch of the size it asks for.
  template <typename F>
  void cub(F algorithm)
  {
    size_t bytes = 0;
    check(algorithm(nullptr, bytes));
    check(algorithm(take<unsigned char>(bytes), bytes));
  }
  static int items(size_t n)
  {
    if (n >= static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw std::length_error{"GPU tracking: too many items"};
    }
    return static_cast<int>(n);
  }
  // Exclusive scan of counts[0, n) into offsets[0, n]: offsets[n] is the
  // total. counts holds n + 1 elements.
  template <typename T>
  void scan(T* counts, T* offsets, size_t n)
  {
    zero(counts + n, sizeof(T));
    cub([&](void* scratch, size_t& bytes) { return tracking_cub::DeviceScan::ExclusiveSum(scratch, bytes, counts, offsets, items(n + 1), stream); });
  }
  // Stably sorts items[0, n) by the low bits of key(item); returns the
  // sorted keys.
  template <typename T, typename Key>
  uint64_t* sortByKey(T* data, size_t n, Key key, int bits)
  {
    auto* keys = take<uint64_t>(n);
    auto* sortedKeys = take<uint64_t>(n);
    auto* order = take<int>(n);
    auto* sortedOrder = take<int>(n);
    auto* sorted = take<T>(n);
    launch(n, [=] __device__(size_t i) {
      keys[i] = key(data[i]);
      order[i] = static_cast<int>(i);
    });
    cub([&](void* scratch, size_t& bytes) {
      return tracking_cub::DeviceRadixSort::SortPairs(scratch, bytes, keys, sortedKeys, order, sortedOrder, items(n), 0, bits, stream);
    });
    launch(n, [=] __device__(size_t i) { sorted[i] = data[sortedOrder[i]]; });
    check(cudaMemcpyAsync(data, sorted, n * sizeof(T), cudaMemcpyDeviceToDevice, stream));
    return sortedKeys;
  }
  // The lookup table [0, nKeys] of items sorted by key(item): the first item
  // of every key, plus the total.
  template <typename T, typename Key>
  void buildLookup(const T* data, size_t n, int* lookup, size_t nKeys, Key key)
  {
    auto* counts = take<int>(nKeys + 1);
    zero(counts, (nKeys + 1) * sizeof(int));
    launch(n, [=] __device__(size_t i) { atomicAdd(counts + key(data[i]), 1); });
    scan(counts, lookup, nKeys);
  }

  Workspace& w;
  size_t index;
  cudaStream_t stream;
  int exceptions = std::uncaught_exceptions();
};

// Evaluates road jobs into output, as the legacy processNeighboursHandler;
// returns the number of results.
size_t runRoads(Lease& lease, DeviceVector<RoadResult>& output, const RoadDeviceJob* jobs, size_t nJobs, const SeedInput* initial,
                const TrackSeed* seeds, const RoadResult* previous, const RoadTarget* targets, float bz, float maxChi2)
{
  auto* counts = lease.take<size_t>(nJobs + 1);
  auto* offsets = lease.take<size_t>(nJobs + 1);
  lease.launch(nJobs, [=] __device__(size_t i) {
    size_t count = 0;
    runRoadJob(jobs[i], initial, seeds, previous, targets, bz, maxChi2, [&](const RoadSeedEmission&) { ++count; });
    counts[i] = count;
  });
  lease.scan(counts, offsets, nJobs);
  const size_t total = lease.get(offsets + nJobs);
  auto* results = output.reserve(total);
  lease.launch(total ? nJobs : 0, [=] __device__(size_t i) {
    size_t slot = offsets[i];
    runRoadJob(jobs[i], initial, seeds, previous, targets, bz, maxChi2, [&](const RoadSeedEmission& emission) { results[slot++] = {jobs[i].job, emission}; });
  });
  lease.sync();
  output.size = total;
  return total;
}
} // namespace

int currentDevice()
{
  int device;
  check(cudaGetDevice(&device));
  return device;
}

DeviceBuffer::DeviceBuffer(DeviceBuffer&& other) noexcept
  : data{std::exchange(other.data, nullptr)}, capacity{std::exchange(other.capacity, 0)}, mMemory{other.mMemory}, mScope{other.mScope}, mEpoch{other.mEpoch}, mOwned{std::exchange(other.mOwned, false)} {}

DeviceBuffer& DeviceBuffer::operator=(DeviceBuffer&& other) noexcept
{
  std::swap(data, other.data);
  std::swap(capacity, other.capacity);
  std::swap(mMemory, other.mMemory);
  std::swap(mScope, other.mScope);
  std::swap(mEpoch, other.mEpoch);
  std::swap(mOwned, other.mOwned);
  return *this;
}

DeviceBuffer::~DeviceBuffer() { free(); }

void DeviceBuffer::free() noexcept
{
  if (mOwned) {
    cudaFree(data);
  }
  data = nullptr;
  capacity = 0;
  mOwned = false;
}

bool DeviceBuffer::reserve(size_t bytes)
{
  if (released()) {
    free();
  }
  if (bytes <= capacity) {
    return false;
  }
  size_t next = std::max<size_t>(capacity, 4096);
  while (next < bytes) {
    next = next > std::numeric_limits<size_t>::max() / 2 ? bytes : next * 2;
  }
  void* replacement{};
  auto* external = mMemory ? mMemory->getFrameworkAllocator() : nullptr;
  if (external) {
    // What the buffer held stays in the pool until the framework releases it.
    replacement = external->allocateDevice(next, mScope == DeviceMemory::Scope::Pass);
    if (!replacement) {
      throw std::bad_alloc{};
    }
  } else {
    check(cudaMalloc(&replacement, next));
  }
  free();
  data = replacement;
  capacity = next;
  mOwned = !external;
  if (mMemory) {
    mEpoch = mMemory->epoch(mScope);
  }
  return true;
}

namespace
{
// The legacy tracker's stack tags: "ITSITER" and the iteration's digit.
uint64_t iterationTag(int iteration)
{
  const char tag[] = {'I', 'T', 'S', 'I', 'T', 'E', 'R', static_cast<char>('0' + iteration % 10)};
  uint64_t value = 0;
  for (size_t i = 0; i < sizeof(tag); ++i) {
    value |= static_cast<uint64_t>(static_cast<unsigned char>(tag[i])) << (i * 8);
  }
  return value;
}
} // namespace

TimeFrameGPU::TimeFrameGPU(int device) : mDevice{device}, mWorkspace{std::make_unique<Workspace>(device, &mMemory)} {}

void TimeFrameGPU::pushMemoryStack(int iteration)
{
  if (auto* external = mMemory.getFrameworkAllocator()) {
    external->pushTagOnStack(iterationTag(iteration));
  }
}

void TimeFrameGPU::popMemoryStack(int iteration)
{
  if (auto* external = mMemory.getFrameworkAllocator()) {
    // Nothing may still be running on what the pass allocated.
    check(cudaSetDevice(mDevice));
    check(cudaDeviceSynchronize());
    mMemory.release(DeviceMemory::Scope::Pass);
    external->popTagOffStack(iterationTag(iteration));
  }
}

TimeFrameGPU::~TimeFrameGPU()
{
  mWorkspace.reset(); // drains the handlers' streams
  int previous = mDevice;
  cudaGetDevice(&previous);
  cudaSetDevice(mDevice);
  for (auto& layer : mLayers) {
    if (layer.stream) {
      cudaStreamSynchronize(toStream(layer.stream));
      cudaStreamDestroy(toStream(layer.stream));
      cudaEventDestroy(toEvent(layer.ready));
    }
  }
  if (mStream) {
    cudaStreamSynchronize(toStream(mStream));
    cudaStreamDestroy(toStream(mStream));
  }
  cudaSetDevice(previous);
}

WorkspaceStatistics TimeFrameGPU::statistics() const
{
  return {mWorkspace->allocations, mWorkspace->reservedBytes, mWorkspace->uploadedBytes};
}

void* TimeFrameGPU::stream() const
{
  check(cudaSetDevice(mDevice));
  std::call_once(mStreamCreated, [this] {
    cudaStream_t stream;
    check(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    mStream = stream;
  });
  return mStream;
}

void TimeFrameGPU::copy(void* to, const void* from, size_t bytes, bool toDevice) const
{
  if (!bytes) {
    return;
  }
  auto* copyStream = toStream(stream());
  check(cudaMemcpyAsync(to, from, bytes, toDevice ? cudaMemcpyHostToDevice : cudaMemcpyDeviceToHost, copyStream));
  check(cudaStreamSynchronize(copyStream));
  if (toDevice) {
    mWorkspace->uploadedBytes += bytes;
  }
}

void TimeFrameGPU::uploadAsync(void* stream, void* to, const void* from, size_t bytes)
{
  if (bytes) {
    check(cudaMemcpyAsync(to, from, bytes, cudaMemcpyHostToDevice, toStream(stream)));
    mWorkspace->uploadedBytes += bytes;
  }
}

void TimeFrameGPU::publishView()
{
  mDeviceView.reserve(1);
  copy(mDeviceView.data(), &mView, sizeof(FrameView), true);
}

void TimeFrameGPU::initialise(size_t nLayers)
{
  if (nLayers > MaxLayoutSurfaces) {
    throw std::out_of_range{"GPU timeframe layers"};
  }
  check(cudaSetDevice(mDevice));
  mLoaded = false;
  if (mMemory.hasFrameworkAllocator()) { // the framework cleared the previous timeframe's memory
    mMemory.release(DeviceMemory::Scope::TimeFrame);
    mMemory.release(DeviceMemory::Scope::Pass);
  }
  while (mLayers.size() < nLayers) {
    mLayers.emplace_back(&mMemory);
  }
  for (auto& layer : mLayers) {
    if (!layer.stream) {
      cudaStream_t stream;
      check(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
      layer.stream = stream;
      cudaEvent_t event;
      check(cudaEventCreateWithFlags(&event, cudaEventDisableTiming));
      layer.ready = event;
    }
  }
  mView = {};
  mView.nLayers = nLayers;
}

void TimeFrameGPU::finishLoading()
{
  publishView();
  mLoaded = true;
}

void TimeFrameGPU::loadLayer(int index, const LayerHostData& host)
{
  if (index < 0 || static_cast<size_t>(index) >= mView.nLayers) {
    throw std::out_of_range{"GPU timeframe layer"};
  }
  if (host.rofClusters.size() == 1 || (!host.rofClusters.empty() && host.rofClusters.back() != static_cast<int>(host.clusters.size()))) {
    throw std::invalid_argument{"GPU timeframe ROF boundaries"};
  }
  const bool tables = !host.tableBuilt.empty();
  const int nROFs = host.rofClusters.empty() ? 0 : static_cast<int>(host.rofClusters.size()) - 1;
  if (tables && (host.tableBuilt.size() != static_cast<size_t>(nROFs) || host.clusterBins.size() != host.clusters.size() || !host.tableSize)) {
    throw std::invalid_argument{"GPU timeframe index table dimensions"};
  }
  check(cudaSetDevice(mDevice));
  auto& layer = mLayers[index];
  auto* stream = layer.stream;
  layer.clusters.reserve(host.clusters.size_bytes());
  layer.rofClusters.reserve(host.rofClusters.size_bytes());
  layer.measurements.reserve(host.surfaceMeasurements.size_bytes());
  uploadAsync(stream, layer.clusters.data, host.clusters.data(), host.clusters.size_bytes());
  uploadAsync(stream, layer.rofClusters.data, host.rofClusters.data(), host.rofClusters.size_bytes());
  uploadAsync(stream, layer.measurements.data, host.surfaceMeasurements.data(), host.surfaceMeasurements.size_bytes());
  auto& view = mView.layers[index];
  view.clusters = static_cast<const GlobalMeasurement*>(layer.clusters.data);
  view.nClusters = host.clusters.size();
  view.measurements = static_cast<const SurfaceMeasurement*>(layer.measurements.data);
  view.nMeasurements = host.surfaceMeasurements.size();
  view.rofClusters = static_cast<const int*>(layer.rofClusters.data);
  view.nROFs = nROFs;
  view.indexTables = nullptr;
  view.tableSize = 0;
  if (tables) {
    // As TimeFrame::prepareClusters: within a built ROF clusters are sorted
    // by bin, so entry b is the number of clusters with bin < b; tables of
    // ROFs that were not built are zero.
    const size_t stride = host.tableSize;
    layer.bins.reserve(host.clusterBins.size_bytes());
    layer.built.reserve(host.tableBuilt.size_bytes());
    layer.indexTables.reserve(nROFs * stride * sizeof(int));
    uploadAsync(stream, layer.bins.data, host.clusterBins.data(), host.clusterBins.size_bytes());
    uploadAsync(stream, layer.built.data, host.tableBuilt.data(), host.tableBuilt.size_bytes());
    auto* table = static_cast<int*>(layer.indexTables.data);
    const auto* bins = static_cast<const int*>(layer.bins.data);
    const auto* built = static_cast<const uint8_t*>(layer.built.data);
    const auto* boundaries = view.rofClusters;
    launch(toStream(stream), nROFs * stride, [=] __device__(size_t i) {
      const size_t rof = i / stride;
      if (!built[rof]) {
        table[i] = 0;
        return;
      }
      const int bin = static_cast<int>(i - rof * stride);
      int lo = boundaries[rof], hi = boundaries[rof + 1];
      while (lo < hi) {
        const int mid = lo + (hi - lo) / 2;
        if (bins[mid] < bin) {
          lo = mid + 1;
        } else {
          hi = mid;
        }
      }
      table[i] = lo - boundaries[rof];
    });
    view.indexTables = table;
    view.tableSize = stride;
  }
  check(cudaEventRecord(toEvent(layer.ready), toStream(stream)));
  publishView();
}

void TimeFrameGPU::waitLayer(void* stream, int layer) const
{
  check(cudaStreamWaitEvent(toStream(stream), toEvent(mLayers[layer].ready), 0));
}

void TimeFrameGPU::loadLayers(DeviceBuffer& buffer, std::span<const std::span<const uint8_t>> layers, const uint8_t* LayerView::*flags)
{
  if (layers.size() > mView.nLayers) {
    throw std::out_of_range{"GPU timeframe flag layers"};
  }
  size_t total = 0;
  for (const auto& layer : layers) {
    total += layer.size();
  }
  check(cudaSetDevice(mDevice));
  buffer.reserve(total);
  auto* data = static_cast<uint8_t*>(buffer.data);
  for (size_t layer = 0; layer < mView.nLayers; ++layer) {
    const auto flagsOfLayer = layer < layers.size() ? layers[layer] : std::span<const uint8_t>{};
    uploadAsync(stream(), data, flagsOfLayer.data(), flagsOfLayer.size());
    mView.layers[layer].*flags = flagsOfLayer.empty() ? nullptr : data;
    if (flags == &LayerView::used) {
      mView.layers[layer].nUsed = flagsOfLayer.size();
    } else if (!flagsOfLayer.empty() && flagsOfLayer.size() != static_cast<size_t>(mView.layers[layer].nROFs)) {
      throw std::invalid_argument{"GPU timeframe ROF mask dimensions"};
    }
    data += flagsOfLayer.size();
  }
  publishView(); // synchronises the uploads
}

void TimeFrameGPU::loadUsedClusters(std::span<const std::span<const uint8_t>> used) { loadLayers(mUsed, used, &LayerView::used); }

void TimeFrameGPU::loadROFEnabled(std::span<const std::span<const uint8_t>> enabled) { loadLayers(mEnabled, enabled, &LayerView::rofEnabled); }

void TimeFrameGPU::loadVertices(std::span<const Vertex> vertices)
{
  check(cudaSetDevice(mDevice));
  mVertices.reserve(vertices.size_bytes());
  uploadAsync(stream(), mVertices.data, vertices.data(), vertices.size_bytes());
  mView.vertices = static_cast<const Vertex*>(mVertices.data);
  publishView();
}

void TimeFrameGPU::loadTimeFrame(TimeFrame& frame, size_t nTableLayers)
{
  const size_t nLayers = frame.getNMeasurementSurfaces();
  initialise(nLayers);
  for (size_t index = 0; index < nLayers; ++index) {
    const int layer = static_cast<int>(index);
    LayerHostData host;
    const auto& clusters = frame.getClusters()[index];
    host.clusters = {clusters.data(), clusters.size()};
    const auto boundaries = frame.getROFrameClusters(layer);
    host.rofClusters = {boundaries.data(), boundaries.size()};
    // Layers whose clusters were never sorted have no tables; the tracklet
    // stage refuses to search them.
    if (index < nTableLayers) {
      const auto bins = frame.getClusterBins(layer);
      const auto built = frame.getIndexTableBuilt(layer);
      if (bins.size() == clusters.size() && built.size() + 1 == boundaries.size()) {
        const auto& utils = frame.getIndexTableUtils(layer);
        host.clusterBins = {bins.data(), bins.size()};
        host.tableBuilt = {built.data(), built.size()};
        host.tableSize = static_cast<size_t>(utils.getNrowBins()) * utils.getNcolBins() + 1;
      }
    }
    const auto measurements = frame.getSurfaceMeasurements(LayerId{static_cast<uint16_t>(index)});
    host.surfaceMeasurements = {measurements.data(), measurements.size()};
    loadLayer(layer, host);
  }
  loadUsedClusters(frame);
  const auto& vertices = frame.getPrimaryVertices();
  loadVertices({vertices.data(), vertices.size()});
  finishLoading();
}

void TimeFrameGPU::loadUsedClusters(TimeFrame& frame)
{
  std::vector<std::span<const uint8_t>> used(frame.getNMeasurementSurfaces());
  for (size_t layer = 0; layer < used.size(); ++layer) {
    const auto flags = frame.getUsedClusters(static_cast<int>(layer));
    used[layer] = {flags.data(), flags.size()};
  }
  loadUsedClusters(used);
}

void TimeFrameGPU::loadROFEnabled(const TimeFrame& frame)
{
  const size_t nLayers = frame.getNMeasurementSurfaces();
  std::vector<std::vector<uint8_t>> masks(nLayers);
  std::vector<std::span<const uint8_t>> spans(nLayers);
  for (size_t layer = 0; layer < nLayers; ++layer) {
    masks[layer].resize(frame.getNrof(static_cast<int>(layer)));
    for (size_t rof = 0; rof < masks[layer].size(); ++rof) {
      masks[layer][rof] = frame.isROFEnabled(static_cast<int>(layer), static_cast<int>(rof));
    }
    spans[layer] = masks[layer];
  }
  loadROFEnabled(spans);
}

int computeTracklets(TimeFrameGPU& frame, const TrackletSearch& hostSearch, int edge)
{
  auto search = hostSearch;
  Lease lease{frame};
  requireLoaded(frame);
  const auto& view = frame.view();
  const int from = search.edgeCache.fromLayer;
  const int to = search.edgeCache.toLayer;
  if (from < 0 || to < 0 || size_t(std::max(from, to)) >= view.nLayers) {
    throw std::out_of_range{"GPU tracklet layers"};
  }
  const auto& source = view.layers[from];
  const auto& target = view.layers[to];
  const size_t nSources = source.nClusters;
  if (nSources && (!source.nROFs || !source.rofEnabled || !target.rofEnabled || search.nOverlapRanges != static_cast<size_t>(source.nROFs) || !target.tableSize ||
                   (!search.useDiamond && search.nVertexRanges != static_cast<size_t>(source.nROFs)))) {
    throw std::out_of_range{"GPU tracklet navigation table dimensions"};
  }
  auto& results = frame.tracklets(edge);
  results.items.size = 0;
  int* lookup = results.lookup.reserve(Lease::items(nSources + 1));
  results.lookup.size = nSources + 1;
  lease.zero(lookup, (nSources + 1) * sizeof(int));
  if (!nSources) {
    lease.sync();
    return 0;
  }
  frame.waitLayer(lease.stream, from);
  frame.waitLayer(lease.stream, to);
  search.overlapRanges = lease.upload(search.overlapRanges, search.nOverlapRanges);
  if (!search.useDiamond) {
    search.vertexRanges = lease.upload(search.vertexRanges, search.nVertexRanges);
  }
  auto* counter = lease.take<int>(1);
  lease.zero(counter, sizeof(int));
  const int capacity = static_cast<int>(std::min<size_t>(results.capacity, std::numeric_limits<int>::max()));
  const Appender<Tracklet> append{results.items.data(), capacity, counter};
  const auto* deviceView = frame.deviceView();
  lease.launch(nSources, [=] __device__(size_t i) { forEachTracklet(*deviceView, search, static_cast<int>(i), append); });
  const int emitted = lease.get(counter);
  if (emitted > capacity || !emitted) {
    return emitted;
  }
  // The same pair can be found from several vertices.
  auto* tracklets = results.items.data();
  const auto* keys = lease.sortByKey(tracklets, emitted, [] __device__(const Tracklet& tracklet) {
    return (uint64_t(uint32_t(tracklet.firstClusterIndex)) << 32) | uint32_t(tracklet.secondClusterIndex);
  }, 64);
  auto* uniqueKeys = lease.take<uint64_t>(emitted);
  auto* unique = lease.take<Tracklet>(emitted);
  auto* nUnique = lease.take<int>(1);
  lease.cub([&](void* scratch, size_t& bytes) {
    return tracking_cub::DeviceSelect::UniqueByKey(scratch, bytes, keys, tracklets, uniqueKeys, unique, nUnique, emitted, lease.stream);
  });
  const int n = lease.get(nUnique);
  check(cudaMemcpyAsync(tracklets, unique, n * sizeof(Tracklet), cudaMemcpyDeviceToDevice, lease.stream));
  lease.buildLookup(tracklets, n, lookup, nSources, [] __device__(const Tracklet& tracklet) { return tracklet.firstClusterIndex; });
  lease.sync();
  results.items.size = n;
  return emitted;
}

int computeCells(TimeFrameGPU& frame, const CellSearch& hostSearch, int firstEdge, int secondEdge, int path)
{
  auto search = hostSearch;
  Lease lease{frame};
  requireLoaded(frame);
  const auto& first = frame.tracklets(firstEdge);
  const auto& second = frame.tracklets(secondEdge);
  auto& results = frame.cells(path);
  results.items.size = 0;
  const size_t nFirst = first.items.size;
  if (!nFirst || !second.items.size) {
    return 0;
  }
  for (const int layer : search.layers) {
    if (layer < 0 || size_t(layer) >= frame.view().nLayers) {
      throw std::out_of_range{"GPU cell layers"};
    }
    frame.waitLayer(lease.stream, layer);
  }
  search.first = first.items.data();
  search.second = second.items.data();
  search.secondLookup = second.lookup.data();
  int* lookup = results.lookup.reserve(Lease::items(nFirst + 1));
  results.lookup.size = nFirst + 1;
  auto* counter = lease.take<int>(1);
  lease.zero(counter, sizeof(int));
  const int capacity = static_cast<int>(std::min<size_t>(results.capacity, std::numeric_limits<int>::max()));
  const Appender<Triplet> append{results.items.data(), capacity, counter};
  const auto* deviceView = frame.deviceView();
  lease.launch(nFirst, [=] __device__(size_t i) { forEachCell(*deviceView, search, static_cast<int>(i), append); });
  const int emitted = lease.get(counter);
  if (emitted > capacity) {
    return emitted;
  }
  auto* cells = results.items.data();
  lease.sortByKey(cells, emitted, [] __device__(const Triplet& cell) {
    return (uint64_t(uint32_t(cell.getFirstTrackletIndex())) << 32) | uint32_t(cell.getSecondTrackletIndex());
  }, 64);
  lease.buildLookup(cells, emitted, lookup, nFirst, [] __device__(const Triplet& cell) { return cell.getFirstTrackletIndex(); });
  lease.sync();
  results.items.size = emitted;
  return emitted;
}

int computeNeighbours(TimeFrameGPU& frame, int targetPath, std::span<const NeighbourSearch> sources)
{
  if (sources.empty()) {
    throw std::invalid_argument{"GPU neighbour search without sources"};
  }
  for (const auto& search : sources) {
    if (search.sourcePath < 0 || targetPath < 0 || size_t(std::max(search.sourcePath, targetPath)) >= frame.getNPaths() || search.sourcePath == targetPath) {
      throw std::out_of_range{"GPU neighbour paths"};
    }
  }
  Lease lease{frame};
  requireLoaded(frame);
  auto& target = frame.cells(targetPath);
  auto& results = frame.neighbours(targetPath);
  const size_t nTargets = target.items.size;
  auto* counter = lease.take<int>(1);
  lease.zero(counter, sizeof(int));
  const int capacity = static_cast<int>(std::min<size_t>(results.capacity, std::numeric_limits<int>::max()));
  const Appender<CellNeighbour> append{results.items.data(), capacity, counter};
  const auto* deviceView = frame.deviceView();
  Triplet* targetCells = target.items.data();
  for (auto search : sources) {
    const auto& cells = frame.cells(search.sourcePath);
    if (!cells.items.size || !nTargets || !target.lookup.size) {
      continue;
    }
    search.sources = cells.items.data();
    search.targets = targetCells;
    search.targetLookup = target.lookup.data();
    search.lookupSize = target.lookup.size;
    lease.launch(cells.items.size, [=] __device__(size_t i) {
      const int level = search.sources[i].getLevel() + 1;
      forEachNeighbour(*deviceView, search, static_cast<int>(i), [&](const CellNeighbour& neighbour) {
        append(neighbour);
        atomicMax(targetCells[neighbour.nextCell].getLevelPtr(), level);
      });
    });
  }
  const int emitted = lease.get(counter);
  if (emitted > capacity) {
    return emitted;
  }
  auto* neighbours = results.items.data();
  int* lookup = results.lookup.reserve(Lease::items(nTargets + 1));
  results.lookup.size = nTargets + 1;
  lease.sortByKey(neighbours, emitted, [] __device__(const CellNeighbour& neighbour) {
    return (uint64_t(uint32_t(neighbour.cellPath)) << 32) | uint32_t(neighbour.cell);
  }, 64);
  lease.sortByKey(neighbours, emitted, [] __device__(const CellNeighbour& neighbour) { return uint64_t(uint32_t(neighbour.nextCell)); }, 32);
  lease.buildLookup(neighbours, emitted, lookup, nTargets, [] __device__(const CellNeighbour& neighbour) { return neighbour.nextCell; });
  lease.sync();
  results.items.size = emitted;
  return emitted;
}

void bindRoadGraph(TimeFrameGPU& frame)
{
  std::lock_guard lock{frame.roadMutex};
  requireLoaded(frame);
  std::vector<RoadGraphPath> graph(frame.getNPaths());
  for (size_t path = 0; path < graph.size(); ++path) {
    const auto& cells = frame.cells(path);
    const auto& neighbours = frame.neighbours(path);
    const bool linked = neighbours.items.size;
    graph[path] = {cells.items.data(), cells.items.size, linked ? neighbours.lookup.data() : nullptr, linked ? neighbours.lookup.size : 0,
                   neighbours.items.data(), neighbours.items.size};
  }
  frame.upload(frame.roadGraph, std::span<const RoadGraphPath>{graph});
  frame.roadGraphGeneration = frame.getCellGeneration();
}

size_t extendRoads(TimeFrameGPU& frame, const RoadExtension& extension)
{
  Lease lease{frame};
  std::lock_guard lock{frame.roadMutex};
  requireLoaded(frame);
  if (frame.roadGraphGeneration != frame.getCellGeneration() || frame.roadGraph.size != frame.getNPaths()) {
    throw std::runtime_error{"GPU road: road graph missing or stale"};
  }
  const bool reuse = extension.reusePreviousSeeds;
  if (!reuse && (extension.startPath < 0 || size_t(extension.startPath) >= frame.getNPaths())) {
    throw std::out_of_range{"GPU road start path"};
  }
  const RoadResult* previous = reuse ? frame.roads[frame.lastRoads].data() : nullptr;
  const size_t n = reuse ? frame.roads[frame.lastRoads].size : frame.cells(extension.startPath).items.size;
  frame.lastRoads = 1 - frame.lastRoads;
  auto& output = frame.roads[frame.lastRoads];
  output.size = 0;
  if (!n) {
    return 0;
  }
  for (size_t layer = 0; layer < frame.view().nLayers; ++layer) {
    frame.waitLayer(lease.stream, static_cast<int>(layer));
  }
  const auto* view = frame.deviceView();
  const RoadStep step{frame.roadGraph.data(), frame.roadGraph.size, extension.nLayers, extension.currentLevel,
                      {lease.upload(extension.surfaces.data(), extension.surfaces.size()), static_cast<uint32_t>(extension.surfaces.size())}};
  const int startPath = extension.startPath;
  auto* ranges = lease.take<RoadRange>(n);
  auto* targetCounts = lease.take<size_t>(n + 1);
  auto* targetOffsets = lease.take<size_t>(n + 1);
  auto* jobCounts = lease.take<size_t>(n + 1);
  auto* jobOffsets = lease.take<size_t>(n + 1);
  auto* error = lease.take<unsigned>(1);
  lease.zero(error, sizeof(unsigned));
  lease.launch(n, [=] __device__(size_t i) {
    auto& range = ranges[i];
    unsigned bits;
    if (previous) {
      const auto& emission = previous[i].emission;
      range.path = emission.cellPathId;
      bits = roadRange(*view, step, range.path, emission.cellId, emission.seed.getLevel(), false, range.begin, range.end);
    } else {
      range.path = startPath;
      bits = roadRange(*view, step, range.path, static_cast<int>(i), step.paths[startPath].cells[i].getLevel(), true, range.begin, range.end);
    }
    if (bits) {
      atomicOr(error, bits);
    }
    targetCounts[i] = range.end - range.begin;
    jobCounts[i] = range.end != range.begin;
  });
  lease.scan(targetCounts, targetOffsets, n);
  lease.scan(jobCounts, jobOffsets, n);
  throwRoadError(lease.get(error));
  const size_t nTargets = lease.get(targetOffsets + n);
  const size_t nJobs = lease.get(jobOffsets + n);
  if (!nJobs) {
    return 0;
  }
  auto* jobs = lease.take<RoadDeviceJob>(nJobs);
  auto* initial = reuse ? nullptr : lease.take<SeedInput>(nJobs);
  auto* targets = lease.take<RoadTarget>(nTargets);
  lease.launch(n, [=] __device__(size_t i) {
    if (jobOffsets[i] == jobOffsets[i + 1]) {
      return;
    }
    const auto range = ranges[i];
    const size_t job = jobOffsets[i];
    const size_t first = targetOffsets[i];
    jobs[job] = {i, first, first + (range.end - range.begin), job, !previous};
    const auto& graph = step.paths[range.path];
    unsigned bits = previous ? 0 : makeSeedInput(*view, step.catalog, resolveCell(*view, graph.cells[i]), initial + job);
    for (int neighbour = range.begin; neighbour < range.end; ++neighbour) {
      bits |= makeRoadTarget(*view, step, graph.neighbours[neighbour], targets[first + neighbour - range.begin]);
    }
    if (bits) {
      atomicOr(error, bits);
    }
  });
  const size_t total = runRoads(lease, output, jobs, nJobs, initial, nullptr, previous, targets, extension.bz, extension.maxChi2);
  if (const unsigned bits = lease.get(error)) {
    output.size = 0;
    throwRoadError(bits);
  }
  return total;
}

void computeRoads(TimeFrameGPU& frame, const RoadInput& input, bounded_vector<RoadResult>& output)
{
  output.clear();
  const bool reuse = input.reusePreviousSeeds;
  if (reuse && frame.roads[frame.lastRoads].size != input.nSources) {
    throw std::invalid_argument{"GPU road predecessor count mismatch"};
  }
  // Construct before the lease: on exception the stream drains before these
  // pageable upload sources are destroyed.
  auto* resource = output.get_allocator().resource();
  bounded_vector<RoadDeviceJob> jobs{resource};
  bounded_vector<SeedInput> initial{resource};
  bounded_vector<TrackSeed> seeds{resource};
  for (const auto& job : input.jobs) {
    if (job.endTarget < job.firstTarget || job.endTarget > input.targets.size()) {
      throw std::out_of_range{"GPU road target range"};
    }
    if (job.firstTarget == job.endTarget) {
      continue;
    }
    if (reuse && job.initialize) {
      throw std::invalid_argument{"GPU road: previous seeds cannot be initialised"};
    }
    if (reuse && job.source >= input.nSources) {
      throw std::out_of_range{"GPU road predecessor index"};
    }
    jobs.push_back({job.source, job.firstTarget, job.endTarget, job.initialize ? initial.size() : seeds.size(), job.initialize});
    if (job.initialize) {
      initial.push_back(job.initial);
    } else if (!reuse) {
      seeds.push_back(job.seed);
    }
  }
  Lease lease{frame};
  std::lock_guard lock{frame.roadMutex};
  const RoadResult* previous = reuse ? frame.roads[frame.lastRoads].data() : nullptr;
  frame.lastRoads = 1 - frame.lastRoads;
  auto& roads = frame.roads[frame.lastRoads];
  roads.size = 0;
  if (jobs.empty()) {
    return;
  }
  runRoads(lease, roads, lease.upload(jobs.data(), jobs.size()), jobs.size(), lease.upload(initial.data(), initial.size()),
           lease.upload(seeds.data(), seeds.size()), previous, lease.upload(input.targets.data(), input.targets.size()), input.bz, input.maxChi2);
  frame.download(roads, output);
}

size_t selectRoadSeeds(TimeFrameGPU& frame, const RoadSeedSelector& selector)
{
  Lease lease{frame};
  std::lock_guard lock{frame.roadMutex};
  requireLoaded(frame);
  const auto& roads = frame.roads[frame.lastRoads];
  const size_t n = roads.size;
  auto& seeds = frame.trackSeeds;
  const size_t first = seeds.size;
  if (!n) {
    return first;
  }
  const auto* results = roads.data();
  auto* flags = lease.take<size_t>(n + 1);
  auto* offsets = lease.take<size_t>(n + 1);
  lease.launch(n, [=] __device__(size_t i) { flags[i] = selector(results[i].emission.seed); });
  lease.scan(flags, offsets, n);
  const size_t selected = lease.get(offsets + n);
  if (!selected) {
    return first;
  }
  if ((first + selected) * sizeof(TrackSeed) > seeds.buffer.usable()) { // grow keeping the seeds
    auto grown = seeds.buffer.sibling();
    grown.reserve(std::max((first + selected) * sizeof(TrackSeed), 2 * seeds.buffer.usable()));
    if (first) {
      check(cudaMemcpyAsync(grown.data, seeds.data(), first * sizeof(TrackSeed), cudaMemcpyDeviceToDevice, lease.stream));
      lease.sync();
    }
    seeds.buffer = std::move(grown);
  }
  auto* output = seeds.data() + first;
  lease.launch(n, [=] __device__(size_t i) {
    if (offsets[i + 1] != offsets[i]) {
      output[offsets[i]] = results[i].emission.seed;
    }
  });
  lease.sync();
  seeds.size = first + selected;
  return seeds.size;
}

size_t refitTrackSeeds(TimeFrameGPU& frame, std::span<const SurfaceDescriptor> surfaces, const RefitParameters& parameters,
                       std::span<const float> minPt)
{
  if (surfaces.size() > MaxLayoutSurfaces || minPt.size() > MaxLayoutSurfaces) {
    throw std::out_of_range{"GPU refit surface count"};
  }
  Lease lease{frame};
  std::lock_guard lock{frame.roadMutex};
  requireLoaded(frame);
  const size_t n = frame.trackSeeds.size;
  frame.tracks.size = 0;
  if (!n) {
    return 0;
  }
  const SurfaceCatalogView catalog{lease.upload(surfaces.data(), surfaces.size()), static_cast<uint32_t>(surfaces.size())};
  const auto* deviceMinPt = lease.upload(minPt.data(), minPt.size());
  const size_t nMinPt = minPt.size();
  const auto* view = frame.deviceView();
  const auto* seeds = frame.trackSeeds.data();
  auto* results = frame.refitResults.reserve(n);
  frame.refitResults.size = n;
  lease.launch(n, [=] __device__(size_t i) {
    RefitJob job;
    prepareRefitJob(*view, seeds[i], deviceMinPt, nMinPt, job);
    results[i] = evaluateRefit(job, catalog, parameters);
  }, 64);
  auto* flags = lease.take<size_t>(n + 1);
  auto* offsets = lease.take<size_t>(n + 1);
  lease.launch(n, [=] __device__(size_t i) { flags[i] = results[i].accepted; });
  lease.scan(flags, offsets, n);
  const size_t accepted = lease.get(offsets + n);
  if (!accepted) {
    return 0;
  }
  auto* tracks = frame.tracks.reserve(accepted);
  lease.launch(n, [=] __device__(size_t i) {
    if (offsets[i + 1] != offsets[i]) {
      TrackingCandidate track{seeds[i]};
      track.track.innerState = results[i].inner;
      track.track.outerState = results[i].outer;
      track.track.chi2 = results[i].chi2;
      tracks[offsets[i]] = track;
    }
  });
  // More clusters first, then lower chi2 (non-negative floats order by their bits).
  lease.sortByKey(tracks, accepted, [] __device__(const TrackingCandidate& track) {
    const float chi2 = track.track.chi2 > 0.f ? track.track.chi2 : 0.f;
    uint32_t bits;
    memcpy(&bits, &chi2, sizeof(bits));
    return (static_cast<uint64_t>(MaxLayoutSurfaces - track.seed.getActiveLayerCount()) << 32) | bits;
  }, 64);
  lease.sync();
  frame.tracks.size = accepted;
  return accepted;
}

void refitTracks(TimeFrameGPU& frame, const RefitInput& input, bounded_vector<RefitResult>& output)
{
  output.clear();
  if (input.surfaces.size() > MaxLayoutSurfaces) {
    throw std::out_of_range{"GPU refit surface count"};
  }
  for (const auto& job : input.jobs) {
    if (job.nSlots > MaxLayoutSurfaces || job.nPoints > MaxLayoutSurfaces) {
      throw std::out_of_range{"GPU refit measurement count"};
    }
  }
  if (input.jobs.empty()) {
    return;
  }
  output.reserve(input.jobs.size());
  Lease lease{frame};
  const auto* jobs = lease.upload(input.jobs.data(), input.jobs.size());
  const SurfaceCatalogView catalog{lease.upload(input.surfaces.data(), input.surfaces.size()), static_cast<uint32_t>(input.surfaces.size())};
  const auto parameters = input.parameters;
  auto* results = frame.refitResults.reserve(input.jobs.size());
  frame.refitResults.size = input.jobs.size();
  lease.launch(input.jobs.size(), [=] __device__(size_t i) { results[i] = evaluateRefit(jobs[i], catalog, parameters); }, 64);
  lease.sync();
  frame.download(frame.refitResults, output);
}

} // namespace ITSMFT_TRACKING_GPU_BACKEND
} // namespace o2::itsmft::tracking::gpu

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// SPDX-License-Identifier: GPL-3.0-only
#ifndef GTRK_TOF_STUDY_PREDICTIONCACHE_H_
#define GTRK_TOF_STUDY_PREDICTIONCACHE_H_

// Local, same-build prediction shards. Propagation is independent for each
// track/cluster pair; arbitration is performed once after ALL shards are read.
#include "GenericTracking/TOF/Matcher.h"
#include "GenericTracking/Tracking/Propagator.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace o2::gtrk::tof::study
{
struct Fingerprint {
  std::uint64_t value = 14695981039346656037ULL;
  template <typename T>
  void add(const T& object)
  {
    static_assert(std::is_trivially_copyable_v<T>);
    const auto* bytes = reinterpret_cast<const unsigned char*>(&object);
    for (std::size_t i = 0; i < sizeof(T); ++i) {
      value = (value ^ bytes[i]) * 1099511628211ULL;
    }
  }
  void addSurface(const TrackingSurface& s)
  {
    add(s.origin);
    add(s.axis0);
    add(s.axis1);
    add(s.normal);
    add(s.bounds.halfX);
    add(s.bounds.halfY);
    add(s.id);
    add(s.layer);
    add(s.radius);
    add(s.surfacePhi);
    add(s.depth);
  }
};
inline std::uint64_t fingerprint(std::span<const Track> tracks, std::span<const Hit> hits, Propagator prop, const Parameters& p, double margin)
{
  Fingerprint f;
  f.add(prop.bzKiloGauss);
  f.add(prop.massGeV);
  f.add(p.maxSpatialChi2);
  f.add(p.maxTimeChi2);
  f.add(p.timeMarginNS);
  f.add(p.maxStepCM);
  f.add(p.maxPathCM);
  f.add(margin);
  for (const auto& t : tracks) {
    f.add(t.valid);
    if (!t.valid) {
      continue;
    }
    f.add(t.state);
    f.add(t.outwardState);
    f.addSurface(t.surface);
    f.addSurface(t.outwardSurface);
    f.add(t.boxStart.bc);
    f.add(t.boxStart.fineNS);
    f.add(t.boxWidthNS);
  }
  for (const auto& h : hits) {
    f.addSurface(h.surface);
    f.add(h.time.bc);
    f.add(h.time.fineNS);
    f.add(h.loc0);
    f.add(h.loc1);
    f.add(h.varLoc0);
    f.add(h.varLoc1);
    f.add(h.sigmaTimeNS);
    f.add(h.index);
    f.add(h.layer);
  }
  return f.value;
}

class TimeIndex
{
 public:
  explicit TimeIndex(std::span<const Hit> hits) : mHits(hits), mRows(hits.size())
  {
    std::iota(mRows.begin(), mRows.end(), 0);
    std::sort(mRows.begin(), mRows.end(), [&](auto a, auto b) {
      const auto dt = differenceNS(hits[a].time, hits[b].time);
      return dt == 0. ? a < b : dt < 0.;
    });
  }
  std::vector<std::size_t> rows(const Track& track, double margin) const
  {
    if (!track.valid || !(track.boxWidthNS > 0.)) {
      return {};
    }
    const auto start = std::lower_bound(mRows.begin(), mRows.end(), -margin, [&](auto row, double limit) { return differenceNS(mHits[row].time, track.boxStart) < limit; });
    const auto end = std::upper_bound(start, mRows.end(), track.boxWidthNS + margin, [&](double limit, auto row) { return limit < differenceNS(mHits[row].time, track.boxStart); });
    std::vector<std::size_t> result(start, end);
    std::sort(result.begin(), result.end());
    return result;
  }

 private:
  std::span<const Hit> mHits;
  std::vector<std::size_t> mRows;
};

struct CachedPrediction {
  std::uint64_t track = 0, hit = 0;
  Prediction prediction;
};
struct CacheHeader {
  std::uint64_t magic = 0x4754524b544f4631ULL;
  std::uint64_t recordSize = sizeof(CachedPrediction);
  std::uint64_t fingerprint = 0, tracks = 0, hits = 0, shard = 0, shards = 0, records = 0, exactCalls = 0;
};
static_assert(std::is_trivially_copyable_v<CachedPrediction> && std::is_trivially_copyable_v<CacheHeader>);

template <typename T>
void write(std::ostream& out, const T& object)
{
  out.write(reinterpret_cast<const char*>(&object), sizeof(T));
  if (!out) {
    throw std::runtime_error("cannot write prediction cache");
  }
}
template <typename T>
void read(std::istream& in, T& object)
{
  if (!in.read(reinterpret_cast<char*>(&object), sizeof(T))) {
    throw std::runtime_error("truncated prediction cache");
  }
}

class PredictionCache
{
 public:
  PredictionCache(std::size_t tracks, std::size_t hits, std::uint64_t fingerprint, const std::vector<std::string>& paths) : mValues(tracks), mRows(tracks)
  {
    std::vector<bool> seen;
    for (const auto& path : paths) {
      std::ifstream file(path, std::ios::binary);
      CacheHeader header;
      read(file, header);
      if (header.magic != CacheHeader{}.magic || header.recordSize != sizeof(CachedPrediction) || header.fingerprint != fingerprint || header.tracks != tracks || header.hits != hits || header.shards != paths.size() || header.shard >= paths.size()) {
        throw std::runtime_error("prediction shard does not describe these inputs: " + path);
      }
      if (seen.empty()) {
        seen.resize(paths.size());
      }
      if (seen[header.shard]) {
        throw std::runtime_error("duplicate prediction shard");
      }
      seen[header.shard] = true;
      exactCalls += header.exactCalls;
      for (std::uint64_t i = 0; i < header.records; ++i) {
        CachedPrediction record;
        read(file, record);
        if (record.track >= tracks || record.hit >= hits || record.track % header.shards != header.shard) {
          throw std::runtime_error("prediction cache row outside its shard");
        }
        mValues[record.track].push_back(record);
      }
      if (file.peek() != std::char_traits<char>::eof()) {
        throw std::runtime_error("trailing prediction cache data");
      }
    }
    for (std::size_t t = 0; t < tracks; ++t) {
      auto& values = mValues[t];
      std::sort(values.begin(), values.end(), [](const auto& a, const auto& b) { return a.hit < b.hit; });
      for (const auto& value : values) {
        if (!mRows[t].empty() && mRows[t].back() == value.hit) {
          throw std::runtime_error("duplicate prediction pair");
        }
        mRows[t].push_back(value.hit);
      }
    }
  }
  std::span<const std::size_t> rows(std::size_t track) const { return mRows[track]; }
  Prediction prediction(std::size_t track, std::size_t hit) const
  {
    const auto& values = mValues[track];
    const auto at = std::lower_bound(values.begin(), values.end(), hit, [](const auto& a, auto b) { return a.hit < b; });
    return at != values.end() && at->hit == hit ? at->prediction : Prediction{};
  }
  std::uint64_t exactCalls = 0;

 private:
  std::vector<std::vector<CachedPrediction>> mValues;
  std::vector<std::vector<std::size_t>> mRows;
};
} // namespace o2::gtrk::tof::study
#endif

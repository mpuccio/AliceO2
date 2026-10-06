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

// Compare complete serialized branch payloads, including floating-point bits.
// ROOT file headers (timestamps, compression, etc.) are deliberately excluded.
#include <cstring>
#include <stdexcept>
#include <TBufferFile.h>
#include <TClass.h>

template <typename T>
bool sameTrackingBranch(const T& a, const T& b)
{
  TBufferFile left{TBuffer::kWrite}, right{TBuffer::kWrite};
  auto* type = TClass::GetClass(typeid(T));
  if (!type) {
    throw std::runtime_error{"Missing branch dictionary"};
  }
  left.WriteObjectAny(&a, type);
  right.WriteObjectAny(&b, type);
  return left.Length() == right.Length() && !std::memcmp(left.Buffer(), right.Buffer(), left.Length());
}

#include <algorithm>
#include <cmath>
#include <map>
#include <iostream>
#include <vector>

// Validation only: abs + relative tolerances in the stored parameter units.
// Reject NaNs/infinities; do not let them satisfy the tolerance accidentally.
inline bool closeTrackingFloat(float a, float b, float absolute, float relative)
{
  return std::isfinite(a) && std::isfinite(b) &&
         std::abs(a - b) <= absolute + relative * std::max(std::abs(a), std::abs(b));
}

template <typename State>
bool equivalentTrackingState(const State& a, State& normalized)
{
  if (!closeTrackingFloat(a.getX(), normalized.getX(), 1.e-5f, 1.e-4f) ||
      !closeTrackingFloat(a.getAlpha(), normalized.getAlpha(), 1.e-6f, 1.e-4f)) {
    return false;
  }
  for (int i = 0; i < 5; ++i) {
    if (!closeTrackingFloat(a.getParam(i), normalized.getParam(i), 1.e-5f, 2.e-3f)) {
      return false;
    }
  }
  for (int i = 0; i < 15; ++i) {
    if (!closeTrackingFloat(a.getCov()[i], normalized.getCov()[i], 1.e-8f, 2.e-2f)) {
      return false;
    }
  }
  normalized.setX(a.getX());
  normalized.setAlpha(a.getAlpha());
  for (int i = 0; i < 5; ++i) {
    normalized.setParam(a.getParam(i), i);
  }
  normalized.setCov(a.getCov());
  return true;
}

template <typename Track>
bool equivalentTrackingTrack(const Track& a, const Track& b)
{
  if (!closeTrackingFloat(a.getChi2(), b.getChi2(), 1.e-2f, 5.e-3f)) {
    return false;
  }
  Track normalized = b;
  if (!equivalentTrackingState(a, normalized) || !equivalentTrackingState(a.getParamOut(), normalized.getParamOut())) {
    return false;
  }
  normalized.setChi2(a.getChi2());
  normalized.setFirstClusterEntry(a.getFirstClusterEntry());
  // All remaining metadata, including timestamps, patterns, flags, PID and
  // charge, must match exactly. Normalize only already-validated floats.
  return sameTrackingBranch(a, normalized);
}

template <typename Tracks, typename Indices, typename ROFs, typename Labels>
bool equivalentTrackingEvent(const Tracks& cpu, const Indices& cpuIndices, const ROFs& cpuROFs, const Labels& cpuLabels,
                             const Tracks& gpu, const Indices& gpuIndices, const ROFs& gpuROFs, const Labels& gpuLabels)
{
  if (cpu.size() != gpu.size() || cpuIndices.size() != gpuIndices.size() ||
      cpuLabels.size() != cpu.size() || gpuLabels.size() != gpu.size() || !sameTrackingBranch(cpuROFs, gpuROFs)) {
    throw std::runtime_error{"Track/cluster/label counts or ROF metadata differ"};
  }
  const auto key = [](const auto& track, const auto& indices) {
    const int first = track.getFirstClusterEntry();
    const int count = track.getNumberOfClusters();
    if (first < 0 || count <= 0 || size_t(first) + count > indices.size()) {
      throw std::runtime_error{"Invalid track cluster range"};
    }
    return std::vector<int>(indices.begin() + first, indices.begin() + first + count);
  };
  size_t processed = 0;
  size_t mismatches = 0, missingIdentities = 0;
  size_t clusters = 0;
  for (const auto& rof : cpuROFs) {
    const size_t start = rof.getFirstEntry(), end = start + rof.getNEntries();
    if (start != processed || end > cpu.size()) {
      throw std::runtime_error{"Invalid ROF track range"};
    }
    std::map<std::vector<int>, std::vector<size_t>> byClusters;
    for (size_t i = start; i < end; ++i) {
      byClusters[key(gpu[i], gpuIndices)].push_back(i);
    }
    for (size_t i = start; i < end; ++i) {
      const auto clusterKey = key(cpu[i], cpuIndices);
      auto& candidates = byClusters[clusterKey];
      const auto match = std::find_if(candidates.begin(), candidates.end(), [&](size_t j) {
        return equivalentTrackingTrack(cpu[i], gpu[j]) && sameTrackingBranch(cpuLabels[i], gpuLabels[j]);
      });
      clusters += clusterKey.size();
      ++processed;
      if (match == candidates.end()) {
        ++mismatches;
        missingIdentities += candidates.empty();
        if (!candidates.empty() && mismatches <= 3) {
          const auto& a = cpu[i];
          const auto& b = gpu[candidates.front()];
          std::cerr << "Track " << i << " chi2 " << a.getChi2() << "/" << b.getChi2() << "\n";
          const auto diagnose = [](const auto& x, const auto& y) {
            for (int k = 0; k < 5; ++k) {
              std::cerr << " p" << k << " " << x.getParam(k) << "/" << y.getParam(k);
            }
            for (int k = 0; k < 15; ++k) {
              if (!closeTrackingFloat(x.getCov()[k], y.getCov()[k], 1.e-8f, 2.e-2f)) {
                std::cerr << " c" << k << " " << x.getCov()[k] << "/" << y.getCov()[k];
              }
            }
            std::cerr << "\n";
          };
          diagnose(a, b);
          diagnose(a.getParamOut(), b.getParamOut());
        }
        continue;
      }
      candidates.erase(match);
    }
  }
  if (processed != cpu.size() || clusters != cpuIndices.size()) {
    throw std::runtime_error{"Unaccounted tracks or cluster references"};
  }
  if (mismatches) {
    throw std::runtime_error{std::to_string(missingIdentities) + " unmatched cluster identities; " +
                             std::to_string(mismatches - missingIdentities) + " metadata, label or numerical mismatches"};
  }
  return true;
}

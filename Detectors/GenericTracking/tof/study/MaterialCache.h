// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// SPDX-License-Identifier: GPL-3.0-only
#ifndef GTRK_TOF_STUDY_MATERIALCACHE_H_
#define GTRK_TOF_STUDY_MATERIALCACHE_H_

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <utility>

namespace o2::gtrk::tof::study
{
/// Reuse bit-identical material segments along the common bounded-step prefixes
/// of different target-plane transports. No coordinate quantization or material
/// approximation. Bound the memory per track; excess keys use the provider.
template <typename Material>
class MaterialCache
{
  using Key = std::array<std::uint32_t, 6>;
  struct Hash {
    std::size_t operator()(const Key& key) const
    {
      std::uint64_t value = 14695981039346656037ULL;
      for (const auto word : key) {
        value = (value ^ word) * 1099511628211ULL;
      }
      return value;
    }
  };
  using Budget = decltype(std::declval<const Material&>().getMatBudget(0.F, 0.F, 0.F, 0.F, 0.F, 0.F));

 public:
  explicit MaterialCache(Material material = {}, std::size_t limit = 65536) : mMaterial(std::move(material)), mLimit(limit) {}
  Budget getMatBudget(float x0, float y0, float z0, float x1, float y1, float z1) const
  {
    ++queries;
    const Key key{std::bit_cast<std::uint32_t>(x0), std::bit_cast<std::uint32_t>(y0), std::bit_cast<std::uint32_t>(z0), std::bit_cast<std::uint32_t>(x1), std::bit_cast<std::uint32_t>(y1), std::bit_cast<std::uint32_t>(z1)};
    if (const auto at = mCache.find(key); at != mCache.end()) {
      ++hits;
      return at->second;
    }
    const auto result = mMaterial.getMatBudget(x0, y0, z0, x1, y1, z1);
    if (mCache.size() < mLimit) {
      mCache.emplace(key, result);
    }
    return result;
  }
  mutable std::uint64_t queries = 0, hits = 0;

 private:
  Material mMaterial;
  std::size_t mLimit;
  mutable std::unordered_map<Key, Budget, Hash> mCache;
};
} // namespace o2::gtrk::tof::study
#endif

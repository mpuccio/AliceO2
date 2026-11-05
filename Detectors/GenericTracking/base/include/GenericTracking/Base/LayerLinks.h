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

/// \file LayerLinks.h
/// \brief What follows at compile time from a detector's link table.
///
/// `Det::layerLinks` says which layer a track reaches next. The grid and the
/// CA read it; neither encodes
/// adjacency again. Everything here is a constexpr function of that table:
///
///  - the links leaving and entering each layer, as index lists;
///  - a topological order of the layers -- a table with a cycle does not
///    compile;
///  - each layer's depth, and the longest chain, which bounds the level;
///  - the first link joining two layers (`linkOf`).
///
/// The index lists do not depend on the order the table is written in. The
/// tool emits it sorted by `(from, to)` only so that it reads well.

#ifndef GTRK_BASE_LAYERLINKS_H_
#define GTRK_BASE_LAYERLINKS_H_

#include "GenericTracking/Base/LayerSpec.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <limits>

namespace o2::gtrk
{

template <DetectorLayout Det>
inline constexpr std::size_t nLinks = std::size(Det::layerLinks);

/// Link indices grouped by one of their ends: the links of layer `l` are
/// `links[offsets[l]]` up to `links[offsets[l + 1]]`, in table order.
template <DetectorLayout Det>
struct LinksByLayer {
  std::array<std::size_t, static_cast<std::size_t>(Detector<Det>::nLayers) + 1> offsets{};
  std::array<std::size_t, nLinks<Det>> links{};

  [[nodiscard]] constexpr std::size_t begin(LayerID layer) const { return offsets[layer]; }
  [[nodiscard]] constexpr std::size_t end(LayerID layer) const { return offsets[layer + 1U]; }
  [[nodiscard]] constexpr std::size_t count(LayerID layer) const { return end(layer) - begin(layer); }
};

namespace detail
{

/// A counting sort of the link indices by their source or target layer.
template <DetectorLayout Det, bool ByTarget>
[[nodiscard]] constexpr LinksByLayer<Det> groupLinks()
{
  const auto key = [](const LayerLink& link) { return static_cast<std::size_t>(ByTarget ? link.to : link.from); };
  LinksByLayer<Det> result{};
  for (const auto& link : Det::layerLinks) {
    ++result.offsets[key(link) + 1U];
  }
  for (std::size_t layer = 0; layer + 1U < result.offsets.size(); ++layer) {
    result.offsets[layer + 1U] += result.offsets[layer];
  }
  auto next = result.offsets;
  for (std::size_t i = 0; i < nLinks<Det>; ++i) {
    result.links[next[key(Det::layerLinks[i])]++] = i;
  }
  return result;
}

} // namespace detail

/// The links leaving each layer.
template <DetectorLayout Det>
inline constexpr LinksByLayer<Det> outLinks = detail::groupLinks<Det, false>();

/// The links entering each layer.
template <DetectorLayout Det>
inline constexpr LinksByLayer<Det> inLinks = detail::groupLinks<Det, true>();

namespace detail
{

template <DetectorLayout Det>
struct TopologicalSort {
  std::array<LayerID, static_cast<std::size_t>(Detector<Det>::nLayers)> order{};
  std::size_t sorted = 0;
};

/// Kahn's algorithm, with `order` doubling as the queue. Layers with nothing
/// entering them start it in index order, so a chain of barrel layers comes
/// out as 0, 1, 2, ... A cycle leaves its layers unsorted.
template <DetectorLayout Det>
[[nodiscard]] constexpr TopologicalSort<Det> topologicalSort()
{
  constexpr auto kLayers = static_cast<std::size_t>(Detector<Det>::nLayers);
  TopologicalSort<Det> result{};
  std::array<std::size_t, kLayers> pending{};
  for (const auto& link : Det::layerLinks) {
    ++pending[link.to];
  }
  for (std::size_t layer = 0; layer < kLayers; ++layer) {
    if (pending[layer] == 0) {
      result.order[result.sorted++] = static_cast<LayerID>(layer);
    }
  }
  for (std::size_t head = 0; head < result.sorted; ++head) {
    const auto layer = result.order[head];
    for (auto i = outLinks<Det>.begin(layer); i < outLinks<Det>.end(layer); ++i) {
      const auto to = Det::layerLinks[outLinks<Det>.links[i]].to;
      if (--pending[to] == 0) {
        result.order[result.sorted++] = to;
      }
    }
  }
  return result;
}

template <DetectorLayout Det>
consteval std::array<LayerID, static_cast<std::size_t>(Detector<Det>::nLayers)> checkedLayerOrder()
{
  constexpr auto sort = topologicalSort<Det>();
  static_assert(sort.sorted == static_cast<std::size_t>(Detector<Det>::nLayers), "Det::layerLinks contains a cycle");
  return sort.order;
}

} // namespace detail

/// The layers in an order where every link points forward, so that one pass
/// over them sees every link's source layer before its target.
template <DetectorLayout Det>
inline constexpr auto layerOrder = detail::checkedLayerOrder<Det>();

/// Each layer's depth in the link DAG: 0 for a layer no link enters, else one
/// more than the deepest source of a link into it.
template <DetectorLayout Det>
inline constexpr auto layerDepth = [] {
  std::array<std::size_t, static_cast<std::size_t>(Detector<Det>::nLayers)> depth{};
  for (const auto layer : layerOrder<Det>) {
    for (auto i = outLinks<Det>.begin(layer); i < outLinks<Det>.end(layer); ++i) {
      const auto to = Det::layerLinks[outLinks<Det>.links[i]].to;
      depth[to] = std::max(depth[to], depth[layer] + 1U);
    }
  }
  return depth;
}();

/// The layers by depth (`layerDepth`), a topological order too, and where each
/// depth's run begins in it: `order[groups[d]]` up to `order[groups[d + 1]]`
/// are the layers of depth `d`, `count` depths in all. The layers of one depth
/// have no link between them and every link into them from a smaller depth,
/// so an ordered pass over the layers may take a depth's layers together.
template <DetectorLayout Det>
struct LayerDepthGroups {
  static constexpr std::size_t kLayers = static_cast<std::size_t>(Detector<Det>::nLayers);
  std::array<LayerID, kLayers> order{};
  std::array<std::size_t, kLayers + 1U> groups{};
  std::size_t count = 0U;
};

template <DetectorLayout Det>
inline constexpr LayerDepthGroups<Det> layerDepthGroups = [] {
  LayerDepthGroups<Det> result{};
  std::size_t next = 0U;
  for (std::size_t depth = 0U; next < LayerDepthGroups<Det>::kLayers; ++depth) {
    result.groups[depth] = next;
    for (std::size_t layer = 0U; layer < LayerDepthGroups<Det>::kLayers; ++layer) {
      if (layerDepth<Det>[layer] == depth) {
        result.order[next++] = static_cast<LayerID>(layer);
      }
    }
    result.count = depth + 1U;
  }
  result.groups[result.count] = LayerDepthGroups<Det>::kLayers;
  return result;
}();

/// No link of the layout.
inline constexpr std::size_t kNoLink = std::numeric_limits<std::size_t>::max();

/// Per (from, to) layer pair, the first link of Det's table joining them;
/// kNoLink for none.
template <DetectorLayout Det>
inline constexpr auto kLinkIndex = [] {
  constexpr auto kLayers = static_cast<std::size_t>(Detector<Det>::nLayers);
  std::array<std::size_t, kLayers * kLayers> result{};
  result.fill(kNoLink);
  for (std::size_t link = 0; link < nLinks<Det>; ++link) {
    auto& at = result[(static_cast<std::size_t>(Det::layerLinks[link].from) * kLayers) + Det::layerLinks[link].to];
    if (at == kNoLink) {
      at = link;
    }
  }
  return result;
}();

/// The first link of Det's table from layer \a from to layer \a to; kNoLink
/// for none or a layer outside the layout.
template <DetectorLayout Det>
[[nodiscard]] constexpr std::size_t linkOf(LayerID from, LayerID to)
{
  constexpr auto kLayers = static_cast<std::size_t>(Detector<Det>::nLayers);
  const auto fromIndex = static_cast<std::size_t>(from);
  const auto toIndex = static_cast<std::size_t>(to);
  return fromIndex < kLayers && toIndex < kLayers ? kLinkIndex<Det>[(fromIndex * kLayers) + toIndex] : kNoLink;
}

/// The most layers any path through the links visits, the deepest layer's
/// depth + 1, which bounds a cell's level at `maxChainLayers - 2`.
template <DetectorLayout Det>
inline constexpr std::size_t maxChainLayers = layerDepthGroups<Det>.count;

} // namespace o2::gtrk

#endif // GTRK_BASE_LAYERLINKS_H_

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

/// \file GridSizing.h
/// \brief The grid's cells sized by the queries the configured iterations make.
///
/// A layer's cells follow from
/// the windows of the queries that search it, and only those. Two kinds of
/// query search a layer: the tracklet stage's, along a link the iteration
/// takes and with that link's pair window at the iteration's minimum pT,
/// into the link's target layer; and in an iteration that follows its tracks,
/// the follow's, along every link of the layout, outward into the target
/// layer and inward into the source layer, sized by the link's pair window.
/// A layer's windows are the widest and the narrowest of them per axis;
/// a layer no query enters is not searched and
/// has one cell per axis.

#ifndef GTRK_TRACKING_GRIDSIZING_H_
#define GTRK_TRACKING_GRIDSIZING_H_

#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LinkHoles.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"
#include "GenericTracking/Tracking/SearchWindow.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace o2::gtrk
{

/// Per layer, the widest and the narrowest window of the queries \a iterations
/// make into it.
///
/// The window is the material-uniform pair window (`layerPairWindow` without a
/// scattering table): the producer, which bins the rows, has no material map,
/// and only the cells' sizes follow from it, never which rows a query
/// accepts. Barrel <-> disc links do not size a layer: their window is per
/// seed, from a few percent of the target to
/// all of it, so no cell width follows from it.
template <DetectorLayout Det>
[[nodiscard]] std::array<GridLayerWindow, static_cast<std::size_t>(Detector<Det>::nLayers)> gridLayerWindows(std::span<const IterationSettings> iterations,
                                                                                                             const VertexConstraint& beamSpot)
{
  std::array<GridLayerWindow, static_cast<std::size_t>(Detector<Det>::nLayers)> result{};
  const auto widen = [&](LayerID layer, const LayerPairWindow& window) {
    auto& entry = result[layer];
    const auto phi = static_cast<double>(window.phiHalfWidth);
    const auto second = static_cast<double>(window.secondHalfWidth);
    entry.narrowPhiHalfWidth = entry.queried ? std::min(entry.narrowPhiHalfWidth, phi) : phi;
    entry.narrowSecondHalfWidth = entry.queried ? std::min(entry.narrowSecondHalfWidth, second) : second;
    entry.phiHalfWidth = std::max(entry.phiHalfWidth, phi);
    entry.secondHalfWidth = std::max(entry.secondHalfWidth, second);
    entry.queried = true;
  };
  constexpr auto kLayers = static_cast<std::size_t>(Detector<Det>::nLayers);
  // A link's window once per distinct window setting: the iterations with
  // equal settings share it, as do an iteration's tracklet and follow sizing.
  std::vector<std::pair<SearchWindowParams, std::array<std::optional<LayerPairWindow>, nLinks<Det>>>> derived;
  const auto windowOf = [&](const SearchWindowParams& params, std::size_t link) {
    auto entry = std::ranges::find_if(derived, [&](const auto& candidate) { return candidate.first == params; });
    if (entry == derived.end()) {
      derived.emplace_back(params, std::array<std::optional<LayerPairWindow>, nLinks<Det>>{});
      entry = std::prev(derived.end());
    }
    auto& window = entry->second[link];
    if (!window.has_value()) {
      const auto& pair = Det::layerLinks[link];
      window = layerPairWindow<Det>(pair.from, pair.to, params, beamSpot);
    }
    return *window;
  };
  for (const auto& settings : iterations) {
    // A displaced iteration's windows hold what its pair test admits, up to
    // half a turn and a target's whole length: sized by them the cells would be the prompt iterations' loss. Its
    // queries visit more cells instead; the cells never decide a pair.
    if (settings.windows.impactMax > 0.) {
      continue;
    }
    const auto taken = takenLinks<Det>(settings);
    const auto sizing = [&](std::size_t link) {
      const auto& pair = Det::layerLinks[link];
      return taken[link] && Detector<Det>::isBarrelLayer(pair.from) == Detector<Det>::isBarrelLayer(pair.to);
    };
    // Without skips, a link whose pairs can cross another layer is queried
    // only from the sources whose pairs can avoid it, at that layer's ends
    // (LinkHoles.h): it sizes a layer only where no link without holes
    // queries it into the target.
    std::array<bool, nLinks<Det>> edgeOnly{};
    std::array<bool, kLayers> holeFreeInto{};
    for (std::size_t link = 0; link < nLinks<Det>; ++link) {
      const auto& pair = Det::layerLinks[link];
      edgeOnly[link] = !settings.skipLinks && linkHoleLayers<Det>(pair.from, pair.to) != 0U;
      if (sizing(link) && !edgeOnly[link]) {
        holeFreeInto[pair.to] = true;
      }
    }
    for (std::size_t link = 0; link < nLinks<Det>; ++link) {
      if (!sizing(link)) {
        continue;
      }
      const auto& pair = Det::layerLinks[link];
      const auto window = windowOf(settings.windows, link);
      if (!(edgeOnly[link] && holeFreeInto[pair.to])) {
        widen(pair.to, window);
      }
    }
    // The follow searches both layers of every link of the layout, from
    // either end; its window is its prediction's, the pair window its size.
    if (settings.follow || settings.followEnd) {
      for (std::size_t link = 0; link < nLinks<Det>; ++link) {
        const auto& pair = Det::layerLinks[link];
        if (Detector<Det>::isBarrelLayer(pair.from) != Detector<Det>::isBarrelLayer(pair.to)) {
          continue;
        }
        const auto window = windowOf(settings.windows, link);
        widen(pair.from, window);
        widen(pair.to, window);
      }
    }
  }
  return result;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GRIDSIZING_H_

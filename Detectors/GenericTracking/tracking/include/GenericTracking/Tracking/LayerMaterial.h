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

/// \file LayerMaterial.h
/// \brief Each layer's average material budget, queried from the LUT once.
///
/// Host-side setup, once per run: the averages the pair cuts and search
/// windows are set up from; the fit and the follow query the LUT itself
/// (TrackFitter.h, Follow.h). A layer's budget is the
/// average `x/X0` of segments crossing it at normal incidence, from half-way
/// to the layer inside it to half-way to the layer outside, so that all
/// material between layers is attributed to the nearer one -- except a
/// passive scatterer's band, which is its own (`passiveMaterial`).

#ifndef GTRK_TRACKING_LAYERMATERIAL_H_
#define GTRK_TRACKING_LAYERMATERIAL_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/SearchWindow.h"

#include "DetectorsBase/MatLayerCylSet.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <utility>
#include <vector>

namespace o2::gtrk
{

template <DetectorLayout Det>
using LayerBudget = std::array<double, static_cast<std::size_t>(Detector<Det>::nLayers)>;

/// The LUT's `x/X0` across a band [\a low, \a high], averaged over azimuth
/// and over [\a axisMin, \a axisMax] of the other coordinate: radial segments
/// at each z for a cylinder; for a disc at \a position, segments from \a low
/// to \a high in z along the ray from the nominal vertex through each radius
/// at \a position, scaled to normal incidence by the cosine of their angle to
/// z. A segment parallel to z would do, but the LUT is a stack of cylindrical
/// layers traced by radius, and a segment that never changes radius reads as
/// empty.
///
/// On a cylinder the band is cut at the LUT's own layer boundaries, and each
/// piece is read once per cell of its LUT layer -- at each azimuth cell's
/// centre, and along z at each z cell's, weighted by how much of the cell the
/// range covers -- so that this is the LUT's own average. A fixed grid is not:
/// thin radial walls, as between the ALICE 3 VD's petals, are each averaged
/// into one cell, and a grid with the walls' symmetry hits every such cell or
/// none; a grid finer than a layer's cells weighs some of them twice.
///
/// A disc's average is the LUT's own as well, in closed form. A ray through
/// radius rho at the disc is at r = rho t, t = z / position, and at normal
/// incidence its x/X0 is the integral of the LUT's 1/X0 over z. At one t the
/// rays of [axisMin, axisMax] cover [axisMin t, axisMax t], which overlaps a
/// LUT layer [rMin, rMax] over a length linear in t between the points where
/// the ends cross. In a cell the LUT is constant, its azimuths' mean taken,
/// so the cell adds that mean times the integral of overlap(t) / t, a
/// logarithm and a line: a fixed grid misreads the supports' sharp radial
/// edges.
[[nodiscard]] inline double sampledBudget(const o2::base::MatLayerCylSet& lut, bool barrel, double low, double high, double axisMin, double axisMax,
                                          double position)
{
  if (barrel) {
    const double from = std::min(low, high);
    const double to = std::max(low, high);
    double total = 0.;
    for (int i = 0; i < lut.getNLayers(); ++i) {
      const auto& layer = lut.getLayer(i);
      const double pieceLow = std::max(from, static_cast<double>(layer.getRMin()));
      const double pieceHigh = std::min(to, static_cast<double>(layer.getRMax()));
      if (!(pieceHigh > pieceLow)) {
        continue;
      }
      // Along z: the layer's cells within [axisMin, axisMax], and what lies
      // outside its extent as one piece each side, by length.
      const auto zMin = static_cast<double>(layer.getZMin());
      const auto zMax = static_cast<double>(layer.getZMax());
      const auto dz = static_cast<double>(layer.getDZ());
      std::vector<std::pair<double, double>> along;
      for (double at = axisMin; at < axisMax;) {
        double next = axisMax;
        if (at < zMin) {
          next = std::min(axisMax, zMin);
        } else if (at < zMax) {
          next = std::min(axisMax, zMin + (dz * (std::floor((at - zMin) / dz) + 1.)));
        }
        along.emplace_back(0.5 * (at + next), next - at);
        at = next;
      }
      const int phiCells = layer.getNPhiBins();
      double sum = 0.;
      for (int iPhi = 0; iPhi < phiCells; ++iPhi) {
        const double phi = (2 * std::numbers::pi_v<double>)*(static_cast<double>(iPhi) + 0.5) / phiCells;
        const double c = std::cos(phi);
        const double s = std::sin(phi);
        for (const auto& [z, weight] : along) {
          const auto segment = lut.getMatBudget(static_cast<float>(pieceLow * c), static_cast<float>(pieceLow * s), static_cast<float>(z),
                                                static_cast<float>(pieceHigh * c), static_cast<float>(pieceHigh * s), static_cast<float>(z));
          sum += weight * static_cast<double>(segment.meanX2X0);
        }
      }
      total += sum / (static_cast<double>(phiCells) * (axisMax - axisMin));
    }
    return total;
  }
  const double tA = low / position;
  const double tB = high / position;
  const double tFrom = std::min(tA, tB);
  const double tTo = std::max(tA, tB);
  double total = 0.;
  for (int i = 0; i < lut.getNLayers(); ++i) {
    const auto& layer = lut.getLayer(i);
    const auto rMin = static_cast<double>(layer.getRMin());
    const auto rMax = static_cast<double>(layer.getRMax());
    // Where the ends of [axisMin t, axisMax t] cross the layer's radii.
    std::vector<double> crossings{rMin / axisMax, rMax / axisMax};
    if (axisMin > 0.) {
      crossings.push_back(rMin / axisMin);
      crossings.push_back(rMax / axisMin);
    }
    for (int iz = 0; iz < layer.getNZBins(); ++iz) {
      const double zA = static_cast<double>(layer.getZBinMin(iz)) / position;
      const double zB = static_cast<double>(layer.getZBinMax(iz)) / position;
      const double from = std::max(tFrom, std::min(zA, zB));
      const double to = std::min(tTo, std::max(zA, zB));
      if (!(to > from) || axisMax * to <= rMin || axisMin * from >= rMax) {
        continue;
      }
      std::vector<double> edges{from, to};
      for (const double t : crossings) {
        if (t > from && t < to) {
          edges.push_back(t);
        }
      }
      std::ranges::sort(edges);
      double integral = 0.;
      for (std::size_t k = 0; k + 1 < edges.size(); ++k) {
        const double u = edges[k];
        const double v = edges[k + 1];
        const double mid = 0.5 * (u + v);
        // overlap(t) = upper - lower = alpha + beta t on [u, v]
        const bool upperIsLayer = rMax < axisMax * mid;
        const bool lowerIsLayer = rMin > axisMin * mid;
        const double alpha = (upperIsLayer ? rMax : 0.) - (lowerIsLayer ? rMin : 0.);
        const double beta = (upperIsLayer ? 0. : axisMax) - (lowerIsLayer ? 0. : axisMin);
        if (alpha + (beta * mid) <= 0.) {
          continue;
        }
        integral += beta * (v - u);
        if (!exactlyEqual(alpha, 0.)) {
          integral += alpha * std::log(v / u);
        }
      }
      if (!(integral > 0.)) {
        continue;
      }
      double meanX2X0 = 0.;
      for (int iPhi = 0; iPhi < layer.getNPhiBins(); ++iPhi) {
        meanX2X0 += static_cast<double>(layer.getCellPhiBin(iPhi, iz).meanX2X0);
      }
      total += integral * meanX2X0 / layer.getNPhiBins();
    }
  }
  return std::abs(position) * total / (axisMax - axisMin);
}

/// Each layer's average normal-incidence `x/X0` from the LUT.
///
/// Sampled as `sampledBudget` does, over the stave length or the disc's radii.
/// A disc's neighbours are the discs of the same endcap; the innermost takes
/// the same half-gap on both sides.
template <DetectorLayout Det>
[[nodiscard]] inline LayerBudget<Det> layerMaterial(const o2::base::MatLayerCylSet& lut)
{
  constexpr auto kLayers = static_cast<std::size_t>(Detector<Det>::nLayers);

  LayerBudget<Det> budget{};
  for (std::size_t layer = 0; layer < kLayers; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    const bool barrel = Detector<Det>::isBarrelLayer(id);

    // The half-gaps to the neighbours of the same kind and, for discs, side.
    const auto neighbour = [&](int step) -> double {
      const auto other = static_cast<int>(layer) + step;
      if (other < 0 || std::cmp_greater_equal(other, kLayers)) {
        return 0.;
      }
      const auto otherId = static_cast<LayerID>(other);
      if (Detector<Det>::isBarrelLayer(otherId) != barrel) {
        return 0.;
      }
      if (!barrel && discZ<Det>(otherId) * discZ<Det>(id) < 0.) {
        return 0.;
      }
      return barrel ? layerDistance<Det>(otherId) : discZ<Det>(otherId);
    };
    const double here = barrel ? layerDistance<Det>(id) : discZ<Det>(id);
    const double below = neighbour(-1);
    const double above = neighbour(+1);
    double low = !exactlyEqual(below, 0.) ? 0.5 * (below + here) : 0.;
    double high = !exactlyEqual(above, 0.) ? 0.5 * (above + here) : 0.;
    if (barrel) {
      low = !exactlyEqual(below, 0.) ? low : 0.5 * here;
      high = !exactlyEqual(above, 0.) ? high : here + (here - low);
    } else {
      // Discs are numbered by increasing z, so "below" is the smaller z.
      if (exactlyEqual(below, 0.)) {
        low = here - (high - here);
      }
      if (exactlyEqual(above, 0.)) {
        high = here + (here - low);
      }
    }

    // A passive scatterer's band is its own (PassiveSpec): the slab stops at
    // it, so its material is counted once, at its own lever arm.
    if constexpr (Detector<Det>::nPassives > 0) {
      for (std::size_t i = 0; i < static_cast<std::size_t>(Detector<Det>::nPassives); ++i) {
        const auto& passive = Detector<Det>::passive(i);
        if ((passive.kind == PassiveSpec::Kind::kBarrel) != barrel) {
          continue;
        }
        const double slabLow = std::min(low, high);
        const double slabHigh = std::max(low, high);
        if (passive.bandHigh <= here && passive.bandHigh > slabLow) {
          low = passive.bandHigh;
          high = slabHigh;
        } else if (passive.bandLow >= here && passive.bandLow < slabHigh) {
          low = slabLow;
          high = passive.bandLow;
        }
      }
      // The outermost barrel layer's open side mirrors its inner side as the
      // passives left it: where its gap inward is a passive's, its slab is its
      // own silicon, and does not reach material beyond the last layer.
      if (barrel && exactlyEqual(above, 0.)) {
        high = here + (here - low);
      }
    }
    budget[layer] = sampledBudget(lut, barrel, low, high, secondMin<Det>(id), secondMax<Det>(id), here);
  }
  return budget;
}

template <DetectorLayout Det>
using PassiveBudget = std::array<double, static_cast<std::size_t>(Detector<Det>::nPassives)>;

/// Each passive scatterer's nominal `x/X0` (PassiveSpec::x2X0), for a run
/// without a material LUT.
template <DetectorLayout Det>
[[nodiscard]] inline PassiveBudget<Det> nominalPassiveMaterial()
{
  PassiveBudget<Det> budget{};
  if constexpr (Detector<Det>::nPassives > 0) {
    for (std::size_t i = 0; i < budget.size(); ++i) {
      budget[i] = Detector<Det>::passive(i).x2X0;
    }
  }
  return budget;
}

/// Each passive scatterer's average normal-incidence `x/X0` from the LUT:
/// segments across its band, sampled over
/// azimuth and its extent -- radial for a cylinder; for a disc along the ray
/// from the nominal vertex, scaled to normal incidence, as `layerMaterial`
/// reads a disc layer.
template <DetectorLayout Det>
[[nodiscard]] inline PassiveBudget<Det> passiveMaterial(const o2::base::MatLayerCylSet& lut)
{
  PassiveBudget<Det> budget{};
  if constexpr (Detector<Det>::nPassives > 0) {
    for (std::size_t i = 0; i < budget.size(); ++i) {
      const auto& passive = Detector<Det>::passive(i);
      budget[i] =
        sampledBudget(lut, passive.kind == PassiveSpec::Kind::kBarrel, passive.bandLow, passive.bandHigh, passive.low, passive.high, passive.position);
    }
  }
  return budget;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_LAYERMATERIAL_H_

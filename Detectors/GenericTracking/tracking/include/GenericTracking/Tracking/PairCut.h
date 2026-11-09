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

/// \file PairCut.h
/// \brief The exact spatial cuts a tracklet passes.
///
/// Two tests on the measured positions, in the frame of the corrected beam
/// line: a common helix pitch from the vertex to both measurements, and the
/// curvature of the circle through the beam point and both measurements. One
/// expression serves every kind of layer pair, and each measurement's own
/// position and errors enter rather than a layer's nominal ones.
///
/// The cheap transverse test runs first. Only its survivors pay for the
/// division and arcsines of the helix-pitch test. The loop over layers that
/// multiple scattering needs happens once per source and vertex, in
/// `sourceScattering`.
///
/// The iteration's cut and the vertex are separate arguments: one source is
/// tested against every vertex seed compatible with it in time, so the vertex
/// changes within one source's work while the cut stays the iteration's.

#ifndef GTRK_TRACKING_PAIRCUT_H_
#define GTRK_TRACKING_PAIRCUT_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Tracking/SearchWindow.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include "GPUCommonMath.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>

namespace o2::gtrk
{

/// A layer as multiple scattering sees it: where it is, how far it extends,
/// and how much it scatters.
///
/// `k` is the Highland factor at normal incidence, `(0.0136 GeV)^2 x (1 +
/// 0.038 ln x)^2` for the layer's average `x = x/X0`, so that a crossing's
/// `theta0^2` is `k` times the path factor over `p^2`.
struct ScatteringLayer {
  Scalar position = 0.F; ///< barrel radius, or disc z [cm]
  Scalar low = 0.F;      ///< extent along the other axis: stave z, or disc r [cm]
  Scalar high = 0.F;
  Scalar k = 0.F; ///< [GeV^2]
};

/// A passive scatterer as multiple scattering sees it (`PassiveSpec`): a
/// layer's fields, and its kind, which device code cannot read off the
/// layout's constexpr array with a run-time index.
struct PassiveScatterer {
  ScatteringLayer at;
  bool barrel = true;
};

/// Every layer and passive scatterer of `Det`, built once at initialisation
/// from their material budgets; the stages read only this.
template <DetectorLayout Det>
struct ScatteringTable {
  static constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);
  static constexpr std::size_t kPassiveCount = static_cast<std::size_t>(Detector<Det>::nPassives);
  ScatteringLayer layers[kLayerCount]{};
  /// One slot even without passives: a C array cannot be empty.
  PassiveScatterer passives[kPassiveCount > 0 ? kPassiveCount : 1]{};
  static_assert(kPassiveCount <= 64U, "a set of passives is a bit each of a 64-bit mask");
  /// Every passive, a bit each: the set `betweenScattering` takes by default.
  static constexpr std::uint64_t kAllPassives = kPassiveCount == 64U ? ~std::uint64_t{0U} : (std::uint64_t{1U} << kPassiveCount) - 1U;
  /// Per layer, the passives no segment of a layout link into it, or out of
  /// it, can cross, a bit each (`passiveCanLieAcross`):
  /// the ones `visitPassivesBetween` need not try for a triplet's or a pair of
  /// triplets' segments. None, so every passive tried, unless
  /// `scatteringTable` found them.
  std::uint64_t notInto[kLayerCount]{};
  std::uint64_t notOutOf[kLayerCount]{};
  /// The passives a segment of a link into \a layer can cross.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint64_t into(LayerID layer) const { return kAllPassives & ~notInto[layer]; }
  /// The passives a segment of a link out of \a layer can cross.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint64_t outOf(LayerID layer) const
  {
    return kAllPassives & ~notOutOf[layer]; // NOLINT(clang-analyzer-security.ArrayBound): layer < kLayerCount, the geometry's
  }
};

namespace detail
{

/// `(0.0136 GeV)^2 x (1 + 0.038 ln x)^2`: Highland's theta0^2 p^2 at normal
/// incidence for \a x = x/X0, 0 where the logarithm turns it negative.
[[nodiscard]] inline double highlandFactor(double x)
{
  if (!(x > 0.)) {
    return 0.;
  }
  const double highland = 1. + (0.038 * std::log(x));
  return highland > 0. ? 0.0136 * 0.0136 * x * highland * highland : 0.;
}

} // namespace detail

/// Whether a straight segment between a point of layer \a a and one of layer
/// \a b can cross \a passive: a cylinder strictly between the layers' least
/// and largest radius, a disc strictly between their least and largest z --
/// what `visitPassivesBetween` asks of the points' own radius or z.
template <DetectorLayout Det>
[[nodiscard]] inline bool passiveCanLieAcross(const PassiveScatterer& passive, LayerID a, LayerID b)
{
  return positionBetweenLayers<Det>(passive.barrel, static_cast<double>(passive.at.position), a, b);
}

/// The table for per-layer average normal-incidence `x/X0` values and, for the
/// layout's passive scatterers, \a passiveX2X0 (`passiveMaterial`) or, if
/// empty, their nominal values (`nominalPassiveMaterial`).
template <DetectorLayout Det>
[[nodiscard]] inline ScatteringTable<Det> scatteringTable(std::span<const double> x2X0, std::span<const double> passiveX2X0 = {})
{
  ScatteringTable<Det> table;
  for (std::size_t layer = 0; layer < ScatteringTable<Det>::kLayerCount; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    table.layers[layer] = {
      .position = static_cast<Scalar>(Detector<Det>::isBarrelLayer(id) ? layerDistance<Det>(id) : discZ<Det>(id)),
      .low = static_cast<Scalar>(secondMin<Det>(id)),
      .high = static_cast<Scalar>(secondMax<Det>(id)),
      .k = static_cast<Scalar>(detail::highlandFactor(x2X0[layer])),
    };
  }
  if constexpr (Detector<Det>::nPassives > 0) {
    if (!passiveX2X0.empty() && passiveX2X0.size() != ScatteringTable<Det>::kPassiveCount) {
      throw std::invalid_argument("scatteringTable: one x/X0 per passive scatterer");
    }
    for (std::size_t i = 0; i < ScatteringTable<Det>::kPassiveCount; ++i) {
      const auto& spec = Detector<Det>::passive(i);
      table.passives[i] = {
        .at =
          {
            .position = static_cast<Scalar>(spec.position),
            .low = static_cast<Scalar>(spec.low),
            .high = static_cast<Scalar>(spec.high),
            .k = static_cast<Scalar>(detail::highlandFactor(passiveX2X0.empty() ? spec.x2X0 : passiveX2X0[i])),
          },
        .barrel = spec.kind == PassiveSpec::Kind::kBarrel,
      };
    }
    // The passives no segment of a layout link into, or out of, a layer can cross.
    for (std::size_t i = 0; i < ScatteringTable<Det>::kPassiveCount; ++i) {
      std::array<bool, ScatteringTable<Det>::kLayerCount> into{};
      std::array<bool, ScatteringTable<Det>::kLayerCount> outOf{};
      for (const auto& link : Det::layerLinks) {
        if (passiveCanLieAcross<Det>(table.passives[i], link.from, link.to)) {
          into[link.to] = true;
          outOf[link.from] = true;
        }
      }
      for (std::size_t layer = 0; layer < ScatteringTable<Det>::kLayerCount; ++layer) {
        table.notInto[layer] |= into[layer] ? 0U : std::uint64_t{1U} << i;
        table.notOutOf[layer] |= outOf[layer] ? 0U : std::uint64_t{1U} << i;
      }
    }
  }
  return table;
}

/// The cut parameters of one iteration, the same for every layer pair and
/// every vertex.
struct PairCut {
  Scalar twoRMin = 0.F;          ///< twice the softest track's transverse radius [cm]; 0 with no field
  Scalar nSigma = 3.F;           ///< how many sigma of the Gaussian terms are admitted
  Scalar invBetaMomentum2 = 0.F; ///< 1 / (beta pT)^2 at pT_min [GeV^-2], for the scattering sum (`invBetaMomentumMin2`)
  /// The largest transverse curvature above pT_min [1/cm]: zero with no field,
  /// unbounded with no pT cut.
  Scalar maxCurvature = std::numeric_limits<Scalar>::infinity();
  /// A track's momentum times its 3D curvature, |B2C Bz| [GeV/cm]: what turns
  /// a fitted curvature into the momentum its scattering is taken at; 0 with
  /// no field.
  Scalar momentumCurvature = 0.F;
  Scalar mass2 = 0.F; ///< the mass hypothesis squared [GeV^2], for beta
};

[[nodiscard]] inline PairCut pairCut(const SearchWindowParams& params)
{
  const double radius = curvatureRadius(params.minPtGeV, params.bzKiloGauss);
  return {
    .twoRMin = static_cast<Scalar>(2. * radius),
    .nSigma = static_cast<Scalar>(params.nSigma),
    .invBetaMomentum2 = params.minPtGeV > 0. ? static_cast<Scalar>(invBetaMomentumMin2(params)) : 0.F,
    .maxCurvature = static_cast<Scalar>(maxCurvature(params)),
    .momentumCurvature = static_cast<Scalar>(std::abs(static_cast<double>(o2::constants::math::B2C) * params.bzKiloGauss)),
    .mass2 = static_cast<Scalar>(params.massGeV * params.massGeV),
  };
}

/// A candidate origin of the tracks in one collision: the beam line, one
/// entry for the whole timeframe, or a vertex from a list sorted by time.
///
/// Transverse positions are in the frame of the corrected beam line, so a
/// vertex contributes only its transverse spread there.
struct VertexSeed {
  /// The time box `[timePS, timePS + widthPS)` [ps from the timeframe's start]:
  /// a point for an MC vertex, the span of
  /// its contributors' boxes for a found one, every time for the beam line.
  std::uint32_t timePS = 0U;
  std::uint32_t widthPS = 0U;
  Scalar z = 0.F;             ///< [cm]
  Scalar varZ = 0.F;          ///< [cm^2]
  Scalar varTransverse = 0.F; ///< transverse spread squared [cm^2]
};

/// A seed from a spatial constraint alone, compatible with every time: the
/// beam line, from the beam spot's z and widths.
[[nodiscard]] inline VertexSeed vertexSeed(const VertexConstraint& constraint)
{
  return {
    .timePS = 0U,
    .widthPS = std::numeric_limits<std::uint32_t>::max(),
    .z = static_cast<Scalar>(constraint.z),
    .varZ = static_cast<Scalar>(constraint.sigmaZ * constraint.sigmaZ),
    .varTransverse = static_cast<Scalar>(constraint.sigmaTransverse * constraint.sigmaTransverse),
  };
}

/// What the tests read of one measurement: its shipped position and its
/// variances along z, r and r-phi.
struct PairPoint {
  Scalar x = 0.F;
  Scalar y = 0.F;
  Scalar z = 0.F;
  Scalar r = 0.F;
  Scalar varZ = 0.F;
  Scalar varR = 0.F;
  Scalar varRPhi = 0.F;
};

template <typename Row>
[[nodiscard]] GPUhdi() GPUAlwaysInline() PairPoint pairPoint(const Row& row)
{
  return {.x = row.x, .y = row.y, .z = row.z, .r = row.r, .varZ = row.varZ, .varR = row.varR, .varRPhi = row.varRPhi};
}

/// The same for row \a row of \a columns, reading only these seven columns: a
/// whole row through a view's `operator[]` is a call that reads all of them.
template <typename Columns>
[[nodiscard]] GPUhdi() GPUAlwaysInline() PairPoint pairPointAt(const Columns& columns, std::uint32_t row)
{
  return {
    .x = columns.x[row],
    .y = columns.y[row],
    .z = columns.z[row],
    .r = columns.r[row],
    .varZ = columns.varZ[row],
    .varR = columns.varR[row],
    .varRPhi = columns.varRPhi[row],
  };
}

namespace detail
{

/// Whether the line from \a vertex through \a source, of slope \a slope in
/// r and \a dz from the vertex to the source, crosses \a scatterer before
/// the source -- the source's own layer, \a own, always -- and if so its
/// term of `sourceScattering`'s sum in \a term.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool crossedBeforeSource(const ScatteringLayer& scatterer, bool barrel, bool own, const VertexSeed& vertex, const PairPoint& source, Scalar dz, Scalar slope, Scalar absDz, Scalar discFactor, Scalar& term)
{
  if (barrel) {
    if (!own) {
      const Scalar zAt = vertex.z + (slope * scatterer.position);
      if (!(scatterer.position < source.r) || zAt < scatterer.low || zAt > scatterer.high) {
        return false;
      }
    }
    term = scatterer.k * scatterer.position * scatterer.position;
    return true;
  }
  const Scalar dzLayer = scatterer.position - vertex.z;
  if (!own) {
    if (!(dzLayer * dz > 0.F) || !(o2::gpu::GPUCommonMath::Abs(dzLayer) < absDz)) {
      return false;
    }
    const Scalar rAt = source.r * dzLayer / dz;
    if (rAt < scatterer.low || rAt > scatterer.high) {
      return false;
    }
  }
  term = scatterer.k * dzLayer * dzLayer * discFactor;
  return true;
}

} // namespace detail

/// `S = sum_k theta_k^2 s_k^2` over the layers the source's line from the
/// vertex crosses up to and including the source's own, for the iteration's
/// minimum pT.
///
/// A kink at distance `s_k` from the vertex moves the collinearity by
/// `theta_k s_k (s_j - s_i)` and the curvature cross product by
/// `theta_k s_k d`, so this one sum serves both tests. The source's direction
/// from the vertex sets the path through each layer, the lever arm and the
/// momentum `p = pT g`; the Highland logarithm stays at normal incidence.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar sourceScattering(const ScatteringTable<Det>& table, const PairCut& cut, const VertexSeed& vertex, LayerID sourceLayer, const PairPoint& source)
{
  if (!(source.r > 0.F)) {
    return 0.F;
  }
  const Scalar dz = source.z - vertex.z;
  const Scalar slope = dz / source.r;
  const Scalar g = o2::gpu::GPUCommonMath::Sqrt(1.F + (slope * slope));
  const Scalar absDz = o2::gpu::GPUCommonMath::Abs(dz);
  const Scalar invAbsDz = absDz > 0.F ? 1.F / absDz : 0.F;
  const Scalar discFactor = source.r * source.r * source.r * invAbsDz * invAbsDz * invAbsDz;

  // The layers in order, then the passive scatterers inside the source, as
  // layers that are not the source's own: the sum's order is the output's.
  Scalar sum = 0.F;
  Scalar term = 0.F;
  for (std::size_t layer = 0; layer < ScatteringTable<Det>::kLayerCount; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    if (detail::crossedBeforeSource(table.layers[layer], Detector<Det>::isBarrelLayer(id), id == sourceLayer, vertex, source, dz, slope, absDz, discFactor, term)) {
      sum += term;
    }
  }
  if constexpr (ScatteringTable<Det>::kPassiveCount > 0) {
    for (std::size_t i = 0; i < ScatteringTable<Det>::kPassiveCount; ++i) {
      const auto& passive = table.passives[i];
      if (detail::crossedBeforeSource(passive.at, passive.barrel, false, vertex, source, dz, slope, absDz, discFactor, term)) {
        sum += term;
      }
    }
  }
  return sum * g * cut.invBetaMomentum2;
}

/// `T = sum_p theta_p^2 (s_j - s_p)^2` over the passive scatterers the line
/// from the vertex through the source crosses between the source and the
/// target, for the iteration's minimum pT:
/// a kink there moves the target alone, by `theta_p (s_j - s_p)`.
///
/// As in `sourceScattering`, the source's direction from the vertex sets the
/// path through the scatterer, its distance `s_p` and the momentum `p = pT
/// g`; `s_j` is the target's own distance from the vertex. A layout's layers
/// never enter: their material is the target's, which does not matter, or
/// the source's, which `sourceScattering` has.
///
/// \a passives: the ones to try, a bit each, in index order; every passive
/// by default. A link's `passivesBetween` holds every passive any of its pairs
/// can cross, so the terms and their order, and
/// the sum, are the same as with all of them.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() Scalar betweenScattering(const ScatteringTable<Det>& table, const PairCut& cut, const VertexSeed& vertex, const PairPoint& source, const PairPoint& target, std::uint64_t passives = ScatteringTable<Det>::kAllPassives)
{
  if constexpr (ScatteringTable<Det>::kPassiveCount == 0) {
    return 0.F;
  } else {
    if (!(source.r > 0.F)) {
      return 0.F;
    }
    const Scalar dz = source.z - vertex.z;
    const Scalar slope = dz / source.r;
    const Scalar g = o2::gpu::GPUCommonMath::Sqrt(1.F + (slope * slope));
    const Scalar absDz = o2::gpu::GPUCommonMath::Abs(dz);
    const Scalar dzTarget = target.z - vertex.z;
    const Scalar sTarget = o2::gpu::GPUCommonMath::Sqrt((target.r * target.r) + (dzTarget * dzTarget));
    Scalar sum = 0.F;
    for (std::uint64_t left = passives & ScatteringTable<Det>::kAllPassives; left != 0U; left &= left - 1U) {
      const auto& passive = table.passives[o2::gpu::CAMath::Ctz64(left)];
      const auto& scatterer = passive.at;
      Scalar pathFactor = 0.F;
      Scalar sPassive = 0.F;
      if (passive.barrel) {
        const Scalar zAt = vertex.z + (slope * scatterer.position);
        if (scatterer.position <= source.r || scatterer.position >= target.r || zAt < scatterer.low || zAt > scatterer.high) {
          continue;
        }
        pathFactor = g;
        sPassive = scatterer.position * g;
      } else {
        const Scalar dzPassive = scatterer.position - vertex.z;
        if (!(dzPassive * dz > 0.F) || !(o2::gpu::GPUCommonMath::Abs(dzPassive) > absDz) || !(o2::gpu::GPUCommonMath::Abs(dzPassive) < o2::gpu::GPUCommonMath::Abs(dzTarget)) || !(absDz > 0.F)) {
          continue;
        }
        const Scalar rAt = source.r * dzPassive / dz;
        if (rAt < scatterer.low || rAt > scatterer.high) {
          continue;
        }
        pathFactor = source.r * g / absDz;
        sPassive = o2::gpu::GPUCommonMath::Abs(dzPassive) * pathFactor;
      }
      // theta_p^2 pT^2 = k (path factor) / g^2, with p = pT g.
      const Scalar lever = sTarget - sPassive;
      sum += scatterer.k * pathFactor / (g * g) * lever * lever;
    }
    return sum * cut.invBetaMomentum2;
  }
}

/// One standard deviation of a source's transverse direction on \a from
/// [rad], at its largest over the layer and the vertex range: `sqrt(S)/r_i`,
/// with `S` the `sourceScattering` sum at pT_min, and the beam width over the
/// radius. The transverse test takes both times the pair's chord.
///
/// Each term of the sum is bounded where the source can give it its largest
/// value. A barrel layer inside the source contributes `k (r_l/r_i)^2 g` and
/// takes the steepest path factor `g`. A disc crossed on the way contributes
/// `k dz_l^2 sqrt(r_i^2 + dz^2) / |dz|^3`, which falls with the source's
/// distance `|dz| > |dz_l|` from the vertex and is largest just beyond the
/// disc; the source's own disc takes its nearest and farthest offsets.
template <DetectorLayout Det>
[[nodiscard]] inline double sourceKinkBound(const ScatteringTable<Det>& table, const SearchWindowParams& params, LayerID from, const VertexConstraint& vertex)
{
  const double rMin = radiusMin<Det>(from);
  const double rMax = radiusMax<Det>(from);
  if (!(rMin > 0.) || !(params.minPtGeV > 0.)) {
    return std::numeric_limits<double>::infinity();
  }
  const Interval vertexZ = vertexReach(params, vertex);
  const Interval sourceZ = zRange<Det>(from);
  // The distance from the vertex range to a position: the nearest, zero inside it.
  const auto nearest = [&](double z) { return nearestApart(vertexZ, {.low = z, .high = z}); };
  const double sourceFar = farthestApart(vertexZ, sourceZ);
  const double g = std::sqrt(1. + (sourceFar * sourceFar / (rMin * rMin)));

  double sum = 0.;
  // A scatterer that is not the source's own layer, at its largest: false
  // where that is unbounded.
  const auto addOther = [&](double k, double position, bool barrel) {
    if (barrel) {
      if (position < rMax) {
        const double ratio = std::min(position / rMin, 1.);
        sum += k * ratio * ratio * g;
      }
      return true;
    }
    // Crossed only between the vertex and the source.
    const bool crossed = (position > vertexZ.low && position < sourceZ.high) || (position < vertexZ.high && position > sourceZ.low);
    if (!crossed) {
      return true;
    }
    const double near = nearest(position);
    if (!(near > 0.)) {
      return false;
    }
    sum += k * std::sqrt((rMax * rMax) + (near * near)) / near;
    return true;
  };
  for (std::size_t layer = 0; layer < ScatteringTable<Det>::kLayerCount; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    const auto k = static_cast<double>(table.layers[layer].k);
    const auto position = static_cast<double>(table.layers[layer].position);
    const bool barrel = Detector<Det>::isBarrelLayer(id);
    if (id != from) {
      if (!addOther(k, position, barrel)) {
        return std::numeric_limits<double>::infinity();
      }
      continue;
    }
    if (barrel) {
      const double ratio = position / rMin;
      sum += k * ratio * ratio * g;
      continue;
    }
    const double sourceNear = nearestApart(vertexZ, sourceZ);
    if (!(sourceNear > 0.)) {
      return std::numeric_limits<double>::infinity();
    }
    const double offset = farthestApart(vertexZ, {.low = position, .high = position});
    sum += k * offset * offset * std::sqrt((rMax * rMax) + (sourceNear * sourceNear)) / (sourceNear * sourceNear * sourceNear);
  }
  // Passive scatterers that can lie inside a source on this layer, bounded as
  // a layer that is not the source's own.
  if constexpr (ScatteringTable<Det>::kPassiveCount > 0) {
    for (std::size_t i = 0; i < ScatteringTable<Det>::kPassiveCount; ++i) {
      const auto& passive = table.passives[i];
      if (!addOther(static_cast<double>(passive.at.k), static_cast<double>(passive.at.position), passive.barrel)) {
        return std::numeric_limits<double>::infinity();
      }
    }
  }
  const double beam = vertex.sigmaTransverse / rMin;
  return std::sqrt((sum / (params.minPtGeV * params.minPtGeV) * betaMomentumFactor(params)) + (beam * beam));
}

/// Whether \a passive can lie between a pair on the link \a from -> \a to:
/// a cylinder between the source layer's least and the target's largest
/// radius, a disc between the layers' z on the target's side.
template <DetectorLayout Det>
[[nodiscard]] inline bool passiveCanLieBetween(const PassiveScatterer& passive, LayerID from, LayerID to)
{
  const auto position = static_cast<double>(passive.at.position);
  if (passive.barrel) {
    return position > radiusMin<Det>(from) && position < radiusMax<Det>(to);
  }
  // Beyond some of the source layer's z and short of some of the target's, on one side.
  const double fromNear = std::min(std::abs(zMin<Det>(from)), std::abs(zMax<Det>(from)));
  const double toFar = std::max(std::abs(zMin<Det>(to)), std::abs(zMax<Det>(to)));
  return std::abs(position) > fromNear && std::abs(position) < toFar && !(position * zMax<Det>(to) < 0.);
}

/// The passives a pair on the link \a from -> \a to can cross, a bit each:
/// those `passiveCanLieBetween` the layers, the test `betweenKinkBound` bounds
/// over. `betweenScattering` on them adds the terms it adds on all.
template <DetectorLayout Det>
[[nodiscard]] inline std::uint64_t passivesBetween(const ScatteringTable<Det>& table, LayerID from, LayerID to)
{
  std::uint64_t bits = 0U;
  for (std::size_t i = 0; i < ScatteringTable<Det>::kPassiveCount; ++i) {
    if (passiveCanLieBetween<Det>(table.passives[i], from, to)) {
      bits |= std::uint64_t{1U} << i;
    }
  }
  return bits;
}

/// A bound on `sqrt(T)` [cm] (`betweenScattering`) for every pair on the link
/// \a from -> \a to: the passive scatterers that can lie between the two
/// layers, each at its largest -- `theta_p^2 pT^2 <= k_p` for both kinds, the
/// path factor over `g^2` never exceeding one -- over the longest lever from
/// the scatterer to the target. Zero if none can, which is how an iteration
/// knows to skip the term (`IterationParameters::passiveBetween`).
template <DetectorLayout Det>
[[nodiscard]] inline double betweenKinkBound(const ScatteringTable<Det>& table, const SearchWindowParams& params, LayerID from, LayerID to, const VertexConstraint& vertex)
{
  if constexpr (ScatteringTable<Det>::kPassiveCount == 0) {
    return 0.;
  } else {
    if (!(params.minPtGeV > 0.)) {
      return 0.;
    }
    const Interval vertexZ = vertexReach(params, vertex);
    const auto nearest = [&](double z) { return nearestApart(vertexZ, {.low = z, .high = z}); };
    const double targetFar = std::hypot(radiusMax<Det>(to), farthestApart(vertexZ, zRange<Det>(to)));
    double sum = 0.;
    for (std::size_t i = 0; i < ScatteringTable<Det>::kPassiveCount; ++i) {
      const auto& passive = table.passives[i];
      if (!passiveCanLieBetween<Det>(passive, from, to)) {
        continue;
      }
      const auto position = static_cast<double>(passive.at.position);
      const double sPassive = passive.barrel ? position : nearest(position);
      const double lever = std::max(targetFar - sPassive, 0.);
      sum += static_cast<double>(passive.at.k) * lever * lever;
    }
    return std::sqrt(sum * betaMomentumFactor(params)) / params.minPtGeV;
  }
}

/// What `betweenScattering` computes of one source and vertex alone: the
/// passives the source's side admits, each one's weight `k pathFactor / g^2`
/// and `sPassive`, and the bound a target has to pass (a barrel's radius, a
/// disc's |z - z_v|). A pair then adds only its own side, the same terms in
/// the same order with the same arithmetic.
template <DetectorLayout Det>
struct BetweenSource {
  static constexpr std::size_t kCount = ScatteringTable<Det>::kPassiveCount > 0 ? ScatteringTable<Det>::kPassiveCount : 1;
  std::uint64_t admitted = 0U;
  std::uint64_t barrels = 0U;
  Scalar weight[kCount]{};
  Scalar sPassive[kCount]{};
  Scalar bound[kCount]{};
};

/// The source's part of `betweenScattering` through \a vertex.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() BetweenSource<Det> betweenSource(const ScatteringTable<Det>& table, const VertexSeed& vertex, const PairPoint& source)
{
  BetweenSource<Det> result;
  if constexpr (ScatteringTable<Det>::kPassiveCount > 0) {
    if (!(source.r > 0.F)) {
      return result;
    }
    const Scalar dz = source.z - vertex.z;
    const Scalar slope = dz / source.r;
    const Scalar g = o2::gpu::GPUCommonMath::Sqrt(1.F + (slope * slope));
    const Scalar absDz = o2::gpu::GPUCommonMath::Abs(dz);
    for (std::size_t i = 0; i < ScatteringTable<Det>::kPassiveCount; ++i) {
      const auto& passive = table.passives[i];
      const auto& scatterer = passive.at;
      Scalar pathFactor = 0.F;
      Scalar sPassive = 0.F;
      if (passive.barrel) {
        const Scalar zAt = vertex.z + (slope * scatterer.position);
        if (scatterer.position <= source.r || zAt < scatterer.low || zAt > scatterer.high) {
          continue;
        }
        pathFactor = g;
        sPassive = scatterer.position * g;
        result.bound[i] = scatterer.position;
        result.barrels |= std::uint64_t{1U} << i;
      } else {
        const Scalar dzPassive = scatterer.position - vertex.z;
        if (!(dzPassive * dz > 0.F) || !(o2::gpu::GPUCommonMath::Abs(dzPassive) > absDz) || !(absDz > 0.F)) {
          continue;
        }
        const Scalar rAt = source.r * dzPassive / dz;
        if (rAt < scatterer.low || rAt > scatterer.high) {
          continue;
        }
        pathFactor = source.r * g / absDz;
        sPassive = o2::gpu::GPUCommonMath::Abs(dzPassive) * pathFactor;
        result.bound[i] = o2::gpu::GPUCommonMath::Abs(dzPassive);
      }
      result.admitted |= std::uint64_t{1U} << i;
      result.weight[i] = scatterer.k * pathFactor / (g * g);
      result.sPassive[i] = sPassive;
    }
  }
  return result;
}

/// `betweenScattering` of the pair from \a from, its source's part, to
/// \a target, over \a passives (`IterationParameters::betweenPassives`).
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar betweenScattering(const BetweenSource<Det>& from, const PairCut& cut, const VertexSeed& vertex, const PairPoint& target, std::uint64_t passives)
{
  if constexpr (ScatteringTable<Det>::kPassiveCount == 0) {
    return 0.F;
  } else {
    const Scalar dzTarget = target.z - vertex.z;
    const Scalar absDzTarget = o2::gpu::GPUCommonMath::Abs(dzTarget);
    const Scalar sTarget = o2::gpu::GPUCommonMath::Sqrt((target.r * target.r) + (dzTarget * dzTarget));
    Scalar sum = 0.F;
    for (std::uint64_t left = passives & from.admitted; left != 0U; left &= left - 1U) {
      const std::uint32_t i = o2::gpu::CAMath::Ctz64(left);
      // A barrel's radius below the target's, a disc's |z - z_v| below the target's.
      if (((from.barrels >> i) & 1U) != 0U ? from.bound[i] >= target.r : !(from.bound[i] < absDzTarget)) {
        continue;
      }
      const Scalar lever = sTarget - from.sPassive[i];
      sum += from.weight[i] * lever * lever;
    }
    return sum * cut.invBetaMomentum2;
  }
}

/// The window a pair query uses: its phi and second-coordinate half-widths
/// hold what the pair tests accept with \a table's scattering, not only the
/// uniform material the grid's cells are sized with.
template <DetectorLayout Det>
[[nodiscard]] inline LayerPairWindow layerPairWindow(LayerID from, LayerID to, const SearchWindowParams& params, const VertexConstraint& vertex, const ScatteringTable<Det>& table)
{
  return layerPairWindow<Det>(from, to, secondMax<Det>(from), params, vertex, sourceKinkBound<Det>(table, params, from, vertex), betweenKinkBound<Det>(table, params, from, to, vertex));
}

/// The two tests' ingredients for one source, vertex and candidate.
struct PairTerms {
  bool beforeTurn = false;      ///< the circle's turning point is not between the two
  Scalar curvature = 0.F;       ///< 2 R_min |P_i x P_j|
  Scalar bound = 0.F;           ///< r_i r_j d: the curvature bound, in the same units
  Scalar varCurvature = 0.F;    ///< (2 R_min sigma_X)^2
  Scalar collinearity2 = 0.F;   ///< D^2
  Scalar varCollinearity = 0.F; ///< sigma_D^2
};

namespace detail
{

struct ArcLengthDerivatives {
  Scalar value;
  Scalar chord;
  Scalar curvature;
};

/// Arc length and its derivatives with respect to chord and signed curvature.
[[nodiscard]] GPUhdi() GPUAlwaysInline() ArcLengthDerivatives arcLengthDerivatives(Scalar curvature, Scalar chord)
{
  const Scalar halfTurnSine = o2::gpu::GPUCommonMath::Clamp(0.5F * curvature * chord, -1.F, 1.F);
  const auto quotient = asinOverArgumentWithDerivative(halfTurnSine);
  const Scalar ratio = quotient.value;
  const Scalar derivative = quotient.derivative;
  return {.value = chord * ratio, .chord = ratio + (halfTurnSine * derivative), .curvature = 0.5F * chord * chord * derivative};
}

struct ArcPairGeometry {
  Scalar curvature;
  Scalar curvatureScale;
  Scalar chord;
  bool bounded; ///< the curvature is the cut's bound, which the positions do not move
  ArcLengthDerivatives source;
  ArcLengthDerivatives target;
};

/// The circle through the beam point and the pair, with its curvature held to
/// what a track above pT_min can have. The positions can put it further out
/// -- up to 2/r when the chord is unresolved -- and the arc length there, and
/// its derivatives, would describe a track the cut does not admit. The bound
/// never exceeds a curvature the measured one already reached, so the arcsine
/// argument stays within [-1, 1]. Arc length is even in the curvature, so the
/// bound's sign does not matter.
[[nodiscard]] GPUhdi() GPUAlwaysInline() ArcPairGeometry arcPairGeometry(const PairPoint& source, const PairPoint& target, Scalar cross, Scalar chord, Scalar maxCurvature)
{
  const Scalar curvatureScale = 2.F / (source.r * target.r * chord);
  const Scalar measured = curvatureScale * cross;
  // Also false for the NaN of a zero chord.
  const bool bounded = !(o2::gpu::GPUCommonMath::Abs(measured) <= maxCurvature);
  const Scalar curvature = bounded ? maxCurvature : measured;
  return {
    .curvature = curvature,
    .curvatureScale = curvatureScale,
    .chord = chord,
    .bounded = bounded,
    .source = arcLengthDerivatives(curvature, source.r),
    .target = arcLengthDerivatives(curvature, target.r),
  };
}

} // namespace detail

/// What the pair tests compute of the two positions alone: the transverse
/// chord and the cross product, once a pair for all of them.
struct PairChord {
  Scalar dx;
  Scalar dy;
  Scalar chord2;
  Scalar chord;
  Scalar cross;
};

[[nodiscard]] GPUhdi() GPUAlwaysInline() PairChord pairChord(const PairPoint& source, const PairPoint& target)
{
  const Scalar dx = target.x - source.x;
  const Scalar dy = target.y - source.y;
  const Scalar chord2 = (dx * dx) + (dy * dy);
  return {.dx = dx, .dy = dy, .chord2 = chord2, .chord = o2::gpu::GPUCommonMath::Sqrt(chord2), .cross = (source.x * target.y) - (target.x * source.y)};
}

/// \a betweenScatter is `betweenScattering` for the pair, where a passive
/// scatterer can lie between its layers, and zero otherwise.
///
/// Always inlined, as `pairCompatible`: the tracklet stage's grid walk runs
/// the transverse test on every candidate, and GCC keeps the two out of
/// line once the walk grows past its inlining limits -- a call each per
/// candidate.
[[nodiscard]] GPUhdi() GPUAlwaysInline() PairTerms transversePairTerms(const PairCut& cut, const VertexSeed& vertex, const PairPoint& source, Scalar sourceScatter, const PairPoint& target, const PairChord& pair, Scalar betweenScatter = 0.F)
{
  const Scalar chord2 = pair.chord2;
  const Scalar chord = pair.chord;
  const Scalar cross = pair.cross;
  // A kink between the two moves P_j against the side of length r_i.
  const Scalar varCross = (source.r * source.r * (target.varRPhi + betweenScatter)) + (target.r * target.r * source.varRPhi) + (chord2 * (vertex.varTransverse + sourceScatter));

  // The circle through the beam point turns where it is farthest from it,
  // opposite the beam point across its centre. The pair lies on one side of
  // that diameter -- no turning point between them -- exactly when
  // `P_i . P_j` lies between `r_i^2` and `r_j^2`: the step from one to the
  // other leads outwards at both ends. That also keeps them in one hemisphere.
  const Scalar dot = (source.x * target.x) + (source.y * target.y);
  const Scalar sourceR2 = source.r * source.r;
  const Scalar targetR2 = target.r * target.r;
  const Scalar lowR2 = sourceR2 < targetR2 ? sourceR2 : targetR2;
  const Scalar highR2 = sourceR2 < targetR2 ? targetR2 : sourceR2;

  return {
    .beforeTurn = dot >= lowR2 && dot <= highR2,
    .curvature = cut.twoRMin * o2::gpu::GPUCommonMath::Abs(cross),
    .bound = source.r * target.r * chord,
    .varCurvature = cut.twoRMin * cut.twoRMin * varCross,
  };
}

GPUhdi() void addLongitudinalPairTerms(PairTerms& terms, const PairCut& cut, const VertexSeed& vertex, const PairPoint& source, Scalar sourceScatter, const PairPoint& target, const PairChord& pair, Scalar betweenScatter = 0.F)
{
  const Scalar dx = pair.dx;
  const Scalar dy = pair.dy;
  const Scalar chord = pair.chord;
  const Scalar cross = pair.cross;
  const auto geometry = detail::arcPairGeometry(source, target, cross, chord, cut.maxCurvature);
  const Scalar dzSource = source.z - vertex.z;
  const Scalar dzTarget = target.z - vertex.z;
  const Scalar collinearity = (geometry.source.value * dzTarget) - (geometry.target.value * dzSource);

  const Scalar invSourceR = 1.F / source.r;
  const Scalar invTargetR = 1.F / target.r;
  // A bounded curvature does not follow the positions: its terms drop out.
  const Scalar curvatureD = geometry.bounded ? 0.F : (dzTarget * geometry.source.curvature) - (dzSource * geometry.target.curvature);
  const Scalar sourceRadiusD = (dzTarget * geometry.source.chord) - (curvatureD * geometry.curvature * invSourceR);
  const Scalar targetRadiusD = (-dzSource * geometry.target.chord) - (curvatureD * geometry.curvature * invTargetR);
  const Scalar crossD = geometry.bounded ? 0.F : curvatureD * geometry.curvatureScale;
  const Scalar chordD = geometry.bounded ? 0.F : curvatureD * geometry.curvature / (geometry.chord * geometry.chord);

  const Scalar sourceX = (sourceRadiusD * source.x * invSourceR) + (crossD * target.y) + (chordD * dx);
  const Scalar sourceY = (sourceRadiusD * source.y * invSourceR) - (crossD * target.x) + (chordD * dy);
  const Scalar targetX = (targetRadiusD * target.x * invTargetR) - (crossD * source.y) - (chordD * dx);
  const Scalar targetY = (targetRadiusD * target.y * invTargetR) + (crossD * source.x) - (chordD * dy);
  const Scalar sourceRadialD = (sourceX * source.x + sourceY * source.y) * invSourceR;
  const Scalar sourceTangentD = (-sourceX * source.y + sourceY * source.x) * invSourceR;
  const Scalar targetRadialD = (targetX * target.x + targetY * target.y) * invTargetR;
  const Scalar targetTangentD = (-targetX * target.y + targetY * target.x) * invTargetR;
  const Scalar vertexXD = -(sourceX + targetX);
  const Scalar vertexYD = -(sourceY + targetY);
  const Scalar arcGap = geometry.target.value - geometry.source.value;
  const Scalar dz = target.z - source.z;

  terms.collinearity2 = collinearity * collinearity;
  terms.varCollinearity = (source.varR * o2::gpu::CAMath::Square(sourceRadialD)) + (o2::gpu::CAMath::Square(sourceTangentD) * source.varRPhi) + (o2::gpu::CAMath::Square(targetRadialD) * target.varR) + (o2::gpu::CAMath::Square(targetTangentD) * target.varRPhi) + ((o2::gpu::CAMath::Square(vertexXD) + o2::gpu::CAMath::Square(vertexYD)) * vertex.varTransverse) + (o2::gpu::CAMath::Square(geometry.target.value) * source.varZ) + (o2::gpu::CAMath::Square(geometry.source.value) * target.varZ) + (o2::gpu::CAMath::Square(arcGap) * vertex.varZ) + ((o2::gpu::CAMath::Square(arcGap) + o2::gpu::CAMath::Square(dz)) * sourceScatter) +
                          // A kink between the two moves the target's (s, z) corner against (s_i, z_i - z_v).
                          ((o2::gpu::CAMath::Square(geometry.source.value) + o2::gpu::CAMath::Square(dzSource)) * betweenScatter);
}

/// Curvature at most the iteration's bound plus n sigma, with no turning point
/// between the two: a looper's legs are found separately. With no field
/// `twoRMin` is zero and every such pair passes.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool transverseCompatible(const PairCut& cut, const PairTerms& terms)
{
  const Scalar excess = terms.curvature - terms.bound;
  return terms.beforeTurn && (!(excess > 0.F) || excess * excess <= cut.nSigma * cut.nSigma * terms.varCurvature);
}

/// Collinear with the vertex within n sigma.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool longitudinalCompatible(const PairCut& cut, const PairTerms& terms) { return terms.collinearity2 <= cut.nSigma * cut.nSigma * terms.varCollinearity; }

/// Always inlined (see `transversePairTerms`); the longitudinal half stays a
/// call.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool pairCompatible(const PairCut& cut, const VertexSeed& vertex, const PairPoint& source, Scalar sourceScatter, const PairPoint& target, const PairChord& pair, Scalar betweenScatter = 0.F)
{
  countWork(WorkCounter::kPairTests);
  auto terms = transversePairTerms(cut, vertex, source, sourceScatter, target, pair, betweenScatter);
  if (!transverseCompatible(cut, terms)) {
    return false;
  }
  countWork(WorkCounter::kPairTransverse);
  addLongitudinalPairTerms(terms, cut, vertex, source, sourceScatter, target, pair, betweenScatter);
  const bool compatible = longitudinalCompatible(cut, terms);
  if (compatible) {
    countWork(WorkCounter::kPairLongitudinal);
  }
  return compatible;
}

/// The same three without a computed chord, one pair each: the harness's and
/// the tests' (StageFunnel, the sink).
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool pairCompatible(const PairCut& cut, const VertexSeed& vertex, const PairPoint& source, Scalar sourceScatter, const PairPoint& target, Scalar betweenScatter = 0.F) { return pairCompatible(cut, vertex, source, sourceScatter, target, pairChord(source, target), betweenScatter); }

[[nodiscard]] GPUhdi() GPUAlwaysInline() PairTerms transversePairTerms(const PairCut& cut, const VertexSeed& vertex, const PairPoint& source, Scalar sourceScatter, const PairPoint& target, Scalar betweenScatter = 0.F) { return transversePairTerms(cut, vertex, source, sourceScatter, target, pairChord(source, target), betweenScatter); }

GPUhdi() GPUAlwaysInline() void addLongitudinalPairTerms(PairTerms& terms, const PairCut& cut, const VertexSeed& vertex, const PairPoint& source, Scalar sourceScatter, const PairPoint& target, Scalar betweenScatter = 0.F) { addLongitudinalPairTerms(terms, cut, vertex, source, sourceScatter, target, pairChord(source, target), betweenScatter); }

} // namespace o2::gtrk

#endif // GTRK_TRACKING_PAIRCUT_H_

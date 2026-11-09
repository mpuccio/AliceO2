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

/// \file LayerSpec.h
/// \brief The compile-time half of a detector description.
///
/// A detector is a type, not a file: a struct of static constexpr members
/// satisfying the DetectorLayout concept below.  The layer count, the
/// barrel/disc split, each layer's sensor count and sensor shape are fixed at
/// compile time, which is what lets `if constexpr (Detector<Det>::hasDiscs)`
/// remove the forward tracking paths entirely from a barrel-only detector.
/// Where each sensor sits is runtime data, the experiment's geometry -- see
/// Geometry.h.
///
/// Barrel and disc layers are kept in two separate arrays rather than one
/// tagged array on purpose.  It makes a detector without discs express that as
/// an empty array, and it makes the loops over each kind naturally separable,
/// so the eliminated branch does not have to be found by the optimiser.
///
/// Units are centimetres and radians, following O2 and TGeo, and every
/// geometric quantity is double.  Geometry is built once, so the precision is
/// free here, and it keeps alignment clear of rounding that has to be reasoned
/// about.  The tracking narrows to float at its own boundary, where bandwidth
/// and registers actually matter; O2 draws the same line, with DetMatrixCache
/// holding double Transform3D but float Rotation2D.

#ifndef GTRK_BASE_LAYERSPEC_H_
#define GTRK_BASE_LAYERSPEC_H_

#include "GenericTracking/Base/Definitions.h"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <utility>

namespace o2::gtrk
{

/// Outline of a sensor in its own plane -- or, for a curved sensor, in the
/// plane its surface unrolls into.
enum class ShapeType : std::uint8_t {
  kRectangle, ///< halfX applies along the whole length
  kCylinder   ///< a rectangle rolled onto a cylinder; see SensorSpec::radius
};

/// The sensor's local frame is right-handed with
///   x across the sensor (pixel rows),
///   y along it (pixel columns),
///   z normal to it (thickness).
/// Every sensor of a layer -- or of a disc ring -- shares one of these.
struct SensorSpec {
  ShapeType shape;
  double halfX;         ///< half width
  double halfY;         ///< half length
  double halfThickness; ///< half the active silicon thickness
  double pitchX;        ///< pixel pitch across the sensor
  double pitchY;        ///< pixel pitch along the sensor
  /// kCylinder only: the radius of the sensor's inner surface about its local z
  /// axis, as O2 measures a curved sensor's flat coordinate on it
  /// (TRKBase/SegmentationChip.h, `curvedToFlat`). The unrolled x is the arc
  /// length at this radius from the local +x axis, counter-clockwise, and y is
  /// the local z; the silicon runs from here outwards over 2 halfThickness.
  double radius = 0.;

  [[nodiscard]] constexpr bool curved() const { return shape == ShapeType::kCylinder; }

  /// A point of the sensor in its own frame, from its unrolled coordinates:
  /// \a x across, \a y along, \a z off the surface. A flat sensor's frame is
  /// those coordinates; a curved one's is the cylinder's, the point at angle
  /// x / radius from the local +x axis and radius + z from the local z axis.
  void unrolledToFrame(double x, double y, double z, double frame[3]) const
  {
    if (!curved()) {
      frame[0] = x;
      frame[1] = y;
      frame[2] = z;
      return;
    }
    const double angle = x / radius;
    frame[0] = (radius + z) * std::cos(angle);
    frame[1] = (radius + z) * std::sin(angle);
    frame[2] = y;
  }

  /// The range of the unrolled z the silicon spans: +-halfThickness about a
  /// flat sensor's plane, 2 halfThickness outwards from a curved one's radius.
  [[nodiscard]] constexpr double surfaceMin() const { return curved() ? 0. : -halfThickness; }
  [[nodiscard]] constexpr double surfaceMax() const { return curved() ? 2. * halfThickness : halfThickness; }
};

/// Where a layer's silicon is, as the measurement grid, the search windows
/// and the material sampling ask for it: the nominal distance from the
/// interaction point and the hull of the active silicon in r and z.
///
/// Measured from the geometry the layout was loaded with (LayerEnvelopes.h):
/// where its chips sit is run-time data, aligned as the geometry is, and
/// nothing compiled states or checks it.
struct LayerEnvelope {
  double position; ///< barrel: nominal radius; disc: signed z [cm]
  double rMin;     ///< [cm]
  double rMax;
  double zMin;
  double zMax;
};

/// A layer whose sensors are placed from outside, from an experiment's own
/// geometry, as O2's GeometryTGeo places the ALICE ITS and MFT chips.
///
/// Only topology is compile-time: the sensor count, which sizes the geometry
/// and tells a wrong file apart, and the sensor shape. The layer's envelope is
/// measured when the geometry is loaded (LayerEnvelopes.h).
///
/// A layout uses this type for both its barrel and its disc arrays. Barrel
/// and disc still differ in which second coordinate a layer is binned in, and
/// that is all the tracking asks of the kind.
struct PlacedLayerSpec {
  int sensorCount;
  SensorSpec sensor;

  [[nodiscard]] constexpr int nSensors() const { return sensorCount; }
};

/// One directed link of a detector's layer graph: a track from the luminous
/// region crosses `from` and then `to`, with no active layer in between.
///
/// The links are part of the detector's specification, like the ring radii
/// they follow from. `o2-gtrk-layer-links` traces them from the layout once
/// and prints the table to paste in; nothing derives them at run time.
///
/// The source interval is the range of the seed's second coordinate -- z on a
/// barrel layer, r on a disc -- from which the trace reaches `to`, padded by
/// the margin the tool was given. A seed outside it does not query the link.
/// `skipsLayer` is traced too: true when every traced trajectory on the link
/// crosses another layer's envelope between its ends, so no pair on it can be
/// direct. On the other links each pair's own segment decides
/// (Tracking/LinkHoles.h).
struct LayerLink {
  LayerID from;
  LayerID to;
  double sourceMin; ///< [cm]
  double sourceMax; ///< [cm]
  bool skipsLayer;
};

/// Passive material that scatters but measures nothing -- the beam pipe first
/// -- placed by its own position among the measurement layers, which may lie
/// on both sides of it.
///
/// A cylinder or a disc. `position` is where its kink is taken, the radius or
/// the signed z; `bandLow`/`bandHigh` are the stretch of that axis whose
/// material it owns -- a layer's material slab stops there, so nothing is
/// counted twice; `low`/`high` its extent along the other axis. `x2X0` is the
/// nominal normal-incidence budget, for a run without a material LUT; with
/// one, the LUT's across the band stands in.
struct PassiveSpec {
  enum class Kind : std::uint8_t { kBarrel,
                                   kDisc };
  Kind kind;
  double position; ///< barrel: radius; disc: signed z [cm]
  double bandLow;  ///< [cm], along the same axis as position
  double bandHigh;
  double low; ///< extent along the other axis: z for a barrel, r for a disc [cm]
  double high;
  double x2X0; ///< nominal normal-incidence x/X0
};

/// O2's beryllium beam pipe, as the ALICE ITS and MFT see it (O2's
/// Detectors/Passive/src/Pipe.cxx): outer radius 1.9 cm, 0.8 mm thick, over
/// |z| < 44.4 cm; 0.08 cm over beryllium's X0 of 35.28 cm.
inline constexpr PassiveSpec kAliceBerylliumPipe{.kind = PassiveSpec::Kind::kBarrel, .position = 1.86, .bandLow = 1.8, .bandHigh = 1.95, .low = -44.4, .high = 44.4, .x2X0 = 0.08 / 35.28};

/// What a detector type has to provide.
///
/// Derived quantities -- layer counts, hasBarrel, hasDiscs -- are deliberately
/// not required members: they are computed below, so a detector cannot state
/// one of them inconsistently with its own arrays. What follows from the link
/// table is in LayerLinks.h.
template <typename Det>
concept DetectorLayout = requires {
  { Det::name } -> std::convertible_to<const char*>;
  { Det::nSensors } -> std::convertible_to<int>;
  { Det::layerLinks[0] } -> std::convertible_to<const LayerLink&>;
  { Det::barrelLayers[0] } -> std::convertible_to<const PlacedLayerSpec&>;
  { Det::discLayers[0] } -> std::convertible_to<const PlacedLayerSpec&>;
};

/// Compile-time queries derived from one detector layout.
///
/// Keeping them on one traits type gives tracking code a detector-level API
/// without requiring each layout to repeat facts that follow from its arrays.
template <DetectorLayout Layout>
struct Detector {
  static constexpr int nBarrelLayers = static_cast<int>(std::size(Layout::barrelLayers));
  static constexpr int nDiscLayers = static_cast<int>(std::size(Layout::discLayers));
  static constexpr int nLayers = nBarrelLayers + nDiscLayers;
  static constexpr bool hasBarrel = nBarrelLayers > 0;
  static constexpr bool hasDiscs = nDiscLayers > 0;
  /// Whether any layer's sensors are curved (SensorSpec::curved).
  static constexpr bool hasCurvedSensors = [] {
    const auto curved = [](const auto& layer) { return layer.sensor.curved(); };
    return std::ranges::any_of(Layout::barrelLayers, curved) || std::ranges::any_of(Layout::discLayers, curved);
  }();
  /// Passive scatterers (PassiveSpec), if the layout declares any.
  static constexpr int nPassives = [] {
    if constexpr (requires { Layout::passives; }) {
      return static_cast<int>(std::size(Layout::passives));
    } else {
      return 0;
    }
  }();

  /// Layers are numbered barrel first, then discs.
  [[nodiscard]] static constexpr bool isBarrelLayer(LayerID layer) { return hasBarrel && std::cmp_less(layer, nBarrelLayers); }

  /// Position of a disc layer within Layout::discLayers.
  [[nodiscard]] static constexpr std::size_t discIndex(LayerID layer) { return static_cast<std::size_t>(layer) - static_cast<std::size_t>(nBarrelLayers); }

  /// Passive scatterer \a index of the layout.
  [[nodiscard]] static constexpr const PassiveSpec& passive(std::size_t index)
    requires(nPassives > 0)
  {
    return Layout::passives[index];
  }

  /// A layer's stated sensor count.
  [[nodiscard]] static constexpr int sensorCount(LayerID layer) { return isBarrelLayer(layer) ? Layout::barrelLayers[layer].sensorCount : Layout::discLayers[discIndex(layer)].sensorCount; }
};

/// Total sensor count of a detector, for the detector struct to use when
/// declaring its own nSensors -- keeping that number derived rather than typed
/// in a second time.
template <typename BarrelArray, typename DiscArray>
[[nodiscard]] constexpr int countSensors(const BarrelArray& barrel, const DiscArray& discs)
{
  int n = 0;
  for (const auto& layer : barrel) {
    n += layer.nSensors();
  }
  for (const auto& layer : discs) {
    n += layer.nSensors();
  }
  return n;
}

} // namespace o2::gtrk

#endif // GTRK_BASE_LAYERSPEC_H_

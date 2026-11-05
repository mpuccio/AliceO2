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

/// \file TRK.h
/// \brief ALICE 3's tracker, TRK's barrel and FT3's discs, as ours: its
///        cluster format, its topology dictionaries' depths, its readout
///        frames and its field.
///
/// **Clusters.** `o2::trkft3::TRKCluster` and `FT3Cluster` come one stream per
/// layer, each with its own frames, and every cluster carries its pattern:
/// `trkLocalClusters` gives it its entry's position and widths in the layer's
/// topology dictionary, and
/// each stream goes through `toClusters` on its own (Clusters.h).
///
/// **Topology dictionaries.** O2's `o2::trkft3::TopologyDictionary` per
/// detector, built by its CheckTopologyDictionary.C, CCDB conditions
/// (`TRK/Calib/ClusterDictionary`, `FT3/Calib/ClusterDictionary`;
/// ConditionsTRK.h). A dictionary's layers are O2's: TRK's numbered over the
/// barrel (VD, then ML/OT), FT3's discs; ours are the same in that order, the
/// discs after the barrel.
///
/// **Readout frames.** `TRKAlmiraParam`, one length, bias and delay per layer,
/// as the digitiser took it from its configuration: a CCDB condition of the
/// device (`TRK/Config/AlmiraParam`, ConditionsTRK.h), as the ITS's
/// `ITS/Config/AlpideParam`.
///
/// **Field.** An ALICE 3 run records no field in its GRP (the L3 current is
/// 0): O2's transport takes it from `o2::field::ALICE3MagneticField`, which
/// compiles the macro `ALICE3_MAGFIELD_MACRO` names at run time. The
/// tracking reads the same object, at run time, and never a compiled value.

#ifndef GTRK_ALICE_TRK_H_
#define GTRK_ALICE_TRK_H_

#include "GenericTracking/Alice/Clusters.h"
#include "GenericTracking/Base/FrameTiming.h"

#include "DataFormatsTRKFT3/Cluster.h"
#include "DetectorsCommonDataFormats/DetID.h"

#include <cstdint>
#include <span>
#include <vector>

namespace o2::trkft3
{
class ROFRecord;
template <int DetID>
class TopologyDictionary;
} // namespace o2::trkft3

namespace o2::gtrk::alice
{

using TrkDictionary = o2::trkft3::TopologyDictionary<o2::detectors::DetID::TRK>;
using Ft3Dictionary = o2::trkft3::TopologyDictionary<o2::detectors::DetID::FT3>;

/// Every TRK cluster's local position, in order, reading its pattern as it
/// goes, in our sensor's frame -- O2's x and z, and on the VD the flat x,
/// which is our unrolled x, not curved --
/// with its spans and pixel count: its \a dictionary entry's position (the
/// centre of gravity, plus the entry's median offsets where the dictionary
/// applies them) and its widths across and along as its errors. Throws
/// std::runtime_error if the patterns run out, or a cluster names no TRK
/// sub-detector or no layer of the dictionary.
[[nodiscard]] std::vector<LocalCluster> trkLocalClusters(std::span<const o2::trkft3::TRKCluster> clusters, std::span<const unsigned char> patterns,
                                                         const TrkDictionary& dictionary);

/// The same for FT3's clusters, ML/OT's chip (sub-detector 2), each chip
/// numbered as our sensor: TRKFT3::kFirstDiscSensor plus FT3's chip ID.
[[nodiscard]] std::vector<LocalCluster> ft3LocalClusters(std::span<const o2::trkft3::FT3Cluster> clusters, std::span<const unsigned char> patterns,
                                                         const Ft3Dictionary& dictionary);

/// Clusters.h's `toClusters` for one layer's stream of TRK or FT3 clusters.
[[nodiscard]] std::vector<Cluster> toClusters(std::span<const LocalCluster> local, std::span<const o2::trkft3::ROFRecord> rofs, ReadoutFrames frames,
                                              std::uint32_t tfFirstOrbit);

/// Clusters.h's `readoutFrames` for one layer's stream of TRK or FT3 clusters.
[[nodiscard]] ReadoutFrames readoutFrames(std::span<const o2::trkft3::ROFRecord> rofs);

/// Our layers' measurement depths [cm] (Geometry::measurementDepth), barrel
/// then discs: each dictionary layer's response depth, O2's local y, in our
/// frame, whose normal is O2's -y on a flat sensor and the distance from the
/// axis minus the radius on a curved one (Placement.cxx, LayerSpec.h): the
/// opposite sign on every layer.
[[nodiscard]] std::vector<double> measurementDepthsCM(const TrkDictionary& trk, const Ft3Dictionary& ft3);

/// The TRKFT3 layers' readout frames from O2's `TRKAlmiraParam` as it stands
/// -- a device's from the CCDB (`TRK/Config/AlmiraParam`, ConditionsTRK.h), a
/// tool's from the digitiser's configuration --: each layer's length (O2's
/// default, 18 BC, for an unset one, as the digitiser takes it) and first
/// start, the bias and the delay together (`TRKDigitizerSpec`).
[[nodiscard]] std::vector<FrameTiming> trkFrameTiming();

/// Bz at the origin [kG] of O2's ALICE 3 field. Throws std::runtime_error if
/// `ALICE3_MAGFIELD_MACRO` is unset or names no file: O2 would then fall back
/// to a toy -10 kG field with only an error in its log.
[[nodiscard]] double alice3BzKiloGauss();

} // namespace o2::gtrk::alice

#endif // GTRK_ALICE_TRK_H_

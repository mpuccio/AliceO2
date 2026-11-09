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

/// \file TrackingInput.h
/// \brief What the sink's tracker reads for one timeframe, on disk.
///
/// The measurement columns, the grid offsets and the vertex list exactly as
/// the sink received them, so that `o2-gtrk-track-bench` runs the tracking
/// loop on the same input without DPL, the producer or the QA around it.
/// The configuration the sink ran with
/// goes beside it as `<file>.ini`, which `ConfigurableParam::updateFromFile`
/// reads back to the bit.
///
/// The file is raw and native-endian: a header with the blocks' counts and the
/// beam line, then the readout frames, every column in schema order, the
/// offsets, the vertices and the measurement depths, each block with its
/// element count and its `Digest`, which the reader checks.
///
/// With MC, the truth the comparison reads goes beside it as `<file>.truth`:
/// each row's particles and layer, and
/// the particles' kinematics, laid out the same way, so that the bench scores
/// its tracks as the sink does.

#ifndef GTRK_STUDY_TRACKINGINPUT_H_
#define GTRK_STUDY_TRACKINGINPUT_H_

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Workflow/BeamLine.h"
#include "GenericTracking/Study/TrackComparison.h"

#include <cstdint>
#include <format>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace o2::gtrk::study
{

/// One timeframe's tracking input, owned.
struct TrackingInput {
  MeasurementStore columns;
  std::vector<std::uint32_t> offsets; ///< the grid's bucket prefix offsets, then the row count
  std::vector<VertexSeed> vertices;   ///< time-sorted
  workflow::BeamLine beamLine;        ///< the one the rows were placed and binned about
  /// Each layer's readout frames, which the grid's time blocks and the slices'
  /// reach follow.
  std::vector<FrameTiming> frames;
  /// Each layer's measurement depth [cm] (Geometry::measurementDepth), which
  /// the rows were placed at.
  std::vector<double> depths;
};

/// Write \a view's \a rows rows, \a offsets, \a vertices, \a beamLine,
/// \a frames and \a depths to \a path; the caller writes the parameters to
/// `<path>.ini` with `writeParameters`. Throws std::runtime_error.
void writeTrackingInput(const std::string& path, const MeasurementView& view, std::size_t rows, std::span<const std::uint32_t> offsets, std::span<const VertexSeed> vertices, const workflow::BeamLine& beamLine, std::span<const FrameTiming> frames, std::span<const double> depths);

/// Every configurable parameter to \a path as INI, as
/// `ConfigurableParam::writeINI` writes it, but GTRKReco<DET>'s floating-point
/// values at full precision: ROOT prints them with `%g`, and a field of
/// -5.0066791 kG would read back as -5.00668.
void writeParameters(const std::string& path);

/// Read what `writeTrackingInput` wrote to \a path; the parameters are not
/// touched. Throws std::runtime_error on a short, foreign or corrupt file.
[[nodiscard]] TrackingInput readTrackingInput(const std::string& path);

/// \a input's readout frames, one per layer of \a Det. Throws
/// std::runtime_error for another count.
template <DetectorLayout Det>
[[nodiscard]] std::vector<FrameTiming> layerFrames(const TrackingInput& input)
{
  constexpr auto kLayers = static_cast<std::size_t>(Detector<Det>::nLayers);
  if (input.frames.size() != kLayers) {
    throw std::runtime_error(std::format("tracking input: readout frames for {} layers, {} has {}", input.frames.size(), Det::name, kLayers));
  }
  return input.frames;
}

/// One timeframe's MC truth, as `compareTrackSets` reads it.
struct TrackingTruth {
  std::vector<TruthRow> rows; ///< one per row of the tracking input, in its order
  ParticleTable particles;
};

/// Write \a rows -- one per row of the tracking input, whose cluster indices
/// are \a clusterIndex -- and \a particles to \a path. Throws std::runtime_error.
void writeTrackingTruth(const std::string& path, std::span<const std::uint32_t> clusterIndex, std::span<const TruthRow> rows, const ParticleTable& particles);

/// Read what `writeTrackingTruth` wrote to \a path for \a input. Throws
/// std::runtime_error on a short, foreign or corrupt file, or one whose rows
/// are not \a input's by count and cluster index.
[[nodiscard]] TrackingTruth readTrackingTruth(const std::string& path, const TrackingInput& input);

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_TRACKINGINPUT_H_

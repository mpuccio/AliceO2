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

/// \file testTrackingTruth.cxx
/// \brief The MC truth beside a tracking input (TrackingInput.h): it reads
///        back as written, and only for its own tracking input.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK tracking truth
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Study/TrackingInput.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace o2::gtrk::study
{

namespace
{

/// A file in the temporary directory, removed with the object.
struct TemporaryFile {
  std::string path = (std::filesystem::temp_directory_path() / ("gtrk-test-truth-" + std::to_string(getpid()))).string();
  TemporaryFile() = default;
  TemporaryFile(const TemporaryFile&) = delete;
  TemporaryFile& operator=(const TemporaryFile&) = delete;
  TemporaryFile(TemporaryFile&&) = delete;
  TemporaryFile& operator=(TemporaryFile&&) = delete;
  ~TemporaryFile() { std::filesystem::remove(path); }
};

/// A tracking input of \a clusterIndex's rows, only their cluster indices set.
TrackingInput inputOf(const std::vector<std::uint32_t>& clusterIndex)
{
  TrackingInput input;
  resizeMeasurementColumns(input.columns, clusterIndex.size());
  std::ranges::copy(clusterIndex, input.columns.clusterIndex.begin());
  return input;
}

const std::vector<std::uint32_t> kClusters{7U, 3U, 11U};

std::vector<TruthRow> someRows()
{
  return {{{.particles = {5U, 9U}, .count = 2, .firstLabelCounts = true}, 0},
          {{.particles = {9U}, .count = 1, .firstLabelCounts = false}, 1},
          {{.count = 0, .firstLabelCounts = false}, 6}}; // noise
}

const ParticleTable kParticles{{5U, {.pt = 0.25F, .primary = true, .eta = -0.5F, .phi = 1.5F, .origin = ParticleOrigin::kPhysicalPrimary, .pdg = 211}},
                               {9U,
                                {.pt = 3.F,
                                 .primary = false,
                                 .eta = 0.75F,
                                 .phi = -2.F,
                                 .origin = ParticleOrigin::kDecay,
                                 .pdg = -2212,
                                 .motherPdg = -3122,
                                 .process = 4,
                                 .vx = 1.5F,
                                 .vy = -2.25F,
                                 .vz = 7.F}}};

} // namespace

BOOST_AUTO_TEST_CASE(ItReadsBackAsWritten)
{
  const TemporaryFile file;
  const auto rows = someRows();
  writeTrackingTruth(file.path, kClusters, rows, kParticles);
  const auto truth = readTrackingTruth(file.path, inputOf(kClusters));

  BOOST_REQUIRE_EQUAL(truth.rows.size(), rows.size());
  for (std::size_t row = 0; row < rows.size(); ++row) {
    BOOST_TEST(truth.rows[row].particles == rows[row].particles);
    BOOST_TEST(truth.rows[row].count == rows[row].count);
    BOOST_TEST(truth.rows[row].layer == rows[row].layer);
    BOOST_TEST(truth.rows[row].firstLabelCounts == rows[row].firstLabelCounts);
  }
  BOOST_REQUIRE_EQUAL(truth.particles.size(), kParticles.size());
  for (const auto& [id, particle] : kParticles) {
    const auto& read = truth.particles.at(id);
    BOOST_TEST(read.pt == particle.pt);
    BOOST_TEST(read.primary == particle.primary);
    BOOST_TEST(read.eta == particle.eta);
    BOOST_TEST(read.phi == particle.phi);
    BOOST_TEST((read.origin == particle.origin));
    BOOST_TEST(read.pdg == particle.pdg);
    BOOST_TEST(read.motherPdg == particle.motherPdg);
    BOOST_TEST(read.process == particle.process);
    BOOST_TEST(read.vx == particle.vx);
    BOOST_TEST(read.vy == particle.vy);
    BOOST_TEST(read.vz == particle.vz);
  }
}

/// A tracking input keeps each layer's readout frames and the boxes its
/// producer made -- widened by the margin it applied -- whatever margin the
/// reader is configured with: the slices' reach reads the boxes against the
/// frames; and the depths its
/// rows were placed at.
BOOST_AUTO_TEST_CASE(TheInputKeepsEachLayersFramesAndItsBoxes)
{
  const TemporaryFile file;
  auto input = inputOf(kClusters);
  const std::vector<FrameTiming> frames{{.lengthBC = 22, .firstStartBC = 0}, {.lengthBC = 22, .firstStartBC = 0}, {.lengthBC = 81, .firstStartBC = 3}};
  const std::array<TimeEstPS, 3> boxes{frameBox(154, 176, 0), frameBox(154, 176, 30), frameBox(165, 246, 0)};
  for (std::size_t row = 0; row < boxes.size(); ++row) {
    input.columns.timePS[row] = boxes[row].start;
    input.columns.widthPS[row] = boxes[row].width;
  }
  const std::vector<std::uint32_t> offsets{0U, 3U};
  const std::vector<double> depths{8.25e-4, -41.25e-4, -40.73e-4};
  writeTrackingInput(file.path, measurementView(input.columns), kClusters.size(), offsets, {}, workflow::nominalBeamLine(), frames, depths);
  const auto read = readTrackingInput(file.path);
  BOOST_TEST(read.frames == frames);
  BOOST_TEST(read.depths == depths);
  for (std::size_t row = 0; row < boxes.size(); ++row) {
    BOOST_TEST(read.columns.timePS[row] == boxes[row].start);
    BOOST_TEST(read.columns.widthPS[row] == boxes[row].width);
  }
}

/// Frames for another count than the layout's layers are refused.
BOOST_AUTO_TEST_CASE(ItReadsOnlyForItsOwnInput)
{
  const TemporaryFile file;
  writeTrackingTruth(file.path, kClusters, someRows(), kParticles);
  BOOST_CHECK_THROW((void)readTrackingTruth(file.path, inputOf({7U, 3U})), std::runtime_error);
  BOOST_CHECK_THROW((void)readTrackingTruth(file.path, inputOf({7U, 11U, 3U})), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ItRefusesAShortOrForeignFile)
{
  const TemporaryFile file;
  BOOST_CHECK_THROW(writeTrackingTruth(file.path, std::span{kClusters}.first(2), someRows(), kParticles), std::runtime_error);

  writeTrackingTruth(file.path, kClusters, someRows(), kParticles);
  const auto size = std::filesystem::file_size(file.path);
  std::filesystem::resize_file(file.path, size - 1U);
  BOOST_CHECK_THROW((void)readTrackingTruth(file.path, inputOf(kClusters)), std::runtime_error);

  {
    std::ofstream out{file.path, std::ios::binary | std::ios::trunc};
    out << "GTRKTIN1 not a truth";
  }
  BOOST_CHECK_THROW((void)readTrackingTruth(file.path, inputOf(kClusters)), std::runtime_error);
}

} // namespace o2::gtrk::study

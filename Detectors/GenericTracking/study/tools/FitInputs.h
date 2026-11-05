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

/// \file FitInputs.h
/// \brief The triplet and neighbour fits' inputs, as o2-gtrk-track-bench
///        writes them for o2-gtrk-fit-bench.
///
/// One record per accepted triplet and per accepted neighbour of a run: the
/// measurements the fit reads and the scattering variances already taken from
/// the iteration's table, so that the fit needs no detector. A triplet carries
/// the curvature and its variance the stage stored, which the fit bench checks
/// bit for bit, and each record the chi2 bound it passed. The file is raw and
/// native-endian, like the tracking input, and rests on the records' layouts
/// (`PairPoint`'s among them): the magic changes with any of them, and a new
/// run of the bench is the only migration, as for the dumps.

#ifndef GTRK_STUDY_FITINPUTS_H_
#define GTRK_STUDY_FITINPUTS_H_

#include "GenericTracking/Tracking/PairCut.h"

#include <array>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace o2::gtrk::study
{

struct TripletFitInput {
  PairPoint first;
  PairPoint middle;
  PairPoint third;
  Scalar variance = 0.F; ///< the space-angle scattering variance at the middle (`tripletScattering`)
  Scalar chi2Max = 0.F;  ///< the stage's bound on the fit's chi2
  Scalar curvature = 0.F;
  Scalar varCurvature = 0.F;
  std::uint32_t iteration = 0U;
};

struct NeighbourFitInput {
  PairPoint a;
  PairPoint b;
  PairPoint c;
  PairPoint d;
  Scalar varianceSecond = 0.F; ///< at `b`, from `a -> c`
  Scalar varianceThird = 0.F;  ///< at `c`, from `b -> d`
  Scalar covariance = 0.F;     ///< of the two, from a passive scatterer between `b` and `c`
  Scalar chi2Max = 0.F;        ///< the stage's bound on the fit's chi2
  std::uint32_t iteration = 0U;
};

struct FitInputs {
  std::vector<TripletFitInput> triplets;
  std::vector<NeighbourFitInput> neighbours;
};

namespace fit_inputs_detail
{
constexpr std::array<char, 8> kFitInputsMagic{'G', 'T', 'R', 'K', 'F', 'I', 'T', '2'};
static_assert(std::is_trivially_copyable_v<TripletFitInput> && std::is_trivially_copyable_v<NeighbourFitInput>);

/// A stream's view of an object's bytes.
[[nodiscard]] inline const char* bytesOf(const void* data) { return static_cast<const char*>(data); }
[[nodiscard]] inline char* bytesOf(void* data) { return static_cast<char*>(data); }

template <typename T>
void writeVector(std::ofstream& out, const std::vector<T>& values)
{
  const auto count = static_cast<std::uint64_t>(values.size());
  out.write(bytesOf(&count), sizeof(count));
  out.write(bytesOf(values.data()), static_cast<std::streamsize>(values.size() * sizeof(T)));
}

template <typename T>
void readVector(std::ifstream& in, std::vector<T>& values)
{
  std::uint64_t count = 0U;
  in.read(bytesOf(&count), sizeof(count));
  values.resize(count);
  in.read(bytesOf(values.data()), static_cast<std::streamsize>(count * sizeof(T)));
}
} // namespace fit_inputs_detail

inline void writeFitInputs(const std::string& path, const FitInputs& inputs)
{
  std::ofstream out{path, std::ios::binary | std::ios::trunc};
  out.write(fit_inputs_detail::kFitInputsMagic.data(), fit_inputs_detail::kFitInputsMagic.size());
  fit_inputs_detail::writeVector(out, inputs.triplets);
  fit_inputs_detail::writeVector(out, inputs.neighbours);
  if (!out) {
    throw std::runtime_error("fit inputs: cannot write " + path);
  }
}

[[nodiscard]] inline FitInputs readFitInputs(const std::string& path)
{
  std::ifstream in{path, std::ios::binary};
  std::array<char, 8> magic{};
  in.read(magic.data(), magic.size());
  if (!in || magic != fit_inputs_detail::kFitInputsMagic) {
    throw std::runtime_error("fit inputs: " + path + " is not a fit-input file");
  }
  FitInputs inputs;
  fit_inputs_detail::readVector(in, inputs.triplets);
  fit_inputs_detail::readVector(in, inputs.neighbours);
  if (!in) {
    throw std::runtime_error("fit inputs: " + path + " is short");
  }
  return inputs;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_FITINPUTS_H_

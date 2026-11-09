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

/// \file fitBench.cxx
/// \brief o2-gtrk-fit-bench: the triplet and the adjacent-triplet fit alone.
///
/// Reads what `o2-gtrk-track-bench
/// --fit-inputs` wrote -- every accepted triplet and neighbour of a run -- and
/// runs the two fits over all of them, again and again, on one thread: the
/// time per call, and for the triplet fit the part in the kinks and their
/// gradient. It checks that every triplet's curvature and
/// variance are the stage's, bit for bit, and that every record passes its
/// chi2 bound, and prints a digest of every result: a faster fit must leave
/// the digests alone.

#include "FitInputs.h"
#include <boost/program_options.hpp>

#include "GenericTracking/Study/Digest.h"
#include "GenericTracking/Tracking/TripletFit.h"
#include "GenericTracking/Study/Validation.h"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <format>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace o2::gtrk::study
{
namespace
{

namespace bpo = boost::program_options;
using study::NeighbourFitInput;
using study::TripletFitInput;

/// Where the timed loops' results go, so that none of their work is left out.
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
volatile Scalar gSink = 0.F;

/// The shortest of \a repeat runs of \a body, in seconds; \a body returns a
/// value that depends on all its work.
template <typename Body>
double shortest(unsigned repeat, const Body& body)
{
  double best = std::numeric_limits<double>::infinity();
  Scalar sink = 0.F;
  for (unsigned run = 0; run < repeat; ++run) {
    const auto start = std::chrono::steady_clock::now();
    sink += body();
    best = std::min(best, std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
  }
  gSink = sink;
  return best;
}

[[nodiscard]] bool fitTripletInput(const TripletFitInput& input, TripletFitResult& result)
{
  TripletFitFactor factor;
  return makeTripletFitFactor(input.first, input.middle, input.third, input.variance, factor) && fitTriplet(factor, result);
}

[[nodiscard]] bool fitNeighbourInput(const NeighbourFitInput& input, TripletFitResult& result)
{
  KinkTerms first;
  KinkTerms second;
  return makeKinkTerms(input.a, input.b, input.c, first) && makeKinkTerms(input.b, input.c, input.d, second) && fitAdjacentTriplets(first, second, input.a, input.b, input.c, input.d, input.varianceSecond, input.varianceThird, result, input.covariance);
}

void addResult(Digest& digest, bool fitted, const TripletFitResult& result)
{
  digest.add(static_cast<std::uint32_t>(fitted));
  digest.add(result.curvature);
  digest.add(result.varCurvature);
  digest.add(result.chi2);
}

double perCall(double seconds, std::size_t calls) { return calls == 0U ? 0. : seconds * 1.e9 / static_cast<double>(calls); }

void benchTriplets(const std::vector<TripletFitInput>& inputs, unsigned repeat)
{
  Digest digest;
  std::size_t same = 0U;
  std::size_t passed = 0U;
  for (const auto& input : inputs) {
    TripletFitResult result;
    const bool fitted = fitTripletInput(input, result);
    addResult(digest, fitted, result);
    same += static_cast<std::size_t>(fitted && std::bit_cast<std::uint32_t>(result.curvature) == std::bit_cast<std::uint32_t>(input.curvature) && std::bit_cast<std::uint32_t>(result.varCurvature) == std::bit_cast<std::uint32_t>(input.varCurvature));
    passed += static_cast<std::size_t>(fitted && result.chi2 <= input.chi2Max);
  }

  const double fit = shortest(repeat, [&] {
    Scalar sum = 0.F;
    for (const auto& input : inputs) {
      TripletFitResult result;
      if (fitTripletInput(input, result)) {
        sum += result.curvature;
      }
    }
    return sum;
  });
  const double gradient = shortest(repeat, [&] {
    Scalar sum = 0.F;
    for (const auto& input : inputs) {
      sum += kinkGradient(input.first, input.middle, input.third).theta[0];
    }
    return sum;
  });

  const double fitNS = perCall(fit, inputs.size());
  const double gradientNS = perCall(gradient, inputs.size());
  std::cout << std::format("triplets: {} records; fit {:.1f} ns/call: the kinks and their gradient {:.1f} ns ({:.0f} %), the rest {:.1f} ns\n", inputs.size(), fitNS, gradientNS, fitNS > 0. ? 100. * gradientNS / fitNS : 0., fitNS - gradientNS);
  std::cout << std::format("  digest {:016x}; the stage's curvature and variance {}/{}, within the chi2 bound {}/{}\n", digest.result(), same, inputs.size(), passed, inputs.size());
}

void benchNeighbours(const std::vector<NeighbourFitInput>& inputs, unsigned repeat)
{
  Digest digest;
  std::size_t passed = 0U;
  for (const auto& input : inputs) {
    TripletFitResult result;
    const bool fitted = fitNeighbourInput(input, result);
    addResult(digest, fitted, result);
    passed += static_cast<std::size_t>(fitted && result.chi2 <= input.chi2Max);
  }

  const double fit = shortest(repeat, [&] {
    Scalar sum = 0.F;
    for (const auto& input : inputs) {
      TripletFitResult result;
      if (fitNeighbourInput(input, result)) {
        sum += result.chi2;
      }
    }
    return sum;
  });

  // From the two triplets' kink terms, computed here beforehand: the fit's
  // time alone, without making them (the stage makes them from the points).
  // A record whose terms cannot be made fails here as in the fit from the
  // points.
  std::vector<std::array<KinkTerms, 2>> terms(inputs.size());
  std::vector<std::uint8_t> made(inputs.size(), 0U);
  for (std::size_t record = 0; record < inputs.size(); ++record) {
    const auto& input = inputs[record];
    made[record] = makeKinkTerms(input.a, input.b, input.c, terms[record][0]) && makeKinkTerms(input.b, input.c, input.d, terms[record][1]) ? 1U : 0U;
  }
  const auto fromTerms = [&](std::size_t record, TripletFitResult& result) {
    const auto& input = inputs[record];
    return made[record] != 0U && fitAdjacentTriplets(terms[record][0], terms[record][1], input.a, input.b, input.c, input.d, input.varianceSecond, input.varianceThird, result, input.covariance);
  };
  Digest storedDigest;
  for (std::size_t record = 0; record < inputs.size(); ++record) {
    TripletFitResult result;
    const bool fitted = fromTerms(record, result);
    addResult(storedDigest, fitted, result);
  }
  const double stored = shortest(repeat, [&] {
    Scalar sum = 0.F;
    for (std::size_t record = 0; record < inputs.size(); ++record) {
      TripletFitResult result;
      if (fromTerms(record, result)) {
        sum += result.chi2;
      }
    }
    return sum;
  });
  std::cout << std::format("neighbours: {} records; fit {:.1f} ns/call, from the triplets' stored terms {:.1f} ns/call\n", inputs.size(), perCall(fit, inputs.size()), perCall(stored, inputs.size()));
  std::cout << std::format("  digest {:016x}, from the stored terms {:016x}; within the chi2 bound {}/{}\n", digest.result(), storedDigest.result(), passed, inputs.size());
}

int run(int argc, char** argv)
{
  std::string input;
  unsigned repeat = 5U;
  std::string which;
  bpo::options_description description("Run the triplet and the adjacent-triplet fit alone");
  description.add_options()                                                                                     //
    ("help,h", "this message")                                                                                  //
    ("input,i", bpo::value(&input)->default_value("fit-inputs.bin"), "o2-gtrk-track-bench --fit-inputs's file") //
    ("repeat,n", bpo::value(&repeat)->default_value(5U), "runs of each loop; the shortest counts")              //
    ("fits", bpo::value(&which)->default_value("all"), "triplets, neighbours or all");

  bpo::variables_map vm;
  try {
    bpo::store(bpo::parse_command_line(argc, argv, description), vm);
    bpo::notify(vm);
  } catch (const bpo::error& error) {
    std::cerr << error.what() << "\n\n"
              << description << '\n';
    return 1;
  }
  if (vm.contains("help")) {
    std::cout << description << '\n';
    return 0;
  }
  if (which != "all" && which != "triplets" && which != "neighbours") {
    std::cerr << "--fits: triplets, neighbours or all, not " << which << "\n\n"
              << description << '\n';
    return 1;
  }

  const auto inputs = study::readFitInputs(input);
  if (which == "all" || which == "triplets") {
    benchTriplets(inputs.triplets, std::max(repeat, 1U));
  }
  if (which == "all" || which == "neighbours") {
    benchNeighbours(inputs.neighbours, std::max(repeat, 1U));
  }
  return 0;
}

} // namespace
} // namespace o2::gtrk::study

int main(int argc, char** argv)
{
  try {
    return o2::gtrk::study::run(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}

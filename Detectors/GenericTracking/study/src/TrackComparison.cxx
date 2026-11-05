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

#include "GenericTracking/Study/TrackComparison.h"
#include "GenericTracking/Study/Statistics.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <format>
#include <functional>
#include <iterator>
#include <map>
#include <numbers>
#include <unordered_map>
#include <utility>
#include <vector>

namespace o2::gtrk::study
{

namespace
{

struct Particle {
  float pt = 0.F;
  bool primary = false;
  float eta = 0.F;
  float phi = 0.F;
  std::uint64_t layers = 0U;      ///< a bit per layer with a measurement of it
  std::uint64_t firstLayers = 0U; ///< the same, from the first labels O2's QA counts
  ParticleLayers chain;           ///< its layers and longest chain (particleLayers)
  ParticleOrigin origin = ParticleOrigin::kMaterial;
  std::int32_t species = 0; ///< |PDG code|
  DecayFamily family = DecayFamily::kOther;
  std::size_t production = 0U; ///< its production bin (productionBin)
};

/// Where one tracker's fakes go wrong: every measurement of a fake that does
/// not carry its label, by layer, by whether the label's particle had a
/// measurement of its own on that layer, and by whose it is instead.
struct FakeAnatomy {
  static constexpr std::size_t kWrongSlots = 5U; ///< 1, 2, 3, 4 or more wrong; slot 0: no label at all
  static constexpr std::size_t kLayerSlots = 64U;
  std::array<std::uint64_t, kWrongSlots> byWrong{};
  std::array<std::uint64_t, kLayerSlots> wrongByLayer{};
  std::uint64_t wrong = 0U;
  std::uint64_t swaps = 0U; ///< the label had its own measurement on the layer
  std::uint64_t holes = 0U; ///< it had none there
  std::uint64_t noise = 0U; ///< the wrong measurement carries no particle
  std::uint64_t primaryStranger = 0U;
  std::uint64_t secondaryStranger = 0U;
  std::uint64_t labelPrimary = 0U;
  std::uint64_t labelSecondary = 0U;
  std::vector<std::uint64_t> labels; ///< every fake's label, for the other tracker's view of it
};

/// One tracker's tracks against the particles.
struct Tally {
  std::unordered_map<std::uint64_t, std::uint32_t> good;       ///< good tracks per particle
  std::unordered_map<std::uint64_t, std::uint64_t> goodLayers; ///< per particle, the layers of its good track with the most
  std::unordered_map<std::uint64_t, bool> fakeLabel;           ///< particles some fake is labelled with
  std::array<std::uint64_t, kScoreLengthSlots> byLength{};
  std::array<std::uint64_t, kScoreLengthSlots> fakesByLength{};
  /// Tracks short of every layer, and their fakes, by their layers' bits.
  std::map<std::uint64_t, std::array<std::uint64_t, 2>> byPattern;
  std::uint64_t tracks = 0U;
  std::uint64_t fakes = 0U;
  FakeAnatomy anatomy;
};

/// Each track's particle by `majorityParticle`, and whether every
/// measurement carries it.
Tally tally(std::span<const TruthRow> rows, const TrackRows& tracks, const std::unordered_map<std::uint64_t, Particle>& particles, int nLayers)
{
  Tally result;
  auto& anatomy = result.anatomy;
  const auto primary = [&](std::uint64_t id) {
    const auto found = particles.find(id);
    return found != particles.end() && found->second.primary;
  };
  for (const auto& track : tracks) {
    const auto [label, most] = majorityParticle(rows, track);
    const auto length = std::min<std::size_t>(track.size(), kScoreLengthSlots - 1U);
    std::uint64_t pattern = 0U;
    for (const auto row : track) {
      pattern |= std::uint64_t{1} << rows[row].layer;
    }
    // Every layer of the detector: the comparison with O2 runs on the ITS
    // and MFT, whose class has no mask.
    const bool partial = std::cmp_less(std::popcount(pattern), nLayers);
    ++result.tracks;
    ++result.byLength[length];
    if (partial) {
      ++result.byPattern[pattern][0];
    }
    if (label != kNoTrackParticle && most == track.size()) {
      ++result.good[label];
      auto& layers = result.goodLayers[label];
      if (std::popcount(pattern) > std::popcount(layers)) {
        layers = pattern;
      }
      continue;
    }
    ++result.fakes;
    ++result.fakesByLength[length];
    if (partial) {
      ++result.byPattern[pattern][1];
    }
    if (label == kNoTrackParticle) {
      ++anatomy.byWrong[0];
      continue;
    }
    result.fakeLabel[label] = true;
    anatomy.labels.push_back(label);
    ++(primary(label) ? anatomy.labelPrimary : anatomy.labelSecondary);
    ++anatomy.byWrong[std::min<std::size_t>(track.size() - most, FakeAnatomy::kWrongSlots - 1U)];
    const auto own = particles.find(label);
    const std::uint64_t ownLayers = own == particles.end() ? 0U : own->second.layers;
    for (const auto row : track) {
      const auto& truth = rows[row];
      if (truth.carries(label)) {
        continue;
      }
      ++anatomy.wrong;
      ++anatomy.wrongByLayer[std::min<std::size_t>(truth.layer, FakeAnatomy::kLayerSlots - 1U)];
      ++((ownLayers & (std::uint64_t{1} << truth.layer)) != 0U ? anatomy.swaps : anatomy.holes);
      if (truth.count == 0U) {
        ++anatomy.noise;
      } else if (std::any_of(truth.particles.begin(), truth.particles.begin() + truth.count, primary)) {
        ++anatomy.primaryStranger;
      } else {
        ++anatomy.secondaryStranger;
      }
    }
  }
  return result;
}

/// \a mine's fakes, and what \a other made of their particles.
void describeFakes(std::vector<std::string>& report, const std::string& name, const Tally& mine, const Tally& other, int nLayers)
{
  const auto& a = mine.anatomy;
  report.push_back(
    std::format("  {:>6} fakes: {}; labelled with a primary {:.1f} %, a secondary {:.1f} %, nothing {:.1f} %; by wrong measurements 1 / 2 / 3 / "
                "4+: {:.1f} / {:.1f} / {:.1f} / {:.1f} %",
                name, mine.fakes, percent(a.labelPrimary, mine.fakes), percent(a.labelSecondary, mine.fakes), percent(a.byWrong[0], mine.fakes),
                percent(a.byWrong[1], mine.fakes), percent(a.byWrong[2], mine.fakes), percent(a.byWrong[3], mine.fakes), percent(a.byWrong[4], mine.fakes)));
  std::string layers;
  for (int layer = 0; layer < nLayers; ++layer) {
    layers += std::format(" {}: {:.1f}", layer, percent(a.wrongByLayer[static_cast<std::size_t>(layer)], a.wrong));
  }
  report.push_back(std::format("         {} wrong measurements, by layer (%):{}", a.wrong, layers));
  report.push_back(
    std::format("         the label's own measurement on that layer passed over {:.1f} %, none of its own there {:.1f} %; the wrong one "
                "noise {:.1f} %, another primary's {:.1f} %, a secondary's {:.1f} %",
                percent(a.swaps, a.wrong), percent(a.holes, a.wrong), percent(a.noise, a.wrong), percent(a.primaryStranger, a.wrong),
                percent(a.secondaryStranger, a.wrong)));
  std::uint64_t cleanElsewhere = 0U;
  for (const auto label : a.labels) {
    cleanElsewhere += other.good.contains(label) ? 1U : 0U;
  }
  report.push_back(std::format("         their particles found with a good track by the other tracker: {:.1f} %", percent(cleanElsewhere, a.labels.size())));
}

/// The tracks short of every layer by their layers, most common first: every
/// pattern with at least 0.5 % of them, the rest summed.
std::string describePatterns(const std::string& name, const Tally& tally)
{
  std::vector<std::pair<std::uint64_t, std::array<std::uint64_t, 2>>> patterns(tally.byPattern.begin(), tally.byPattern.end());
  std::ranges::stable_sort(patterns, std::greater{}, [](const auto& entry) { return entry.second[0]; });
  std::uint64_t total = 0U;
  for (const auto& [pattern, counts] : patterns) {
    total += counts[0];
  }
  std::string line = std::format("  {:>6}: {} tracks short of a layer, by layers (fake %):", name, total);
  std::array<std::uint64_t, 2> rest{};
  for (const auto& [pattern, counts] : patterns) {
    if (200U * counts[0] < total) {
      rest[0] += counts[0];
      rest[1] += counts[1];
      continue;
    }
    std::string layers;
    for (std::uint64_t bits = pattern; bits != 0U; bits &= bits - 1U) {
      const auto layer = std::countr_zero(bits);
      layers += layer < 10 ? std::format("{}", layer) : std::format("[{}]", layer);
    }
    line += std::format(" {}: {} ({:.1f})", layers, counts[0], percent(counts[1], counts[0]));
  }
  if (rest[0] != 0U) {
    line += std::format(" other: {} ({:.1f})", rest[0], percent(rest[1], rest[0]));
  }
  return line;
}

std::size_t ptBin(float pt)
{
  const auto bin = static_cast<std::size_t>(std::ranges::upper_bound(kScorePtEdges, static_cast<double>(pt)) - kScorePtEdges.begin());
  return bin == 0U ? 0U : std::min(bin - 1U, kScorePtBins - 1U);
}

using EtaEdges = std::array<double, kScoreEtaBins + 1U>;

std::size_t etaBin(float eta, const EtaEdges& edges)
{
  const auto bin = static_cast<std::size_t>(std::ranges::upper_bound(edges, static_cast<double>(eta)) - edges.begin());
  return std::clamp<std::size_t>(bin, 1U, kScoreEtaBins) - 1U;
}

/// In double, unlike `wrapPhi` (TrackingQA.h).
std::size_t phiBin(float phi)
{
  constexpr double kTwoPi = 2. * std::numbers::pi;
  const double wrapped = std::fmod(std::fmod(static_cast<double>(phi), kTwoPi) + kTwoPi, kTwoPi);
  return std::min(static_cast<std::size_t>(wrapped / kTwoPi * static_cast<double>(kScorePhiBins)), kScorePhiBins - 1U);
}

/// Efficiency, fake-only and clone rates of one tracker for one findable class.
FindableScore rates(const std::unordered_map<std::uint64_t, Particle>& particles, const Tally& tally, const EtaEdges& etaEdges, auto findable)
{
  FindableScore result;
  for (const auto& [id, particle] : particles) {
    if (!findable(particle)) {
      continue;
    }
    const auto bin = ptBin(particle.pt);
    const auto eta = etaBin(particle.eta, etaEdges);
    const auto phi = phiBin(particle.phi);
    ++result.findable;
    ++result.findableByPt[bin];
    ++result.findableByEta[eta];
    ++result.findableByPhi[phi];
    auto& species = result.bySpecies[particle.species];
    ++species[0];
    const auto good = tally.good.find(id);
    const std::uint32_t count = good == tally.good.end() ? 0U : good->second;
    result.withFake += tally.fakeLabel.contains(id) ? 1U : 0U;
    if (count != 0U) {
      ++result.found;
      ++result.foundByPt[bin];
      ++result.foundByEta[eta];
      ++result.foundByPhi[phi];
      ++species[1];
      result.cloned += count > 1U ? 1U : 0U;
      result.cloneTracks += count - 1U;
      const auto layers = tally.goodLayers.find(id);
      const std::uint64_t covered = layers == tally.goodLayers.end() ? 0U : (layers->second & particle.layers);
      result.layersFound += static_cast<std::uint64_t>(std::popcount(covered));
      result.layersOf += static_cast<std::uint64_t>(std::popcount(particle.layers));
      const double fraction = static_cast<double>(std::popcount(covered)) / static_cast<double>(std::popcount(particle.layers));
      result.layerFraction += fraction;
      result.layerFractionByPt[bin] += fraction;
      if (covered == particle.layers) {
        ++result.complete;
        ++result.completeByPt[bin];
      }
    } else if (tally.fakeLabel.contains(id)) {
      ++result.fakeOnly;
    }
  }
  return result;
}

std::string describe(const std::string& name, const FindableScore& r)
{
  std::string line = std::format(
    "  {:>6}: {:6.2f} % found ({:.2f} % complete, {:.1f} % of their layers), {:5.2f} % only with a fake, {:5.2f} % cloned; by pT, found (complete):", name,
    percent(r.found, r.findable), percent(r.complete, r.findable), percent(r.layersFound, r.layersOf), percent(r.fakeOnly, r.findable),
    percent(r.cloned, r.findable));
  for (std::size_t bin = 0; bin < kScorePtBins; ++bin) {
    line += std::format(" [{:g},{:g}) {:.1f} ({:.1f})", kScorePtEdges[bin], kScorePtEdges[bin + 1U], percent(r.foundByPt[bin], r.findableByPt[bin]),
                        percent(r.completeByPt[bin], r.findableByPt[bin]));
  }
  return line;
}

/// Found by eta and by phi, one line.
std::string describeAngles(const std::string& name, const FindableScore& r, const EtaEdges& etaEdges)
{
  std::string line = std::format("  {:>6} by eta:", name);
  for (std::size_t bin = 0; bin < kScoreEtaBins; ++bin) {
    line += std::format(" [{:.3g},{:.3g}) {:.1f}", etaEdges[bin], etaEdges[bin + 1U], percent(r.foundByEta[bin], r.findableByEta[bin]));
  }
  line += "; by phi, octants from 0:";
  for (std::size_t bin = 0; bin < kScorePhiBins; ++bin) {
    line += std::format(" {:.1f}", percent(r.foundByPhi[bin], r.findableByPhi[bin]));
  }
  return line;
}

/// Every particle the rows carry, with the layers it has measurements on
/// and its longest chain.
std::unordered_map<std::uint64_t, Particle> particlesOf(std::span<const TruthRow> rows, const ParticleTable& truth, const ComparisonSettings& settings)
{
  std::unordered_map<std::uint64_t, Particle> particles;
  for (const auto& row : rows) {
    for (std::uint8_t i = 0; i < row.count; ++i) {
      const auto found = truth.find(row.particles[i]);
      if (found == truth.end()) {
        continue;
      }
      auto& particle = particles[row.particles[i]];
      particle.pt = found->second.pt;
      particle.primary = found->second.primary;
      particle.eta = found->second.eta;
      particle.phi = found->second.phi;
      particle.origin = found->second.origin;
      particle.species = std::abs(found->second.pdg);
      particle.family = decayFamily(found->second.motherPdg);
      particle.production = productionBin(found->second, settings);
      particle.layers |= std::uint64_t{1} << row.layer;
      if (i == 0U && row.firstLabelCounts) {
        particle.firstLayers |= std::uint64_t{1} << row.layer;
      }
    }
  }
  for (const auto& [id, layers] : particleLayers(rows, settings)) {
    const auto particle = particles.find(id);
    if (particle != particles.end()) {
      particle->second.chain = layers;
    }
  }
  return particles;
}

/// The findable classes.
bool onAllLayers(const Particle& p, const ComparisonSettings& settings)
{
  return p.primary && static_cast<double>(p.pt) >= settings.minPtGeV && onEveryLayer(p.layers, settings);
}

bool onMinLayers(const Particle& p, const ComparisonSettings& settings)
{
  return p.primary && static_cast<double>(p.pt) >= settings.minPtGeV && std::popcount(p.layers) >= settings.minLayers;
}

/// O2's own QA definition (CheckTracksCA.C): first labels on every layer, a
/// primary, any pT; found with a good track, its further good tracks clones.
bool o2Findable(const Particle& p, const ComparisonSettings& settings)
{
  const std::uint64_t everyLayer = settings.nLayers >= 64 ? ~std::uint64_t{0} : (std::uint64_t{1} << settings.nLayers) - 1U;
  return p.primary && p.firstLayers == everyLayer;
}

/// The decay daughters by family and production bin.
std::array<std::vector<DecayCounts>, kDecayFamilies> decayCounts(const std::unordered_map<std::uint64_t, Particle>& particles, const Tally& tally,
                                                                 const ComparisonSettings& settings)
{
  std::array<std::vector<DecayCounts>, kDecayFamilies> result;
  for (auto& bins : result) {
    bins.resize(settings.productionEdges.size() + 1U);
  }
  for (const auto& [id, particle] : particles) {
    if (particle.origin != ParticleOrigin::kDecay || static_cast<double>(particle.pt) < settings.minPtGeV ||
        std::popcount(particle.layers) < settings.minLayers) {
      continue;
    }
    auto& counts = result[static_cast<std::size_t>(particle.family)][particle.production];
    const bool found = tally.good.contains(id);
    ++counts.onMinLayers;
    counts.foundOnMinLayers += found ? 1U : 0U;
    if (isFindable(particle.pt, particle.chain, settings)) {
      ++counts.findable;
      counts.found += found ? 1U : 0U;
    }
  }
  return result;
}

/// Whether \a p is a findable physical primary with (\a onDisc) or without a
/// measurement on a disc: a layer that is not a barrel layer.
bool primaryInRegion(const Particle& p, bool onDisc, const ComparisonSettings& settings)
{
  const bool disc = (p.layers & ~settings.barrelLayers) != 0U;
  return p.origin == ParticleOrigin::kPhysicalPrimary && disc == onDisc && isFindable(p.pt, p.chain, settings);
}

/// \a score of each origin.
std::array<FindableScore, kParticleOrigins> byOrigin(auto score)
{
  std::array<FindableScore, kParticleOrigins> result;
  for (std::size_t origin = 0; origin < kParticleOrigins; ++origin) {
    result[origin] = score(static_cast<ParticleOrigin>(origin));
  }
  return result;
}

TrackScore scoreOf(const std::unordered_map<std::uint64_t, Particle>& particles, const Tally& tally, const ComparisonSettings& settings)
{
  return {.tracks = tally.tracks,
          .fakes = tally.fakes,
          .tracksByLength = tally.byLength,
          .fakesByLength = tally.fakesByLength,
          .allLayers = rates(particles, tally, settings.etaEdges, [&](const Particle& p) { return onAllLayers(p, settings); }),
          .minLayers = rates(particles, tally, settings.etaEdges, [&](const Particle& p) { return onMinLayers(p, settings); }),
          .o2 = rates(particles, tally, settings.etaEdges, [&](const Particle& p) { return o2Findable(p, settings); }),
          .originAllLayers = byOrigin([&](ParticleOrigin origin) {
            return rates(particles, tally, settings.etaEdges, [&](const Particle& p) {
              return p.origin == origin && static_cast<double>(p.pt) >= settings.minPtGeV && onEveryLayer(p.layers, settings);
            });
          }),
          .originMinLayers = byOrigin([&](ParticleOrigin origin) {
            return rates(particles, tally, settings.etaEdges, [&](const Particle& p) {
              return p.origin == origin && static_cast<double>(p.pt) >= settings.minPtGeV && std::popcount(p.layers) >= settings.minLayers;
            });
          }),
          .originFindable = byOrigin([&](ParticleOrigin origin) {
            return rates(particles, tally, settings.etaEdges, [&](const Particle& p) { return p.origin == origin && isFindable(p.pt, p.chain, settings); });
          }),
          .primaryRegion = {rates(particles, tally, settings.etaEdges, [&](const Particle& p) { return primaryInRegion(p, false, settings); }),
                            rates(particles, tally, settings.etaEdges, [&](const Particle& p) { return primaryInRegion(p, true, settings); })},
          .decays = decayCounts(particles, tally, settings)};
}

/// The lines of a score, each for one tracker or for the findable they share.
std::string describeTracks(const std::string& name, const TrackScore& score)
{
  std::string line = std::format("  {:>6}: {} tracks, {:.2f} % fake; by clusters (fake %):", name, score.tracks, percent(score.fakes, score.tracks));
  for (std::size_t length = 0; length < kScoreLengthSlots; ++length) {
    if (score.tracksByLength[length] != 0U) {
      line += std::format(" {}: {} ({:.1f})", length, score.tracksByLength[length], percent(score.fakesByLength[length], score.tracksByLength[length]));
    }
  }
  return line;
}

std::string allLayersHeader(const ComparisonSettings& settings, std::uint64_t findable)
{
  if (settings.allLayersMask != 0U) {
    return std::format("  findable, primaries above {} GeV on every layer of mask {}: {}", settings.minPtGeV, settings.allLayersMask, findable);
  }
  return std::format("  findable, primaries above {} GeV on all {} layers: {}", settings.minPtGeV, settings.nLayers, findable);
}

std::string minLayersHeader(const ComparisonSettings& settings, std::uint64_t findable)
{
  return std::format("  findable, primaries above {} GeV on {} layers or more: {}", settings.minPtGeV, settings.minLayers, findable);
}

std::string o2Header(const ComparisonSettings& settings, std::uint64_t findable)
{
  return std::format("  O2's QA definition (CheckTracksCA.C), first labels on all {} layers, any pT: {} findable", settings.nLayers, findable);
}

std::string describeO2(const std::string& name, const FindableScore& r)
{
  return std::format("  {:>6}: found {} ({:.2f} %), clone tracks {}, with a fake track {}", name, r.found, percent(r.found, r.findable), r.cloneTracks,
                     r.withFake);
}

constexpr std::array<const char*, kParticleOrigins> kOriginNames{"physical primaries", "decay daughters", "material"};

std::string originHeader(const ComparisonSettings& settings, const TrackScore& score)
{
  std::string line = std::format("  findable by origin above {} GeV, on all {} layers / on {} or more / on a chain of {}:", settings.minPtGeV, settings.nLayers,
                                 settings.minLayers, settings.minLayers);
  for (std::size_t origin = 0; origin < kParticleOrigins; ++origin) {
    line += std::format("{} {} {} / {} / {}", origin == 0U ? "" : ",", kOriginNames[origin], score.originAllLayers[origin].findable,
                        score.originMinLayers[origin].findable, score.originFindable[origin].findable);
  }
  return line;
}

std::string describeOrigins(const std::string& name, const ComparisonSettings& settings, const TrackScore& score)
{
  std::string line = std::format("  {:>6} by origin, found % on all layers / {}+ / a chain of {}:", name, settings.minLayers, settings.minLayers);
  for (std::size_t origin = 0; origin < kParticleOrigins; ++origin) {
    const auto& all = score.originAllLayers[origin];
    const auto& some = score.originMinLayers[origin];
    const auto& adjacent = score.originFindable[origin];
    line += std::format("{} {} {:.2f} / {:.2f} / {:.2f}", origin == 0U ? "" : ",", kOriginNames[origin], percent(all.found, all.findable),
                        percent(some.found, some.findable), percent(adjacent.found, adjacent.findable));
  }
  return line;
}

/// One origin's particles on `minLayers` or more by species, most findable
/// first: every species with at least 0.5 % of them, the rest summed.
std::string describeSpecies(const std::string& name, const ComparisonSettings& settings, const TrackScore& score, std::size_t origin)
{
  const auto& r = score.originMinLayers[origin];
  std::vector<std::pair<std::int32_t, std::array<std::uint64_t, 2>>> species(r.bySpecies.begin(), r.bySpecies.end());
  std::ranges::stable_sort(species, std::greater{}, [](const auto& entry) { return entry.second[0]; });
  std::string line = std::format("  {:>6} {} on {}+ layers by species, found % (findable):", name, kOriginNames[origin], settings.minLayers);
  std::array<std::uint64_t, 2> rest{};
  for (const auto& [pdg, counts] : species) {
    if (200U * counts[0] < r.findable) {
      rest[0] += counts[0];
      rest[1] += counts[1];
      continue;
    }
    line += std::format(" {} {:.1f} ({})", speciesName(pdg), percent(counts[1], counts[0]), counts[0]);
  }
  if (rest[0] != 0U) {
    line += std::format(" other {:.1f} ({})", percent(rest[1], rest[0]), rest[0]);
  }
  return line;
}

constexpr std::array<const char*, kDecayFamilies> kFamilyNames{"K0S", "Lambda", "Sigma+-", "Xi, Omega", "K0L", "pi", "K", "other"};

DecayCounts operator+(DecayCounts a, const DecayCounts& b)
{
  return {.onMinLayers = a.onMinLayers + b.onMinLayers,
          .foundOnMinLayers = a.foundOnMinLayers + b.foundOnMinLayers,
          .findable = a.findable + b.findable,
          .found = a.found + b.found};
}

std::string decaysHeader(const ComparisonSettings& settings)
{
  std::string edges = "0";
  for (const auto edge : settings.productionEdges) {
    edges += std::format(", {:.1f}", edge);
  }
  return std::format(
    "  decay daughters above {} GeV by mother and production {} [cm] (bins from {}): found % of the findable, on a chain of {} "
    "layers (their count / all on {}+ layers)",
    settings.minPtGeV, settings.productionAlongZ ? "|z|" : "radius", edges, settings.minLayers, settings.minLayers);
}

/// One line of decay daughters: the total, then bin by bin.
std::string describeDecays(const std::string& name, const std::string& family, std::span<const DecayCounts> bins)
{
  DecayCounts total;
  for (const auto& bin : bins) {
    total = total + bin;
  }
  std::string line = std::format("  {:>6} {:<16} {:6.2f} ({}/{}), of all {:.2f} %; by production:", name, family, percent(total.found, total.findable),
                                 total.findable, total.onMinLayers, percent(total.foundOnMinLayers, total.onMinLayers));
  for (const auto& bin : bins) {
    line += std::format(" {:.1f} ({}/{})", percent(bin.found, bin.findable), bin.findable, bin.onMinLayers);
  }
  return line;
}

/// The decay daughters' lines: the V0s and hyperons and the decays in
/// flight summed, then each family.
std::vector<std::string> describeDecays(const std::string& name, const TrackScore& score)
{
  const auto bins = score.decays[0].size();
  std::vector<DecayCounts> hyperons(bins);
  std::vector<DecayCounts> flight(bins);
  for (std::size_t family = 0; family < kDecayFamilies; ++family) {
    const auto kind = static_cast<DecayFamily>(family);
    if (kind == DecayFamily::kOther) {
      continue;
    }
    auto& group = inFlight(kind) ? flight : hyperons;
    for (std::size_t bin = 0; bin < bins; ++bin) {
      group[bin] = group[bin] + score.decays[family][bin];
    }
  }
  std::vector<std::string> lines{describeDecays(name, "V0s and hyperons", hyperons), describeDecays(name, "in flight", flight)};
  for (std::size_t family = 0; family < kDecayFamilies; ++family) {
    lines.push_back(describeDecays(name, kFamilyNames[family], score.decays[family]));
  }
  return lines;
}

} // namespace

std::string speciesName(std::int32_t species)
{
  switch (species) {
    case 11:
      return "e";
    case 13:
      return "mu";
    case 211:
      return "pi";
    case 321:
      return "K";
    case 2212:
      return "p";
    case 3112:
      return "Sigma-";
    case 3222:
      return "Sigma+";
    case 3312:
      return "Xi";
    case 3334:
      return "Omega";
    case 1000010020:
      return "d";
    case 1000010030:
      return "t";
    case 1000020030:
      return "He3";
    case 1000020040:
      return "alpha";
    default:
      return std::format("{}", species);
  }
}

std::size_t productionBin(const ParticleTruth& particle, const ComparisonSettings& settings)
{
  const auto production =
    settings.productionAlongZ ? std::abs(static_cast<double>(particle.vz)) : std::hypot(static_cast<double>(particle.vx), static_cast<double>(particle.vy));
  return static_cast<std::size_t>(std::ranges::upper_bound(settings.productionEdges, production) - settings.productionEdges.begin());
}

DecayFamily decayFamily(std::int32_t motherPdg)
{
  switch (std::abs(motherPdg)) {
    case 310:
      return DecayFamily::kK0S;
    case 3122:
      return DecayFamily::kLambda;
    case 3112:
    case 3222:
      return DecayFamily::kSigma;
    case 3312:
    case 3322:
    case 3334:
      return DecayFamily::kCascade;
    case 130:
      return DecayFamily::kK0L;
    case 211:
      return DecayFamily::kPion;
    case 321:
      return DecayFamily::kKaon;
    default:
      return DecayFamily::kOther;
  }
}

std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> rowsOfParticles(std::span<const TruthRow> rows)
{
  std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> rowsOf;
  for (std::uint32_t row = 0; row < rows.size(); ++row) {
    for (std::uint8_t i = 0; i < rows[row].count; ++i) {
      rowsOf[rows[row].particles[i]].push_back(row);
    }
  }
  return rowsOf;
}

std::unordered_map<std::uint64_t, ParticleLayers> particleLayers(std::span<const TruthRow> rows, const ComparisonSettings& settings)
{
  const auto rowsOf = rowsOfParticles(rows);
  std::unordered_map<std::uint64_t, ParticleLayers> result;
  result.reserve(rowsOf.size());
  std::array<std::uint32_t, std::numeric_limits<std::uint64_t>::digits> depth{};
  // A step joins one of the particle's measurements on each of its layers
  // without crossing another layer.
  const auto direct = [&](std::span<const std::uint32_t> own, const ChainLink& link) {
    for (const auto a : own) {
      for (const auto b : own) {
        const auto& from = rows[a];
        const auto& to = rows[b];
        if (from.layer == link.from && to.layer == link.to &&
            !segmentLeavesHole(link.holes, settings.barrelLayers, settings.layerSpans.data(), from.r, from.z, to.r, to.z)) {
          return true;
        }
      }
    }
    return false;
  };
  for (const auto& [id, own] : rowsOf) {
    ParticleLayers particle;
    depth.fill(0U);
    for (const auto row : own) {
      particle.layers |= std::uint64_t{1} << rows[row].layer;
      depth[rows[row].layer] = 1U;
    }
    // Links in topological order: a link's `from` is complete when it comes.
    for (const auto& link : settings.chainLinks) {
      if (depth[link.from] != 0U && depth[link.to] != 0U && depth[link.from] + 1U > depth[link.to] && direct(own, link)) {
        depth[link.to] = depth[link.from] + 1U;
      }
    }
    particle.chain = *std::ranges::max_element(depth);
    result.emplace(id, particle);
  }
  return result;
}

TrackScore scoreTracks(std::span<const TruthRow> rows, const ParticleTable& truth, const TrackRows& tracks, const ComparisonSettings& settings)
{
  const auto particles = particlesOf(rows, truth, settings);
  return scoreOf(particles, tally(rows, tracks, particles, settings.nLayers), settings);
}

std::vector<std::string> describeScore(const TrackScore& score, const ComparisonSettings& settings, const std::string& name)
{
  std::vector<std::string> lines{describeTracks(name, score),
                                 allLayersHeader(settings, score.allLayers.findable),
                                 describe(name, score.allLayers),
                                 describeAngles(name, score.allLayers, settings.etaEdges),
                                 minLayersHeader(settings, score.minLayers.findable),
                                 describe(name, score.minLayers),
                                 o2Header(settings, score.o2.findable),
                                 describeO2(name, score.o2),
                                 originHeader(settings, score),
                                 describeOrigins(name, settings, score),
                                 describe(name + " primaries", score.originFindable[0]),
                                 describeSpecies(name, settings, score, 0U),
                                 describeSpecies(name, settings, score, 1U),
                                 describeSpecies(name, settings, score, 2U),
                                 decaysHeader(settings)};
  std::ranges::move(describeDecays(name, score), std::back_inserter(lines));
  return lines;
}

std::vector<std::string> compareTrackSets(std::span<const TruthRow> rows, const ParticleTable& truth, const TrackRows& ours, const TrackRows& theirs,
                                          const ComparisonSettings& settings, const std::string& ourName, const std::string& theirName)
{
  const auto particles = particlesOf(rows, truth, settings);
  std::vector<std::string> report;
  const auto ourTally = tally(rows, ours, particles, settings.nLayers);
  const auto theirTally = tally(rows, theirs, particles, settings.nLayers);
  const auto ourScore = scoreOf(particles, ourTally, settings);
  const auto theirScore = scoreOf(particles, theirTally, settings);
  report.push_back(std::format("comparison with {}: one set of truth rules for both", theirName));
  report.push_back(describeTracks(ourName, ourScore));
  report.push_back(describeTracks(theirName, theirScore));
  report.push_back(describePatterns(ourName, ourTally));
  report.push_back(describePatterns(theirName, theirTally));
  report.push_back(allLayersHeader(settings, ourScore.allLayers.findable));
  report.push_back(describe(ourName, ourScore.allLayers));
  report.push_back(describe(theirName, theirScore.allLayers));
  report.push_back(describeAngles(ourName, ourScore.allLayers, settings.etaEdges));
  report.push_back(describeAngles(theirName, theirScore.allLayers, settings.etaEdges));
  report.push_back(minLayersHeader(settings, ourScore.minLayers.findable));
  report.push_back(describe(ourName, ourScore.minLayers));
  report.push_back(describe(theirName, theirScore.minLayers));
  report.push_back(o2Header(settings, ourScore.o2.findable));
  report.push_back(describeO2(ourName, ourScore.o2));
  report.push_back(describeO2(theirName, theirScore.o2));
  report.push_back(originHeader(settings, ourScore));
  report.push_back(describeOrigins(ourName, settings, ourScore));
  report.push_back(describeOrigins(theirName, settings, theirScore));
  report.push_back(describe(ourName + " primaries", ourScore.originFindable[0]));
  report.push_back(describe(theirName + " primaries", theirScore.originFindable[0]));
  for (std::size_t origin = 0; origin < kParticleOrigins; ++origin) {
    report.push_back(describeSpecies(ourName, settings, ourScore, origin));
    report.push_back(describeSpecies(theirName, settings, theirScore, origin));
  }
  report.push_back(decaysHeader(settings));
  const auto ourDecays = describeDecays(ourName, ourScore);
  const auto theirDecays = describeDecays(theirName, theirScore);
  for (std::size_t line = 0; line < ourDecays.size(); ++line) {
    report.push_back(ourDecays[line]);
    report.push_back(theirDecays[line]);
  }

  describeFakes(report, ourName, ourTally, theirTally, settings.nLayers);
  describeFakes(report, theirName, theirTally, ourTally, settings.nLayers);

  // Track by track, by measurement sets.
  std::unordered_map<std::uint32_t, std::vector<std::size_t>> ourTracksOfRow;
  std::vector<std::vector<std::uint32_t>> ourSorted(ours.size());
  for (std::size_t track = 0; track < ours.size(); ++track) {
    ourSorted[track] = ours[track];
    std::ranges::sort(ourSorted[track]);
    for (const auto row : ours[track]) {
      ourTracksOfRow[row].push_back(track);
    }
  }
  enum Relation : std::uint8_t { kIdentical,
                                 kTheirsInside,
                                 kOursInside,
                                 kOverlap,
                                 kDisjoint,
                                 kRelations };
  std::array<std::uint64_t, kRelations> relations{};
  std::vector<std::uint32_t> sorted;
  std::map<std::size_t, std::size_t> shared;
  for (const auto& track : theirs) {
    sorted = track;
    std::ranges::sort(sorted);
    shared.clear();
    for (const auto row : sorted) {
      const auto found = ourTracksOfRow.find(row);
      if (found != ourTracksOfRow.end()) {
        for (const auto our : found->second) {
          ++shared[our];
        }
      }
    }
    if (shared.empty()) {
      ++relations[kDisjoint];
      continue;
    }
    const auto best = std::ranges::max_element(shared, {}, [](const auto& entry) { return entry.second; });
    const auto& our = ourSorted[best->first];
    if (best->second == sorted.size() && our.size() == sorted.size()) {
      ++relations[kIdentical];
    } else if (best->second == sorted.size()) {
      ++relations[kTheirsInside];
    } else if (best->second == our.size()) {
      ++relations[kOursInside];
    } else {
      ++relations[kOverlap];
    }
  }
  report.push_back(
    std::format("  {}'s {} tracks against ours, by clusters: identical {:.2f} %, inside one of ours {:.2f} %, containing one of ours "
                "{:.2f} %, overlapping {:.2f} %, sharing none {:.2f} %",
                theirName, theirs.size(), percent(relations[kIdentical], theirs.size()), percent(relations[kTheirsInside], theirs.size()),
                percent(relations[kOursInside], theirs.size()), percent(relations[kOverlap], theirs.size()), percent(relations[kDisjoint], theirs.size())));

  std::array<std::uint64_t, 4> who{}; // both, ours only, theirs only, neither
  for (const auto& [id, particle] : particles) {
    if (!onAllLayers(particle, settings)) {
      continue;
    }
    const bool mine = ourTally.good.contains(id);
    const bool other = theirTally.good.contains(id);
    std::size_t slot = 3U;
    if (mine && other) {
      slot = 0U;
    } else if (mine) {
      slot = 1U;
    } else if (other) {
      slot = 2U;
    }
    ++who[slot];
  }
  report.push_back(
    std::format("  the all-layer findable: found by both {}, by {} only {}, by {} only {}, by neither {}", who[0], ourName, who[1], theirName, who[2], who[3]));
  return report;
}

} // namespace o2::gtrk::study

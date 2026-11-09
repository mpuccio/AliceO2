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

/// \file TrackComparison.h
/// \brief Two trackers' tracks on the same measurements, under one set of
///        truth rules.
///
/// Both track sets are lists of measurement rows, so the rules cannot depend
/// on either tracker:
///
/// - a particle is **findable** if it is a primary above the pT floor with a
///   measurement on every layer (the ITS standard: 7 of 7), and, as a second
///   class, on a chain of at least `minLayers` layers (`isFindable`);
/// - a track takes the particle most of its measurements carry among their
///   labels, and is **good** if every one of them carries it and a **fake**
///   otherwise -- one wrong measurement is enough. This follows O2's rule
///   (`its::Tracker::computeTracksMClabels`): a measurement of two particles,
///   which pile-up makes common, counts for both, and one carrying no
///   particle -- noise -- for none;
/// - a findable particle is **found** with a good track, **found with a fake**
///   if its only track is a fake labelled with it, and **cloned** if it has
///   more than one good track; it is found **complete** when its good track
///   with the most layers has one on every layer the particle has a
///   measurement on (a 4-cluster stub of a 7-layer particle is found, not
///   complete).
///
/// The tracks short of a measurement on every layer are counted by their
/// layers, which shows holes and where roads end. The all-layer class is also
/// split by eta and phi, and counted once more under O2's own QA definition
/// (`CheckTracksCA.C`): a primary whose clusters' first labels name it on
/// every layer, at any pT.
///
/// Then track by track, by their measurement sets: identical, one a strict
/// subset of the other, overlapping, or sharing nothing; and particle by
/// particle, which of the findable each tracker found.
///
/// The decay daughters are split once more:
/// by their mother's family and by where they were produced, in bins at the
/// layers' positions, and counted on `minLayers` layers or more and, the
/// findable among them, on a chain of `minLayers` layers (`isFindable`).

#ifndef GTRK_STUDY_TRACKCOMPARISON_H_
#define GTRK_STUDY_TRACKCOMPARISON_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/LayerEnvelopes.h"
#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/LinkHoles.h"
#include "GenericTracking/Workflow/TrackLabel.h"

#include <algorithm>
#include <cmath>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace o2::gtrk::study
{

/// What the comparison knows about one measurement: the particles its labels
/// carry (the label rule's, TrackLabel.h), its layer and where it is.
struct TruthRow : workflow::RowParticles {
  LayerID layer = 0;
  /// The measurement's global radius and z [cm], which a findable chain's
  /// steps are tested with (LinkHoles.h).
  Scalar r = 0.F;
  Scalar z = 0.F;
};

// The label rule is the workflow's, shared with the O2-format adaptors.
using workflow::goodParticle;
using workflow::kNoTrackParticle;
using workflow::majorityParticle;
using workflow::MajorityParticle;

/// Where a particle comes from.
enum class ParticleOrigin : std::uint8_t {
  kPhysicalPrimary, ///< ALICE's definition (O2's MCTrackNavigator::isPhysicalPrimary)
  kDecay,           ///< otherwise a decay's: in transport, or feed-down the generator decayed
  kMaterial,        ///< otherwise the detector's: hadronic interactions, conversions, delta rays, ...
  kCount
};
inline constexpr std::size_t kParticleOrigins = static_cast<std::size_t>(ParticleOrigin::kCount);

/// A particle's MC truth.
struct ParticleTruth {
  float pt = 0.F;       ///< [GeV]
  bool primary = false; ///< the generator's (MCTrack::isPrimary), which O2's QA takes
  float eta = 0.F;
  float phi = 0.F; ///< [rad]
  ParticleOrigin origin = ParticleOrigin::kMaterial;
  std::int32_t pdg = 0;
  std::int32_t motherPdg = 0; ///< 0 without a mother
  std::uint8_t process = 0;   ///< how it was produced: a TMCProcess
  float vx = 0.F;             ///< the production vertex [cm]
  float vy = 0.F;
  float vz = 0.F;
};

using ParticleTable = std::unordered_map<std::uint64_t, ParticleTruth>;

/// Tracks as the measurement rows they hold, in any order.
using TrackRows = std::vector<std::vector<std::uint32_t>>;

/// \a tracks (a tracker's tracks) as the measurement rows each holds, in the
/// tracks' and their clusters' order.
template <typename Tracks>
[[nodiscard]] TrackRows trackRows(const Tracks& tracks)
{
  TrackRows rows(tracks.size());
  const auto view = tracks.view();
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    rows[track].assign(std::begin(view.clusters[track].index), std::begin(view.clusters[track].index) + view.count[track]);
  }
  return rows;
}

/// The bins of the efficiencies: pT [GeV], eta, and phi in octants from 0.
inline constexpr std::array<double, 9> kScorePtEdges{0.05, 0.1, 0.2, 0.3, 0.5, 1., 2., 5., std::numeric_limits<double>::infinity()};
inline constexpr std::size_t kScorePtBins = kScorePtEdges.size() - 1U;
/// The eta bins of a barrel; a disc-only layout takes its own acceptance
/// (`comparisonSettings`), the same number of bins.
inline constexpr std::array<double, 7> kScoreEtaEdges{-1.5, -1., -0.5, 0., 0.5, 1., 1.5};
inline constexpr std::size_t kScoreEtaBins = kScoreEtaEdges.size() - 1U;
inline constexpr std::size_t kScorePhiBins = 8U;
/// Track lengths in measurements; the last slot holds the longer ones too.
inline constexpr std::size_t kScoreLengthSlots = 16U;

/// A link of a findable chain: its layers, and the layers a step on it can
/// cross (`linkHoleLayers`), a bit per layer.
struct ChainLink {
  LayerID from = 0;
  LayerID to = 0;
  std::uint64_t holes = 0U;
};

struct ComparisonSettings {
  int nLayers = 0;       ///< the detector's
  int minLayers = 4;     ///< the second findable class
  double minPtGeV = 0.1; ///< the findable floor
  /// The QA's second findable rule (TrackingQA.h, `outwardLayers`): the
  /// layers its MC hits reach on their way out, at the same floor.
  int minOutwardLayers = 7;
  /// The inner edges of the decay daughters' production bins after the first,
  /// which starts at 0 [cm]: the layers' radii for a barrel, their |z| for
  /// discs. Empty, one bin.
  std::vector<double> productionEdges{}; // NOLINT(readability-redundant-member-init): GCC's -Wmissing-field-initializers wants it
  bool productionAlongZ = false;         ///< production measured as |vz|, not as the radius
  /// The efficiency's eta bins; outside them a particle counts in the nearest.
  std::array<double, kScoreEtaBins + 1U> etaEdges = kScoreEtaEdges;
  /// The links a findable chain follows, in topological source order: the
  /// table's links that do not skip a layer, each with the layers a step on
  /// it can cross. A step counts when its segment crosses none of them, as a
  /// tracklet of an iteration without skips.
  std::vector<ChainLink> chainLinks{}; // NOLINT(readability-redundant-member-init): GCC's -Wmissing-field-initializers wants it
  /// The barrel layers, a bit per layer, and each layer's span, for that test.
  std::uint64_t barrelLayers = 0U;
  std::vector<LayerSpan> layerSpans{}; // NOLINT(readability-redundant-member-init): GCC's -Wmissing-field-initializers wants it
  /// The layers the "on every layer" classes ask a measurement on, a bit per
  /// layer; 0 for every layer of the detector (e.g. the TRK's eleven barrel
  /// layers).
  std::uint64_t allLayersMask = 0U;
};

/// Whether a particle measured on \a layers is on every layer of \a settings.
[[nodiscard]] inline bool onEveryLayer(std::uint64_t layers, const ComparisonSettings& settings)
{
  if (settings.allLayersMask != 0U) {
    return (layers & settings.allLayersMask) == settings.allLayersMask;
  }
  return std::popcount(layers) == settings.nLayers;
}

/// The settings for detector \a Det: the decay daughters' production bins at
/// its layers' positions -- for a placed layout of barrel layers or of discs;
/// one bin otherwise. A placed disc-only layout bins eta over its own
/// acceptance, from the outer edge of its first disc to the inner edge of its
/// last, seen from the origin, rounded outwards to 0.1 -- the MFT's
/// [-3.8, -2.0] rather than the barrel's [-1.5, 1.5], which would hold every
/// one of its tracks in one bin.
template <DetectorLayout Det>
[[nodiscard]] ComparisonSettings comparisonSettings(int minLayers = 4, double minPtGeV = 0.1)
{
  ComparisonSettings settings{.nLayers = Detector<Det>::nLayers, .minLayers = minLayers, .minPtGeV = minPtGeV, .productionAlongZ = !Detector<Det>::hasBarrel};
  for (const auto from : layerOrder<Det>) {
    for (auto at = outLinks<Det>.begin(from); at < outLinks<Det>.end(from); ++at) {
      const auto link = outLinks<Det>.links[at];
      if (!Det::layerLinks[link].skipsLayer) {
        const auto to = Det::layerLinks[link].to;
        settings.chainLinks.push_back({.from = from, .to = to, .holes = linkHoleLayers<Det>(from, to)});
      }
    }
  }
  settings.barrelLayers = barrelLayerBits<Det>();
  for (int layer = 0; layer < Detector<Det>::nLayers; ++layer) {
    settings.layerSpans.push_back(layerSpan<Det>(static_cast<LayerID>(layer)));
  }
  if constexpr (!(Detector<Det>::hasBarrel && Detector<Det>::hasDiscs)) {
    for (int layer = 0; layer < Detector<Det>::nLayers; ++layer) {
      settings.productionEdges.push_back(std::abs(static_cast<double>(PlacedEnvelopes<Det>::get(static_cast<LayerID>(layer)).position)));
    }
  }
  if constexpr (!Detector<Det>::hasBarrel && Detector<Det>::hasDiscs) {
    const auto eta = [](double r, double z) { return std::asinh(z / r); };
    const auto& first = PlacedEnvelopes<Det>::get(static_cast<LayerID>(0));
    const auto& last = PlacedEnvelopes<Det>::get(static_cast<LayerID>(Detector<Det>::nLayers - 1));
    const double wide = eta(static_cast<double>(first.rMax), static_cast<double>(first.position));
    const double narrow = eta(static_cast<double>(last.rMin), static_cast<double>(last.position));
    // Rounded outwards to 0.1, then split evenly.
    const double low = std::floor(std::min(wide, narrow) * 10.) / 10.;
    const double high = std::ceil(std::max(wide, narrow) * 10.) / 10.;
    for (std::size_t bin = 0; bin <= kScoreEtaBins; ++bin) {
      settings.etaEdges[bin] = low + (high - low) * static_cast<double>(bin) / static_cast<double>(kScoreEtaBins);
    }
  }
  return settings;
}

/// A decay daughter's mother, as the score splits them: the V0s and
/// hyperons, then the decays in flight, then the rest.
enum class DecayFamily : std::uint8_t {
  kK0S,
  kLambda,
  kSigma,   ///< Sigma+ and Sigma-
  kCascade, ///< Xi and Omega, charged or not
  kK0L,
  kPion, ///< in flight
  kKaon, ///< in flight
  kOther,
  kCount
};
inline constexpr std::size_t kDecayFamilies = static_cast<std::size_t>(DecayFamily::kCount);
[[nodiscard]] DecayFamily decayFamily(std::int32_t motherPdg);
[[nodiscard]] constexpr bool inFlight(DecayFamily family) { return family == DecayFamily::kPion || family == DecayFamily::kKaon; }

/// The particles a study takes: the physical primaries, the decay daughters,
/// or those of a K0S or a hyperon (a K0S, Lambda, Sigma or cascade mother).
enum class Population : std::uint8_t { kPrimaries,
                                       kDecays,
                                       kV0s };
inline constexpr std::array<std::string_view, 3> kPopulationNames{"primaries", "decay daughters", "V0 daughters"};

/// The population an option names: primary, decay or v0; none for another.
[[nodiscard]] constexpr std::optional<Population> populationOf(std::string_view name)
{
  if (name == "primary") {
    return Population::kPrimaries;
  }
  if (name == "decay") {
    return Population::kDecays;
  }
  if (name == "v0") {
    return Population::kV0s;
  }
  return std::nullopt;
}

/// \a particle belongs to \a population.
[[nodiscard]] inline bool inPopulation(const ParticleTruth& particle, Population population)
{
  switch (population) {
    case Population::kPrimaries:
      return particle.origin == ParticleOrigin::kPhysicalPrimary;
    case Population::kDecays:
      return particle.origin == ParticleOrigin::kDecay;
    case Population::kV0s: {
      const auto family = decayFamily(particle.motherPdg);
      return particle.origin == ParticleOrigin::kDecay && (family == DecayFamily::kK0S || family == DecayFamily::kLambda || family == DecayFamily::kSigma || family == DecayFamily::kCascade);
    }
  }
  return false;
}

/// A particle's layers: a bit per layer with a measurement carrying it, and
/// its longest chain, the most layers a path of `chainLinks` visits whose
/// every step joins two of its measurements without crossing another layer.
struct ParticleLayers {
  std::uint64_t layers = 0U;
  std::uint32_t chain = 0U;
};

/// Whether a particle of transverse momentum \a pt [GeV] whose longest chain
/// has \a chain layers is findable: at or above
/// the pT floor, on a chain of `minLayers`.
[[nodiscard]] constexpr bool isFindable(float pt, const ParticleLayers& layers, const ComparisonSettings& settings) { return static_cast<double>(pt) >= settings.minPtGeV && settings.minLayers > 0 && std::cmp_greater_equal(layers.chain, settings.minLayers); }

/// Each particle the rows carry, with its rows in ascending order: every row
/// among whose labels it is.
[[nodiscard]] std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> rowsOfParticles(std::span<const TruthRow> rows);

/// Each particle the rows carry, with its layers and longest chain.
[[nodiscard]] std::unordered_map<std::uint64_t, ParticleLayers> particleLayers(std::span<const TruthRow> rows, const ComparisonSettings& settings);

/// A species, |PDG code|, as the score prints it: e, mu, pi, K, p, ...
[[nodiscard]] std::string speciesName(std::int32_t species);

/// \a particle's production bin under \a settings: 0 before the first edge.
[[nodiscard]] std::size_t productionBin(const ParticleTruth& particle, const ComparisonSettings& settings);

/// One family's decay daughters above the pT floor in one production bin.
struct DecayCounts {
  std::uint64_t onMinLayers = 0U; ///< with measurements on `minLayers` layers or more
  std::uint64_t foundOnMinLayers = 0U;
  std::uint64_t findable = 0U; ///< of them, on a chain of `minLayers` layers (`isFindable`)
  std::uint64_t found = 0U;
};

/// One tracker's particles of one findable class: found with a good track,
/// found only with a fake, cloned.
struct FindableScore {
  std::uint64_t findable = 0U;
  std::uint64_t found = 0U;
  std::uint64_t fakeOnly = 0U;
  std::uint64_t cloned = 0U;
  std::array<std::uint64_t, kScorePtBins> findableByPt{};
  std::array<std::uint64_t, kScorePtBins> foundByPt{};
  std::array<std::uint64_t, kScoreEtaBins> findableByEta{};
  std::array<std::uint64_t, kScoreEtaBins> foundByEta{};
  std::array<std::uint64_t, kScorePhiBins> findableByPhi{};
  std::array<std::uint64_t, kScorePhiBins> foundByPhi{};
  std::uint64_t cloneTracks = 0U; ///< good tracks beyond each particle's first
  std::uint64_t withFake = 0U;    ///< particles some fake is labelled with, found or not
  /// Found complete: the good track with the most layers covers every layer
  /// the particle has a measurement on. The layers that track covers and the
  /// layers the found particles have, for the fraction attached.
  std::uint64_t complete = 0U;
  std::array<std::uint64_t, kScorePtBins> completeByPt{};
  std::uint64_t layersFound = 0U;
  std::uint64_t layersOf = 0U;
  /// The layer fraction, summed over the findable: each found particle adds
  /// the fraction of its layers its good track with the most layers covers, an unfound one 0, clones and fakes nothing;
  /// divided by `findable` it is the mean a stub counts in part towards.
  double layerFraction = 0.;
  std::array<double, kScorePtBins> layerFractionByPt{};
  /// By species, |PDG code|: findable and found.
  std::map<std::int32_t, std::array<std::uint64_t, 2>> bySpecies;
};

/// `TrackScore::primaryRegion`'s classes: on no disc, on a disc too.
inline constexpr std::size_t kPrimaryRegions = 2;

/// One tracker's tracks under the truth rules, as numbers:
/// what `compareTrackSets` reports of it.
struct TrackScore {
  std::uint64_t tracks = 0U;
  std::uint64_t fakes = 0U;
  std::array<std::uint64_t, kScoreLengthSlots> tracksByLength{};
  std::array<std::uint64_t, kScoreLengthSlots> fakesByLength{};
  FindableScore allLayers; ///< primaries above the pT floor on every layer
  FindableScore minLayers; ///< the same on `minLayers` layers or more
  FindableScore o2;        ///< O2's QA definition (`CheckTracksCA.C`)
  /// Every particle above the pT floor by its origin, primary or not, on
  /// every layer and on `minLayers` or more.
  std::array<FindableScore, kParticleOrigins> originAllLayers;
  std::array<FindableScore, kParticleOrigins> originMinLayers;
  /// The same on a chain of `minLayers` layers: the findable (`isFindable`).
  std::array<FindableScore, kParticleOrigins> originFindable;
  /// The physical primaries of `originFindable` by where they are measured:
  /// [0] on no disc, [1] on a disc too. The
  /// class is the layers the particle crosses, not an eta cut.
  std::array<FindableScore, kPrimaryRegions> primaryRegion;
  /// The decay daughters by their mother's family and production bin, one
  /// more bin than `productionEdges`.
  std::array<std::vector<DecayCounts>, kDecayFamilies> decays;
};

/// \a tracks scored against \a truth.
[[nodiscard]] TrackScore scoreTracks(std::span<const TruthRow> rows, const ParticleTable& truth, const TrackRows& tracks, const ComparisonSettings& settings);

/// \a score as the lines `compareTrackSets` prints for the tracker \a name,
/// in its order.
[[nodiscard]] std::vector<std::string> describeScore(const TrackScore& score, const ComparisonSettings& settings, const std::string& name);

/// The comparison's report, one line per entry.
[[nodiscard]] std::vector<std::string> compareTrackSets(std::span<const TruthRow> rows, const ParticleTable& truth, const TrackRows& ours, const TrackRows& theirs, const ComparisonSettings& settings, const std::string& ourName, const std::string& theirName);

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_TRACKCOMPARISON_H_

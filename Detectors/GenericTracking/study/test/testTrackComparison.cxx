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

/// \file testTrackComparison.cxx
/// \brief The truth rules two trackers are compared under (TrackComparison.h).

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK track comparison
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Study/TrackComparison.h"
#include "GenericTracking/Study/TrackPulls.h"
#include "ToyBarrel3.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace o2::gtrk::study
{

namespace
{

/// Numbering differs from path order; an alternate direct link need not skip.
struct ReorderedDirect : detectors::ToyBarrel3 {
  static constexpr std::array<LayerLink, 3> layerLinks{{
    {.from = 0, .to = 2, .sourceMin = -1., .sourceMax = 1., .skipsLayer = false},
    {.from = 2, .to = 1, .sourceMin = -1., .sourceMax = 1., .skipsLayer = false},
    {.from = 0, .to = 1, .sourceMin = -1., .sourceMax = 1., .skipsLayer = true},
  }};
};

/// ToyBarrel3 with a link 0 -> 2 beside 0 -> 1 -> 2 that skips no layer:
/// whether a step on it leaves a hole is the pair's.
struct PastTheEnd : detectors::ToyBarrel3 {
  static constexpr std::array<LayerLink, 3> layerLinks{{
    {.from = 0, .to = 1, .sourceMin = -1.5, .sourceMax = 1.5, .skipsLayer = false},
    {.from = 0, .to = 2, .sourceMin = -1.5, .sourceMax = 1.5, .skipsLayer = false},
    {.from = 1, .to = 2, .sourceMin = -1.5, .sourceMax = 1.5, .skipsLayer = false},
  }};
};

[[maybe_unused]] const bool kReorderedDirectPublished = publishToy<ReorderedDirect>();
[[maybe_unused]] const bool kPastTheEndPublished = publishToy<PastTheEnd>();

/// Whether some report line contains \a text.
bool reports(const std::vector<std::string>& lines, const std::string& text)
{
  return std::ranges::any_of(lines, [&](const std::string& line) { return line.find(text) != std::string::npos; });
}

TruthRow row(LayerID layer, std::initializer_list<std::uint64_t> particles)
{
  TruthRow result;
  result.layer = layer;
  for (const auto particle : particles) {
    result.particles[result.count++] = particle;
  }
  return result;
}

/// A measurement of particle \a particle on \a layer at radius \a r and \a z.
TruthRow at(LayerID layer, std::uint64_t particle, float r, float z)
{
  auto result = row(layer, {particle});
  result.r = r;
  result.z = z;
  return result;
}

} // namespace

BOOST_AUTO_TEST_CASE(FindabilityFollowsTheLinksThatSkipNoLayer)
{
  const auto settings = comparisonSettings<ReorderedDirect>(2);
  // Particle 1 on 0 and 2, joined by 0 -> 2 despite the numbering gap;
  // particle 2 on 0 and 1, joined only by the skipping 0 -> 1; particle 3 on
  // all three, 0 -> 2 -> 1.
  const std::vector<TruthRow> rows{row(0, {1, 3}), row(2, {1, 3}), row(0, {2}), row(1, {2, 3})};
  const auto layers = particleLayers(rows, settings);
  BOOST_TEST(layers.at(1).chain == 2U);
  BOOST_TEST(layers.at(2).chain == 1U);
  BOOST_TEST(layers.at(3).chain == 3U);
  BOOST_TEST(layers.at(3).layers == 0b111U);
  BOOST_TEST(isFindable(1.F, layers.at(1), settings));
  BOOST_TEST(!isFindable(1.F, layers.at(2), settings));
  BOOST_TEST(!isFindable(0.05F, layers.at(1), settings));
  // The score and decay-family denominator use exactly the same predicate.
  const std::vector<TruthRow> decayRows{row(0, {1}), row(2, {1}), row(0, {2}), row(1, {2})};
  const ParticleTable particles{{1, {.pt = 1.F, .origin = ParticleOrigin::kDecay, .motherPdg = 310}}, {2, {.pt = 1.F, .origin = ParticleOrigin::kDecay, .motherPdg = 310}}};
  const auto score = scoreTracks(decayRows, particles, {{0, 1}, {2, 3}}, settings);
  BOOST_TEST(score.originFindable[static_cast<std::size_t>(ParticleOrigin::kDecay)].findable == 1U);
  BOOST_TEST(score.originFindable[static_cast<std::size_t>(ParticleOrigin::kDecay)].found == 1U);
  BOOST_TEST(score.decays[static_cast<std::size_t>(DecayFamily::kK0S)][0].onMinLayers == 2U);
  BOOST_TEST(score.decays[static_cast<std::size_t>(DecayFamily::kK0S)][0].findable == 1U);
}

BOOST_AUTO_TEST_CASE(AChainStepCrossingAnotherLayerIsAHole)
{
  // ToyBarrel3 with 0 -> 2 as a link that skips nothing: past layer 1's ends
  // (|z| 1 cm at r 4) a step 0 -> 2 is direct; through its length a hole.
  const auto settings = comparisonSettings<PastTheEnd>(2);
  BOOST_TEST(settings.chainLinks.size() == 3U);
  const auto step = std::ranges::find_if(settings.chainLinks, [](const ChainLink& link) { return link.from == 0 && link.to == 2; });
  BOOST_TEST((step != settings.chainLinks.end() && step->holes == 0b010U));
  // Particle 1 at z 0 on both: the segment crosses layer 1 at z 0. Particle 2
  // from z 0.9 to 3.9: at r 4 it is at z 1.9, past the end.
  const std::vector<TruthRow> rows{at(0, 1, 2.F, 0.F), at(2, 1, 8.F, 0.F), at(0, 2, 2.F, 0.9F), at(2, 2, 8.F, 3.9F)};
  const auto layers = particleLayers(rows, settings);
  BOOST_TEST(layers.at(1).chain == 1U);
  BOOST_TEST(layers.at(2).chain == 2U);
}

BOOST_AUTO_TEST_CASE(DirectChainsIncludeLayersAboveThirtyOne)
{
  const ComparisonSettings settings{.nLayers = 36, .minLayers = 3, .chainLinks = {{.from = 32, .to = 34}, {.from = 34, .to = 35}}};
  const std::vector<TruthRow> rows{row(32, {1}), row(34, {1}), row(35, {1})};
  const auto layers = particleLayers(rows, settings);
  BOOST_TEST(layers.at(1).chain == 3U);
  const ParticleTable particles{{1, {.pt = 1.F, .primary = true, .origin = ParticleOrigin::kPhysicalPrimary}}};
  const auto score = scoreTracks(rows, particles, {{0, 1, 2}}, settings);
  BOOST_TEST(score.originFindable[0].findable == 1U);
  BOOST_TEST(score.originFindable[0].complete == 1U);
}

BOOST_AUTO_TEST_CASE(ASharedMeasurementCountsForEveryParticleItCarries)
{
  // Particle 1 on layers 0-2, particle 2 on layers 0-2; the layer-1
  // measurement carries both. Row 6 is noise.
  const std::vector<TruthRow> rows{row(0, {1}), row(1, {1, 2}), row(2, {1}), row(0, {2}), row(2, {2}), row(1, {3}), row(2, {})};
  const ParticleTable particles{{1, {.pt = 1.F, .primary = true}}, {2, {.pt = 1.F, .primary = true}}, {3, {.pt = 1.F, .primary = true}}};
  const ComparisonSettings settings{.nLayers = 3, .minLayers = 2, .minPtGeV = 0.1, .chainLinks = {{.from = 0, .to = 1}, {.from = 1, .to = 2}}};

  // Ours: particle 1 through the shared measurement, good; particle 2 through
  // it too, good; particle 1 again with the noise row, a fake labelled 1.
  const TrackRows ours{{0, 1, 2}, {3, 1, 4}, {0, 5, 6}};
  // Theirs: particle 1 without the shared measurement's alternative -- the
  // same set as ours -- and particle 2 with one wrong measurement.
  const TrackRows theirs{{2, 1, 0}, {3, 5, 4}};

  const auto lines = compareTrackSets(rows, particles, ours, theirs, settings, "A", "B");
  BOOST_TEST(reports(lines, "A: 3 tracks, 33.33 % fake"));
  BOOST_TEST(reports(lines, "B: 2 tracks, 50.00 % fake"));
  // All three layers: particles 1 and 2. A finds both, B particle 1 only, and
  // particle 2 only with its fake.
  BOOST_TEST(reports(lines, "on all 3 layers: 2"));
  BOOST_TEST(reports(lines, "A: 100.00 % found (100.00 % complete, 100.0 % of their layers),  0.00 % only with a fake,  0.00 % cloned"));
  BOOST_TEST(reports(lines, "B:  50.00 % found (50.00 % complete, 100.0 % of their layers), 50.00 % only with a fake"));
  // Theirs track by track: {0,1,2} is identical to ours; {3,5,4} overlaps
  // {3,1,4} in two of three.
  BOOST_TEST(reports(lines, "identical 50.00 %"));
  BOOST_TEST(reports(lines, "overlapping 50.00 %"));
  BOOST_TEST(reports(lines, "found by both 1, by A only 1, by B only 0, by neither 0"));
}

BOOST_AUTO_TEST_CASE(TheScoreIsWhatTheComparisonReports)
{
  // As above: A finds particles 1 and 2, the second through the measurement
  // they share, and makes a fake labelled 1 with a noise row.
  const std::vector<TruthRow> rows{row(0, {1}), row(1, {1, 2}), row(2, {1}), row(0, {2}), row(2, {2}), row(1, {3}), row(2, {})};
  const ParticleTable particles{{1, {.pt = 1.F, .primary = true, .origin = ParticleOrigin::kPhysicalPrimary, .pdg = 211}}, {2, {.pt = 0.15F, .primary = true, .origin = ParticleOrigin::kDecay, .pdg = -2212}}, {3, {.pt = 1.F, .primary = true, .origin = ParticleOrigin::kMaterial, .pdg = 11}}};
  const ComparisonSettings settings{.nLayers = 3, .minLayers = 1, .minPtGeV = 0.1, .chainLinks = {{.from = 0, .to = 1}, {.from = 1, .to = 2}}};
  const TrackRows ours{{0, 1, 2}, {3, 1, 4}, {0, 5, 6}};
  const TrackRows theirs{{2, 1, 0}};

  const auto score = scoreTracks(rows, particles, ours, settings);
  BOOST_TEST(score.tracks == 3U);
  BOOST_TEST(score.fakes == 1U);
  BOOST_TEST(score.tracksByLength[3] == 3U);
  BOOST_TEST(score.fakesByLength[3] == 1U);
  BOOST_TEST(score.allLayers.findable == 2U);
  BOOST_TEST(score.allLayers.found == 2U);
  BOOST_TEST(score.allLayers.withFake == 1U); // particle 1, found too
  BOOST_TEST(score.allLayers.fakeOnly == 0U);
  BOOST_TEST(score.allLayers.foundByPt[1] == 1U); // particle 2, [0.1, 0.2)
  BOOST_TEST(score.allLayers.foundByPt[5] == 1U); // particle 1, [1, 2)
  // Particle 3, on layer 1 alone, is findable on one layer or more and not found.
  BOOST_TEST(score.minLayers.findable == 3U);
  BOOST_TEST(score.minLayers.found == 2U);
  // No row's first label counts for O2's QA here.
  BOOST_TEST(score.o2.findable == 0U);

  // By origin, primary or not: particle 1 a physical primary, 2 a decay's
  // (an antiproton), 3 the material's (an electron) and on one layer only.
  BOOST_TEST(score.originAllLayers[0].findable == 1U);
  BOOST_TEST(score.originAllLayers[1].found == 1U);
  BOOST_TEST(score.originAllLayers[2].findable == 0U);
  BOOST_TEST(score.originMinLayers[2].findable == 1U);
  BOOST_TEST(score.originMinLayers[2].found == 0U);
  BOOST_TEST((score.originMinLayers[1].bySpecies.at(2212) == std::array<std::uint64_t, 2>{1U, 1U}));
  BOOST_TEST((score.minLayers.bySpecies.at(11) == std::array<std::uint64_t, 2>{1U, 0U}));

  // Every line of the score is one the comparison prints, for the same tracker.
  const auto lines = compareTrackSets(rows, particles, ours, theirs, settings, "A", "B");
  for (const auto& line : describeScore(score, settings, "A")) {
    BOOST_TEST(std::ranges::count(lines, line) == 1, line);
  }
}

BOOST_AUTO_TEST_CASE(AStubIsFoundButNotComplete)
{
  // Particle 1 on three layers: track A holds two of them, found; track B all
  // three, complete. Particle 2 on three layers, found with a stub only.
  const std::vector<TruthRow> rows{row(0, {1}), row(1, {1}), row(2, {1}), row(0, {2}), row(1, {2}), row(2, {2})};
  const ParticleTable particles{{1, {.pt = 1.F, .primary = true, .origin = ParticleOrigin::kPhysicalPrimary, .pdg = 211}}, {2, {.pt = 0.15F, .primary = true, .origin = ParticleOrigin::kPhysicalPrimary, .pdg = 211}}};
  const ComparisonSettings settings{.nLayers = 3, .minLayers = 2, .minPtGeV = 0.1, .chainLinks = {{.from = 0, .to = 1}, {.from = 1, .to = 2}}};

  const auto stubs = scoreTracks(rows, particles, {{0, 1}, {3, 4}}, settings);
  BOOST_TEST(stubs.allLayers.findable == 2U);
  BOOST_TEST(stubs.allLayers.found == 2U);
  BOOST_TEST(stubs.allLayers.complete == 0U);
  BOOST_TEST(stubs.allLayers.layersFound == 4U);
  BOOST_TEST(stubs.allLayers.layersOf == 6U);
  BOOST_TEST(stubs.allLayers.layerFraction == 4. / 3., boost::test_tools::tolerance(1e-9)); // two stubs of 2 of 3 layers
  BOOST_TEST(stubs.allLayers.layerFractionByPt[1] == 2. / 3., boost::test_tools::tolerance(1e-9));

  const auto mixed = scoreTracks(rows, particles, {{0, 1}, {0, 1, 2}, {3, 4}}, settings);
  BOOST_TEST(mixed.allLayers.found == 2U);
  BOOST_TEST(mixed.allLayers.cloned == 1U);
  BOOST_TEST(mixed.allLayers.complete == 1U); // particle 1, by its longer track
  BOOST_TEST(mixed.allLayers.completeByPt[5] == 1U);
  BOOST_TEST(mixed.allLayers.completeByPt[1] == 0U);
  BOOST_TEST(mixed.allLayers.layersFound == 5U);
  BOOST_TEST(mixed.allLayers.layerFraction == 5. / 3., boost::test_tools::tolerance(1e-9)); // the complete one 1, the stub 2 of 3; the clone nothing
  BOOST_TEST(reports(describeScore(mixed, settings, "A"), "100.00 % found (50.00 % complete, 83.3 % of their layers)"));
}

BOOST_AUTO_TEST_CASE(SecondariesAndSoftParticlesAreNotFindable)
{
  const std::vector<TruthRow> rows{row(0, {1}), row(1, {1}), row(0, {2}), row(1, {2})};
  const ParticleTable particles{{1, {.pt = 1.F, .primary = false}}, {2, {.pt = 0.05F, .primary = true}}};
  const TrackRows tracks{{0, 1}, {2, 3}};
  const auto lines = compareTrackSets(rows, particles, tracks, tracks, {.nLayers = 2, .minLayers = 2, .minPtGeV = 0.1, .chainLinks = {{.from = 0, .to = 1}}}, "A", "B");
  BOOST_TEST(reports(lines, "on all 2 layers: 0"));
  BOOST_TEST(reports(lines, "A: 2 tracks, 0.00 % fake"));
}

BOOST_AUTO_TEST_CASE(TwoGoodTracksOfOneParticleAreAClone)
{
  const std::vector<TruthRow> rows{row(0, {1}), row(1, {1}), row(0, {1}), row(1, {1})};
  const ParticleTable particles{{1, {.pt = 0.4F, .primary = true}}};
  const TrackRows ours{{0, 1}, {2, 3}};
  const TrackRows theirs{{0, 1}};
  const auto lines = compareTrackSets(rows, particles, ours, theirs, {.nLayers = 2, .minLayers = 2, .minPtGeV = 0.1, .chainLinks = {{.from = 0, .to = 1}}}, "A", "B");
  BOOST_TEST(reports(lines, "A: 100.00 % found (100.00 % complete, 100.0 % of their layers),  0.00 % only with a fake, 100.00 % cloned"));
  BOOST_TEST(reports(lines, "B: 100.00 % found (100.00 % complete, 100.0 % of their layers),  0.00 % only with a fake,  0.00 % cloned"));
}

BOOST_AUTO_TEST_CASE(OneWrongMeasurementIsAFakeAndItsLayerIsCounted)
{
  // Particle 1 on layers 0-2, with its own layer-1 measurement passed over for
  // particle 2's.
  const std::vector<TruthRow> rows{row(0, {1}), row(1, {1}), row(2, {1}), row(1, {2})};
  const ParticleTable particles{{1, {.pt = 1.F, .primary = true}}, {2, {.pt = 1.F, .primary = false}}};
  const TrackRows ours{{0, 3, 2}};
  const auto lines = compareTrackSets(rows, particles, ours, {}, {.nLayers = 3, .minLayers = 2, .minPtGeV = 0.1, .chainLinks = {{.from = 0, .to = 1}, {.from = 1, .to = 2}}}, "A", "B");
  BOOST_TEST(reports(lines, "A fakes: 1; labelled with a primary 100.0 %"));
  BOOST_TEST(reports(lines, "by wrong measurements 1 / 2 / 3 / 4+: 100.0 / 0.0 / 0.0 / 0.0 %"));
  BOOST_TEST(reports(lines, "1 wrong measurements, by layer (%): 0: 0.0 1: 100.0 2: 0.0"));
  BOOST_TEST(reports(lines, "passed over 100.0 %, none of its own there 0.0 %"));
  BOOST_TEST(reports(lines, "a secondary's 100.0 %"));
  BOOST_TEST(goodParticle(rows, std::vector<std::uint32_t>{0, 1, 2}) == 1U);
  BOOST_TEST(goodParticle(rows, std::vector<std::uint32_t>{0, 3, 2}) == kNoTrackParticle);
}

BOOST_AUTO_TEST_CASE(ShortTracksAreCountedByTheirLayers)
{
  // Particles 1 and 2 on layers 0-2, particle 3 on layers 0, 1 and 3; a
  // four-layer detector, so every track is short of a layer.
  const std::vector<TruthRow> rows{row(0, {1}), row(1, {1}), row(2, {1}), row(0, {2}), row(1, {2}), row(2, {2}), row(0, {3}), row(1, {3}), row(3, {3})};
  const ParticleTable particles{{1, {.pt = 1.F, .primary = true}}, {2, {.pt = 1.F, .primary = true}}, {3, {.pt = 1.F, .primary = true}}};
  // Two good tracks on 0-2, one fake with a hole on layer 2, and the one good
  // track of the full set: none.
  const TrackRows ours{{0, 1, 2}, {3, 4, 5}, {0, 7, 8}};
  const auto lines = compareTrackSets(rows, particles, ours, {}, {.nLayers = 4, .minLayers = 3, .minPtGeV = 0.1, .chainLinks = {{.from = 0, .to = 1}, {.from = 1, .to = 2}, {.from = 2, .to = 3}}}, "A", "B");
  BOOST_TEST(reports(lines, "A: 3 tracks short of a layer, by layers (fake %): 012: 2 (0.0) 013: 1 (100.0)"));
  BOOST_TEST(reports(lines, "B: 0 tracks short of a layer, by layers (fake %):"));
}

BOOST_AUTO_TEST_CASE(O2sQADefinitionReadsOnlyTheFirstLabel)
{
  // Particles 1 and 2 on both layers of a two-layer detector, at 50 MeV --
  // below our findable floor, not O2's. On layer 1 particle 2's measurement
  // is shared and its first label is particle 3's, so O2's QA does not count
  // it findable; ours, reading every label, would above the floor.
  auto first = [](TruthRow row) {
    row.firstLabelCounts = true;
    return row;
  };
  const std::vector<TruthRow> rows{first(row(0, {1})), first(row(1, {1})), first(row(0, {2})), first(row(1, {3, 2}))};
  const ParticleTable particles{{1, {.pt = 0.05F, .primary = true, .eta = 0.2F, .phi = 1.F}}, {2, {.pt = 0.05F, .primary = true, .eta = -0.7F, .phi = 4.F}}, {3, {.pt = 1.F, .primary = false}}};
  const TrackRows ours{{0, 1}, {0, 1}, {2, 3}};
  const auto lines = compareTrackSets(rows, particles, ours, {}, {.nLayers = 2, .minLayers = 2, .minPtGeV = 0.1, .chainLinks = {{.from = 0, .to = 1}}}, "A", "B");
  BOOST_TEST(reports(lines, "first labels on all 2 layers, any pT: 1 findable"));
  BOOST_TEST(reports(lines, "A: found 1 (100.00 %), clone tracks 1, with a fake track 0"));
  BOOST_TEST(reports(lines, "on all 2 layers: 0"));
}

BOOST_AUTO_TEST_CASE(DecayDaughtersAreSplitByMotherAndProduction)
{
  // A five-layer detector, production bins [0,1), [1,5), [5,inf) cm.
  // Particle 1, a K0S's at 3 cm, on layers 2-4: findable, found. Particle 2,
  // an antilambda's at 0.5 cm, on layers 0, 1, 3, 4: on three layers or more
  // but no three adjacent, found. Particle 3, a pion's in flight at |z| 2 and
  // radius 10 cm, on layers 0-2: findable, not found. Particle 4, a primary,
  // and particle 5, a K0S's below the pT floor, are not decay daughters here.
  const std::vector<TruthRow> rows{row(2, {1}), row(3, {1}), row(4, {1}), row(0, {2}), row(1, {2}), row(3, {2}), row(4, {2}), row(0, {3}), row(1, {3}), row(2, {3}), row(0, {4}), row(1, {4}), row(2, {4}), row(2, {5}), row(3, {5}), row(4, {5})};
  const ParticleTable particles{{1, {.pt = 0.5F, .origin = ParticleOrigin::kDecay, .pdg = 211, .motherPdg = 310, .vx = 3.F}}, {2, {.pt = 0.5F, .origin = ParticleOrigin::kDecay, .pdg = -2212, .motherPdg = -3122, .vx = 0.3F, .vy = 0.4F}}, {3, {.pt = 0.5F, .origin = ParticleOrigin::kDecay, .pdg = 13, .motherPdg = 211, .vx = 10.F, .vz = -2.F}}, {4, {.pt = 0.5F, .primary = true, .origin = ParticleOrigin::kPhysicalPrimary, .pdg = 211}}, {5, {.pt = 0.05F, .origin = ParticleOrigin::kDecay, .pdg = 211, .motherPdg = 310, .vx = 3.F}}};
  const TrackRows tracks{{0, 1, 2}, {3, 4, 5, 6}, {10, 11, 12}};
  ComparisonSettings settings{
    .nLayers = 5,
    .minLayers = 3,
    .minPtGeV = 0.1,
    .productionEdges = {1., 5.},
    .chainLinks = {{.from = 0, .to = 1}, {.from = 1, .to = 2}, {.from = 2, .to = 3}, {.from = 3, .to = 4}},
  };

  const auto score = scoreTracks(rows, particles, tracks, settings);
  const auto counts = [&](DecayFamily family, std::size_t bin) {
    const auto& c = score.decays[static_cast<std::size_t>(family)][bin];
    return std::array<std::uint64_t, 4>{c.onMinLayers, c.foundOnMinLayers, c.findable, c.found};
  };
  BOOST_TEST(score.decays[0].size() == 3U);
  BOOST_TEST((counts(DecayFamily::kK0S, 1) == std::array<std::uint64_t, 4>{1U, 1U, 1U, 1U}));
  BOOST_TEST((counts(DecayFamily::kLambda, 0) == std::array<std::uint64_t, 4>{1U, 1U, 0U, 0U}));
  BOOST_TEST((counts(DecayFamily::kPion, 2) == std::array<std::uint64_t, 4>{1U, 0U, 1U, 0U}));
  BOOST_TEST((counts(DecayFamily::kOther, 0) == std::array<std::uint64_t, 4>{}));

  const auto lines = compareTrackSets(rows, particles, tracks, {}, settings, "A", "B");
  BOOST_TEST(reports(lines, "production radius [cm] (bins from 0, 1.0, 5.0)"));
  BOOST_TEST(reports(lines, "A V0s and hyperons 100.00 (1/2), of all 100.00 %; by production: 0.0 (0/1) 100.0 (1/1) 0.0 (0/0)"));
  BOOST_TEST(reports(lines, "A in flight          0.00 (1/1), of all 0.00 %; by production: 0.0 (0/0) 0.0 (0/0) 0.0 (1/1)"));
  BOOST_TEST(reports(lines, "B K0S                0.00 (1/1)"));

  // Along z, particle 3 moves to the middle bin.
  settings.productionAlongZ = true;
  const auto alongZ = scoreTracks(rows, particles, tracks, settings);
  BOOST_TEST(alongZ.decays[static_cast<std::size_t>(DecayFamily::kPion)][1].findable == 1U);
}

BOOST_AUTO_TEST_CASE(PullWidthsAreRobust)
{
  PullTable table{"t", {{.name = "x", .unit = "um", .scale = 1e4}}};
  // Pulls -2 .. 2 in one bin: median 0, median absolute deviation 1; a
  // thousand-sigma outlier does not move either much.
  for (const double residual : {-2., -1., 0., 1., 2., 1000.}) {
    table.add(0.5, std::array<double, 1>{residual * 1e-4}, std::array<double, 1>{1e-4});
  }
  table.add(0.5, std::array<double, 1>{0.}, std::array<double, 1>{1e-4});
  const auto lines = table.report();
  BOOST_TEST(reports(lines, "7 residual median     0.00 sigma     1.48; pull median  0.000 sigma  1.483 beyond 5 14.29 %"));
}

} // namespace o2::gtrk::study

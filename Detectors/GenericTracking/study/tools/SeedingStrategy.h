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

/// \file SeedingStrategy.h
/// \brief The truth-removal seeding study: the configured iterations as
///        seeding passes, each on the rows the earlier passes' seeded
///        particles left, scored against the MC truth.
///
/// A pass is an iteration's graph alone (`Tracker::runGraph`): tracklets to
/// the numbered roads, with no road fit, selection, refit or follow. A road
/// seeds every particle all of its rows carry among their labels (the one
/// truth rule, as the funnel's `forEachOwner`); a road no particle owns is
/// impure. After the pass every row of every particle it seeded is held
/// exclusively -- a perfect follow -- whatever the particle's origin, and the
/// next pass runs on what is left. Without removal every pass reads the whole
/// frame (the check against the real run's first iteration).
///
/// The denominator is the score's `originFindable` population of physical
/// primaries (`isFindable`: a chain of `minLayers` layers, the score's pT
/// floor), counted at each of the study's pT steps -- or, by the study's
/// population, the findable decay daughters, or only those of a K0S or a
/// hyperon. Coverage is cumulative:
/// the share of the findable at or above a step that some pass so far seeded.
/// Per pass: the graph's sizes (tracklets, triplets, neighbours, roads; the
/// triplet and neighbour fits in a counters build), its wall time and pool
/// peak, the impure roads, the owners per seeded particle (the multiplicity),
/// the particles seeded by origin, the pass's own gain and the cumulative
/// coverage by class at its floor -- the barrel, the discs on one side by the
/// sign of their z, both sides -- the layer patterns of the findable at or
/// above that floor still unseeded, and what the rows left free are: noise,
/// the findable's, other particles' at or above the lowest step,
/// or only particles below it.
///
/// A shadow pass runs on the same free rows as the
/// pass after it and holds nothing: what it seeds and the next pass does not,
/// and the reverse, are the two settings' differences on identical rows, each
/// with the other's mechanism -- the particle's rows inside a road with
/// foreign rows (at least the pass's minimum of its own), only partly in a
/// road, in none, or (with the fit) in a pure road the fit refused. With the
/// fit, a particle is seeded by an accepted road
/// fit record all of whose rows carry it (`Tracker::runGraph`'s fit), and the
/// findable with a pure road but no pure record are counted by the length of
/// their longest pure road, and by why the fit refuses that road
/// (`roadFitRefusal`, the road fit stage's tests replayed).
///
/// After a real run's first iterations (`startAfter`): their
/// good tracks seed their particles before the first pass, and every pass
/// starts from the used mask they left, so that the passes are judged on the
/// rows such a run leaves to an appended iteration.
///
/// A real run's iterations are scored the same way (`observeReal`): each
/// iteration's good tracks (O2's rule, `goodParticle`) seed their particles, the used mask after it is the held rows, and the
/// gap to truth removal is counted -- the rows of the particles found so far
/// the iteration left free (an incomplete follow), and the held rows that
/// carry no found particle (fakes' rows).

#ifndef GTRK_STUDY_SEEDINGSTRATEGY_H_
#define GTRK_STUDY_SEEDINGSTRATEGY_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/RoadCount.h"
#include "GenericTracking/Study/Statistics.h"
#include "GenericTracking/Study/TrackComparison.h"
#include "GenericTracking/Study/TrackingInput.h"
#include "GenericTracking/Study/Validation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <functional>
#include <map>
#include <optional>
#include <ostream>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace o2::gtrk::study
{

/// The particles the study scores.
enum class SeedingPopulation : std::uint8_t { kPrimaries,
                                              kDecays,
                                              kV0s };
inline constexpr std::array<std::string_view, 3> kSeedingPopulationNames{"primaries", "decay daughters", "V0 daughters"};

/// \a particle belongs to \a population: a physical primary; a decay's
/// daughter; one whose mother is a K0S or a hyperon (`decayFamily`: a K0S,
/// Lambda, Sigma or cascade).
[[nodiscard]] inline bool inPopulation(const study::ParticleTruth& particle, SeedingPopulation population)
{
  switch (population) {
    case SeedingPopulation::kPrimaries:
      return particle.origin == study::ParticleOrigin::kPhysicalPrimary;
    case SeedingPopulation::kDecays:
      return particle.origin == study::ParticleOrigin::kDecay;
    case SeedingPopulation::kV0s: {
      const auto family = study::decayFamily(particle.motherPdg);
      return particle.origin == study::ParticleOrigin::kDecay &&
             (family == study::DecayFamily::kK0S || family == study::DecayFamily::kLambda || family == study::DecayFamily::kSigma ||
              family == study::DecayFamily::kCascade);
    }
  }
  return false;
}

template <DetectorLayout Det>
class SeedingStrategy
{
 public:
  /// The road fit stage replayed on one road's rows, innermost first.
  using Refit = std::function<FitRefusal(std::span<const std::uint32_t>)>;

  /// What a pass's graph cost, beyond the sizes `observe` counts.
  struct PassCost {
    double seconds = 0.;            ///< the graph's wall time, the truth's bookkeeping left out
    std::size_t poolPeakBytes = 0U; ///< the pool's peak inside the pass, from a released pool
    bool counted = false;           ///< the fits below were counted (a counters build)
    std::uint64_t tripletFits = 0U;
    std::uint64_t neighbourFits = 0U;
    double fitSeconds = 0.; ///< with the fit: its wall time, inside `seconds`
  };

  /// The particles of \a truth, the findable of \a comparison, coverage at
  /// \a steps [GeV]; with \a remove, each pass's seeded particles' rows held;
  /// with \a fit, seeded by the road fit's records rather than the roads;
  /// the findable of \a population.
  SeedingStrategy(const study::TrackingTruth& truth, const study::ComparisonSettings& comparison, std::vector<double> steps, bool remove, bool fit,
                  SeedingPopulation population = SeedingPopulation::kPrimaries)
    : mTruth(truth),
      mComparison(comparison),
      mSteps(std::move(steps)),
      mRemove(remove),
      mFit(fit),
      mPopulation(kSeedingPopulationNames[static_cast<std::size_t>(population)])
  {
    const auto& rows = truth.rows;
    // Every particle the rows carry, numbered, and the rows of each.
    std::vector<std::size_t> counts;
    for (const auto& row : rows) {
      for (std::uint8_t at = 0U; at < row.count; ++at) {
        const auto [found, inserted] = mIndex.try_emplace(row.particles[at], static_cast<std::uint32_t>(mIds.size()));
        if (inserted) {
          mIds.push_back(row.particles[at]);
          counts.push_back(0U);
        }
        ++counts[found->second];
      }
    }
    mParticleOffsets.assign(mIds.size() + 1U, 0U);
    for (std::size_t particle = 0; particle < mIds.size(); ++particle) {
      mParticleOffsets[particle + 1U] = mParticleOffsets[particle] + counts[particle];
    }
    mParticleRows.resize(mParticleOffsets.back());
    std::vector<std::size_t> filled(mParticleOffsets.begin(), mParticleOffsets.end() - 1);
    for (std::uint32_t row = 0; row < rows.size(); ++row) {
      for (std::uint8_t at = 0U; at < rows[row].count; ++at) {
        mParticleRows[filled[mIndex.at(rows[row].particles[at])]++] = row;
      }
    }
    mRowOffsets.assign(rows.size() + 1U, 0U);
    for (std::size_t row = 0; row < rows.size(); ++row) {
      mRowOffsets[row + 1U] = mRowOffsets[row] + rows[row].count;
    }
    mRowParticles.resize(mRowOffsets.back());
    for (std::size_t row = 0; row < rows.size(); ++row) {
      for (std::uint8_t at = 0U; at < rows[row].count; ++at) {
        mRowParticles[mRowOffsets[row] + at] = mIndex.at(rows[row].particles[at]);
      }
    }
    mHeld.assign(rows.size(), std::uint8_t{0});
    mPt.assign(mIds.size(), 0.F);
    for (std::size_t particle = 0; particle < mIds.size(); ++particle) {
      const auto found = truth.particles.find(mIds[particle]);
      mPt[particle] = found == truth.particles.end() ? 0.F : found->second.pt;
    }
    mIsFindable.assign(mIds.size(), std::uint8_t{0});
    mSeededBy.assign(mIds.size(), kUnseeded);
    mSeededIn.assign(mIds.size(), kUnseeded);
    mShadowIn.assign(mIds.size(), kUnseeded);
    mShadowBy.assign(mIds.size(), kUnseeded);
    mRoadIn.assign(mIds.size(), kUnseeded);
    mRoadBy.assign(mIds.size(), kUnseeded);
    mEver.assign(mIds.size(), Mechanism::kNone);
    mRefusal.assign(mIds.size(), std::int8_t{-1});
    mRefusals.assign(mIds.size(), std::uint8_t{0});
    for (std::size_t layer = 0; layer < comparison.layerSpans.size(); ++layer) {
      if (((comparison.barrelLayers >> layer) & 1U) == 0U) {
        (comparison.layerSpans[layer].position > 0.F ? mPositiveDiscs : mNegativeDiscs) |= std::uint64_t{1} << layer;
      }
    }
    for (const auto& [id, layers] : study::particleLayers(rows, comparison)) {
      const auto found = truth.particles.find(id);
      if (found == truth.particles.end() || !inPopulation(found->second, population) || !study::isFindable(found->second.pt, layers, comparison)) {
        continue;
      }
      const bool positive = (layers.layers & mPositiveDiscs) != 0U;
      const bool negative = (layers.layers & mNegativeDiscs) != 0U;
      const auto kind = positive && negative ? Class::kDiscBoth : (positive ? Class::kDiscPositive : (negative ? Class::kDiscNegative : Class::kBarrel));
      mFindable.push_back({.particle = mIndex.at(id), .pt = found->second.pt, .layers = layers.layers, .kind = kind});
      mIsFindable[mIndex.at(id)] = 1U;
    }
  }

  /// A pass starts: its counts at zero; a \a shadow pass holds nothing.
  void beginPass(bool shadow)
  {
    mPass = {};
    mSeeded.clear();
    mShadow = shadow;
    mBest.assign(mIds.size(), std::uint8_t{0});
    mPureLength.assign(mIds.size(), std::uint8_t{0});
    mPureRows.clear();
  }

  /// One slice of the pass: its roads, walked over its graph, and with the
  /// fit its \a records.
  template <typename Result, typename Roads, typename Records>
  void observe(const Result& result, const Roads& roads, const Records* records)
  {
    constexpr std::size_t kMax = maxChainLayers<Det>;
    const auto tracklets = result.tracklets->view();
    const auto triplets = result.triplets->view();
    mPass.tracklets += result.tracklets->size();
    mPass.triplets += result.triplets->size();
    mPass.neighbours += result.neighbours->size();
    const auto walk = roads.walk(*result.triplets, *result.neighbours);
    std::array<std::uint32_t, kMax - 2U> chain{};
    std::array<std::uint32_t, kMax> rows{};
    for (std::uint32_t road = 0; road < roads.roadCount(); ++road) {
      const auto length = walk(road, chain.data());
      rows[0] = tracklets.source[triplets.first[chain[0]]];
      rows[1] = tracklets.target[triplets.first[chain[0]]];
      for (std::uint32_t i = 0; i < length; ++i) {
        rows[i + 2U] = tracklets.target[triplets.second[chain[i]]];
      }
      const std::span<const std::uint32_t> held{rows.data(), length + 2U};
      ++mPass.roads;
      noteOwnRows(held);
      const bool owned = forEachOwner(held, [&](std::uint32_t particle) {
        ++mPass.owners;
        if (mFit && !mShadow && mIsFindable[particle] != 0U && held.size() > mPureLength[particle]) {
          mPureRows[particle].assign(held.begin(), held.end());
        }
        mPureLength[particle] = std::max(mPureLength[particle], static_cast<std::uint8_t>(held.size()));
        if (!mShadow) {
          mRoadIn[particle] = mPassIndex;
          mRoadBy[particle] = mRoadBy[particle] == kUnseeded ? mPassIndex : mRoadBy[particle];
        }
        if (!mFit) {
          seed(particle);
        }
      });
      mPass.pureRoads += owned ? 1U : 0U;
      mPass.impureRoads += owned ? 0U : 1U;
    }
    if (records == nullptr) {
      return;
    }
    const auto view = records->view();
    mPass.records += records->size();
    for (std::size_t record = 0; record < records->size(); ++record) {
      const std::span<const std::uint32_t> held{view.clusters[record].index, static_cast<std::size_t>(view.count[record])};
      const bool owned = forEachOwner(held, [&](std::uint32_t particle) { seed(particle); });
      mPass.pureRecords += owned ? 1U : 0U;
    }
  }

  /// Before the first pass: a real run's first \a iterations iterations left
  /// \a tracks and the \a used mask. Their good tracks seed
  /// their particles (`seededBy` -2), the rows the mask holds exclusively are
  /// held, and the passes are numbered on from \a iterations; reported to
  /// \a out.
  template <typename Tracks>
  void startAfter(std::size_t iterations, const Tracks& tracks, std::span<const std::uint8_t> used, std::ostream& out)
  {
    std::uint64_t fakes = 0U;
    for (const auto& track : tracks) {
      const auto id = study::goodParticle(mTruth.rows, std::span<const std::uint32_t>{track.clusters.index, static_cast<std::size_t>(track.count)});
      if (id == study::kNoTrackParticle) {
        ++fakes;
        continue;
      }
      const auto particle = mIndex.at(id);
      mSeededBy[particle] = kByRun;
      mSeededIn[particle] = kByRun;
    }
    std::uint64_t held = 0U;
    for (std::size_t row = 0; row < mHeld.size(); ++row) {
      mHeld[row] = used[row] == kExclusive ? 1U : 0U;
      held += mHeld[row];
    }
    mPassIndex = static_cast<std::uint32_t>(iterations);
    std::string line =
      std::format("after the run's first {} iterations: {} tracks, {} fake; {} rows held ({:.2f} %); found of the findable {} at or above:", iterations,
                  tracks.size(), fakes, held, study::percent(held, mHeld.size()), mPopulation);
    for (const auto step : mSteps) {
      std::uint64_t findable = 0U;
      std::uint64_t found = 0U;
      for (const auto& particle : mFindable) {
        if (static_cast<double>(particle.pt) >= step) {
          ++findable;
          found += mSeededBy[particle.particle] != kUnseeded ? 1U : 0U;
        }
      }
      line += std::format(" {} {:.2f}", step, study::percent(found, findable));
    }
    out << line << '\n';
  }

  /// Iteration \a pass's tracks of a real run, \a first to \a end of \a view,
  /// with \a used the mask after it and \a roads its roads: the pass as
  /// `observe` would have it, its good tracks seeding; `endPass` reports it
  /// with the real-run counts (a study made without removal).
  template <typename TrackView>
  void observeReal(const TrackView& view, std::size_t first, std::size_t end, std::span<const std::uint8_t> used, std::uint64_t roads)
  {
    beginPass(false);
    mPass.roads = roads;
    mReal = RealCounts{};
    for (std::size_t track = first; track < end; ++track) {
      const std::span<const std::uint32_t> rows{view.clusters[track].index, static_cast<std::size_t>(view.count[track])};
      ++mReal->tracks;
      const auto id = study::goodParticle(mTruth.rows, rows);
      if (id == study::kNoTrackParticle) {
        ++mReal->fakes;
        continue;
      }
      seed(mIndex.at(id));
    }
    for (std::size_t row = 0; row < mHeld.size(); ++row) {
      mHeld[row] = used[row] == kExclusive ? 1U : 0U;
      mReal->heldRows += mHeld[row];
    }
    for (std::size_t particle = 0; particle < mIds.size(); ++particle) {
      if (mSeededBy[particle] == kUnseeded) {
        continue;
      }
      std::uint64_t free = 0U;
      for (auto at = mParticleOffsets[particle]; at < mParticleOffsets[particle + 1U]; ++at) {
        free += mHeld[mParticleRows[at]] == 0U ? 1U : 0U;
      }
      mReal->foundRows += mParticleOffsets[particle + 1U] - mParticleOffsets[particle];
      mReal->foundRowsFree += free;
      mReal->foundWithFree += free != 0U ? 1U : 0U;
      mReal->findableWithFree += free != 0U && mIsFindable[particle] != 0U ? 1U : 0U;
    }
    for (std::size_t row = 0; row < mHeld.size(); ++row) {
      if (mHeld[row] == 0U) {
        continue;
      }
      bool found = false;
      for (auto at = mRowOffsets[row]; at < mRowOffsets[row + 1U] && !found; ++at) {
        found = mSeededBy[mRowParticles[at]] != kUnseeded;
      }
      mReal->wrongHeld += found ? 0U : 1U;
    }
  }

  /// The pass ends: its report to \a out and \a json, and the rows the next
  /// pass must not read -- every row of the particles it seeded, with removal.
  /// With the fit, \a refit says why it refused the pure roads of the
  /// findable without a pure record.
  [[nodiscard]] std::vector<std::uint32_t> endPass(const IterationSettings& settings, const PassCost& cost, std::ostream& out, std::ostream* json,
                                                   const Refit& refit = {})
  {
    std::vector<std::uint32_t> hold;
    if (mRemove && !mShadow) {
      for (const auto particle : mSeeded) {
        const std::span<const std::uint32_t> rows{mParticleRows};
        hold.insert(hold.end(), rows.begin() + static_cast<std::ptrdiff_t>(mParticleOffsets[particle]),
                    rows.begin() + static_cast<std::ptrdiff_t>(mParticleOffsets[particle + 1U]));
      }
      std::ranges::sort(hold);
      const auto [last, end] = std::ranges::unique(hold);
      hold.erase(last, end);
      for (const auto row : hold) {
        mHeld[row] = 1U;
      }
    }
    const double floor = settings.windows.minPtGeV;
    // The particles the pass seeded, by origin: the findable among them
    // below.
    std::array<std::uint64_t, study::kParticleOrigins> seededByOrigin{};
    for (const auto particle : mSeeded) {
      const auto found = mTruth.particles.find(mIds[particle]);
      ++seededByOrigin[static_cast<std::size_t>(found == mTruth.particles.end() ? study::ParticleOrigin::kMaterial : found->second.origin)];
    }
    std::uint64_t seededFindable = 0U;
    std::uint64_t firstFindable = 0U;
    std::vector<std::uint64_t> covered(mSteps.size(), 0U);
    std::vector<std::uint64_t> findable(mSteps.size(), 0U);
    std::array<std::uint64_t, kClasses> classFindable{};
    std::array<std::uint64_t, kClasses> classCovered{};
    std::array<std::uint64_t, kClasses> classGained{};
    std::map<std::uint64_t, std::uint64_t> unseeded;
    const auto& seededIn = mShadow ? mShadowIn : mSeededIn;
    // The two settings' differences on identical rows: this pass against the
    // shadow pass just before it, each difference by the other's mechanism.
    const bool versus = !mShadow && mShadowPass != kUnseeded && mShadowPass + 1U == mPassIndex;
    std::array<std::uint64_t, kMechanisms> mainOnly{};
    std::array<std::uint64_t, kMechanisms> shadowOnly{};
    std::uint64_t bothSeeded = 0U;
    // With the fit: the findable with a pure road and no pure record, by the
    // length of their longest pure road.
    std::vector<std::uint64_t> lostByLength(maxChainLayers<Det> + 1U, 0U);
    std::array<std::uint64_t, kClasses> lostByClass{};
    std::array<std::uint64_t, kFitRefusals> lostBy{};
    std::uint64_t roadFindable = 0U;
    std::array<std::uint64_t, kMechanisms> unseededBy{};
    for (const auto& particle : mFindable) {
      const auto number = particle.particle;
      if (!mShadow && mSeededIn[number] != mPassIndex) {
        const auto mechanism = mechanismOf(mPureLength[number], mBest[number], settings.minClusters);
        mEver[number] = rank(mechanism) < rank(mEver[number]) ? mechanism : mEver[number];
      }
      if (versus) {
        const bool main = mSeededIn[number] == mPassIndex;
        const bool shadow = mShadowIn[number] == mShadowPass;
        bothSeeded += main && shadow ? 1U : 0U;
        if (main && !shadow) {
          ++mainOnly[static_cast<std::size_t>(mechanismOf(mShadowPure[number], mShadowBest[number], mShadowMinClusters))];
        } else if (shadow && !main) {
          ++shadowOnly[static_cast<std::size_t>(mechanismOf(mPureLength[number], mBest[number], settings.minClusters))];
        }
      }
      if (mFit && !mShadow && mRoadIn[number] == mPassIndex) {
        ++roadFindable;
        if (mSeededIn[number] != mPassIndex) {
          ++lostByLength[mPureLength[number]];
          ++lostByClass[static_cast<std::size_t>(particle.kind)];
          if (refit) {
            const auto why = refit(mPureRows.at(number));
            ++lostBy[static_cast<std::size_t>(why)];
            mRefusal[number] = static_cast<std::int8_t>(why);
            ++mRefusals[number];
          }
        }
      }
    }
    for (const auto& particle : mFindable) {
      const bool seeded = mSeededBy[particle.particle] != kUnseeded;
      const bool first = mSeededBy[particle.particle] == mPassIndex;
      seededFindable += seededIn[particle.particle] == mPassIndex ? 1U : 0U;
      firstFindable += first ? 1U : 0U;
      for (std::size_t step = 0; step < mSteps.size(); ++step) {
        if (static_cast<double>(particle.pt) >= mSteps[step]) {
          ++findable[step];
          covered[step] += seeded ? 1U : 0U;
        }
      }
      if (static_cast<double>(particle.pt) >= floor) {
        ++classFindable[static_cast<std::size_t>(particle.kind)];
        classCovered[static_cast<std::size_t>(particle.kind)] += seeded ? 1U : 0U;
        classGained[static_cast<std::size_t>(particle.kind)] += first ? 1U : 0U;
        if (!seeded) {
          ++unseeded[particle.layers];
          ++unseededBy[static_cast<std::size_t>(mEver[particle.particle])];
        }
      }
    }
    const double multiplicity = mSeeded.empty() ? 0. : static_cast<double>(mPass.owners) / static_cast<double>(mSeeded.size());
    out << std::format(
             "pass {}{}: floor {} GeV, layers {} (mask {}), starts {}, {} clusters, kink sigma {}; {} tracklets, {} triplets{}, {} neighbours{}, {} roads "
             "({:.2f} % impure, {:.2f} owners per seeded particle); {:.3f} s, pool peak {:.2f} GB\n",
             mPassIndex, mShadow ? " (shadow: holds nothing)" : "", floor, layerList(settings.layerMask), settings.layerMask,
             settings.roadStartLayers == 0U ? "any" : layerList(settings.roadStartLayers), settings.minClusters,
             settings.graphKinkSigma > 0. ? std::format("{} (graph {})", settings.kinkSigma, settings.graphKinkSigma) : std::format("{}", settings.kinkSigma),
             mPass.tracklets, mPass.triplets, cost.counted ? std::format(" ({} fits)", cost.tripletFits) : std::string{}, mPass.neighbours,
             cost.counted ? std::format(" ({} fits)", cost.neighbourFits) : std::string{}, mPass.roads, study::percent(mPass.impureRoads, mPass.roads), multiplicity,
             cost.seconds, static_cast<double>(cost.poolPeakBytes) / 1.e9)
        << std::format("  seeded {} particles: {} findable {} ({} first seeded here); by origin {} primaries, {} decay daughters, {} material; {} rows held\n",
                       mSeeded.size(), seededFindable, mPopulation, firstFindable,
                       seededByOrigin[static_cast<std::size_t>(study::ParticleOrigin::kPhysicalPrimary)],
                       seededByOrigin[static_cast<std::size_t>(study::ParticleOrigin::kDecay)],
                       seededByOrigin[static_cast<std::size_t>(study::ParticleOrigin::kMaterial)], hold.size());
    std::uint64_t lost = 0U;
    for (const auto count : lostByLength) {
      lost += count;
    }
    if (mFit) {
      std::string line =
        std::format("  road fit: {} records ({} pure) in {:.3f} s; findable {} with a pure road {}, without a pure record {} ({:.2f} %), by pure road length:",
                    mPass.records, mPass.pureRecords, cost.fitSeconds, mPopulation, roadFindable, lost, study::percent(lost, roadFindable));
      for (std::size_t length = 0; length < lostByLength.size(); ++length) {
        line += lostByLength[length] == 0U ? std::string{} : std::format(" {}: {}", length, lostByLength[length]);
      }
      line += "; refused by";
      for (std::size_t why = 0; why < kFitRefusals; ++why) {
        line += std::format(" {} {}", kFitRefusalNames[why], lostBy[why]);
      }
      out << line << '\n';
    }
    if (mReal) {
      out << std::format(
        "  real run: {} tracks, {} good, {} fake; held rows {}, of them {} carry no found particle; the found particles' rows {}, {} of them "
        "free ({:.2f} %), on {} particles ({} findable {})\n",
        mReal->tracks, mReal->tracks - mReal->fakes, mReal->fakes, mReal->heldRows, mReal->wrongHeld, mReal->foundRows, mReal->foundRowsFree,
        study::percent(mReal->foundRowsFree, mReal->foundRows), mReal->foundWithFree, mReal->findableWithFree, mPopulation);
    }
    if (versus) {
      out << std::format("  against shadow pass {}: both seed {} findable {}; only this pass {} (in the shadow: {}), only the shadow {} (here: {})\n",
                         mShadowPass, bothSeeded, mPopulation, mainOnly[0] + mainOnly[1] + mainOnly[2] + mainOnly[3], mechanismLine(mainOnly),
                         shadowOnly[0] + shadowOnly[1] + shadowOnly[2] + shadowOnly[3], mechanismLine(shadowOnly));
    }
    std::string line = std::format("  coverage of the findable {} at or above:", mPopulation);
    for (std::size_t step = 0; step < mSteps.size(); ++step) {
      line += std::format(" {} {:.2f}", mSteps[step], study::percent(covered[step], findable[step]));
    }
    out << line << '\n';
    line = std::format("  at the floor, by class, covered % (this pass's gain %) of the findable:");
    for (std::size_t kind = 0; kind < kClasses; ++kind) {
      line += std::format(" {} {:.2f} ({:+.2f}) of {}", kClassNames[kind], study::percent(classCovered[kind], classFindable[kind]),
                          study::percent(classGained[kind], classFindable[kind]), classFindable[kind]);
    }
    out << line << '\n';
    // The rows left free: noise, a findable particle's (unseeded), another
    // particle's at or above the lowest step, or only particles below it.
    const double lowest = *std::ranges::min_element(mSteps);
    std::array<std::uint64_t, kFreeKinds> free{};
    for (std::size_t row = 0; row < mHeld.size(); ++row) {
      if (mHeld[row] != 0U) {
        continue;
      }
      const auto& truthRow = mTruth.rows[row];
      bool findableRow = false;
      bool above = false;
      for (std::uint8_t at = 0U; at < truthRow.count; ++at) {
        const auto particle = mIndex.at(truthRow.particles[at]);
        findableRow = findableRow || mIsFindable[particle] != 0U;
        above = above || static_cast<double>(mPt[particle]) >= lowest;
      }
      ++free[truthRow.count == 0U ? 0U : (findableRow ? 1U : (above ? 2U : 3U))];
    }
    std::uint64_t freeTotal = 0U;
    for (const auto count : free) {
      freeTotal += count;
    }
    out << std::format(
      "  rows free after the pass: {:.2f} % of all; of them noise {:.2f} %, findable {} {:.2f} %, other particles >= {} GeV {:.2f} %, "
      "only particles below {:.2f} %\n",
      study::percent(freeTotal, mHeld.size()), study::percent(free[0], freeTotal), mPopulation, study::percent(free[1], freeTotal), lowest, study::percent(free[2], freeTotal),
      study::percent(free[3], freeTotal));
    std::vector<std::pair<std::uint64_t, std::uint64_t>> patterns(unseeded.begin(), unseeded.end());
    std::ranges::sort(patterns, [](const auto& a, const auto& b) { return a.second != b.second ? a.second > b.second : a.first < b.first; });
    std::uint64_t unseededTotal = 0U;
    for (const auto& [mask, count] : patterns) {
      unseededTotal += count;
    }
    out << std::format("  unseeded at the floor: {} (over the passes so far: {}); the most common layer patterns:\n", unseededTotal, mechanismLine(unseededBy));
    for (std::size_t at = 0; at < std::min<std::size_t>(patterns.size(), kPatterns); ++at) {
      out << std::format("    {:7d} ({:5.2f} %)  {} (mask {})\n", patterns[at].second, study::percent(patterns[at].second, unseededTotal),
                         layerList(patterns[at].first), patterns[at].first);
    }
    if (json != nullptr) {
      std::string coverage;
      for (std::size_t step = 0; step < mSteps.size(); ++step) {
        coverage += std::format("{}[{},{},{}]", step == 0U ? "" : ",", mSteps[step], covered[step], findable[step]);
      }
      std::string classes;
      for (std::size_t kind = 0; kind < kClasses; ++kind) {
        classes += std::format(R"({}"{}":[{},{},{}])", kind == 0U ? "" : ",", kClassNames[kind], classCovered[kind], classGained[kind], classFindable[kind]);
      }
      std::string top;
      for (std::size_t at = 0; at < std::min<std::size_t>(patterns.size(), kPatterns); ++at) {
        top += std::format("{}[{},{}]", at == 0U ? "" : ",", patterns[at].first, patterns[at].second);
      }
      std::string origins;
      for (std::size_t origin = 0; origin < study::kParticleOrigins; ++origin) {
        origins += std::format("{}{}", origin == 0U ? "" : ",", seededByOrigin[origin]);
      }
      const auto fits = [&](std::uint64_t count) { return cost.counted ? std::to_string(count) : std::string{"null"}; };
      const auto byMechanism = [](const std::array<std::uint64_t, kMechanisms>& counts) {
        return std::format("[{},{},{},{}]", counts[0], counts[1], counts[2], counts[3]);
      };
      std::string lostLengths;
      for (std::size_t length = 0; length < lostByLength.size(); ++length) {
        if (lostByLength[length] != 0U) {
          lostLengths += std::format("{}[{},{}]", lostLengths.empty() ? "" : ",", length, lostByLength[length]);
        }
      }
      const std::string fitJson =
        mFit && !mShadow
          ? std::format(R"({{"records":{},"pureRecords":{},"seconds":{:.6f},"roadFindable":{},"lost":{},"lostByLength":[{}],"lostByClass":[{},{},{},{}],)"
                        R"("lostBy":[{},{},{},{},{},{}]}})",
                        mPass.records, mPass.pureRecords, cost.fitSeconds, roadFindable, lost, lostLengths, lostByClass[0], lostByClass[1], lostByClass[2],
                        lostByClass[3], lostBy[0], lostBy[1], lostBy[2], lostBy[3], lostBy[4], lostBy[5])
          : std::string{"null"};
      const std::string realJson =
        mReal
          ? std::format(R"({{"tracks":{},"fakes":{},"heldRows":{},"wrongHeld":{},"foundRows":{},"foundRowsFree":{},"foundWithFree":{},"findableWithFree":{}}})",
                        mReal->tracks, mReal->fakes, mReal->heldRows, mReal->wrongHeld, mReal->foundRows, mReal->foundRowsFree, mReal->foundWithFree,
                        mReal->findableWithFree)
          : std::string{"null"};
      const std::string versusJson =
        versus ? std::format(R"({{"shadowPass":{},"both":{},"mainOnly":{},"shadowOnly":{}}})", mShadowPass, bothSeeded, byMechanism(mainOnly), byMechanism(shadowOnly))
               : std::string{"null"};
      *json << std::format(
                 R"({{"pass":{},"floor":{},"layerMask":{},"startLayers":{},"minClusters":{},"remove":{},"tracklets":{},"triplets":{},"neighbours":{},)"
                 R"("tripletFits":{},"neighbourFits":{},"roads":{},"impureRoads":{},"pureRoads":{},"owners":{},"seeded":{},"seededByOrigin":[{}],)"
                 R"("seededFindable":{},"firstFindable":{},"rowsHeld":{},"seconds":{:.6f},"poolPeakBytes":{},)"
                 R"("coverage":[{}],"classes":{{{}}},"unseeded":{},"unseededBy":{},"patterns":[{}],"free":[{},{},{},{},{}],)"
                 R"("shadow":{},"kinkSigma":{},"graphKinkSigma":{},"fit":{},"versus":{},"real":{}}})",
                 mPassIndex, floor, settings.layerMask, settings.roadStartLayers, settings.minClusters, mRemove, mPass.tracklets, mPass.triplets,
                 mPass.neighbours, fits(cost.tripletFits), fits(cost.neighbourFits), mPass.roads, mPass.impureRoads, mPass.pureRoads, mPass.owners,
                 mSeeded.size(), origins, seededFindable, firstFindable, hold.size(), cost.seconds, cost.poolPeakBytes, coverage, classes, unseededTotal,
                 byMechanism(unseededBy), top, freeTotal, free[0], free[1], free[2], free[3], mShadow, settings.kinkSigma, settings.graphKinkSigma, fitJson, versusJson, realJson)
            << '\n';
    }
    if (mShadow) {
      mShadowPass = mPassIndex;
      mShadowBest = mBest;
      mShadowPure = mPureLength;
      mShadowMinClusters = settings.minClusters;
    }
    ++mPassIndex;
    return hold;
  }

  /// Per findable, one CSV line: its id, pT, class and layers, the
  /// first pass that seeded it, the first that gave it a pure road, the first
  /// shadow pass that seeded it, and its best mechanism over the passes
  /// (0 inside a road, 1 partly, 2 in none, 3 a pure road the fit refused,
  /// -1 seeded); -1 for never, -2 for found by the real run before the
  /// passes (`startAfter`); its PDG
  /// code, eta and production vertex (radius, z [cm]); with the fit, the passes whose fit refused its pure road
  /// and the last refusal (`FitRefusal`, -1 none); its mother's PDG code (0 none).
  void writeParticles(std::ostream& out) const
  {
    const auto pass = [](std::uint32_t at) { return at == kUnseeded ? std::string{"-1"} : (at == kByRun ? std::string{"-2"} : std::to_string(at)); };
    out << "id,pt,class,layers,seededBy,roadBy,shadowBy,mechanism,pdg,eta,vr,vz,refusals,refusal,motherPdg\n";
    for (const auto& particle : mFindable) {
      const auto number = particle.particle;
      const auto found = mTruth.particles.find(mIds[number]);
      const auto truth = found == mTruth.particles.end() ? study::ParticleTruth{} : found->second;
      out << std::format("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n", mIds[number], particle.pt, static_cast<int>(particle.kind), particle.layers,
                         pass(mSeededBy[number]), pass(mRoadBy[number]), pass(mShadowBy[number]),
                         mSeededBy[number] == kUnseeded ? static_cast<int>(mEver[number]) : -1, truth.pdg, truth.eta, std::hypot(truth.vx, truth.vy), truth.vz,
                         mRefusals[number], mRefusal[number], truth.motherPdg);
    }
  }

 private:
  static constexpr std::uint32_t kUnseeded = ~0U;
  static constexpr std::uint32_t kByRun = kUnseeded - 1U; ///< found by the real run before the passes (`startAfter`)
  /// A real run's iteration (`observeReal`).
  struct RealCounts {
    std::uint64_t tracks = 0U;
    std::uint64_t fakes = 0U;
    std::uint64_t heldRows = 0U;         ///< held exclusively after the iteration
    std::uint64_t wrongHeld = 0U;        ///< of them, carrying no particle found so far
    std::uint64_t foundRows = 0U;        ///< the rows of every particle found so far
    std::uint64_t foundRowsFree = 0U;    ///< of them, not held
    std::uint64_t foundWithFree = 0U;    ///< the found particles with a row not held
    std::uint64_t findableWithFree = 0U; ///< of them, findable
  };
  std::optional<RealCounts> mReal;
  static constexpr std::size_t kPatterns = 12U;
  enum class Class : std::uint8_t { kBarrel,
                                    kDiscPositive,
                                    kDiscNegative,
                                    kDiscBoth };
  /// Where an unseeded particle's rows were in a pass: inside a road with
  /// foreign rows (at least the pass's minimum of its own), only partly in a
  /// road, in no road, or -- with the fit -- in a pure road the fit refused
  /// (the codes are the CSV's; `rank` orders them).
  enum class Mechanism : std::uint8_t { kInside,
                                        kPartly,
                                        kNone,
                                        kRefused };
  static constexpr std::size_t kMechanisms = 4U;
  static constexpr std::size_t kClasses = 4U;
  static constexpr std::array<std::string_view, kClasses> kClassNames{"barrel", "discs z>0", "discs z<0", "discs both"};

  struct Findable {
    std::uint32_t particle = 0U;
    float pt = 0.F;
    std::uint64_t layers = 0U;
    Class kind = Class::kBarrel;
  };
  struct PassCounts {
    std::uint64_t tracklets = 0U;
    std::uint64_t triplets = 0U;
    std::uint64_t neighbours = 0U;
    std::uint64_t roads = 0U;
    std::uint64_t impureRoads = 0U;
    std::uint64_t pureRoads = 0U;   ///< roads some particle owns
    std::uint64_t owners = 0U;      ///< the owners of the pure roads, a road of two particles twice
    std::uint64_t records = 0U;     ///< with the fit: the accepted roads' records
    std::uint64_t pureRecords = 0U; ///< of them, those some particle owns
  };

  /// The mechanism of a particle not seeded in a pass where its longest pure
  /// road had \a pure rows (0: none) and its most rows in one road were \a best.
  [[nodiscard]] static Mechanism mechanismOf(std::uint8_t pure, std::uint8_t best, int minClusters)
  {
    if (pure > 0U) {
      return Mechanism::kRefused;
    }
    return best >= minClusters ? Mechanism::kInside : (best > 0U ? Mechanism::kPartly : Mechanism::kNone);
  }

  /// Mechanisms best first: refused, inside, partly, none.
  [[nodiscard]] static int rank(Mechanism mechanism) { return mechanism == Mechanism::kRefused ? -1 : static_cast<int>(mechanism); }

  [[nodiscard]] static std::string mechanismLine(const std::array<std::uint64_t, kMechanisms>& counts)
  {
    return std::format("{} inside an impure road, {} partly in one, {} in none, {} in a pure road the fit refused", counts[0], counts[1], counts[2], counts[3]);
  }

  /// \a particle seeded by this pass: a shadow pass's own record, else
  /// the study's, and held after the pass.
  void seed(std::uint32_t particle)
  {
    auto& in = mShadow ? mShadowIn : mSeededIn;
    if (in[particle] == mPassIndex) {
      return;
    }
    in[particle] = mPassIndex;
    mSeeded.push_back(particle);
    auto& by = mShadow ? mShadowBy : mSeededBy;
    by[particle] = by[particle] == kUnseeded ? mPassIndex : by[particle];
  }

  /// Call \a body with the number of every particle all of \a held carries
  /// (each of the first row's particles the others carry too: the one truth
  /// rule); whether there was one.
  template <typename Body>
  bool forEachOwner(std::span<const std::uint32_t> held, Body&& body)
  {
    return forEachCommonParticle(mTruth.rows, held, [&](std::uint64_t id) { body(mIndex.at(id)); });
  }

  /// Each particle's most rows in one road, over the pass: \a rows' labels
  /// counted.
  void noteOwnRows(std::span<const std::uint32_t> rows)
  {
    std::array<std::uint32_t, maxChainLayers<Det> * study::TruthRow::kMaxParticles> labels{};
    std::size_t count = 0U;
    for (const auto row : rows) {
      for (auto at = mRowOffsets[row]; at < mRowOffsets[row + 1U]; ++at) {
        labels[count++] = mRowParticles[at];
      }
    }
    std::sort(labels.begin(), labels.begin() + static_cast<std::ptrdiff_t>(count));
    for (std::size_t at = 0; at < count;) {
      std::size_t end = at + 1U;
      while (end < count && labels[end] == labels[at]) {
        ++end;
      }
      mBest[labels[at]] = std::max(mBest[labels[at]], static_cast<std::uint8_t>(end - at));
      at = end;
    }
  }

  /// A layer mask as its layers, the barrel's before a bar and the discs'.
  [[nodiscard]] std::string layerList(std::uint64_t mask) const
  {
    std::string barrel;
    std::string discs;
    for (int layer = 0; layer < Detector<Det>::nLayers; ++layer) {
      if (((mask >> layer) & 1U) == 0U) {
        continue;
      }
      auto& part = ((mComparison.barrelLayers >> layer) & 1U) != 0U ? barrel : discs;
      part += std::format("{}{}", part.empty() ? "" : " ", layer);
    }
    return discs.empty() ? barrel : std::format("{} | {}", barrel, discs);
  }

  const study::TrackingTruth& mTruth;
  const study::ComparisonSettings& mComparison;
  std::vector<double> mSteps;
  bool mRemove = true;
  bool mFit = false;
  std::string_view mPopulation;          ///< the scored particles' name
  bool mShadow = false;                  ///< the pass is a shadow pass
  std::uint32_t mShadowPass = kUnseeded; ///< the last shadow pass
  int mShadowMinClusters = 0;
  std::vector<std::uint8_t> mBest;                                         ///< per number, its most rows in one road of the pass
  std::vector<std::uint8_t> mShadowBest;                                   ///< `mBest` of the last shadow pass
  std::vector<std::uint8_t> mPureLength;                                   ///< per number, its longest pure road in the pass
  std::vector<std::uint8_t> mShadowPure;                                   ///< `mPureLength` of the last shadow pass
  std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> mPureRows; ///< with the fit: per findable, that road's rows
  std::vector<std::uint32_t> mShadowIn;                                    ///< per number, the last shadow pass that seeded it
  std::vector<std::uint32_t> mShadowBy;                                    ///< per number, the first shadow pass that seeded it
  std::vector<std::uint32_t> mRoadIn;                                      ///< per number, the last pass a pure road of which it owned
  std::vector<std::uint32_t> mRoadBy;                                      ///< per number, the first such pass
  std::vector<Mechanism> mEver;                                            ///< per number, its best mechanism over the passes that did not seed it
  std::vector<std::int8_t> mRefusal;                                       ///< with the fit: per number, the last refusal of its pure road, -1 none
  std::vector<std::uint8_t> mRefusals;                                     ///< with the fit: per number, the passes that refused it
  std::vector<std::size_t> mRowOffsets;                                    ///< per row, its first number in `mRowParticles`
  std::vector<std::uint32_t> mRowParticles;                                ///< the rows' particles, as numbers
  std::unordered_map<std::uint64_t, std::uint32_t> mIndex;                 ///< particle id -> number
  std::vector<std::uint64_t> mIds;                                         ///< number -> particle id
  std::vector<std::size_t> mParticleOffsets;                               ///< per number, its first row in `mParticleRows`
  std::vector<std::uint32_t> mParticleRows;
  std::vector<std::uint32_t> mSeededBy; ///< per number, the pass that first seeded it
  std::vector<std::uint32_t> mSeededIn; ///< per number, the last pass that seeded it
  std::uint64_t mPositiveDiscs = 0U;    ///< the disc layers at z > 0, a bit per layer
  std::uint64_t mNegativeDiscs = 0U;
  static constexpr std::size_t kFreeKinds = 4U;
  std::vector<std::uint8_t> mHeld; ///< per row, held by a seeded particle
  std::vector<float> mPt;          ///< per number, the particle's pT [GeV]; 0 without its truth
  std::vector<std::uint8_t> mIsFindable;
  std::vector<Findable> mFindable;
  std::vector<std::uint32_t> mSeeded; ///< the particles the pass seeded
  std::uint32_t mPassIndex = 0U;
  PassCounts mPass;
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_SEEDINGSTRATEGY_H_

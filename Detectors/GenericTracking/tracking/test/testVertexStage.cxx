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

/// \file testVertexStage.cxx
/// \brief The vertexer on synthetic triplets: separated and close vertices,
///        frames, the contributor floor, the cuts, a fake tracklet shared by
///        many triplets, and a list that does not depend on the threads or on
///        the triplets' order.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK VertexStage
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/Parallel.h"
#include "GenericTracking/Tracking/VertexStage.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <bit>
#include <numbers>
#include <random>
#include <vector>

namespace o2::gtrk
{

namespace
{

/// The ITS's readout frame [BC].
constexpr std::int64_t kFrameBC = 198;
constexpr double kSigmaZ = 5.e-4; // [cm]
constexpr std::array<double, 3> kRadii{2.34, 3.15, 3.93};

/// The vertexer's candidate errors from the points alone: the r-z line's.
struct LineErrors {
  [[nodiscard]] Scalar operator()(const PairPoint& /*first*/, const VertexCandidateFit& fit) const { return fit.varZ0Line; }
};

/// The columns `pairPointAt` and the vertexer read, as plain vectors.
struct Measurements {
  std::vector<Scalar> x, y, z, r, varZ, varR, varRPhi;
  std::vector<std::uint32_t> timePS, widthPS;
};
struct Tracklets {
  std::vector<std::uint32_t> source, target;
};
struct Triplets {
  std::vector<std::uint32_t> first, second;
};

/// A timeframe of triplets, built track by track.
class Event
{
 public:
  explicit Event(std::uint64_t seed = 1U) : mRandom(seed) {}

  std::uint32_t point(double x, double y, double z, const TimeEstPS& box)
  {
    const auto row = static_cast<std::uint32_t>(m.x.size());
    m.x.push_back(static_cast<Scalar>(x));
    m.y.push_back(static_cast<Scalar>(y));
    m.z.push_back(static_cast<Scalar>(z));
    m.r.push_back(static_cast<Scalar>(std::hypot(x, y)));
    m.varZ.push_back(static_cast<Scalar>(kSigmaZ * kSigmaZ));
    m.varR.push_back(0.F);
    m.varRPhi.push_back(static_cast<Scalar>(kSigmaZ * kSigmaZ));
    m.timePS.push_back(box.start);
    m.widthPS.push_back(box.width);
    return row;
  }

  std::uint32_t tracklet(std::uint32_t source, std::uint32_t target)
  {
    t.source.push_back(source);
    t.target.push_back(target);
    return static_cast<std::uint32_t>(t.source.size() - 1U);
  }

  void triplet(std::uint32_t first, std::uint32_t second)
  {
    p.first.push_back(first);
    p.second.push_back(second);
  }

  /// A track from (0, 0, \a zVertex) at azimuth \a phi, with cot theta
  /// \a cotTheta and transverse radius of curvature \a radius [cm] (bending
  /// by \a charge), in frame \a frame; z measured with kSigmaZ. \a offset
  /// shifts its points sideways [cm], away from the beam axis.
  void track(double zVertex, double phi, double cotTheta, double radius, int charge, std::int64_t frame, double offset = 0., double kink = 0.)
  {
    std::normal_distribution<double> noise{0., kSigmaZ};
    std::array<std::uint32_t, 3> rows{};
    for (std::size_t layer = 0; layer < kRadii.size(); ++layer) {
      const double r = kRadii[layer];
      const double half = std::asin(r / (2. * radius));
      const double angle = phi + (charge * half);
      const double s = 2. * radius * half;
      const double z = zVertex + (s * cotTheta) + (layer == 2U ? kink : 0.) + noise(mRandom);
      rows[layer] = point((r * std::cos(angle)) - (offset * std::sin(phi)), (r * std::sin(angle)) + (offset * std::cos(phi)), z, frameTime(frame));
    }
    triplet(tracklet(rows[0], rows[1]), tracklet(rows[1], rows[2]));
  }

  /// \a count tracks from \a zVertex in \a frame, spread in azimuth and dip.
  void vertex(double zVertex, std::size_t count, std::int64_t frame)
  {
    std::uniform_real_distribution<double> phi{0., 2. * std::numbers::pi};
    std::uniform_real_distribution<double> cot{-1.2, 1.2};
    std::uniform_real_distribution<double> radius{300., 1500.};
    for (std::size_t i = 0; i < count; ++i) {
      track(zVertex, phi(mRandom), cot(mRandom), radius(mRandom), (i % 2U) == 0U ? 1 : -1, frame);
    }
  }

  /// \a count triplets of unrelated points in \a frame.
  void noise(std::size_t count, std::int64_t frame)
  {
    std::uniform_real_distribution<double> phi{0., 2. * std::numbers::pi};
    std::uniform_real_distribution<double> z{-15., 15.};
    for (std::size_t i = 0; i < count; ++i) {
      std::array<std::uint32_t, 3> rows{};
      for (std::size_t layer = 0; layer < kRadii.size(); ++layer) {
        const double angle = phi(mRandom);
        rows[layer] = point(kRadii[layer] * std::cos(angle), kRadii[layer] * std::sin(angle), z(mRandom), frameTime(frame));
      }
      triplet(tracklet(rows[0], rows[1]), tracklet(rows[1], rows[2]));
    }
  }

  /// Move every point by (\a dx, \a dy): the tracks then come from a beam
  /// point off the axis the columns are about.
  void shift(double dx, double dy)
  {
    for (std::size_t i = 0; i < m.x.size(); ++i) {
      const double x = static_cast<double>(m.x[i]) + dx;
      const double y = static_cast<double>(m.y[i]) + dy;
      m.x[i] = static_cast<Scalar>(x);
      m.y[i] = static_cast<Scalar>(y);
      m.r[i] = static_cast<Scalar>(std::hypot(x, y));
    }
  }

  static TimeEstPS frameTime(std::int64_t frame) { return frameBox(frame * kFrameBC, (frame + 1) * kFrameBC, 0); }

  [[nodiscard]] VertexerResult find(const VertexerSettings& settings, detail::TaskArena* arena = nullptr) const { return Vertexer{settings, VertexConstraint{}}.find(m, t, p, p.first.size(), arena, true, LineErrors{}); }

  Measurements m;
  Tracklets t;
  Triplets p;

 private:
  std::mt19937_64 mRandom;
};

bool sameVertices(const VertexerResult& a, const VertexerResult& b)
{
  if (a.vertices.size() != b.vertices.size() || a.contributors != b.contributors) {
    return false;
  }
  // Bit for bit: every field widened to double exactly, then its bits.
  const auto same = [](auto x, auto y) { return std::bit_cast<std::uint64_t>(static_cast<double>(x)) == std::bit_cast<std::uint64_t>(static_cast<double>(y)); };
  for (std::size_t v = 0; v < a.vertices.size(); ++v) {
    const auto& x = a.vertices[v];
    const auto& y = b.vertices[v];
    if (x.timePS != y.timePS || x.widthPS != y.widthPS || !same(x.z, y.z) || !same(x.varZ, y.varZ) || !same(x.varTransverse, y.varTransverse)) {
      return false;
    }
  }
  return true;
}

} // namespace

BOOST_AUTO_TEST_CASE(EmptyTimeframe)
{
  const Event event;
  const auto result = event.find({});
  BOOST_TEST(result.vertices.empty());
  BOOST_TEST(result.candidates == 0U);
}

BOOST_AUTO_TEST_CASE(TwoVerticesFarApartInOneFrame)
{
  Event event;
  event.vertex(-3.2, 60U, 7);
  event.vertex(1.8, 40U, 7);
  const auto result = event.find({});
  BOOST_TEST_REQUIRE(result.vertices.size() == 2U);
  BOOST_TEST(std::abs(result.vertices[0].z - -3.2F) < 1.e-3F);
  BOOST_TEST(std::abs(result.vertices[1].z - 1.8F) < 1.e-3F);
  // The chi2 cut takes the tail of the z noise, a few in a hundred; every
  // candidate left belongs to one of the two.
  BOOST_TEST(result.contributors[0] + result.contributors[1] == result.candidates);
  BOOST_TEST(result.contributors[0] >= 54U);
  BOOST_TEST(result.contributors[1] >= 36U);
  // The vertex's box spans its triplets', here all in frame 7.
  BOOST_TEST(result.vertices[0].timePS == Event::frameTime(7).start);
  BOOST_TEST(result.vertices[0].widthPS == Event::frameTime(7).width);
  // Its z error is the mean's from the core's spread, no floor:
  // a few micrometres for sixty candidates.
  BOOST_TEST(result.vertices[0].varZ > 0.F);
  BOOST_TEST(std::sqrt(result.vertices[0].varZ) < 1.e-3F);
  BOOST_TEST(result.rounds == 3U);
}

BOOST_AUTO_TEST_CASE(CloseVerticesMergeWithinTheClaim)
{
  Event event;
  event.vertex(0.50, 60U, 3);
  event.vertex(0.60, 30U, 3);
  VertexerSettings settings;
  // 1 mm apart: the 2 mm claim takes both.
  const auto merged = event.find(settings);
  BOOST_TEST_REQUIRE(merged.vertices.size() == 1U);
  BOOST_TEST(merged.contributors[0] == merged.candidates);
  // A 500 um claim separates them.
  settings.claimUM = 500;
  const auto separate = event.find(settings);
  BOOST_TEST_REQUIRE(separate.vertices.size() == 2U);
  BOOST_TEST(std::abs(separate.vertices[0].z - 0.50F) < 1.e-3F);
  BOOST_TEST(std::abs(separate.vertices[1].z - 0.60F) < 1.e-3F);
}

BOOST_AUTO_TEST_CASE(FramesStaySeparateAndInTimeOrder)
{
  Event event;
  event.vertex(2.0, 20U, 9);
  event.vertex(2.0, 20U, 5);
  const auto result = event.find({});
  BOOST_TEST_REQUIRE(result.vertices.size() == 2U);
  BOOST_TEST(result.vertices[0].timePS == Event::frameTime(5).start);
  BOOST_TEST(result.vertices[1].timePS == Event::frameTime(9).start);
}

BOOST_AUTO_TEST_CASE(TheContributorFloor)
{
  const std::uint32_t floor = VertexerSettings{}.minContributors;
  Event below;
  below.vertex(1.0, floor - 1U, 1);
  BOOST_TEST(below.find({}).vertices.empty());
  Event at;
  at.vertex(1.0, floor, 1);
  BOOST_TEST(at.find({}).vertices.size() == 1U);
}

/// A frame whose rounds find no vertex runs them again at
/// aloneMinContributors; a frame with a vertex does not, so a small
/// collision beside a large one stays below the floor.
BOOST_AUTO_TEST_CASE(AFrameWithoutAVertexRerunsAtTheLowerFloor)
{
  Event event;
  event.vertex(2.0, 20U, 3);
  event.vertex(-6.0, 2U, 3);
  event.vertex(1.0, 2U, 6);
  BOOST_TEST(event.find({}).vertices.size() == 1U);
  const auto result = event.find({.aloneMinContributors = 2U});
  BOOST_TEST_REQUIRE(result.vertices.size() == 2U);
  BOOST_TEST(result.vertices[0].timePS == Event::frameTime(3).start);
  BOOST_TEST(std::abs(result.vertices[0].z - 2.F) < 0.05F);
  BOOST_TEST(result.vertices[1].timePS == Event::frameTime(6).start);
  BOOST_TEST(std::abs(result.vertices[1].z - 1.F) < 0.05F);
  BOOST_TEST(result.contributors[1] == 2U);
}

/// A dense collision's measurements recombined into triplets that point
/// elsewhere make no vertex: the dense vertex owns the measurements of its
/// frame's triplets within its claim, and the free candidates with one are
/// dropped. A small collision beside it, on measurements of
/// its own, keeps its vertex.
BOOST_AUTO_TEST_CASE(AVertexOwnsItsMeasurements)
{
  Event event;
  event.vertex(-2.0, 80U, 4);
  // Twelve fakes at z = 3 cm: each a dense track's innermost measurement
  // with two new points on its transverse path, on a line to z = 3 cm.
  const auto denseTriplets = event.p.first.size();
  for (std::size_t i = 0; i < 12U; ++i) {
    const auto a = event.t.source[event.p.first[i]];
    const auto b = event.t.target[event.p.first[i]];
    const auto c = event.t.target[event.p.second[i]];
    const double slope = (static_cast<double>(event.m.z[a]) - 3.) / static_cast<double>(event.m.r[a]);
    const auto box = Event::frameTime(4);
    const auto fakeB = event.point(event.m.x[b], event.m.y[b], 3. + (slope * static_cast<double>(event.m.r[b])), box);
    const auto fakeC = event.point(event.m.x[c], event.m.y[c], 3. + (slope * static_cast<double>(event.m.r[c])) + 2.e-4, box);
    event.triplet(event.tracklet(a, fakeB), event.tracklet(fakeB, fakeC));
  }
  event.vertex(6.0, 8U, 4);
  const auto result = event.find({});
  BOOST_TEST_REQUIRE(result.vertices.size() == 2U);
  BOOST_TEST(std::abs(result.vertices[0].z - -2.0F) < 1.e-3F);
  BOOST_TEST(std::abs(result.vertices[1].z - 6.0F) < 1.e-3F);
  BOOST_TEST(result.contributors[1] == 8U);
  // Every fake is free: kept by the deduplication or not, none claimed.
  for (std::size_t t = denseTriplets; t < denseTriplets + 12U; ++t) {
    BOOST_TEST(result.tripletVertex[t] < 0);
  }
}

BOOST_AUTO_TEST_CASE(TheCuts)
{
  Event event;
  event.track(1.0, 0.3, 0.5, 800., 1, 2, 0.1);      // 1 mm from the axis
  event.track(1.0, 1.3, 0.5, 800., 1, 2, 0., 0.05); // a kink of 500 um in z
  event.track(1.0, 2.3, 0.5, 800., 1, 2);
  const auto result = event.find({});
  BOOST_TEST(result.accepted == 1U);
  BOOST_TEST(result.candidates == 1U);
  BOOST_TEST(result.tripletVertex[0] == VertexerResult::kNone);
  BOOST_TEST(result.tripletVertex[1] == VertexerResult::kNone);
  BOOST_TEST(result.tripletVertex[2] == VertexerResult::kFree);
}

BOOST_AUTO_TEST_CASE(OneTrackletSharedByManyTripletsCountsOnce)
{
  // One tracklet 0 -> 1 pointing at z = 4, with 12 third points along its
  // line: 12 triplets, nearly one z -- one candidate once they are made
  // cluster-disjoint, so no vertex.
  Event event;
  const auto a = event.point(2.34, 0., 4. + (2.34 * 0.5), Event::frameTime(4));
  const auto b = event.point(3.15, 0., 4. + (3.15 * 0.5), Event::frameTime(4));
  const auto shared = event.tracklet(a, b);
  for (int i = 0; i < 12; ++i) {
    const double angle = 1.e-4 * i;
    const auto c = event.point(3.93 * std::cos(angle), 3.93 * std::sin(angle), 4. + (3.93 * 0.5) + (2.e-4 * i), Event::frameTime(4));
    event.triplet(shared, event.tracklet(b, c));
  }
  VertexerSettings settings;
  settings.maxChi2RZ = 0.;
  settings.maxDcaCM = 0.;
  const auto result = event.find(settings);
  BOOST_TEST(result.accepted == 12U);
  BOOST_TEST(result.candidates == 1U);
  BOOST_TEST(result.vertices.empty());
}

/// A particle with two clusters on a layer, where sensors overlap, makes two
/// triplets; it counts once, through the better-fitting one.
BOOST_AUTO_TEST_CASE(AParticleWithTwoTripletsCountsOnce)
{
  Event event;
  constexpr int kTracks = 8;
  for (int i = 0; i < kTracks; ++i) {
    const double phi = 0.7 * i;
    const auto at = [&](double r, double dz) { return event.point(r * std::cos(phi), r * std::sin(phi), 1. + (0.5 * r) + dz, Event::frameTime(2)); };
    const auto a = at(2.34, 0.);
    const auto b = at(3.15, 0.);
    const auto bOverlap = at(3.16, 4.e-4);
    const auto c = at(3.93, 0.);
    event.triplet(event.tracklet(a, b), event.tracklet(b, c));
    event.triplet(event.tracklet(a, bOverlap), event.tracklet(bOverlap, c));
  }
  const auto result = event.find({});
  BOOST_TEST(result.accepted == 2U * kTracks);
  BOOST_TEST(result.candidates == static_cast<std::size_t>(kTracks));
  BOOST_REQUIRE(result.vertices.size() == 1U);
  BOOST_TEST(result.contributors[0] == static_cast<std::uint32_t>(kTracks));
  BOOST_TEST(std::abs(result.vertices[0].z - 1.F) < 1.e-4F);
}

BOOST_AUTO_TEST_CASE(NoiseMakesNoVertex)
{
  Event event{7U};
  event.noise(3000U, 1);
  event.vertex(-6.0, 200U, 1);
  const auto result = event.find({});
  BOOST_TEST_REQUIRE(result.vertices.size() == 1U);
  BOOST_TEST(std::abs(result.vertices[0].z - -6.0F) < 1.e-3F);
}

BOOST_AUTO_TEST_CASE(TheListDoesNotDependOnThreadsOrOrder)
{
  Event event{11U};
  for (std::int64_t frame = 1; frame < 40; frame += 2) {
    event.vertex(-10. + (0.5 * static_cast<double>(frame)), 30U + static_cast<std::size_t>(frame), frame);
    event.vertex(4. - (0.3 * static_cast<double>(frame)), 12U, frame);
    event.noise(400U, frame);
  }
  const auto serial = event.find({});
  BOOST_TEST(serial.vertices.size() >= 40U);
  for (const int threads : {1, 4, 16}) {
    detail::TaskArena arena{threads};
    const auto parallel = event.find({}, &arena);
    BOOST_TEST(sameVertices(serial, parallel));
    BOOST_TEST(serial.tripletVertex == parallel.tripletVertex);
  }
  // The triplets in reverse: the same vertices, the claims carried along.
  Event reversed = event;
  std::ranges::reverse(reversed.p.first);
  std::ranges::reverse(reversed.p.second);
  const auto backwards = reversed.find({});
  BOOST_TEST(sameVertices(serial, backwards));
  auto claims = backwards.tripletVertex;
  std::ranges::reverse(claims);
  BOOST_TEST(claims == serial.tripletVertex);
}

} // namespace o2::gtrk

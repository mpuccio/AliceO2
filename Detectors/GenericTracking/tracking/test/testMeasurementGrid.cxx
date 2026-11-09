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

/// \file testMeasurementGrid.cxx
/// \brief Grid construction, ordering and candidate lookup contracts.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK MeasurementGrid
#include <boost/test/unit_test.hpp>
#include <numbers>

#include "GenericTracking/Base/Geometry.h"
#include "ToyBarrel4Disc8.h"
#include "ToyBarrel3.h"
#include "GenericTracking/Tracking/GridSizing.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"
#include "GenericTracking/Tracking/MeasurementGridBuilder.h"

#include "GridTestSupport.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <span>
#include <vector>

namespace o2::gtrk
{

namespace
{

using Toy = detectors::ToyBarrel3;
using Mixed = detectors::ToyBarrel4Disc8;
using GeometryGrid = GridGeometry<Toy>;
using Grid = MeasurementGrid<Toy>;

/// Every layer \a phi x \a second cells, whatever the windows.
GridParams fixedCellParams(std::int64_t blockWidthPS, std::int32_t phi, std::int32_t second)
{
  auto params = test::withToyFrames(GridParams{});
  params.blockWidth = blockWidthPS;
  for (std::size_t layer = 0; layer < GridParams::kMaxLayers; ++layer) {
    params.phiBins[layer] = phi;
    params.secondBins[layer] = second;
  }
  return params;
}

GridParams oneCellParams(std::int64_t blockWidthPS) { return fixedCellParams(blockWidthPS, 1, 1); }

/// The windows of one iteration with the default window physics.
auto defaultWindows() { return test::oneIterationWindows<Toy>(SearchWindowParams{}, VertexConstraint{}); }

/// Grid one row per point, all in the box [0, 1 ps), and return the positions in
/// \a points of the rows a query on \a target visits.
template <DetectorLayout Det, std::size_t N>
std::vector<std::uint32_t> visitedPoints(const GridGeometry<Det>& geometry, LayerID target, const std::array<GridPoint, N>& points, const CandidateWindow& window)
{
  std::vector<std::uint32_t> offsets(static_cast<std::size_t>(geometry.bucketCount()) + 1U);
  MeasurementStore columns;
  resizeMeasurementColumns(columns, N);
  test::buildGridRows(
    N, offsets, [&](std::size_t i) { return geometry.bucketAt(points[i].layer, 0U, points[i].phi, points[i].second); },
    [&](std::size_t i, std::uint32_t row) {
      columns.widthPS[row] = 1U;
      columns.clusterIndex[row] = static_cast<std::uint32_t>(i);
    });

  const MeasurementGrid<Det> grid{measurementView(columns), geometry, std::span<const std::uint32_t>{offsets}};
  std::vector<std::uint32_t> visited;
  grid.visitCandidates(target, window, geometry.allBlocks(), [&](std::uint32_t row) { visited.push_back(columns.clusterIndex[row]); });
  std::ranges::sort(visited);
  return visited;
}

} // namespace

BOOST_AUTO_TEST_CASE(VertexConstraintMovesAndWidensThePairWindow)
{
  SearchWindowParams params;
  params.minPtGeV = 1.;
  params.materialX2X0 = 0.;
  params.positionSigma = 0.;
  params.nSigma = 1.;

  const VertexConstraint pointVertex{.z = 1.5, .sigmaZ = 0.};
  const auto pointWindow = layerPairWindow<Toy>(0, 1, params, pointVertex);

  BOOST_TEST(pointWindow.secondScale == 2.F, boost::test_tools::tolerance(1e-6F));
  BOOST_TEST(pointWindow.secondOffset == -1.5F, boost::test_tools::tolerance(1e-6F));
  // With every Gaussian term zero, the measured-radius spread and the
  // difference between radial and arc-length extrapolation remain.
  const auto expectedHalfWidth = conservativeScalarBound(secondCentreSpread<Toy>(0, 1, pointVertex.z) + helixSecondSpread<Toy>(0, 1, params, pointVertex));
  BOOST_TEST(pointWindow.secondHalfWidth == expectedHalfWidth);
  BOOST_TEST(secondCentreSpread<Toy>(0, 1, pointVertex.z) > 0.);

  const auto candidate = candidateWindow(pointWindow, test::toyTimeBox(100.), (2 * std::numbers::pi_v<Scalar>)-0.1F, 2.F, 2.F, 41'000U);
  BOOST_TEST(candidate.phiCentre == (2 * std::numbers::pi_v<Scalar>)-0.1F, boost::test_tools::tolerance(1e-6F));
  BOOST_TEST(candidate.secondCentre == 2.5F, boost::test_tools::tolerance(1e-6F));

  const VertexConstraint broadVertex{.z = 1.5, .sigmaZ = 1.};
  const auto broadWindow = layerPairWindow<Toy>(0, 1, params, broadVertex);
  BOOST_TEST(broadWindow.secondHalfWidth > pointWindow.secondHalfWidth);

  constexpr LayerID kFirstPositiveDisc = 8;
  constexpr LayerID kSecondPositiveDisc = 9;
  const auto centralDiscWindow = layerPairWindow<Mixed>(kFirstPositiveDisc, kSecondPositiveDisc, params, VertexConstraint{.z = 0., .sigmaZ = 0.});
  const auto displacedDiscWindow = layerPairWindow<Mixed>(kFirstPositiveDisc, kSecondPositiveDisc, params, VertexConstraint{.z = 10., .sigmaZ = 0.});
  // The lever arm between the discs' measured positions (26 / 22 cm nominal,
  // moved by their rings' offsets), from the vertex.
  const double first = discZ<Mixed>(kFirstPositiveDisc);
  const double second = discZ<Mixed>(kSecondPositiveDisc);
  BOOST_TEST(centralDiscWindow.secondScale == static_cast<Scalar>(second / first), boost::test_tools::tolerance(1e-6F));
  BOOST_TEST(displacedDiscWindow.secondScale == static_cast<Scalar>((second - 10.) / (first - 10.)), boost::test_tools::tolerance(1e-6F));
  BOOST_TEST(displacedDiscWindow.secondScale > centralDiscWindow.secondScale);
}

/// A barrel <-> disc window is the hull of the lines from the ends of the
/// vertex range through the seed, evaluated per seed.
BOOST_AUTO_TEST_CASE(CrossingWindowFollowsTheSeed)
{
  SearchWindowParams params;
  params.minPtGeV = 1.;
  params.materialX2X0 = 0.;
  params.positionSigma = 0.;
  params.nSigma = 3.;
  // A vertex range of +-3 cm.
  const VertexConstraint beamSpot{.z = 0., .sigmaZ = 1.};
  const auto seed = test::toyTimeBox(0.);
  constexpr LayerID kInnermostBarrel = 0;
  constexpr LayerID kOuterBarrel = 3;
  constexpr LayerID kFirstPositiveDisc = 8;

  // Barrel 0 -> the +22 cm disc: r_j = r_i (Z - z_v) / (z_i - z_v).
  const auto toDisc = layerPairWindow<Mixed>(kInnermostBarrel, kFirstPositiveDisc, params, beamSpot);
  BOOST_TEST((toDisc.kind == PairKind::kBarrelToDisc));
  const double zLow = zMin<Mixed>(kFirstPositiveDisc);
  const double zHigh = zMax<Mixed>(kFirstPositiveDisc);
  // With the Gaussian terms zero, the padding is the helix's reach alone.
  const double helix = helixSecondSpread<Mixed>(kInnermostBarrel, kFirstPositiveDisc, params, beamSpot);
  BOOST_TEST(helix > 0.);
  const auto far = candidateWindow(toDisc, seed, 0.F, 3.F, 10.F, seed.width);
  const double expectedLow = (3. * (zLow + 3.) / (10. + 3.)) - helix;
  const double expectedHigh = (3. * (zHigh - 3.) / (10. - 3.)) + helix;
  BOOST_TEST(static_cast<double>(far.secondCentre - far.secondHalfWidth) == expectedLow, boost::test_tools::tolerance(1e-5));
  BOOST_TEST(static_cast<double>(far.secondCentre + far.secondHalfWidth) == expectedHigh, boost::test_tools::tolerance(1e-5));
  // a part of the disc, not all of it
  BOOST_TEST(static_cast<double>(far.secondHalfWidth) < 0.1 * (secondMax<Mixed>(kFirstPositiveDisc) - secondMin<Mixed>(kFirstPositiveDisc)));

  // A seed level with the vertex range sends tracks at any radius: the window
  // opens to the disc's outer edge.
  const auto near = candidateWindow(toDisc, seed, 0.F, 3.F, 2.9F, seed.width);
  BOOST_TEST(static_cast<double>(near.secondCentre + near.secondHalfWidth) == secondMax<Mixed>(kFirstPositiveDisc) + helix, boost::test_tools::tolerance(1e-5));
  BOOST_TEST(static_cast<double>(near.secondCentre - near.secondHalfWidth) == (3. * (zLow + 3.) / (static_cast<double>(2.9F) + 3.)) - helix, boost::test_tools::tolerance(1e-5));

  // The +22 cm disc -> barrel layer 3: z_j = z_v + (z_i - z_v) R / r_i.
  const auto toBarrel = layerPairWindow<Mixed>(kFirstPositiveDisc, kOuterBarrel, params, beamSpot);
  BOOST_TEST((toBarrel.kind == PairKind::kDiscToBarrel));
  const auto back = candidateWindow(toBarrel, seed, 0.F, 14.F, 22.F, seed.width);
  const double rLow = radiusMin<Mixed>(kOuterBarrel);
  const double rHigh = radiusMax<Mixed>(kOuterBarrel);
  const double backHelix = helixSecondSpread<Mixed>(kFirstPositiveDisc, kOuterBarrel, params, beamSpot);
  BOOST_TEST(static_cast<double>(back.secondCentre - back.secondHalfWidth) == 3. + (19. * rLow / 14.) - backHelix, boost::test_tools::tolerance(1e-5));
  BOOST_TEST(static_cast<double>(back.secondCentre + back.secondHalfWidth) == -3. + (25. * rHigh / 14.) + backHelix, boost::test_tools::tolerance(1e-5));
}

BOOST_AUTO_TEST_CASE(GridPointUsesThePlacedSensor)
{
  const auto geometry = toyGeometry<Toy>();
  constexpr SensorID kSensor = 5;
  constexpr double kLoc0 = 0.2;
  constexpr double kLoc1 = -0.4;

  const auto point = gridPoint<Toy>(geometry, kSensor, kLoc0, kLoc1);
  const double local[3] = {kLoc0, kLoc1, 0.};
  double global[3] = {};
  geometry.toGlobal(kSensor, local, global);
  double expectedPhi = std::atan2(global[1], global[0]);
  if (expectedPhi < 0.) {
    expectedPhi += 2 * std::numbers::pi_v<double>;
  }

  BOOST_TEST(point.layer == geometry.sensor(kSensor).layer);
  BOOST_TEST(point.phi == static_cast<Scalar>(expectedPhi), boost::test_tools::tolerance(1e-6F));
  BOOST_TEST(point.second == static_cast<Scalar>(global[2]), boost::test_tools::tolerance(1e-6F));

  // The shipped position is the double transform narrowed once, and the grid
  // coordinates come from those floats, exactly as the tracker derives them.
  BOOST_TEST(std::bit_cast<std::uint32_t>(point.x) == std::bit_cast<std::uint32_t>(static_cast<Scalar>(global[0])));
  BOOST_TEST(std::bit_cast<std::uint32_t>(point.y) == std::bit_cast<std::uint32_t>(static_cast<Scalar>(global[1])));
  BOOST_TEST(std::bit_cast<std::uint32_t>(point.z) == std::bit_cast<std::uint32_t>(static_cast<Scalar>(global[2])));
  BOOST_TEST(std::bit_cast<std::uint32_t>(point.phi) == std::bit_cast<std::uint32_t>(azimuth(point.x, point.y)));
  BOOST_TEST(std::bit_cast<std::uint32_t>(point.second) == std::bit_cast<std::uint32_t>(secondCoordinate<Toy>(point.layer, point.r, point.z)));
  BOOST_TEST(std::bit_cast<std::uint32_t>(point.r) == std::bit_cast<std::uint32_t>(o2::gpu::CAMath::Hypot(point.x, point.y)));
}

/// A bin comes from a range check, then the conversion: a NaN or an infinite
/// window edge must land in an edge cell on every platform, not wherever an
/// out-of-range float-to-int conversion happens to go (INT_MIN on x86).
BOOST_AUTO_TEST_CASE(BinsClampBeforeConverting)
{
  constexpr std::int32_t kCount = 10;
  BOOST_TEST(clampedBin(std::numeric_limits<float>::quiet_NaN(), kCount) == 0);
  BOOST_TEST(clampedBin(-std::numeric_limits<float>::infinity(), kCount) == 0);
  BOOST_TEST(clampedBin(-0.5F, kCount) == 0);
  BOOST_TEST(clampedBin(0.F, kCount) == 0);
  BOOST_TEST(clampedBin(3.7F, kCount) == 3);
  BOOST_TEST(clampedBin(9.999F, kCount) == 9);
  BOOST_TEST(clampedBin(10.F, kCount) == 9);
  BOOST_TEST(clampedBin(3.e9F, kCount) == 9);
  BOOST_TEST(clampedBin(std::numeric_limits<float>::infinity(), kCount) == 9);
  BOOST_TEST(clampedBin(std::numeric_limits<double>::quiet_NaN(), kCount) == 0);
  BOOST_TEST(clampedBin(1.e300, kCount) == 9);
}

BOOST_AUTO_TEST_CASE(GridGeometryMapsBlocksBucketsAndSlices)
{
  auto params = test::withToyFrames(GridParams{});
  params.blockWidth = 10'000;
  const auto geometry = test::cappedGrid<Toy>(25., params, SearchWindowParams{}, VertexConstraint{}, 8, 4);

  BOOST_TEST(geometry.blockCount() == 3);
  BOOST_TEST(geometry.blockOfPS(-1'000) == 0);
  BOOST_TEST(geometry.blockOfPS(0) == 0);
  BOOST_TEST(geometry.blockOfPS(9'999) == 0);
  BOOST_TEST(geometry.blockOfPS(10'000) == 1);
  BOOST_TEST(geometry.blockOfPS(25'000) == 2);
  BOOST_TEST(geometry.blocksPerSlice(0) == geometry.blockCount());
  BOOST_TEST(geometry.blocksPerSlice(10'010) == 2);

  const auto firstLayerCells = geometry.axes(0).cellCount();
  BOOST_TEST(geometry.layerBucketBase(0) == 0U);
  BOOST_TEST(geometry.layerBucketBase(1) == static_cast<std::uint32_t>(geometry.blockCount() * firstLayerCells));
  BOOST_TEST(geometry.bucket(0, 1, 0, 0) == static_cast<std::uint32_t>(firstLayerCells));

  const auto slice = geometry.slice({.begin = 1, .end = 2}, 10'010);
  BOOST_TEST(slice.core.begin == 1);
  BOOST_TEST(slice.core.end == 2);
  BOOST_TEST(slice.halo.begin == 0);
  BOOST_TEST(slice.halo.end == 3);
}

BOOST_AUTO_TEST_CASE(PerLayerCountsOverrideTheDerivedOnes)
{
  // A layer's own counts, 0 for the derived.
  auto params = test::withToyFrames(GridParams{});
  params.blockWidth = 10'000;
  const GeometryGrid derived{25'000, params, defaultWindows()};
  params.phiBins[1] = 3;
  params.secondBins[1] = 5;
  params.secondBins[2] = 135;
  const GeometryGrid overridden{25'000, params, defaultWindows()};

  BOOST_TEST(overridden.axes(0).nPhi == derived.axes(0).nPhi);
  BOOST_TEST(overridden.axes(0).nSecond == derived.axes(0).nSecond);
  BOOST_TEST(overridden.axes(1).nPhi == 3);
  BOOST_TEST(overridden.axes(1).nSecond == 5);
  BOOST_TEST(overridden.axes(2).nPhi == derived.axes(2).nPhi);
  BOOST_TEST(overridden.axes(2).nSecond == 135);
  BOOST_TEST(overridden.layerBucketBase(2) == static_cast<std::uint32_t>(overridden.blockCount() * (overridden.axes(0).cellCount() + 15)));
  BOOST_TEST(!overridden.overflows());

  // More buckets than a 32-bit bucket number: flagged, never wrapped.
  params.phiBins[0] = 1 << 16;
  params.secondBins[0] = 1 << 16;
  const GeometryGrid huge{25'000, params, defaultWindows()};
  BOOST_TEST(huge.overflows());
}

namespace
{
/// Per link, the most layers a path through it visits, every path walked
/// one by one: the predicate's two sweeps checked by brute force.
template <DetectorLayout Det>
std::array<int, nLinks<Det>> longestPathsThrough()
{
  std::array<int, nLinks<Det>> longest{};
  std::vector<std::vector<std::size_t>> paths;
  paths.reserve(nLinks<Det>);
  for (std::size_t link = 0; link < nLinks<Det>; ++link) {
    paths.push_back({link});
  }
  while (!paths.empty()) {
    const auto path = std::move(paths.back());
    paths.pop_back();
    for (const auto link : path) {
      longest[link] = std::max(longest[link], static_cast<int>(path.size()) + 1);
    }
    for (std::size_t next = 0; next < nLinks<Det>; ++next) {
      if (Det::layerLinks[next].from == Det::layerLinks[path.back()].to) {
        auto longer = path;
        longer.push_back(next);
        paths.push_back(std::move(longer));
      }
    }
  }
  return longest;
}
} // namespace

BOOST_AUTO_TEST_CASE(TheLinkPredicateFollowsSkipsMaskAndPathLength)
{
  // One predicate for the stages' link lists and the grid's sizing.
  const auto all = takenLinks<Mixed>(IterationSettings{.windows = {}, .minClusters = 2});
  const auto noSkips = takenLinks<Mixed>(IterationSettings{.windows = {}, .minClusters = 2, .skipLinks = false});
  const auto inner = takenLinks<Mixed>(IterationSettings{.windows = {}, .minClusters = 2, .layerMask = 0b111U});
  const auto none = takenLinks<Mixed>(IterationSettings{.windows = {}, .minClusters = 2, .layerMask = 0b1U}); // a single layer: no link has both ends in it
  for (std::size_t link = 0; link < nLinks<Mixed>; ++link) {
    const auto& pair = Mixed::layerLinks[link];
    BOOST_TEST_CONTEXT("link " << link)
    {
      BOOST_TEST(all[link]);
      BOOST_TEST(noSkips[link] == !pair.skipsLayer);
      BOOST_TEST(inner[link] == (pair.from < 3U && pair.to < 3U));
      BOOST_TEST(!none[link]);
    }
  }
  // A link is taken only on a path of minClusters layers or more:
  // every length from 2 to one past the
  // longest path, which takes none.
  const auto longest = longestPathsThrough<Mixed>();
  for (int length = 2; std::cmp_less_equal(length, maxChainLayers<Mixed> + 1U); ++length) {
    const auto taken = takenLinks<Mixed>(IterationSettings{.windows = {}, .minClusters = length});
    for (std::size_t link = 0; link < nLinks<Mixed>; ++link) {
      BOOST_TEST_CONTEXT("length " << length << ", link " << link) { BOOST_TEST(taken[link] == (longest[link] >= length)); }
    }
  }
  // Within a mask, the paths are the mask's: the three innermost layers
  // hold a triplet's three on the one path 0 -> 1 -> 2, and no more; the
  // link 0 -> 2 past layer 1's end is on none.
  const auto triplets = takenLinks<Mixed>(IterationSettings{.windows = {}, .minClusters = 3, .layerMask = 0b111U});
  for (std::size_t link = 0; link < nLinks<Mixed>; ++link) {
    const auto& pair = Mixed::layerLinks[link];
    BOOST_TEST_CONTEXT("link " << link) { BOOST_TEST(triplets[link] == (inner[link] && pair.to == pair.from + 1U)); }
  }
  BOOST_TEST(std::ranges::none_of(takenLinks<Mixed>(IterationSettings{.windows = {}, .minClusters = 4, .layerMask = 0b111U}), std::identity{}));
}

BOOST_AUTO_TEST_CASE(ALayerIsSizedByTheQueriesIntoIt)
{
  const VertexConstraint beam{.z = 0., .sigmaZ = 6., .sigmaTransverse = 0.002};
  const SearchWindowParams low{.minPtGeV = 0.1, .bzKiloGauss = 5.};
  const SearchWindowParams high{.minPtGeV = 1., .bzKiloGauss = 5.};
  const auto pair = [&](LayerID from, LayerID to, const SearchWindowParams& windows) { return layerPairWindow<Toy>(from, to, windows, beam); };
  const auto sizes = [&](std::initializer_list<IterationSettings> iterations) {
    const std::vector<IterationSettings> list{iterations};
    return gridLayerWindows<Toy>(list, beam);
  };

  // Tracklet queries only: into each link's target; layer 0 is not searched
  // and has one cell per axis.
  const auto tracklets = sizes({IterationSettings{.windows = low, .minClusters = 3, .follow = false}});
  BOOST_TEST(!tracklets[0].queried);
  BOOST_TEST(tracklets[1].queried);
  BOOST_TEST(tracklets[1].phiHalfWidth == static_cast<double>(pair(0, 1, low).phiHalfWidth));
  BOOST_TEST(tracklets[2].secondHalfWidth == static_cast<double>(pair(1, 2, low).secondHalfWidth));
  const GeometryGrid grid{25'000, test::withToyFrames(GridParams{.blockWidth = 10'000}), tracklets};
  BOOST_TEST(grid.axes(0).nPhi == 1);
  BOOST_TEST(grid.axes(0).nSecond == 1);
  BOOST_TEST(grid.axes(1).nPhi == derivedBinCount(2 * std::numbers::pi_v<double>, 2. * tracklets[1].phiHalfWidth, kPhiCellsPerWindow));
  BOOST_TEST(!grid.overflows());

  // An iteration that follows also searches inward: layer 0 along 0 -> 1,
  // and layer 1 by the wider of its two queries.
  const auto extending = sizes({IterationSettings{.windows = low, .minClusters = 3, .follow = true}});
  BOOST_TEST(extending[0].queried);
  BOOST_TEST(extending[0].phiHalfWidth == static_cast<double>(pair(0, 1, low).phiHalfWidth));
  BOOST_TEST(extending[1].phiHalfWidth == std::max(static_cast<double>(pair(0, 1, low).phiHalfWidth), static_cast<double>(pair(1, 2, low).phiHalfWidth)));

  // A link only a high-pT iteration takes is sized at that pT; the low-pT one
  // takes only 0 -> 1 (its mask), so layer 1 is sized at its pT.
  const auto mixed = sizes({IterationSettings{.windows = high, .minClusters = 3, .follow = false}, IterationSettings{.windows = low, .minClusters = 2, .layerMask = 0b011U, .follow = false}});
  BOOST_TEST(mixed[2].phiHalfWidth == static_cast<double>(pair(1, 2, high).phiHalfWidth));
  BOOST_TEST(mixed[1].phiHalfWidth == static_cast<double>(pair(0, 1, low).phiHalfWidth));
  BOOST_TEST(mixed[1].phiHalfWidth > static_cast<double>(pair(0, 1, high).phiHalfWidth));

  // A layer the mask excludes from every link: not searched, one cell.
  const auto masked = sizes({IterationSettings{.windows = low, .minClusters = 2, .layerMask = 0b011U, .follow = false}});
  BOOST_TEST(!masked[2].queried);
  const GeometryGrid maskedGrid{25'000, test::withToyFrames(GridParams{.blockWidth = 10'000}), masked};
  BOOST_TEST(maskedGrid.axes(2).cellCount() == 1);

  // The same roads with the follow, which walks every link of the layout:
  // layer 2 is searched, in the 1 -> 2 pair window.
  const auto followed = sizes({IterationSettings{.windows = low, .minClusters = 2, .layerMask = 0b011U, .follow = true}});
  BOOST_TEST(followed[2].queried);
  BOOST_TEST(followed[2].phiHalfWidth == static_cast<double>(pair(1, 2, low).phiHalfWidth));
  BOOST_TEST(followed[0].phiHalfWidth == static_cast<double>(pair(0, 1, low).phiHalfWidth));
}

/// Queries of several widths: phi cells
/// balanced between the narrowest and the widest, the convention itself for
/// one width; the second coordinate's by the narrowest.
BOOST_AUTO_TEST_CASE(CellsServeTheNarrowestAndTheWidestQuery)
{
  // One width: the convention.
  BOOST_TEST(balancedBinCount(2 * std::numbers::pi_v<double>, 0.2, 0.2, kPhiCellsPerWindow) == derivedBinCount(2 * std::numbers::pi_v<double>, 0.2, kPhiCellsPerWindow));
  BOOST_TEST(balancedBinCount(2 * std::numbers::pi_v<double>, 0., 0.2, kPhiCellsPerWindow) == derivedBinCount(2 * std::numbers::pi_v<double>, 0.2, kPhiCellsPerWindow));
  // Two widths: between their conventions, each exceeding it by the same factor.
  const double a = 0.01;
  const double b = 0.4;
  const double c = kPhiCellsPerWindow;
  const auto count = balancedBinCount(2 * std::numbers::pi_v<double>, a, b, c);
  BOOST_TEST(count > derivedBinCount(2 * std::numbers::pi_v<double>, b, c));
  BOOST_TEST(count < derivedBinCount(2 * std::numbers::pi_v<double>, a, c));
  const double w = (2 * std::numbers::pi_v<double>) / count;
  const double rows = (1. + (w / a)) / (1. + (1. / c));
  const double cells = ((b / w) + 1.) / (c + 1.);
  BOOST_TEST(std::abs(rows - cells) < 0.1 * rows);
  // With c = 1, the geometric mean.
  BOOST_TEST(balancedBinCount(1., 0.01, 0.04, 1.) == derivedBinCount(1., 0.02, 1.));

  std::array<GridLayerWindow, 3> windows{};
  for (auto& window : windows) {
    window = {.phiHalfWidth = 0.2, .secondHalfWidth = 2., .queried = true, .narrowPhiHalfWidth = 0.005, .narrowSecondHalfWidth = 0.5};
  }
  const GeometryGrid grid{25'000, test::withToyFrames(GridParams{.blockWidth = 10'000}), windows};
  BOOST_TEST(grid.axes(1).nPhi == balancedBinCount(2 * std::numbers::pi_v<double>, 0.01, 0.4, kPhiCellsPerWindow));
  const double extent = secondMax<Toy>(1) - secondMin<Toy>(1);
  BOOST_TEST(grid.axes(1).nSecond == derivedBinCount(extent, 1., kSecondCellsPerWindow));
}

BOOST_AUTO_TEST_CASE(CountsPastTheNumberingAreRefusedNotWrapped)
{
  BOOST_TEST(derivedBinCount(1., 0.5, 2.) == 4);
  BOOST_TEST(derivedBinCount(1., 1.e-12, 1.) == 0); // 10^12 cells: past int32, refused
  // A derived count past int32: the layout is refused, every axis still one
  // cell or more.
  std::array<GridLayerWindow, 3> tiny{};
  for (auto& window : tiny) {
    window = {.phiHalfWidth = 1.e-12, .secondHalfWidth = 1., .queried = true};
  }
  const GeometryGrid refused{25'000, test::withToyFrames(GridParams{.blockWidth = 10'000}), tiny};
  BOOST_TEST(refused.overflows());
  BOOST_TEST(refused.axes(0).nPhi >= 1);
  // A block's cells past a 32-bit signed offset, from counts that each fit.
  auto params = test::withToyFrames(GridParams{});
  params.blockWidth = 10'000;
  params.phiBins[0] = INT32_MAX;
  params.secondBins[0] = 2;
  const GeometryGrid cells{25'000, params, defaultWindows()};
  BOOST_TEST(cells.overflows());
}

BOOST_AUTO_TEST_CASE(BuilderProducesPrefixOffsetsAndStablePlacements)
{
  std::array<std::uint32_t, 6> offsets{};
  MeasurementGridBuilder builder{offsets};
  constexpr std::array<std::uint32_t, 5> buckets{2U, 0U, 2U, 1U, 2U};
  for (const auto bucket : buckets) {
    builder.count(bucket);
  }
  builder.prefixSum();

  std::array<std::uint32_t, buckets.size()> placements{};
  for (std::size_t i = 0; i < buckets.size(); ++i) {
    placements[i] = builder.place(buckets[i]);
  }
  builder.finish();

  constexpr std::array<std::uint32_t, buckets.size()> expectedPlacements{2U, 0U, 3U, 1U, 4U};
  constexpr std::array<std::uint32_t, 6> expectedOffsets{0U, 1U, 2U, 5U, 5U, 5U};
  for (std::size_t i = 0; i < placements.size(); ++i) {
    BOOST_TEST(placements[i] == expectedPlacements[i]);
  }
  for (std::size_t i = 0; i < offsets.size(); ++i) {
    BOOST_TEST(offsets[i] == expectedOffsets[i]);
  }
  BOOST_TEST(offsets.back() == static_cast<std::uint32_t>(buckets.size()));
}

BOOST_AUTO_TEST_CASE(QueryWrapsPhiAcrossZero)
{
  auto params = test::withToyFrames(GridParams{});
  params.blockWidth = 100'000;
  for (std::size_t layer = 0; layer < GridParams::kMaxLayers; ++layer) {
    params.phiBins[layer] = 8;
    params.secondBins[layer] = 1;
  }
  const GeometryGrid geometry{100'000, params, defaultWindows()};
  BOOST_REQUIRE(geometry.axes(1).nPhi == 8);

  constexpr std::array<Scalar, 3> phis{0.01F, (2 * std::numbers::pi_v<Scalar>)-0.01F, std::numbers::pi_v<Scalar>};
  constexpr std::array<std::uint32_t, 3> ids{10U, 11U, 12U};
  std::vector<std::uint32_t> offsets(static_cast<std::size_t>(geometry.bucketCount()) + 1U);
  MeasurementStore columns;
  resizeMeasurementColumns(columns, phis.size());
  test::buildGridRows(
    phis.size(), offsets, [&](std::size_t i) { return geometry.bucketAt(1, 0U, phis[i], 0.F); },
    [&](std::size_t i, std::uint32_t row) {
      columns.widthPS[row] = 1U;
      columns.clusterIndex[row] = ids[i];
    });

  const Grid grid{measurementView(columns), geometry, std::span<const std::uint32_t>{offsets}};
  const CandidateWindow window{.phiCentre = 0.F, .phiHalfWidth = 0.1F, .secondCentre = 0.F, .secondHalfWidth = 2.F};
  BOOST_TEST(spatiallyCompatible(window, (2 * std::numbers::pi_v<Scalar>)-0.01F, 0.F));
  BOOST_TEST(!spatiallyCompatible(window, std::numbers::pi_v<Scalar>, 0.F));
  BOOST_TEST(!spatiallyCompatible(window, 0.F, 3.F));
  std::vector<std::uint32_t> visited;
  grid.visitCandidates(1, window, geometry.allBlocks(), [&](std::uint32_t row) { visited.push_back(columns.clusterIndex[row]); });
  std::ranges::sort(visited);

  BOOST_REQUIRE(visited.size() == 2U);
  BOOST_TEST(visited[0] == 10U);
  BOOST_TEST(visited[1] == 11U);
}

BOOST_AUTO_TEST_CASE(QueryCoversEveryPhiCellWhenAWrappingWindowEndsInItsStartCell)
{
  auto params = test::withToyFrames(GridParams{});
  params.blockWidth = 100'000;
  for (std::size_t layer = 0; layer < GridParams::kMaxLayers; ++layer) {
    params.phiBins[layer] = 3;
    params.secondBins[layer] = 1;
  }
  const GeometryGrid geometry{100'000, params, defaultWindows()};
  BOOST_REQUIRE(geometry.axes(1).nPhi == 3);

  // Cells are 2pi/3 wide. A window of +-3 rad about 0.5 leaves out only
  // (3.5, 3.78), so both of its ends fall in cell 1 -- yet it spans all
  // three cells, including the seed's own direction in cell 0.
  const CandidateWindow window{.phiCentre = 0.5F, .phiHalfWidth = 3.F, .secondCentre = 0.F, .secondHalfWidth = 100.F};
  constexpr std::array<GridPoint, 3> points{{{.layer = 1, .phi = 0.5F, .second = 0.F}, {.layer = 1, .phi = 3.F, .second = 0.F}, {.layer = 1, .phi = 5.F, .second = 0.F}}};
  for (const auto& point : points) {
    BOOST_TEST(spatiallyCompatible(window, point.phi, point.second));
  }

  const auto visited = visitedPoints(geometry, 1, points, window);
  BOOST_TEST(visited == (std::vector<std::uint32_t>{0U, 1U, 2U}), boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(QueryReachesDiscSensorCornersBeyondTheRingEdge)
{
  const auto placement = toyGeometry<Mixed>();
  constexpr std::size_t kPositiveDisc = 4;
  constexpr auto kDisc = static_cast<LayerID>(Detector<Mixed>::nBarrelLayers + kPositiveDisc);
  const auto& ring = Mixed::discSpecs[kPositiveDisc].rings[1];
  const auto& sensor = Mixed::discSpecs[kPositiveDisc].sensor;
  const SensorID chip = firstSensorOfLayer<Mixed>(kDisc) + Mixed::discSpecs[kPositiveDisc].rings[0].nModules;
  BOOST_REQUIRE(placement.sensor(chip).layer == kDisc);
  BOOST_REQUIRE(std::abs(std::hypot(placement.sensor(chip).trans[0], placement.sensor(chip).trans[1]) - ring.r) < 1e-9);

  // Just inside the outer ring's outer edge, at a corner: further out than the
  // ring's nominal r + halfY, so the axis must not end there.
  const double loc1 = sensor.halfY - 0.01;
  const double loc0 = sensor.halfX - 0.01;
  const auto point = gridPoint<Mixed>(placement, chip, loc0, loc1);
  BOOST_REQUIRE(point.layer == kDisc);
  BOOST_REQUIRE(static_cast<double>(point.second) > ring.r + sensor.halfY + 0.05);
  BOOST_TEST(static_cast<double>(point.second) <= secondMax<Mixed>(kDisc));

  // A window this narrow lies wholly beyond the nominal edge.
  const CandidateWindow window{.phiCentre = point.phi, .phiHalfWidth = 0.01F, .secondCentre = point.second, .secondHalfWidth = 0.02F};
  const GridGeometry<Mixed> geometry{100'000, test::withToyFrames(GridParams{.blockWidth = 100'000}), test::oneIterationWindows<Mixed>(SearchWindowParams{}, VertexConstraint{})};
  const auto visited = visitedPoints(geometry, kDisc, std::array<GridPoint, 1>{point}, window);
  BOOST_TEST(visited.size() == 1U);
}

BOOST_AUTO_TEST_CASE(QueryVisitsTheOverlappingBoxesInsideTheSearchBracket)
{
  const auto params = oneCellParams(100'000);
  const GeometryGrid geometry{100'000, params, defaultWindows()};

  // Boxes of 10 ps starting at 0, 5, 10 and 20; the seed's is [10, 20).
  constexpr std::array<std::uint32_t, 4> starts{0U, 5U, 10U, 20U};
  constexpr std::uint32_t kWidth = 10U;
  std::vector<std::uint32_t> offsets(static_cast<std::size_t>(geometry.bucketCount()) + 1U);
  const auto bucket = geometry.bucketAt(1, 0U, 0.F, 0.F);
  MeasurementStore columns;
  resizeMeasurementColumns(columns, starts.size());
  test::buildGridRows(
    starts.size(), offsets, [&](std::size_t) { return bucket; },
    [&](std::size_t i, std::uint32_t row) {
      columns.timePS[row] = starts[i];
      columns.widthPS[row] = kWidth;
      columns.clusterIndex[row] = static_cast<std::uint32_t>(i);
    });

  const Grid grid{measurementView(columns), geometry, std::span<const std::uint32_t>{offsets}};
  const LayerPairWindow pair{.phiHalfWidth = std::numbers::pi_v<Scalar>, .secondHalfWidth = 2.F};
  const auto window = candidateWindow(pair, TimeEstPS{.start = 10U, .width = kWidth}, 0.F, 0.F, 0.F, kWidth);
  BOOST_TEST(window.earliestPS == 1);
  BOOST_TEST(window.latestPS == 19);
  std::vector<std::uint32_t> visited;
  grid.visitCandidates(1, window, geometry.allBlocks(), [&](std::uint32_t row) { visited.push_back(columns.clusterIndex[row]); });

  // [0, 10) ends where the seed starts; [20, 30) starts where it ends.
  BOOST_TEST(visited == (std::vector<std::uint32_t>{1U, 2U}), boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(LanesShareAQueryAndVisitEachRowOnce)
{
  // A cell denser than the scan's limit (its bracket searched) and
  // one under it (scanned); boxes of 10 ps every 5 ps.
  const auto params = oneCellParams(100'000);
  const GeometryGrid geometry{100'000, params, defaultWindows()};
  for (const std::uint32_t rows : {40U, 5U}) {
    std::vector<std::uint32_t> offsets(static_cast<std::size_t>(geometry.bucketCount()) + 1U);
    const auto bucket = geometry.bucketAt(1, 0U, 0.F, 0.F);
    MeasurementStore columns;
    resizeMeasurementColumns(columns, rows);
    test::buildGridRows(
      rows, offsets, [&](std::size_t) { return bucket; },
      [&](std::size_t i, std::uint32_t row) {
        columns.timePS[row] = 5U * i;
        columns.widthPS[row] = 10U;
        columns.clusterIndex[row] = i;
      });

    const Grid grid{measurementView(columns), geometry, std::span<const std::uint32_t>{offsets}};
    const LayerPairWindow pair{.phiHalfWidth = std::numbers::pi_v<Scalar>, .secondHalfWidth = 2.F};
    const auto window = candidateWindow(pair, TimeEstPS{.start = 10U, .width = 60U}, 0.F, 0.F, 0.F, 10U);
    std::vector<std::uint32_t> single;
    grid.visitCandidates(1, window, geometry.allBlocks(), [&](std::uint32_t row) { single.push_back(columns.clusterIndex[row]); });
    BOOST_TEST_REQUIRE(single.size() >= 2U);
    for (const std::uint32_t count : {1U, 2U, 3U, 8U, 32U}) {
      std::vector<std::uint32_t> shared;
      for (std::uint32_t lane = 0; lane < count; ++lane) {
        grid.visitCandidates(1, window, geometry.allBlocks(), [&](std::uint32_t row) { shared.push_back(columns.clusterIndex[row]); }, Lanes{.lane = lane, .count = count});
      }
      std::ranges::sort(shared);
      BOOST_TEST(shared == single, boost::test_tools::per_element());
    }
  }
}

/// One lane of a group walking in step, on the host: the walk itself calls
/// no collective, so a lane's calls can be recorded alone.
template <std::uint32_t L>
struct LaneOf {
  static constexpr std::uint32_t kLanes = L;
  std::uint32_t index = 0U;
  [[nodiscard]] std::uint32_t lane() const { return index; }
};

template <std::uint32_t L>
static void checkGroupInStep(const Grid& grid, const MeasurementStore& columns, const CandidateWindow& window, const GeometryGrid& geometry)
{
  std::size_t steps = 0U;
  for (std::uint32_t lane = 0; lane < L; ++lane) {
    std::vector<std::uint32_t> own;
    grid.visitCandidates(1, window, geometry.allBlocks(), [&](std::uint32_t row) { own.push_back(columns.clusterIndex[row]); }, Lanes{.lane = lane, .count = L});
    std::vector<std::uint32_t> inStep;
    std::size_t calls = 0U;
    grid.visitCandidates(
      1, window, geometry.allBlocks(),
      [&](std::uint32_t row, bool has) {
        ++calls;
        if (has) {
          inStep.push_back(columns.clusterIndex[row]);
        }
      },
      LaneOf<L>{.index = lane});
    BOOST_TEST_CONTEXT("lanes " << L << ", lane " << lane)
    {
      BOOST_TEST(inStep == own, boost::test_tools::per_element());
      if (lane == 0U) {
        steps = calls;
      }
      BOOST_TEST(calls == steps);
    }
  }
}

BOOST_AUTO_TEST_CASE(AGroupInStepGivesEachLaneItsLanesRowsAtEveryStep)
{
  // Phi cells 6, 7 and then 0, 1 after the wrap in the window, cell
  // 3 outside it, cell 7 empty; rows at the time bound's ends and a box
  // ending exactly at afterPS. The window: earliest 1, latest 69, after 10.
  const auto params = fixedCellParams(100'000, 8, 1);
  const GeometryGrid geometry{100'000, params, defaultWindows()};
  const Scalar cellWidth = 2.F * std::numbers::pi_v<Scalar> / 8.F;
  // Per cell its rows' box starts (ascending) and the row count: a sparse
  // set (the range scanned) and a dense one (each cell searched).
  struct Row {
    std::int32_t cell = 0;
    std::uint32_t start = 0U;
    std::uint32_t width = 10U;
  };
  const auto rowsOf = [](std::uint32_t perCell) {
    std::vector<Row> rows;
    for (const std::int32_t cell : {0, 1, 3, 6}) {
      rows.push_back({.cell = cell, .start = 1U, .width = 9U}); // box ends at afterPS: not a candidate
      for (std::uint32_t i = 0; i < perCell; ++i) {
        rows.push_back({.cell = cell, .start = 2U + ((67U * i) / perCell)});
      }
      rows.push_back({.cell = cell, .start = 69U}); // the bound's last start
      rows.push_back({.cell = cell, .start = 70U}); // past it
    }
    return rows;
  };
  for (const std::uint32_t perCell : {2U, 37U}) {
    const auto rows = rowsOf(perCell);
    std::vector<std::uint32_t> offsets(static_cast<std::size_t>(geometry.bucketCount()) + 1U);
    MeasurementStore columns;
    resizeMeasurementColumns(columns, rows.size());
    const auto bucketOf = [&](const Row& row) { return geometry.bucketAt(1, row.start, (static_cast<Scalar>(row.cell) + 0.5F) * cellWidth, 0.F); };
    test::buildGridRows(
      rows.size(), offsets, [&](std::size_t i) { return bucketOf(rows[i]); },
      [&](std::size_t i, std::uint32_t at) {
        columns.timePS[at] = rows[i].start;
        columns.widthPS[at] = rows[i].width;
        columns.clusterIndex[at] = static_cast<std::uint32_t>(i);
      });

    const Grid grid{measurementView(columns), geometry, std::span<const std::uint32_t>{offsets}};
    const LayerPairWindow pair{.phiHalfWidth = 1.5F * cellWidth, .secondHalfWidth = 2.F};
    const auto window = candidateWindow(pair, TimeEstPS{.start = 10U, .width = 60U}, 0.1F, 0.F, 0.F, 10U);
    BOOST_TEST_REQUIRE(window.earliestPS == 1);
    BOOST_TEST_REQUIRE(window.latestPS == 69);
    std::vector<std::uint32_t> single;
    grid.visitCandidates(1, window, geometry.allBlocks(), [&](std::uint32_t row) { single.push_back(columns.clusterIndex[row]); });
    BOOST_TEST_REQUIRE(single.size() == 3U * (perCell + 1U));
    BOOST_TEST_CONTEXT("rows per cell " << perCell)
    {
      checkGroupInStep<2>(grid, columns, window, geometry);
      checkGroupInStep<3>(grid, columns, window, geometry);
      checkGroupInStep<8>(grid, columns, window, geometry);
      checkGroupInStep<32>(grid, columns, window, geometry);
      checkGroupInStep<64>(grid, columns, window, geometry);
    }
  }
}

BOOST_AUTO_TEST_CASE(QueryVisitsOrderedSecondRunsWithoutDuplicates)
{
  auto params = test::withToyFrames(GridParams{});
  params.blockWidth = 100'000;
  for (std::size_t layer = 0; layer < GridParams::kMaxLayers; ++layer) {
    params.phiBins[layer] = 1;
    params.secondBins[layer] = 8;
  }
  const GeometryGrid geometry{100'000, params, defaultWindows()};
  BOOST_REQUIRE(geometry.axes(1).nSecond == 8);

  constexpr std::size_t kBins = 8;
  const auto& axes = geometry.axes(1);
  const Scalar width = 1.F / axes.secondInvWidth;
  std::array<GridPoint, kBins> points{};
  for (std::size_t i = 0; i < points.size(); ++i) {
    points[i] = {.layer = 1, .phi = 0.F, .second = axes.secondMin + ((static_cast<Scalar>(i) + 0.5F) * width)};
  }

  std::vector<std::uint32_t> offsets(static_cast<std::size_t>(geometry.bucketCount()) + 1U);
  MeasurementStore columns;
  resizeMeasurementColumns(columns, kBins);
  test::buildGridRows(
    points.size(), offsets, [&](std::size_t i) { return geometry.bucketAt(points[i].layer, 0U, points[i].phi, points[i].second); },
    [&](std::size_t i, std::uint32_t row) {
      columns.widthPS[row] = 1U;
      columns.clusterIndex[row] = static_cast<std::uint32_t>(i);
    });

  const Grid grid{measurementView(columns), geometry, std::span<const std::uint32_t>{offsets}};
  const Scalar first = points[1].second;
  const std::array<SecondRun, 3> runs{{
    {.low = first - (0.2F * width), .high = first - (0.1F * width)},
    {.low = first + (0.1F * width), .high = first + (0.2F * width)},
    {.low = points[5].second - (0.1F * width), .high = points[5].second + (0.1F * width)},
  }};
  const CandidateWindow window{.phiCentre = 0.F, .phiHalfWidth = std::numbers::pi_v<Scalar>};
  std::vector<std::uint32_t> visited;
  grid.visitCandidates(1, window, runs.data(), static_cast<std::uint32_t>(runs.size()), geometry.allBlocks(), [&](std::uint32_t row) { visited.push_back(columns.clusterIndex[row]); });

  BOOST_TEST(visited == (std::vector<std::uint32_t>{1U, 5U}), boost::test_tools::per_element());
}

} // namespace o2::gtrk

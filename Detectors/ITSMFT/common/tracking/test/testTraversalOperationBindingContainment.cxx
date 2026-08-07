// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// M5c source-level guard: the four shared tracklet/cell/neighbour/road
// hot-loop entry points (TrackerTraits::computeLayerTracklets/
// computeLayerCells/findCellsNeighbours/findRoads, TrackerTraits.cxx) used to
// each re-run their own runtime StateFamily/SurfaceKind dispatch to pick which
// leaf orchestration to invoke. The binding replaces all four with a
// direct call into the already-bound TraversalOperationBinding
// (mTraversalOperation, established once per iteration by
// bindTraversalOperation() -- see TrackerTraits.h/.cxx). This test proves
// that conversion by construction, not just by behavior: it extracts each of
// the four methods' own source text (from its `void
// TrackerTraits::<name>(...)` signature to the matching closing
// brace, via simple brace counting -- these four bodies contain no nested
// class/lambda braces that would confuse that count), strips `//` line
// comments (so explanatory prose about what these methods *used to* do,
// referencing the same tokens, cannot produce a false failure -- only actual
// code matters here), and greps what remains for the tokens that would mean
// a SurfaceKind/StateFamily selection branch crept back in: SurfaceKind,
// StateFamily, and `if constexpr`.
// The corresponding public-header containment test already proves the
// analogous claim for public *headers*; this file is the .cxx/hot-loop
// counterpart M5c adds.

#define BOOST_TEST_MODULE ITSMFT TraversalOperationBindingContainment
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace
{

std::string readFile(const std::string& path)
{
  std::ifstream input{path};
  BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path);
  std::ostringstream buffer;
  buffer << input.rdbuf();
  return buffer.str();
}

std::string readTrackerTraitsSource()
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  return readFile(testDirectory + "/../src/TrackerTraits.cxx");
}

std::string readTrackerTraitsHeader()
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  return readFile(testDirectory + "/../include/ITSMFTTracking/TrackerTraits.h");
}

/// Extracts the body (inclusive of the opening/closing braces) of the first
/// `void TrackerTraits::<methodName>(` definition found in `source`,
/// by counting braces from that signature's opening `{` to its matching
/// closing `}`. `methodName` must be immediately followed by `(` in the
/// match (a plain string search on the qualified name, not a regex), so
/// `computeLayerTracklets` never matches `computeLayerTrackletsForKind`.
std::string extractMethodBody(const std::string& source, const std::string& methodName)
{
  const std::string signature = "TrackerTraits::" + methodName + "(";
  const auto signaturePos = source.find(signature);
  BOOST_REQUIRE_MESSAGE(signaturePos != std::string::npos, "cannot find " << signature << " in TrackerTraits.cxx");
  const auto openBrace = source.find('{', signaturePos);
  BOOST_REQUIRE_MESSAGE(openBrace != std::string::npos, "cannot find opening brace for " << methodName);

  int depth = 0;
  size_t pos = openBrace;
  for (; pos < source.size(); ++pos) {
    if (source[pos] == '{') {
      ++depth;
    } else if (source[pos] == '}') {
      --depth;
      if (depth == 0) {
        break;
      }
    }
  }
  BOOST_REQUIRE_MESSAGE(depth == 0, "unbalanced braces while extracting " << methodName);
  return source.substr(openBrace, pos - openBrace + 1);
}

/// Strips `//`-style line comments (this file's own convention throughout
/// TrackerTraits.cxx -- no `/* */` block comments are used there), so a
/// token search below only ever matches actual code, never explanatory prose
/// about what a method used to do.
std::string stripLineComments(const std::string& body)
{
  std::istringstream lines{body};
  std::ostringstream stripped;
  std::string line;
  while (std::getline(lines, line)) {
    const auto commentPos = line.find("//");
    stripped << (commentPos == std::string::npos ? line : line.substr(0, commentPos)) << '\n';
  }
  return stripped.str();
}

bool mentionsToken(const std::string& body, const std::string& token)
{
  const std::regex tokenRegex{"\\b" + token + "\\b"};
  return std::regex_search(body, tokenRegex);
}

} // namespace

BOOST_AUTO_TEST_CASE(SharedHotLoopEntryPointsContainNoKindOrFamilyDispatchBranch)
{
  const auto source = readTrackerTraitsSource();
  const std::array<std::string, 4> hotLoopMethods{
    "computeLayerTracklets", "computeLayerCells", "findCellsNeighbours", "findRoads"};
  const std::array<std::string, 3> forbiddenTokens{"SurfaceKind", "StateFamily", "constexpr"};

  for (const auto& method : hotLoopMethods) {
    const auto body = extractMethodBody(source, method);
    BOOST_REQUIRE_GT(body.size(), 0u);
    const auto code = stripLineComments(body);
    for (const auto& token : forbiddenTokens) {
      BOOST_CHECK_MESSAGE(!mentionsToken(code, token),
                          "TrackerTraits::" << method << "() still mentions " << token
                                            << " in code -- a SurfaceKind/StateFamily dispatch branch leaked back"
                                            << " into the shared hot loop instead of staying inside"
                                            << " bindTraversalOperation()");
    }
    // Every one of these methods must instead route through the bound
    // operation (TraversalOperationBinding, TrackerTraits.h) established by
    // bindTraversalOperation() in initialiseTimeFrame().
    BOOST_CHECK_MESSAGE(code.find("mTraversalOperation") != std::string::npos,
                        "TrackerTraits::" << method << "() no longer invokes mTraversalOperation");
  }
}

BOOST_AUTO_TEST_CASE(BindTraversalOperationIsTheSoleProducerOfTheBoundOperation)
{
  // The inverse of the check above: mTraversalOperation's callables must be
  // assigned in exactly one place (bindTraversalOperation()), never
  // piecemeal from the hot loops themselves -- otherwise the "established
  // once per iteration" contract (TraversalOperationBinding's own doc) would
  // not actually hold.
  const auto source = readTrackerTraitsSource();
  const std::string signature = "TrackerTraits::bindTraversalOperation(";
  const auto signaturePos = source.find(signature);
  BOOST_REQUIRE_MESSAGE(signaturePos != std::string::npos, "cannot find " << signature);
  const auto bindBody = extractMethodBody(source, "bindTraversalOperation");
  const auto bindBodyStart = source.find(bindBody, signaturePos);
  BOOST_REQUIRE_MESSAGE(bindBodyStart != std::string::npos, "cannot relocate bindTraversalOperation() body");
  const auto bindBodyEnd = bindBodyStart + bindBody.size();

  const std::regex assignment{"mTraversalOperation\\.(computeTracklets|computeCells|findNeighbours|findRoads)\\s*="};
  size_t count = 0;
  for (auto it = std::sregex_iterator(source.begin(), source.end(), assignment); it != std::sregex_iterator(); ++it) {
    ++count;
    const auto matchPos = static_cast<size_t>(it->position());
    BOOST_CHECK_MESSAGE(matchPos >= bindBodyStart && matchPos < bindBodyEnd,
                        "mTraversalOperation callable assignment found outside bindTraversalOperation() at offset "
                          << matchPos << ": " << it->str());
  }
  // Four callables (computeTracklets/computeCells/findNeighbours/findRoads),
  // each assigned exactly twice (once per SurfaceKind branch inside
  // bindTraversalOperation()'s own `if constexpr`, which is explicitly
  // allowed -- see that method's own doc).
  BOOST_CHECK_EQUAL(count, 8u);
}

BOOST_AUTO_TEST_CASE(TraversalOperationBindingUsesNonAllocatingMemberFunctionPointers)
{
  // M5c corrective slice: TraversalOperationBinding's four members must be
  // plain pointers-to-member (TrackerTraits::*), never std::function or any
  // other type-erasing/potentially-heap-allocating callable wrapper. Proven
  // two ways: the header never mentions std::function/<functional> at all
  // (a coarse but sufficient guard, since nothing else in this header has a
  // legitimate reason to need either), and each of the four member
  // declarations is spelled as a pointer-to-member-function type.
  const auto header = readTrackerTraitsHeader();
  BOOST_CHECK_MESSAGE(header.find("std::function") == std::string::npos,
                      "TrackerTraits.h mentions std::function -- TraversalOperationBinding must use "
                      "non-allocating member-function pointers instead");
  BOOST_CHECK_MESSAGE(header.find("<functional>") == std::string::npos,
                      "TrackerTraits.h includes <functional> -- no longer needed once "
                      "TraversalOperationBinding uses member-function pointers");

  const std::regex trackerTraitsMemberPointer{"\\(TrackerTraits::\\*\\)"};
  auto matches = std::distance(std::sregex_iterator(header.begin(), header.end(), trackerTraitsMemberPointer), std::sregex_iterator());
  // ComputeTrackletsFn/ComputeCellsFn/FindNeighboursFn/FindRoadsFn: one
  // pointer-to-member-function type alias each, inside TraversalOperationBinding.
  BOOST_CHECK_EQUAL(matches, 4);
}

BOOST_AUTO_TEST_CASE(EightNonTemplateWrapperTargetsForwardToKindSpecificLeafImplementations)
{
  // The eight (operation, SurfaceKind) wrapper targets a bound
  // pointer-to-member may point to (TrackerTraits.h) must each be a thin,
  // non-template forwarder to the corresponding pre-existing *ForKind<SurfaceKind>
  // leaf implementation -- never a reimplementation, and never itself
  // templated (a template member function cannot be the unique target of a
  // plain, non-template pointer-to-member-function type).
  const auto source = readTrackerTraitsSource();
  const std::array<std::pair<std::string, std::string>, 8> wrapperToLeaf{{
    {"computeLayerTrackletsCylinder", "computeLayerTrackletsForKind<SurfaceKind::Cylinder>"},
    {"computeLayerTrackletsDisk", "computeLayerTrackletsForKind<SurfaceKind::Disk>"},
    {"computeLayerCellsCylinder", "computeLayerCellsForKind<SurfaceKind::Cylinder>"},
    {"computeLayerCellsDisk", "computeLayerCellsForKind<SurfaceKind::Disk>"},
    {"findCellsNeighboursCylinder", "findCellsNeighboursForKind<SurfaceKind::Cylinder>"},
    {"findCellsNeighboursDisk", "findCellsNeighboursForKind<SurfaceKind::Disk>"},
    {"findRoadsCylinder", "findRoadsForKind<SurfaceKind::Cylinder>"},
    {"findRoadsDisk", "findRoadsForKind<SurfaceKind::Disk>"},
  }};
  for (const auto& [wrapper, leafCall] : wrapperToLeaf) {
    const auto body = extractMethodBody(source, wrapper);
    BOOST_REQUIRE_GT(body.size(), 0u);
    BOOST_CHECK_MESSAGE(body.find(leafCall) != std::string::npos,
                        "TrackerTraits::" << wrapper << "() does not forward to " << leafCall);
  }
}

BOOST_AUTO_TEST_CASE(ProcessNeighboursKeepsSurfaceKindSelectionCompileTime)
{
  const auto source = readTrackerTraitsSource();
  BOOST_CHECK(source.find("template <SurfaceKind Kind, typename InputSeed>") != std::string::npos);
  BOOST_CHECK(source.find("params.kind") == std::string::npos);
  BOOST_CHECK(source.find("processNeighbours<Kind>(") != std::string::npos);
}

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// M4 (GenericTrackingEngineMigration.md; ADR 0007 decisions 7-8) compile/grep
// boundary guard: TransitionPolicyTag is a temporary legacy hot-loop-dispatch
// implementation detail, confined to ITSMFTTracking/detail/ (headers whose
// consumers are, by construction, only other ITSMFTTracking translation
// units -- never a DPL workflow spec or another detector). No public,
// non-detail header under ITSMFTTracking/ may name TransitionPolicyTag.
//
// M4b correction: StateFamily (ITSMFTTracking/StateFamily.h) is a permanent
// public concept, but only as a state-representation-local classification
// (SurfaceKinematicState::family, SurfaceLinearizationReference::family, and
// the barrel::/forward:: operations that validate them). It must never
// reappear as a substitute transition-policy/dispatch key in a public
// topology, traversal, scheduling, binding-count, or adapter-facing
// observability API -- those must derive everything from
// SurfaceDescriptor/SurfaceKind directly. This file therefore also scans for
// that narrower misuse, with its own, much smaller allowlist of genuine
// state-representation headers.
//
// Both scans read every *.h file directly under include/ITSMFTTracking/ (not
// recursing into detail/, which is expected to mention the tag) so the guard
// also catches a *future* header that starts misusing either concept -- not
// just the files known at the time this test was written.

#define BOOST_TEST_MODULE ITSMFT TransitionPolicyTagContainment
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <fstream>
#include <regex>
#include <set>
#include <string>

namespace fs = std::filesystem;

namespace
{

// Explicit, reasoned allowlist: every one of these is confirmed (repo-wide
// grep, M4 audit) to have zero consumers outside this same library -- never
// included by a DPL workflow spec, another detector, or any other component
// -- so the tag-templated internals they still declare (all unaffected leaf
// hot-loop machinery this milestone leaves behind for M5, per
// GenericTrackingEngineMigration.md's M4 entry) are not part of any
// adapter-facing contract, even though the file itself is not physically
// under detail/.
const std::set<std::string>& transitionPolicyTagAllowlist()
{
  static const std::set<std::string> files = {
    // TrackerTraits<NLayers> is the legacy CA participant implementation
    // (M2/M3): every tag-templated member left here is private
    // (processNeighbours<Tag> included, moved private by M4); the public
    // surface no longer mentions the tag at all (M4b removed
    // getPolicyBindingCount(StateFamily), the last public trace of it).
    "TrackerTraits.h",
    // Tag-templated legacy index-table binder (bindIndexTableConfiguration<Tag>);
    // declaration-only, host-only, called exclusively from TrackerTraits.cxx.
    "IndexTableConfiguration.h",
  };
  return files;
}

// M4b: much narrower than the tag allowlist above -- every file here is
// genuinely about a SurfaceKinematicState/SurfaceLinearizationReference
// instance's own family, never about which topology transition, traversal
// schedule, or policy binding is active.
const std::set<std::string>& stateFamilyAllowlist()
{
  static const std::set<std::string> files = {
    // Defines StateFamily itself and its SurfaceKind bridge.
    "StateFamily.h",
    // SurfaceKinematicState::family / SurfaceLinearizationReference::family:
    // the state's own representation, not a dispatch key.
    "SurfaceKinematicState.h",
    "SurfaceKinematicStateLegacyAdapters.h",
    "SurfaceLinearizationReference.h",
    // barrel::/forward:: operations validate the family of the
    // SurfaceKinematicState/SurfaceLinearizationReference they are handed
    // before doing family-specific physics -- a precondition check on the
    // state, not a topology/scheduling decision.
    "BarrelSurfaceStateOperations.h",
    "ForwardSurfaceStateOperations.h",
    // Adapter-edge refit conversion inspects the representation family of
    // the generic state before exporting detector compatibility data.
    "DetectorTrackingOperationAdapterSupport.h",
    // M5d: Propagator::convertFamily/propagateToMeasurement route on the
    // target SurfaceDescriptor::kind-derived StateFamily and, on mismatch,
    // perform a real state-representation conversion between families --
    // exactly the same "operate on the state's own family classification"
    // role BarrelSurfaceStateOperations.h/ForwardSurfaceStateOperations.h
    // already have below, never a topology/scheduling/dispatch key.
    "Propagator.h",
    "NativeRefitDriver.h",
  };
  return files;
}

// Word-boundary match: `token` must appear as a standalone identifier, not
// as a prefix of a longer, unrelated one (e.g. `StateFamily` must not match
// `TraversalFailureReason::StateFamilyMismatch`, an enumerator name that
// documents a failure mode in English but names no StateFamily type).
bool lineMentionsToken(const std::string& line, const std::regex& tokenRegex)
{
  return std::regex_search(line, tokenRegex);
}

// Scans every *.h file directly under include/ITSMFTTracking/ (never
// recursing into detail/) and fails if any non-allowlisted file's text
// mentions `token`. Returns the number of files scanned, so callers can
// assert it against a sanity floor.
size_t scanNonDetailHeadersForToken(const fs::path& includeDir, const std::string& token,
                                    const std::set<std::string>& allowlist)
{
  const std::regex tokenRegex{"\\b" + token + "\\b"};
  size_t scanned = 0;
  for (const auto& entry : fs::directory_iterator(includeDir)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".h") {
      continue;
    }
    const auto name = entry.path().filename().string();
    ++scanned;
    if (allowlist.count(name) != 0) {
      continue;
    }
    std::ifstream input{entry.path()};
    BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << entry.path().string());
    std::string line;
    size_t lineNumber = 0;
    while (std::getline(input, line)) {
      ++lineNumber;
      BOOST_CHECK_MESSAGE(!lineMentionsToken(line, tokenRegex),
                          name << ":" << lineNumber << " (non-detail, non-allowlisted) mentions " << token << ": " << line);
    }
  }
  return scanned;
}

} // namespace

BOOST_AUTO_TEST_CASE(NoPublicNonDetailHeaderExposesTheLegacyTransitionPolicyTag)
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  const fs::path includeDir = fs::path(testDirectory) / ".." / "include" / "ITSMFTTracking";
  BOOST_REQUIRE_MESSAGE(fs::is_directory(includeDir), "cannot find " << includeDir.string());

  // detail/ is expected to still declare/define TransitionPolicyTag and its
  // dispatch/grouping/tag-keyed-traits machinery -- that is exactly what
  // "confined to detail/" means -- so it is deliberately not scanned here.
  const fs::path detailDir = includeDir / "detail";
  BOOST_REQUIRE_MESSAGE(fs::is_directory(detailDir), "detail/ boundary directory is missing: " << detailDir.string());

  const auto scanned = scanNonDetailHeadersForToken(includeDir, "TransitionPolicyTag", transitionPolicyTagAllowlist());
  // Sanity: this must actually have scanned a realistic number of headers --
  // an empty/near-empty loop (e.g. from a resolved-but-wrong includeDir)
  // would make every BOOST_CHECK_MESSAGE above vacuously pass.
  BOOST_CHECK_GT(scanned, 40u);
}

BOOST_AUTO_TEST_CASE(NoPublicTopologyOrSchedulingHeaderExposesStateFamilyAsATagSubstitute)
{
  // M4b: StateFamily leaked into TrackerTraits::getPolicyBindingCount(StateFamily)
  // (a test-only introspection seam, since deleted) and into Cell.h's
  // stateFamilyFromNLayers<NLayers>() (a legacy NLayers-to-family dispatch
  // bridge, since removed with the temporary parameter layer). Neither
  // public-facing use derived anything from SurfaceDescriptor/SurfaceKind;
  // both used StateFamily as a stand-in transition-policy key. This proves
  // neither -- nor any future equivalent -- reappears outside the narrow,
  // reasoned state-representation allowlist above.
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  const fs::path includeDir = fs::path(testDirectory) / ".." / "include" / "ITSMFTTracking";
  BOOST_REQUIRE_MESSAGE(fs::is_directory(includeDir), "cannot find " << includeDir.string());

  const auto scanned = scanNonDetailHeadersForToken(includeDir, "StateFamily", stateFamilyAllowlist());
  BOOST_CHECK_GT(scanned, 40u);
}

BOOST_AUTO_TEST_CASE(TrackerTraitsNoLongerDeclaresAFamilyOrTagKeyedBindingCountAccessor)
{
  // Explicit deletion proof for getPolicyBindingCount(...) (M4b): rather than
  // relying solely on the (necessarily coarser) StateFamily token scan above,
  // grep TrackerTraits.h directly for the removed accessor's name -- guards
  // against a future re-introduction under either a TransitionPolicyTag or a
  // StateFamily parameter.
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  const fs::path header = fs::path(testDirectory) / ".." / "include" / "ITSMFTTracking" / "TrackerTraits.h";
  std::ifstream input{header};
  BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << header.string());
  std::string line;
  size_t lineNumber = 0;
  while (std::getline(input, line)) {
    ++lineNumber;
    BOOST_CHECK_MESSAGE(line.find("getPolicyBindingCount") == std::string::npos,
                        "TrackerTraits.h:" << lineNumber << " still declares getPolicyBindingCount: " << line);
  }
}

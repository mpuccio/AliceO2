// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// GenericTrackingEngineMigration.md M1 dependency-boundary proof (ADR 0007
// decisions 2, 5, 7): the public engine contract
// (ITSMFTTracking/TrackingEngine.h, TrackingParticipant.h, ParticipantId.h)
// must not depend on ITS/MFT tracking implementation headers, DPL/workflow
// headers, output-writer headers, or TransitionPolicyTag.
//
// Two independent halves, same convention as
// testForwardSurfaceStateOperations.cxx's
// NewProductionFilesHaveNoLegacyForwardDependency:
//  - compile half (this file's own include list is restricted to
//    ITSMFTTracking/TrackingEngine.h plus ordinary standard headers --
//    nothing else project-related -- and still names every public
//    engine-contract type, including TrackingEngine::executeEvent()'s
//    exact member-function type, which mentions TimeFrame only as an
//    incomplete forward-declared type);
//  - grep half (reads the three production headers themselves and asserts
//    none contains any forbidden token).

#define BOOST_TEST_MODULE ITSMFT TrackingEngineDependencyBoundary
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>

#include "ITSMFTTracking/TrackingEngine.h"

using namespace o2::itsmft::tracking;

// --- Compile half ------------------------------------------------------

BOOST_AUTO_TEST_CASE(EngineHeaderAloneNamesTheCompleteEngineContract)
{
  // TrackingEngine.h (via TrackingParticipant.h) only forward-declares
  // TimeFrame: naming TrackingEngine::executeEvent()'s member-function type
  // does not require a complete TimeFrame, so this line compiles from this
  // file's include list alone.
  using ExecuteEventType = EventResult (TrackingEngine::*)(TimeFrame&, gsl::span<TrackingParticipant* const>, const o2::InteractionRecord&);
  static_assert(std::is_same_v<decltype(&TrackingEngine::executeEvent), ExecuteEventType>);

  static_assert(std::is_default_constructible_v<TrackingEngine>);
  static_assert(std::is_default_constructible_v<ParticipantId>);
  static_assert(std::is_default_constructible_v<EventResult>);
  static_assert(std::is_default_constructible_v<ParticipantEventResult>);
  static_assert(std::is_default_constructible_v<ParticipantLoadResult>);
  static_assert(std::is_default_constructible_v<ParticipantTrackingResult>);
  static_assert(std::is_default_constructible_v<ParticipantPublicationExport>);
  static_assert(std::is_abstract_v<TrackingParticipant>);
  static_assert(std::has_virtual_destructor_v<TrackingParticipant>);
}

// --- Grep half -----------------------------------------------------------

BOOST_AUTO_TEST_CASE(EngineAndParticipantHeadersExcludeForbiddenDependencies)
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  const std::array<std::string, 3> productionHeaders = {
    testDirectory + "/../include/ITSMFTTracking/TrackingEngine.h",
    testDirectory + "/../include/ITSMFTTracking/TrackingParticipant.h",
    testDirectory + "/../include/ITSMFTTracking/ParticipantId.h"};

  // ITS/MFT tracking implementation headers:
  const std::vector<std::string> forbidden = {
    "CATracker.h", "LegacyTrackerScratch.h", "TrackerTraits.h",
    // TransitionPolicyTag machinery (ADR 0007 decision 7):
    "TransitionPolicyDispatch.h", "TransitionPolicyOperations.h",
    "TransitionPolicyBinding.h", "TransitionPolicyState.h", "TransitionPolicyTag",
    // Detector identity / current-output-type / workflow-facing headers:
    "DetID.h", "TrackingInterface.h", "CombinedTimeFrameCoordinator.h",
    "CommonTrackOutputAdapter.h",
    // DPL/workflow headers:
    "DataProcessorSpec.h", "WorkflowSpec.h", "ConfigParamSpec.h", "DeviceSpec.h",
    // Output-writer headers:
    "RootTreeWriter.h", "MakeRootTreeWriterSpec.h", "Writer.h"};

  // Scans only #include directive lines, not the surrounding prose -- this
  // file's own doc comments name several of these forbidden headers
  // precisely to explain why they are *not* included, so a whole-file
  // substring scan would misfire on that explanation.
  for (const auto& file : productionHeaders) {
    std::ifstream input{file};
    BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << file);
    std::string line;
    while (std::getline(input, line)) {
      const auto hashPos = line.find_first_not_of(" \t");
      if (hashPos == std::string::npos || line[hashPos] != '#') {
        continue;
      }
      if (line.find("include", hashPos + 1) == std::string::npos) {
        continue;
      }
      for (const auto& token : forbidden) {
        BOOST_CHECK_MESSAGE(line.find(token) == std::string::npos, file << " includes forbidden header via: " << line);
      }
    }
  }
}

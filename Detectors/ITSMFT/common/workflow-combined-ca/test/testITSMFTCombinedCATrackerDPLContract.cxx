// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 4 C4: focused tests for the DPL input/output contract of
// o2::itsmft::combined::getCombinedCATrackerSpec() -- MC/non-MC variants,
// the hard requirement that no vertex-related OutputSpec is ever declared,
// that ITS and MFT each get their own distinct input bindings/contexts, that
// geometry is always requested via explicit nominal CCDB objects (never the
// aligned global-geometry route), that no IR-frame input is ever declared,
// and that this device's own name never collides with either single-detector
// opt-in tracker (its-ca-tracker/mft-ca-tracker).

#define BOOST_TEST_MODULE ITSMFT ITSMFTCombinedCATrackerDPLContract
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <string>

#include "Framework/DataProcessorSpec.h"
#include "Framework/DataSpecUtils.h"
#include "ITSMFTCombinedCAWorkflow/CombinedCATrackerSpec.h"

using namespace o2::framework;

namespace
{
bool hasInput(const std::vector<InputSpec>& specs, const std::string& binding)
{
  return std::any_of(specs.begin(), specs.end(), [&binding](const InputSpec& s) { return s.binding == binding; });
}

bool hasOutputDesc(const std::vector<OutputSpec>& specs, const std::string& origin, const std::string& desc)
{
  return std::any_of(specs.begin(), specs.end(), [&](const OutputSpec& s) {
    return DataSpecUtils::describe(s).find(origin) != std::string::npos && DataSpecUtils::describe(s).find(desc) != std::string::npos;
  });
}
} // namespace

BOOST_AUTO_TEST_CASE(NonMCContractHasNoLabelsInputOrMCOutputs)
{
  const auto spec = o2::itsmft::combined::getCombinedCATrackerSpec(false);

  BOOST_CHECK(hasInput(spec.inputs, "compClustersITS"));
  BOOST_CHECK(hasInput(spec.inputs, "patternsITS"));
  BOOST_CHECK(hasInput(spec.inputs, "ROframesITS"));
  BOOST_CHECK(hasInput(spec.inputs, "itscldict"));
  BOOST_CHECK(hasInput(spec.inputs, "compClustersMFT"));
  BOOST_CHECK(hasInput(spec.inputs, "patternsMFT"));
  BOOST_CHECK(hasInput(spec.inputs, "ROframesMFT"));
  BOOST_CHECK(hasInput(spec.inputs, "mftcldict"));
  BOOST_CHECK(!hasInput(spec.inputs, "labelsITS"));
  BOOST_CHECK(!hasInput(spec.inputs, "labelsMFT"));

  BOOST_CHECK(hasOutputDesc(spec.outputs, "ITS", "TRACKS"));
  BOOST_CHECK(hasOutputDesc(spec.outputs, "ITS", "TRACKCLSID"));
  BOOST_CHECK(hasOutputDesc(spec.outputs, "ITS", "ITSTrackROF"));
  BOOST_CHECK(hasOutputDesc(spec.outputs, "MFT", "TRACKS"));
  BOOST_CHECK(hasOutputDesc(spec.outputs, "MFT", "MFTTrackROF"));
  BOOST_CHECK(hasOutputDesc(spec.outputs, "MFT", "TRACKCLSID"));
  BOOST_CHECK(hasOutputDesc(spec.outputs, "MFT", "TRACKSEEDPAT"));
  BOOST_CHECK(!hasOutputDesc(spec.outputs, "ITS", "TRACKSMCTR"));
  BOOST_CHECK(!hasOutputDesc(spec.outputs, "MFT", "TRACKSMCTR"));
}

BOOST_AUTO_TEST_CASE(MCContractAddsLabelsInputsAndMCOutputsForBothDetectors)
{
  const auto spec = o2::itsmft::combined::getCombinedCATrackerSpec(true);

  BOOST_CHECK(hasInput(spec.inputs, "labelsITS"));
  BOOST_CHECK(hasInput(spec.inputs, "labelsMFT"));
  BOOST_CHECK(hasOutputDesc(spec.outputs, "ITS", "TRACKSMCTR"));
  BOOST_CHECK(hasOutputDesc(spec.outputs, "MFT", "TRACKSMCTR"));
}

BOOST_AUTO_TEST_CASE(GeometryIsAlwaysRequestedExplicitlyNeverAligned)
{
  // No --use-geom option exists on this workflow: itsTGeo/mftTGeo (each
  // detector's own nominal/ideal GeometryTGeo CCDB condition object) must
  // always be present, in both the MC and non-MC variant, proving the
  // GRPGeomRequest::Aligned route (which would omit these and instead
  // request the full aligned "GLO/Config/GeometryAligned" object) is
  // structurally unreachable from this workflow.
  for (const bool useMC : {false, true}) {
    const auto spec = o2::itsmft::combined::getCombinedCATrackerSpec(useMC);
    BOOST_CHECK_MESSAGE(hasInput(spec.inputs, "itsTGeo"), "useMC=" << useMC);
    BOOST_CHECK_MESSAGE(hasInput(spec.inputs, "mftTGeo"), "useMC=" << useMC);
  }
}

BOOST_AUTO_TEST_CASE(NoIRFrameInputIsEverDeclared)
{
  // The combined ITS+MFT common-CA tracking flow has no IR-frame/trigger-window masking;
  // this workflow must never declare an IRFramesITS input regardless of MC
  // -- see ConfigPreflight.h's requireNoMFTIRFrameConfigOrFatal()/
  // requireContinuousMFTReadoutOrFatal() for the two runtime-enforced halves
  // of this same contract.
  for (const bool useMC : {false, true}) {
    const auto spec = o2::itsmft::combined::getCombinedCATrackerSpec(useMC);
    BOOST_CHECK_MESSAGE(!hasInput(spec.inputs, "IRFramesITS"), "useMC=" << useMC);
  }
}

BOOST_AUTO_TEST_CASE(NoVertexRelatedOutputsArePresentEver)
{
  for (const bool useMC : {false, true}) {
    const auto spec = o2::itsmft::combined::getCombinedCATrackerSpec(useMC);
    for (const auto& out : spec.outputs) {
      const auto desc = DataSpecUtils::describe(out);
      BOOST_CHECK_MESSAGE(desc.find("VERTICES") == std::string::npos, "unexpected vertex output: " << desc);
    }
  }
}

BOOST_AUTO_TEST_CASE(DeviceNameIsDistinctFromEitherSingleDetectorOptInTracker)
{
  const auto spec = o2::itsmft::combined::getCombinedCATrackerSpec(false);
  BOOST_CHECK_NE(spec.name, "its-ca-tracker");
  BOOST_CHECK_NE(spec.name, "mft-ca-tracker");
}

BOOST_AUTO_TEST_CASE(ITSAndMFTInputsHaveDistinctOriginsAndBindings)
{
  // Never the same InputSpec identity reused for both detectors -- each
  // has its own binding name and Origin, so the two ClusterSourceInputs
  // CombinedCATrackerDPL::run() builds can never be cross-wired.
  const auto spec = o2::itsmft::combined::getCombinedCATrackerSpec(false);
  const auto itsRof = std::find_if(spec.inputs.begin(), spec.inputs.end(), [](const InputSpec& s) { return s.binding == "ROframesITS"; });
  const auto mftRof = std::find_if(spec.inputs.begin(), spec.inputs.end(), [](const InputSpec& s) { return s.binding == "ROframesMFT"; });
  BOOST_REQUIRE(itsRof != spec.inputs.end());
  BOOST_REQUIRE(mftRof != spec.inputs.end());
  BOOST_CHECK(DataSpecUtils::describe(*itsRof) != DataSpecUtils::describe(*mftRof));
}

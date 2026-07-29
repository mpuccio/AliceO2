// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 3 ITS common-CA config-ownership correction: IOUtils.cxx's
// per-cluster systematic-error correction must not consult
// TrackerParamRef<ITS>::get() (o2::its::TrackerParamConfig, the frozen
// legacy "ITSCATrackerParam" namespace)'s sysErrY2/sysErrZ2 fields. The
// dedicated ITSCommonCATrackerParam configuration does not support a
// systematic-error feature at all, so common ITS normalized decoding must
// treat it as an explicit no-op. MFT is unchanged: it keeps reading its own
// live TrackerParamConfig<MFT> ("MFTCATrackerParam") sysErr2Row/sysErr2Col
// exactly as before.
//
// o2::itsmft::ioutils::detail::shouldApplySysErrors<DetId>()/addSysErrors<DetId>()
// (IOUtils.h) are the smallest seam that lets this be tested directly,
// without a geometry/dictionary/CompClusterExt fixture: they are the exact
// functions IOUtils.cxx's decodeCluster()/decodeClusterBounded()/
// convertCompactClusters() call, moved from an anonymous namespace in
// IOUtils.cxx into IOUtils.h's pre-existing `detail` namespace (which
// already held isSensorInGeometry()/isLayerInDetector() for the same
// reason: internal, not broad public API).

#define BOOST_TEST_MODULE ITSMFT IOUtilsSysErrorOwnership
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITStracking/TrackingConfigParam.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;
using namespace o2::itsmft::ioutils::detail;

namespace
{

o2::its::TrackerParamConfig& mutableLegacyITSTrackerParamConfig()
{
  return const_cast<o2::its::TrackerParamConfig&>(o2::its::TrackerParamConfig::Instance());
}

TrackerParamConfig<o2::detectors::DetID::MFT>& mutableMFTTrackerParamConfig()
{
  return const_cast<TrackerParamConfig<o2::detectors::DetID::MFT>&>(
    TrackerParamConfig<o2::detectors::DetID::MFT>::Instance());
}

// Restores every singleton field this test file stages, on scope exit
// (including a thrown/failed check), so no test case here can leak state
// into another test case in this binary.
struct ScopedSysErrParams {
  ScopedSysErrParams()
  {
    auto& itsLegacy = mutableLegacyITSTrackerParamConfig();
    for (int i = 0; i < ITSNLayers; ++i) {
      originalLegacyITSSysErrY2[i] = itsLegacy.sysErrY2[i];
      originalLegacyITSSysErrZ2[i] = itsLegacy.sysErrZ2[i];
    }
    auto& mft = mutableMFTTrackerParamConfig();
    for (int i = 0; i < MFTNLayers; ++i) {
      originalMFTSysErr2Row[i] = mft.sysErr2Row[i];
      originalMFTSysErr2Col[i] = mft.sysErr2Col[i];
    }
  }
  ~ScopedSysErrParams()
  {
    auto& itsLegacy = mutableLegacyITSTrackerParamConfig();
    for (int i = 0; i < ITSNLayers; ++i) {
      itsLegacy.sysErrY2[i] = originalLegacyITSSysErrY2[i];
      itsLegacy.sysErrZ2[i] = originalLegacyITSSysErrZ2[i];
    }
    auto& mft = mutableMFTTrackerParamConfig();
    for (int i = 0; i < MFTNLayers; ++i) {
      mft.sysErr2Row[i] = originalMFTSysErr2Row[i];
      mft.sysErr2Col[i] = originalMFTSysErr2Col[i];
    }
  }
  float originalLegacyITSSysErrY2[ITSNLayers];
  float originalLegacyITSSysErrZ2[ITSNLayers];
  float originalMFTSysErr2Row[MFTNLayers];
  float originalMFTSysErr2Col[MFTNLayers];
};

} // namespace

// --- ITS: a legacy sysErrY2/sysErrZ2 override must not reach common ITS
// normalized decoding: it is an explicit no-op, not merely a config with no
// effect. ---

BOOST_FIXTURE_TEST_CASE(ITSDoesNotApplyLegacySysErrorsByDefault, ScopedSysErrParams)
{
  BOOST_CHECK_EQUAL(shouldApplySysErrors<o2::detectors::DetID::ITS>(), false);
}

BOOST_FIXTURE_TEST_CASE(ITSIgnoresLegacySysErrY2AndSysErrZ2, ScopedSysErrParams)
{
  auto& itsLegacy = mutableLegacyITSTrackerParamConfig();
  for (int i = 0; i < ITSNLayers; ++i) {
    itsLegacy.sysErrY2[i] = 5.e-4f;
    itsLegacy.sysErrZ2[i] = 7.e-4f;
  }

  // If the legacy field still leaked in, this would be true, and the
  // addSysErrors() call below would perturb sigma2Row/sigma2Col.
  BOOST_CHECK_EQUAL(shouldApplySysErrors<o2::detectors::DetID::ITS>(), false);

  float sigma2Row = 1.5f;
  float sigma2Col = 2.5f;
  addSysErrors<o2::detectors::DetID::ITS>(0, sigma2Row, sigma2Col);
  BOOST_CHECK_CLOSE(sigma2Row, 1.5f, 1e-4);
  BOOST_CHECK_CLOSE(sigma2Col, 2.5f, 1e-4);
}

// --- MFT: no regression. TrackerParamConfig<MFT>::sysErr2Row/sysErr2Col
// ("MFTCATrackerParam", MFT's own live configuration, not a legacy
// namespace) must remain authoritative for MFT exactly as before this
// correction. ---

BOOST_FIXTURE_TEST_CASE(MFTDoesNotApplySysErrorsByDefault, ScopedSysErrParams)
{
  BOOST_CHECK_EQUAL(shouldApplySysErrors<o2::detectors::DetID::MFT>(), false);
}

BOOST_FIXTURE_TEST_CASE(MFTStillAppliesItsOwnLiveSysErrorConfig, ScopedSysErrParams)
{
  auto& mft = mutableMFTTrackerParamConfig();
  mft.sysErr2Row[3] = 5.e-4f;
  mft.sysErr2Col[3] = 7.e-4f;

  BOOST_CHECK_EQUAL(shouldApplySysErrors<o2::detectors::DetID::MFT>(), true);

  float sigma2Row = 1.5f;
  float sigma2Col = 2.5f;
  addSysErrors<o2::detectors::DetID::MFT>(3, sigma2Row, sigma2Col);
  BOOST_CHECK_CLOSE(sigma2Row, 1.5f + 5.e-4f, 1e-4);
  BOOST_CHECK_CLOSE(sigma2Col, 2.5f + 7.e-4f, 1e-4);
}

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

/// \file MeasurementSinkSpec.h
/// \brief The study's sink: the tracker over the producer's measurements, with
///        the truth reports, the QA trees and the comparison with O2's tracks.

#ifndef GTRK_STUDY_MEASUREMENTSINKSPEC_H_
#define GTRK_STUDY_MEASUREMENTSINKSPEC_H_

#include "GenericTracking/Base/LayerSpec.h"

#include "Framework/DataProcessorSpec.h"

#include <cstdint>

namespace o2::gtrk::study
{

/// The sink for \a Det. Its tracker seeds each iteration from the
/// source GTRKReco<DET>.iterationVertexSource configures, the beam line by default.
/// \a withVertices adds the vertex list input, which a vertex-list iteration
/// needs. With \a products other than kNone it also takes the tracking
/// device's products (workflow/TrackerSpec.h), every timeframe. Defined for
/// ITS, MFT and TRKFT3, each by its own compilation of MeasurementSinkSpec.cxx.
enum class TrackerProducts : std::uint8_t {
  kNone,
  kCompare, ///< --compare-tracks: against its own tracker's, every track record field by field, every vertex and contributor count; the run fails on a difference
  kDigests, ///< --tracker-digests: their digests logged, as the bench prints a device's
};
template <DetectorLayout Det>
o2::framework::DataProcessorSpec getMeasurementSinkSpec(bool useMC, bool withVertices, TrackerProducts products);

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_MEASUREMENTSINKSPEC_H_

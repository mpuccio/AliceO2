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

/// \file TrackerSpec.h
/// \brief The tracking device, gtrk-<det>-tracker: the measurements in, the
///        tracks and the own vertexer's vertices out, in GenericTracking's
///        own formats (DataSpecs.h's kTracks, kVertexSeeds,
///        kVertexContributors). One device per detector, its own
///        instantiation.

#ifndef GTRK_WORKFLOW_TRACKERSPEC_H_
#define GTRK_WORKFLOW_TRACKERSPEC_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/TrackRefit.h"

#include "Framework/DataProcessorSpec.h"

namespace o2::gtrk::workflow
{

/// A track as the device publishes it: the tracker's own record, its refit's
/// `status` with it (kOk, or why the refit failed).
template <DetectorLayout Det>
using TrackRecord = typename TrackRefit<Det>::Tracks::Record;

/// The device; with \a withVertices the native vertex seeds (kVertices) are
/// an input.
template <DetectorLayout Det>
o2::framework::DataProcessorSpec getTrackerSpec(bool withVertices);

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_TRACKERSPEC_H_

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

/// \file BeamLine.h
/// \brief The run's beam line, from the MeanVertex calibration object, as
///        the producer and the trackers read it.
///
/// A calibration input, not a reconstruction parameter:
/// offline it is the run's `GLO/Calib/MeanVertex`. The producer subtracts its
/// x and y from every measurement's global position, so that r, phi and the
/// beam-frame variances the finding reads are about the corrected beam line,
/// as O2's ITS tracker does (TimeFrame.cxx, the index tables); the fit works
/// in the local frames and the detector frame, untouched. Its z and widths
/// are the trackers' VertexConstraint. Without an object (the TRK, whose beam
/// line is the nominal one: ConditionsTRK.h) the beam line is the nominal one:
/// the axis, 20 um wide, 6 cm long, as the MC's Diamond gives.

#ifndef GTRK_WORKFLOW_BEAMLINE_H_
#define GTRK_WORKFLOW_BEAMLINE_H_

#include "GenericTracking/Tracking/SearchWindow.h"

#include "DataFormatsCalibration/MeanVertexObject.h"

#include <algorithm>
#include <format>
#include <string>

namespace o2::gtrk::workflow
{

struct BeamLine {
  double x = 0.;         ///< [cm]
  double y = 0.;         ///< [cm]
  double z = 0.;         ///< [cm]
  double sigmaX = 0.002; ///< [cm]
  double sigmaY = 0.002; ///< [cm]
  double sigmaZ = 6.;    ///< [cm]

  bool operator==(const BeamLine&) const = default;
};

/// The nominal axis, as the MC's Diamond places the beam.
[[nodiscard]] inline BeamLine nominalBeamLine() { return {}; }

/// The run's, from its MeanVertex. The slopes are left out, as O2 leaves them:
/// 2e-4 over the ITS's half length is 30 um, below the width.
[[nodiscard]] inline BeamLine beamLineOf(const o2::dataformats::MeanVertexObject& object) { return {.x = object.getX(), .y = object.getY(), .z = object.getZ(), .sigmaX = object.getSigmaX(), .sigmaY = object.getSigmaY(), .sigmaZ = object.getSigmaZ()}; }

/// The trackers' constraint: z and its width, and across the beam the larger
/// of the two widths, since the pair test's beam term is a bound. The windows
/// and cuts admit \a windowNSigma of the constraint's z width; that width is
/// the beam's taken \a zReachNSigma / \a windowNSigma times, so that they
/// reach \a zReachNSigma of the beam's.
[[nodiscard]] inline VertexConstraint beamSpotConstraint(const BeamLine& beam, double zReachNSigma, double windowNSigma)
{
  return {
    .z = beam.z,
    .sigmaZ = beam.sigmaZ * zReachNSigma / windowNSigma,
    .sigmaTransverse = std::max(beam.sigmaX, beam.sigmaY),
    .x = beam.x,
    .y = beam.y,
  };
}

[[nodiscard]] inline std::string describe(const BeamLine& beam) { return std::format("({:.4f}, {:.4f}, {:.4f}) cm, widths ({:.4f}, {:.4f}, {:.3f}) cm", beam.x, beam.y, beam.z, beam.sigmaX, beam.sigmaY, beam.sigmaZ); }

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_BEAMLINE_H_

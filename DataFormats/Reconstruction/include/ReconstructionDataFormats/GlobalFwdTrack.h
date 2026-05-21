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

/// \file GlobalFwdTrack.h
/// \brief Global Forward Muon tracks

#ifndef ALICEO2_TRACKGLOBALFWD_H
#define ALICEO2_TRACKGLOBALFWD_H

#include <array>
#include "ReconstructionDataFormats/TrackFwd.h"
#include "ReconstructionDataFormats/MatchInfoFwd.h"

namespace o2
{
namespace dataformats
{
using FwdResiduals = std::array<double, 5>;

class GlobalFwdTrack : public o2::track::TrackParCovFwd, public o2::dataformats::MatchInfoFwd
{
 public:
  GlobalFwdTrack() = default;
  GlobalFwdTrack(const GlobalFwdTrack& t) = default;
  GlobalFwdTrack(o2::track::TrackParCovFwd const& t) { *this = t; }
  ~GlobalFwdTrack() = default;

  FwdResiduals computeResiduals2Cov(const o2::track::TrackParCovFwd& t) const
  {
    FwdResiduals Residuals2Cov;

    Residuals2Cov[0] = (getX() - t.getX()) / o2::math_utils::sqrtd(getSigma2X() + t.getSigma2X());
    Residuals2Cov[1] = (getY() - t.getY()) / o2::math_utils::sqrtd(getSigma2Y() + t.getSigma2Y());
    Residuals2Cov[2] = (getPhi() - t.getPhi()) / o2::math_utils::sqrtd(getSigma2Phi() + t.getSigma2Phi());
    Residuals2Cov[3] = (getTanl() - t.getTanl()) / o2::math_utils::sqrtd(getSigma2Tanl() + t.getSigma2Tanl());
    Residuals2Cov[4] = (getInvQPt() - t.getInvQPt()) / o2::math_utils::sqrtd(getSigma2InvQPt() + t.getSigma2InvQPt());
    return Residuals2Cov;
  }

 private:
  ClassDefNV(GlobalFwdTrack, 2);
};

} // namespace dataformats

namespace framework
{
template <typename T>
struct is_messageable;
template <>
struct is_messageable<o2::dataformats::GlobalFwdTrack> : std::true_type {
};
} // namespace framework

} // namespace o2

#endif

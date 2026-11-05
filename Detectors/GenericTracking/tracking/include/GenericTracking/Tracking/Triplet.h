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

/// \file Triplet.h
/// \brief A triplet as SoA columns: two tracklets, their local curvature, and
///        what the neighbour fit reads of them.
///
/// Adapted from O2 ITSMFTTracking/Triplet.h. GTRK stores only the fields read
/// by later CA stages; measurement rows and layers follow from the tracklets.

#ifndef GTRK_TRACKING_TRIPLET_H_
#define GTRK_TRACKING_TRIPLET_H_

#include "GenericTracking/Tracking/Definitions.h"

#include <memlayout/wrapper.h>

#include <cstdint>

namespace o2::gtrk
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
template <template <class> class C>
struct TripletT {
  MEMLAYOUT_APPLY(TripletT, first, second, curvature, varCurvature, level, kinked)
  C<std::uint32_t> first, second; ///< indices of tracklets a -> b and b -> c
  C<Scalar> curvature;            ///< fitted signed 3D curvature, prop. to q/p [1/cm]; transverse is curvature / sin(theta)
  C<Scalar> varCurvature;         ///< curvature variance [1/cm^2]
  C<std::uint32_t>
    level;                ///< CA level: 1 from the triplet stage, set by the neighbour stage stratum by stratum; 32 bits for the tests' fixed-point sweep, an atomic max
  C<std::uint8_t> kinked; ///< 1 if it passed only with the iteration's kink allowance
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)

using Triplet = TripletT<memlayout::value>;
using TripletConstView = memlayout::Wrapper<TripletT, memlayout::const_pointer, memlayout::Layout::soa>;

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRIPLET_H_

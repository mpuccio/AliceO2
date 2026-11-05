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

/// \file Neighbour.h
/// \brief A neighbour as SoA columns: the predecessor of the triplet that
///        emitted it.
///
/// The neighbour stage runs one element per triplet `T2` and emits the
/// triplets `T1` it continues, so the target is the element and only the
/// predecessor is stored. The road walk
/// (`RoadCount.h`) reads exactly this: each triplet's stored predecessors.

#ifndef GTRK_TRACKING_NEIGHBOUR_H_
#define GTRK_TRACKING_NEIGHBOUR_H_

#include <memlayout/wrapper.h>

#include <cstdint>

namespace o2::gtrk
{

// As for TrackletT: upstream's `apply` macro also uses the member names as
// enum constants.
// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
template <template <class> class C>
struct NeighbourT {
  MEMLAYOUT_APPLY(NeighbourT, from)
  C<std::uint32_t> from; ///< index of the predecessor triplet `T1`, whose second tracklet is this triplet's first
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)

using Neighbour = NeighbourT<memlayout::value>;

} // namespace o2::gtrk

#endif // GTRK_TRACKING_NEIGHBOUR_H_

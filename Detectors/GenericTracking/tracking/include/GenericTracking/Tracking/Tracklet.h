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

/// \file Tracklet.h
/// \brief A tracklet as SoA columns: two measurement rows, their link and the
///        chord's direction.
///
/// A tracklet joins a *source* -- the inner measurement whose position and
/// time open the grid query -- to a *target*, the outer measurement the query
/// found and the pair cut accepted. Both are global measurement rows of the
/// timeframe; `Measurement::clusterIndex` leads from a row back to its cluster.
///
/// `link` indexes `Det::layerLinks`. The rows imply it, but only through a
/// search over the layer boundaries; the triplet stage needs it for every
/// combination it forms, because one source's range mixes all of its layer's
/// outgoing links and the tolerances and scattering terms are per link.
///
/// `lambda` and `phi` are the direction of the chord from source to target.
/// They are computed once per accepted tracklet and cached for the triplet
/// stage's pre-cut: a middle row's outgoing tracklets are contiguous, so
/// comparing cached directions is a contiguous read where each outer position
/// would be a gather from the measurement columns.
///
/// The dip angle is stored rather than tan(lambda). Between two discs a
/// forward track's transverse chord is short while dz is the disc spacing, so
/// dz over the chord grows without bound and is dominated by the position
/// resolution; atan2(dz, chord) stays in [-pi/2, pi/2] and well conditioned
/// for every kind of link. The chord's azimuth has no such cure: its error
/// grows as the r-phi resolution over the chord, which the triplet stage's
/// tolerance has to carry.

#ifndef GTRK_TRACKING_TRACKLET_H_
#define GTRK_TRACKING_TRACKLET_H_

#include "GenericTracking/Tracking/Definitions.h"

#include <memlayout/wrapper.h>

#include <cstdint>

namespace o2::gtrk
{

// As for MeasurementT: an SoA element proxy holds references, and `apply` is
// upstream's macro, which also uses the camelBack member names as enum
// constants.
// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
template <template <class> class C>
struct TrackletT {
  MEMLAYOUT_APPLY(TrackletT, source, target, link, lambda, phi)
  C<std::uint32_t> source, target; ///< global measurement rows
  C<std::uint16_t> link;           ///< index into Det::layerLinks
  C<Scalar> lambda;                ///< atan2(dz, transverse chord), source to target, in [-pi/2, pi/2]
  C<Scalar> phi;                   ///< atan2(dy, dx) of the chord, source to target, in [-pi, pi]
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)

/// One tracklet record, used by emitters and host-side validation.
using Tracklet = TrackletT<memlayout::value>;

/// The columns in a stage's arena, as the tracklet stage writes them.
using TrackletView = memlayout::Wrapper<TrackletT, memlayout::pointer, memlayout::Layout::soa>;

/// The same columns, as every later stage reads them.
using TrackletConstView = memlayout::Wrapper<TrackletT, memlayout::const_pointer, memlayout::Layout::soa>;

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRACKLET_H_

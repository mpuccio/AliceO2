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

/// \file ProfileReport.h
/// \brief The host profile's reports: a profile (Tracking/Profile.h) as an
///        indented tree for people and as JSON for the scripts, and one
///        iteration's summary line, the sink's and the bench's.

#ifndef GTRK_STUDY_PROFILEREPORT_H_
#define GTRK_STUDY_PROFILEREPORT_H_

#include "GenericTracking/Tracking/Profile.h"

#include <cstddef>
#include <format>
#include <ostream>
#include <string>
#include <vector>

namespace o2::gtrk::study
{

/// A profile as an indented tree: each phase with its share of its parent,
/// and an `(other)` line for the time its children do not cover.
void report(std::ostream& out, const profile::Profile& profile, profile::Phase root, int threads);

/// A profile as one JSON object, phase name to stats, for the scripts.
void writeJson(std::ostream& out, const profile::Profile& profile);

/// One iteration's summary (a tracker's IterationSummary) as the lines the
/// sink logs and the bench prints: the counts and
/// the phases' times, and the earlier iterations' tracks followed when there
/// are any.
template <typename Summary>
[[nodiscard]] std::vector<std::string> describeIteration(const Summary& summary, std::size_t iteration)
{
  std::vector<std::string> lines;
  lines.push_back(
    std::format("  iteration {}: {} slices, halo {} ps; {} roads, {} records ({:.1f} MB), {} duplicates, {} mismatches, {} rounds, {} tracks; "
                "{} refit failed; {} followed (+{} outward, +{} inward), {} dropped; {:.3f} s slices, {:.3f} s roads, {:.3f} s selection, {:.3f} s refit, "
                "{:.3f} s follow",
                iteration, summary.slices, summary.haloPS, summary.roads, summary.records, static_cast<double>(summary.recordBytes) / 1.e6, summary.duplicates, summary.mismatches, summary.rounds, summary.tracks, summary.refitFailures, summary.tracksFollowed, summary.clustersFollowedOutward, summary.clustersFollowedInward, summary.tracksDropped, summary.seconds(profile::Phase::kSlices), summary.seconds(profile::Phase::kRoads), summary.seconds(profile::Phase::kSelection), summary.seconds(profile::Phase::kRefit), summary.seconds(profile::Phase::kFollow)));
  if (summary.tracksFollowedEarlier != 0U || summary.seconds(profile::Phase::kFollowEarlier) > 0.) {
    lines.push_back(std::format("    earlier iterations' tracks followed: {} (+{} outward, +{} inward) in {:.3f} s", summary.tracksFollowedEarlier, summary.clustersFollowedEarlierOutward, summary.clustersFollowedEarlierInward, summary.seconds(profile::Phase::kFollowEarlier)));
  }
  return lines;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_PROFILEREPORT_H_

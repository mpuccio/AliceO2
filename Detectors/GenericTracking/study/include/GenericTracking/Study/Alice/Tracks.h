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

/// \file Tracks.h
/// \brief O2's ITS and MFT tracks, as the clusters they hold, their MC labels'
///        fake flags, fitted kinematics and ITS innermost-cluster states.
///
/// All a comparison with our tracks needs, because both trackers read the
/// same clusters in the same order (Clusters.h): a track is its list of O2
/// cluster indices, which are our measurements' `clusterIndex`.

#ifndef GTRK_STUDY_ALICE_TRACKS_H_
#define GTRK_STUDY_ALICE_TRACKS_H_

#include "GenericTracking/Study/Alice/Hits.h"

#include "ReconstructionDataFormats/Track.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace o2::gtrk::study
{

struct TrackKinematics;

/// One timeframe's tracks, each as its cluster indices.
using TrackClusters = std::vector<std::vector<std::uint32_t>>;

/// Reads `o2trac_its.root` or `mfttracks.root`, one tree entry per timeframe.
class TrackFile
{
 public:
  /// Throws std::runtime_error if \a path has no o2sim tree.
  TrackFile(const std::string& path, Chips chips);
  ~TrackFile();
  TrackFile(const TrackFile&) = delete;
  TrackFile& operator=(const TrackFile&) = delete;
  TrackFile(TrackFile&&) = delete;
  TrackFile& operator=(TrackFile&&) = delete;

  /// The tracks of entry \a entry; empty past the last. Throws
  /// std::runtime_error if a track points past the cluster index array.
  [[nodiscard]] TrackClusters read(long long entry);

  /// Whether each track's own MC label is flagged fake, as O2's tracker set it.
  [[nodiscard]] std::vector<bool> fakeFlags(long long entry);

  /// Fitted inward kinematics in the same order as read(). Native outputs
  /// contain accepted tracks; fitValid checks finite, positive momentum.
  [[nodiscard]] std::vector<TrackKinematics> kinematics(long long entry);

  /// Every ITS track's state at its innermost cluster. Throws
  /// std::runtime_error for the MFT, whose tracks are parameterised in z.
  [[nodiscard]] std::vector<o2::track::TrackParCov> innerStates(long long entry);

 private:
  Chips mChips;
  std::unique_ptr<TFile> mFile;
  TTree* mTree = nullptr;
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_ALICE_TRACKS_H_

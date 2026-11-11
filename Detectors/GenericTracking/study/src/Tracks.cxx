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

#include "GenericTracking/Study/Alice/Tracks.h"
#include "GenericTracking/Study/TrackingQA.h"

#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "SimulationDataFormat/MCCompLabel.h"

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <format>
#include <stdexcept>
#include <type_traits>

namespace o2::gtrk::study
{

namespace
{

/// The cluster indices of every track in \a tracks, whose own range into
/// \a indices \a range gives as (first, count).
template <typename Track, typename Range>
TrackClusters clustersOf(const std::vector<Track>& tracks, const std::vector<int>& indices, Range range)
{
  TrackClusters result;
  result.reserve(tracks.size());
  for (const auto& track : tracks) {
    const auto [first, count] = range(track);
    if (first < 0 || count < 0 || static_cast<std::size_t>(first) + static_cast<std::size_t>(count) > indices.size()) {
      throw std::runtime_error(std::format("a track's clusters [{}, {}) run past the {} cluster indices", first, first + count, indices.size()));
    }
    auto& clusters = result.emplace_back();
    clusters.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
      clusters.push_back(static_cast<std::uint32_t>(indices[static_cast<std::size_t>(first) + static_cast<std::size_t>(i)]));
    }
  }
  return result;
}

template <typename Track>
TrackClusters readEntry(TTree& tree, long long entry, const char* trackBranch, const char* indexBranch, auto range)
{
  std::vector<Track> tracks;
  std::vector<int> indices;
  auto* trackPointer = &tracks;
  auto* indexPointer = &indices;
  tree.SetBranchAddress(trackBranch, &trackPointer);
  tree.SetBranchAddress(indexBranch, &indexPointer);
  tree.GetEntry(entry);
  auto result = clustersOf(tracks, indices, range);
  tree.ResetBranchAddresses();
  return result;
}

template <typename Track>
std::vector<TrackKinematics> readKinematics(TTree& tree, long long entry, const char* branch)
{
  std::vector<Track> tracks;
  auto* pointer = &tracks;
  tree.SetBranchAddress(branch, &pointer);
  tree.GetEntry(entry);
  tree.ResetBranchAddresses();
  std::vector<TrackKinematics> result;
  result.reserve(tracks.size());
  for (const auto& track : tracks) {
    const auto pt = static_cast<float>(track.getPt());
    const auto eta = static_cast<float>(track.getEta());
    const auto phi = static_cast<float>(track.getPhi());
    float chi2 = 0.F;
    int ndf = 0;
    if constexpr (std::is_same_v<Track, o2::its::TrackITS>) {
      chi2 = track.getChi2();
      ndf = 2 * track.getNumberOfClusters() - 5;
    } else {
      chi2 = static_cast<float>(track.getTrackChi2());
      ndf = track.getNDF();
    }
    result.push_back({.pt = pt, .eta = eta, .phi = phi, .chi2 = chi2, .ndf = ndf, .fitValid = std::isfinite(pt) && pt > 0.F && std::isfinite(eta) && std::isfinite(phi)});
  }
  return result;
}

} // namespace

TrackFile::TrackFile(const std::string& path, Chips chips) : mChips(chips), mFile(TFile::Open(path.c_str()))
{
  if (!mFile || mFile->IsZombie()) {
    throw std::runtime_error("cannot open " + path);
  }
  mTree = mFile->Get<TTree>("o2sim");
  if (mTree == nullptr) {
    throw std::runtime_error("no o2sim tree in " + path);
  }
}

TrackFile::~TrackFile() = default;

TrackClusters TrackFile::read(long long entry)
{
  if (entry >= mTree->GetEntries()) {
    return {};
  }
  if (mChips == Chips::kITS) {
    return readEntry<o2::its::TrackITS>(*mTree, entry, "ITSTrack", "ITSTrackClusIdx", [](const o2::its::TrackITS& track) { return std::pair{track.getFirstClusterEntry(), track.getNumberOfClusters()}; });
  }
  return readEntry<o2::mft::TrackMFT>(*mTree, entry, "MFTTrack", "MFTTrackClusIdx", [](const o2::mft::TrackMFT& track) { return std::pair{track.getExternalClusterIndexOffset(), track.getNumberOfPoints()}; });
}

std::vector<bool> TrackFile::fakeFlags(long long entry)
{
  if (entry >= mTree->GetEntries()) {
    return {};
  }
  std::vector<o2::MCCompLabel> labels;
  auto* pointer = &labels;
  mTree->SetBranchAddress(mChips == Chips::kITS ? "ITSTrackMCTruth" : "MFTTrackMCTruth", &pointer);
  mTree->GetEntry(entry);
  mTree->ResetBranchAddresses();
  std::vector<bool> result;
  result.reserve(labels.size());
  for (const auto& label : labels) {
    result.push_back(label.isFake());
  }
  return result;
}

std::vector<TrackKinematics> TrackFile::kinematics(long long entry)
{
  if (entry >= mTree->GetEntries()) {
    return {};
  }
  return mChips == Chips::kITS
           ? readKinematics<o2::its::TrackITS>(*mTree, entry, "ITSTrack")
           : readKinematics<o2::mft::TrackMFT>(*mTree, entry, "MFTTrack");
}

std::vector<o2::track::TrackParCov> TrackFile::innerStates(long long entry)
{
  if (mChips != Chips::kITS) {
    throw std::runtime_error("innerStates: only the ITS's tracks are parameterised in the tracking frame");
  }
  if (entry >= mTree->GetEntries()) {
    return {};
  }
  std::vector<o2::its::TrackITS> tracks;
  auto* pointer = &tracks;
  mTree->SetBranchAddress("ITSTrack", &pointer);
  mTree->GetEntry(entry);
  mTree->ResetBranchAddresses();
  std::vector<o2::track::TrackParCov> result;
  result.reserve(tracks.size());
  for (const auto& track : tracks) {
    result.push_back(track);
  }
  return result;
}

} // namespace o2::gtrk::study

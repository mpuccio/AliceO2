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

#include "GenericTracking/Study/Alice/Hits.h"
#include "GenericTracking/Detectors/TRKFT3.h"

#include "DataFormatsTRKFT3/Hit.h"
#include "ITSMFTSimulation/Hit.h"

#include <TFile.h>
#include <TTree.h>

#include <format>
#include <stdexcept>

namespace o2::gtrk::study
{

namespace
{
/// Every event's hits of type \a HitType from \a branch, as segments, appended to its event's.
template <typename HitType>
void readSegments(TTree& tree, const char* branch, std::vector<std::vector<HitSegment>>& events, std::uint32_t firstSensor = 0U)
{
  std::vector<HitType> hits;
  auto* pointer = &hits;
  tree.SetBranchAddress(branch, &pointer);
  for (std::size_t event = 0; event < events.size(); ++event) {
    tree.GetEntry(static_cast<Long64_t>(event));
    auto& segments = events[event];
    segments.reserve(segments.size() + hits.size());
    for (const auto& hit : hits) {
      segments.push_back({.track = hit.GetTrackID(),
                          .chip = firstSensor + hit.GetDetectorID(),
                          .start = {hit.GetStartX(), hit.GetStartY(), hit.GetStartZ()},
                          .end = {hit.GetX(), hit.GetY(), hit.GetZ()},
                          .momentum = {hit.GetPx(), hit.GetPy(), hit.GetPz()},
                          .time = hit.GetTime()});
    }
  }
  tree.ResetBranchAddresses();
}

TTree* openTree(const std::string& path, std::unique_ptr<TFile>& file)
{
  file.reset(TFile::Open(path.c_str()));
  if (!file || file->IsZombie()) {
    throw std::runtime_error("cannot open " + path);
  }
  auto* tree = file->Get<TTree>("o2sim");
  if (tree == nullptr) {
    throw std::runtime_error("no o2sim tree in " + path);
  }
  return tree;
}
} // namespace

HitFile::HitFile(const std::string& path, Chips chips, const std::string& discPath) : mChips(chips), mTree(openTree(path, mFile))
{
  if (!discPath.empty()) {
    if (chips != Chips::kTRK) {
      throw std::invalid_argument("HitFile: a second file of disc hits is FT3's, for the TRK only");
    }
    mDiscTree = openTree(discPath, mDiscFile);
    if (mDiscTree->GetEntries() != mTree->GetEntries()) {
      throw std::runtime_error(std::format("{} has {} events, {} has {}", path, mTree->GetEntries(), discPath, mDiscTree->GetEntries()));
    }
  }
}

HitFile::~HitFile() = default;

int HitFile::nEvents() const { return static_cast<int>(mTree->GetEntries()); }

const std::vector<HitSegment>& HitFile::event(int event)
{
  // Every event at the first call, in the file's order: the callers ask in their particles' order, and
  // a basket holds many events, so entry by entry would read the file many times over.
  if (!mRead) {
    mRead = true;
    mEvents.resize(static_cast<std::size_t>(mTree->GetEntries()));
    switch (mChips) {
      case Chips::kITS:
        readSegments<o2::itsmft::Hit>(*mTree, "ITSHit", mEvents);
        break;
      case Chips::kMFT:
        readSegments<o2::itsmft::Hit>(*mTree, "MFTHit", mEvents);
        break;
      case Chips::kTRK:
        readSegments<o2::trkft3::Hit>(*mTree, "TRKHit", mEvents);
        if (mDiscTree != nullptr) {
          readSegments<o2::trkft3::Hit>(*mDiscTree, "FT3Hit", mEvents, static_cast<std::uint32_t>(detectors::TRKFT3::kFirstDiscSensor));
        }
        break;
    }
  }
  static const std::vector<HitSegment> kNone;
  return event < 0 || static_cast<std::size_t>(event) >= mEvents.size() ? kNone : mEvents[static_cast<std::size_t>(event)];
}

} // namespace o2::gtrk::study

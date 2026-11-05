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

/// \file Hits.h
/// \brief O2's simulated ITS, MFT and TRK hits, as the segments a particle crossed
///        each chip along.

#ifndef GTRK_STUDY_ALICE_HITS_H_
#define GTRK_STUDY_ALICE_HITS_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"

#include <array>
#include <concepts>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace o2::gtrk::study
{

/// Whose MC files: the ITS's, the MFT's, or the ALICE 3 tracker's (TRK's barrel, then FT3's discs).
enum class Chips : std::uint8_t { kITS,
                                  kMFT,
                                  kTRK };

/// Whose MC files a layout's are.
template <DetectorLayout Det>
  requires std::same_as<Det, detectors::ITS> || std::same_as<Det, detectors::MFT> || std::same_as<Det, detectors::TRKFT3>
[[nodiscard]] constexpr Chips chipsOf()
{
  if constexpr (std::same_as<Det, detectors::ITS>) {
    return Chips::kITS;
  } else if constexpr (std::same_as<Det, detectors::MFT>) {
    return Chips::kMFT;
  } else {
    return Chips::kTRK;
  }
}

/// One particle's crossing of one chip, in the global frame [cm, GeV].
struct HitSegment {
  int track = -1;          ///< MC track ID within its event
  std::uint32_t chip = 0U; ///< our sensor ID: O2's chip ID; FT3's after TRK's (TRKFT3::kFirstDiscSensor)
  std::array<double, 3> start{};
  std::array<double, 3> end{};
  std::array<double, 3> momentum{}; ///< at the start
  double time = 0.;                 ///< the time of flight at the end: the order along the particle's path
};

/// `sgn_HitsITS.root`, `sgn_HitsMFT.root` or `sgn_HitsTRK.root`, one tree entry per MC event;
/// for the TRK, with `sgn_HitsFT3.root` beside it, FT3's hits too.
class HitFile
{
 public:
  /// Throws std::runtime_error if \a path, or a non-empty \a discPath, has
  /// no o2sim tree. \a discPath, FT3's hits, only for Chips::kTRK.
  HitFile(const std::string& path, Chips chips, const std::string& discPath = {});
  ~HitFile();
  HitFile(const HitFile&) = delete;
  HitFile& operator=(const HitFile&) = delete;
  HitFile(HitFile&&) = delete;
  HitFile& operator=(HitFile&&) = delete;

  /// The segments of event \a event; every event's read at the first call and kept.
  [[nodiscard]] const std::vector<HitSegment>& event(int event);

  /// How many MC events the file holds.
  [[nodiscard]] int nEvents() const;

 private:
  Chips mChips;
  std::unique_ptr<TFile> mFile;
  TTree* mTree = nullptr;
  std::unique_ptr<TFile> mDiscFile; ///< FT3's, for the TRK
  TTree* mDiscTree = nullptr;
  bool mRead = false;
  std::vector<std::vector<HitSegment>> mEvents; ///< by event
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_ALICE_HITS_H_

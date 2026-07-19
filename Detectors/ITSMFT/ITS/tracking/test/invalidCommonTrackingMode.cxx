// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include <cstdlib>

#include "ITSCommonTracking/CommonTrackingParameters.h"

int main(int argc, char** argv)
{
  if (argc != 2) {
    return 2;
  }
  const auto mode = static_cast<o2::its::TrackingMode::Type>(std::atoi(argv[1]));
  (void)o2::its::commontracking::makeTrackingParameters(mode);
  return 0;
}

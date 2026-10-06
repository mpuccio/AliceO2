#!/usr/bin/env python3
# Copyright 2019-2026 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

"""Check the tolerant comparator's rejection rules against a real output file."""
import sys
from pathlib import Path
import ROOT

ROOT.gROOT.LoadMacro(str(Path(__file__).with_name("compareITSTracking.C")))
ROOT.gInterpreter.Declare(r'''
template <typename Tracks, typename Indices, typename ROFs, typename Labels>
void testTrackingComparison(const Tracks& tracks, const Indices& indices, const ROFs& rofs, const Labels& labels)
{
  const auto require = [](bool condition) {
    if (!condition) { throw std::runtime_error{"Comparison self-test failed"}; }
  };
  const auto rejected = [&](const auto& changed, const auto& changedIndices, const auto& changedLabels) {
    try { equivalentTrackingEvent(tracks, indices, rofs, labels, changed, changedIndices, rofs, changedLabels); }
    catch (const std::runtime_error&) { return true; }
    return false;
  };
  require(!tracks.empty());
  require(equivalentTrackingEvent(tracks, indices, rofs, labels, tracks, indices, rofs, labels));
  auto changed = tracks;
  changed[0].setChi2(tracks[0].getChi2() + 1.e-4f);
  require(!rejected(changed, indices, labels));
  changed[0].setChi2(tracks[0].getChi2() + 1.f + std::abs(tracks[0].getChi2()));
  require(rejected(changed, indices, labels));
  changed = tracks;
  changed[0].setParam(tracks[0].getParam(0) + 1.f + std::abs(tracks[0].getParam(0)), 0);
  require(rejected(changed, indices, labels));
  changed = tracks;
  changed[0].getParamOut().setCov(tracks[0].getParamOut().getCov()[0] + 1.f, 0);
  require(rejected(changed, indices, labels));
  changed = tracks;
  changed[0].setParam(std::numeric_limits<float>::quiet_NaN(), 0);
  require(rejected(changed, indices, labels));
  changed = tracks;
  changed[0].setPattern(tracks[0].getPattern() ^ 1u);
  require(rejected(changed, indices, labels));
  changed = tracks;
  auto changedIndices = indices;
  changedIndices[0] ^= 1;
  require(rejected(changed, changedIndices, labels));
  auto changedLabels = labels;
  changedLabels[0].setFakeFlag(!labels[0].isFake());
  require(rejected(changed, indices, changedLabels));
  changedLabels = labels;
  bool reordered = false;
  for (const auto& rof : rofs) {
    if (rof.getNEntries() < 2) { continue; }
    const size_t i = rof.getFirstEntry();
    std::swap(changed[i], changed[i + 1]);
    std::swap(changedLabels[i], changedLabels[i + 1]);
    reordered = true;
    break;
  }
  require(reordered);
  require(!rejected(changed, indices, changedLabels));
}
''')
f = ROOT.TFile.Open(sys.argv[1])
if not f or f.IsZombie():
    raise RuntimeError("Cannot open test input")
t = f.Get("o2sim")
t.GetEntry(0)
ROOT.testTrackingComparison(t.ITSTrack, t.ITSTrackClusIdx, t.ITSTracksROF, t.ITSTrackMCTruth)
f.Close()
print("Tracking comparison self-test passed")

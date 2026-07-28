// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Static-diamond timing correctness (ITSCommonCATrackerParam.useDiamond /
// TrackerTraits::computeLayerTrackletsForPolicy, TrackerTraits.cxx). A
// diamond vertex has no genuine per-event timing, but Vertex::getTimeStamp()
// is a fixed-width TimeEstBC (error field is only 16 bits -- see
// DataFormatsITS/TimeEstBC.h) that cannot literally span a whole real
// TimeFrame's BC range (order 1e5-1e6 BC). TrackerTraits.cxx's
// diamondVertexForROF() resolves this by deriving the diamond's timestamp
// fresh for each (layer, rofId) it is tested against, from that ROF's own
// real interval envelope (ROFOverlapTableView::getLayer(layer)
// .getROFTimeBounds(rofId, true) -- the exact function every real-vertex
// ROF timing check in this file already uses). This is a focused,
// standalone proof of that mechanism's core mathematical property, using
// the same underlying legacy o2::its types (ROFOverlapTable,
// ROFVertexLookupTable, Vertex, TimeEstBC) production code calls, without
// needing a full TrackerTraits/TimeFrame fixture.

#define BOOST_TEST_MODULE ITSMFT StaticDiamondROFCompatibility
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cstdint>

#include "DataFormatsITS/TimeEstBC.h"
#include "DataFormatsITS/Vertex.h"
#include "ITStracking/ROFLookupTables.h"

using namespace o2::its;

namespace
{
constexpr int NLayers = 7;

// Mirrors o2::itsmft::tracking::diamondVertexForROF() (TrackerTraits.cxx):
// derive a diamond vertex's timestamp from one specific (layer, rofId)'s own
// real interval envelope. Kept independent (not #include-ing the anonymous-
// namespace production helper) so this test exercises the same public
// building blocks the production code is built from, not a re-export of it.
template <typename View>
Vertex diamondVertexForROF(const Vertex& base, const View& view, int layer, int rofId)
{
  Vertex v = base;
  v.setTimeStamp(view.getLayer(layer).getROFTimeBounds(rofId, true));
  return v;
}

// Builds an NLayers-uniform ROFOverlapTable/ROFVertexLookupTable pair for
// one explicit ROFTiming, exactly as ITSMFTTrackingInterface::
// configureROFLookupTables() does for a real, uniform-across-layers
// TimeFrame, then asserts the diamond derived for every (layer, rofId) in
// [0, nROFsTF) is compatible with that very ROF via the real, unconditional
// isVertexCompatible() -- the same function TrackerTraits.cxx now calls
// with no useDiamond-gated bypass.
void checkEveryROFCompatible(uint32_t nROFsTF, uint32_t rofLength, uint32_t rofDelay, uint32_t rofBias, uint32_t addTimeErr)
{
  LayerTiming timing{};
  timing.mNROFsTF = nROFsTF;
  timing.mROFLength = rofLength;
  timing.mROFDelay = rofDelay;
  timing.mROFBias = rofBias;
  timing.mROFAddTimeErr = addTimeErr;

  ROFOverlapTable<NLayers> rofTable;
  ROFVertexLookupTable<NLayers> vtxTable;
  for (int layer = 0; layer < NLayers; ++layer) {
    rofTable.defineLayer(layer, timing);
    vtxTable.defineLayer(layer, timing);
  }
  rofTable.init();
  vtxTable.init();

  const auto rofView = rofTable.getView();
  const auto vtxView = vtxTable.getView();

  const float pos[3] = {0.f, 0.f, 0.f};
  const float cov[6] = {25.e-6f, 0.f, 0.f, 25.e-6f, 0.f, 36.f};
  const Vertex base(pos, cov, 1, 1.f);

  for (int layer = 0; layer < NLayers; ++layer) {
    for (uint32_t rof = 0; rof < nROFsTF; ++rof) {
      const auto diamond = diamondVertexForROF(base, rofView, layer, static_cast<int>(rof));
      BOOST_CHECK_MESSAGE(vtxView.isVertexCompatible(layer, rof, diamond),
                           "layer=" << layer << " rof=" << rof << " nROFsTF=" << nROFsTF
                                    << " rofLength=" << rofLength << " rofDelay=" << rofDelay
                                    << " rofBias=" << rofBias << " addTimeErr=" << addTimeErr);
    }
  }
}
} // namespace

// A realistic ITS continuous-readout fixture: nROFsTF=2304 matches the
// Gate 3 Slice 3 characterization fixture's "Input ROFs" count
// (gate3-slice3-its-ca-validation/characterization_summary.md), rofLength
// 198 BC == LHCMaxBunches/18. Total TF span (2304*198 ~ 456k BC) vastly
// exceeds TimeEstBC's 16-bit error-field width (65535): this is exactly the
// case a single, TF-wide static TimeEstBC could not represent, which
// per-ROF derivation sidesteps entirely.
BOOST_AUTO_TEST_CASE(StaticDiamondCompatibleWithEveryROF_ContinuousReadout)
{
  checkEveryROFCompatible(2304, 198, 0, 0, 0);
}

// Non-zero delay/bias/error margin (a triggered-like or offset-clock
// configuration), and a much smaller ROF count.
BOOST_AUTO_TEST_CASE(StaticDiamondCompatibleWithEveryROF_DelayBiasAndErrorMargin)
{
  checkEveryROFCompatible(500, 198, 37, 5, 10);
}

// Degenerate but valid (nROFsTF == 1): the first and last ROF are the same
// ROF, exercising the tfBegin == tfEnd-anchor edge case.
BOOST_AUTO_TEST_CASE(StaticDiamondCompatibleWithEveryROF_SingleROF)
{
  checkEveryROFCompatible(1, 40, 0, 0, 0);
}

// A single ROF's own tracklet-time window (the intersection of two
// overlapping ROFs' bounds, TrackerTraits.cxx's `ts`) is always a subset of
// either source ROF's own bounds -- so it is always TimeEstBC::isCompatible
// with a diamond derived from that ROF, exactly the property
// TrackerTraits.cxx's second (post-bypass-removal) check relies on. This
// checks that subset property directly, independent of the vertex-lookup
// table above.
BOOST_AUTO_TEST_CASE(TrackletTimeWindowIsAlwaysSubsetOfItsSourceROFWindow)
{
  LayerTiming timing{};
  timing.mNROFsTF = 200;
  timing.mROFLength = 198;
  timing.mROFAddTimeErr = 5;

  ROFOverlapTable<NLayers> rofTable;
  for (int layer = 0; layer < NLayers; ++layer) {
    rofTable.defineLayer(layer, timing);
  }
  rofTable.init();
  const auto rofView = rofTable.getView();

  constexpr int fromLayer = 0;
  constexpr int toLayer = 1;
  for (uint32_t pivotROF = 0; pivotROF + 1 < timing.mNROFsTF; ++pivotROF) {
    if (!rofView.doROFsOverlap(fromLayer, pivotROF, toLayer, pivotROF)) {
      continue;
    }
    const auto ts = rofView.getTimeStamp(fromLayer, pivotROF, toLayer, pivotROF);
    const auto t0 = rofView.getLayer(fromLayer).getROFTimeBounds(pivotROF, true);
    BOOST_CHECK_MESSAGE(ts.isCompatible(t0), "pivotROF=" << pivotROF << " ts not compatible with its own fromLayer/pivotROF window");
  }
}

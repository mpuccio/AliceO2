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

/// \file RecoParam.h
/// \brief Each detector's tracking settings as a ConfigurableParam.
///
/// A ConfigurableParam, where `--configKeyValues "GTRKReco<DET>.<field>=<value>"`
/// reaches them and every device of the chain sees the same value.
///
/// The fields are one struct, `RecoSettings`; each detector has its own
/// parameter of them, `RecoParam<DetID>` (as `DPLAlpideParam<N>`):
/// `GTRKRecoITS`, `GTRKRecoMFT`, `GTRKRecoTRK`, so that one workflow carries
/// every detector's. A detector layout's is `recoParam<Det>()`.

#ifndef GTRK_WORKFLOW_RECOPARAM_H_
#define GTRK_WORKFLOW_RECOPARAM_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DetectorsCommonDataFormats/DetID.h"

// ENABLE_UPGRADES comes with DetID.h (UpgradesStatus.h).
#ifdef ENABLE_UPGRADES
#include "GenericTracking/Detectors/TRKFT3.h"
#endif

#include <cstdint>
#include <string_view>
#include <type_traits>
#include <vector>

namespace o2::gtrk
{

struct RecoSettings {
  /// O2's `addTimeError` for the ITS and the MFT [BC]: every readout frame's
  /// time box widened by this much at both ends, for a hit recorded in the
  /// neighbouring frame. One entry per layer, or
  /// none for 0 on every layer, O2's default.
  std::vector<int> frameMarginBC;

  /// Width of the outer time block in the measurement grid [ps], 0 for one
  /// block. This is the minimum slice quantum and is chosen from the worker
  /// memory budget; a detector with readout frames takes the whole number of
  /// frames nearest it.
  std::int64_t blockWidthPS = 10'000'000;

  // The grid's cells are derived from the iterations' queries, with no key of
  // their own (GridSizing.h).
  /// The mass hypothesis [GeV], positive (`checkRecoParam`): beta for the
  /// scattering, in the windows at pT_min and in the fit, and the fit's energy
  /// loss. The charged pion's by default.
  double massHypothesisGeV = o2::constants::physics::MassPionCharged;
  /// The uniform material per layer [X/X0] that sizes the windows and the
  /// grid, and the scattering when no material LUT is given; with a LUT the
  /// tracking takes its budgets.
  double materialX2X0 = 0.01;
  /// The windows' position resolution [cm].
  double positionSigma = 0.005;
  /// Added in quadrature to every cluster's local errors [cm], for what the
  /// cluster shape does not know: the alignment's residuals on real data. O2
  /// adds 10 um per layer in its async reconstruction (ITSCATrackerParam.sysErrY2,
  /// sysErrZ2 = 100e-8 cm^2); 0, the default, is the MC's dictionary error alone.
  double clusterSystematicSigmaCM = 0.;
  /// The windows' and the cuts' nSigma, the run's (every iteration's).
  double windowNSigma = 3.;
  /// How many sigma of the beam's z width the tracking reaches about the
  /// beam's z (`beamSpotConstraint`): the z range it is to be efficient in,
  /// for the beam-line windows and pair cut, the vertexing pass, the grid and
  /// the vertex-list iterations alike.
  double beamZReachNSigma = 3.;

  /// The own vertexer: before the tracking, a pass at
  /// the beam line on the three innermost layers, one slice over the
  /// timeframe, tracklets and triplets at `vertexerMinPtGeV`; its triplets'
  /// z at the beam line, frame by frame, give the vertex list the
  /// vertex-list iterations read (`iterationVertexSource`). Off by default.
  /// The rest are the vertexer's knobs (VertexStage.h, VertexerSettings).
  bool vertexer = false;
  double vertexerMinPtGeV = 0.5;
  double vertexerMaxChi2RZ = 5.;   ///< a triplet's r-z line chi2 cut, 0 for none
  double vertexerMaxDcaCM = 0.02;  ///< its circle's transverse distance to the beam axis [cm], 0 for none
  int vertexerWindowUM = 150;      ///< a peak's window [um]
  int vertexerClaimUM = 2000;      ///< the half-width a vertex claims around its peak [um]
  int vertexerMinContributors = 4; ///< the least candidates a peak's window holds
  /// The least in a frame where none holds vertexerMinContributors: that
  /// frame's rounds run again at this minimum; 0, none, by default.
  int vertexerAloneMinContributors = 0;
  double vertexerSignificance = 5.; ///< a peak's excess over the sidebands, in their Poisson sigmas
  /// Whether the vertexing pass tests its pairs and triplets against their
  /// measurements' cluster shapes (`iterationShapeTest`); 0, off, by default.
  int vertexerShapeTest = 0;

  /// Tracking iterations, run one after another: one per entry of
  /// `iterationMinPtGeV`. Each other per-iteration vector holds one entry per
  /// iteration, or none for its default in every iteration (GridConfig.h,
  /// `iterationCount`). They are set whole, as O2's parameters take vectors:
  /// `GTRKRecoITS.iterationMinPtGeV=[1,0.5];GTRKRecoITS.iterationMinClusters=[7,4]`.
  /// What an iteration loosens is its minimum pT; nSigma (`windowNSigma`) is
  /// the run's, the same in every iteration.
  ///
  /// Three iterations at 1.0, 0.5 and 0.3 GeV by default; a schedule sets its
  /// own.
  std::vector<double> iterationMinPtGeV{1., 0.5, 0.3};
  /// Where each iteration's vertex seeds come from: 0 the beam line (the
  /// default), 1 the timeframe's vertex list.
  std::vector<int> iterationVertexSource;
  /// The shortest road each iteration fits, in clusters, 4 by default; a road
  /// is at least 4 clusters, two triplets, whatever is set. An iteration takes
  /// only the links on a path of this many layers (Iteration.h, `takenLinks`).
  /// ITS asks 7 of its first two iterations and 4 of the third.
  std::vector<int> iterationMinClusters;
  /// The kinks each iteration's fit may take per pass;
  /// 0 for none, 1 by default.
  std::vector<int> iterationMaxKinks;
  /// The kink allowance of each iteration's triplet and neighbour cuts: a
  /// triplet failing the nominal cuts is kept, flagged kinked, if it passes with the middle
  /// layer's scattering at pT_min scaled by this factor; 1 (the default) for
  /// none.
  std::vector<double> iterationKinkSigma;
  /// The triplet and neighbour cuts' own kink allowance, each iteration's:
  /// above 0 it replaces
  /// iterationKinkSigma there, which then sets only the road fit's and the
  /// follow's; 0, the default, iterationKinkSigma everywhere.
  std::vector<double> iterationGraphKinkSigma;
  /// A displaced iteration's bounds [cm]:
  /// how far its tracks' circles may pass from the beam line, and how far
  /// the z of that closest approach from the beam spot's centre. 0, the
  /// default, a prompt iteration; above, on the beam line and with both set,
  /// the displaced pair test and its windows.
  std::vector<double> iterationImpactMaxCm;
  std::vector<double> iterationImpactZMaxCm;
  /// Whether each iteration makes tracklets that leave a hole, whose segment
  /// crosses another layer's envelope (the default); 0 for none, so the
  /// tracks past a shorter layer's end are kept and those through a gap
  /// are not.
  std::vector<int> iterationSkipLinks;
  /// The layers each iteration's links join, a bit per layer; 0 for every
  /// layer (the default). The vertexing pass's 7 takes layers 0-2 alone.
  std::vector<long> iterationLayerMask;
  /// The layers each iteration's roads may end on, a bit per layer; 0 for any
  /// (the default). O2's ITS async asks layer 3 or 6 of its third,
  /// `(1 << 6) + (1 << 3)` = 72.
  std::vector<long> iterationRoadEndLayers;
  /// The layers each iteration's roads may start on, their innermost
  /// measurement's, a bit per layer; 0 for any (the default).
  std::vector<long> iterationRoadStartLayers;
  /// The lowest pT of a track by its cluster count [GeV], 0 for none -- the
  /// final track's: after the follow where the iteration follows, else the
  /// road's. Iteration
  /// `i`'s floor for `n` clusters is entry `16 i + n`, so 16 entries per
  /// iteration, or none for no floor at all. O2's ITS async asks 1/12 of 7
  /// clusters in every iteration, and 1/5, 1 and 1/6 of 6, 5 and 4 in its
  /// third (entries 7, 23, 39, 38, 37, 36).
  static constexpr int kClusterSlots = 16;
  std::vector<double> iterationClusterMinPtGeV;
  /// How many of a track's innermost clusters each iteration lets it share
  /// with another; 0, strict, by default.
  std::vector<int> iterationSharedInnermost;
  /// Whether each iteration follows its selected tracks after their refit,
  /// outward and then inward, as a CKF does, before the next iteration;
  /// 1, on, by default.
  std::vector<int> iterationFollow;
  /// The holes each iteration's follow lets a track cross in each direction,
  /// a sensor crossed without a row; 0 by default.
  std::vector<int> iterationFollowHoles;
  /// Whether each iteration decides its roads' outermost row again with the
  /// follow: fitted without it, the follow's first row in its place where it
  /// takes another; without iterationFollow,
  /// that row alone. 0, off, by default; not with a kink allowance.
  std::vector<int> iterationFollowEnd;
  /// The innermost clusters each iteration's road fit may drop from a road it
  /// refuses, one more each time, down to the iteration's minimum, the first
  /// it takes the road's record; 0, off, by
  /// default.
  std::vector<int> iterationInnerTrim;
  /// Which of a barrel plane's two crossings each iteration's follow refit
  /// takes (`CrossingRule`): 0 the first ahead
  /// (the default), 1 of those ahead the one nearer the plane's centre, 2 one
  /// inside the sensor first, else nearer the centre.
  std::vector<int> iterationCrossingRule;
  /// Whether each iteration's follow refit, where it refuses the rows the
  /// follow added, keeps the longest prefix of them it accepts instead of
  /// none; 0, off, by default.
  std::vector<int> iterationFollowPrefix;
  /// Whether each iteration's follow takes, where more than one row passes a
  /// step, the smallest chi2 instead of the lookahead's confirmed row; 0, off, by
  /// default.
  std::vector<int> iterationFollowGreedy;
  /// The chi2 below which each iteration's follow keeps a direction's
  /// unconfirmed last row, which the unconfirmed drop removes otherwise; 0, none kept, by
  /// default.
  std::vector<double> iterationFollowKeepChi2;
  /// Whether each iteration refuses its inward walk's last row, where a
  /// nominal step took it, whose chi2 against that step's prediction updated
  /// with the beam line is above the gate; needs the iteration's follow and
  /// beam constraint; 0, off, by default.
  std::vector<int> iterationFollowBeamEnd;
  /// Whether each iteration tests its pairs and triplets against their
  /// measurements' cluster shapes: a candidate's local slope at each of its measurements -- its line's
  /// |dz| over its length along the sensor's normal -- must lie in the
  /// measurement's slope interval, which the measurement producer makes from
  /// its detector's calibrated band (ClusterShapeBand.h); 0, off, by default.
  std::vector<int> iterationShapeTest;
  /// Each iteration's follow's end decisions on the last row's own update:
  /// inward, the sqrt|V| [cm^2]
  /// above which a last row at the layout's end is refused; outward, the
  /// density [cm^-2] above which the unconfirmed drop keeps a last row; 0, off, by default.
  std::vector<double> iterationFollowEndMaxSqrtDet;
  std::vector<double> iterationFollowKeepDensity;
  /// The Delta chi2 below which each iteration's follow trims a track's
  /// last outward row where it is one of two rows or more passing the gate
  /// and the track keeps the cuts after the follow without it: the other
  /// row's chi2 less its own; 0, off, by default; above the gate, every such row.
  std::vector<double> iterationFollowTrimMargin;
  /// The clusters each iteration's tracks must hold after the follow; 0 for
  /// none (the default). A shorter track is dropped and its rows given back.
  std::vector<int> iterationMinClustersAfterFollow;
  /// Whether each iteration updates its roads with the beam line on their
  /// innermost surface and refuses them above the update's gate, to start
  /// the follow inward from the constrained state; 0, off, by default.
  std::vector<int> iterationBeamConstraint;
  /// Whether each iteration follows every earlier iteration's tracks too,
  /// with its own settings, their rows kept, as a recovery pass does; 0, off,
  /// by default. Needs the iteration's
  /// iterationFollow; iterationMinClustersAfterFollow applies to the
  /// iteration's own tracks only.
  std::vector<int> iterationFollowEarlier;
  /// Whether each iteration takes the previous one's graph instead of
  /// building its own, one graph for both rounds; 0, off, by default. Needs
  /// the previous iteration's graph settings, and its graph kink allowance 1
  /// and this one's above.
  std::vector<int> iterationSharedGraph;

  /// Requested width of one time slice [ps]. The actual width is rounded up
  /// to whole grid blocks, and zero means one slice over the timeframe.
  std::int64_t sliceWidthPS = 10'000'000;

  /// How far a slice reads beyond its core on each side [ps], rounded up to
  /// whole blocks. The workflow raises this to `requiredHaloPS` when
  /// necessary, so the output does not depend on the slicing.
  std::int64_t sliceHaloPS = 100'000;

  /// Consecutive source elements assigned to one CPU task.
  int trackletChunkSize = 256;

  /// Maximum CPU workers; zero uses TBB's default.
  int threads = 0;

  /// The host tracker's memory [GB]: a stage that would hold more fails with
  /// kOutOfMemory, and the timeframe's tracking with it, instead of growing.
  /// Zero, the default, for no limit.
  double memoryLimitGB = 0.;

  /// A timeframe whose tracking runs out of memory (memoryLimitGB) gives
  /// empty products and is counted, instead of failing the run. Off, the
  /// default, as ITS's DropTFUponFailure: the run fails.
  bool dropTFUponFailure = false;
};

/// Detector \a N's settings (`o2::detectors::DetID`): `GTRKReco<DET>`.
// O2ParamDef deletes the copy constructor and supplies the rest of the
// boilerplate, so the rule about declaring all five has nothing to act on here.
template <int N>
// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions)
struct RecoParam : public o2::conf::ConfigurableParamPromoter<RecoParam<N>, RecoSettings> {
  static constexpr std::string_view getParamName() { return N == o2::detectors::DetID::ITS ? "GTRKRecoITS" : (N == o2::detectors::DetID::MFT ? "GTRKRecoMFT" : "GTRKRecoTRK"); }
  O2ParamDef(RecoParam, getParamName().data());
};

template <int N>
RecoParam<N> RecoParam<N>::sInstance;

/// A detector layout's O2 detector: its parameter's, its data origin's
/// (DataSpecs.h's `detectorOrigin`).
template <DetectorLayout Det>
struct DetIDOf;
template <>
struct DetIDOf<detectors::ITS> : std::integral_constant<o2::detectors::DetID::ID, o2::detectors::DetID::ITS> {
};
template <>
struct DetIDOf<detectors::MFT> : std::integral_constant<o2::detectors::DetID::ID, o2::detectors::DetID::MFT> {
};
#ifdef ENABLE_UPGRADES
template <>
struct DetIDOf<detectors::TRKFT3> : std::integral_constant<o2::detectors::DetID::ID, o2::detectors::DetID::TRK> {
};
#endif

/// A detector layout's parameter.
template <DetectorLayout Det>
struct RecoParamOf {
  using type = RecoParam<DetIDOf<Det>::value>; // NOLINT(readability-identifier-naming): a type trait's member
};

/// The settings of \a Det's parameter.
template <DetectorLayout Det>
[[nodiscard]] const RecoSettings& recoParam()
{
  return RecoParamOf<Det>::type::Instance();
}

} // namespace o2::gtrk

#endif // GTRK_WORKFLOW_RECOPARAM_H_

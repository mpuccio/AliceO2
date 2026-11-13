// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// Distributed under the GNU GPL v3. CERN privileges and immunities are preserved.

#include "StrangenessTracking/StrangenessTracker.h"
#include "StrangenessTracking/ITSShortTrackFit.h"
#include "ITStracking/Configuration.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsBase/GlobalParams.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "ITStracking/IndexTableUtils.h"
#include "ITSMFTTracking/ROFLookupTables.h"
#include "ITSMFTTracking/ITSTrackingConfigParam.h"
#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TrkClusRef.h"
#include <algorithm>
#include <bit>
#include <cmath>
#include <map>
#include <limits>
#include <set>
#include <tuple>

namespace o2::strangeness_tracking
{
namespace
{
using Track = o2::track::TrackParCovF;
// Fit/search bookkeeping is transient. Public device outputs reuse StrangeTrack
// and TrackITS rather than persisting this implementation-specific state.
struct ClusterCandidate {
  Track mTrack, mOuter;
  bool isShort = false;
  o2::dataformats::StrangeTrack source{};
  o2::dataformats::kPartType mPartType{};
  int mPVRef = -1, mDecayRef = -1;
  std::array<int, 7> mClusterRefs{{-1, -1, -1, -1, -1, -1, -1}};
  std::array<int, 7> mROFRefs{{-1, -1, -1, -1, -1, -1, -1}};
  std::array<uint8_t, 7> mClusterSizes{};
  double mTimeMinBC = 0., mTimeMaxBC = 0.;
  float mHitChi2 = 0.f, mPVChi2 = 0.f;
  int mNDF = 0;
  uint8_t mMissingLayerMask = 0;
  int getNClusters() const
  {
    return std::count_if(mClusterRefs.begin(), mClusterRefs.end(), [](int id) { return id >= 0; });
  }
};
struct Hit {
  int ref = -1, rof = -1, clock = -1, layer = -1;
  float x = 0, y = 0, z = 0;
  double lower = 0, upper = 0;
};
// Reject invalid states before they can participate in the global ordering.
bool valid(const Track& t)
{
  if (std::abs(t.getSnp()) >= 1.f) {
    return false;
  }
  // A Cholesky check also catches finite but non-positive covariance matrices.
  double l[5][5]{};
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      double v = t.getCov()[i * (i + 1) / 2 + j];
      for (int k = 0; k < j; ++k) {
        v -= l[i][k] * l[j][k];
      }
      if (i == j) {
        if (v <= 0.) {
          return false;
        }
        l[i][j] = std::sqrt(v);
      } else {
        l[i][j] = v / l[j][j];
      }
    }
  }
  return true;
}
int holes(const ClusterCandidate& t) { return std::popcount(unsigned(t.mMissingLayerMask)); }
bool better(const ClusterCandidate& a, const ClusterCandidate& b)
{
  if (a.getNClusters() != b.getNClusters()) {
    return a.getNClusters() > b.getNClusters();
  }
  if (holes(a) != holes(b)) {
    return holes(a) < holes(b);
  }
  const auto qa = a.mHitChi2 / std::max(1, a.getNClusters());
  const auto qb = b.mHitChi2 / std::max(1, b.getNClusters());
  if (qa != qb) {
    return qa < qb;
  }
  return a.mClusterRefs < b.mClusterRefs;
}
} // namespace

void StrangenessTracker::remapExtensionSeeds(int thread, const std::vector<int>& v0, const std::vector<int>& casc, const std::vector<int>& three)
{
  mExtensionSeedsRemapped = true;
  for (auto& seed : mExtensionSeeds[thread]) {
    for (const auto& accepted : mStrangeTrackVec[thread]) {
      seed.matched |= accepted.mPartType == seed.source.mPartType && accepted.mDecayRef == seed.source.mDecayRef;
    }
    const auto& map = seed.source.mPartType == o2::dataformats::kStrkV0 ? v0 : seed.source.mPartType == o2::dataformats::kStrkCascade ? casc
                                                                                                                                      : three;
    seed.source.mDecayRef = map.at(seed.source.mDecayRef);
  }
}

void StrangenessTracker::finishClusterTracking()
{
  if (!clusterTrackingEnabled() || mClusterTrackingFinished) {
    return;
  }
  mClusterTrackingFinished = true;
  if (!mClusterTrackingData) {
    throw std::runtime_error("ITS extension/short tracking has no current reconstruction inputs");
  }
  const auto* grp = o2::base::GRPGeomHelper::instance().getGRPECS();
  if (!grp || !grp->isDetContinuousReadOut(o2::detectors::DetID::ITS)) {
    throw std::runtime_error("ITS extension/short tracking currently requires continuous ITS readout and GRPECS");
  }
#ifdef ENABLE_UPGRADES
  if (o2::GlobalParams::Instance().withITS3) {
    throw std::runtime_error("ITS extension/short tracking currently supports ITS2 geometry only");
  }
#endif
  const auto& data = *mClusterTrackingData;
  const auto& cfg = *mStrParams;
  if (cfg.extensionMinHits < 2 || cfg.extensionMinHits > 7 || cfg.extensionMaxMissing < 0 || cfg.extensionMaxMissing > 1 ||
      cfg.extensionMaxBranches < 1 || cfg.shortTrackMaxPerPV < 1 || cfg.shortTrackMaxPerTF < 1 ||
      mPVTimeMargin < 0 || !(cfg.shortTrackMinPt > 0) || !(cfg.clusterSearchNSigma > 0) || !(cfg.extensionMaxHitChi2 > 0) ||
      !(cfg.shortTrackMaxPointingChi2 > 0) || !(cfg.shortTrackMaxChi2NDF > 0)) {
    throw std::runtime_error("Invalid ITS extension/short tracking configuration");
  }
  if (!cfg.itsOnly && !data.hasITSABClusterUsage()) {
    throw std::runtime_error("ITS extension/short tracking requires afterburner references (or explicit itsOnly profile)");
  }
  mStrangenessTrackingStats.bz = mBz;
  mStrangenessTrackingStats.itsOnly = cfg.itsOnly;
  mStrangenessTrackingStats.clustersPerLayer = data.getITSPerLayer();
  mStrangenessTrackingStats.decayFollowingEnabled = cfg.enableDecayFollowing;
  mStrangenessTrackingStats.tripletsEnabled = cfg.enableShortIB;
  mStrangenessTrackingStats.tripletsDeferredZeroField = cfg.enableShortIB && std::abs(mBz) < 0.01f;
  auto* geom = o2::its::GeometryTGeo::Instance();
  const o2::its::TrackingParameters itsParams;
  o2::its::IndexTableUtils<7> utils;
  utils.setTrackingParameters(itsParams);
  const auto& alpide = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>::Instance();
  const auto& itsCfg = o2::its::TrackerParamConfig::Instance();
  std::array<o2::itsmft::tracking::LayerTiming, 7> timing;
  std::vector<Hit> hits(mInputITSclusters.size());
  // Sparse ITS phi-z grid per physical layer and timing-clock ROF. The stored
  // value remains a dense transient index; persistent refs keep the source slot.
  using Bin = std::tuple<int, int, int>; // clock, phi bin, z bin
  std::array<std::map<Bin, std::vector<int>>, 7> grid;
  std::array<float, 7> radiusMin, radiusMax, maxErrY{}, maxErrZ{};
  for (int layer = 0; layer < 7; ++layer) {
    radiusMin[layer] = radiusMax[layer] = itsParams.LayerRadii[layer];
    if (alpide.getROFLengthInBC(layer) <= 0 || alpide.getROFBiasInBC(layer) < 0 || alpide.getROFDelayInBC(layer) < 0 || itsCfg.addTimeError[layer] < 0) {
      throw std::runtime_error("Invalid configured ITS timing for cluster tracking");
    }
    const auto length = alpide.getROFLengthInBC(layer);
    if (o2::constants::lhc::LHCMaxBunches % length != 0) {
      throw std::runtime_error("ITS extension/short tracking ROF length must divide an orbit");
    }
    const auto nROFs = (o2::constants::lhc::LHCMaxBunches / length) * o2::base::GRPGeomHelper::getNHBFPerTF();
    if (nROFs <= 0) {
      throw std::runtime_error("ITS extension/short tracking requires the configured TF orbit count");
    }
    timing[layer] = {uint32_t(nROFs), uint32_t(alpide.getROFLengthInBC(layer)), uint32_t(alpide.getROFDelayInBC(layer)),
                     uint32_t(alpide.getROFBiasInBC(layer)), uint32_t(itsCfg.addTimeError[layer])};
  }
  const int nSlots = mInputITSclusters.getNLayers();
  if (nSlots == 1) {
    for (int layer = 1; layer < 7; ++layer) {
      if (timing[layer].mROFLength != timing[0].mROFLength || timing[layer].mROFBias != timing[0].mROFBias || timing[layer].mROFDelay != timing[0].mROFDelay) {
        throw std::runtime_error("Staggered ITS extension/short tracking timing requires layer-separated cluster inputs");
      }
    }
  }
  for (int slot = 0; slot < nSlots; ++slot) {
    const auto rofs = data.getITSClustersROFRecords(slot);
    if (rofs.size() != timing[slot].mNROFsTF) {
      throw std::runtime_error("ITS extension/short tracking ROF count disagrees with configured TF timing");
    }
    for (int ir = 0; ir < int(rofs.size()); ++ir) {
      const auto& rof = rofs[ir];
      for (int local = rof.getFirstEntry(); local < rof.getFirstEntry() + rof.getNEntries(); ++local) {
        if (local < 0 || local >= mInputITSclusters.getNClusters(slot)) {
          throw std::runtime_error("Invalid extension/short tracking ROF cluster range");
        }
        const int ref = o2::itsmft::composeClusID(slot, local);
        const int flat = mInputITSclusters.flatIndex(ref);
        const auto& cl = mInputITSclusters[ref];
        const int layer = geom->getLayer(cl.getSensorID());
        if (layer < 0 || layer >= 7 || (nSlots > 1 && layer != slot)) {
          throw std::runtime_error("Invalid ITS extension/short tracking layer mapping");
        }
        // ITS tracking uses the supplied ROF ordinal. Effective bias/delay is
        // applied by LayerTiming, not subtracted from ROFRecord acquisition IR.
        const int clock = ir;
        auto& h = hits[flat];
        if (h.ref >= 0) {
          throw std::runtime_error("Overlapping ITS extension/short tracking ROF ranges");
        }
        const auto p = geom->getMatrixT2G(cl.getSensorID())(o2::math_utils::Point3D<float>{cl.getX(), cl.getY(), cl.getZ()});
        const auto bounds = timing[layer].getROFTimeBounds(clock, true);
        h = {ref, ir, clock, layer, p.X(), p.Y(), p.Z(), double(bounds.lower()), double(bounds.upper())};
        const float r = std::hypot(h.x, h.y);
        radiusMin[layer] = std::min(radiusMin[layer], r);
        radiusMax[layer] = std::max(radiusMax[layer], r);
        maxErrY[layer] = std::max(maxErrY[layer], cl.getSigmaY2());
        maxErrZ[layer] = std::max(maxErrZ[layer], cl.getSigmaZ2());
        const float phi = o2::its::math_utils::getNormalizedPhi(std::atan2(h.y, h.x));
        const int pb = std::clamp(utils.getPhiBinIndex(phi), 0, utils.getNphiBins() - 1);
        const int zb = std::clamp(utils.getZBinIndex(layer, h.z), 0, utils.getNzBins() - 1);
        grid[layer][{clock, pb, zb}].push_back(flat);
      }
    }
  }
  for (const auto& h : hits) {
    if (h.ref < 0) {
      throw std::runtime_error("ITS extension/short tracking input cluster has no ROF");
    }
  }
  std::vector<uint8_t> used(hits.size(), 0);
  auto mark = [&](int ref) {
    const int slot = o2::itsmft::clusID2Layer(ref), index = o2::itsmft::clusID2Index(ref);
    if (ref < 0 || slot >= nSlots || index >= mInputITSclusters.getNClusters(slot)) {
      throw std::runtime_error("Invalid upstream ITS extension/short tracking cluster reference");
    }
    used[mInputITSclusters.flatIndex(ref)] = 1;
  };
  for (const auto& tr : mInputITStracks) {
    if (tr.getFirstClusterEntry() < 0 || tr.getFirstClusterEntry() + tr.getNumberOfClusters() > mInputITSidxs.size()) {
      throw std::runtime_error("Invalid ITS track cluster span");
    }
    for (int i = 0; i < tr.getNumberOfClusters(); ++i) {
      mark(mInputITSidxs[tr.getFirstClusterEntry() + i]);
    }
  }
  if (data.hasITSABClusterUsage()) {
    const auto refs = data.getITSABClusterRefs();
    for (const auto& ab : data.getITSABRefs()) {
      if (ab.getFirstEntry() + ab.getNClusters() > refs.size()) {
        throw std::runtime_error("Invalid afterburner cluster span");
      }
      for (int i = 0; i < ab.getNClusters(); ++i) {
        mark(refs[ab.getFirstEntry() + i]);
      }
    }
  }
  const auto vertices = data.getPrimaryVertices();
  std::vector<std::pair<double, double>> pvWindows(vertices.size(), {0., 0.});
  for (int pv = 0; pv < vertices.size(); ++pv) {
    const auto& v = vertices[pv];
    if (v.getIRMin().isDummy() || v.getIRMax().isDummy() || v.getIRMax() < v.getIRMin()) {
      ++mStrangenessTrackingStats.nInvalidPVTime;
      continue;
    }
    const double margin = mPVTimeMargin / o2::constants::lhc::LHCBunchSpacingMUS;
    pvWindows[pv] = {v.getIRMin().differenceInBC(data.startIR) - 0.5 - margin,
                     v.getIRMax().differenceInBC(data.startIR) + 0.5 + margin};
  }
  auto pvTime = [&](int pv) { return pvWindows[pv]; };
  auto query = [&](int layer, float phi, float z, float dphi, float dz, double lo, double hi) {
    std::vector<int> result;
    if (!(hi > lo) || dphi < 0 || dz < 0) {
      return result;
    }
    const auto rr = timing[layer].getROFRange(float(lo), float(hi));
    auto rect = o2::its::getBinsRect(layer, phi, z, dz, dphi, utils);
    if (dphi >= float(M_PI)) {
      rect.y = 0;
      rect.w = utils.getNphiBins() - 1;
    }
    if (rect.x < 0) {
      return result;
    }
    int p = rect.y;
    do {
      for (int clock = rr.getFirstEntry(); clock < rr.getFirstEntry() + rr.getEntries(); ++clock) {
        auto it = grid[layer].lower_bound({clock, p, rect.x});
        const auto end = grid[layer].upper_bound({clock, p, rect.z});
        for (; it != end; ++it) {
          for (int flat : it->second) {
            if (!used[flat] && hits[flat].lower < hi && hits[flat].upper > lo) {
              result.push_back(flat);
            }
          }
        }
      }
      if (p == rect.w) {
        break;
      }
      p = (p + 1) % utils.getNphiBins();
    } while (true);
    std::sort(result.begin(), result.end(), [&](int a, int b) { return hits[a].ref < hits[b].ref; });
    return result;
  };
  auto attach = [&](ClusterCandidate& out, int flat) {
    const auto& h = hits[flat];
    out.mClusterRefs[h.layer] = h.ref;
    out.mROFRefs[h.layer] = h.rof;
    out.mClusterSizes[h.layer] = std::clamp(mInputClusterSizes[flat], 0, 255);
    out.mTimeMinBC = std::max(out.mTimeMinBC, h.lower);
    out.mTimeMaxBC = std::min(out.mTimeMaxBC, h.upper);
  };
  auto advance = [&](Track& tr, float x) {
    return o2::base::Propagator::Instance()->propagateToX(tr, x, mBz, o2::base::PropagatorF::MAX_SIN_PHI, o2::base::PropagatorF::MAX_STEP, mCorrType);
  };
  struct Proposal {
    ClusterCandidate track;
    std::array<uint32_t, 4> key;
  };
  std::vector<Proposal> proposals;
  std::set<std::array<uint32_t, 4>> matched;
  std::vector<ExtensionSeed> seeds;
  for (int thread = 0; thread < mNThreads; ++thread) {
    for (auto seed : mExtensionSeeds[thread]) {
      // Standalone IDs are already final. Embedded IDs were remapped and the
      // accepted-association result captured before the legacy output merge.
      if (!seed.matched && !mExtensionSeedsRemapped) {
        for (const auto& accepted : mStrangeTrackVec[thread]) {
          seed.matched |= accepted.mPartType == seed.source.mPartType && accepted.mDecayRef == seed.source.mDecayRef;
        }
      }
      if (seed.matched) {
        matched.insert(seed.key);
      }
      seeds.push_back(seed);
    }
  }
  mStrangenessTrackingStats.nSeeds = seeds.size();
  std::sort(seeds.begin(), seeds.end(), [](const auto& a, const auto& b) { return std::tie(a.key, a.pv, a.source.mDecayRef) < std::tie(b.key, b.pv, b.source.mDecayRef); });
  struct FollowBranch : ClusterCandidate {
    Track lastMeasured;
  };
  for (const auto& seed : seeds) {
    if (matched.count(seed.key) || seed.pv < 0 || seed.pv >= vertices.size() || !valid(seed.source.mMother)) {
      continue;
    }
    ++mStrangenessTrackingStats.nUnmatched;
    const auto [lo, hi] = pvTime(seed.pv);
    if (!(hi > lo)) {
      continue;
    }
    FollowBranch initial;
    initial.source = seed.source;
    initial.mTrack = seed.source.mMother;
    initial.lastMeasured = seed.source.mMother;
    initial.mPVRef = seed.pv;
    initial.mDecayRef = seed.source.mDecayRef;
    initial.mPartType = seed.source.mPartType;
    initial.mTimeMinBC = lo;
    initial.mTimeMaxBC = hi;
    std::vector<FollowBranch> branches{initial};
    const auto xyz = seed.source.mMother.getXYZGlo();
    const float decayR = std::hypot(xyz.X(), xyz.Y());
    for (int layer = 6; layer >= 0 && !branches.empty(); --layer) {
      if (radiusMin[layer] >= decayR) {
        continue;
      }
      std::vector<FollowBranch> next;
      for (const auto& branch : branches) {
        auto projected = branch.mTrack;
        float x;
        const auto here = projected.getXYZGlo();
        const float targetR = std::min(itsParams.LayerRadii[layer], std::hypot(here.X(), here.Y()) - 1.e-4f);
        if (!projected.getXatLabR(targetR, x, mBz, o2::track::DirInward) || !advance(projected, x)) {
          ++mStrangenessTrackingStats.nPropagationFailures;
          continue; // propagation failure is not a missing hit
        }
        if (!valid(projected)) {
          ++mStrangenessTrackingStats.nPropagationFailures;
          continue;
        }
        const auto p = projected.getXYZGlo();
        if (std::abs(p.Z()) > itsParams.LayerZ[layer]) {
          next.push_back(branch); // no active layer crossing
          continue;
        }
        const float r = std::hypot(p.X(), p.Y());
        if (!(r > 0)) {
          ++mStrangenessTrackingStats.nPropagationFailures;
          continue;
        }
        const float dr = radiusMax[layer] - radiusMin[layer];
        const float csp = std::sqrt(std::max(1.e-6f, 1.f - projected.getSnp() * projected.getSnp()));
        const float dy = cfg.clusterSearchNSigma * std::sqrt(projected.getSigmaY2() + maxErrY[layer]) + dr * (1.f + std::abs(projected.getSnp()) / csp);
        float dz = cfg.clusterSearchNSigma * std::sqrt(projected.getSigmaZ2() + maxErrZ[layer]) + dr * std::abs(projected.getTgl()) / csp;
        float dphi = std::min(float(M_PI), dy / r);
        const float ca = std::cos(projected.getAlpha()), sa = std::sin(projected.getAlpha());
        const float radialCos = std::abs((p.X() * (ca * csp - sa * projected.getSnp()) +
                                          p.Y() * (sa * csp + ca * projected.getSnp())) /
                                         r);
        // Local-frame csp is not the incidence cosine at a cylinder. Bound the
        // radius-to-sensor displacement and covariance projection with the
        // radial incidence instead. Sums of standard deviations also cover
        // unknown YZ correlations. Near a turning/grazing crossing, query the
        // full layer rather than let a singular linearization exclude hits.
        const float radialVariation = dr * (std::abs(projected.getCurvature(mBz)) + 1.f / radiusMin[layer]);
        if (radialCos < 1.e-3f || radialCos * radialCos <= 4.f * radialVariation) {
          dphi = float(M_PI);
          dz = std::max(dz, itsParams.LayerZ[layer] + std::abs(p.Z()));
        } else {
          const float sigmaY = std::sqrt(projected.getSigmaY2()) + std::sqrt(maxErrY[layer]);
          const float sigmaZ = std::sqrt(projected.getSigmaZ2()) + std::sqrt(maxErrZ[layer]);
          const float transverseSpan = (cfg.clusterSearchNSigma * sigmaY + dr) / radialCos;
          dphi = std::min(float(M_PI), std::max(dphi, transverseSpan / r));
          dz = std::max(dz, cfg.clusterSearchNSigma * sigmaZ + std::abs(projected.getTgl()) * transverseSpan);
        }
        const auto candidates = query(layer, std::atan2(p.Y(), p.X()), p.Z(), dphi, dz, branch.mTimeMinBC, branch.mTimeMaxBC);
        for (int flat : candidates) {
          auto candidate = branch;
          auto& tr = candidate.mTrack;
          const auto& cl = mInputITSclusters[hits[flat].ref];
          if (!tr.rotate(geom->getSensorRefAlpha(cl.getSensorID()))) {
            ++mStrangenessTrackingStats.nPropagationFailures;
            continue;
          }
          if (cl.getX() > tr.getX()) {
            continue;
          } // only measurements crossed inward of the branch
          if (!advance(tr, cl.getX())) {
            ++mStrangenessTrackingStats.nPropagationFailures;
            continue;
          }
          if (mCorrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE &&
              !tr.correctForMaterial(itsParams.LayerxX0[layer], itsParams.LayerxX0[layer] * o2::its::constants::Radl * o2::its::constants::Rho, true)) {
            ++mStrangenessTrackingStats.nPropagationFailures;
            continue;
          }
          const std::array<float, 2> pos{cl.getY(), cl.getZ()};
          const std::array<float, 3> cov{cl.getSigmaY2(), cl.getSigmaYZ(), cl.getSigmaZ2()};
          const float chi2 = tr.getPredictedChi2(pos, cov);
          if (chi2 < 0 || chi2 >= cfg.extensionMaxHitChi2 || !tr.update(pos, cov) || !valid(tr)) {
            continue;
          }
          attach(candidate, flat);
          candidate.mHitChi2 += chi2;
          candidate.lastMeasured = tr;
          next.push_back(candidate);
        }
        // Retain a missing-hit alternative even when a hit alternative exists.
        if (holes(branch) < cfg.extensionMaxMissing) {
          auto missing = branch;
          if (mCorrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE &&
              !projected.correctForMaterial(itsParams.LayerxX0[layer], itsParams.LayerxX0[layer] * o2::its::constants::Radl * o2::its::constants::Rho, true)) {
            ++mStrangenessTrackingStats.nPropagationFailures;
            continue;
          }
          if (!valid(projected)) {
            ++mStrangenessTrackingStats.nPropagationFailures;
            continue;
          }
          missing.mTrack = projected;
          missing.mMissingLayerMask |= 1 << layer;
          next.push_back(missing);
          ++mStrangenessTrackingStats.nMissingHits;
        }
      }
      std::sort(next.begin(), next.end(), better);
      if (next.size() > cfg.extensionMaxBranches) {
        mStrangenessTrackingStats.nMotherTruncated += next.size() - cfg.extensionMaxBranches;
        next.resize(cfg.extensionMaxBranches);
      }
      branches.swap(next);
    }
    for (auto& tr : branches) {
      if (tr.getNClusters() < cfg.extensionMinHits || !(tr.mTimeMaxBC > tr.mTimeMinBC)) {
        continue;
      }
      tr.mTrack = tr.lastMeasured;     // persist the innermost successful measurement surface
      tr.mNDF = 2 * tr.getNClusters(); // innovations conditioned on the decay prior
      proposals.push_back({tr, seed.key});
    }
  }
  auto commit = [&](const ClusterCandidate& tr) {
    for (int ref : tr.mClusterRefs) {
      if (ref >= 0) {
        mark(ref);
      }
    }
    o2::MCCompLabel best;
    if (mMCTruthON) {
      std::map<o2::MCCompLabel, int> counts;
      for (int ref : tr.mClusterRefs) {
        if (ref < 0) {
          continue;
        }
        const auto* labels = data.getITSClustersMCLabels(o2::itsmft::clusID2Layer(ref));
        if (!labels) {
          throw std::runtime_error("Missing per-stream extension/short tracking cluster MC labels");
        }
        std::set<o2::MCCompLabel> unique;
        for (auto label : labels->getLabels(o2::itsmft::clusID2Index(ref))) {
          if (label.isValid()) {
            label.setFakeFlag(false);
            unique.insert(label);
          }
        }
        for (const auto& label : unique) {
          ++counts[label];
        }
      }
      int n = 0;
      for (const auto& [label, count] : counts) {
        if (count > n) {
          n = count;
          best = label;
        }
      }
      if (best.isValid()) {
        best.setFakeFlag(n != tr.getNClusters());
      }
    }
    if (!tr.isShort) {
      auto output = tr.source;
      output.mMother = tr.mTrack;
      output.mITSRef = -1u;
      output.mMatchChi2 = -1.f;
      output.mTopoChi2 = tr.mHitChi2;
      output.mClusterSizes = 0;
      for (int layer = 0; layer < 7; ++layer) {
        output.setClusterSize(layer, tr.mClusterSizes[layer]);
      }
      mExtendedStrangeTracks.push_back(output);
      mExtendedStrangeTimes.emplace_back(float(0.5 * (tr.mTimeMinBC + tr.mTimeMaxBC)), float(0.5 * (tr.mTimeMaxBC - tr.mTimeMinBC)));
      if (mMCTruthON) {
        mExtendedStrangeLabels.push_back(best);
      }
    } else {
      TrackITS output(tr.mTrack, tr.mHitChi2, tr.mOuter);
      output.setClusterRefs(mITSShortClusterRefs.size(), 3);
      output.setPattern(0x7);
      output.getTimeStamp() = o2::its::TimeStamp{float(0.5 * (tr.mTimeMinBC + tr.mTimeMaxBC)), float(0.5 * (tr.mTimeMaxBC - tr.mTimeMinBC))};
      for (int layer = 0; layer < 3; ++layer) {
        output.setClusterSize(layer, tr.mClusterSizes[layer]);
        mITSShortClusterRefs.push_back(tr.mClusterRefs[layer]);
        if (mMCTruthON && best.isValid()) {
          const int ref = tr.mClusterRefs[layer];
          const auto* labels = data.getITSClustersMCLabels(o2::itsmft::clusID2Layer(ref));
          for (auto label : labels->getLabels(o2::itsmft::clusID2Index(ref))) {
            label.setFakeFlag(false);
            auto selected = best;
            selected.setFakeFlag(false);
            if (label == selected) {
              output.setPattern(output.getPattern() | (1u << (16 + layer)));
              break;
            }
          }
        }
      }
      mITSShortTracks.push_back(output);
      mITSShortPVRefs.push_back(tr.mPVRef);
      if (mMCTruthON) {
        mITSShortLabels.push_back(best);
      }
    }
  };
  auto free = [&](const ClusterCandidate& tr) {
    for (int ref : tr.mClusterRefs) {
      if (ref >= 0 && used[mInputITSclusters.flatIndex(ref)]) {
        return false;
      }
    }
    return true;
  };
  std::sort(proposals.begin(), proposals.end(), [](const auto& a, const auto& b) {
    if (a.track.getNClusters() != b.track.getNClusters()) {
      return a.track.getNClusters() > b.track.getNClusters();
    }
    if (holes(a.track) != holes(b.track)) {
      return holes(a.track) < holes(b.track);
    }
    const float qa = a.track.mHitChi2 / a.track.getNClusters(), qb = b.track.mHitChi2 / b.track.getNClusters();
    if (qa != qb) {
      return qa < qb;
    }
    return std::tie(a.key, a.track.mPVRef, a.track.mClusterRefs, a.track.mPartType, a.track.mDecayRef) < std::tie(b.key, b.track.mPVRef, b.track.mClusterRefs, b.track.mPartType, b.track.mDecayRef);
  });
  std::set<std::array<uint32_t, 4>> selected;
  for (const auto& proposal : proposals) {
    if (!selected.count(proposal.key) && free(proposal.track)) {
      commit(proposal.track);
      selected.insert(proposal.key);
    }
  }
  if (!cfg.enableShortIB || mStrangenessTrackingStats.tripletsDeferredZeroField) {
    return;
  }
  std::vector<ClusterCandidate> triplets;
  const std::array<float, 3> material{itsParams.LayerxX0[0], itsParams.LayerxX0[1], itsParams.LayerxX0[2]};
  // Build spatially bounded PV-pointing pairs before fitting exactly 0-1-2.
  // The low pT floor bounds curvature only; momentum never orders proposals.
  for (int pv = 0; pv < vertices.size(); ++pv) {
    const auto [lo, hi] = pvTime(pv);
    if (!(hi > lo)) {
      continue;
    }
    int perPV = 0;
    bool capped = false;
    const auto inner = query(0, 0.f, 0.f, float(M_PI), itsParams.LayerZ[0], lo, hi);
    for (int flat0 : inner) {
      if (capped) {
        break;
      }
      const auto& h0 = hits[flat0];
      if (used[flat0] || h0.layer != 0 || h0.upper <= lo || h0.lower >= hi) {
        continue;
      }
      const auto& vertex = vertices[pv];
      const float dx = h0.x - vertex.getX(), dy = h0.y - vertex.getY();
      const float r0 = std::hypot(dx, dy);
      if (r0 < 1.e-6f) {
        continue;
      }
      const float ux = dx / r0, uy = dy / r0;
      const float dot = vertex.getX() * ux + vertex.getY() * uy;
      const float pvR2 = vertex.getX() * vertex.getX() + vertex.getY() * vertex.getY();
      const float slope = (h0.z - vertex.getZ()) / r0;
      std::array<std::vector<int>, 2> outer;
      for (int layer = 1; layer <= 2; ++layer) {
        const float r = itsParams.LayerRadii[layer];
        const float discriminant = dot * dot + r * r - pvR2;
        if (discriminant <= 0) {
          continue;
        }
        const float distance = -dot + std::sqrt(discriminant);
        if (distance <= r0) {
          continue;
        }
        const float phi = std::atan2(vertex.getY() + distance * uy, vertex.getX() + distance * ux);
        const float curvature = std::abs(mBz * o2::constants::math::B2C) / cfg.shortTrackMinPt;
        float dphi = std::min(float(M_PI), curvature * (r + radiusMax[0]) + cfg.clusterSearchNSigma * std::sqrt(maxErrY[0] + maxErrY[layer] + vertex.getSigmaX2() + vertex.getSigmaY2()) / r0);
        float dz = cfg.clusterSearchNSigma * std::sqrt(maxErrZ[layer] + maxErrZ[0] * r * r / (r0 * r0) + vertex.getSigmaZ2()) +
                   std::abs(slope) * (radiusMax[layer] - radiusMin[layer] + curvature * curvature * r * r * r);
        const float ratio = distance / r0;
        const float pvWeight = std::abs(1.f - ratio);
        const float rayCos = std::abs((dot + distance) / r);
        // Extrapolation is PV + ratio*(L0-PV), not r/r0 for a displaced PV.
        // PV XY uncertainty changes the ray length and therefore Z at nonzero
        // slope. Sum-of-SD bounds include that term and arbitrary PV/cluster
        // coordinate correlations; these widths only enlarge candidate lookup.
        if (rayCos < 1.e-3f) {
          dphi = float(M_PI);
          dz = std::max(dz, itsParams.LayerZ[layer] + std::abs(vertex.getZ() + slope * distance));
        } else {
          const float sigmaPVXY = std::sqrt(vertex.getSigmaX2() + vertex.getSigmaY2());
          const float transverseSigma = std::sqrt(maxErrY[layer]) + ratio * std::sqrt(maxErrY[0]) + pvWeight * sigmaPVXY;
          const float radialSpan = radiusMax[layer] - radiusMin[layer];
          dphi = std::min(float(M_PI), std::max(dphi, curvature * (r + radiusMax[0]) +
                                                        (cfg.clusterSearchNSigma * transverseSigma + radialSpan) / (r * rayCos)));
          const float sigmaZ = std::sqrt(maxErrZ[layer]) + ratio * std::sqrt(maxErrZ[0]) + pvWeight * std::sqrt(vertex.getSigmaZ2()) +
                               std::abs(slope) * transverseSigma / rayCos;
          dz = std::max(dz, cfg.clusterSearchNSigma * sigmaZ + std::abs(slope) *
                                                                 (radialSpan / rayCos + curvature * curvature * r * r * r));
        }
        outer[layer - 1] = query(layer, phi, vertex.getZ() + slope * distance, dphi, dz,
                                 std::max(lo, h0.lower), std::min(hi, h0.upper));
      }
      for (int flat1 : outer[0]) {
        if (capped) {
          break;
        }
        for (int flat2 : outer[1]) {
          const auto& h1 = hits[flat1];
          const auto& h2 = hits[flat2];
          const double lower = std::max({lo, h0.lower, h1.lower, h2.lower});
          const double upper = std::min({hi, h0.upper, h1.upper, h2.upper});
          if (!(upper > lower)) {
            continue;
          }
          if (perPV >= cfg.shortTrackMaxPerPV || mStrangenessTrackingStats.nTripletProposals >= cfg.shortTrackMaxPerTF) {
            ++mStrangenessTrackingStats.nTripletTruncated;
            capped = true;
            break;
          }
          ++perPV;
          ++mStrangenessTrackingStats.nTripletProposals;
          ClusterCandidate tr;
          tr.isShort = true;
          tr.mPVRef = pv;
          tr.mTimeMinBC = lower;
          tr.mTimeMaxBC = upper;
          const std::array<ITSCluster, 3> clusters{mInputITSclusters[h0.ref], mInputITSclusters[h1.ref], mInputITSclusters[h2.ref]};
          std::array<float, 2> preDCA;
          if (!fitITSShortTrack(clusters, vertex, mBz, mCorrType, material, cfg.shortTrackMinPt, cfg.shortTrackMaxPointingChi2,
                                cfg.shortTrackMaxChi2NDF, tr.mTrack, tr.mHitChi2, tr.mPVChi2, tr.mNDF, &preDCA, &tr.mOuter)) {
            continue;
          }
          if (!valid(tr.mTrack)) {
            continue;
          }
          attach(tr, flat0);
          attach(tr, flat1);
          attach(tr, flat2);
          triplets.push_back(tr);
        }
      }
    }
  }
  std::sort(triplets.begin(), triplets.end(), [](const auto& a, const auto& b) {
    const auto qa = a.mHitChi2 / a.mNDF, qb = b.mHitChi2 / b.mNDF;
    if (qa != qb) {
      return qa < qb;
    }
    if (a.mPVChi2 != b.mPVChi2) {
      return a.mPVChi2 < b.mPVChi2;
    }
    return std::tie(a.mClusterRefs, a.mPVRef) < std::tie(b.mClusterRefs, b.mPVRef);
  });
  for (const auto& tr : triplets) {
    if (free(tr)) {
      commit(tr);
    }
  }
}
} // namespace o2::strangeness_tracking

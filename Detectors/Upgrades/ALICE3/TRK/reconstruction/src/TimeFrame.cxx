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
///
/// \file TimeFrame.cxx
/// \brief TRK TimeFrame implementation
///

#include "TRKReconstruction/TimeFrame.h"
#include "TRKReconstruction/Clusterer.h"
#include "TRKSimulation/Hit.h"
#include "TRKBase/GeometryTGeo.h"
#include "TRKBase/SegmentationChip.h"
#include "Framework/Logger.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/DigitizationContext.h"
#include "Steer/MCKinematicsReader.h"
#include <TTree.h>
#include <TRandom3.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <array>

namespace o2::trk
{

template <int nLayers>
int TimeFrame<nLayers>::loadROFsFromHitTree(TTree* hitsTree, GeometryTGeo* gman, const nlohmann::json& config)
{
  constexpr std::array<int, 2> startLayer{0, 3};
  const Long64_t nEvents = hitsTree->GetEntries();

  gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L) | o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));

  std::vector<o2::trk::Hit>* trkHit = nullptr;
  hitsTree->SetBranchAddress("TRKHit", &trkHit);

  const int inROFpileup{config.contains("inROFpileup") ? config["inROFpileup"].get<int>() : 1};

  // Calculate number of ROFs and initialize data structures
  this->mNrof = (nEvents + inROFpileup - 1) / inROFpileup;

  // Reset and prepare ROF data structures
  for (int iLayer{0}; iLayer < nLayers; ++iLayer) {
    this->mMinR[iLayer] = std::numeric_limits<float>::max();
    this->mMaxR[iLayer] = std::numeric_limits<float>::lowest();
    this->mROFramesClusters[iLayer].clear();
    this->mROFramesClusters[iLayer].resize(this->mNrof + 1, 0);
    this->mUnsortedClusters[iLayer].clear();
    this->mTrackingFrameInfo[iLayer].clear();
    this->mClusterExternalIndices[iLayer].clear();
  }

  // Pre-count hits to reserve memory efficiently
  int totalNHits{0};
  std::array<int, nLayers> clusterCountPerLayer{};
  for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) {
    hitsTree->GetEntry(iEvent);
    for (const auto& hit : *trkHit) {
      if (gman->getDisk(hit.GetDetectorID()) != -1) {
        continue; // skip non-barrel hits
      }
      int subDetID = gman->getSubDetID(hit.GetDetectorID());
      const int layer = startLayer[subDetID] + gman->getLayer(hit.GetDetectorID());
      if (layer >= nLayers) {
        continue;
      }
      ++clusterCountPerLayer[layer];
      totalNHits++;
    }
  }

  // Reserve memory for all layers
  for (int iLayer{0}; iLayer < nLayers; ++iLayer) {
    this->mUnsortedClusters[iLayer].reserve(clusterCountPerLayer[iLayer]);
    this->mTrackingFrameInfo[iLayer].reserve(clusterCountPerLayer[iLayer]);
    this->mClusterExternalIndices[iLayer].reserve(clusterCountPerLayer[iLayer]);
  }
  clearResizeBoundedVector(this->mClusterSize, totalNHits, this->mMemoryPool.get());

  std::array<float, 11> resolution{0.001, 0.001, 0.001, 0.001, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004};
  if (config["geometry"]["pitch"].size() == nLayers) {
    for (int iLayer{0}; iLayer < config["geometry"]["pitch"].size(); ++iLayer) {
      LOGP(info, "Setting resolution for layer {} from config", iLayer);
      LOGP(info, "Layer {} pitch {} cm", iLayer, config["geometry"]["pitch"][iLayer].get<float>());
      resolution[iLayer] = config["geometry"]["pitch"][iLayer].get<float>() / std::sqrt(12.f);
    }
  }
  LOGP(info, "Number of active parts in VD: {}", gman->getNumberOfActivePartsVD());

  int hitCounter{0};
  auto labels = new dataformats::MCTruthContainer<MCCompLabel>();

  int iRof{0}; // Current ROF index
  for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) {
    hitsTree->GetEntry(iEvent);

    for (auto& hit : *trkHit) {
      if (gman->getDisk(hit.GetDetectorID()) != -1) {
        continue; // skip non-barrel hits for this test
      }
      int subDetID = gman->getSubDetID(hit.GetDetectorID());
      const int layer = startLayer[subDetID] + gman->getLayer(hit.GetDetectorID());

      float alpha{0.f};
      o2::math_utils::Point3D<float> gloXYZ;
      o2::math_utils::Point3D<float> trkXYZ;
      float r{0.f};
      if (layer >= nLayers) {
        continue;
      }
      if (layer >= 3) {
        int chipID = hit.GetDetectorID();
        alpha = gman->getSensorRefAlphaMLOT(chipID);
        const o2::math_utils::Transform3D& l2g = gman->getMatrixL2G(chipID);
        auto locXYZ = l2g ^ (hit.GetPos());
        locXYZ.SetX(locXYZ.X() + gRandom->Gaus(0.0, resolution[layer]));
        locXYZ.SetZ(locXYZ.Z() + gRandom->Gaus(0.0, resolution[layer]));
        gloXYZ = gman->getMatrixL2G(chipID) * locXYZ;
        trkXYZ = gman->getMatrixT2L(chipID - gman->getNumberOfActivePartsVD()) ^ locXYZ;
        r = std::hypot(gloXYZ.X(), gloXYZ.Y());
      } else {
        const auto& hitPos = hit.GetPos();
        r = std::hypot(hitPos.X(), hitPos.Y());
        alpha = std::atan2(hitPos.Y(), hitPos.X()) + gRandom->Gaus(0.0, resolution[layer] / r);
        o2::math_utils::bringTo02Pi(alpha);
        gloXYZ.SetX(r * std::cos(alpha));
        gloXYZ.SetY(r * std::sin(alpha));
        gloXYZ.SetZ(hitPos.Z() + gRandom->Gaus(0.0, resolution[layer]));
        trkXYZ.SetX(r);
        trkXYZ.SetY(0.f);
        trkXYZ.SetZ(gloXYZ.Z());
      }
      this->mMinR[layer] = std::min(this->mMinR[layer], r);
      this->mMaxR[layer] = std::max(this->mMaxR[layer], r);
      this->addTrackingFrameInfoToLayer(layer, gloXYZ.x(), gloXYZ.y(), gloXYZ.z(), trkXYZ.x(), alpha,
                                        std::array<float, 2>{trkXYZ.y(), trkXYZ.z()},
                                        std::array<float, 3>{resolution[layer] * resolution[layer], 0., resolution[layer] * resolution[layer]});
      /// Rotate to the global frame
      this->addClusterToLayer(layer, gloXYZ.x(), gloXYZ.y(), gloXYZ.z(), this->mUnsortedClusters[layer].size());
      this->addClusterExternalIndexToLayer(layer, hitCounter);
      MCCompLabel label{hit.GetTrackID(), static_cast<int>(iEvent), 0};
      labels->addElement(hitCounter, label);
      this->mClusterSize[hitCounter] = 1; // For compatibility with cluster-based tracking, set cluster size to 1 for hits
      hitCounter++;
    }
    trkHit->clear();

    // Update ROF structure when we complete an ROF or reach the last event
    if ((iEvent + 1) % inROFpileup == 0 || iEvent == nEvents - 1) {
      iRof++;
      for (unsigned int iLayer{0}; iLayer < this->mUnsortedClusters.size(); ++iLayer) {
        this->mROFramesClusters[iLayer][iRof] = this->mUnsortedClusters[iLayer].size(); // effectively calculating an exclusive sum
      }
      // Update primary vertices ROF structure
    }
    this->mClusterLabels = labels;
  }
  return this->mNrof;
}

template <int nLayers>
int TimeFrame<nLayers>::loadROFrameData(gsl::span<const o2::trk::ROFRecord> rofs,
                                        gsl::span<const o2::trk::Cluster> clusters,
                                        gsl::span<const unsigned char> patterns,
                                        const dataformats::MCTruthContainer<MCCompLabel>* mcLabels,
                                        float yPlaneMLOT)
{
  constexpr std::array<int, 2> startLayer{0, 3};
  GeometryTGeo* geom = GeometryTGeo::Instance();
  geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L) | o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));

  this->mNrof = rofs.size();

  for (int iLayer{0}; iLayer < nLayers; ++iLayer) {
    this->mMinR[iLayer] = std::numeric_limits<float>::max();
    this->mMaxR[iLayer] = std::numeric_limits<float>::lowest();
    this->mROFramesClusters[iLayer].clear();
    this->mROFramesClusters[iLayer].resize(this->mNrof + 1, 0);
    this->mUnsortedClusters[iLayer].clear();
    this->mTrackingFrameInfo[iLayer].clear();
    this->mClusterExternalIndices[iLayer].clear();
  }

  std::array<int, nLayers> clusterCountPerLayer{};
  for (const auto& c : clusters) {
    if (c.subDetID < 0 || c.subDetID > 1) {
      continue;
    }
    if (c.disk != -1) {
      continue; // skip non-barrel clusters for now
    }
    const int layer = startLayer[c.subDetID] + c.layer;
    if (layer < 0 || layer >= nLayers) {
      continue;
    }
    ++clusterCountPerLayer[layer];
  }

  for (int iLayer{0}; iLayer < nLayers; ++iLayer) {
    this->mUnsortedClusters[iLayer].reserve(clusterCountPerLayer[iLayer]);
    this->mTrackingFrameInfo[iLayer].reserve(clusterCountPerLayer[iLayer]);
    this->mClusterExternalIndices[iLayer].reserve(clusterCountPerLayer[iLayer]);
  }
  clearResizeBoundedVector(this->mClusterSize, clusters.size(), this->mMemoryPool.get());

  const uint8_t* pattPtr = patterns.data();
  const uint8_t* pattEnd = pattPtr + patterns.size();

  for (size_t iRof{0}; iRof < rofs.size(); ++iRof) {
    const auto& rof = rofs[iRof];
    const int first = rof.getFirstEntry();
    const int last = first + rof.getNEntries();

    for (int clusterId{first}; clusterId < last; ++clusterId) {
      if (clusterId < 0 || clusterId >= static_cast<int>(clusters.size())) {
        LOGP(warning, "Skipping out-of-range cluster id {} for ROF {}", clusterId, iRof);
        continue;
      }

      const auto& c = clusters[clusterId];
      if (c.subDetID < 0 || c.subDetID > 1 || c.disk != -1) {
        continue;
      }

      if (pattPtr + 2 > pattEnd) {
        LOGP(error, "Pattern stream exhausted while decoding cluster {}", clusterId);
        break;
      }
      const uint8_t* pattForCluster = pattPtr;
      const int nBytes = (pattForCluster[0] * pattForCluster[1] + 7) / 8;
      if (pattPtr + 2 + nBytes > pattEnd) {
        LOGP(error, "Pattern stream truncated for cluster {}", clusterId);
        break;
      }

      const int layer = startLayer[c.subDetID] + c.layer;
      if (layer < 0 || layer >= nLayers) {
        LOGP(error, "Skipping cluster with invalid layer {} (subDetID {}, layer {})", layer, c.subDetID, c.layer);
        pattPtr += 2 + nBytes;
        continue;
      }

      auto locXYZ = Clusterer::getClusterLocalCoordinates(c, pattForCluster, yPlaneMLOT);
      pattPtr += 2 + nBytes;

      const auto gloXYZ = geom->getMatrixL2G(c.chipID) * locXYZ;

      float alpha{0.f};
      o2::math_utils::Point3D<float> trkXYZ;
      if (c.subDetID == 1) {
        alpha = geom->getSensorRefAlphaMLOT(c.chipID);
        trkXYZ = geom->getMatrixT2L(c.chipID - geom->getNumberOfActivePartsVD()) ^ locXYZ;
      } else {
        const float r = std::hypot(gloXYZ.X(), gloXYZ.Y());
        alpha = std::atan2(gloXYZ.Y(), gloXYZ.X());
        o2::math_utils::bringTo02Pi(alpha);
        trkXYZ.SetX(r);
        trkXYZ.SetY(0.f);
        trkXYZ.SetZ(gloXYZ.Z());
      }

      const float r = std::hypot(gloXYZ.X(), gloXYZ.Y());
      this->mMinR[layer] = std::min(this->mMinR[layer], r);
      this->mMaxR[layer] = std::max(this->mMaxR[layer], r);

      const float sigmaY2 = (c.subDetID == 0)
                              ? 0.25f * SegmentationChip::PitchRowVD * SegmentationChip::PitchRowVD
                              : 0.25f * SegmentationChip::PitchRowMLOT * SegmentationChip::PitchRowMLOT;
      const float sigmaZ2 = (c.subDetID == 0)
                              ? 0.25f * SegmentationChip::PitchColVD * SegmentationChip::PitchColVD
                              : 0.25f * SegmentationChip::PitchColMLOT * SegmentationChip::PitchColMLOT;

      this->addTrackingFrameInfoToLayer(layer, gloXYZ.x(), gloXYZ.y(), gloXYZ.z(), trkXYZ.x(), alpha,
                                        std::array<float, 2>{trkXYZ.y(), trkXYZ.z()},
                                        std::array<float, 3>{sigmaY2, 0.f, sigmaZ2});
      this->addClusterToLayer(layer, gloXYZ.x(), gloXYZ.y(), gloXYZ.z(), this->mUnsortedClusters[layer].size());
      this->addClusterExternalIndexToLayer(layer, clusterId);
      this->mClusterSize[clusterId] = std::clamp(static_cast<unsigned int>(c.size), 0u, 255u);
    }

    for (unsigned int iL{0}; iL < this->mUnsortedClusters.size(); ++iL) {
      this->mROFramesClusters[iL][iRof + 1] = this->mUnsortedClusters[iL].size();
    }
  }

  for (auto i = 0; i < this->mNTrackletsPerCluster.size(); ++i) {
    this->mNTrackletsPerCluster[i].resize(this->mUnsortedClusters[1].size());
    this->mNTrackletsPerClusterSum[i].resize(this->mUnsortedClusters[1].size() + 1);
  }

  if (mcLabels != nullptr) {
    this->mClusterLabels = mcLabels;
  }

  return this->mNrof;
}

template <int nLayers>
void TimeFrame<nLayers>::getPrimaryVerticesFromMC(TTree* mcHeaderTree, int nRofs, Long64_t nEvents, int inROFpileup)
{
  auto mcheader = new o2::dataformats::MCEventHeader;
  mcHeaderTree->SetBranchAddress("MCEventHeader.", &mcheader);

  this->mROFramesPV.clear();
  this->mROFramesPV.resize(nRofs + 1, 0);
  this->mPrimaryVertices.clear();

  int iRof{0};
  for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) {
    mcHeaderTree->GetEntry(iEvent);
    o2::its::Vertex vertex;
    vertex.setXYZ(mcheader->GetX(), mcheader->GetY(), mcheader->GetZ());
    vertex.setNContributors(30);
    vertex.setChi2(0.f);
    LOGP(debug, "ROF {}: Added primary vertex at ({}, {}, {})", iRof, mcheader->GetX(), mcheader->GetY(), mcheader->GetZ());
    this->mPrimaryVertices.push_back(vertex);
    if ((iEvent + 1) % inROFpileup == 0 || iEvent == nEvents - 1) {
      iRof++;
      this->mROFramesPV[iRof] = this->mPrimaryVertices.size(); // effectively calculating an exclusive sum
    }
  }
  this->mMultiplicityCutMask.resize(nRofs, true); /// all ROFs are valid with MC primary vertices.
}



// Explicit template instantiation for TRK with 11 layers

template <int nLayers>
void TimeFrame<nLayers>::addTruthSeedingVertices(gsl::span<const o2::trk::ROFRecord> rofs)
{
  LOGP(info, "TRK: using truth seeds as vertices from DigitizationContext");
  this->resetRofPV();

  const auto dc = o2::steer::DigitizationContext::loadFromFile("collisioncontext.root");
  const auto irs = dc->getEventRecords();
  o2::steer::MCKinematicsReader mcReader(dc);

  // Pre-compute ROF start BC (as absolute long) for binary search
  std::vector<int64_t> rofStartBC(rofs.size());
  for (size_t i = 0; i < rofs.size(); ++i) {
    rofStartBC[i] = rofs[i].getBCData().toLong();
  }

  using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;
  struct VertInfo {
    std::pmr::vector<Vertex> vertices;
    std::pmr::vector<int> srcs;
    std::pmr::vector<int> events;
  };
  std::map<int, VertInfo> vertMap;

  const int iSrc = 0; // primary collision generator source
  auto eveId2colId = dc->getCollisionIndicesForSource(iSrc);
  for (int iEve{0}; iEve < mcReader.getNEvents(iSrc); ++iEve) {
    const auto& ir = irs[eveId2colId[iEve]];
    if (!ir.isDummy()) {
      const auto& eve = mcReader.getMCEventHeader(iSrc, iEve);
      const int64_t evBC = ir.toLong();
      // Find ROF: last ROF whose start BC <= evBC
      auto it = std::upper_bound(rofStartBC.begin(), rofStartBC.end(), evBC);
      if (it != rofStartBC.begin()) {
        --it;
        int rofId = static_cast<int>(std::distance(rofStartBC.begin(), it));
        auto* mr = this->mMemoryPool.get();
        if (!vertMap.contains(rofId)) {
          vertMap[rofId] = {
            .vertices = std::pmr::vector<Vertex>(mr),
            .srcs = std::pmr::vector<int>(mr),
            .events = std::pmr::vector<int>(mr),
          };
        }
        Vertex vert;
        vert.setTimeStamp(rofId);
        vert.setNContributors(std::max(1L, std::ranges::count_if(
                                            mcReader.getTracks(iSrc, iEve),
                                            [](const auto& trk) {
                                              return trk.isPrimary() && trk.GetPt() > 0.05 && std::abs(trk.GetEta()) < 1.1;
                                            })));
        vert.setXYZ((float)eve.GetX(), (float)eve.GetY(), (float)eve.GetZ());
        vert.setChi2(1);
        constexpr float cov = 50e-9f;
        vert.setCov(cov, cov, cov, cov, cov, cov);
        vertMap[rofId].vertices.push_back(vert);
        vertMap[rofId].srcs.push_back(iSrc);
        vertMap[rofId].events.push_back(iEve);
      }
    }
    mcReader.releaseTracksForSourceAndEvent(iSrc, iEve);
  }

  size_t nVerts{0};
  auto* mr = this->mMemoryPool.get();
  for (int iROF{0}; iROF < this->getNrof(); ++iROF) {
    std::pmr::vector<Vertex> verts(mr);
    std::pmr::vector<std::pair<o2::MCCompLabel, float>> polls(mr);
    if (vertMap.contains(iROF)) {
      const auto& info = vertMap[iROF];
      verts = info.vertices;
      nVerts += verts.size();
      for (size_t i{0}; i < verts.size(); ++i) {
        o2::MCCompLabel lbl(o2::MCCompLabel::maxTrackID(), info.events[i], info.srcs[i], false);
        polls.emplace_back(lbl, 1.f);
      }
    } else {
      this->getNoVertexROF()++;
    }
    this->addPrimaryVertices(verts, 0);
    this->addPrimaryVerticesLabels(polls);
  }
  this->mMultiplicityCutMask.resize(this->getNrof(), true);
  LOGP(info, "TRK truth seeding: {}/{} ROFs with {} vertices -> <NV>={:.2f}",
       vertMap.size(), this->getNrof(), nVerts,
       vertMap.size() > 0 ? (float)nVerts / (float)vertMap.size() : 0.f);
}

template class TimeFrame<11>;

} // namespace o2::trk

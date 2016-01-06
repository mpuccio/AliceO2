//-*- Mode: C++ -*-
// **************************************************************************
// This file is property of and copyright by the ALICE ITSU Project         *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Author: Maximiliano Puccio <maximiliano.puccio@cern.ch>          *
//                 for the ITS Upgrade project                              *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************

#include "CATracker.h"

// STD
#include <algorithm>
#include <cassert>
// ROOT
#include <TBranch.h>
#include <TMath.h>
#include <TTree.h>
#include <Riostream.h>
// ALIROOT
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliLog.h"
// ALIROOT ITSU
#include "Aux.h"
#include "AliITSUClusterPix.h"
#include "AliITSURecoDet.h"
#include "AliITSUReconstructor.h"
#include "AliITSUTrackCooked.h"

using AliceO2::ITS::CA;

using TMath::Abs;
using TMath::Sort;
using TMath::Sqrt;
using std::sort;
using std::cout;
using std::endl;
using std::flush;

// tolerance for layer on-surface check
const Double_t Tracker::mkChi2Cut =  600.f;
const int Tracker::mkNumberOfIterations =  2;
const float Tracker::mkR[7] = {2.33959,3.14076,3.91924,19.6213,24.5597,34.388,39.3329};
//
const float kmaxDCAxy[5] = {0.05f,0.04f,0.05f,0.2f,0.4f};
const float kmaxDCAz[5] = {0.2f,0.4f,0.5f,0.6f,3.f};
const float kmaxDN[4] = {0.002f,0.009f,0.002f,0.005f};
const float kmaxDP[4] = {0.008f,0.0025f,0.003f,0.0035f};
const float kmaxDZ[6] = {0.1f,0.1f,0.3f,0.3f,0.3f,0.3f};
const float kDoublTanL = 0.025;
const float kDoublPhi = 0.14;

const float kmaxDCAxy1[5] = /*{1.f,0.5,0.5,1.7,3.};/*/{1.f,0.4f,0.4f,1.5f,3.f};
const float kmaxDCAz1[5] = /*{2.f,0.8,0.8,3.,5.};/*/{1.f,0.4f,0.4f,1.5f,3.f};
const float kmaxDN1[4] = /*{0.006f,0.0045f,0.01f,0.04f};/*/{0.005f,0.0035f,0.009f,0.03f};
const float kmaxDP1[4] = /*{0.04f,0.01f,0.012f,0.014f};/*/{0.02f,0.005f,0.006f,0.007f};
const float kmaxDZ1[6] = /*{1.5f,1.5f,2.f,2.f,2.f,2.f};/*/{1.f,1.f,1.5f,1.5f,1.5f,1.5f};
const float kDoublTanL1 = /*0.12f;/*/0.05f;
const float kDoublPhi1 = /*0.4f;/*/0.2f;

Tracker::Tracker(AliITSUReconstructor* rec)
:mLayer()
,mUsedClusters()
,mChi2Cut(mkChi2Cut)
,mPhiCut(1)
,mZCut(0.5f)
,mCandidates()
,mSAonly(kTRUE)
,mCPhi()
,mCDTanL()
,mCDPhi()
,mCZ()
,mCDCAz()
,mCDCAxy()
,mCDN()
,mCDP()
,mCDZ() {
  // This default constructor needs to be provided
  for (int i = 0; i < 4; ++i) {
    mCandidates[i] = new TClonesArray("AliITSUTrackCooked",100000);
  }

#ifdef _TUNING_
  for (int i = 0; i < 6; ++i)
  {
    fGDZ[i] = new TH1F(Form("DGZ%i",i),Form("DZ%i;#Deltaz;Entries",i),500,0.f,5.f);
    fGDXY[i] = new TH1F(Form("DGPhi%i",i),Form("#Delta#Phi%i;#DeltaPhi;Entries",i),500,0.f,TMath::TwoPi());
    fFDZ[i] = new TH1F(Form("DFZ%i",i),Form("DZ%i;#Deltaz;Entries",i),500,0.f,5.f);
    fFDXY[i] = new TH1F(Form("DFPhi%i",i),Form("#Delta#Phi%i;#Delta#Phi;Entries",i),500,0.f,TMath::TwoPi());
  }
  for (int i = 0; i < 5; ++i)
  {
    fGDCAZ[i] = new TH1F(Form("DCAGZ%i",i),Form("DCAZ%i;#Deltaz;Entries",i),500,0.f,5.f);
    fGDCAXY[i] = new TH1F(Form("DCAGXY%i",i),Form("DCAXY%i;#Deltar;Entries",i),500,0.f,5.f);
    fFDCAZ[i] = new TH1F(Form("DCAFZ%i",i),Form("DCAZ%i;#Deltaz;Entries",i),500,0.f,5.f);
    fFDCAXY[i] = new TH1F(Form("DCAFXY%i",i),Form("DCAXY%i;#Deltar;Entries",i),500,0.f,5.f);
  }
  for(int i = 0; i < 4; ++i)
  {
    fGoodCombChi2[i] = new TH1F(Form("goodcombchi2%i",i),Form("%i;#chi^{2};Entries",i),2000,0,0.1);
    fFakeCombChi2[i] = new TH1F(Form("fakecombchi2%i",i),Form("%i;#chi^{2};Entries",i),2000,0,0.1);
    fGoodCombN[i] = new TH1F(Form("goodcombn%i",i),Form("%i;#Deltan;Entries",i),300,0.f,0.03f);
    fFakeCombN[i] = new TH1F(Form("fakecombn%i",i),Form("%i;#Deltan;Entries",i),300,0.f,0.03f);
  }
  fGoodCombChi2[4] = new TH1F("goodcomb4",";#chi^{2};Entries",200,0,500);
  fFakeCombChi2[4] = new TH1F("fakecomb4",";#chi^{2};Entries",200,0,500);
  fTan = new TH1F("tan","tan",2500,0,0.5);
  fTanF = new TH1F("tanF","tanF",2500,0,0.5);
  fPhi = new TH1F("phi","phi",2500,0,TMath::Pi());
  fPhiF = new TH1F("phi","phiF",2500,0,TMath::Pi());
  fNEntries = new TH1F("nentries","nentries",2001,-0.5,2000.5);
#endif
}

Tracker::~Tracker() {
  // d-tor
   for (int i = 0; i < 4; ++i) {
    if (mCandidates[i])
      delete mCandidates[i];
  }
}

bool Tracker::CellParams(int l, ClsInfo_t* c1, ClsInfo_t* c2, ClsInfo_t* c3,
                                  float &curv, float n[3]) {
  // Calculation of cell params and filtering using a DCA cut wrt beam line position.
  // The hit are mapped on a paraboloid space: there a circle is described as plane.
  // The parameter n of the cells is the normal vector to the plane describing the circle in the
  // paraboloid.

  // Mapping the hits
  const float mHit0[3] = {c1->x, c1->y, c1->r * c1->r};
  const float mHit1[3] = {c2->x, c2->y, c2->r * c2->r};
  const float mHit2[3] = {c3->x, c3->y, c3->r * c3->r};
  // Computing the deltas
  const float mD10[3] = {mHit1[0] - mHit0[0],mHit1[1] - mHit0[1],mHit1[2] - mHit0[2]};
  const float mD20[3] = {mHit2[0] - mHit0[0],mHit2[1] - mHit0[1],mHit2[2] - mHit0[2]};
  // External product of the deltas -> n
  n[0] = (mD10[1] * mD20[2]) - (mD10[2] * mD20[1]);
  n[1] = (mD10[2] * mD20[0]) - (mD10[0] * mD20[2]);
  n[2] = (mD10[0] * mD20[1]) - (mD10[1] * mD20[0]);
  // Normalisation
  float norm = TMath::Sqrt((n[0] * n[0]) + (n[1] * n[1]) + (n[2] * n[2]));
  if (norm < 1e-20f || fabs(n[2]) < 1e-20f)
    return false;
  else
    norm = 1.f / norm;
  n[0] *= norm;
  n[1] *= norm;
  n[2] *= norm;
  // Center of the circle
  const float c[2] = {-0.5f * n[0] / n[2], -0.5f * n[1] / n[2]};
  // Constant
  const float k = - n[0] * mHit1[0] - n[1] * mHit1[1] - n[2] * mHit1[2];
  // Radius of the circle
  curv = TMath::Sqrt((1.f - n[2] * n[2] - 4.f * k * n[2]) / (4.f * n[2] * n[2]));
  // Distance of closest approach to the beam line
  const float dca = fabs(curv - sqrt(c[0] * c[0] + c[1] * c[1]));
  // Cut on the DCA
  if (dca > mCDCAxy[l]) {
    return false;
  }
#ifdef _TUNING_
  if (fGood) {
    fGDCAXY[l]->Fill(dca);
  } else {
    fFDCAXY[l]->Fill(dca);
  }
#endif

  curv = 1.f / curv;
  return true;
}

void Tracker::CellsTreeTraversal(vector<CARoad> &roads,
                                          const int &iD, const int &doubl) {

  // Each cells contains a list of neighbours. Each neighbour has presumably other neighbours.
  // This chain of neighbours is, as a matter of fact, a tree and this function goes recursively
  // through it. This function implements a depth first tree browsing.

  // End of the road
  if (doubl < 0) return;

  // [1] add current cell to current cell
  roads.back().AddElement(doubl,iD);
  // We want the right number of elements in the roads also in the case of multiple neighbours
  const int currentN = roads.back().N;

  // [2] loop on the neighbours of the current cell
  for (size_t iN = 0; iN < mCells[doubl][iD].NumberOfNeighbours(); ++iN) {
    const int currD = doubl - 1;
    const int neigh = mCells[doubl][iD](iN);

    // [3] for each neighbour one road
    if (iN > 0) {
      roads.push_back(roads.back());
      roads.back().N = currentN;
    }
    // [4] play this game again until the end of the road
    CellsTreeTraversal(roads,neigh,currD);
  }

  mCells[doubl][iD].SetLevel(0u); // Level = -1
}

Int_t Tracker::Clusters2Tracks(AliESDEvent *event) {
  // This is the main tracking function
  // The clusters must already be loaded

  if (!mSAonly) {
    AliITSUTrackerGlo::Clusters2Tracks(event);
    for (int iL = 0; iL < 7; ++iL) {
      for (int iC = 0; iC < mLayer[iL].GetNClusters(); ++iC)
        if (mLayer[iL].GetClusterUnSorted(iC)->IsClusterUsed())
          mUsedClusters[iL][iC] = true;
    }
  }

  int ntrk = 0, ngood = 0;
  for (int iteration = 0; iteration < mkNumberOfIterations; ++iteration) {

    mCandidates[0]->Clear();
    mCandidates[1]->Clear();
    mCandidates[2]->Clear();
    mCandidates[3]->Clear();

    MakeCells(iteration);
    FindTracksCA(iteration);

    for (int iL = 3; iL >= 0; --iL) {
      const int nCand = mCandidates[iL]->GetEntries();
      int index[nCand];
      float chi2[nCand];

      for (int iC = 0; iC < nCand; ++iC) {
        AliITSUTrackCooked *tr = (AliITSUTrackCooked*)mCandidates[iL]->At(iC);
        chi2[iC] = tr->GetChi2();//(RefitTrack(tr,clInfo,0.,-1));
        index[iC] = iC;
        CookLabel(tr,0.f);
      }

      TMath::Sort(nCand,chi2,index,false);

      for (int iUC = 0; iUC < nCand; ++iUC) {
        const int iC = index[iUC];
        if (chi2[iC] < 0.f) {
          continue;
        }

        AliITSUTrackCooked *tr = (AliITSUTrackCooked*)mCandidates[iL]->At(iC);

        bool sharingCluster = false;
        for (int k = 0; k < tr->GetNumberOmClusters(); ++k) {
          const int layer = (tr->GetClusterIndex(k) & 0xf0000000) >> 28;
          const int idx = (tr->GetClusterIndex(k) & 0x0fffffff);
          if (mUsedClusters[layer][idx]) {
            sharingCluster = true;
            break;
          }
        }

        if (sharingCluster)
          continue;

        for (int k = 0; k < tr->GetNumberOmClusters(); ++k) {
          const int layer = (tr->GetClusterIndex(k) & 0xf0000000) >> 28;
          const int idx = (tr->GetClusterIndex(k) & 0x0fffffff);
          mUsedClusters[layer][idx] = true;
        }

        AliESDtrack outTrack;
        CookLabel(tr,0.f);
        ntrk++;
        if(tr->GetLabel() >= 0) {
          ngood++;
#ifdef _TUNING_
          fGoodCombChi2[4]->Fill(chi2[iC] / (4 + iL));
        } else {
          fFakeCombChi2[4]->Fill(chi2[iC] / (4 + iL));
#endif
        }

        outTrack.UpdateTrackParams(tr,AliESDtrack::kITSin);
        outTrack.SetLabel(tr->GetLabel());
        if (mSAonly) outTrack.SetStatus(AliESDtrack::kITSpureSA);
	event->AddTrack(&outTrack);
      }
    }
  }
  Info("Clusters2Tracks","Reconstructed tracks: %d",ntrk);
  if (ntrk)
    Info("Clusters2Tracks","Good tracks/reconstructed: %f",Float_t(ngood)/ntrk);
  //
  return 0;
}

void Tracker::FindTracksCA(int iteration) {
  // Main pattern recognition routine. It has 4 steps (planning to split in different methods)
  // 1. Tracklet finding (using vertex position)
  // 2. Tracklet association, cells building
  // 3. Handshake between neighbour cells
  // 4. Candidates ("roads") finding and fitting

  // Road finding and fitting. The routine starts from cells with level 5 since they are the head
  // of a full road (candidate with 7 points). Minimum level is 2, so candidates at the end have
  // at least 4 points.
  const int itLevelLimit[3] = {4, 4, 1};
  for (int level = 5; level > itLevelLimit[iteration]; --level) {
    vector<CARoad> roads;
    // Road finding. For each cell at level $(level) a loop on their neighbours to start building
    // the roads.
    for (int iCL = 4; iCL >= level - 1; --iCL) {
      for (size_t iCell = 0; iCell < mCells[iCL].size(); ++iCell) {
        if (mCells[iCL][iCell].GetLevel() != level)
        {
          continue;
        }
        // [1] Add current cell to road
        roads.push_back(CARoad(iCL,iCell));
        // [2] Loop on current cell neighbours
        for(size_t iN = 0; iN < mCells[iCL][iCell].NumberOfNeighbours(); ++iN) {
          const int currD = iCL - 1;
          const int neigh = mCells[iCL][iCell](iN);
          // [3] if more than one neighbour => more than one road, one road for each neighbour
          if(iN > 0)
          {
            roads.push_back(CARoad(iCL,iCell));
          }
          // [4] Essentially the neighbour became the current cell and then go to [1]
          CellsTreeTraversal(roads,neigh,currD);
        }
        mCells[iCL][iCell].SetLevel(0u); // Level = -1
      }
    }

    // Roads fitting
    for (size_t iR = 0; iR < roads.size(); ++iR) {
      if (roads[iR].N != level)
        continue;
      AliITSUTrackCooked tr;
      int first = -1,last = -1;
      ClsInfo_t *cl0 = 0x0,*cl1 = 0x0,*cl2 = 0x0;
      for(int i = 0; i < 5; ++i) {
        if (roads[iR][i] < 0)
          continue;

        if (first < 0) {
          cl0 = mLayer[i][mCells[i][roads[iR][i]].x()];
          tr.SetClusterIndex(i,mLayer[i][mCells[i][roads[iR][i]].x()]->index);
          tr.SetClusterIndex(i + 1,mLayer[i + 1][mCells[i][roads[iR][i]].y()]->index);
          first = i;
        }
        tr.SetClusterIndex(i + 2,mLayer[i + 2][mCells[i][roads[iR][i]].z()]->index);
        last = i;
      }
      AliITSUClusterPix* c = mLayer[last + 2].GetClusterSorted(mCells[last][roads[iR][last]].z());
      cl2 = mLayer[last + 2][mCells[last][roads[iR][last]].z()];
      first = (last + first) / 2;
      cl1 = mLayer[first + 1][mCells[first][roads[iR][first]].y()];
      // Init track parameters
      double cv = Curvature(cl0->x,cl0->y,cl1->x,cl1->y,cl2->x,cl2->y);
      double tgl = TanLambda(cl0->x,cl0->y,cl2->x,cl2->y,cl0->z,cl2->z);

      CATrackingStation::ITSDetInfo_t det = mLayer[last + 2].GetDetInfo(cl2->detid);
      double x = det.xTF + c->GetX(); // I'd like to avoit using AliITSUClusterPix...
      double alp = det.phiTF;
      double par[5] = {c->GetY(),c->GetZ(),0,tgl,cv};
      double cov[15] = {
        5*5,
        0,  5*5,
        0,  0  , 0.7*0.7,
        0,  0,   0,       0.7*0.7,
        0,  0,   0,       0,       10
      };
      tr.Set(x,alp,par,cov);
      AliITSUTrackCooked tt(tr);
      tt.ResetClusters();
      if (RefitAt(2.1, &tt, &tr))
        new((*mCandidates[level - 2])[fCandidates[level - 2]->GetEntriesFast()]) AliITSUTrackCooked(tt);
    }
  }
}

inline AliCluster* Tracker::GetCluster(Int_t index) const {
	const Int_t l=(index & 0xf0000000) >> 28;
	const Int_t c=(index & 0x0fffffff);
	return (AliCluster*)mLayer[l].GetClusterUnSorted(c) ;
}

Double_t Tracker::GetMaterialBudget(const double* pnt0,const double* pnt1, double& x2x0,
	                                           double& rhol) const {
  double par[7];
  if (fUseMatLUT && fMatLUT) {
    double d = fMatLUT->GetMatBudget(pnt0,pnt1,par);
    x2x0 = par[AliITSUMatLUT::kParX2X0];
    rhol = par[AliITSUMatLUT::kParRhoL];
    return d;
  } else {
    MeanMaterialBudget(pnt0,pnt1,par);
    x2x0 = par[1];
    rhol = par[0]*par[4];
    return par[4];
  }
}

Int_t Tracker::LoadClusters(TTree *cluTree) {
  // This function reads the ITSU clusters from the tree,
  // sort them, distribute over the internal tracker arrays, etc

  AliITSUTrackerGlo::LoadClusters(cluTree); // === fITS->LoadClusters(cluTree);
  if (mSAonly) fITS->ProcessClusters();

  // I consider a single vertex event for the moment.
  //TODO: pile-up (trivial here), fast reco of primary vertices (not so trivial)
  double vertex[3] = {GetX(),GetY(),GetZ()};
  for (int iL = 0; iL < 7; ++iL) {
    mLayer[iL].Init(fITS->GetLayerActive(iL),fITS->GetGeom());
    AliVertex v(vertex,1,1);
    mLayer[iL].SortClusters(&v);
    mUsedClusters[iL].resize(mLayer[iL].GetNClusters(),false);
  }
  return 0;
}

void Tracker::MakeCells(int iteration) {
#ifdef _TUNING_
  unsigned int numberOfGoodDoublets = 0, totalNumberOmDoublets = 0;
  unsigned int numberOfGoodCells = 0, totalNumberOmCells = 0;
  unsigned int cellsCombiningSuccesses = 0, cellsWrongCombinations = 0;
#endif

  SetCuts(iteration);
  if (iteration >= 1) {
#ifdef _TUNING_
    ResetHistos();
#endif
    for (int i = 0; i < 5; ++i)
      vector<CACell>().swap(mCells[i]);
    for (int i = 0; i < 6; ++i)
      vector<Doublets>().swap(mDoublets[i]);
  }

  // Trick to speed up the navigation of the doublets array. The lookup table is build like:
  // dLUT[l][i] = n;
  // where n is the index inside mDoublets[l+1] of the first doublets that uses the point
  // mLayer[l+1][i]
  vector<int> dLUT[5];
  for (int iL = 0; iL < 6; ++iL) {
    if (mLayer[iL].GetNClusters() == 0) continue;
    if (iL < 5)
      dLUT[iL].resize(mLayer[iL + 1].GetNClusters(),-1);
    if (dLUT[iL - 1].size() == 0u)
      continue;
    for (int iC = 0; iC < mLayer[iL].GetNClusters(); ++iC) {
      ClsInfo_t* cls = mLayer[iL].GetClusterInfo(iC);
      if (mUsedClusters[iL][cls->index]) {
        continue;
      }
      const float tanL = (cls->z - GetZ()) / cls->r;
      const float extz = tanL * (mkR[iL + 1] - cls->r) + cls->z;
      const int nClust = mLayer[iL + 1].SelectClusters(extz - 2 * mCZ, extz + 2 * fCZ,
                                                       cls->phi - mCPhi, cls->phi + fCPhi);
      bool first = true;

      for (int iC2 = 0; iC2 < nClust; ++iC2) {
        const int iD2 = mLayer[iL + 1].GetNextClusterInfoID();
        ClsInfo_t* cls2 = mLayer[iL + 1].GetClusterInfo(iD2);
        if (mUsedClusters[iL + 1][cls2->index]) {
          continue;
        }
        const float dz = tanL * (cls2->r - cls->r) + cls->z - cls2->z;
        if (TMath::Abs(dz) < mCDZ[iL] && CompareAngles(cls->phi, cls2->phi, mCPhi)) {
          if (first && iL > 0) {
            dLUT[iL - 1][iC] = mDoublets[iL].size();
            first = false;
          }
          const float dTanL = (cls->z - cls2->z) / (cls->r - cls2->r);
          const float phi = TMath::ATan2(cls->y - cls2->y, cls->x - cls2->x);
          mDoublets[iL].push_back(Doublets(iC,iD2,dTanL,phi));
#ifdef _TUNING_
          if (mLayer[iL].GetClusterSorted(iC)->GetLabel(0) ==
              mLayer[iL + 1].GetClusterSorted(iD2)->GetLabel(0) &&
              mLayer[iL].GetClusterSorted(iC)->GetLabel(0) > 0) {
            numberOfGoodDoublets++;
            fGDZ[iL]->Fill(dz);
            fGDXY[iL]->Fill(fabs(cls->phi - cls2->phi));
          } else {
            fFDZ[iL]->Fill(dz);
            fFDXY[iL]->Fill(fabs(cls->phi - cls2->phi));
          }
          totalNumberOmDoublets++;
#endif
        }
      }
      mLayer[iL + 1].ResetFoundIterator();
    }
  }

  // Trick to speed up the navigation of the cells array. The lookup table is build like:
  // tLUT[l][i] = n;
  // where n is the index inside mCells[l+1] of the first cells that uses the doublet
  // mDoublets[l+1][i]
  vector<int> tLUT[4];
  tLUT[0].resize(mDoublets[1].size(),-1);
  tLUT[1].resize(mDoublets[2].size(),-1);
  tLUT[2].resize(mDoublets[3].size(),-1);
  tLUT[3].resize(mDoublets[4].size(),-1);

  for (int iD = 0; iD < 5; ++iD) {
    if (mDoublets[iD + 1].size() == 0u || fDoublets[iD].size() == 0u) continue;

    for (size_t iD0 = 0; iD0 < mDoublets[iD].size(); ++iD0) {
      const int idx = mDoublets[iD][iD0].y;
      bool first = true;
      if (dLUT[iD][idx] == -1) continue;
      for (size_t iD1 = dLUT[iD][idx]; iD1 < mDoublets[iD + 1].size(); ++iD1 {
        if (idx != mDoublets[iD + 1][iD1].x) break;
        if (TMath::Abs(mDoublets[iD][iD0].tanL - fDoublets[iD + 1][iD1].tanL) < mCDTanL &&
            TMath::Abs(mDoublets[iD][iD0].phi - fDoublets[iD + 1][iD1].phi) < mCDPhi) {
          const float tan = 0.5f * (mDoublets[iD][iD0].tanL + fDoublets[iD + 1][iD1].tanL);
          const float extz = -tan * mLayer[iD][mDoublets[iD][iD0].x]->r +
                              mLayer[iD][mDoublets[iD][iD0].x]->z;
          if (fabs(extz - GetZ()) < mCDCAz[iD]) {
#ifdef _TUNING_
            fGood = (mLayer[iD].GetClusterSorted(mDoublets[iD][iD0].x)->GetLabel(0) ==
                     mLayer[iD + 1].GetClusterSorted(mDoublets[iD][iD0].y)->GetLabel(0) &&
                     mLayer[iD].GetClusterSorted(mDoublets[iD][iD0].x)->GetLabel(0) ==
                     mLayer[iD + 2].GetClusterSorted(mDoublets[iD + 1][iD1].y)->GetLabel(0) &&
                     mLayer[iD].GetClusterSorted(mDoublets[iD][iD0].x)->GetLabel(0) > 0);
#endif
            float curv, n[3];
            if (CellParams(iD,mLayer[iD][mDoublets[iD][iD0].x],fLayer[iD + 1][fDoublets[iD][iD0].y],
                           mLayer[iD + 2][mDoublets[iD + 1][iD1].y],curv,n)) {
              if (first && iD > 0) {
                tLUT[iD - 1][iD0] = mCells[iD].size();
                first = false;
              }
              mCells[iD].push_back(CACell(mDoublets[iD][iD0].x,fDoublets[iD][iD0].y,
                                                 mDoublets[iD + 1][iD1].y,iD0,iD1,curv,n));
#ifdef _TUNING_
              if (fGood) {
                fTan->Fill(TMath::Abs(mDoublets[iD][iD0].tanL - fDoublets[iD + 1][iD1].tanL));
                fPhi->Fill(TMath::Abs(mDoublets[iD][iD0].phi - fDoublets[iD + 1][iD1].phi));
                fGDCAZ[iD]->Fill(fabs(extz-GetZ()));
                numberOfGoodCells++;
              } else {
                fTanF->Fill(TMath::Abs(mDoublets[iD][iD0].tanL - fDoublets[iD + 1][iD1].tanL));
                fPhiF->Fill(TMath::Abs(mDoublets[iD][iD0].phi - fDoublets[iD + 1][iD1].phi));
                fFDCAZ[iD]->Fill(fabs(extz - GetZ()));
              }
              totalNumberOmCells++;
#endif
            }
          }
        }
      }
    }
  }

  // Adjacent cells: cells that share 2 points. In the following code adjacent cells are combined.
  // If they meet some requirements (~ same curvature, ~ same n) the innermost cell id is added
  // to the list of neighbours of the outermost cell. When the cell is added to the neighbours of
  // the outermost cell the "level" of the latter is set to the level of the innermost one + 1.
  // ( only if $(level of the innermost) + 1 > $(level of the outermost) )
  for (int iD = 0; iD < 4; ++iD) {
    if (mCells[iD + 1].size() == 0u || tLUT[iD].size() == 0u) continue; // TODO: dealing with holes
    for (size_t c0 = 0; c0 < mCells[iD].size(); ++c0) {
      const int idx = mCells[iD][c0].d1();
      if (tLUT[iD][idx] == -1) continue;
      for (size_t c1 = tLUT[iD][idx]; c1 < mCells[iD + 1].size(); ++c1) {
        if (idx != mCells[iD + 1][c1].d0()) break;
#ifdef _TUNING_
        fGood = (mLayer[iD].GetClusterSorted(mCells[iD][c0].x())->GetLabel(0) ==
                 mLayer[iD + 1].GetClusterSorted(mCells[iD][c0].y())->GetLabel(0) &&
                 mLayer[iD + 1].GetClusterSorted(mCells[iD][c0].y())->GetLabel(0) ==
                 mLayer[iD + 2].GetClusterSorted(mCells[iD][c0].z())->GetLabel(0) &&
                 mLayer[iD + 2].GetClusterSorted(mCells[iD][c0].z())->GetLabel(0) ==
                 mLayer[iD + 3].GetClusterSorted(mCells[iD + 1][c1].z())->GetLabel(0) &&
                 mLayer[iD].GetClusterSorted(mCells[iD][c0].x())->GetLabel(0) > 0);
#endif
        float* n0 = mCells[iD][c0].GetN();
        float* n1 = mCells[iD + 1][c1].GetN();
        const float dn2 = ((n0[0] - n1[0]) * (n0[0] - n1[0]) + (n0[1] - n1[1]) * (n0[1] - n1[1]) +
                           (n0[2] - n1[2]) * (n0[2] - n1[2]));
        const float dp = fabs(mCells[iD][c0].GetCurvature() - fCells[iD + 1][c1].GetCurvature());
        if (dn2 < mCDN[iD] && dp < fCDP[iD]) {
#ifdef _TUNING_
          assert(mCells[iD + 1][c1].Combine(fCells[iD][c0], c0));
          if (fGood) {
            fGoodCombChi2[iD]->Fill(dp);
            fGoodCombN[iD]->Fill(dn2);
            cellsCombiningSuccesses++;
          } else {
            fFakeCombChi2[iD]->Fill(dp);
            fFakeCombN[iD]->Fill(dn2);
            cellsWrongCombinations++;
          }
#else
          mCells[iD + 1][c1].Combine(fCells[iD][c0], c0);
#endif
        }
      }
    }
  }
#ifdef _TUNING_
  Info("MakeCells","Good doublets: %d",numberOfGoodDoublets);
  Info("MakeCells","Number of doublets: %d",totalNumberOmDoublets);
  Info("MakeCells","Good cells: %d",numberOfGoodCells);
  Info("MakeCells","Number of cells: %d",totalNumberOmCells);
  Info("MakeCells","Cells combining successes: %d",cellsCombiningSuccesses);
  Info("MakeCells","Cells wrong combinations: %d",cellsWrongCombinations);
#endif
}

Int_t Tracker::PropagateBack(AliESDEvent * event) {

  Int_t n=event->GetNumberOfTracks();
  Int_t ntrk=0;
  Int_t ngood=0;
  for (Int_t i=0; i<n; i++) {
    AliESDtrack *esdTrack=event->GetTrack(i);

    if ((esdTrack->GetStatus()&AliESDtrack::kITSin)==0) continue;
    if (esdTrack->IsOn(AliESDtrack::kTPCin)) continue; //skip a TPC+ITS track

    AliITSUTrackCooked track(*esdTrack);
    AliITSUTrackCooked temp(*esdTrack);

    temp.ResetCovariance(10.);
    temp.ResetClusters();

    if (RefitAt(40., &temp, &track)) {

      CookLabel(&temp, 0.); //For comparison only
      Int_t label = temp.GetLabel();
      if (label > 0) ngood++;

      esdTrack->UpdateTrackParams(&temp,AliESDtrack::kITSout);
      ntrk++;
    }
  }

  Info("PropagateBack","Back propagated tracks: %d",ntrk);
  if (ntrk)
    Info("PropagateBack","Good tracks/back propagated: %f",Float_t(ngood)/ntrk);

  if (!mSAonly) return AliITSUTrackerGlo::PropagateBack(event);
  return 0;
}

Bool_t Tracker::RefitAt(Double_t xx, AliITSUTrackCooked *t, const AliITSUTrackCooked *c) {
  // This function refits the track "t" at the position "x" using
  // the clusters from "c"

  const int nLayers = 7;
  Int_t index[nLayers];
  Int_t k;
  for (k = 0; k < nLayers; k++) index[k] = -1;
  Int_t nc = c->GetNumberOmClusters();
  for (k = 0; k < nc; k++) {
    Int_t idx = c->GetClusterIndex(k), nl = (idx&0xf0000000)>>28;
    index[nl] = idx;
  }

  Int_t from, to, step;
  if (xx > t->GetX()) {
    from = 0;
    to = nLayers;
    step = +1;
  } else {
    from = nLayers - 1;
    to = -1;
    step = -1;
  }

  for (Int_t i = from; i != to; i += step) {
    Int_t idx = index[i];
    if (idx >= 0) {
      const AliCluster *cl = GetCluster(idx);
      Float_t xr,ar;
      cl->GetXAlphaRefPlane(xr, ar);
      if (!t->Propagate(Double_t(ar), Double_t(xr), GetBz())) {
        return kFALSE;
      }
      Double_t chi2 = t->GetPredictedChi2(cl);
      //      if (chi2 < 100)
      t->Update(cl, chi2, idx);
    } else {
      Double_t r = mkR[i];
      Double_t phi,z;
      if (!t->GetPhiZat(r,phi,z)) {
        return kFALSE;
      }
      if (!t->Propagate(phi, r, GetBz())) {
        return kFALSE;
      }
    }
    Double_t xx0 = (i > 2) ? 0.008 : 0.003;  // Rough layer thickness
    Double_t x0  = 9.36; // Radiation length of Si [cm]
    Double_t rho = 2.33; // Density of Si [g/cm^3]
    Double_t mass = t->GetMass();
    t->CorrectForMeanMaterial(xx0, - step * xx0 * x0 * rho, mass, kTRUE);
  }

  if (!t->PropagateTo(xx,0.,0.)) return kFALSE;
  return kTRUE;
}

Int_t Tracker::RefitInward(AliESDEvent * event) {
  // Some final refit, after the outliers get removed by the smoother ?
  // The clusters must be loaded

  Int_t n = event->GetNumberOfTracks();
  Int_t ntrk = 0;
  Int_t ngood = 0;
  for (Int_t i = 0; i < n; i++) {
    AliESDtrack *esdTrack = event->GetTrack(i);

    if ((esdTrack->GetStatus() & AliESDtrack::kITSout) == 0) continue;
    if (esdTrack->IsOn(AliESDtrack::kTPCin)) continue; //skip a TPC+ITS track
    AliITSUTrackCooked track(*esdTrack);
    AliITSUTrackCooked temp(*esdTrack);

    temp.ResetCovariance(10.);
    temp.ResetClusters();

    if (!RefitAt(2.1, &temp, &track)) continue;
    //Cross the beam pipe
    if (!temp.PropagateTo(1.8, 2.27e-3, 35.28 * 1.848)) continue;

    CookLabel(&temp, 0.); //For comparison only
    Int_t label = temp.GetLabel();
    if (label > 0) ngood++;

    esdTrack->UpdateTrackParams(&temp,AliESDtrack::kITSrefit);
    ntrk++;
  }

  Info("RefitInward","Refitted tracks: %d",ntrk);
  if (ntrk)
    Info("RefitInward","Good tracks/refitted: %f",Float_t(ngood)/ntrk);

  if (!mSAonly) return AliITSUTrackerGlo::RefitInward(event);
  return 0;
}

void Tracker::SetCuts(int it) {
  switch (it) {
    case 0:
      mCPhi = mPhiCut;
      mCDTanL = kDoublTanL;
      mCDPhi = kDoublPhi;
      mCZ = mZCut;
      for (int i = 0; i < 5; ++i) {
        mCDCAxy[i] = kmaxDCAxy[i];
        mCDCAz[i] = kmaxDCAz[i];
      }
      for (int i = 0; i < 4; ++i) {
        mCDN[i] = kmaxDN[i];
        mCDP[i] = kmaxDP[i];
      }
      for (int i = 0; i < 6; ++i) {
        mCDZ[i] = kmaxDZ[i];
      }
      break;

    default:
      mCPhi = 3.f * mPhiCut;
      mCDTanL = kDoublTanL1;
      mCDPhi = kDoublPhi1;
      mCZ = mZCut;
      for (int i = 0; i < 5; ++i) {
        mCDCAxy[i] = kmaxDCAxy1[i];
        mCDCAz[i] = kmaxDCAz1[i];
      }
      for (int i = 0; i < 4; ++i) {
        mCDN[i] = kmaxDN1[i];
        mCDP[i] = kmaxDP1[i];
      }
      for (int i = 0; i < 6; ++i) {
        mCDZ[i] = kmaxDZ1[i];
      }

      break;
  }
}

void Tracker::UnloadClusters() {
  // This function unloads ITSU clusters from the memory
  for (int i = 0;i < 7;++i) {
    mLayer[i].Clear();
    mUsedClusters[i].clear();
  }
  for (int i = 0; i < 6; ++i) {
    mDoublets[i].clear();
  }
  for (int i = 0; i < 5; ++i) {
    mCells[i].clear();
  }
  for (int i = 0; i < 4; ++i) {
    mCandidates[i]->Clear("C");
  }
  AliITSUTrackerGlo::UnloadClusters();
}

#ifdef _TUNING_

void Tracker::ResetHistos() {
  for (int i = 0; i < 6; ++i) {
    fGDZ[i]->Reset();
    fGDXY[i]->Reset();
    fFDZ[i]->Reset();
    fFDXY[i]->Reset();
  }
  for (int i = 0; i < 5; ++i) {
    fGoodCombChi2[i]->Reset();
    fFakeCombChi2[i]->Reset();
    fGDCAZ[i]->Reset();
    fGDCAXY[i]->Reset();
    fFDCAZ[i]->Reset();
    fFDCAXY[i]->Reset();
  }
  for (int i = 0; i < 4; ++i) {
    fGoodCombN[i]->Reset();
    fFakeCombN[i]->Reset();
  }
  fTan->Reset();
  fTanF->Reset();
  fPhi->Reset();
  fPhiF->Reset();
  fNEntries->Reset();
}
#endif

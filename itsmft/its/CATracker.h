//-*- Mode: C++ -*-
// ************************************************************************
// This file is property of and copyright by the ALICE ITSU Project       *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

#ifndef ALIITSUCATRACKER_H
#define ALIITSUCATRACKER_H

#include <vector>

#include "CAaux.h"
#include "CATrackingStation.h"

using std::vector;
class TrackPC;

namespace AliceO2 {
  namespace ITS {
    namespace CA {
      class Tracker {
        public:
          Tracker();
          // These functions must be implemented
          int Clusters2Tracks();
          int PropagateBack();
          int RefitInward();
          int LoadClusters(void *ct);
          void UnloadClusters();
          // Possibly, other public functions
          float GetMaterialBudget(const double* p0, const double* p1, double& x2x0, double& rhol) const;
          bool     GetSAonly() const { return mSAonly; }
          void     SetChi2Cut(float cut) { mChi2Cut = cut; }
          void     SetPhiCut(float cut) { mPhiCut = cut; }
          void     SetSAonly(bool sa = kTRUE) { mSAonly = sa; }
          void     SetZCut(float cut) { mZCut = cut; }
          //
        private:
          Tracker(const Tracker&);
          Tracker &operator=(const Tracker &tr);
          //
          bool   CellParams(int l, Cluster* c1, Cluster* c2, Cluster* c3, float &curv, float np[3]);
          void   CellsTreeTraversal(vector<Road> &roads, const int &iD, const int &doubl);
          void   FindTracksCA(int iteration);
          void   MakeCells(int iteration);
          bool   RefitAt(float xx, TrackPC* t, const TrackPC* c);
          void   SetCuts(int it);
          void   SetLabel(TrackPC &t, float wrong);
          //
          TrackingStation       mLayer[7];
          vector<bool>          mUsedClusters[7];
          float                 mChi2Cut;
          float                 mPhiCut;
          float                 mZCut;
          vector<Doublets>      mDoublets[6];
          vector<Cell>          mCells[5];
          TClonesArray         *mCandidates[4];
          bool                  mSAonly;             // kTRUE if the standalone tracking only
          // Cuts
          float mCPhi;
          float mCDTanL;
          float mCDPhi;
          float mCZ;
          float mCDCAz[5];
          float mCDCAxy[5];
          float mCDN[4];
          float mCDP[4];
          float mCDZ[6];
          //
          static const float              mkChi2Cut;      // chi2 cut during track merging
          static const int                mkNumberOfIterations;
          static const float              mkR[7];
          //
      };
    } // namespace CA
  } // namespace ITS
} // namespace AliceO2

#endif // ALIITSUCATRACKER_H

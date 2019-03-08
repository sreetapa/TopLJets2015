#ifndef ForestPFCands_h
#define ForestPFCands_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "vector"

class ForestPFCands {
public :
  ForestPFCands(TChain *t)
    {
      t->SetBranchAddress("pfId", &pfId);
      t->SetBranchAddress("pfPt", &pfPt);
      t->SetBranchAddress("pfEta", &pfEta);
      t->SetBranchAddress("pfPhi", &pfPhi);
      t->SetBranchAddress("pfM", &pfM);
      t->SetBranchAddress("trkAlgo", &trkAlgo);
      t->SetBranchAddress("trkPtError", &trkPtError);
      t->SetBranchAddress("trkNHit", &trkNHit);
      t->SetBranchAddress("trkChi2", &trkChi2);
      t->SetBranchAddress("trkNdof", &trkNdof);
    }
  ~ForestPFCands() {}

   vector<int>     *pfId;
   vector<float>   *pfPt;
   vector<float>   *pfEta;
   vector<float>   *pfPhi;
   vector<float>   *pfM;
   vector<int>     *trkAlgo;
   vector<float>   *trkPtError;
   vector<float>   *trkNHit;
   vector<float>   *trkChi2;
   vector<float>   *trkNdof;
};



#endif 

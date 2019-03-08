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
      t->SetBranchAddress("pfId", &pfId, &b_pfId);
      t->SetBranchAddress("pfPt", &pfPt, &b_pfPt);
      t->SetBranchAddress("pfEta", &pfEta, &b_pfEta);
      t->SetBranchAddress("pfPhi", &pfPhi, &b_pfPhi);
      t->SetBranchAddress("pfM", &pfM, &b_pfM);
      t->SetBranchAddress("trkAlgo", &trkAlgo, &b_trkAlgo);
      t->SetBranchAddress("trkPtError", &trkPtError, &b_trkPtError);
      t->SetBranchAddress("trkNHit", &trkNHit, &b_trkNHit);
      t->SetBranchAddress("trkChi2", &trkChi2, &b_trkChi2);
      t->SetBranchAddress("trkNdof", &trkNdof, &b_trkNdof);     
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

#ifndef ForestPFCands_h
#define ForestPFCands_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>

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

   std::vector<int>     *pfId;
   std::vector<float>   *pfPt;
   std::vector<float>   *pfEta;
   std::vector<float>   *pfPhi;
   std::vector<float>   *pfM;
   std::vector<int>     *trkAlgo;
   std::vector<float>   *trkPtError;
   std::vector<float>   *trkNHit;
   std::vector<float>   *trkChi2;
   std::vector<float>   *trkNdof;
};



#endif 

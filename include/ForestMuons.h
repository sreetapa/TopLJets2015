#ifndef ForestMuons_h
#define ForestMuons_h

#define maxForestMuons 10

#include <iostream>
#include <vector>
#include "TBranch.h"

class ForestMuons {
public :
   ForestMuons(TChain *t)
     {
       t->SetBranchStatus("mu*", 1);
       t->SetBranchAddress("muPt", &muPt);
       t->SetBranchAddress("muPhi", &muPhi);
       t->SetBranchAddress("muEta", &muEta);
       t->SetBranchAddress("muCharge", &muCharge);
       t->SetBranchAddress("muType", &muType);
       t->SetBranchAddress("muIsGood", &muIsGood);
       t->SetBranchAddress("muD0", &muD0);
       t->SetBranchAddress("muDz", &muDz);
       t->SetBranchAddress("muChi2NDF", &muChi2NDF);
       t->SetBranchAddress("muInnerD0", &muInnerD0);
       t->SetBranchAddress("muInnerDz", &muInnerDz);
       t->SetBranchAddress("muMuonHits", &muMuonHits);
       t->SetBranchAddress("muStations", &muStations);
       t->SetBranchAddress("muTrkLayers", &muTrkLayers);
       t->SetBranchAddress("muPixelHits", &muPixelHits);    
       t->SetBranchAddress("muMuonHits", &muMuonHits);
       t->SetBranchAddress("muTrkQuality", &muTrkQuality);
       t->SetBranchAddress("muStations", &muStations);
       t->SetBranchAddress("muIsoTrk", &muIsoTrk);
       t->SetBranchAddress("muPFChIso", &muPFChIso);
       t->SetBranchAddress("muPFPhoIso", &muPFPhoIso);
       t->SetBranchAddress("muPFNeuIso", &muPFNeuIso);
       t->SetBranchAddress("muPFPUIso", &muPFPUIso);
     };
   ~ForestMuons(){};

   // RecoMuon info
   std::vector<float>   *muPt = 0;
   std::vector<float>   *muEta = 0;
   std::vector<float>   *muPhi = 0;
   std::vector<int>     *muCharge = 0;
   std::vector<int>     *muType = 0;
   std::vector<int>     *muIsGood = 0;
   std::vector<float>   *muD0 = 0;
   std::vector<float>   *muDz = 0;
   std::vector<float>   *muChi2NDF = 0;
   std::vector<float>   *muInnerD0 = 0;
   std::vector<float>   *muInnerDz = 0;
   std::vector<int>     *muTrkLayers = 0;
   std::vector<int>     *muPixelLayers = 0;
   std::vector<int>     *muPixelHits = 0;
   std::vector<int>     *muMuonHits = 0;
   std::vector<int>     *muTrkQuality = 0;
   std::vector<int>     *muStations = 0;
   std::vector<float>   *muIsoTrk = 0;
   std::vector<float>   *muPFChIso = 0;
   std::vector<float>   *muPFPhoIso = 0;
   std::vector<float>   *muPFNeuIso = 0;
   std::vector<float>   *muPFPUIso = 0;
 };
#endif

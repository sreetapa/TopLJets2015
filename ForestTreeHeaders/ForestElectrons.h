#ifndef ForestElectrons_h
#define ForestElectrons_h

#include <iostream>
#include <vector>
#include "TBranch.h"

class ForestElectrons {
public :
   ForestElectrons(TChain *t)
     {
       t->SetBranchStatus("ele*", 1);    
       t->SetBranchAddress("elePt", &elePt);
       t->SetBranchAddress("elePhi", &elePhi);
       t->SetBranchAddress("eleEta", &eleEta);
       t->SetBranchAddress("eleCharge", &eleCharge);
       t->SetBranchAddress("eleSigmaIEtaIEta", &eleSigmaIEtaIEta);
       t->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx);
       t->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx);
       t->SetBranchAddress("eleHoverE", &eleHoverE);
       t->SetBranchAddress("eleD0", &eleD0);
       t->SetBranchAddress("eleDz", &eleDz);
       t->SetBranchAddress("eleEoverPInv", &eleEoverPInv);
       t->SetBranchAddress("eleSCPhi", &eleSCPhi);
       t->SetBranchAddress("eleSCEta", &eleSCEta);
       t->SetBranchAddress("eleSigmaIEtaIEta_2012", &eleSigmaIEtaIEta_2012);
       t->SetBranchAddress("eleSigmaIPhiIPhi", &eleSigmaIPhiIPhi);
       t->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth);
       t->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth);
       t->SetBranchAddress("eleMissHits", &eleMissHits);
       t->SetBranchAddress("elePFChIso03", &elePFChIso03);
       t->SetBranchAddress("elePFPhoIso03", &elePFPhoIso03);
       t->SetBranchAddress("elePFNeuIso03", &elePFNeuIso03);
     };
   ~ForestElectrons(){};

   // RecoElectron info
   Int_t           nEle;
   std::vector<Int_t>     *eleCharge=0;
   std::vector<Int_t>     *eleChargeConsistent=0;
   std::vector<Float_t>   *eleEn=0;
   std::vector<Float_t>   *eleD0=0;
   std::vector<Float_t>   *eleDz=0;
   std::vector<Float_t>   *eleD0Err=0;
   std::vector<Float_t>   *eleDzErr=0;
   std::vector<Float_t>   *eleTrkPt=0;
   std::vector<Float_t>   *eleTrkEta=0;
   std::vector<Float_t>   *eleTrkPhi=0;
   std::vector<Int_t>     *eleTrkCharge=0;
   std::vector<Float_t>   *eleTrkChi2=0;
   std::vector<Float_t>   *eleTrkNdof=0;
   std::vector<Float_t>   *eleTrkNormalizedChi2=0;
   std::vector<Int_t>     *eleTrkValidHits=0;
   std::vector<Int_t>     *eleTrkLayers=0;
   std::vector<Float_t>   *elePt=0;
   std::vector<Float_t>   *eleEta=0;
   std::vector<Float_t>   *elePhi=0;
   std::vector<Float_t>   *eleSCEn=0;
   std::vector<Float_t>   *eleESEn=0;
   std::vector<Float_t>   *eleSCEta=0;
   std::vector<Float_t>   *eleSCPhi=0;
   std::vector<Float_t>   *eleSCRawEn=0;
   std::vector<Float_t>   *eleSCEtaWidth=0;
   std::vector<Float_t>   *eleSCPhiWidth=0;
   std::vector<Float_t>   *eleHoverE=0;
   std::vector<Float_t>   *eleEoverP=0;
   std::vector<Float_t>   *eleEoverPInv=0;
   std::vector<Float_t>   *eleBrem=0;
   std::vector<Float_t>   *eledEtaAtVtx=0;
   std::vector<Float_t>   *eledPhiAtVtx=0;
   std::vector<Float_t>   *eleSigmaIEtaIEta=0;
   std::vector<Float_t>   *eleSigmaIEtaIEta_2012=0;
   std::vector<Float_t>   *eleSigmaIPhiIPhi=0;
   std::vector<Int_t>     *eleMissHits=0;
   std::vector<Float_t>   *eleESEffSigmaRR=0;
   std::vector<Float_t>   *elePFChIso=0;
   std::vector<Float_t>   *elePFPhoIso=0;
   std::vector<Float_t>   *elePFNeuIso=0;
   std::vector<Float_t>   *elePFPUIso=0;
   std::vector<Float_t>   *elePFChIso03=0;
   std::vector<Float_t>   *elePFPhoIso03=0;
   std::vector<Float_t>   *elePFNeuIso03=0;
   std::vector<Float_t>   *elePFChIso04=0;
   std::vector<Float_t>   *elePFPhoIso04=0;
   std::vector<Float_t>   *elePFNeuIso04=0;
   std::vector<Float_t>   *eleBC1E=0;
   std::vector<Float_t>   *eleBC1Eta=0;
   std::vector<Float_t>   *eleBC2E=0;
   std::vector<Float_t>   *eleBC2Eta=0;
};
#endif

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
       t->SetBranchAddress("elePt", &fForestEle.elePt);
       t->SetBranchAddress("elePhi", &fForestEle.elePhi);
       t->SetBranchAddress("eleEta", &fForestEle.eleEta);
       t->SetBranchAddress("eleCharge", &fForestEle.eleCharge);
       t->SetBranchAddress("eleSigmaIEtaIEta", &fForestEle.eleSigmaIEtaIEta);
       t->SetBranchAddress("eledEtaAtVtx", &fForestEle.eledEtaAtVtx);
       t->SetBranchAddress("eledPhiAtVtx", &fForestEle.eledPhiAtVtx);
       t->SetBranchAddress("eleHoverE", &fForestEle.eleHoverE);
       t->SetBranchAddress("eleD0", &fForestEle.eleD0);
       t->SetBranchAddress("eleDz", &fForestEle.eleDz);
       t->SetBranchAddress("eleEoverPInv", &fForestEle.eleEoverPInv);
       t->SetBranchAddress("eleSCPhi", &fForestEle.eleSCPhi);
       t->SetBranchAddress("eleSCEta", &fForestEle.eleSCEta);
       t->SetBranchAddress("eleSigmaIEtaIEta_2012", &fForestEle.eleSigmaIEtaIEta_2012);
       t->SetBranchAddress("eleSigmaIPhiIPhi", &fForestEle.eleSigmaIPhiIPhi);
       t->SetBranchAddress("eleSCEtaWidth", &fForestEle.eleSCEtaWidth);
       t->SetBranchAddress("eleSCPhiWidth", &fForestEle.eleSCPhiWidth);
       t->SetBranchAddress("eleMissHits", &fForestEle.eleMissHits);
       t->SetBranchAddress("elePFChIso03", &fForestEle.elePFChIso03);
       t->SetBranchAddress("elePFPhoIso03", &fForestEle.elePFPhoIso03);
       t->SetBranchAddress("elePFNeuIso03", &fForestEle.elePFNeuIso03);
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

#ifndef HiTree_h
#define HiTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


class HiTree {
public :
 HiTree(TChain *t) : weight(1.0),  ttbar_w(0)
    {
      t->SetBranchAddress("run", &run);
      t->SetBranchAddress("evt", &evt);
      t->SetBranchAddress("lumi", &lumi);
      t->SetBranchAddress("vx", &vx);
      t->SetBranchAddress("vy", &vy);
      t->SetBranchAddress("vz", &vz);
      t->SetBranchAddress("hiBin", &hiBin);
      t->SetBranchAddress("weight", &weight);
      t->SetBranchAddress("ttbar_w", &ttbar_w);
      t->SetBranchAddress("hiHF", &hiHF);
      t->SetBranchAddress("hiHFplus", &hiHFplus);
      t->SetBranchAddress("hiHFminus", &hiHFminus);
      t->SetBranchAddress("hiHFECut", &hiHFECut);
      t->SetBranchAddress("hiHFECutPlus", &hiHFECutPlus);
      t->SetBranchAddress("hiHFECutMinus", &hiHFECutMinus);
      t->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4);
      t->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4);
      t->SetBranchAddress("hiZDC", &hiZDC);
      t->SetBranchAddress("hiZDCplus", &hiZDCplus);
      t->SetBranchAddress("hiZDCminus", &hiZDCminus);
      t->SetBranchAddress("hiHFhit", &hiHFhit);
      t->SetBranchAddress("hiHFhitPlus", &hiHFhitPlus);
      t->SetBranchAddress("hiHFhitMinus", &hiHFhitMinus);
      t->SetBranchAddress("hiET", &hiET);
      t->SetBranchAddress("hiEE", &hiEE);
      t->SetBranchAddress("hiEB", &hiEB);
      t->SetBranchAddress("hiEEplus", &hiEEplus);
      t->SetBranchAddress("hiEEminus", &hiEEminus);
      t->SetBranchAddress("hiNpix", &hiNpix);
      t->SetBranchAddress("hiNpixPlus", &hiNpixPlus);
      t->SetBranchAddress("hiNpixMinus", &hiNpixMinus);
      t->SetBranchAddress("hiNpixelTracks", &hiNpixelTracks);
      t->SetBranchAddress("hiNpixelTracksPlus", &hiNpixelTracksPlus);
      t->SetBranchAddress("hiNpixelTracksMinus", &hiNpixelTracksMinus);
      t->SetBranchAddress("hiNtracks", &hiNtracks);
      t->SetBranchAddress("hiNtracksPtCut", &hiNtracksPtCut);
      t->SetBranchAddress("hiNtracksEtaCut", &hiNtracksEtaCut);
      t->SetBranchAddress("hiNtracksEtaPtCut", &hiNtracksEtaPtCut);
      t->SetBranchAddress("hiNevtPlane", &hiNevtPlane);
      t->SetBranchAddress("hiEvtPlanes", &hiEvtPlanes);
    }
  ~HiTree(){}
   // Declaration of leaf types
   UInt_t          run;
   ULong64_t       evt;
   UInt_t          lumi;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Int_t           hiBin;
   Float_t         weight;
   Float_t         hiHF;
   Float_t         hiHFplus;
   Float_t         hiHFminus;
   Float_t         hiHFECut;
   Float_t         hiHFECutPlus;
   Float_t         hiHFECutMinus;
   Float_t         hiHFplusEta4;
   Float_t         hiHFminusEta4;
   Float_t         hiZDC;
   Float_t         hiZDCplus;
   Float_t         hiZDCminus;
   Float_t         hiHFhit;
   Float_t         hiHFhitPlus;
   Float_t         hiHFhitMinus;
   Float_t         hiET;
   Float_t         hiEE;
   Float_t         hiEB;
   Float_t         hiEEplus;
   Float_t         hiEEminus;
   Int_t           hiNpix;
   Int_t           hiNpixPlus;
   Int_t           hiNpixMinus;
   Int_t           hiNpixelTracks;
   Int_t           hiNpixelTracksPlus;
   Int_t           hiNpixelTracksMinus;
   Int_t           hiNtracks;
   Int_t           hiNtracksPtCut;
   Int_t           hiNtracksEtaCut;
   Int_t           hiNtracksEtaPtCut;
   Int_t           hiNevtPlane;
   Float_t         hiEvtPlanes[1];   //[hiNevtPlane]
   std::vector<float> *ttbar_w;
};

#endif

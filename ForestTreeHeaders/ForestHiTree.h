#ifndef HiTree_h
#define HiTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


class HiTree {
public :
  HiTree(TChain *t) 
    {
      t->SetBranchAddress("run", &run, &b_run);
      t->SetBranchAddress("evt", &evt, &b_evt);
      t->SetBranchAddress("lumi", &lumi, &b_lumi);
      t->SetBranchAddress("vx", &vx, &b_vx);
      t->SetBranchAddress("vy", &vy, &b_vy);
      t->SetBranchAddress("vz", &vz, &b_vz);
      t->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
      t->SetBranchAddress("hiHF", &hiHF, &b_hiHF);
      t->SetBranchAddress("hiHFplus", &hiHFplus, &b_hiHFplus);
      t->SetBranchAddress("hiHFminus", &hiHFminus, &b_hiHFminus);
      t->SetBranchAddress("hiHFECut", &hiHFECut, &b_hiHFECut);
      t->SetBranchAddress("hiHFECutPlus", &hiHFECutPlus, &b_hiHFECutPlus);
      t->SetBranchAddress("hiHFECutMinus", &hiHFECutMinus, &b_hiHFECutMinus);
      t->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4, &b_hiHFplusEta4);
      t->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4, &b_hiHFminusEta4);
      t->SetBranchAddress("hiZDC", &hiZDC, &b_hiZDC);
      t->SetBranchAddress("hiZDCplus", &hiZDCplus, &b_hiZDCplus);
      t->SetBranchAddress("hiZDCminus", &hiZDCminus, &b_hiZDCminus);
      t->SetBranchAddress("hiHFhit", &hiHFhit, &b_hiHFhit);
      t->SetBranchAddress("hiHFhitPlus", &hiHFhitPlus, &b_hiHFhitPlus);
      t->SetBranchAddress("hiHFhitMinus", &hiHFhitMinus, &b_hiHFhitMinus);
      t->SetBranchAddress("hiET", &hiET, &b_hiET);
      t->SetBranchAddress("hiEE", &hiEE, &b_hiEE);
      t->SetBranchAddress("hiEB", &hiEB, &b_hiEB);
      t->SetBranchAddress("hiEEplus", &hiEEplus, &b_hiEEplus);
      t->SetBranchAddress("hiEEminus", &hiEEminus, &b_hiEEminus);
      t->SetBranchAddress("hiNpix", &hiNpix, &b_hiNpix);
      t->SetBranchAddress("hiNpixPlus", &hiNpixPlus, &b_hiNpixPlus);
      t->SetBranchAddress("hiNpixMinus", &hiNpixMinus, &b_hiNpixMinus);
      t->SetBranchAddress("hiNpixelTracks", &hiNpixelTracks, &b_hiNpixelTracks);
      t->SetBranchAddress("hiNpixelTracksPlus", &hiNpixelTracksPlus, &b_hiNpixelTracksPlus);
      t->SetBranchAddress("hiNpixelTracksMinus", &hiNpixelTracksMinus, &b_hiNpixelTracksMinus);
      t->SetBranchAddress("hiNtracks", &hiNtracks, &b_hiNtracks);
      t->SetBranchAddress("hiNtracksPtCut", &hiNtracksPtCut, &b_hiNtracksPtCut);
      t->SetBranchAddress("hiNtracksEtaCut", &hiNtracksEtaCut, &b_hiNtracksEtaCut);
      t->SetBranchAddress("hiNtracksEtaPtCut", &hiNtracksEtaPtCut, &b_hiNtracksEtaPtCut);
      t->SetBranchAddress("hiNevtPlane", &hiNevtPlane, &b_hiNevtPlane);
      t->SetBranchAddress("hiEvtPlanes", &hiEvtPlanes, &b_hiEvtPlanes);      
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
};

#endif

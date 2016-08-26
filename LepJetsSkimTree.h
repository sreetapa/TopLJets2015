#ifndef LepJetsSkimTree_h
#define LepJetsSkimTree_h
#include "TTree.h"
#include "TH2F.h"
#include <iostream>

TTree* skimTree_p = 0;
const Int_t nLep = 2; // perguntar
const Int_t eleID = 11;
const Int_t muID = 13;
const Float_t eleM = .000510998910;
const Float_t muM = .1056583715;
UInt_t run_, lumi_;
ULong64_t evt_;
Int_t hiBin_;
Float_t vz_;
Float_t weight_;
Int_t nLep_;
Int_t lepID_[nLep];
Float_t lepPt_[nLep];
Float_t lepPhi_[nLep];
Float_t lepEta_[nLep];
Int_t lepChg_[nLep];
Float_t lepIso_[nLep];
Float_t lepInnerDz_[nLep];
Int_t EventCount_;
const Int_t nMaxJets = 500;

Int_t nJt_;
Float_t jtPt_[nMaxJets]; //*
Int_t Index_;

Float_t IsoPt_; 
Float_t NIsoPt_;
Float_t csvV2IsoPt_; 
Float_t csvV2NIsoPt_;  

Float_t jtPhi_[nMaxJets];
Float_t jtEta_[nMaxJets];
Float_t jtM_[nMaxJets];
Float_t discr_csvV1_[nMaxJets];
Float_t discr_csvV2_[nMaxJets];  //*
Float_t discr_tcHighEff_[nMaxJets];
Float_t discr_tcHighPur_[nMaxJets];
Int_t refparton_flavorForB_[nMaxJets];
const Int_t nMaxGen = 5000;
Int_t nGen_;
Int_t genPdg_[nMaxGen];
Float_t genPt_[nMaxGen];
Float_t genPhi_[nMaxGen];
Float_t genEta_[nMaxGen];
Float_t genChg_[nMaxGen];

Float_t deltaPhi_[nMaxJets]; //*
Float_t deltaEta_[nMaxJets]; //*
Float_t ptMuJet_[nMaxJets]; //*

Float_t svtxM_[nMaxJets]; //*
Float_t svtxPt_[nMaxJets]; //*
Float_t discr_Prob_[nMaxJets]; //*
Float_t discr_ProbforB_[nMaxJets];
Float_t discr_ProbforC_[nMaxJets];
Float_t discr_ProbforBack_[nMaxJets];

Float_t ptJP_[nMaxJets]; //*
Float_t R_[nMaxJets]; //*
Float_t deltaphiJJ_[nMaxJets]; //*
//Float_t CSVJP_[nMaxJets]; //*
Float_t svtxptJP_[nMaxJets]; //*
Float_t discrProbJP_[nMaxJets]; //*
Float_t discrcsvV2JP_[nMaxJets]; //*
Float_t Pthat_;

// List of branches
TBranch        *b_run;   //!
TBranch        *b_evt;   //!
TBranch        *b_lumi;   //!
TBranch        *b_hiBin;   //!
TBranch        *b_vz;   //!
TBranch        *b_weight; //!
TBranch        *b_nLep;   //!
TBranch        *b_lepID;   //!
TBranch        *b_lepPt;   //!
TBranch        *b_lepPhi;   //!
TBranch        *b_lepEta;   //!
TBranch        *b_lepChg;   //!
TBranch        *b_lepInnerDz; //!
TBranch        *b_nJt;   //!
TBranch        *b_jtPt;   //! IsoPt
TBranch        *b_EventCount;   //!
TBranch        *b_IsoPt;
TBranch        *b_NIsoPt; 
TBranch        *b_csvV2IsoPt;
TBranch        *b_Index;
TBranch        *b_Pthat;


TBranch        *b_jtPhi;   //!
TBranch        *b_jtEta;   //!
TBranch        *b_jtM;   //!
TBranch        *b_discr_csvV1;//!
TBranch        *b_discr_csvV2;//!
TBranch        *b_discr_tcHighEff;//!
TBranch        *b_discr_tcHighPur;//!
TBranch        *b_refparton_flavorForB;//!
TBranch        *b_nGen;//!
TBranch        *b_genPdg;//!
TBranch        *b_genPt;//!
TBranch        *b_genPhi;//!
TBranch        *b_genEta;//!
TBranch        *b_genChg;//!
TBranch        *b_deltaPhi;   //!
TBranch        *b_deltaEta;  
TBranch        *b_ptmujet;   //!

TBranch        *b_svtxM;   //!
TBranch        *b_svtxPt;  
TBranch        *b_discr_Prob;   //!
TBranch        *b_discr_ProbforB;
TBranch        *b_discr_ProbforC;
TBranch        *b_discr_ProbforBack;


TBranch        *b_ptJP;   //!
TBranch        *b_R;  
TBranch        *b_deltaphiJJ;   //!
//TBranch        *b_CSVJP;  
TBranch        *b_lepIso;   //!
TBranch        *b_svtxptJP;   //!
TBranch        *b_discrProbJP;  
TBranch        *b_discrcsvV2JP;   //!

void BookTree()
{
  if(skimTree_p == NULL){
    std::cout << "BOOKTREE error; skimTree_p is NULL. return" << std::endl;
    return;
  }
  skimTree_p->Branch("run", &run_, "run/i");
  skimTree_p->Branch("evt", &evt_, "evt/l");
  skimTree_p->Branch("lumi", &lumi_, "lumi/i");
  skimTree_p->Branch("hiBin", &hiBin_, "hiBin/I");
  skimTree_p->Branch("vz", &vz_, "vz/F");
  skimTree_p->Branch("weight", &weight_, "weight/F");
  skimTree_p->Branch("nLep", &nLep_, "nLep/I");
  skimTree_p->Branch("lepID", lepID_, Form("lepID[%d]/I",nLep));
  skimTree_p->Branch("lepPt", lepPt_, Form("lepPt[%d]/F", nLep));
  skimTree_p->Branch("lepPhi", lepPhi_, Form("lepPhi[%d]/F", nLep));
  skimTree_p->Branch("lepEta", lepEta_, Form("lepEta[%d]/F", nLep));
  skimTree_p->Branch("lepChg", lepChg_, Form("lepChg[%d]/I", nLep));
  skimTree_p->Branch("pthat", &Pthat_, "pthat/F");
  skimTree_p->Branch("lepIso", lepIso_, Form("lepIso[%d]/F", nLep));

  skimTree_p->Branch("lepInnerDz", lepInnerDz_, Form("lepInnerDz[%d]/F", nLep)); 
  skimTree_p->Branch("EventCount", &EventCount_, "EventCount/I");
  skimTree_p->Branch("Index", &Index_, "Index/I");
  skimTree_p->Branch("nJt", &nJt_, "nJt/I");
  skimTree_p->Branch("jtPt", jtPt_, "jtPt[nJt]/F");

  skimTree_p->Branch("IsoPt", &IsoPt_, "IsoPt/F");
  skimTree_p->Branch("NIsoPt", &NIsoPt_, "NIsoPt/F");
  skimTree_p->Branch("csvV2IsoPt", &csvV2IsoPt_, "csvV2IsoPt/F");
  skimTree_p->Branch("csvV2NIsoPt", &csvV2NIsoPt_, "csvV2NIsoPt/F");

  skimTree_p->Branch("jtPhi", jtPhi_, "jtPhi[nJt]/F");
  skimTree_p->Branch("jtEta", jtEta_, "jtEta[nJt]/F");
  skimTree_p->Branch("jtM", jtM_, "jtM[nJt]/F");
  skimTree_p->Branch("discr_csvV1", discr_csvV1_, "discr_csvV1[nJt]/F");
  skimTree_p->Branch("discr_csvV2", discr_csvV2_, "discr_csvV2[nJt]/F");

  skimTree_p->Branch("deltaPhi",deltaPhi_,"deltaPhi[nJt]/F");
  skimTree_p->Branch("deltaEta",deltaEta_,"deltaEta[nJt]/F");
  skimTree_p->Branch("ptMuJet",ptMuJet_,"ptMuJet[nJt]/F");
  skimTree_p->Branch("svtxM",svtxM_,"svtxM[nJt]/F");
  skimTree_p->Branch("svtxPt",svtxPt_,"svtxPt[nJt]/F");
  skimTree_p->Branch("discr_Prob",discr_Prob_,"discr_Prob[nJt]/F");
  skimTree_p->Branch("discr_ProbforB",discr_ProbforB_,"discr_ProbforB[nJt]/F");
  skimTree_p->Branch("discr_ProbforC",discr_ProbforC_,"discr_ProbforC[nJt]/F");
  skimTree_p->Branch("discr_ProbforBack",discr_ProbforBack_,"discr_ProbforBack[nJt]/F");
  skimTree_p->Branch("ptJP",ptJP_,"ptJP[nJt]/F");
  skimTree_p->Branch("R",R_,"R[nJt]/F");
  skimTree_p->Branch("deltaphiJJ",deltaphiJJ_,"deltaphiJJ[nJt]/F");
  //skimTree_p->Branch("CSVJP",CSVJP_,"CSVJP[nJt]/F");
  skimTree_p->Branch("svtxptJP",svtxptJP_,"svtxptJP[nJt]/F");
  skimTree_p->Branch("discrProbJP",discrProbJP_,"discrProbJP[nJt]/F");
  skimTree_p->Branch("discrcsvV2JP",discrcsvV2JP_,"discrcsvV2JP[nJt]/F");

  skimTree_p->Branch("discr_tcHighEff", discr_tcHighEff_, "discr_tcHighEff[nJt]/F");
  skimTree_p->Branch("discr_tcHighPur", discr_tcHighPur_, "discr_tcHighPur[nJt]/F");
  skimTree_p->Branch("refparton_flavorForB", refparton_flavorForB_, "refparton_flavorForB[nJt]/I");
  skimTree_p->Branch("nGen",&nGen_,"nGen/I");
  skimTree_p->Branch("genPdg",genPdg_,"genPdg[nGen]/I");
  skimTree_p->Branch("genPt",genPt_,"genPt[nGen]/F");
  skimTree_p->Branch("genPhi",genPhi_,"genPhi[nGen]/F");
  skimTree_p->Branch("genEta",genEta_,"genEta[nGen]/F");
  skimTree_p->Branch("genChg",genChg_,"genChg[nGen]/F");

  return;
}

void ReadTree(){
  if(skimTree_p == NULL){
    std::cout << "READTREE error; skimTree_p is NULL. return" << std::endl;
    return;
  }
  return;
}

#endif

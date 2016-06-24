#ifndef TopEMuTree_h
#define TopEMuTree_h

#include "TTree.h"

#include <iostream>

TTree* topEMuTree_p = 0;

UInt_t run_, lumi_;
ULong64_t evt_;
Int_t hiBin_;
Float_t vz_;

const int nDilepMax = 500;
Int_t nDilep_;
Float_t dilepPt_[nDilepMax];
Float_t dilepPhi_[nDilepMax];
Float_t dilepEta_[nDilepMax];
Float_t dilepM_[nDilepMax];
Float_t dilepDeltaPhi_[nDilepMax];
Float_t muIso_[nDilepMax];
Float_t eleIso_[nDilepMax];

const Int_t nMaxJets = 500;
Int_t nJt_;
Float_t jtPt_[nMaxJets];
Float_t jtPhi_[nMaxJets];
Float_t jtEta_[nMaxJets];
Float_t jtM_[nMaxJets];
Float_t discr_csvV1_[nMaxJets];

// List of branches
TBranch        *b_run;   //!
TBranch        *b_evt;   //!
TBranch        *b_lumi;   //!
TBranch        *b_hiBin;   //!
TBranch        *b_vz;   //!

TBranch        *b_nDilep;   //!
TBranch        *b_dilepPt;   //!
TBranch        *b_dilepPhi;   //!
TBranch        *b_dilepEta;   //!
TBranch        *b_dilepM;   //!
TBranch        *b_dilepDeltaPhi;   //!
TBranch        *b_muIso; //!
TBranch        *b_eleIso; //!

TBranch        *b_nJt;   //!
TBranch        *b_jtPt;   //!
TBranch        *b_jtPhi;   //!
TBranch        *b_jtEta;   //!
TBranch        *b_jtM;   //!
TBranch        *b_discr_csvV1;//!

void BookTree()
{
  if(topEMuTree_p == NULL){
    std::cout << "BOOKTREE error; topEMuTree_p is NULL. return" << std::endl;
    return;
  }

  topEMuTree_p->Branch("run", &run_, "run/i");
  topEMuTree_p->Branch("evt", &evt_, "evt/l");
  topEMuTree_p->Branch("lumi", &lumi_, "lumi/i");
  topEMuTree_p->Branch("hiBin", &hiBin_, "hiBin/I");
  topEMuTree_p->Branch("vz", &vz_, "vz/F");

  topEMuTree_p->Branch("nDilep", &nDilep_, "nDilep/I");
  topEMuTree_p->Branch("dilepPt", dilepPt_, "dilepPt[nDilep]/F");
  topEMuTree_p->Branch("dilepPhi", dilepPhi_, "dilepPhi[nDilep]/F");
  topEMuTree_p->Branch("dilepEta", dilepEta_, "dilepEta[nDilep]/F");
  topEMuTree_p->Branch("dilepM", dilepM_, "dilepM[nDilep]/F");
  topEMuTree_p->Branch("dilepDeltaPhi", dilepDeltaPhi_, "dilepDeltaPhi[nDilep]/F");
  topEMuTree_p->Branch("muIso",muIso_,"muIso[nDilep]/F");
  topEMuTree_p->Branch("eleIso",eleIso_,"eleIso[nDilep]/F");

  topEMuTree_p->Branch("nJt", &nJt_, "nJt/I");
  topEMuTree_p->Branch("jtPt", jtPt_, "jtPt[nJt]/F");
  topEMuTree_p->Branch("jtPhi", jtPhi_, "jtPhi[nJt]/F");
  topEMuTree_p->Branch("jtEta", jtEta_, "jtEta[nJt]/F");
  topEMuTree_p->Branch("jtM", jtM_, "jtM[nJt]/F");
  topEMuTree_p->Branch("discr_csvV1", discr_csvV1_, "discr_csvV1[nJt]/F");
 
  return;
}


void ReadTree()
{
  if(topEMuTree_p == NULL){
    std::cout << "READTREE error; topEMuTree_p is NULL. return" << std::endl;
    return;
  }

  return;
}

#endif

#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "TopEMuTree.h"
#include <string>
#include <vector>

const Float_t eleM = .000510998910;
const Float_t muM = .1056583715;

void anaSkimTree(TString strIn = "emuskim_0.root", const std::string outFileName = "topEmuSkim_0.root") {

  const int iDebug = 0;
  
  TFile *f = TFile::Open(strIn.Data());
  TTree *tr = dynamic_cast<TTree*>(f->Get("skimTree"));

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  topEMuTree_p = new TTree("topEMuTree", "topEMuTree");
  BookTree();

  //input tree variables
  // unsigned int run, lumi;
  // ULong64_t evt;
  int hiBin;
  float vz;
  
  int nLep;
  const int nLepMax = 4;
  int lepID[nLepMax];
  float lepPt[nLepMax];
  float lepPhi[nLepMax];
  float lepEta[nLepMax];
  int lepChg[nLepMax];
  float lepIso[nLepMax];

  const int nMaxJets = 500;
  int nJt;
  float jtPt[nMaxJets];
  float jtPhi[nMaxJets];
  float jtEta[nMaxJets];
  float jtM[nMaxJets];
  float discr_csvV1[nMaxJets];

  tr->SetBranchAddress("run", &run_);
  tr->SetBranchAddress("evt", &evt_);
  tr->SetBranchAddress("lumi", &lumi_);
  tr->SetBranchAddress("hiBin", &hiBin);
  tr->SetBranchAddress("vz", &vz);
  tr->SetBranchAddress("nLep", &nLep);
  tr->SetBranchAddress("lepID", lepID);
  tr->SetBranchAddress("lepPt", lepPt);
  tr->SetBranchAddress("lepPhi", lepPhi);
  tr->SetBranchAddress("lepEta", lepEta);
  tr->SetBranchAddress("lepChg", lepChg);
  tr->SetBranchAddress("lepIso", lepIso);
  tr->SetBranchAddress("nJt", &nJt);
  tr->SetBranchAddress("jtPt", jtPt);
  tr->SetBranchAddress("jtPhi", jtPhi);
  tr->SetBranchAddress("jtEta", jtEta);
  tr->SetBranchAddress("jtM", jtM);
  tr->SetBranchAddress("discr_csvV1", discr_csvV1);

  int emuPairs = 0;
  int emuPairsMassCut = 0;
  int emuPairsMassCutNJet2 = 0;
  
  int nEntries = (int)tr->GetEntries();
  for(int entry = 0; entry < nEntries; entry++){
    tr->GetEntry(entry);

    hiBin_ = hiBin;
    vz_    = vz;

    std::vector<int> indexEmu1;
    std::vector<int> indexEmu2;
    for(int ilep = 0; ilep<nLep;  ilep++) {
      for(int jlep = ilep+1; jlep<nLep;  jlep++) {
        if(lepID[ilep]==lepID[jlep]) continue;   //reject same flavor
        if(lepChg[ilep]==lepChg[jlep]) continue; //reject same sign

        //do lepton isolation cut
        if(hiBin<20 && lepIso[ilep]>0.58) continue;
        else if(hiBin>=20 && hiBin<60 && lepIso[ilep]>0.45) continue;
        else if(hiBin>=60 && hiBin<100 && lepIso[ilep]>0.3) continue;
        else if(hiBin>=100 && hiBin<140 && lepIso[ilep]>0.24) continue;
        else if(hiBin>=140 && lepIso[ilep]>0.18) continue;
        
        if(iDebug) Printf("%d event %llu Found emu pair: %d %d",entry,evt_,ilep,jlep);
        if(lepID[ilep] == 11) {//first lepton is electron 
          indexEmu1.push_back(ilep);
          indexEmu2.push_back(jlep);
          if(iDebug) Printf("Ele pt=%f eta=%f phi=%f ch=%d",lepPt[ilep],lepEta[ilep],lepPhi[ilep],lepChg[ilep]);
          if(iDebug) Printf("Muon pt=%f eta=%f phi=%f ch=%d",lepPt[jlep],lepEta[jlep],lepPhi[jlep],lepChg[jlep]);
        } else {//first lepton is muon
          indexEmu1.push_back(jlep);
          indexEmu2.push_back(ilep);
          if(iDebug) Printf("Ele pt=%f eta=%f phi=%f ch=%d",lepPt[jlep],lepEta[jlep],lepPhi[jlep],lepChg[jlep]);
          if(iDebug) Printf("Muon pt=%f eta=%f phi=%f ch=%d",lepPt[ilep],lepEta[ilep],lepPhi[ilep],lepChg[ilep]);
        }
        ++emuPairs;
      }
    }

    //make dilepton pairs
    std::vector<TLorentzVector> dileptons;
    std::vector<TLorentzVector> electrons;
    std::vector<TLorentzVector> muons;
    std::vector<int> index_muons;
    std::vector<int> index_electrons;
    
    int npairs = (int)indexEmu1.size();
    for(int ilep = 0; ilep<npairs; ++ilep) {
      int ie = indexEmu1[ilep];
      int im = indexEmu2[ilep];
      TLorentzVector ele;
      ele.SetPtEtaPhiM(lepPt[ie],lepEta[ie],lepPhi[ie],eleM);

      TLorentzVector muon;
      muon.SetPtEtaPhiM(lepPt[im],lepEta[im],lepPhi[im],muM);

      TLorentzVector dilepton = ele + muon;
      if(dilepton.M()>20.) {
        dileptons.push_back(dilepton);
        electrons.push_back(ele);
        muons.push_back(muon);
        if(iDebug) Printf("dilepton pair %d mass: %f pt: %f",ilep,dilepton.M(),dilepton.Pt());
        index_muons.push_back(im);
        index_electrons.push_back(ie);
        ++emuPairsMassCut;
      }
    }

    if(dileptons.size()<1) continue;

    std::vector<int> indexJets;
    int njets = 0;
    for(int ij = 0; ij<nJt; ++ij) {
      if(jtPt[ij]<30.) continue;
      if(fabs(jtEta[ij])>2.) continue;

      //reject jets close to selected leptons
      bool accept = true;
      for(int ilep=0; ilep<(int)dileptons.size(); ++ilep) {
        double drToMuon = sqrt(pow(acos(cos(jtPhi[ij]-muons[ilep].Phi())),2)+pow(jtEta[ij]-muons[ilep].Eta(),2));
        double drToEle = sqrt(pow(acos(cos(jtPhi[ij]-electrons[ilep].Phi())),2)+pow(jtEta[ij]-electrons[ilep].Eta(),2));
        if(drToMuon<0.3 || drToEle<0.3) accept = false;
      }
      if(accept) {
        indexJets.push_back(ij);
        ++njets;
      }
    }//jet loop
    
    if(iDebug) Printf("njets: %d",njets);
    if(njets>1) {

      nDilep_ = (int)dileptons.size();
      for(int ilep=0; ilep<(int)dileptons.size(); ++ilep) {
        dilepPt_[ilep] = dileptons[ilep].Pt();
        dilepEta_[ilep] = dileptons[ilep].Eta();
        dilepPhi_[ilep] = dileptons[ilep].Phi();
        dilepM_[ilep] = dileptons[ilep].M();

        dilepDeltaPhi_[ilep] = fabs(TVector2::Phi_mpi_pi(electrons[ilep].Phi()-muons[ilep].Phi()))/TMath::Pi();
        muIso_[ilep] = lepIso[index_muons[ilep]];
        eleIso_[ilep] = lepIso[index_electrons[ilep]];
      }
      nJt_ = njets;
      for(int ij = 0; ij<(int)indexJets.size(); ++ij) {
        jtPt_[ij] = jtPt[indexJets[ij]];
        jtEta_[ij] = jtEta[indexJets[ij]];
        jtPhi_[ij] = jtPhi[indexJets[ij]];
        jtM_[ij] = jtM[indexJets[ij]];
        discr_csvV1_[ij] = discr_csvV1[indexJets[ij]];
      }
      
      ++emuPairsMassCutNJet2;
      topEMuTree_p->Fill();
    }
  }

  Printf("Total emu pairs: %d with mass cut: %d with NJet cut: %d",emuPairs,emuPairsMassCut,emuPairsMassCutNJet2);

  outFile_p->cd();
  topEMuTree_p->Write();
  outFile_p->Close();
  delete outFile_p;
  
}

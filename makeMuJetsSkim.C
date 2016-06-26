#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "LepJetsSkimTree.h"
#include <string>
#include <vector>

#include "ForestTreeHeaders/ForestMuons.h"

const bool isDebug = false;

const float jetPtCut = 25.;
const float jetEtaCut = 3.;

const float lepPtCut = 15;

// FIXME: Need to check lepton selections for PbPb
const float muEtaCut = 2.1;
const float muPtCut = 15;

const float muChi2NDFCut   = 10;
const float muInnerD0Cut   = 0.2;
const float muInnerDzCut   = 0.5;
const int   muMuonHitsCut  = 0;
const int   muStationsCut  = 1;
const int   muTrkLayersCut = 5;
const int   muPixelHitsCut = 0;

double calcLeptonIsolation(float lepPt, float lepEta, float lepPhi, std::vector<float> *pfPt, std::vector<float> *pfEta, std::vector<float> *pfPhi);

void makeMuJetsSkim(const std::string outFileName = "", const std::string inFileName = "")
{
  if(!strcmp(inFileName.c_str(), "")){
    std::cout << "No inputs specified. return" << std::endl;
    return;
  }

  if(isDebug) std::cout << __LINE__ << std::endl;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  skimTree_p = new TTree("skimTree", "skimTree");
  BookTree();

  std::vector<std::string>* inFileNames_p = new std::vector<std::string>;
  inFileNames_p->push_back(inFileName);

  TChain *lepTree_p = new TChain("ggHiNtuplizer/EventTree");
  TChain *jetTree_p = new TChain("akCs2PFJetAnalyzer/t");
  TChain *hiTree_p = new TChain("hiEvtAnalyzer/HiTree");
  TChain *hltTree_p = new TChain("hltanalysis/HltTree");
  TChain *pfTree_p = new TChain("pfcandAnalyzerCS/pfTree");
  TChain *skimAnaTree_p = new TChain("skimanalysis/HltTree");
  
  const int nFiles = (int)inFileNames_p->size();

  for(int fileIter = 0; fileIter < nFiles; fileIter++){
    std::cout << "On file: " << fileIter << "/" << nFiles << "  " << inFileNames_p->at(fileIter).c_str()  << std::endl;
    lepTree_p->Add(inFileNames_p->at(fileIter).c_str());
    jetTree_p->Add(inFileNames_p->at(fileIter).c_str());
    hiTree_p->Add(inFileNames_p->at(fileIter).c_str());
    hltTree_p->Add(inFileNames_p->at(fileIter).c_str());
    pfTree_p->Add(inFileNames_p->at(fileIter).c_str());
    skimAnaTree_p->Add(inFileNames_p->at(fileIter).c_str());
  }
  
  ForestMuons fForestMu;
      
  const int maxJets = 5000;
  int           nref;
  float         jtpt[maxJets];   //[nref]
  float         jteta[maxJets];   //[nref]
  float         jtphi[maxJets];   //[nref]
  float         jtm[maxJets];   //[nref]
  float         discr_csvV1[maxJets]; //[nref]

  //pf particles pfId, pfPt, pfEta, pfPhi
  std::vector<int>           *pfId = 0;
  std::vector<float>         *pfPt = 0;
  std::vector<float>         *pfEta = 0;
  std::vector<float>         *pfPhi = 0;
    
  int trig = 0;

  //event selections
  int phfCoincFilter = 1;
  int HBHENoiseFilterResult = 1;
  int pprimaryVertexFilter = 1;
  int pcollisionEventSelection = 1;

  lepTree_p->SetBranchStatus("*", 0);
  lepTree_p->SetBranchStatus("mu*", 1);
      
  lepTree_p->SetBranchAddress("muPt", &fForestMu.muPt);
  lepTree_p->SetBranchAddress("muPhi", &fForestMu.muPhi);
  lepTree_p->SetBranchAddress("muEta", &fForestMu.muEta);
  lepTree_p->SetBranchAddress("muCharge", &fForestMu.muCharge);
  lepTree_p->SetBranchAddress("muChi2NDF", &fForestMu.muChi2NDF);
  lepTree_p->SetBranchAddress("muInnerD0", &fForestMu.muInnerD0);
  lepTree_p->SetBranchAddress("muInnerDz", &fForestMu.muInnerDz);
  lepTree_p->SetBranchAddress("muMuonHits", &fForestMu.muMuonHits);
  lepTree_p->SetBranchAddress("muStations", &fForestMu.muStations);
  lepTree_p->SetBranchAddress("muTrkLayers", &fForestMu.muTrkLayers);
  lepTree_p->SetBranchAddress("muPixelHits", &fForestMu.muPixelHits);    
 
  jetTree_p->SetBranchStatus("*", 0);
  jetTree_p->SetBranchStatus("nref", 1);
  jetTree_p->SetBranchStatus("jtpt", 1);
  jetTree_p->SetBranchStatus("jtphi", 1);
  jetTree_p->SetBranchStatus("jteta", 1);
  jetTree_p->SetBranchStatus("jtm", 1);
  jetTree_p->SetBranchStatus("discr_csvV1", 1);
        
  jetTree_p->SetBranchAddress("nref", &nref);
  jetTree_p->SetBranchAddress("jtpt", jtpt);
  jetTree_p->SetBranchAddress("jtphi", jtphi);
  jetTree_p->SetBranchAddress("jteta", jteta);
  jetTree_p->SetBranchAddress("jtm", jtm);
  jetTree_p->SetBranchAddress("discr_csvV1", discr_csvV1);
    
  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("run", 1);
  hiTree_p->SetBranchStatus("evt", 1);
  hiTree_p->SetBranchStatus("lumi", 1);
  hiTree_p->SetBranchStatus("hiBin", 1);
  hiTree_p->SetBranchStatus("vz", 1);
    
  hiTree_p->SetBranchAddress("run", &run_);
  hiTree_p->SetBranchAddress("evt", &evt_);
  hiTree_p->SetBranchAddress("lumi", &lumi_);
  hiTree_p->SetBranchAddress("hiBin", &hiBin_);
  hiTree_p->SetBranchAddress("vz", &vz_);

  pfTree_p->SetBranchAddress("pfId", &pfId);
  pfTree_p->SetBranchAddress("pfPt", &pfPt);
  pfTree_p->SetBranchAddress("pfEta", &pfEta);
  pfTree_p->SetBranchAddress("pfPhi", &pfPhi);
    
  hltTree_p->SetBranchStatus("HLT_HIL2Mu15_v2",1);
  hltTree_p->SetBranchAddress("HLT_HIL2Mu15_v2",&trig);

  skimAnaTree_p->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
  skimAnaTree_p->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult);
  skimAnaTree_p->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
  skimAnaTree_p->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
    
  if(isDebug) std::cout << __LINE__ << std::endl;
    
  if(isDebug) std::cout << __LINE__ << std::endl;
    
  int nEntries = (int)lepTree_p->GetEntries();
  //nEntries = 100;
  int entryDiv = ((int)(nEntries/20));
    
  if(isDebug) std::cout << __LINE__ << std::endl;
    
  for(int entry = 0; entry < nEntries; entry++){
    if(isDebug) std::cout << __LINE__ << std::endl;
      
    if(entry%entryDiv == 0 && nEntries >= 10000) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
    if(isDebug) std::cout << __LINE__ << std::endl;

    lepTree_p->GetEntry(entry);
    jetTree_p->GetEntry(entry);
    hiTree_p->GetEntry(entry);
    hltTree_p->GetEntry(entry);
    pfTree_p->GetEntry(entry);
    skimAnaTree_p->GetEntry(entry);
      
    if(!trig) continue;
    if(!phfCoincFilter) continue;
    if(!HBHENoiseFilterResult) continue;
    if(!pcollisionEventSelection) continue;
    if(!pprimaryVertexFilter) continue;
      
    if(TMath::Abs(vz_) > 15) continue;
      
    if(isDebug) std::cout << __LINE__ << std::endl;

    float tempMuPt_[nLep];
    float tempMuPhi_[nLep];
    float tempMuEta_[nLep];
    int tempMuChg_[nLep];
    float tempMuIso_[nLep];     

    for(int lepIter = 0; lepIter < nLep; lepIter++){
      lepPt_[lepIter] = -999;
      lepPhi_[lepIter] = -999;
      lepEta_[lepIter] = -999;
      lepChg_[lepIter] = -999;
      lepID_[lepIter] = -999;
      lepIso_[lepIter] = -999;
    }
       
    for(int lepIter = 0; lepIter < 2; lepIter++){
      tempMuPt_[lepIter] = -999;
      tempMuPhi_[lepIter] = -999;
      tempMuEta_[lepIter] = -999;
      tempMuChg_[lepIter] = -999;
      tempMuIso_[lepIter] = -999;
    }
       
    for(int ij = 0; ij<nMaxJets; ++ij) {
      jtPt_[ij] = -999.;
      jtEta_[ij] = -999.;
      jtPhi_[ij] = -999.;
      jtM_[ij] = -999.;
      discr_csvV1_[ij] = -999.;
    }
      
    if(isDebug) std::cout << __LINE__ << std::endl;
      
    //Find two leading muons
    for(unsigned int muIter = 0; muIter < fForestMu.muPt->size(); ++muIter) {
      if(TMath::Abs(fForestMu.muEta->at(muIter)) > muEtaCut) continue;
      if(fForestMu.muPt->at(muIter) < muPtCut) continue;
	
      if(fForestMu.muChi2NDF->at(muIter) >= muChi2NDFCut) continue;
      if(TMath::Abs(fForestMu.muInnerD0->at(muIter)) > muInnerD0Cut) continue;
      if(TMath::Abs(fForestMu.muInnerDz->at(muIter)) > muInnerDzCut) continue;
      if(fForestMu.muMuonHits->at(muIter) <= muMuonHitsCut) continue;
      if(fForestMu.muStations->at(muIter) <= muStationsCut) continue;
      if(fForestMu.muTrkLayers->at(muIter) <= muTrkLayersCut) continue;
      if(fForestMu.muPixelHits->at(muIter) <= muPixelHitsCut) continue;

      if(fForestMu.muPt->at(muIter) > lepPt_[0]){
        tempMuPt_[1]  = tempMuPt_[0];
        tempMuPhi_[1] = tempMuPhi_[0];
        tempMuEta_[1] = tempMuEta_[0];
        tempMuChg_[1] = tempMuChg_[0];
        tempMuIso_[1] = tempMuIso_[0]; 
 
        tempMuPt_[0]  = fForestMu.muPt->at(muIter);
        tempMuPhi_[0] = fForestMu.muPhi->at(muIter);
        tempMuEta_[0] = fForestMu.muEta->at(muIter);
        tempMuChg_[0] = fForestMu.muCharge->at(muIter);
        tempMuIso_[0] = calcLeptonIsolation(fForestMu.muPt->at(muIter),fForestMu.muEta->at(muIter),fForestMu.muPhi->at(muIter),pfPt,pfEta,pfPhi); //iso;
      }
      else if(fForestMu.muPt->at(muIter) > tempMuPt_[1]){
        tempMuPt_[1]  = fForestMu.muPt->at(muIter);
        tempMuPhi_[1] = fForestMu.muPhi->at(muIter);
        tempMuEta_[1] = fForestMu.muEta->at(muIter);
        tempMuChg_[1] = fForestMu.muCharge->at(muIter);
        tempMuIso_[1] = calcLeptonIsolation(fForestMu.muPt->at(muIter),fForestMu.muEta->at(muIter),fForestMu.muPhi->at(muIter),pfPt,pfEta,pfPhi);//iso;
      }
    }

    //store muons in out tree
    int lepIter = 0;
    for(int muIter = 0; muIter < 2; muIter++){
      if(tempMuPt_[muIter]<0.) continue;
      lepPt_[lepIter] = tempMuPt_[muIter];
      lepPhi_[lepIter] = tempMuPhi_[muIter];
      lepEta_[lepIter] = tempMuEta_[muIter];
      lepChg_[lepIter] = tempMuChg_[muIter];
      lepID_[lepIter] = muID;
      lepIso_[lepIter] = tempMuIso_[muIter];
      ++lepIter;
    }
    if(lepIter<1) continue;
    nLep_ = lepIter;

    int njets = 0;
    for(int jetIter = 0; jetIter < nref; jetIter++){
      if(jtpt[jetIter]<jetPtCut) continue;
      if(fabs(jteta[jetIter])>jetEtaCut) continue;
      jtPt_[njets]  = jtpt[jetIter];
      jtEta_[njets] = jteta[jetIter];
      jtPhi_[njets] = jtphi[jetIter];
      jtM_[njets]   = jtm[jetIter]; 
      discr_csvV1_[njets] = discr_csvV1[jetIter];
      ++njets;
    }
    nJt_ = njets;
    if(nJt_<4) continue; //need at least 2 b jets (t->Wb) and 2 light jets (W->qqbar)
    
    skimTree_p->Fill();
    
  }//entries loop

  
  outFile_p->cd();
  TNamed pathStr1("pathStr1", outFileName.c_str());
  pathStr1.Write("", TObject::kOverwrite);
  TNamed pathStr2("pathStr2", inFileName.c_str());
  pathStr2.Write("", TObject::kOverwrite);
  skimTree_p->Write("", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  return;
}

double calcLeptonIsolation(float lepPt, float lepEta, float lepPhi, std::vector<float> *pfPt, std::vector<float> *pfEta, std::vector<float> *pfPhi) {
  //calculate lepton isolation from pf candidates.
  //Isolation cone R=0.3
  //excluding pf candidates at a distance less than 0.03 from lepton
  
  double conePt = 0.;
  for(unsigned int i = 0; i<pfPt->size(); ++i) {

    double deltaR = sqrt(pow(acos(cos(lepPhi-pfPhi->at(i))),2)+pow(lepEta-pfEta->at(i),2));

    if(deltaR<0.03 || deltaR>0.3) continue;

    conePt+=pfPt->at(i);
  }
  double relIso = conePt;
  if(lepPt>0.) relIso = conePt/lepPt;
  
  return relIso;
}

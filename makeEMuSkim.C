#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "EMuSkimTree.h"
#include <string>
#include <vector>

#include "ForestTreeHeaders/ForestElectrons.h"
#include "ForestTreeHeaders/ForestMuons.h"

#include "Helpers/EnergyRegression.h"

const bool isDebug = false;

const float jetPtCut = 25.;

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

//PbPb: https://twiki.cern.ch/twiki/bin/view/CMS/ElectronPbPb5TeV#3_Selection_for_different_centra
//pp:   https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns
//WARNING: this code has PbPb selections activated
const float eleEtaCut = 2.4;
const float elePtCut = 15;

const float barrelEndcapEta = 1.479;
const int nBarrelEndcap = 2; // barrel = 0, endcap = 1
//For PbPb electron ID is centrality dependent
//centrality bins 0-20%, 20-40%, 40-70%, 70-100%
const int nCentEleId = 4;
double centMinEleId[4] = {0.,20.,40.,70.};
double centMaxEleId[4] = {20.,40.,70.,100.};

float eleSigmaIEtaIEta_VetoCut[nBarrelEndcap][nCentEleId];
float eleDEtaIn_VetoCut[nBarrelEndcap][nCentEleId];
float eleDPhiIn_VetoCut[nBarrelEndcap][nCentEleId];
float eleHOverE_VetoCut[nBarrelEndcap][nCentEleId];
float eleRelIsoWithEA_VetoCut[nBarrelEndcap][nCentEleId];
float eleOOEmooP_VetoCut[nBarrelEndcap][nCentEleId];
float eleD0_VetoCut[nBarrelEndcap][nCentEleId];
float eleDZ_VetoCut[nBarrelEndcap][nCentEleId];
float eleMissingInnerHits_VetoCut[nBarrelEndcap][nCentEleId];
float eleEoverPInv_VetoCut[nBarrelEndcap][nCentEleId];

// const float eleSigmaIEtaIEta_VetoCut[nBarrelEndcap] = {0.0114, 0.0352};
// const float eleDEtaIn_VetoCut[nBarrelEndcap] = {0.0152, 0.0113};
// const float eleDPhiIn_VetoCut[nBarrelEndcap] = {0.216, 0.237};
// const float eleHOverE_VetoCut[nBarrelEndcap] = {0.181, 0.116};
// const float eleRelIsoWithEA_VetoCut[nBarrelEndcap] = {0.126, 0.144};
// const float eleOOEmooP_VetoCut[nBarrelEndcap] = {0.207, 0.174};
// const float eleD0_VetoCut[nBarrelEndcap] = {0.0564, 0.222};
// const float eleDZ_VetoCut[nBarrelEndcap] = {0.472, 0.921};
// const float eleMissingInnerHits_VetoCut[nBarrelEndcap] = {2, 3};

int getCentBinEleId(double cent) {

  int centBin = -1;
  for(int i = 0; i<nCentEleId; ++i) {
    if(cent>=centMinEleId[i] && cent<centMaxEleId[i]) centBin = i;
  }
  return centBin;
}

double calcLeptonIsolation(float lepPt, float lepEta, float lepPhi, std::vector<float> *pfPt, std::vector<float> *pfEta, std::vector<float> *pfPhi);

void makeEMuSkim(const std::string outFileName = "", const std::string inFileName = "")
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
  // if(strcmp(inFileName.c_str(), "") != 0) inFileNames_p->push_back(inFileName);

  //  Printf("Create EnergyRegression object");
  EnergyRegression energyRegression;
  //Printf("Create EnergyRegression object done. calling initreader");
  energyRegression.initreader();
  
  //initialize electron Id variables [barrel/endcap][centBin]
  eleSigmaIEtaIEta_VetoCut[0][0] = 0.01325;
  eleSigmaIEtaIEta_VetoCut[1][0] = 0.04272;
  eleSigmaIEtaIEta_VetoCut[0][1] = 0.01098;
  eleSigmaIEtaIEta_VetoCut[1][1] = 0.03569;
  eleSigmaIEtaIEta_VetoCut[0][2] = 0.01078;
  eleSigmaIEtaIEta_VetoCut[1][2] = 0.03155;
  eleSigmaIEtaIEta_VetoCut[0][3] = 0.01038;
  eleSigmaIEtaIEta_VetoCut[1][3] = 0.02955;

  eleDEtaIn_VetoCut[0][0] = 0.04524;
  eleDEtaIn_VetoCut[1][0] = 0.42269;
  eleDEtaIn_VetoCut[0][1] = 0.02799;
  eleDEtaIn_VetoCut[1][1] = 0.01592;
  eleDEtaIn_VetoCut[0][2] = 0.01275;
  eleDEtaIn_VetoCut[1][2] = 0.01074;
  eleDEtaIn_VetoCut[0][3] = 0.00595;
  eleDEtaIn_VetoCut[1][3] = 0.00927;

  eleDPhiIn_VetoCut[0][0] = 0.15133;
  eleDPhiIn_VetoCut[1][0] = 0.33656;
  eleDPhiIn_VetoCut[0][1] = 0.10697;
  eleDPhiIn_VetoCut[1][1] = 0.18786;
  eleDPhiIn_VetoCut[0][2] = 0.11228;
  eleDPhiIn_VetoCut[1][2] = 0.12940;
  eleDPhiIn_VetoCut[0][3] = 0.22180;
  eleDPhiIn_VetoCut[1][3] = 0.20464;

  eleHOverE_VetoCut[0][0] = 0.12879;
  eleHOverE_VetoCut[1][0] = 0.14855;
  eleHOverE_VetoCut[0][1] = 0.09844;
  eleHOverE_VetoCut[1][1] = 0.11125;
  eleHOverE_VetoCut[0][2] = 0.02355;
  eleHOverE_VetoCut[1][2] = 0.05202;
  eleHOverE_VetoCut[0][3] = 0.02997;
  eleHOverE_VetoCut[1][3] = 0.01670;

  eleD0_VetoCut[0][0] = 0.18115;
  eleD0_VetoCut[1][0] = 0.12069;
  eleD0_VetoCut[0][1] = 0.06520;
  eleD0_VetoCut[1][1] = 0.16610;
  eleD0_VetoCut[0][2] = 0.05574;
  eleD0_VetoCut[1][2] = 0.14651;
  eleD0_VetoCut[0][3] = 0.05396;
  eleD0_VetoCut[1][3] = 0.12990;

  eleDZ_VetoCut[0][0] = 0.11313;
  eleDZ_VetoCut[1][0] = 0.28022;
  eleDZ_VetoCut[0][1] = 0.06983;
  eleDZ_VetoCut[1][1] = 0.24015;
  eleDZ_VetoCut[0][2] = 0.02011;
  eleDZ_VetoCut[1][2] = 0.18170;
  eleDZ_VetoCut[0][3] = 0.06513;
  eleDZ_VetoCut[1][3] = 0.24262;

  eleMissingInnerHits_VetoCut[0][0] = 1.00005;
  eleMissingInnerHits_VetoCut[1][0] = 1.00005;
  eleMissingInnerHits_VetoCut[0][1] = 1.00005;
  eleMissingInnerHits_VetoCut[1][1] = 1.00005;
  eleMissingInnerHits_VetoCut[0][2] = 1.00005;
  eleMissingInnerHits_VetoCut[1][2] = 1.00005;
  eleMissingInnerHits_VetoCut[0][3] = 1.00005;
  eleMissingInnerHits_VetoCut[1][3] = 1.00005;

  eleEoverPInv_VetoCut[0][0] = 0.20965;
  eleEoverPInv_VetoCut[1][0] = 0.20556;
  eleEoverPInv_VetoCut[0][1] = 0.12291;
  eleEoverPInv_VetoCut[1][1] = 0.29670;
  eleEoverPInv_VetoCut[0][2] = 0.30592;
  eleEoverPInv_VetoCut[1][2] = 0.20473;
  eleEoverPInv_VetoCut[0][3] = 0.23732;
  eleEoverPInv_VetoCut[1][3] = 0.11148;
  
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
    // TFile *inFile_p = TFile::Open(inFileNames_p->at(fileIter).c_str(), "READ");
    // TTree *lepTree_p = dynamic_cast<TTree*>(inFile_p->Get("ggHiNtuplizer/EventTree"));
    // TTree *jetTree_p = dynamic_cast<TTree*>(inFile_p->Get("akCs2PFJetAnalyzer/t"));
    // TTree *hiTree_p = dynamic_cast<TTree*>(inFile_p->Get("hiEvtAnalyzer/HiTree"));
    // TTree *hltTree_p = dynamic_cast<TTree*>(inFile_p->Get("hltanalysis/HltTree"));
    // TTree *pfTree_p = dynamic_cast<TTree*>(inFile_p->Get("pfcandAnalyzerCS/pfTree"));
    // TTree *skimAnaTree_p = dynamic_cast<TTree*>(inFile_p->Get("skimanalysis/HltTree"));

    ForestMuons fForestMu;
    ForestElectrons fForestEle;
    
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
    
    lepTree_p->SetBranchStatus("mu*", 1);
    // lepTree_p->SetBranchStatus("muPt", 1);
    // lepTree_p->SetBranchStatus("muPhi", 1);
    // lepTree_p->SetBranchStatus("muEta", 1);
    // lepTree_p->SetBranchStatus("muCharge", 1);
    // lepTree_p->SetBranchStatus("muChi2NDF", 1);
    // lepTree_p->SetBranchStatus("muInnerD0", 1);
    // lepTree_p->SetBranchStatus("muInnerDz", 1);
    // lepTree_p->SetBranchStatus("muMuonHits", 1);
    // lepTree_p->SetBranchStatus("muStations", 1);
    // lepTree_p->SetBranchStatus("muTrkLayers", 1);
    // lepTree_p->SetBranchStatus("muPixelHits", 1);
    
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
    
    
    lepTree_p->SetBranchStatus("ele*", 1);
    
    lepTree_p->SetBranchAddress("elePt", &fForestEle.elePt);
    lepTree_p->SetBranchAddress("elePhi", &fForestEle.elePhi);
    lepTree_p->SetBranchAddress("eleEta", &fForestEle.eleEta);
    lepTree_p->SetBranchAddress("eleCharge", &fForestEle.eleCharge);
    lepTree_p->SetBranchAddress("eleSigmaIEtaIEta", &fForestEle.eleSigmaIEtaIEta);
    lepTree_p->SetBranchAddress("eledEtaAtVtx", &fForestEle.eledEtaAtVtx);
    lepTree_p->SetBranchAddress("eledPhiAtVtx", &fForestEle.eledPhiAtVtx);
    lepTree_p->SetBranchAddress("eleHoverE", &fForestEle.eleHoverE);
    lepTree_p->SetBranchAddress("eleD0", &fForestEle.eleD0);
    lepTree_p->SetBranchAddress("eleDz", &fForestEle.eleDz);
    lepTree_p->SetBranchAddress("eleEoverPInv", &fForestEle.eleEoverPInv);

    lepTree_p->SetBranchAddress("eleSCPhi", &fForestEle.eleSCPhi);
    lepTree_p->SetBranchAddress("eleSCEta", &fForestEle.eleSCEta);
    lepTree_p->SetBranchAddress("eleSigmaIEtaIEta_2012", &fForestEle.eleSigmaIEtaIEta_2012);
    lepTree_p->SetBranchAddress("eleSigmaIPhiIPhi", &fForestEle.eleSigmaIPhiIPhi);
    lepTree_p->SetBranchAddress("eleSCEtaWidth", &fForestEle.eleSCEtaWidth);
    lepTree_p->SetBranchAddress("eleSCPhiWidth", &fForestEle.eleSCPhiWidth);
        
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

      int centEleId = getCentBinEleId((double)hiBin_/2.);
            
      float tempMuPt_[nLep];
      float tempMuPhi_[nLep];
      float tempMuEta_[nLep];
      int tempMuChg_[nLep];
      float tempMuIso_[nLep];     
 
      float tempElePt_[nLep];
      float tempElePhi_[nLep];
      float tempEleEta_[nLep];
      int tempEleChg_[nLep];
      float tempEleIso_[nLep];            

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
	
	tempElePt_[lepIter] = -999;
	tempElePhi_[lepIter] = -999;
	tempEleEta_[lepIter] = -999;
	tempEleChg_[lepIter] = -999;
        tempEleIso_[lepIter] = -999;
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

      //Find two leading electrons
      for(unsigned int eleIter = 0; eleIter < fForestEle.elePt->size(); ++eleIter) {
	if(TMath::Abs(fForestEle.eleEta->at(eleIter)) > eleEtaCut) continue;
	if(fForestEle.elePt->at(eleIter) < elePtCut) continue;	
        
	int eleEtaCutPos = 0;
	if(TMath::Abs(fForestEle.eleEta->at(eleIter)) > barrelEndcapEta) eleEtaCutPos = 1;

   	if(fForestEle.eleSigmaIEtaIEta->at(eleIter) > eleSigmaIEtaIEta_VetoCut[eleEtaCutPos][centEleId]) continue;
	if(TMath::Abs(fForestEle.eledEtaAtVtx->at(eleIter)) > eleDEtaIn_VetoCut[eleEtaCutPos][centEleId]) continue;
	if(TMath::Abs(fForestEle.eledPhiAtVtx->at(eleIter)) > eleDPhiIn_VetoCut[eleEtaCutPos][centEleId]) continue;
	if(fForestEle.eleHoverE->at(eleIter) > eleHOverE_VetoCut[eleEtaCutPos][centEleId]) continue;
	if(TMath::Abs(fForestEle.eleD0->at(eleIter)) > eleD0_VetoCut[eleEtaCutPos][centEleId]) continue;
	if(TMath::Abs(fForestEle.eleDz->at(eleIter)) > eleDZ_VetoCut[eleEtaCutPos][centEleId]) continue;
        if(TMath::Abs(fForestEle.eleEoverPInv->at(eleIter)) > eleEoverPInv_VetoCut[eleEtaCutPos][centEleId]) continue;

        if(isDebug) std::cout << __LINE__ << std::endl;
        float ptEle = fForestEle.elePt->at(eleIter);
        //Printf("ptEle: %f",ptEle);
        float eleCorr = energyRegression.ElectronRegressionTMVA(fForestEle.eleSCPhi->at(eleIter),
                                                                fForestEle.eleSCEta->at(eleIter),
                                                                fForestEle.eleSigmaIEtaIEta->at(eleIter),
                                                                fForestEle.eleSigmaIPhiIPhi->at(eleIter),
                                                                fForestEle.eleSCEtaWidth->at(eleIter),
                                                                fForestEle.eleSCPhiWidth->at(eleIter),
                                                                fForestEle.eleHoverE->at(eleIter));

       
        fForestEle.elePt->at(eleIter) = ptEle * eleCorr;
        //Printf("ptEle: %f  eleCorr: %f elePt: %f",ptEle,eleCorr,fForestEle.elePt->at(eleIter));
        
     	if(fForestEle.elePt->at(eleIter) > tempElePt_[0]){
	  tempElePt_[1] = tempElePt_[0];
	  tempElePhi_[1] = tempElePhi_[0];
	  tempEleEta_[1] = tempEleEta_[0];
	  tempEleChg_[1] = tempEleChg_[0];
          tempEleIso_[1] = tempEleIso_[0];	 
 
	  tempElePt_[0] = fForestEle.elePt->at(eleIter);
	  tempElePhi_[0] = fForestEle.elePhi->at(eleIter);
	  tempEleEta_[0] = fForestEle.eleEta->at(eleIter);
	  tempEleChg_[0] = fForestEle.eleCharge->at(eleIter);
          tempEleIso_[0] = calcLeptonIsolation(fForestEle.elePt->at(eleIter),fForestEle.eleEta->at(eleIter),fForestEle.elePhi->at(eleIter),pfPt,pfEta,pfPhi); //iso;
	}
	else if(fForestEle.elePt->at(eleIter) > tempElePt_[1]){
	  tempElePt_[1]  = fForestEle.elePt->at(eleIter);
	  tempElePhi_[1] = fForestEle.elePhi->at(eleIter);
	  tempEleEta_[1] = fForestEle.eleEta->at(eleIter);
	  tempEleChg_[1] = fForestEle.eleCharge->at(eleIter);
          tempEleIso_[1] = calcLeptonIsolation(fForestEle.elePt->at(eleIter),fForestEle.eleEta->at(eleIter),fForestEle.elePhi->at(eleIter),pfPt,pfEta,pfPhi);//iso;
	}
      }

      //store electrons and muons in out tree
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
      for(int eleIter = 0; eleIter < 2; eleIter++){
        if(tempElePt_[eleIter]<0.) continue;
        lepPt_[lepIter] = tempElePt_[eleIter];
        lepPhi_[lepIter] = tempElePhi_[eleIter];
        lepEta_[lepIter] = tempEleEta_[eleIter];
        lepChg_[lepIter] = tempEleChg_[eleIter];
        lepID_[lepIter] = eleID;
        lepIso_[lepIter] = tempEleIso_[eleIter];
        ++lepIter;
      }
      if(lepIter<2) continue;
      nLep_ = lepIter;

      int njets = 0;
      for(int jetIter = 0; jetIter < nref; jetIter++){
        if(jtpt[jetIter]<jetPtCut) continue;
        jtPt_[njets]  = jtpt[jetIter];
        jtEta_[njets] = jteta[jetIter];
        jtPhi_[njets] = jtphi[jetIter];
        jtM_[njets]   = jtm[jetIter]; 
        discr_csvV1_[njets] = discr_csvV1[jetIter];
        ++njets;
      }
      nJt_ = njets;

      skimTree_p->Fill();
      
    }//entries loop

    //    inFile_p->Close();
    // delete inFile_p;    }

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

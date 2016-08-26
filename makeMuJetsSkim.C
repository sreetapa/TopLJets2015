#include "TFile.h"
#include "TObjArray.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "LepJetsSkimTree.h"
//#include "drawMuJetsSkimControlPlots.C"
#include <string>
#include <vector>
#include "ForestTreeHeaders/ForestMuons.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH1.h"

// Jet and lepton selection
float jetPtCut  = 30.;
const float jetEtaCut = 2.4;
const int   minNJets  = 1;   //Note: you will have to change this to make control distributions for events with less than 4 jets
const float muEtaCut = 2.1;
const float muPtCut  = 15.;
const float muChi2NDFCut   = 10;
const float muInnerD0Cut   = 0.2;
const float muInnerDzCut   = 20.;//0.5;
const int   muMuonHitsCut  = 0;
const int   muStationsCut  = 1;
const int   muTrkLayersCut = 5;///phijet-phimuon
const int   muPixelHitsCut = 0;
//"root://eoscms//eos/cms/store/cmst3/group/hintt/mverweij/PP5TeV/data/SingleMuHighPt/crab_FilteredSingleMuHighPt_v3/160425_163333/merge/HiForest_" PP

// "root://eoscms//eos/cms/store/cmst3/group/hintt/mverweij/PbPb5TeV/data/HIEWQExo/crab_FilteredSingleMuHighPt_v2/160421_135925/mergePartial/HiForest_0.root") PBPB

// "root://eoscms//eos/cms/store/cmst3/group/hintt/psilva/PP5TeV/MC/Pythia8_bJet30_pp502_TuneCUETP8M1/crab_MCpp502_bJet_pthat30_pythia8/160824_090828/0000

// MC store

//PP30 Files: 2,3,4,5,6,7
///root://eoscms//eos/cms/store/cmst3/group/hintt/psilva/PP5TeV/MC/Pythia8_bJet30_pp502_TuneCUETP8M1/crab_MCpp502_bJet_pthat30_pythia8/160824_090828/0000/


//PP50 Files: 2,3,4,5,6,7,8,9,10,11
///root://eoscms//eos/cms/store/cmst3/group/hintt/psilva/PP5TeV/MC/Pythia8_bJet50_pp502_TuneCUETP8M1/crab_MCpp502_bJet_pthat50_pythia8/160824_142558/0000/HiForestAOD_


//PP80 Files: 1,2,3,4,5,6,7,8,9,10
///root://eoscms//eos/cms/store/cmst3/group/hintt/psilva/PP5TeV/MC/Pythia8_bJet80_pp502_TuneCUETP8M1/crab_MCpp502_bJet_pthat80_pythia8/160824_142546/0000/

//PP100 Files: 1,2,3,4,5,6,7,8,9
///root://eoscms//eos/cms/store/cmst3/group/hintt/psilva/PP5TeV/MC/Pythia8_bJet100_pp502_TuneCUETP8M1/crab_MCpp502_bJet_pthat100_pythia8/160824_142629/0000/


//PP120 Files: 1,2,3,4,5,6,7,8,9,10,11
///root://eoscms//eos/cms/store/cmst3/group/hintt/psilva/PP5TeV/MC/Pythia8_bJet120_pp502_TuneCUETP8M1/crab_MCpp502_bJet_pthat120_pythia8/160823_153738/0000/

//PB120 Files: Random
//root://eoscms//eos/cms/store/cmst3/group/hintt/psilva/PbPb5TeV/MC/Pythia6_bJet120_pp502_Hydjet_MB/crab_MCPbPb_bJet_pthat120_pythia6/160824_143451/0000/


double calcLeptonIsolation(float lepPt, float lepEta, float lepPhi, std::vector<float> *pfPt, std::vector<float> *pfEta, std::vector<float> *pfPhi);

void makeMuJetsSkim(const string PNGname,
		    float Rlimit,
		    bool PbPb,	
		    int number_of_entries=300,
		    int start_entry=0,
		    bool createpng=false, 
		    bool isMC = false, 
		    const std::string outFileName = "test.root", 
		    TString inFileList ="root://eoscms//eos/cms/store/cmst3/group/hintt/mverweij/PbPb5TeV/data/HIEWQExo/crab_FilteredSingleMuHighPt_v2/160421_135925/mergePartial/HiForest_0.root")
{

  if(PbPb)jetPtCut  = 50.; 


  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  skimTree_p = new TTree("skimTree", "skimTree");
  BookTree();
  skimTree_p->SetAutoSave(1000000);
  skimTree_p->SetDirectory(outFile_p);
  std::vector<std::string>* inFileNames_p = new std::vector<std::string>;
  TChain *lepTree_p = new TChain("ggHiNtuplizer/EventTree");
  TChain *jetTree_p = PbPb ?  new TChain("akCs2PFJetAnalyzer/t"): /*PB-PB*/ new TChain("ak4PFJetAnalyzer/t"); /*P-P*/
  TChain *hiTree_p = new TChain("hiEvtAnalyzer/HiTree");
  TChain *hltTree_p = new TChain("hltanalysis/HltTree");
  TChain *pfTree_p = PbPb ? new TChain("pfcandAnalyzerCS/pfTree"):/* Pb-Pb*/ new TChain("pfcandAnalyzer/pfTree"); // P-P
  TChain *skimAnaTree_p = new TChain("skimanalysis/HltTree");
 
/* int num[]={1,115,124,127,132,152,154,155,156,160,163,164,167,168,17,170,171,173,46, 5, 53, 63,8,94};
 for(int i=0;i<number_of_files;++i){
 inFileNames_p->push_back(inFileName+to_string(num[i])+".root");
 }*/

  TObjArray *tx = inFileList.Tokenize(",");
  for (Int_t i = 0; i < tx->GetEntries(); i++) 
    inFileNames_p->push_back(tx->At(i)->GetName());
     
  const int nFiles = (int)inFileNames_p->size();
  for(int fileIter = 0; fileIter < nFiles; fileIter++){

  std::cout << "On file: " << fileIter+1 << "/" << nFiles << "  " << inFileNames_p->at(fileIter).c_str()  << std::endl;
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
  float         jteta[maxJets];   
  float         jtphi[maxJets];
  float         jtm[maxJets];   
  float         discr_csvV1[maxJets]; 
  float         svtxm[maxJets];  
  float         svtxpt[maxJets];  
  float         discr_prob[maxJets]; 
  float         discr_csvV2[maxJets];
  float         discr_tcHighEff[maxJets];
  float         discr_tcHighPur[maxJets];
  int           refparton_flavorForB[maxJets];


  std::vector<int>     *genpdg = 0;
  std::vector<float>   *genpt = 0;
  std::vector<float>   *geneta = 0;
  std::vector<float>   *genphi = 0;
  std::vector<int>     *genchg = 0;
  std::vector<int>     *pfId = 0;
  std::vector<float>   *pfPt = 0;
  std::vector<float>   *pfEta = 0;
  std::vector<float>   *pfPhi = 0;
    
  int trig = 1;

  //event selections
  int phfCoincFilter = 1;
  int HBHENoiseFilterResult = 1;
  int pprimaryVertexFilter = 1;
  int pcollisionEventSelection = 1;
 
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
  jetTree_p->SetBranchStatus("discr_csvV2", 1);
  jetTree_p->SetBranchStatus("discr_tcHighEff", 1);
  jetTree_p->SetBranchStatus("discr_tcHighPur", 1);
  jetTree_p->SetBranchStatus("refparton_flavorForB", 1);
  
  jetTree_p->SetBranchStatus("svtxm",1);
  jetTree_p->SetBranchStatus("svtxpt",1);
  jetTree_p->SetBranchStatus("discr_prob",1);
        
  jetTree_p->SetBranchAddress("nref", &nref);
  jetTree_p->SetBranchAddress("jtpt", jtpt);
  jetTree_p->SetBranchAddress("jtphi", jtphi);

  jetTree_p->SetBranchAddress("jteta", jteta);
  jetTree_p->SetBranchAddress("jtm", jtm);
  jetTree_p->SetBranchAddress("svtxm", svtxm);
  jetTree_p->SetBranchAddress("svtxpt", svtxpt);
  jetTree_p->SetBranchAddress("discr_prob", discr_prob);
  jetTree_p->SetBranchAddress("discr_csvV1", discr_csvV1);
  jetTree_p->SetBranchAddress("discr_csvV2", discr_csvV2);
  jetTree_p->SetBranchAddress("discr_tcHighEff", discr_tcHighEff);
  jetTree_p->SetBranchAddress("discr_tcHighPur", discr_tcHighPur);
  jetTree_p->SetBranchAddress("refparton_flavorForB", refparton_flavorForB);
 
  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("run", 1);
  hiTree_p->SetBranchStatus("evt", 1);
  hiTree_p->SetBranchStatus("lumi", 1);
  hiTree_p->SetBranchStatus("hiBin", 1);
  hiTree_p->SetBranchStatus("vz", 1);
  hiTree_p->SetBranchStatus("weight", 1);
  if(isMC && !PbPb)
  hiTree_p->SetBranchStatus("pthat", 1);

  hiTree_p->SetBranchAddress("run", &run_);
  hiTree_p->SetBranchAddress("evt", &evt_);
  hiTree_p->SetBranchAddress("lumi", &lumi_);
  hiTree_p->SetBranchAddress("hiBin", &hiBin_);
  hiTree_p->SetBranchAddress("vz", &vz_);
  if(isMC && !PbPb)
  hiTree_p->SetBranchAddress("pthat", &Pthat_);
  hiTree_p->SetBranchAddress("weight", &weight_);
  pfTree_p->SetBranchAddress("pfId", &pfId);
  pfTree_p->SetBranchAddress("pfPt", &pfPt);
  pfTree_p->SetBranchAddress("pfEta", &pfEta);
  pfTree_p->SetBranchAddress("pfPhi", &pfPhi);
  string triggerName;
  const char* trigname;
  if(PbPb==false){
  triggerName  = isMC ? "HLT_HIL2Mu15ForPPRef_v1" : "HLT_HIL2Mu15_v1"; //v2 para PB-PB
  trigname=triggerName.c_str();
  }
  else {triggerName = "HLT_HIL2Mu15_v2";trigname=triggerName.c_str(); }//v2 para PB-PB;

  hltTree_p->SetBranchStatus(trigname,1);
  hltTree_p->SetBranchAddress(trigname,&trig);

  skimAnaTree_p->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult);

  if(PbPb){
  skimAnaTree_p->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
  skimAnaTree_p->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
  skimAnaTree_p->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  }

  if(number_of_entries<0) number_of_entries=(int)lepTree_p->GetEntries();
  int entryDiv = ((int)(number_of_entries/20));

/*Editted by Luuk*/
  //Used to count the percentage of decrease after eacht cut. Are printed at end of the script
  int t1(0), t2(0), t3(0), t4(0), t5(0), t6(0), t7(0), t8(0);
  int m1(0), m2(0), m3(0), m4(0), m5(0), m6(0), m7(0), m8(0), m9(0), m10(0), m11(0), m12(0), m13(0), m14(0);

  for(int entry = start_entry; entry < start_entry+number_of_entries; entry++){
    if(entry%entryDiv == 0 && number_of_entries >= 10000) std::cout << "Entry # " << entry << "/" << number_of_entries << std::endl;
    if(entry%10000==0)cout<<entry<<endl;
    hiTree_p->GetEntry(entry);
    lepTree_p->GetEntry(entry);
    jetTree_p->GetEntry(entry);
    hltTree_p->GetEntry(entry);
    pfTree_p->GetEntry(entry);
    skimAnaTree_p->GetEntry(entry);
    if(!isMC && !trig) continue;
   
    if(PbPb)
    if(!phfCoincFilter) continue;
  
    if(!HBHENoiseFilterResult) continue;

    if(PbPb)
    if(!pcollisionEventSelection) continue;

    if(PbPb)
    if(!pprimaryVertexFilter) continue;
   
    if(TMath::Abs(vz_) > 15) continue;


    float tempMuPt_[nLep];
    float tempMuPhi_[nLep];
    float tempMuEta_[nLep];
    int tempMuChg_[nLep];
    float tempMuIso_[nLep];     
    float tempMuInnerDz_[nLep];


    for(int lepIter = 0; lepIter < nLep; lepIter++){
      lepPt_[lepIter] = -999;
      lepPhi_[lepIter] = -999;
      lepEta_[lepIter] = -999;
      lepChg_[lepIter] = -999;
      lepID_[lepIter] = -999;
      lepIso_[lepIter] = -999;
      lepInnerDz_[lepIter] = -999;
    }

   
       
    for(int lepIter = 0; lepIter < nLep; lepIter++){
      tempMuPt_[lepIter] = -999;
      tempMuPhi_[lepIter] = -999;
      tempMuEta_[lepIter] = -999;
      tempMuChg_[lepIter] = -999;
      tempMuIso_[lepIter] = -999;
      tempMuInnerDz_[lepIter] = -999;
    }
       
    for(int ij = 0; ij<nMaxJets; ++ij) {
      jtPt_[ij] = -999.; 

     // IsoPt_[ij] = -999.;
     // NIsoPt_[ij] = -999.;

      jtEta_[ij] = -999.;
      jtPhi_[ij] = -999.;
      jtM_[ij] = -999.;
      svtxM_[ij] = -999.;
      svtxPt_[ij] = -999.;
      discr_Prob_[ij] = -999.;
      discr_ProbforB_[ij] = -999.;
      discr_ProbforC_[ij] = -999.;
      discr_ProbforBack_[ij] = -999.;
      deltaPhi_[ij] = -999.;
      deltaEta_[ij] = -999.;
      ptMuJet_[ij] = -999.;
      discr_csvV1_[ij] = -999.;
      discr_csvV2_[ij] = -999.;
      discr_tcHighEff_[ij] = -999.;
      discr_tcHighPur_[ij] = -999.;
      refparton_flavorForB_[ij] = -999;

      ptJP_[ij] = -999.;
      R_[ij] = -999.;
      deltaphiJJ_[ij] = -999.;
      //CSVJP_[ij] = -999.;
      svtxptJP_[ij] = -999.;
      discrProbJP_[ij] = -999.;
      discrcsvV2JP_[ij] = -999.;

    }
    for(int ijgen = 0; ijgen < nMaxGen; ++ijgen){
      genPdg_[ijgen] = -999;
      genPt_[ijgen] = -999;
      genEta_[ijgen] = -999;
      genPhi_[ijgen] = -999;
      genChg_[ijgen] = -999;
    }

    //Find two leading muons
    for(unsigned int muIter = 0; muIter < fForestMu.muPt->size(); ++muIter) {
      if(TMath::Abs(fForestMu.muEta->at(muIter)) > muEtaCut) continue;
      m1++;
      if(fForestMu.muPt->at(muIter) < muPtCut) continue;
      m2++;
      if(fForestMu.muChi2NDF->at(muIter) >= muChi2NDFCut) continue;
      m3++;
      if(TMath::Abs(fForestMu.muInnerD0->at(muIter)) > muInnerD0Cut) continue;
      m4++;
      if(TMath::Abs(fForestMu.muInnerDz->at(muIter)) > muInnerDzCut) continue;
      m5++;
      if(fForestMu.muMuonHits->at(muIter) <= muMuonHitsCut) continue;
      m6++;
      if(fForestMu.muStations->at(muIter) <= muStationsCut) continue;
      m7++;
      if(fForestMu.muTrkLayers->at(muIter) <= muTrkLayersCut) continue;
      m8++;
      if(fForestMu.muPixelHits->at(muIter) <= muPixelHitsCut) continue;
      m9++;
      if(fForestMu.muPt->at(muIter) > lepPt_[0]){
        tempMuPt_[1]  = tempMuPt_[0];
        tempMuPhi_[1] = tempMuPhi_[0];
        tempMuEta_[1] = tempMuEta_[0];
        tempMuChg_[1] = tempMuChg_[0];
        tempMuIso_[1] = tempMuIso_[0]; 
        tempMuInnerDz_[1] = tempMuInnerDz_[0];
 
        tempMuPt_[0]  = fForestMu.muPt->at(muIter);
        tempMuPhi_[0] = fForestMu.muPhi->at(muIter);
        tempMuEta_[0] = fForestMu.muEta->at(muIter);
        tempMuChg_[0] = fForestMu.muCharge->at(muIter);
        tempMuIso_[0] = calcLeptonIsolation(fForestMu.muPt->at(muIter),fForestMu.muEta->at(muIter),fForestMu.muPhi->at(muIter),pfPt,pfEta,pfPhi); //iso;
        tempMuInnerDz_[0] = fForestMu.muInnerDz->at(muIter);
      }
      else if(fForestMu.muPt->at(muIter) > tempMuPt_[1]){
        tempMuPt_[1]  = fForestMu.muPt->at(muIter);
        tempMuPhi_[1] = fForestMu.muPhi->at(muIter);
        tempMuEta_[1] = fForestMu.muEta->at(muIter);
        tempMuChg_[1] = fForestMu.muCharge->at(muIter);
        tempMuIso_[1] = calcLeptonIsolation(fForestMu.muPt->at(muIter),fForestMu.muEta->at(muIter),fForestMu.muPhi->at(muIter),pfPt,pfEta,pfPhi);//iso;
        tempMuInnerDz_[1] = fForestMu.muInnerDz->at(muIter);
      }
    }
    //store muons in out tree
    int lepIter = 0;
    for(int muIter = 0; muIter < 2; muIter++){
      if(tempMuPt_[muIter]<0.) continue;
      m11++;
      m12++;
      lepPt_[lepIter] = tempMuPt_[muIter];
      lepPhi_[lepIter] = tempMuPhi_[muIter];
      lepEta_[lepIter] = tempMuEta_[muIter];
      lepChg_[lepIter] = tempMuChg_[muIter];
      lepID_[lepIter] = muID;
      lepIso_[lepIter] = tempMuIso_[muIter];
      lepInnerDz_[lepIter] = tempMuInnerDz_[muIter];
      ++lepIter; 
    }
   // cout<<lepIter <<endl;
    if(lepIter<1) continue;
    t7++;
    m13 += lepIter;
    nLep_ = lepIter;
    int njets = 1;

    int tempindex=0; 

    for(int jetIter = 0; jetIter < nref; jetIter++){
      
      if(jtpt[jetIter]<jetPtCut) continue;
      if(fabs(jteta[jetIter])>jetEtaCut) continue;
      
      float dphilj( TMath::Abs(TVector2::Phi_mpi_pi(jtphi[jetIter]-lepPhi_[0])) );
      float detalj ( TMath::Abs(jteta[jetIter]-lepEta_[0]) );
      float deltaR( TMath::Sqrt(dphilj*dphilj+detalj*detalj) );
      

      if(deltaR < Rlimit){ // Valor normal=0.2
	tempindex=0;
	}
      else if(dphilj>3.*M_PI/4.){
	++njets;
	tempindex=njets-1;	
	} 
      else continue;

      jtPt_[tempindex]  = jtpt[jetIter];
      jtEta_[tempindex] = jteta[jetIter];
      jtPhi_[tempindex] = jtphi[jetIter];
      jtM_[tempindex]   = jtm[jetIter]; 
      svtxM_[tempindex]   = svtxm[jetIter];
      discr_Prob_[tempindex]   = discr_prob[jetIter];
      svtxPt_[tempindex]   = svtxpt[jetIter];
     
      deltaPhi_[tempindex] = dphilj; 
      deltaEta_[tempindex] = detalj;
      ptMuJet_[tempindex] = lepPt_[0]/jtPt_[tempindex];
      discr_csvV1_[tempindex] = discr_csvV1[jetIter];
      discr_csvV2_[tempindex] = discr_csvV2[jetIter];
      discr_tcHighEff_[tempindex] = discr_tcHighEff[jetIter];
      discr_tcHighPur_[tempindex] = discr_tcHighPur[jetIter];
      refparton_flavorForB_[tempindex] = refparton_flavorForB[jetIter];
    }
    
   int index=0;
   float minR=0;
   int Bcount=0,Ccount=0,Backcount=0;
   
   nJt_=njets;
   if(nJt_<2) continue; 

   EventCount_=0;
   
  if(!PbPb){
   for(int i=1;i<njets;++i){

    if(jtPt_[0]==-999){
      R_[i-1] = jtPt_[i]/lepPt_[0];
      deltaphiJJ_[i-1]=TMath::Abs(TVector2::Phi_mpi_pi(jtPhi_[i]-jtPhi_[0])); 
     
    }

    else {
    
    R_[i-1] = jtPt_[i]/jtPt_[0];
    deltaphiJJ_[i-1]=TMath::Abs(TVector2::Phi_mpi_pi(jtPhi_[i]-jtPhi_[0])); 

    }
 
    if(i!=1 && TMath::Abs(R_[i-1]-1)<TMath::Abs(R_[i-2]-1)){
	   minR=TMath::Abs(R_[i-1]-1);
	   index=i; // posicao no vetor R
	} 
    else if(i==1){minR=TMath::Abs(R_[0]-1);index=1;}
 }
}


   else {
   index=1;
   for(int i=1;i<njets;++i){
   if(jtPt_[i+1]> jtPt_[i])index=i; 
       }
   }

   int tagreference= abs(refparton_flavorForB_[index]); 
  
   if(tagreference==5){ 
     discr_ProbforB_[Bcount]=discr_Prob_[index];
     Bcount++;
   }

   else if(tagreference==4){ 
     discr_ProbforC_[Ccount]=discr_Prob_[index];
     Ccount++;
   }

   else if(tagreference==21 || tagreference<4){ 
     discr_ProbforBack_[Backcount]=discr_Prob_[index];
     Backcount++;
   }
  
   Index_=index;
  

    float cut[5]; cut[0]=0.58;cut[1]=0.45; cut[2]=0.3;cut[3]=0.24;cut[4]=0.18;
    float div=1.;

    if(lepIso_[0]<cut[0]/div && hiBin_>0 && hiBin_< 10){   
    IsoPt_=jtPt_[index];
    csvV2IsoPt_=discr_csvV2_[index];
    }

    else if(lepIso_[0]<cut[1]/div && hiBin_>10 && hiBin_< 30){   
    IsoPt_=jtPt_[index];
    csvV2IsoPt_=discr_csvV2_[index];
    }

    else if(lepIso_[0]<cut[2]/div && hiBin_>30 && hiBin_< 50){   
    IsoPt_=jtPt_[index];
    csvV2IsoPt_=discr_csvV2_[index];
    }

    else if(lepIso_[0]<cut[3]/div && hiBin_>50 && hiBin_<70){   
    IsoPt_=jtPt_[index];
    csvV2IsoPt_=discr_csvV2_[index];
    }

    else if(lepIso_[0]<cut[4]/div && hiBin_>70 && hiBin_<100){   
    IsoPt_=jtPt_[index];
    csvV2IsoPt_=discr_csvV2_[index];
    }

    else{
    NIsoPt_=jtPt_[index];
    csvV2NIsoPt_=discr_csvV2_[index];
    }
  
    //cout<<"MinR: "<<R_[index]<<" "<< "Index: "<<index<<" "<<" JtPt(index): "<<jtPt_[index]<<" "<<"JtPt(0): "<<jtPt_[0]<<endl; 

    skimTree_p->Fill();
    
  }//entries loop

  outFile_p->cd();
  TNamed pathStr1("pathStr1", outFileName.c_str());
  pathStr1.Write("", TObject::kOverwrite);
  TNamed pathStr2("pathStr2", inFileList.Data());
  pathStr2.Write("", TObject::kOverwrite);
  skimTree_p->Write("", TObject::kOverwrite);
  outFile_p->Close();
  delete outFile_p;

  //  drawMuJetsSkimControlPlots(PNGname,outFileName,createpng);
}

double calcLeptonIsolation(float lepPt, float lepEta, float lepPhi,
                           std::vector<float> *pfPt, std::vector<float> *pfEta, std::vector<float> *pfPhi) {
  //calculate lepton isolation from pf candidates.
  //Isolation cone R=0.3
  //excluding pf candidates at a distance less than 0.03 from lepton
 
  double conePt = 0.;
  for(unsigned int i = 0; i<pfPt->size(); ++i) {
 
    float deta=TMath::Abs(lepEta-pfEta->at(i));
    if(deta>0.3) continue;
    float dphi=TMath::Abs(TVector2::Phi_mpi_pi(lepPhi-pfPhi->at(i)));
    if(dphi>0.3) continue;
    float dr(sqrt(dphi*dphi+deta*deta));
    if(dr>0.3 || dr<0.03) continue;
    conePt+=pfPt->at(i);
  }
  double relIso = conePt;
  if(lepPt>0.) relIso = conePt/lepPt;
 
  return relIso;
}


//FIXME: R; lepIso; weight; 

//FIXME: PbPb: Ptjet cut =50; escolher jato com mais pt em vez de o mais proximo; mas no hemisf. oposto 

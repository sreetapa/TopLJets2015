#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TSystem.h"

#include <string>
#include <vector>

#include "HeavyIonsAnalysis/topskim/include/ForestHiTree.h"
#include "HeavyIonsAnalysis/topskim/include/ForestElectrons.h"
#include "HeavyIonsAnalysis/topskim/include/ForestMuons.h"
#include "HeavyIonsAnalysis/topskim/include/ForestPFCands.h"
#include "HeavyIonsAnalysis/topskim/include/ForestJets.h"
#include "HeavyIonsAnalysis/topskim/include/LumiRun.h"
#include "HeavyIonsAnalysis/topskim/include/HistTool.h"
#include "HeavyIonsAnalysis/topskim/include/PFAnalysis.h"
#include "HeavyIonsAnalysis/topskim/include/LeptonSummary.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "../scripts/functions.cc"


const bool isDebug = true;
const float lepPtCut  = 20.;
const float lepEtaCut = 2.1;
//see https://indico.cern.ch/event/803679/contributions/3342407/attachments/1808912/2953435/egm-minipog-190308.pdf
const float eeScaleShift = 6.8182E-2/5.9097E-2;
const int firstEEScaleShiftRun = 327402; 
const float barrelEndcapEta[2]={1.4442,1.5660};
const float csvWP = 0.8838;

using namespace std;
using namespace fastjet;
using namespace TMVA;

//dispersion of rapidities of the final state
std::vector<float> getRapidityMoments(std::vector<TLorentzVector> & coll){
  std::vector<float> mom(3,0.);
  for(size_t i=0; i<coll.size(); i++) {
    TLorentzVector pi(coll[i]);
    mom[0]+=pi.Rapidity();
    mom[1]+=pow(pi.Rapidity(),2);
    for(size_t j=i; j<coll.size(); j++) {
      TLorentzVector pj(coll[j]);
      float dy=fabs(pj.Rapidity()-pi.Rapidity());
      mom[2]=max(dy,mom[2]);
    }
  }
  
  mom[0]=mom[0]/float(mom.size());
  mom[1]=sqrt(mom[1]/float(mom.size())-pow(mom[0],2));
  return mom;
}



// index, ntks in svtx, m svtx, csv
typedef std::tuple<int,int,float,float> BtagInfo_t;
static bool orderByBtagInfo(const BtagInfo_t &a, const BtagInfo_t &b)
{
  //int ntks_a(std::get<1>(a)), ntks_b(std::get<1>(b));
  //if(ntks_a>ntks_b) return true;

  float csv_a(std::get<3>(a)), csv_b(std::get<3>(b));
  if(csv_a>csv_b) return true;
  return false;
}


static bool orderByPt(const LeptonSummary &a, const LeptonSummary &b)
{
  float pt_a(a.p4.Pt()), pt_b(b.p4.Pt());
  if(pt_a>pt_b) return true;
  return false;
}



//
int main(int argc, char* argv[])
{
  bool blind(true);
  TString inURL,outURL;
  bool isMC(false),isPP(false);
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--in")!=string::npos && i+1<argc)       { inURL=TString(argv[i+1]); i++;}
    else if(arg.find("--out")!=string::npos && i+1<argc) { outURL=TString(argv[i+1]); i++;}
    else if(arg.find("--mc")!=string::npos)              { isMC=true;  }
    else if(arg.find("--pp")!=string::npos)              { isPP=true;  }
  }

  bool isSingleMuPD( !isMC && inURL.Contains("SkimMuons"));
  bool isSingleElePD( !isMC && inURL.Contains("SkimElectrons"));
  LumiRun lumiTool;

  if(isPP)
    cout << "Treating as a pp collision file" << endl;
  if(isMC)
    cout << "Treating as a MC file" << endl;

  //book some histograms
  HistTool ht;

  if(!isMC) ht.addHist("ratevsrun",lumiTool.getLumiMonitor());

  //generic histograms
  for(int i=0; i<2; i++) {
    TString pf(Form("l%d",i+1));
    ht.addHist(pf+"pt",             new TH1F(pf+"pt",            ";Lepton transverse momentum [GeV];Events",20,20,200));
    ht.addHist(pf+"eta",            new TH1F(pf+"eta",           ";Lepton pseudo-rapidity;Events",20,0,2.5));

    for(int j=0; j<3; j++) {
      TString comp("ch");
      if(j==1) comp="pho";
      if(j==2) comp="nh";
      ht.addHist(pf+comp+"iso",      new TH1F(pf+comp+"iso",      ";PF "+comp+" isolation;Leptons",50,0,250));
      ht.addHist(pf+comp+"reliso",   new TH1F(pf+comp+"reliso",   ";Relative PF "+comp+" isolation;Leptons",50,0,2.0));
      ht.addHist(pf+comp+"isovscen", new TH2F(pf+comp+"isovscen", ";Centrality bin;PF "+comp+" isolation [GeV];Leptons",10,0,100,50,0,100));
      ht.addHist(pf+comp+"isovsrho", new TH2F(pf+comp+"isovsrho", ";#rho_{"+comp+"};PF "+comp+" isolation [GeV];Leptons",10,0,100,20,0,100));
    }    
  }

  //electron specific
  ht.addHist("esihih",      new TH1F("esihih",      ";#sigma(i#etai#eta);Electrons",       50,0,0.06));
  ht.addHist("edetavtx",    new TH1F("edetavtx",    ";#Delta#eta(vtx);Electrons",          50,0,0.015));
  ht.addHist("edphivtx",    new TH1F("edphivtx",    ";#Delta#phi(vtx) [rad];Electrons",    50,0,0.015));
  ht.addHist("ehoe",        new TH1F("ehoe"    ,    ";h/e;Electrons",                      50,0,0.25));
  ht.addHist("eempinv",     new TH1F("eempinv",     ";|1/E-1/p| [1/GeV];Electrons",        50,0,0.05));
  ht.addHist("ed0",         new TH1F("ed0",         ";d_{0} [cm];Electrons",               50,0,0.05));
  ht.addHist("edz",         new TH1F("edz",         ";d_{z} [cm];Electrons",               50,0,0.05));

  //muon specific
  ht.addHist("mmusta",     new TH1F("mmusta",      ";Muon stations;Muons",            15,0,15));   
  ht.addHist("mtrklay",    new TH1F("mtrklay",     ";Tracker layers;Muons",           25,0,25));
  ht.addHist("mchi2ndf",   new TH1F("mchi2ndf",    ";#chi^2/ndf;Muons",               50,0,15));
  ht.addHist("mmuhits",    new TH1F("mmuhits",     ";Muon hits;Muons",                25,0,25));
  ht.addHist("mpxhits",    new TH1F("mpxhits",     ";Pixel hits;Muons",               15,0,15));
  ht.addHist("md0",        new TH1F("md0",         ";d_{0} [cm];Muons",               50,0,0.5));
  ht.addHist("mdz",        new TH1F("mdz",         ";d_{z} [cm];Muons",               50,0,1.0));
 
  ht.addHist("mll",      new TH1F("mll",      ";Dilepton invariant mass [GeV];Events",40,0,200));
  ht.addHist("ptll",     new TH1F("ptll",     ";Dilepton transverse momentum [GeV];Events",25,0,200));
  ht.addHist("ptsum",    new TH1F("ptsum",    ";p_{T}(l)+p_{T}(l') [GeV];Events",25,0,200));
  ht.addHist("acopl",    new TH1F("acopl" ,   ";1-#Delta#phi(l,l')/#pi;Events",20,0,1.0));
  ht.addHist("detall",   new TH1F("detall",   ";#Delta#eta(l,l');Events",20,0,4));
  ht.addHist("drll",     new TH1F("drll"  ,   ";#DeltaR(l,l');Events",20,0,2*TMath::Pi()));
  ht.addHist("chrho",    new TH1F("chrho",    ";#rho_{ch};Events",25,0,50));
  ht.addHist("phorho",   new TH1F("phorho",   ";#rho_{#gama};Events",25,0,50));
  ht.addHist("nhrho",    new TH1F("nhrho",    ";#rho_{nh};Events",25,0,50));

  ht.addHist("pfrapavg",      new TH1F("pfrapavg",     ";Average rapidity;Events",25,0,2.5));
  ht.addHist("pfraprms",      new TH1F("pfraprms",     ";#sigma(rapidity);Events",50,0,2.5));
  ht.addHist("pfrapmaxspan",  new TH1F("pfrapmaxspan", ";Maximum rapidity span;Events",25,0,5));
  ht.addHist("pfht",          new TH1F("pfht",         ";H_{T} [GeV];Events",25,0,500));
  ht.addHist("pfmht",         new TH1F("pfmht",        ";Missing H_{T} [GeV];Events",25,0,200));
  ht.addHist("npfjets",    new TH1F("npfjets",    ";Jet multiplicity;Events",8,0,8));
  ht.addHist("npfbjets",   new TH1F("npfbjets",   ";b-jet multiplicity;Events",5,0,5));    
  ht.addHist("npfsvtx",    new TH1F("npfsvtx",    ";Secondary vertex multiplicity;Events",5,0,5));
  for(size_t j=1; j<=2; j++){
    TString ppf(j==1 ? "1" : "2");
    ht.addHist("pf"+ppf+"jbalance",    new TH1F("pf"+ppf+"jbalance", ";R = p_{T}(j)/p_{T}(ll);Events",50,0,3));
    ht.addHist("pf"+ppf+"jpt",         new TH1F("pf"+ppf+"jpt",      ";Jet transverse momentum [GeV];Events",30,00,300));
    ht.addHist("pf"+ppf+"jeta",        new TH1F("pf"+ppf+"jeta",     ";Jet pseudo-rapidity;Events",20,0,2.5));
    ht.addHist("pf"+ppf+"jetavsphi",   new TH2F("pf"+ppf+"jetavsphi", ";Jet pseudo-rapidity;Jet azimuthal angle [rad];Events",100,-2.5,2.5,100,-TMath::Pi(),TMath::Pi()));
    ht.addHist("pf"+ppf+"jsvtxm",      new TH1F("pf"+ppf+"jsvtxm",   ";Secondary vertex mass;Events",25,0,6));
    ht.addHist("pf"+ppf+"jsvtxntk",    new TH1F("pf"+ppf+"jsvtxntk", ";Secondary vertex track multiplicity;Events",5,0,5));
    ht.addHist("pf"+ppf+"jcsv",        new TH1F("pf"+ppf+"jcsv",     ";CSVv2;Events",25,0,1));
  }

  //configure leptons
  TChain *lepTree_p     = new TChain(isPP ? "ggHiNtuplizer/EventTree" : "ggHiNtuplizerGED/EventTree");
  lepTree_p->Add(inURL);
  ForestMuons fForestMu(lepTree_p);  
  ForestElectrons fForestEle(lepTree_p);

  //configure PF cands
  TChain *pfCandTree_p  = new TChain("pfcandAnalyzer/pfTree");
  pfCandTree_p->Add(inURL);
  ForestPFCands fForestPF(pfCandTree_p);

  //configure jets
  TChain *jetTree_p     = new TChain(isPP ? "ak4PFJetAnalyzer/t" : "akPu4PFJetAnalyzer/t");
  jetTree_p->Add(inURL);
  ForestJets fForestJets(jetTree_p);

  //global variables
  TChain *hiTree_p      = new TChain("hiEvtAnalyzer/HiTree");
  hiTree_p->Add(inURL);
  HiTree fForestTree(hiTree_p);

  //trigger
  TChain *hltTree_p     = new TChain("hltanalysis/HltTree");
  hltTree_p->Add(inURL);
  int etrig(0),mtrig(0);
  if(isPP){
    TString muTrig("HLT_HIL3Mu20_v1");
    if( !hltTree_p->FindBranch(muTrig) ) muTrig="HLT_HIL3Mu15ForPPRef_v1";
    hltTree_p->SetBranchStatus(muTrig,1);
    hltTree_p->SetBranchAddress(muTrig,&mtrig);
    TString eTrig("HLT_HIEle20_WPLoose_Gsf_v1");
    if( !hltTree_p->FindBranch(eTrig) ) eTrig="HLT_HISinglePhoton20_Eta3p1ForPPRef_v1";
    hltTree_p->SetBranchStatus(eTrig,1);
    hltTree_p->SetBranchAddress(eTrig,&etrig);
    cout << "Using " << muTrig << " " << eTrig << " as MC triggers" << endl;
  }else{
    hltTree_p->SetBranchStatus("HLT_HIL3Mu15_v1",1);
    hltTree_p->SetBranchAddress("HLT_HIL3Mu15_v1",&mtrig);
    hltTree_p->SetBranchStatus("HLT_HIEle20Gsf_v1",1);
    hltTree_p->SetBranchAddress("HLT_HIEle20Gsf_v1",&etrig);    
  }

  TChain *rhoTree_p = new TChain("hiFJRhoAnalyzer/t");
  rhoTree_p->Add(inURL);
  std::vector<Double_t> *t_rho=0,*t_rhom=0;
  if(rhoTree_p){
    rhoTree_p->SetBranchAddress("rho", &t_rho);
    rhoTree_p->SetBranchAddress("rhom", &t_rhom);
  }else{
    std::cout << "[WARN] Can't find rho tree hiFJRhoAnalyzer/t" << std::endl;
  }

  
  // =============================================================
  // marc here make the output tree

  TTree * outTree = new TTree("tree", "tree with 2lepton selection and combined collections");

  // event and trigger variables
  Int_t  t_run, t_lumi, t_etrig, t_mtrig;
  Long_t t_event;
  Float_t t_weight, t_cenbin;
  Float_t t_chrho[3], t_phorho[3], t_nhrho[3];
  outTree->Branch("run"   , &t_run  , "run/I");
  outTree->Branch("lumi"  , &t_lumi , "lumi/I");
  outTree->Branch("event" , &t_event, "event/L");

  outTree->Branch("weight", &t_weight, "weight/F");

  // centrality and different flavors of rho
  outTree->Branch("cenbin", &t_cenbin, "cenbin/F");
  
  outTree->Branch("rho",    &t_rho);
  outTree->Branch("rhom",   &t_rhom);
  outTree->Branch("chrho" , t_chrho , "chrho[3]/F");
  outTree->Branch("nhrho" , t_nhrho , "nhrho[3]/F");
  outTree->Branch("phorho", t_phorho, "phorho[3]/F");

  outTree->Branch("etrig" , &t_etrig , "etrig/I");
  outTree->Branch("mtrig" , &t_mtrig , "mtrig/I");

  // variables per lepton, including iso
  Int_t t_nlep;
  std::vector<Float_t> t_lep_pt, t_lep_eta, t_lep_phi, t_lep_phiso, t_lep_chiso, t_lep_nhiso, t_lep_rho, t_lep_chrho, t_lep_nhrho, t_lep_phrho;
  std::vector<Int_t  > t_lep_pdgId, t_lep_charge;
  outTree->Branch("nlep"      , &t_nlep      , "nlep/I"            );
  outTree->Branch("lep_pt"    , &t_lep_pt);
  outTree->Branch("lep_eta"   , &t_lep_eta);
  outTree->Branch("lep_phi"   , &t_lep_phi);
  outTree->Branch("lep_phiso" , &t_lep_phiso);
  outTree->Branch("lep_chiso" , &t_lep_chiso);
  outTree->Branch("lep_nhiso" , &t_lep_nhiso);
  outTree->Branch("lep_rho"   , &t_lep_rho);
  outTree->Branch("lep_phrho" , &t_lep_phrho);
  outTree->Branch("lep_chrho" , &t_lep_chrho);
  outTree->Branch("lep_nhrho" , &t_lep_nhrho);
  outTree->Branch("lep_pdgId" , &t_lep_pdgId);
  outTree->Branch("lep_charge", &t_lep_charge);

  // variables from dilepton system
  Float_t t_llpt, t_lleta, t_llphi, t_llm, t_dphi, t_deta, t_sumeta;
  outTree->Branch("llpt"   , &t_llpt   , "llpt/F");
  outTree->Branch("lleta"  , &t_lleta  , "lleta/F");
  outTree->Branch("llphi"  , &t_llphi  , "llphi/F");
  outTree->Branch("llm"    , &t_llm    , "llm/F");
  outTree->Branch("dphi"   , &t_dphi   , "dphi/F");
  outTree->Branch("deta"   , &t_deta   , "deta/F");
  outTree->Branch("sumeta" , &t_sumeta , "sumeta/F");

  // variables per jet (jets ordered by pt)
  Int_t t_njet;
  std::vector<Float_t> t_jet_pt, t_jet_eta, t_jet_phi, t_jet_mass, t_jet_csvv2;
  outTree->Branch("njet"      , &t_njet      , "njet/I"            );
  outTree->Branch("jet_pt"    , &t_jet_pt    , "jet_pt[njet]/F"    );
  outTree->Branch("jet_eta"   , &t_jet_eta   , "jet_eta[njet]/F"   );
  outTree->Branch("jet_phi"   , &t_jet_phi   , "jet_phi[njet]/F"   );
  outTree->Branch("jet_mass"  , &t_jet_mass  , "jet_mass[njet]/F" );
  outTree->Branch("jet_csvv2" , &t_jet_csvv2 , "jet_csvv2[njet]/F");

  // variables per bjet (jets ordered by csvv2)
  Int_t t_nbjet;
  std::vector<Float_t> t_bjet_pt, t_bjet_eta, t_bjet_phi, t_bjet_mass, t_bjet_csvv2;
  outTree->Branch("nbjet"      , &t_nbjet      , "nbjet/I"            );
  outTree->Branch("bjet_pt"    , &t_bjet_pt    );
  outTree->Branch("bjet_eta"   , &t_bjet_eta   );
  outTree->Branch("bjet_phi"   , &t_bjet_phi   );
  outTree->Branch("bjet_mass"  , &t_bjet_mass  );
  outTree->Branch("bjet_csvv2" , &t_bjet_csvv2 );

  // constructed variables like ht and stuff
  Float_t t_ht, t_mht, t_apt, t_dphilll2;
  outTree->Branch("ht"     , &t_ht     , "ht/F");
  outTree->Branch("mht"    , &t_mht    , "mht/F");
  outTree->Branch("apt"    , &t_apt    , "apt/F");
  outTree->Branch("dphilll2"    , &t_dphilll2    , "dphilll2/F");

  Float_t t_bdt, t_bdt_rarity, t_fisher2;
  outTree->Branch("bdt", &t_bdt, "bdt/F");
  outTree->Branch("bdtrarity", &t_bdt_rarity, "bdtrarity/F");
  outTree->Branch("fisher2", &t_fisher2, "fisher2/F");
  // =============================================================
  //
  TMVA::Tools::Instance();
  TMVA::Reader *reader        = new TMVA::Reader( "!Color:!Silent" );
  TMVA::Reader *readerFisher2 = new TMVA::Reader( "!Color:!Silent" );

  // make new variables because i'm too lazy to think right now
  Float_t bdt_l1pt, bdt_apt, bdt_abslleta, bdt_dphilll2, bdt_sumabseta, bdt_flavor;

  // these must have the same name as in the training. and the same ordeeeeeeer
  // copy directly from the script that runs the training:
  //dataloader.AddVariable('lep_pt[0]'  , 'p_{T}^{lep1}'     , 'GeV' , 'F')
  //dataloader.AddVariable('apt'        , 'A_{pt}'           , ''    , 'F')
  //dataloader.AddVariable('llpt'       , 'p_{T}^{ll}'       , 'GeV' , 'F')
  //dataloader.AddVariable('abs(lleta)' , '|#eta^{ll}|'      , ''    , 'F')
  //dataloader.AddVariable('dphi'       , '|#Delta #phi|'    , 'rad' , 'F')
  //dataloader.AddVariable('abs(lep_eta[0])+abs(lep_eta[1])' , '#sum |#eta_{i}|', ''    , 'F')
  //dataloader.AddVariable('abs(lep_pdgId[0]*lep_pdgId[1])'  , 'flavor', ''    , 'F')
  //
  reader->AddVariable("lep_pt[0]"  , &bdt_l1pt    );
  reader->AddVariable("apt"        , &bdt_apt     );
  reader->AddVariable("llpt"       , &t_llpt      );
  reader->AddVariable("abs(lleta)" , &bdt_abslleta);
  reader->AddVariable("dphi"       , &t_dphi      );
  reader->AddVariable("abs(lep_eta[0])+abs(lep_eta[1])", &bdt_sumabseta);
  reader->AddVariable("abs(lep_pdgId[0]*lep_pdgId[1])" , &bdt_flavor);

  // for the fisher just take these two
  //dataloader.AddVariable('llpt'       , 'p_{T}^{ll}'       , 'GeV' , 'F')
  //dataloader.AddVariable('dphi'       , '|#Delta #phi|'    , 'rad' , 'F')
  readerFisher2->AddVariable("llpt", &t_llpt);
  readerFisher2->AddVariable("dphi", &t_dphi);

  TString methodName       ("BDTG method");
  TString methodNameFisher2("Fisher method");
  // hard coded path for now ...
  TString weightFile("/afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_10_3_1/src/HeavyIonsAnalysis/topskim/scripts/trainingV2_sevenVars_includeEMuZ/weights/TMVAClassification_BDTG.weights.xml");
  reader->BookMVA( methodName, weightFile);

  TString weightFileFisher2("/afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_10_3_1/src/HeavyIonsAnalysis/topskim/scripts/trainingV2_Fisher2_includeEMuZ/weights/TMVAClassification_Fisher.weights.xml");
  readerFisher2->BookMVA( methodNameFisher2, weightFileFisher2);

    
  Double_t wgtSum(0);
  int nEntries = (int)lepTree_p->GetEntries();  
  int entryDiv = ((int)(nEntries/20));    
  cout << inURL << " has " << nEntries << " events to process" << endl;
  for(int entry = 0; entry < nEntries; entry++){
    
    if(entry%entryDiv == 0) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
    
    lepTree_p->GetEntry(entry);
    pfCandTree_p->GetEntry(entry);
    jetTree_p->GetEntry(entry);    
    hltTree_p->GetEntry(entry);
    hiTree_p->GetEntry(entry);
    if(rhoTree_p) rhoTree_p->GetEntry(entry);

    float evWgt(1.0);
    if(isMC && fForestTree.ttbar_w->size()) { evWgt=fForestTree.ttbar_w->at(0); }
    wgtSum += evWgt;    
    float plotWgt(evWgt);

    //first of all require a trigger
    int trig=etrig+mtrig;
    if(trig==0) continue;
    if(isSingleMuPD) {
      if(mtrig==0) continue;
      if(etrig!=0) continue;
    }
    if(isSingleElePD) {
      if(etrig==0) continue;
    }

    //FIXME: WHAT SHOULD WE DO HERE ?
    //apply global filters
    //if(!isPP){
    //  if(TMath::Abs(fForestTree.vz) > 15) continue;
    // }

    //build jets from different PF candidate collections   
    SlimmedPFCollection_t pfColl;
    for(size_t ipf=0; ipf<fForestPF.pfId->size(); ipf++) {
      int id(abs(fForestPF.pfId->at(ipf)));
      float mass(0.13957);  //pions
      if(id==4) mass=0.;    //photons
      if(id>=5) mass=0.497; //K0L
      pfColl.push_back( getSlimmedPF( id, fForestPF.pfPt->at(ipf),fForestPF.pfEta->at(ipf),fForestPF.pfPhi->at(ipf),mass) );
    }

    Float_t chrho[3]  = { getRho(pfColl,{1,2,3}),  getRho(pfColl,{1,2,3},-1,1.4442), getRho(pfColl,{1,2,3},1.4442,2.4) };
    Float_t phorho[3] = { getRho(pfColl,{4}),      getRho(pfColl,{4},-1,1.4442),     getRho(pfColl,{4},1.4442,2.4) };
    Float_t nhrho[3]  = { getRho(pfColl,{5,6}),    getRho(pfColl,{5,6},-1,1.4442),   getRho(pfColl,{5,6},1.4442,2.4) };

    //monitor trigger and centrality
    float cenBin=0;
    bool isCentralEvent(false);
    if(!isMC){
      cenBin=0.5*fForestTree.hiBin;
      isCentralEvent=(!isPP && cenBin<30);
      Int_t runBin=lumiTool.getRunBin(fForestTree.run);
      Float_t lumi=lumiTool.getLumi(fForestTree.run);
      if(lumi>0.){
        if(etrig>0) ht.fill("ratevsrun",runBin,1./lumi,"e");
        if(mtrig>0) ht.fill("ratevsrun",runBin,1./lumi,"m");
      }
    }

    //the selected leptons
    std::vector<LeptonSummary> selLeptons;

    //select muons
    std::vector<LeptonSummary> noIdMu;
    for(unsigned int muIter = 0; muIter < fForestMu.muPt->size(); ++muIter) {

      //kinematics selection
      TLorentzVector p4(0,0,0,0);
      p4.SetPtEtaPhiM(fForestMu.muPt->at(muIter),fForestMu.muEta->at(muIter),fForestMu.muPhi->at(muIter),0.1057);
      if(TMath::Abs(p4.Eta()) > lepEtaCut) continue;
      if(p4.Pt() < lepPtCut) continue;

      LeptonSummary l(13,p4);
      l.charge  = fForestMu.muCharge->at(muIter);
      l.chiso   = fForestMu.muPFChIso->at(muIter);
      l.nhiso   = fForestMu.muPFNeuIso->at(muIter);
      l.phoiso  = fForestMu.muPFPhoIso->at(muIter);
      l.rho     = getRho(pfColl,{1,2,3,4,5,6},p4.Eta()-0.5,p4.Eta()+0.5);
      l.chrho   = getRho(pfColl,{1,2,3},      p4.Eta()-0.5,p4.Eta()+0.5);
      l.phorho  = getRho(pfColl,{4},          p4.Eta()-0.5,p4.Eta()+0.5);
      l.nhrho   = getRho(pfColl,{5,6},        p4.Eta()-0.5,p4.Eta()+0.5);
      l.origIdx = muIter;
      noIdMu.push_back(l);

      //id (Tight muon requirements)
      int type=fForestMu.muType->at(muIter);
      bool isGlobal( ((type>>1)&0x1) );
      if(!isGlobal) continue;
      bool isPF( ((type>>5)&0x1) );
      if(!isPF) continue;
      bool isGlobalMuonPromptTight(fForestMu.muChi2NDF->at(muIter)<10. && fForestMu.muMuonHits->at(muIter)>0);
      if(!isGlobalMuonPromptTight) continue;
      if(fForestMu.muStations->at(muIter)<=1) continue;
      if(fForestMu.muTrkLayers->at(muIter) <= 5) continue;
      if(fForestMu.muPixelHits->at(muIter) == 0) continue;
      if(TMath::Abs(fForestMu.muInnerD0->at(muIter)) >=0.2 ) continue;
      if(TMath::Abs(fForestMu.muInnerDz->at(muIter)) >=0.5) continue;

      //selected a good muon
      selLeptons.push_back(l);
    }
    std::sort(noIdMu.begin(),noIdMu.end(),orderByPt);


    //monitor muon id variables
    if(noIdMu.size()>1) {
      TLorentzVector p4[2] = {noIdMu[0].p4,      noIdMu[1].p4};
      int midx[2]          = {noIdMu[0].origIdx, noIdMu[1].origIdx};
      int charge(noIdMu[0].charge*noIdMu[1].charge);
      if( fabs((p4[0]+p4[1]).M()-91)<15 ) {
        TString cat("zmmctrl");
        if(charge>0) cat="ss"+cat;
        for(size_t i=0; i<2; i++) {
          ht.fill("mmusta",    fForestMu.muStations->at(midx[i]),            plotWgt,cat);
          ht.fill("mtrklay",   fForestMu.muTrkLayers->at(midx[i]),           plotWgt,cat);
          ht.fill("mchi2ndf",  fForestMu.muChi2NDF->at(midx[i]),             plotWgt,cat);
          ht.fill("mmuhits",   fForestMu.muMuonHits->at(midx[i]),            plotWgt,cat);
          ht.fill("mpxhits",   fForestMu.muPixelHits->at(midx[i]),           plotWgt,cat);
          ht.fill("md0",       TMath::Abs(fForestMu.muInnerD0->at(midx[i])), plotWgt,cat);
          ht.fill("mdz",       TMath::Abs(fForestMu.muInnerDz->at(midx[i])), plotWgt,cat);          
        }
      }
    }
       
    //select electrons
    //cf. https://twiki.cern.ch/twiki/pub/CMS/HiHighPt2019/HIN_electrons2018_followUp.pdf
    std::vector<LeptonSummary> noIdEle;
    for(unsigned int eleIter = 0; eleIter < fForestEle.elePt->size(); ++eleIter) {
      
      //kinematics selection
      TLorentzVector p4(0,0,0,0);
      p4.SetPtEtaPhiM(fForestEle.elePt->at(eleIter),fForestEle.eleEta->at(eleIter),fForestEle.elePhi->at(eleIter),0.000511);

      //apply ad-hoc shift for endcap electrons if needed
      if(!isPP && fForestTree.run<=firstEEScaleShiftRun && TMath::Abs(p4.Eta())>=barrelEndcapEta[1])
        p4 *=eeScaleShift;         

      if(TMath::Abs(p4.Eta()) > lepEtaCut) continue;
      if(TMath::Abs(p4.Eta()) > barrelEndcapEta[0] && TMath::Abs(p4.Eta()) < barrelEndcapEta[1] ) continue;
      if(p4.Pt() < lepPtCut) continue;	      
      
      LeptonSummary l(11,p4);
      l.charge  = fForestEle.eleCharge->at(eleIter);
      l.chiso   = fForestEle.elePFChIso03->at(eleIter);
      l.nhiso   = fForestEle.elePFNeuIso03->at(eleIter);
      l.phoiso  = fForestEle.elePFPhoIso03->at(eleIter);
      l.rho     = getRho(pfColl,{1,2,3,4,5,6},p4.Eta()-0.5,p4.Eta()+0.5);
      l.chrho   = getRho(pfColl,{1,2,3},      p4.Eta()-0.5,p4.Eta()+0.5);
      l.phorho  = getRho(pfColl,{4},          p4.Eta()-0.5,p4.Eta()+0.5);
      l.nhrho   = getRho(pfColl,{5,6},        p4.Eta()-0.5,p4.Eta()+0.5);
      l.origIdx=eleIter;
      noIdEle.push_back(l);
      
      //electron id (separate for EB and EE, depending on centrality)
      //see https://indico.cern.ch/event/811374/contributions/3382348/attachments/1823016/2983559/HIN_electrons2018_centrDep.pdf
      bool isEB(TMath::Abs(p4.Eta()) <= barrelEndcapEta[0]);      
      if(isCentralEvent){
        if(fForestEle.eleSigmaIEtaIEta->at(eleIter)        >=(isEB ? 0.012 : 0.048)) continue;
        if(TMath::Abs(fForestEle.eledEtaAtVtx->at(eleIter))>=(isEB ? 0.006 : 0.008)) continue;
        if(TMath::Abs(fForestEle.eledPhiAtVtx->at(eleIter))>=(isEB ? 0.027 : 0.067)) continue;
        if(fForestEle.eleHoverE->at(eleIter)               >=(isEB ? 0.150 : 0.150)) continue;
        if(fForestEle.eleEoverPInv->at(eleIter)            >=(isEB ? 0.040 : 0.250)) continue;
        if(TMath::Abs(fForestEle.eleD0->at(eleIter))       >=(isEB ? 0.011 : 0.083)) continue;
        if(TMath::Abs(fForestEle.eleDz->at(eleIter))       >=(isEB ? 0.037 : 0.028)) continue;
      }
      else {
        if(fForestEle.eleSigmaIEtaIEta->at(eleIter)        >=(isEB ? 0.010 : 0.041)) continue;
        if(TMath::Abs(fForestEle.eledEtaAtVtx->at(eleIter))>=(isEB ? 0.007 : 0.011)) continue;
        if(TMath::Abs(fForestEle.eledPhiAtVtx->at(eleIter))>=(isEB ? 0.021 : 0.026)) continue;
        if(fForestEle.eleHoverE->at(eleIter)               >=(isEB ? 0.130 : 0.150)) continue;
        if(fForestEle.eleEoverPInv->at(eleIter)            >=(isEB ? 0.020 : 0.330)) continue;
        if(TMath::Abs(fForestEle.eleD0->at(eleIter))       >=(isEB ? 0.008 : 0.015)) continue;
        if(TMath::Abs(fForestEle.eleDz->at(eleIter))       >=(isEB ? 0.020 : 0.045)) continue;
      }

      //id'ed electron
      selLeptons.push_back(l);
    }
    std::sort(noIdEle.begin(),noIdEle.end(),orderByPt);       

    //monitor electron id variables
    if(noIdEle.size()>1) {
      TLorentzVector p4[2] = {noIdEle[0].p4,noIdEle[1].p4};
      int eidx[2]          = {noIdEle[0].origIdx,noIdEle[1].origIdx};
      int charge(noIdEle[0].charge*noIdEle[1].charge);
      if( fabs((p4[0]+p4[1]).M()-91)<15 ) {
        TString basecat("zeectrl");
        if(charge>0) basecat="ss"+basecat;
        for(size_t i=0; i<2; i++) {
          TString cat(basecat);
          cat += (fabs(p4[i].Eta())>=barrelEndcapEta[1] ? "EE" : "EB");
          ht.fill("esihih",  fForestEle.eleSigmaIEtaIEta->at(eidx[i]),         plotWgt,cat);
          ht.fill("edetavtx", TMath::Abs(fForestEle.eledEtaAtVtx->at(eidx[i])), plotWgt,cat);
          ht.fill("edphivtx", TMath::Abs(fForestEle.eledPhiAtVtx->at(eidx[i])), plotWgt,cat);
          ht.fill("ehoe",     fForestEle.eleHoverE->at(eidx[i]),                plotWgt,cat);
          ht.fill("eempinv",  fForestEle.eleEoverPInv->at(eidx[i]),             plotWgt,cat);
          ht.fill("ed0",      TMath::Abs(fForestEle.eleD0->at(eidx[i])),        plotWgt,cat);
          ht.fill("edz",      TMath::Abs(fForestEle.eleDz->at(eidx[i])),        plotWgt,cat);          
        }
      }
    }

    //sort selected electrons by pt
    std::sort(selLeptons.begin(),selLeptons.end(),orderByPt);

    //apply basic preselection
    if(mtrig+etrig==0) continue;
    if(selLeptons.size()<2) continue;
    TLorentzVector ll(selLeptons[0].p4+selLeptons[1].p4);
    t_llpt=ll.Pt();
    t_lleta=ll.Eta();
    t_llphi=ll.Phi();
    t_llm=ll.M();
    t_dphi=selLeptons[0].p4.DeltaPhi(selLeptons[1].p4);
    t_deta=fabs(selLeptons[0].p4.Eta()-selLeptons[1].p4.Eta());
    t_sumeta=selLeptons[0].p4.Eta()+selLeptons[1].p4.Eta();
    int dilCode(selLeptons[0].id*selLeptons[1].id);
    TString dilCat("mm");
    if(dilCode==11*13) dilCat="em";
    if(dilCode==11*11) dilCat="ee";

    if(blind) {
      bool isZ( fabs(t_llm-91)<15);
      int charge(selLeptons[0].charge*selLeptons[1].charge);
      if(!isMC && !isZ && charge<0 && fForestTree.run>=326887) continue;
    }

    //analyze jets
    std::vector<BtagInfo_t> pfJetsIdx,nodr_pfJetsIdx;
    std::vector<TLorentzVector> pfJetsP4,nodr_pfJetsP4;
    int npfjets(0),npfbjets(0); 
    for(int jetIter = 0; jetIter < fForestJets.nref; jetIter++){

      //at least two tracks
      if(fForestJets.trackN[jetIter]<2) continue;

      TLorentzVector jp4(0,0,0,0);
      jp4.SetPtEtaPhiM( fForestJets.jtpt[jetIter],fForestJets.jteta[jetIter],fForestJets.jtphi[jetIter],fForestJets.jtm[jetIter]);

      float csvVal=fForestJets.discr_csvV2[jetIter];
      int nsvtxTk=fForestJets.svtxntrk[jetIter];
      float msvtx=fForestJets.svtxm[jetIter];

      if(jp4.Pt()<30.) continue;
      if(fabs(jp4.Eta())>2.4) continue;
      bool isBTagged(csvVal>csvWP);      

      nodr_pfJetsIdx.push_back( std::make_tuple(nodr_pfJetsP4.size(),nsvtxTk,msvtx,csvVal) );
      nodr_pfJetsP4.push_back(jp4);

      if(jp4.DeltaR(selLeptons[0].p4)<0.4 || jp4.DeltaR(selLeptons[1].p4)<0.4) {
        pfJetsIdx.push_back(std::make_tuple(pfJetsP4.size(),nsvtxTk,msvtx,csvVal));
        pfJetsP4.push_back(jp4);
        npfjets++;
        npfbjets += isBTagged;
      }
    }
    std::sort(pfJetsIdx.begin(),       pfJetsIdx.end(),      orderByBtagInfo);
    std::sort(nodr_pfJetsIdx.begin(),  nodr_pfJetsIdx.end(), orderByBtagInfo);

    //define categories for pre-selection control histograms
    std::vector<TString> categs;
    categs.push_back(dilCat);

    std::vector<TString> addCategs;

    //monitor after run where EE scale shift changed
    if(!isPP){
      addCategs.clear();
      TString pf( fForestTree.run>=firstEEScaleShiftRun ? "after" : "before" );
      for(auto c : categs) {
        addCategs.push_back(c); addCategs.push_back(c+pf); 
      }
      categs=addCategs;
    }

    //monitor according to the b-tagging category
    addCategs.clear();
    TString pfbcat(Form("%dpfb",min(npfbjets,2)));
    for(auto c : categs) { 
      addCategs.push_back(c); 
      if(npfbjets==0) addCategs.push_back(c+"0pfb"); 
      if(npfbjets>0) addCategs.push_back(c+"geq1pfb"); 
    }
    categs=addCategs;

    
    //fill histogams
    for(int i=0; i<2; i++) {
      TString pf(Form("l%d",i+1));
      float pt(selLeptons[i].p4.Pt());
      ht.fill(pf+"pt",             pt,                            plotWgt, categs);
      ht.fill(pf+"eta",            fabs(selLeptons[i].p4.Eta()),  plotWgt, categs);
      ht.fill(pf+"chiso",          selLeptons[i].chiso,           plotWgt, categs);
      ht.fill(pf+"phoiso",         selLeptons[i].phoiso,          plotWgt, categs);
      ht.fill(pf+"nhiso",          selLeptons[i].nhiso,           plotWgt, categs);

    }

    ht.fill( "acopl",     1-fabs(t_dphi)/TMath::Pi(),                   plotWgt, categs);
    ht.fill( "detall",    t_deta,                                       plotWgt, categs);
    ht.fill( "drll",      selLeptons[0].p4.DeltaR(selLeptons[1].p4),    plotWgt, categs);
    ht.fill( "mll",       t_llm,                                        plotWgt, categs);
    ht.fill( "ptll",      t_llpt,                                       plotWgt, categs);
    ht.fill( "ptsum",     selLeptons[0].p4.Pt()+selLeptons[1].p4.Pt(),  plotWgt, categs);

    ht.fill( "chrho",  chrho[0],   plotWgt, categs);
    ht.fill( "phorho", phorho[0],  plotWgt, categs);
    ht.fill( "nhrho",  nhrho[0],   plotWgt, categs);

    //PF jets
    ht.fill( "npfjets",   npfjets,   plotWgt, categs);
    ht.fill( "npfbjets",  npfbjets,  plotWgt, categs);
    std::vector<TLorentzVector> pfFinalState;
    pfFinalState.push_back(selLeptons[0].p4);
    pfFinalState.push_back(selLeptons[1].p4);
    for(size_t ij=0; ij<min(pfJetsIdx.size(),size_t(2)); ij++) {     
      int idx(std::get<0>(pfJetsIdx[ij]));
      int ntks(std::get<1>(pfJetsIdx[ij]));
      float svm(std::get<2>(pfJetsIdx[ij]));
      float csv(std::get<3>(pfJetsIdx[ij]));
      TLorentzVector p4=pfJetsP4[idx];
      if(csv>csvWP) pfFinalState.push_back(p4);
      TString ppf(ij==0 ? "1" : "2");
      ht.fill( "pf"+ppf+"jbalance", p4.Pt()/ll.Pt(), plotWgt, categs);
      ht.fill( "pf"+ppf+"jpt",      p4.Pt(),         plotWgt, categs);
      ht.fill( "pf"+ppf+"jeta",     fabs(p4.Eta()),  plotWgt, categs);
      ht.fill( "pf"+ppf+"jsvtxm",   ntks,            plotWgt, categs);
      ht.fill( "pf"+ppf+"jsvtxntk", svm,             plotWgt, categs);
      ht.fill( "pf"+ppf+"jcsv",     csv,             plotWgt, categs);
      ht.fill2D( "pf"+ppf+"jetavsphi",   p4.Eta(),p4.Phi(),   plotWgt, categs);
    }

    std::vector<float> rapMoments=getRapidityMoments(pfFinalState);
    ht.fill( "pfrapavg",     rapMoments[0], plotWgt, categs);
    ht.fill( "pfraprms",     rapMoments[1], plotWgt, categs);
    ht.fill( "pfrapmaxspan", rapMoments[2], plotWgt, categs);
    float pfht(0.);
    TLorentzVector vis(0,0,0,0);
    for(auto p : pfFinalState) { vis+=p; pfht+=p.Pt(); }
    ht.fill( "pfht",         pfht, plotWgt, categs);
    ht.fill( "pfmht",        vis.Pt(), plotWgt, categs);

    // for tree filling set all the proper variables
    t_run    = fForestTree.run;
    t_lumi   = fForestTree.lumi;
    t_event  = fForestTree.evt;
    t_weight = plotWgt;
    t_cenbin = cenBin;
    t_etrig  = etrig;
    t_mtrig  = mtrig;

    for(size_t i=0; i<3; i++) {
      t_chrho[i]  = chrho[i];
      t_phorho[i] = phorho[i];
      t_nhrho[i]  = nhrho[i];
    }
 
    // fill the leptons ordered by pt
    t_lep_pt    .clear();
    t_lep_eta   .clear();
    t_lep_phi   .clear();
    t_lep_pdgId .clear();
    t_lep_charge.clear();
    t_lep_chiso.clear();    
    t_lep_phiso.clear();
    t_lep_nhiso.clear();
    t_lep_rho.clear();    
    t_lep_chrho.clear();    
    t_lep_phrho.clear();
    t_lep_nhrho.clear();
    t_nlep = selLeptons.size();
    for (int ilep = 0; ilep < t_nlep; ++ilep){
      t_lep_pt .push_back( selLeptons[ilep].p4.Pt()  );
      t_lep_eta.push_back( selLeptons[ilep].p4.Eta() );
      t_lep_phi.push_back( selLeptons[ilep].p4.Phi() );
      t_lep_chiso.push_back( selLeptons[ilep].chiso );
      t_lep_phiso.push_back( selLeptons[ilep].phoiso );
      t_lep_nhiso.push_back( selLeptons[ilep].nhiso );
      t_lep_rho.push_back( selLeptons[ilep].rho );
      t_lep_chrho.push_back( selLeptons[ilep].chrho );
      t_lep_phrho.push_back( selLeptons[ilep].phorho );
      t_lep_nhrho.push_back( selLeptons[ilep].nhrho );
      t_lep_pdgId .push_back( selLeptons[ilep].id );
      t_lep_charge.push_back( selLeptons[ilep].charge );
    }

    // fill the jets ordered by b-tag
    t_bjet_pt   .clear();
    t_bjet_eta  .clear();
    t_bjet_phi  .clear();
    t_bjet_mass .clear();
    t_bjet_csvv2.clear();
    t_nbjet = nodr_pfJetsIdx.size();
    for (int ij = 0; ij < t_nbjet; ij++) {
      int idx = std::get<0>(nodr_pfJetsIdx[ij]);
      t_bjet_pt   .push_back( nodr_pfJetsP4[idx].Pt()  );
      t_bjet_eta  .push_back( nodr_pfJetsP4[idx].Eta() );
      t_bjet_phi  .push_back( nodr_pfJetsP4[idx].Phi() );
      t_bjet_mass .push_back( nodr_pfJetsP4[idx].M()   );
      t_bjet_csvv2.push_back( std::get<3>(nodr_pfJetsIdx[ij])   );
    }


    t_ht  = pfht;
    t_mht = vis.Pt();

    // now set the 4 variables that we added for the tmva reader for the bdt evaluation
    bdt_l1pt      = t_lep_pt[0];
    bdt_apt       = (t_lep_pt[0]-t_lep_pt[1])/(t_lep_pt[0]+t_lep_pt[1]);
    bdt_abslleta  = fabs(t_lleta);
    bdt_dphilll2  = fabs(dphi_2(t_lep_pt[0],t_lep_eta[0],t_lep_phi[0],t_lep_pt[1],t_lep_eta[1],t_lep_phi[1],2)); // this function is in functions.cc in scripts/
    bdt_sumabseta = fabs(t_lep_eta[0])+fabs(t_lep_eta[1]);
    bdt_flavor    = abs(t_lep_pdgId[0]*t_lep_pdgId[1]); //abs should be fine here, it's an int
    t_apt         = bdt_apt;
    t_dphilll2    = bdt_dphilll2;
    t_bdt         = reader->EvaluateMVA( methodName );
    t_bdt_rarity  = reader->GetRarity  ( methodName );
    t_fisher2     = readerFisher2->EvaluateMVA( methodNameFisher2 );

    outTree->Fill();
    break;
  }

  //save histos to file  
  if(outURL!=""){
    TFile *fOut=TFile::Open(outURL,"RECREATE");
    fOut->cd();

    outTree->Write();

    //store the weight sum for posterior normalization
    if(isMC) {
      TH1D *wgtH=new TH1D("wgtsum","wgtsum",1,0,1);
      wgtH->SetBinContent(1,wgtSum);
      wgtH->SetDirectory(fOut);
      wgtH->Write();
    }
    for (auto& it : ht.getPlots())  { 
      if(it.second->GetEntries()==0) continue;
      it.second->SetDirectory(fOut); it.second->Write(); 
    }
    for (auto& it : ht.get2dPlots())  { 
      if(it.second->GetEntries()==0) continue;
      it.second->SetDirectory(fOut); it.second->Write(); 
    }
    fOut->Close();
  }

  return 0;
}

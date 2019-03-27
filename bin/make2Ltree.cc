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

#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "../scripts/functions.cc"


const bool isDebug = true;

enum ElectronIDs{LOOSE,TIGHT};
const int applyEleID=LOOSE;
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

//subtracts the UE pedestal from the charged isolatioin
float subtractedChIso(float chiso,float chrho){
  double ue=92.2-8765.8/(98.83+chrho);
  return max(chiso-ue,0.);
} 

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


typedef std::tuple<int,TLorentzVector> LeptonInfo_t;
static bool orderByPt(const LeptonInfo_t &a, const LeptonInfo_t &b)
{
  float pt_a(std::get<1>(a).Pt()), pt_b(std::get<1>(b).Pt());
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
  for(size_t i=0; i<2; i++) {
    TString pf(i==0 ? "tk" : "pf");

    ht.addHist(pf+"rapavg",      new TH1F(pf+"rapavg",     ";Average rapidity;Events",25,0,2.5));
    ht.addHist(pf+"raprms",      new TH1F(pf+"raprms",     ";#sigma(rapidity);Events",50,0,2.5));
    ht.addHist(pf+"rapmaxspan",  new TH1F(pf+"rapmaxspan", ";Maximum rapidity span;Events",25,0,5));
    ht.addHist(pf+"ht",          new TH1F(pf+"ht",         ";H_{T} [GeV];Events",25,0,500));
    ht.addHist(pf+"mht",         new TH1F(pf+"mht",        ";Missing H_{T} [GeV];Events",25,0,200));

    ht.addHist("n"+pf+"jets",    new TH1F("n"+pf+"jets",    ";Jet multiplicity;Events",8,0,8));
    ht.addHist("n"+pf+"bjets",   new TH1F("n"+pf+"bjets",   ";b-jet multiplicity;Events",5,0,5));    
    ht.addHist("n"+pf+"svtx",    new TH1F("n"+pf+"svtx",    ";Secondary vertex multiplicity;Events",5,0,5));
    for(size_t j=1; j<=2; j++){
      TString ppf(j==1 ? "1" : "2");
      ht.addHist(pf+ppf+"jbalance",    new TH1F(pf+ppf+"jbalance", ";R = p_{T}(j)/p_{T}(ll);Events",50,0,3));
      ht.addHist(pf+ppf+"jpt",         new TH1F(pf+ppf+"jpt",      ";Jet transverse momentum [GeV];Events",30,00,300));
      ht.addHist(pf+ppf+"jeta",        new TH1F(pf+ppf+"jeta",     ";Jet pseudo-rapidity;Events",20,0,2.5));
      ht.addHist(pf+ppf+"jetavsphi",   new TH2F(pf+ppf+"jetavsphi", ";Jet pseudo-rapidity;Jet azimuthal angle [rad];Events",100,-2.5,2.5,100,-TMath::Pi(),TMath::Pi()));
      ht.addHist(pf+ppf+"jsvtxm",      new TH1F(pf+ppf+"jsvtxm",   ";Secondary vertex mass;Events",25,0,6));
      ht.addHist(pf+ppf+"jsvtxntk",    new TH1F(pf+ppf+"jsvtxntk", ";Secondary vertex track multiplicity;Events",5,0,5));
      ht.addHist(pf+ppf+"jcsv",        new TH1F(pf+ppf+"jcsv",     ";CSVv2;Events",25,0,1));
    }
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
    hltTree_p->SetBranchStatus("HLT_HIL3Mu20_v1",1);
    hltTree_p->SetBranchAddress("HLT_HIL3Mu20_v1",&mtrig);
    hltTree_p->SetBranchStatus("HLT_HIEle20_WPLoose_Gsf_v1",1);
    hltTree_p->SetBranchAddress("HLT_HIEle20_WPLoose_Gsf_v1",&etrig);
  }else{
    hltTree_p->SetBranchStatus("HLT_HIL3Mu15_v1",1);
    hltTree_p->SetBranchAddress("HLT_HIL3Mu15_v1",&mtrig);
    hltTree_p->SetBranchStatus("HLT_HIEle20Gsf_v1",1);
    hltTree_p->SetBranchAddress("HLT_HIEle20Gsf_v1",&etrig);    
  }
  
  // =============================================================
  // marc here make the output tree

  TTree * outTree = new TTree("tree", "tree with 2lepton selection and combined collections");

  // event and trigger variables
  Int_t  t_run, t_lumi, t_etrig, t_mtrig;
  Long_t t_event;
  Float_t t_weight, t_cenbin, t_chrho, t_phorho, t_nhrho;
  outTree->Branch("run"   , &t_run  , "run/I");
  outTree->Branch("lumi"  , &t_lumi , "lumi/I");
  outTree->Branch("event" , &t_event, "event/L");

  outTree->Branch("weight", &t_weight, "weight/F");

  // centrality and 9 flavors of rho
  outTree->Branch("cenbin", &t_cenbin, "cenbin/F");

  outTree->Branch("chrho" , &t_chrho , "chrho/F");
  outTree->Branch("nhrho" , &t_nhrho , "nhrho/F");
  outTree->Branch("phorho", &t_phorho, "phorho/F");

  outTree->Branch("etrig" , &t_etrig , "etrig/I");
  outTree->Branch("mtrig" , &t_mtrig , "mtrig/I");

  // variables per lepton, including iso
  Int_t t_nlep;
  std::vector<Float_t> t_lep_pt, t_lep_eta, t_lep_phi, t_lep_phiso, t_lep_chiso, t_lep_chiso;
  std::vector<Int_t  > t_lep_pdgId, t_lep_charge;
  outTree->Branch("nlep"      , &t_nlep      , "nlep/I"            );
  outTree->Branch("lep_pt"    , &t_lep_pt    );
  outTree->Branch("lep_eta"   , &t_lep_eta   );
  outTree->Branch("lep_phi"   , &t_lep_phi   );
  outTree->Branch("lep_phiso" , &t_lep_phiso   );
  outTree->Branch("lep_chiso" , &t_lep_chiso   );
  outTree->Branch("lep_nhiso" , &t_lep_nhiso   );
  outTree->Branch("lep_pdgId" , &t_lep_pdgId );
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

    float evWgt(1.0);
    if(isMC && fForestTree.ttbar_w->size()) { evWgt=fForestTree.ttbar_w->at(0); }
    wgtSum += evWgt;    
    float plotWgt(evWgt);

    //first of all require a trigger
    int trig=etrig+mtrig;
    if(trig==0) continue;

    //apply global filters
    //if(!isPP){
    //  if(TMath::Abs(fForestTree.vz) > 15) continue;
    // }

    //build jets from different PF candidate collections   
    std::vector<PseudoJet> chCands,phoCands,nhCands;
    TLorentzVector p4(0,0,0,0);
    for(size_t ipf=0; ipf<fForestPF.pfId->size(); ipf++) {
      int id(abs(fForestPF.pfId->at(ipf)));

      //consider pions, photons and K0L
      float mass(0.13957);
      bool isCH(true),isPho(false),isNH(false);
      if(id==22) { mass=0.; isCH=false; isPho=true; isNH=false; }
      if(id==130 || id==2112 || id==1 || id==2) { mass=0.497; isCH=false; isPho=false; isNH=true; }

      p4.SetPtEtaPhiM(fForestPF.pfPt->at(ipf),fForestPF.pfEta->at(ipf),fForestPF.pfPhi->at(ipf),mass);

      //some basic kinematic cuts
      if(p4.Pt()<0.5) continue;
      if(fabs(p4.Eta())>2.5) continue;

      PseudoJet ip=PseudoJet(p4.Px(),p4.Py(),p4.Pz(),p4.E());
      ip.set_user_index( ipf );

      if(isCH)  chCands.push_back(ip);
      if(isPho) phoCands.push_back(ip);
      if(isNH)  nhCands.push_back(ip);

    }

    //rho estimators
    JetDefinition jet_def(antikt_algorithm, 0.4);
    ClusterSequence chcs(chCands, jet_def),phocs(phoCands,jet_def),nhcs(nhCands,jet_def);
    Selector sel_rapmax = SelectorAbsRapMax(2.4);
    JetDefinition jet_def_for_rho(kt_algorithm,0.5);
    AreaDefinition area_def(active_area,GhostedAreaSpec(2.4+1.));
    JetMedianBackgroundEstimator chbge(sel_rapmax, jet_def_for_rho, area_def);
    chbge.set_particles(chCands);
    JetMedianBackgroundEstimator phobge(sel_rapmax, jet_def_for_rho, area_def);
    phobge.set_particles(phoCands);
    JetMedianBackgroundEstimator nhbge(sel_rapmax, jet_def_for_rho, area_def);
    nhbge.set_particles(nhCands);
    float chrho=chbge.rho();
    float nhrho=nhbge.rho();
    float phorho=phobge.rho();


    //monitor trigger and centrality
    float cenBin=0;
    if(!isMC){
      cenBin=0.5*fForestTree.hiBin;
      Int_t runBin=lumiTool.getRunBin(fForestTree.run);
      Float_t lumi=lumiTool.getLumi(fForestTree.run);
      if(lumi>0.){
        if(etrig>0) ht.fill("ratevsrun",runBin,1./lumi,"e");
        if(mtrig>0) ht.fill("ratevsrun",runBin,1./lumi,"m");
      }
    }

    //select muons
    std::vector<LeptonInfo_t> mu,noIdMu;
    for(unsigned int muIter = 0; muIter < fForestMu.muPt->size(); ++muIter) {

      //kinematics selection
      TLorentzVector p4(0,0,0,0);
      p4.SetPtEtaPhiM(fForestMu.muPt->at(muIter),fForestMu.muEta->at(muIter),fForestMu.muPhi->at(muIter),0.1057);
      if(TMath::Abs(p4.Eta()) > lepEtaCut) continue;
      if(p4.Pt() < lepPtCut) continue;

      noIdMu.push_back(LeptonInfo_t(muIter,p4));

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
      mu.push_back(LeptonInfo_t(muIter,p4));
    }
    std::sort(noIdMu.begin(),noIdMu.end(),orderByPt);
    std::sort(mu.begin(),mu.end(),orderByPt);


    //monitor muon id variables
    if(noIdMu.size()>1) {
      TLorentzVector p4[2] = {std::get<1>(noIdMu[0]),std::get<1>(noIdMu[1])};
      int midx[2]          = {std::get<0>(noIdMu[0]),std::get<0>(noIdMu[1])};
      int charge(fForestMu.muCharge->at(midx[0])*fForestMu.muCharge->at(midx[1]));
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
    std::vector<LeptonInfo_t> ele,noIdEle;
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

      noIdEle.push_back(LeptonInfo_t(eleIter,p4));
      
      //electron id (separate for EB and EE)
      bool isEB(TMath::Abs(p4.Eta()) <= barrelEndcapEta[0]);
      if(applyEleID==TIGHT) {
        if(fForestEle.eleSigmaIEtaIEta->at(eleIter)>=(isEB?0.009630:0.042415)) continue;
        if(TMath::Abs(fForestEle.eledEtaAtVtx->at(eleIter))>=(isEB?0.00661:0.007572)) continue;
        if(TMath::Abs(fForestEle.eledPhiAtVtx->at(eleIter))>=(isEB?0.019314:0.016572)) continue;
        if(fForestEle.eleHoverE->at(eleIter)>=(isEB?0.142968:0.145809)) continue;
        if(fForestEle.eleEoverPInv->at(eleIter)>=(isEB?0.010806:0.008948)) continue;
        if(TMath::Abs(fForestEle.eleD0->at(eleIter))>=(isEB?0.007551:0.09230)) continue;
        if(TMath::Abs(fForestEle.eleDz->at(eleIter))>=(isEB?0.015361:0.015449)) continue;
      }
      else if(applyEleID==LOOSE) {
        if(fForestEle.eleSigmaIEtaIEta->at(eleIter)>=(isEB?0.011122:0.045489)) continue;
        if(TMath::Abs(fForestEle.eledEtaAtVtx->at(eleIter))>=(isEB?0.007627:0.014374)) continue;
        if(TMath::Abs(fForestEle.eledPhiAtVtx->at(eleIter))>=(isEB?0.028118:0.030909)) continue;
        if(fForestEle.eleHoverE->at(eleIter)>=(isEB?0.149803:0.147656)) continue;
        if(fForestEle.eleEoverPInv->at(eleIter)>=(isEB?0.027230:0.030291)) continue;
        if(TMath::Abs(fForestEle.eleD0->at(eleIter))>=(isEB?0.023967:0.018704)) continue;
        if(TMath::Abs(fForestEle.eleDz->at(eleIter))>=(isEB?0.031008:0.041391)) continue;
      }

      //selected electron
      ele.push_back(LeptonInfo_t(eleIter,p4));      
    }
    std::sort(noIdEle.begin(),noIdEle.end(),orderByPt);
    std::sort(ele.begin(),ele.end(),orderByPt);
       

    //monitor electron id variables
    if(noIdEle.size()>1) {
      TLorentzVector p4[2] = {std::get<1>(noIdEle[0]),std::get<1>(noIdEle[1])};
      int eidx[2]          = {std::get<0>(noIdEle[0]),std::get<0>(noIdEle[1])};
      int charge(fForestEle.eleCharge->at(eidx[0])*fForestEle.eleCharge->at(eidx[1]));
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


    //check if all electrons are barrel
    bool allEleInEB(true);
    for(size_t i=0; i<min(ele.size(),size_t(2)); i++){
      float abseta(fabs(std::get<1>(ele[i]).Eta()));
      if(abseta>barrelEndcapEta[1]) allEleInEB=false;
    }

    int nLep=mu.size()+ele.size();
    if(nLep<2) continue;

    std::vector<TLorentzVector> selLeptons;
    std::vector<int>            selLeptonsPdgIds;
    int dilCode(0);
    int charge(0);    
    TLorentzVector ll(0,0,0,0);
    vector< std::tuple<float,float,float> > liso;
    bool hasIsoElecs(true);
    if(mu.size()>1 && mtrig>0) {

      //muon final states from the muon PD only
      if(!isMC && !isSingleMuPD) continue;
      
      dilCode=13*13;
      selLeptons.push_back(std::get<1>(mu[0]));
      selLeptons.push_back(std::get<1>(mu[1]));
      ll=selLeptons[0]+selLeptons[1];
      int muIdx[2]={std::get<0>(mu[0]),std::get<0>(mu[1])};
      selLeptonsPdgIds.push_back(-13*fForestMu.muCharge->at(muIdx[0]));
      selLeptonsPdgIds.push_back(-13*fForestMu.muCharge->at(muIdx[1]));
      charge=fForestMu.muCharge->at(muIdx[0])*fForestMu.muCharge->at(muIdx[1]);
      for(size_t i=0; i<2; i++){
        float chIso( fForestMu.muPFChIso->at(muIdx[i]) );
        float phoIso( fForestMu.muPFPhoIso->at(muIdx[i]) );
        float neutIso( fForestMu.muPFNeuIso->at(muIdx[i]) );        
        liso.push_back( std::make_tuple(chIso,phoIso,neutIso) );
      }
    }
    else if(mu.size()>0 && ele.size()>0 && (etrig>0 || mtrig>0)) {

      //simultaneous triggering only in the single muon PD
      //to avoid double counting
      if(!isMC && etrig>0 && mtrig>0 && isSingleElePD) continue;
      
      dilCode=11*13;
      selLeptons.push_back(std::get<1>(mu[0]));
      selLeptons.push_back(std::get<1>(ele[0]));
      ll=selLeptons[0]+selLeptons[1];
      int lIdx[2]={std::get<0>(mu[0]),std::get<0>(ele[0])};
      selLeptonsPdgIds.push_back(-13*fForestMu .muCharge ->at(lIdx[0]));
      selLeptonsPdgIds.push_back(-11*fForestEle.eleCharge->at(lIdx[1]));
      charge=fForestMu.muCharge->at(lIdx[0])*fForestEle.eleCharge->at(lIdx[1]);            
      float chIso( fForestMu.muPFChIso->at(lIdx[0]) );
      float phoIso( fForestMu.muPFPhoIso->at(lIdx[0]) );
      float neutIso( fForestMu.muPFNeuIso->at(lIdx[0]) );        

      liso.push_back( std::make_tuple(chIso,phoIso,neutIso) );
      chIso=fForestEle.elePFChIso03->at(lIdx[1]);
      float chIsoP( subtractedChIso(chIso,chrho) );
      chIsoP=subtractedChIso(chIso,chrho);
      //hasIsoElecs &= (chIso<0.95*selLeptons[0].Pt());
      //hasIsoElecs &= (chIsoP<0.25*selLeptons[0].Pt());
      hasIsoElecs &= (chIsoP<8.5);
      phoIso=fForestEle.elePFPhoIso03->at(lIdx[1]);
      neutIso=fForestEle.elePFNeuIso03->at(lIdx[1]);        
      liso.push_back( std::make_tuple(chIso,phoIso,neutIso) );
    }
    else if(ele.size()>1 && etrig>0) {

      //ee final states from the electron PD only 
      if(!isMC && !isSingleElePD) continue;

      dilCode=11*11;
      selLeptons.push_back(std::get<1>(ele[0]));
      selLeptons.push_back(std::get<1>(ele[1]));
      ll=selLeptons[0]+selLeptons[1];
      int eleIdx[2]={std::get<0>(ele[0]),std::get<0>(ele[1])};
      selLeptonsPdgIds.push_back(-11*fForestEle.eleCharge->at(eleIdx[0]));
      selLeptonsPdgIds.push_back(-11*fForestEle.eleCharge->at(eleIdx[1]));
      charge=fForestEle.eleCharge->at(eleIdx[0])*fForestEle.eleCharge->at(eleIdx[1]);
      for(size_t i=0; i<2; i++){
        float chIso=fForestEle.elePFChIso03->at(eleIdx[i]);
        float chIsoP( subtractedChIso(chIso,chrho) );
        //hasIsoElecs &= (chIso<0.95*selLeptons[i].Pt());
        //hasIsoElecs &= (chIsoP<0.25*selLeptons[i].Pt());
        hasIsoElecs &= (chIsoP<8.5);
        float phoIso=fForestEle.elePFPhoIso03->at(eleIdx[i]);
        float neutIso=fForestEle.elePFNeuIso03->at(eleIdx[i]);        
        liso.push_back( std::make_tuple(chIso,phoIso,neutIso) );
      }
    }else{
      continue;
    }

    if(ll.M()<20) continue;
    bool isZ( dilCode!=11*13 && fabs(ll.M()-91)<15 );

    if(blind) {
      if(!isMC && !isZ && charge<0 && fForestTree.run>=326887) continue;
    }

    TString dilCat("em");
    if(dilCode==11*11) { dilCat=isZ ? "zee" : "ee"; }
    if(dilCode==13*13) { dilCat=isZ ? "zmm" : "mm"; }    
    if(charge>0) dilCat="ss"+dilCat;


    //analyze track jets
    std::vector<TLorentzVector> tkJetsP4;
    int ntkjets(0);
    std::vector<PseudoJet> tkjets = sorted_by_pt(chcs.inclusive_jets());
    for(auto j : tkjets) {
      if(j.constituents().size()<2) continue;
      TLorentzVector p4(j.px(),j.py(),j.pz(),j.e());
      if(p4.DeltaR(selLeptons[0])<0.4 || p4.DeltaR(selLeptons[1])<0.4) continue;
      if(fabs(p4.Eta())>2.4) continue;
      if(p4.Pt()<15) continue;
      tkJetsP4.push_back(p4);
      ntkjets ++;
    }
   
    //b-tag jet the track jets by matching in deltaR to PF jets
    std::vector<BtagInfo_t> matchedJetsIdx,pfJetsIdx;
    std::vector<TLorentzVector> pfJetsP4;
    int npfjets(0),npfbjets(0); 
    bool allPFBInEB(true);//,hasAwayPFJet(true);
    for(int jetIter = 0; jetIter < fForestJets.nref; jetIter++){

      //at least two tracks
      if(fForestJets.trackN[jetIter]<2) continue;

      TLorentzVector jp4(0,0,0,0);
      jp4.SetPtEtaPhiM( fForestJets.jtpt[jetIter],fForestJets.jteta[jetIter],fForestJets.jtphi[jetIter],fForestJets.jtm[jetIter]);
      if(jp4.DeltaR(selLeptons[0])<0.4 || jp4.DeltaR(selLeptons[1])<0.4) continue;            

      float csvVal=fForestJets.discr_csvV2[jetIter];
      int nsvtxTk=fForestJets.svtxntrk[jetIter];
      float msvtx=fForestJets.svtxm[jetIter];

      for(size_t ij=0; ij<tkJetsP4.size(); ij++) {
        if(jp4.DeltaR( tkJetsP4[ij] ) >0.4) continue;
        matchedJetsIdx.push_back(std::make_tuple(ij,nsvtxTk,msvtx,csvVal));
        break;
      }

      if(jp4.Pt()<30.) continue;
      if(fabs(jp4.Eta())>2.4) continue;
      bool isBTagged(csvVal>csvWP);      
      pfJetsIdx.push_back(std::make_tuple(pfJetsP4.size(),nsvtxTk,msvtx,csvVal));
      pfJetsP4.push_back(jp4);
      npfjets++;
      npfbjets += isBTagged;
      if(isBTagged && fabs(jp4.Eta())>1.2) allPFBInEB=false;

      /*
      if(npfjets==1){
        float dr2ll(ll.DeltaR(jp4));
        if(fabs(dr2ll)>2*TMath::Pi()/3.) hasAwayPFJet=true;
      }
      */
    }
    std::sort(pfJetsIdx.begin(),      pfJetsIdx.end(),      orderByBtagInfo);

    //finalize analysing track jets
    std::sort(matchedJetsIdx.begin(), matchedJetsIdx.end(), orderByBtagInfo);
    // bool allTkBInEB(true),hasAwayTkJet(false);
    int nbtkjets(0);
    for(size_t ij=0; ij<min(matchedJetsIdx.size(),size_t(2)); ij++) {     
      int idx(std::get<0>(matchedJetsIdx[ij]));
      float csv(std::get<3>(matchedJetsIdx[ij]));
      TLorentzVector p4=tkJetsP4[idx];
      tkJetsP4.push_back(p4);
      bool isBTagged(csv>csvWP);
      nbtkjets += isBTagged;
      //if(isBTagged && fabs(p4.Eta())>1.2) allTkBInEB=false;
      
      /*
      if(ij==0){
        float dr2ll(ll.DeltaR(p4));
        if(fabs(dr2ll)>2*TMath::Pi()/3.) hasAwayTkJet=true;
      }
      */
    }

    //fill control histograms
    std::vector<TString> categs;
    categs.push_back(dilCat);

    if(ll.Pt()>20)
      categs.push_back(dilCat+"hpur");

    //monitor where the electrons are reconstructed
    if(dilCode==11*11 || dilCode==11*13){

      bool l1EE(fabs(selLeptons[0].Eta())>barrelEndcapEta[1]);
      bool l2EE(fabs(selLeptons[1].Eta())>barrelEndcapEta[1]);
      TString etaCateg(l2EE ? "E" : "B");
      if(dilCode==11*11) etaCateg += l1EE ? "E" : "B";
      categs.push_back(dilCat+etaCateg);

      /*
      if(isZ) {
        if(ntkjets==1 && hasAwayTkJet) {
          categs.push_back(dilCat+"awaytkj");
          categs.push_back(dilCat+etaCateg+"awaytkj");
        }
        if(npfjets==1 && hasAwayPFJet) {
          categs.push_back(dilCat+"awaypfj");
          categs.push_back(dilCat+etaCateg+"awaypfj");
        }
      }
      */
    }
    

    std::vector<TString> addCategs;

    //monitor also after run where EE scale shift changed
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
    //TString tkbcat(Form("%dtkb",min(nbtkjets,2)));
    for(auto c : categs) { 
      addCategs.push_back(c); 
      //addCategs.push_back(c+tkbcat); 
      if(npfbjets==0) addCategs.push_back(c+"0pfb"); 
      if(npfbjets>0) addCategs.push_back(c+"geq1pfb"); 
      if(npfbjets==1) addCategs.push_back(c+"eq1pfb"); 
      if(npfbjets>1) addCategs.push_back(c+"geq2pfb");         
    }
    categs=addCategs;

    //monitor according to the centrality of the jets and electrons
    if(allEleInEB) {
      addCategs.clear();
      for(auto c: categs) {
        addCategs.push_back(c);
        //if(allTkBInEB) addCategs.push_back(c+"alltkeb");
        if(allPFBInEB) addCategs.push_back(c+"allpfeb");
      }
      categs=addCategs;
    }

    for(int i=0; i<2; i++) {
      TString pf(Form("l%d",i+1));
      float pt(selLeptons[i].Pt());
      ht.fill(pf+"pt",        pt,                         plotWgt, categs);
      ht.fill(pf+"eta",       fabs(selLeptons[i].Eta()),  plotWgt, categs);
      float chiso(std::get<0>(liso[i]));
      float phoiso(std::get<1>(liso[i]));
      float nhiso(std::get<2>(liso[i]));

      ht.fill(pf+"chiso",          chiso,     plotWgt, categs);
      ht.fill(pf+"chreliso",       chiso/pt,  plotWgt, categs);
      ht.fill2D(pf+"chisovscen",   cenBin,    chiso,   plotWgt, categs);
      ht.fill2D(pf+"chisovsrho",   chrho,     chiso,   plotWgt, categs);

      ht.fill(pf+"phoiso",        phoiso,     plotWgt,  categs);
      ht.fill(pf+"phoreliso",     phoiso/pt,  plotWgt,  categs);
      ht.fill2D(pf+"phoisovscen", cenBin,     phoiso,   plotWgt, categs);
      ht.fill2D(pf+"phoisovsrho", phorho,     phoiso,   plotWgt, categs);

      ht.fill(pf+"nhiso",        nhiso,     plotWgt, categs);
      ht.fill(pf+"nhreliso",     nhiso/pt,  plotWgt, categs);
      ht.fill2D(pf+"nhisovscen", cenBin,    nhiso,   plotWgt, categs);
      ht.fill2D(pf+"nhisovsrho", nhrho,     nhiso,   plotWgt, categs);
    }

    ht.fill( "acopl",     1-fabs(selLeptons[0].DeltaPhi(selLeptons[1]))/TMath::Pi(), plotWgt, categs);
    ht.fill( "detall",    fabs(selLeptons[0].Eta()-selLeptons[1].Eta()), plotWgt, categs);
    ht.fill( "drll",      selLeptons[0].DeltaR(selLeptons[1]),           plotWgt, categs);
    ht.fill( "mll",       ll.M(),                                        plotWgt, categs);
    ht.fill( "ptll",      ll.Pt(),                                       plotWgt, categs);
    ht.fill( "ptsum",     selLeptons[0].Pt()+selLeptons[1].Pt(),         plotWgt, categs);

    ht.fill( "chrho",  chrho,            plotWgt, categs);
    ht.fill( "phorho", phorho,            plotWgt, categs);
    ht.fill( "nhrho",  nhrho,            plotWgt, categs);

    //track jets
    std::vector<TLorentzVector> tkFinalState;
    ht.fill( "ntkjets",   ntkjets,                                     plotWgt, categs);
    ht.fill( "ntkbjets",  nbtkjets,                                    plotWgt, categs);
    tkFinalState.push_back(selLeptons[0]);
    tkFinalState.push_back(selLeptons[1]);
    for(size_t ij=0; ij<min(matchedJetsIdx.size(),size_t(2)); ij++) {     
      int idx(std::get<0>(matchedJetsIdx[ij]));
      int ntks(std::get<1>(matchedJetsIdx[ij]));
      float svm(std::get<2>(matchedJetsIdx[ij]));
      float csv(std::get<3>(matchedJetsIdx[ij]));
      TLorentzVector p4=tkJetsP4[idx];
      if(csv>csvWP) tkFinalState.push_back(p4);

      TString ppf(ij==0 ? "1" : "2");
      ht.fill( "tk"+ppf+"jbalance",  p4.Pt()/ll.Pt(),  plotWgt, categs);
      ht.fill( "tk"+ppf+"jpt",      p4.Pt(),          plotWgt, categs);
      ht.fill( "tk"+ppf+"jeta",     fabs(p4.Eta()),   plotWgt, categs);
      ht.fill2D( "tk"+ppf+"jetavsphi",   p4.Eta(),p4.Phi(),   plotWgt, categs);
      ht.fill( "tk"+ppf+"jsvtxm",   ntks,             plotWgt, categs);
      ht.fill( "tk"+ppf+"jsvtxntk", svm,              plotWgt, categs);
      ht.fill( "tk"+ppf+"jcsv",     csv,              plotWgt, categs);
      
    }    

    std::vector<float> rapMoments=getRapidityMoments(tkFinalState);
    ht.fill( "tkrapavg",     rapMoments[0], plotWgt, categs);
    ht.fill( "tkraprms",     rapMoments[1], plotWgt, categs);
    ht.fill( "tkrapmaxspan", rapMoments[2], plotWgt, categs);
    float tkht(0.);
    TLorentzVector tkvis(0,0,0,0);
    for(auto p : tkFinalState) { tkvis+=p; tkht+=p.Pt(); }
    ht.fill( "tkht",  tkht, plotWgt, categs);
    ht.fill( "tkmht", tkvis.Pt(), plotWgt, categs);


    //PF jets
    ht.fill( "npfjets",   npfjets,                                     plotWgt, categs);
    ht.fill( "npfbjets",  npfbjets,                                    plotWgt, categs);
    std::vector<TLorentzVector> pfFinalState;
    pfFinalState.push_back(selLeptons[0]);
    pfFinalState.push_back(selLeptons[1]);
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
      ht.fill2D( "pf"+ppf+"jetavsphi",   p4.Eta(),p4.Phi(),   plotWgt, categs);
      ht.fill( "pf"+ppf+"jsvtxm",   ntks,            plotWgt, categs);
      ht.fill( "pf"+ppf+"jsvtxntk", svm,             plotWgt, categs);
      ht.fill( "pf"+ppf+"jcsv",     csv,             plotWgt, categs);
    }

    rapMoments=getRapidityMoments(pfFinalState);
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

    t_chrho  = chrho;
    t_phorho = phorho;
    t_nhrho  = nhrho;
    t_etrig = etrig;
    t_mtrig = mtrig;
    t_llpt   = ll.Pt();
    t_lleta  = ll.Eta();
    t_llphi  = ll.Phi();
    t_llm    = ll.M();
    t_dphi   = fabs(selLeptons[0].DeltaPhi(selLeptons[1]));
    t_deta   = fabs(selLeptons[0].Eta()-selLeptons[1].Eta());
    t_sumeta = selLeptons[0].Eta()+selLeptons[1].Eta();

    // fill the leptons ordered by pt
    t_lep_pt    .clear();
    t_lep_eta   .clear();
    t_lep_phi   .clear();
    t_lep_pdgId .clear();
    t_lep_charge.clear();

    t_nlep = selLeptons.size();
    bool revertOrder = selLeptons[1].Pt() > selLeptons[0].Pt() && dilCode == 11*13;
    if (!revertOrder) {
      for (int ilep = 0; ilep < t_nlep; ++ilep){
        t_lep_pt .push_back( selLeptons[ilep].Pt()  );
        t_lep_eta.push_back( selLeptons[ilep].Eta() );
        t_lep_phi.push_back( selLeptons[ilep].Phi() );
        t_lep_phiso.push_back( 10. );
        t_lep_chiso.push_back( 20. );
        t_lep_nhiso.push_back( 30. );
        t_lep_pdgId .push_back( selLeptonsPdgIds[ilep] );
        t_lep_charge.push_back( TMath::Power(-1.,selLeptonsPdgIds[ilep] > 0));
      }
    }
    else {
      for (int ilep = t_nlep-1; ilep >= 0; --ilep){
        t_lep_pt .push_back( selLeptons[ilep].Pt()  );
        t_lep_eta.push_back( selLeptons[ilep].Eta() );
        t_lep_phi.push_back( selLeptons[ilep].Phi() );
        t_lep_phiso.push_back( 10. );
        t_lep_chiso.push_back( 20. );
        t_lep_nhiso.push_back( 30. );
        t_lep_pdgId .push_back( selLeptonsPdgIds[ilep] );
        t_lep_charge.push_back( TMath::Power(-1.,selLeptonsPdgIds[ilep] > 0));
      }
    }

    // fill the jets ordered by b-tag
    t_bjet_pt   .clear();
    t_bjet_eta  .clear();
    t_bjet_phi  .clear();
    t_bjet_mass .clear();
    t_bjet_csvv2.clear();

    t_nbjet = pfJetsIdx.size();
    for (int ij = 0; ij < t_nbjet; ij++) {
      int idx = std::get<0>(pfJetsIdx[ij]);
      t_bjet_pt   .push_back( pfJetsP4[idx].Pt()  );
      t_bjet_eta  .push_back( pfJetsP4[idx].Eta() );
      t_bjet_phi  .push_back( pfJetsP4[idx].Phi() );
      t_bjet_mass .push_back( pfJetsP4[idx].M()   );
      t_bjet_csvv2.push_back( std::get<3>(pfJetsIdx[ij])   );
    }


    // std::vector<Int_t  > t_lep_pdgId, t_lep_charge;
    // t_njet;
    // std::vector<Float_t> t_jet_pt, t_jet_eta, t_jet_phi, t_jet_mass, t_jet_csvv2;
    t_ht  = pfht;
    t_mht = vis.Pt();

    // now set the 4 variables that we added for the tmva reader for the bdt evaluation
    bdt_l1pt = t_lep_pt[0];
    bdt_apt =  (t_lep_pt[0]-t_lep_pt[1])/(t_lep_pt[0]+t_lep_pt[1]);
    bdt_abslleta = fabs(t_lleta);
    bdt_dphilll2 = fabs(dphi_2(t_lep_pt[0],t_lep_eta[0],t_lep_phi[0],t_lep_pt[1],t_lep_eta[1],t_lep_phi[1],2)); // this function is in functions.cc in scripts/

    bdt_sumabseta = fabs(t_lep_eta[0])+fabs(t_lep_eta[1]);
    bdt_flavor = abs(t_lep_pdgId[0]*t_lep_pdgId[1]); //abs should be fine here, it's an int


    t_apt = bdt_apt;
    t_dphilll2 = bdt_dphilll2;
    t_bdt        = reader->EvaluateMVA( methodName );
    t_bdt_rarity = reader->GetRarity  ( methodName );
    t_fisher2    = readerFisher2->EvaluateMVA( methodNameFisher2 );

    outTree->Fill();

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

#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TChain.h"

#include <string>
#include <vector>

#include "simpleTT2L.h"
#include "ForestTreeHeaders/MuonIdCriteria.h"
#include "ForestTreeHeaders/ElectronIdCriteria.h"
#include "ForestTreeHeaders/ForestElectrons.h"
#include "ForestTreeHeaders/ForestMuons.h"
#include "ForestTreeHeaders/HistTool.h"

const bool isDebug = true;

const float jetPtCut  = 30.;
const float jetEtaCut = 2.4;
const float lepPtCut  = 20.;
const float lepEtaCut = 2.1;
const float barrelEndcapEta[2]={1.4442,1.5660};
const float csvWP = 0.8838;

void simpleTT2L(const std::string outFileName = "", const std::string inFileName = "", bool isMC = false, bool isPP=true,bool doSameSign=false)
{
  if(!strcmp(inFileName.c_str(), "")){
    std::cout << "No inputs specified. return" << std::endl;
    return;
  }

  if(isDebug) std::cout << __LINE__ << std::endl;

  HistTool ht;
  ht.addHist("lpt",     new TH1F("lpt",    ";Lepton transverse momentum [GeV];Events",20,20,200));
  ht.addHist("leta",    new TH1F("leta",   ";Lepton pseudo-rapidity;Events",20,0,2.5));
  ht.addHist("mll",     new TH1F("mll",    ";Dilepton invariant mass [GeV];Events",20,20,200));
  ht.addHist("ptll",    new TH1F("ptll",    ";Dilepton transverse momentum [GeV];Events",20,0,200));
  ht.addHist("dphill",  new TH1F("dphill", ";#Delta#phi(l,l');Events",20,0,3.15));
  ht.addHist("njets",   new TH1F("njets",  ";Jet multiplicity;Events",5,0,5));
  ht.addHist("nbjets",  new TH1F("nbjets", ";b-jet multiplicity;Events",5,0,5));
  ht.addHist("jpt",     new TH1F("jpt",    ";Jet transverse momentum [GeV];Events",20,30,200));
  ht.addHist("jeta",    new TH1F("jeta",   ";Jet pseudo-rapidity;Events",20,0,2.5));
  ht.addHist("jcsv",    new TH1F("jcsv",   ";CSVv2;Events",25,0,1));

  //electron id histograms
  ht.addHist("esihih",     new TH1F("esihih",   ";#sigma(i#eta,i#eta);Electrons",20,0,0.03));
  ht.addHist("edetain",    new TH1F("edetain",  ";#Delta#eta(in);Electrons",20,0,0.2));
  ht.addHist("ephiin",     new TH1F("ephiin",   ";#Delta#phi(in) [rad];Electrons",20,0,0.4));
  ht.addHist("ehoe",       new TH1F("ehoe",     ";h/e;Electrons",20,0,0.2));
  ht.addHist("ed0",        new TH1F("ed0",      ";d_{0} [cm];Electrons",20,0,0.2));
  ht.addHist("edz",        new TH1F("edz",      ";d_{z} [cm];Electrons",20,0,0.2));
  ht.addHist("emisshits",  new TH1F("emisshits",";Missing inner hits;Electrons",3,0,3));
  ht.addHist("eeop",       new TH1F("eeop",     ";|1/E-1/p|;Electrons",20,0,0.5));
  ht.addHist("echreliso",  new TH1F("echreliso",";Relative PF charged isolation;Electrons",20,0,2.0));
  ht.addHist("ephoreliso", new TH1F("ephoreliso",";Relative PF photon isolation;Electrons",20,0,1.0));
  ht.addHist("eneureliso", new TH1F("eneureliso",";Relative PF neutral hadron isolation;Electrons",20,0,1.0));
 
  std::vector<std::string>* inFileNames_p = new std::vector<std::string>;
  inFileNames_p->push_back(inFileName);

  setEleIdCuts();
  
  TChain *lepTree_p     = new TChain(isPP ? "ggHiNtuplizer/EventTree" : "ggHiNtuplizerGED/EventTree");
  TChain *jetTree_p     = new TChain(isPP ? "ak4PFJetAnalyzer/t"      : "akPu4CaloJetAnalyzer/t");
  TChain *hiTree_p      = new TChain("hiEvtAnalyzer/HiTree");
  TChain *hltTree_p     = new TChain("hltanalysis/HltTree");
  TChain *skimAnaTree_p = new TChain("skimanalysis/HltTree");
  
  const int nFiles = (int)inFileNames_p->size();

  for(int fileIter = 0; fileIter < nFiles; fileIter++){
    std::cout << "On file: " << fileIter << "/" << nFiles << "  " << inFileNames_p->at(fileIter).c_str()  << std::endl;
    lepTree_p->Add(inFileNames_p->at(fileIter).c_str());
    jetTree_p->Add(inFileNames_p->at(fileIter).c_str());
    hiTree_p->Add(inFileNames_p->at(fileIter).c_str());
    hltTree_p->Add(inFileNames_p->at(fileIter).c_str());
    skimAnaTree_p->Add(inFileNames_p->at(fileIter).c_str());
  }

  ForestMuons fForestMu;
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
  
  ForestElectrons fForestEle;
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
  lepTree_p->SetBranchAddress("eleMissHits", &fForestEle.eleMissHits);
  lepTree_p->SetBranchAddress("elePFChIso03", &fForestEle.elePFChIso03);
  lepTree_p->SetBranchAddress("elePFPhoIso03", &fForestEle.elePFPhoIso03);
  lepTree_p->SetBranchAddress("elePFNeuIso03", &fForestEle.elePFNeuIso03);
    
  const int maxJets = 5000;
  int           nref;
  float         jtpt[maxJets];
  float         jteta[maxJets];
  float         jtphi[maxJets];
  float         jtm[maxJets];
  float         discr_csv[maxJets];
  int           trackN[maxJets];
  jetTree_p->SetBranchStatus("*", 0);
  jetTree_p->SetBranchStatus("nref", 1);
  jetTree_p->SetBranchStatus("jtpt", 1);
  jetTree_p->SetBranchStatus("jtphi", 1);
  jetTree_p->SetBranchStatus("jteta", 1);
  jetTree_p->SetBranchStatus("jtm", 1);
  jetTree_p->SetBranchStatus("discr_csvV2", 1);
  jetTree_p->SetBranchStatus("trackN", 1);
  jetTree_p->SetBranchAddress("nref", &nref);
  jetTree_p->SetBranchAddress("jtpt", jtpt);
  jetTree_p->SetBranchAddress("jtphi", jtphi);
  jetTree_p->SetBranchAddress("jteta", jteta);
  jetTree_p->SetBranchAddress("jtm", jtm);
  jetTree_p->SetBranchAddress("discr_csvV2", discr_csv);
  jetTree_p->SetBranchAddress("trackN", trackN);
  
  UInt_t run, lumi;
  ULong64_t evt;
  Int_t hiBin;
  Float_t hiHF,vz;
  Float_t weight;
  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("run", 1);
  hiTree_p->SetBranchStatus("evt", 1);
  hiTree_p->SetBranchStatus("lumi", 1);
  hiTree_p->SetBranchStatus("hiBin", 1);
  hiTree_p->SetBranchStatus("hiHF", 1);
  hiTree_p->SetBranchStatus("vz", 1);
  hiTree_p->SetBranchStatus("weight", 1);
  hiTree_p->SetBranchAddress("run", &run);
  hiTree_p->SetBranchAddress("evt", &evt);
  hiTree_p->SetBranchAddress("lumi", &lumi);
  hiTree_p->SetBranchAddress("hiBin", &hiBin);
  hiTree_p->SetBranchAddress("hiHF", &hiHF);
  hiTree_p->SetBranchAddress("vz", &vz);
  hiTree_p->SetBranchAddress("weight",&weight);

  int eetrig(0),emtrig(1);
  if(isPP){
    hltTree_p->SetBranchStatus("HLT_HIL3Mu20_v1",1);
    hltTree_p->SetBranchAddress("HLT_HIL3Mu20_v1",&emtrig);
    hltTree_p->SetBranchStatus("HLT_HIEle20_WPLoose_Gsf_v1",1);
    hltTree_p->SetBranchAddress("HLT_HIEle20_WPLoose_Gsf_v1",&eetrig);
  }else{
    eetrig=inFileName.find("Data2018PbPb_ZEE")!=std::string::npos ? 1 : 0;
    emtrig=inFileName.find("Data2018PbPb_ZEE")!=std::string::npos ? 0 : 1;
    // hltTree_p->SetBranchStatus("HLT_HIL1Mu5Eta2p5_Ele20Gsf_v",1);
    // hltTree_p->SetBranchAddress("HLT_HIL1Mu5Eta2p5_Ele20Gsf_v",&emtrig);
    // hltTree_p->SetBranchStatus("HLT_PADoublePhoton15_Eta3p1_Mass50_1000_v",1);
    // hltTree_p->SetBranchAddress("HLT_PADoublePhoton15_Eta3p1_Mass50_1000_v",&eetrig);
  }

  
  //event selections
  int phfCoincFilter = 1;
  int HBHENoiseFilterResult = 1;
  int pprimaryVertexFilter = 1;
  int pcollisionEventSelection = 1;
  skimAnaTree_p->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
  skimAnaTree_p->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult);
  skimAnaTree_p->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
  skimAnaTree_p->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  
  
  Float_t wgtSum(0);
  int nEntries = (int)lepTree_p->GetEntries();  
  int entryDiv = ((int)(nEntries/20));    
  for(int entry = 0; entry < nEntries; entry++){
    
    if(entry%entryDiv == 0 && nEntries >= 10000) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
    
    lepTree_p->GetEntry(entry);
    jetTree_p->GetEntry(entry);    
    hltTree_p->GetEntry(entry);
    hiTree_p->GetEntry(entry);
    skimAnaTree_p->GetEntry(entry);
    wgtSum += weight;

    int trig=emtrig+eetrig;
    if(trig==0) continue;

    if(!isPP){
      //if(!phfCoincFilter) continue;
      //if(!HBHENoiseFilterResult) continue;
      //if(!pcollisionEventSelection) continue;
      //if(!pprimaryVertexFilter) continue;
      if(TMath::Abs(vz) > 15) continue;
    }

    //use cuts corresponding to most peripheral, if PP
    int centEleId = isPP ? 3 : getCentBinEleId((double)hiBin/2.);
    
    std::vector<int> muIdx,eleIdx,noIdEleIdx;    
    std::vector<TLorentzVector> muP4,eP4;
    for(unsigned int muIter = 0; muIter < fForestMu.muPt->size(); ++muIter) {
      if(TMath::Abs(fForestMu.muEta->at(muIter)) > lepEtaCut) continue;
      if(fForestMu.muPt->at(muIter) < lepPtCut) continue;
      if(fForestMu.muChi2NDF->at(muIter) >= muChi2NDFCut) continue;
      if(TMath::Abs(fForestMu.muInnerD0->at(muIter)) > muInnerD0Cut) continue;
      if(TMath::Abs(fForestMu.muInnerDz->at(muIter)) > muInnerDzCut) continue;
      if(fForestMu.muMuonHits->at(muIter) <= muMuonHitsCut) continue;
      if(fForestMu.muStations->at(muIter) <= muStationsCut) continue;
      if(fForestMu.muTrkLayers->at(muIter) <= muTrkLayersCut) continue;
      if(fForestMu.muPixelHits->at(muIter) <= muPixelHitsCut) continue;
      muIdx.push_back(muIter);
      muP4.push_back(TLorentzVector(0,0,0,0) );
      muP4[muP4.size()-1].SetPtEtaPhiM(fForestMu.muPt->at(muIter),
                                       fForestMu.muEta->at(muIter),
                                       fForestMu.muPhi->at(muIter),
                                       muM);
    }
    
    for(unsigned int eleIter = 0; eleIter < fForestEle.elePt->size(); ++eleIter) {
      if(TMath::Abs(fForestEle.eleEta->at(eleIter)) > lepEtaCut) continue;
      if(fForestEle.elePt->at(eleIter) < lepPtCut) continue;	      
      if(TMath::Abs(fForestEle.eleEta->at(eleIter)) > barrelEndcapEta[0]
         && TMath::Abs(fForestEle.eleEta->at(eleIter)) < barrelEndcapEta[1] ) continue;
      noIdEleIdx.push_back(eleIter);

      unsigned int eleEtaCutPos(TMath::Abs(fForestEle.eleEta->at(eleIter))<barrelEndcapEta[0] ? 0 : 1);
      if(fForestEle.eleSigmaIEtaIEta->at(eleIter) > eleSigmaIEtaIEta_VetoCut[eleEtaCutPos][centEleId]) continue;
      if(TMath::Abs(fForestEle.eledEtaAtVtx->at(eleIter)) > eleDEtaIn_VetoCut[eleEtaCutPos][centEleId]) continue;
      if(TMath::Abs(fForestEle.eledPhiAtVtx->at(eleIter)) > eleDPhiIn_VetoCut[eleEtaCutPos][centEleId]) continue;
      if(fForestEle.eleHoverE->at(eleIter) > eleHOverE_VetoCut[eleEtaCutPos][centEleId]) continue;
      if(TMath::Abs(fForestEle.eleD0->at(eleIter)) > eleD0_VetoCut[eleEtaCutPos][centEleId]) continue;
      if(TMath::Abs(fForestEle.eleDz->at(eleIter)) > eleDZ_VetoCut[eleEtaCutPos][centEleId]) continue;
      if(TMath::Abs(fForestEle.eleEoverPInv->at(eleIter)) > eleEoverPInv_VetoCut[eleEtaCutPos][centEleId]) continue;
      eleIdx.push_back(eleIter);
      eP4.push_back(TLorentzVector(0,0,0,0) );
      eP4[eP4.size()-1].SetPtEtaPhiM(fForestEle.elePt->at(eleIter),
                                     fForestEle.eleEta->at(eleIter),
                                     fForestEle.elePhi->at(eleIter),
                                     eleM);
    }

    //z->ee control region for electron id variables
    if(noIdEleIdx.size()>1) {
      TLorentzVector e1(0,0,0,0),e2(0,0,0,0);
      e1.SetPtEtaPhiM(fForestEle.elePt->at(noIdEleIdx[0]),fForestEle.eleEta->at(noIdEleIdx[0]),fForestEle.elePhi->at(noIdEleIdx[0]),eleM);
      e2.SetPtEtaPhiM(fForestEle.elePt->at(noIdEleIdx[1]),fForestEle.eleEta->at(noIdEleIdx[1]),fForestEle.elePhi->at(noIdEleIdx[1]),eleM);
      float charge(fForestEle.eleCharge->at(noIdEleIdx[0])*fForestEle.eleCharge->at(noIdEleIdx[1]));
      if(doSameSign) charge*=-1;
      float mll((e1+e2).M());
      if(fabs(mll-91)<30 && charge<0 && eetrig>0 && fabs(e1.Eta())<1.44 && fabs(e2.Eta())<1.44) {

        std::vector<TString> baseCat(1,"zeectrl");
        if(isPP) {
          baseCat.push_back("zeectrlhighcen");
          baseCat.push_back("zeectrllowcen");
        }
        else {
          int centEleId=getCentBinEleId((double)hiBin/2.);
          baseCat.push_back(TString("zeectrl")+(centEleId<1 ? "highcen" : "lowcen"));
        }
            
        float plotWgt(isMC ? weight : 1.0);
        ht.fill("mll", mll, plotWgt, baseCat);

        for(size_t iele=0; iele<2; iele++) {
          int eleIter=noIdEleIdx[iele];
          bool isEB(TMath::Abs(fForestEle.eleEta->at(eleIter))<barrelEndcapEta[0]);
          std::vector<TString> cat(baseCat);
          for(size_t icat=0; icat<cat.size(); icat++) cat[icat]+=(isEB?"eb":"ee");
          
          ht.fill("esihih",     fForestEle.eleSigmaIEtaIEta->at(eleIter),         plotWgt,cat);
          ht.fill("edetain",    TMath::Abs(fForestEle.eledEtaAtVtx->at(eleIter)), plotWgt,cat);
          ht.fill("ephiin",     TMath::Abs(fForestEle.eledPhiAtVtx->at(eleIter)), plotWgt,cat);
          ht.fill("ehoe",       fForestEle.eleHoverE->at(eleIter),                plotWgt,cat);
          ht.fill("ed0",        TMath::Abs(fForestEle.eleD0->at(eleIter)),        plotWgt,cat);
          ht.fill("edz",        TMath::Abs(fForestEle.eleDz->at(eleIter)),        plotWgt,cat);
          ht.fill("eeop",       TMath::Abs(fForestEle.eleEoverPInv->at(eleIter)), plotWgt,cat);         
          ht.fill("emisshits",  fForestEle.eleMissHits->at(eleIter),              plotWgt,cat);         
          ht.fill("echreliso",  fForestEle.elePFChIso03->at(eleIter)/fForestEle.elePt->at(eleIter),  plotWgt,cat);         
          ht.fill("ephoreliso", fForestEle.elePFPhoIso03->at(eleIter)/fForestEle.elePt->at(eleIter), plotWgt,cat);         
          ht.fill("eneureliso", fForestEle.elePFNeuIso03->at(eleIter)/fForestEle.elePt->at(eleIter), plotWgt,cat);         
        }
      }
    }
    
    int nLep=muIdx.size()+eleIdx.size();
    if(nLep<2) continue;
  
    std::vector<TLorentzVector> selLeps;
    bool isZee(false),isEM(false);
    if(eP4.size()>1) {
      TLorentzVector ee=eP4[0]+eP4[1];
      float charge=fForestEle.eleCharge->at(eleIdx[0])*fForestEle.eleCharge->at(eleIdx[1]);
      if(doSameSign) charge*=-1;
      if( fabs(ee.M()-91)<15 && charge<0 && eetrig>0 && fabs(eP4[0].Eta())<1.44 && fabs(eP4[1].Eta())<1.44) {
        isZee=true; 
        selLeps.push_back(eP4[0]);
        selLeps.push_back(eP4[1]);
      }
    }
    if(eP4.size() && muP4.size()) {
      TLorentzVector em=eP4[0]+muP4[0];
      float charge=fForestEle.eleCharge->at(eleIdx[0])*fForestMu.muCharge->at(muIdx[0]);
      if(doSameSign) charge*=-1;
      if(em.M()>20 && charge<0  && fabs(eP4[0].Eta())<1.44) {
        if(!isPP) {
          if(emtrig>0) isEM=true;
        }
        else
          if(emtrig>0 || eetrig>0) isEM=true;
      }

      if(isEM) {
        selLeps.push_back(eP4[0]);
        selLeps.push_back(muP4[0]);
      }
      
    }    
    if(!isZee && !isEM) continue;

    TString baseCat(isZee ? "zee" : "em");
    float plotWgt(isMC ? weight : 1.0);
    ht.fill( "lpt",     selLeps[0].Pt(),                       plotWgt, baseCat+"l1");
    ht.fill( "lpt",     selLeps[1].Pt(),                       plotWgt, baseCat+"l2");
    ht.fill( "leta",    fabs(selLeps[0].Eta()),                plotWgt, baseCat+"l1");
    ht.fill( "leta",    fabs(selLeps[1].Eta()),                plotWgt, baseCat+"l2");
    ht.fill( "dphill",  fabs(selLeps[0].DeltaPhi(selLeps[1])), plotWgt, baseCat);
    ht.fill( "mll",     (selLeps[0]+selLeps[1]).M(),           plotWgt, baseCat);
    ht.fill( "ptll",    (selLeps[0]+selLeps[1]).Pt(),          plotWgt, baseCat);

    std::vector<TLorentzVector> selJets;
    std::vector<int> selJetsIdx;
    int njets(0),nbjets(0);
    for(int jetIter = 0; jetIter < nref; jetIter++){
      if(jtpt[jetIter]<jetPtCut) continue;
      if(fabs(jteta[jetIter])>jetEtaCut) continue;
      TLorentzVector jp4(0,0,0,0);
      jp4.SetPtEtaPhiM(jtpt[jetIter],jteta[jetIter],jtphi[jetIter],jtm[jetIter]);
      if(jp4.DeltaR(selLeps[0])<0.4 || jp4.DeltaR(selLeps[1])<0.4) continue;      
      if(!isPP && trackN[jetIter]<2) continue;
      selJets.push_back(jp4);
      selJetsIdx.push_back(jetIter);
      ++njets;
      ht.fill( "jpt",   jtpt[jetIter],        plotWgt, baseCat);
      ht.fill( "jeta",  fabs(jteta[jetIter]), plotWgt, baseCat);
      ht.fill( "jcsv",  discr_csv[jetIter],   plotWgt, baseCat);
      if(discr_csv[jetIter]>csvWP) ++nbjets;
    }
    
    ht.fill( "njets",   njets,   plotWgt, baseCat);
    ht.fill( "nbjets",  nbjets,  plotWgt, baseCat);
    if(nbjets>=2 && !isPP && !isMC) {
      cout << "****************" << endl;
      cout << run << ":" << lumi << ":" << evt << endl;
      cout << "e: pt=" << selLeps[0].Pt() << " eta=" << selLeps[0].Eta() << " phi=" << selLeps[0].Phi() << endl;
      cout << "m: pt=" << selLeps[1].Pt() << " eta=" << selLeps[1].Eta() << " phi=" << selLeps[1].Phi() << endl;
      for(size_t ij=0; ij<selJetsIdx.size(); ij++)
        cout << "j" << ij << " pt=" << jtpt[ij] << " eta=" << jteta[ij] << " phi=" << jtphi[ij] << " csvV2=" << discr_csv[ij] << endl;
      cout << "****************" << endl;
    }
  }

  //save histos to file  
  if(outFileName!=""){
    TFile *fOut=TFile::Open(outFileName.c_str(),"RECREATE");
    fOut->cd();
    for (auto& it : ht.getPlots())  { 
      if(it.second->GetEntries()==0) continue;
      if(isMC && wgtSum!=0) it.second->Scale(1./wgtSum);
      it.second->SetDirectory(fOut); it.second->Write(); 
    }
    for (auto& it : ht.get2dPlots())  { 
      if(it.second->GetEntries()==0) continue;
      if(isMC && wgtSum!=0) it.second->Scale(1./wgtSum);
      it.second->SetDirectory(fOut); it.second->Write(); 
    }
    fOut->Close();
  }

  return;
}

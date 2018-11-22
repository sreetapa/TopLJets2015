void doSimpleScanEMu(int run=326392, int lumi=135, int evt=7205380,
                  TString jetType = "akPu4CaloJetAnalyzer",
                  TString selSuf = "&& nEle>0 && nMu>0",
                  TString printSuf = "",
                  TString strIn = "HiForestEMSkim.root") {

  //akPu4CaloJetAnalyzer
  //akPu4PFJetAnalyzer

  TFile *f = new TFile(strIn.Data());
  
  TTree *hiTree = dynamic_cast<TTree*>(f->Get("hiEvtAnalyzer/HiTree"));
  TTree *jetTree = dynamic_cast<TTree*>(f->Get(Form("%s/t",jetType.Data())));
  TTree *ggTree = dynamic_cast<TTree*>(f->Get("ggHiNtuplizerGED/EventTree"));
  
  hiTree->AddFriend(jetTree);
  hiTree->AddFriend(ggTree);

  hiTree->Scan(Form("run:lumi:evt:elePt:eleEta:elePhi:muPt:muEta:muPhi%s",printSuf.Data()),Form("run==%d && lumi==%d && evt==%d %s",run,lumi,evt,selSuf.Data()));

  
}

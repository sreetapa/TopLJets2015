void doSimpleScan(TString jetType = "akPu4CaloJetAnalyzer",
                  TString selSuf = "nEle>0 && nMu>0",
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

  //hiTree->Scan(Form("run:lumi:evt:elePt:eleEta:elePhi:muPt:muEta:muPhi:jtpt:jteta:jtphi:discr_csvV2%s",printSuf.Data()),Form("%s",selSuf.Data()));

  hiTree->Scan(Form("run:lumi:evt:jtpt:jteta:jtphi:discr_csvV2%s",printSuf.Data()),Form("%s",selSuf.Data()));
  
}

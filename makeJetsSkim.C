#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TChain.h"

#include <string>
#include <vector>
#include <iostream>

using namespace std;

//
float ptRel(TLorentzVector &jet, TLorentzVector &part)
{  
  float lj_x = jet.Px();
  float lj_y = jet.Py();
  float lj_z = jet.Pz();
    
  // absolute values squared
  float lj2  = lj_x*lj_x+lj_y*lj_y+lj_z*lj_z;
  float lep2 = part.X()*part.X()+part.Y()*part.Y()+part.Z()*part.Z();
  
  // projection vec(mu) to lepjet axis
  float lepXlj = part.X()*lj_x+part.Y()*lj_y+part.Z()*lj_z;
  
  // absolute value squared and normalized
  float pLrel2 = lepXlj*lepXlj/lj2;
  
  // lep2 = pTrel2 + pLrel2
  float pTrel2 = lep2-pLrel2;
  
  return (pTrel2 > 0) ? std::sqrt(pTrel2) : 0.0;
}


//
void makeJetsSkim(TString outFileName = "", TString inFileName = "", float jetPtCut = 25.,float jetEtaCut = 2.1, bool debug=false)
{
  if(inFileName=="" || outFileName=="")
    {
      cout << "Check inputs in=" << inFileName << " out=" << outFileName << endl;
      return;
    }
  cout << "[makeJetsSkim] in=" << inFileName << " -> out=" << outFileName 
       << " will be analysed for jets with pT>" << jetPtCut << " |eta|<" << jetEtaCut 
       << endl;

  //
  //prepare to read the trees
  //
  TFile *fIn=TFile::Open(inFileName);
  TTree *jetTree=(TTree *)fIn->Get("akCs3PFJetAnalyzer/t");
  int           nref;
  float         rawpt[60],jtpt[60],jteta[60],jtphi[60],jtm[60],jtarea[60];
  float         jttau1[60],jttau2[60],jttau3[60],discr_csvV1[60],discr_csvV2[60];
  int           svtxntrk[60];
  float         svtxdl[60],svtxdls[60],svtxdl2d[60],svtxdls2d[60],svtxm[60],svtxpt[60],svtxeta[60],svtxphi[60],svtxmcorr[60];
  float svtxTrkSumChi2[60];
  int svtxTrkNetCharge[60],svtxNtrkInCone[60];
  int          refparton_flavor[60],refparton_flavorForB[60];
  float  refpt[60];
  Float_t mue[60],mupt[60],mudr[60],muptrel[60];
  jetTree->SetBranchAddress("nref", &nref);
  jetTree->SetBranchAddress("rawpt", rawpt);
  jetTree->SetBranchAddress("jtpt", jtpt);
  jetTree->SetBranchAddress("jteta", jteta);
  jetTree->SetBranchAddress("jtphi", jtphi);
  jetTree->SetBranchAddress("jtm", jtm);
  jetTree->SetBranchAddress("jttau1", jttau1);
  jetTree->SetBranchAddress("jttau2", jttau2);
  jetTree->SetBranchAddress("jttau3", jttau3);
  jetTree->SetBranchAddress("discr_csvV1", discr_csvV1);
  jetTree->SetBranchAddress("discr_csvV2", discr_csvV2);
  jetTree->SetBranchAddress("svtxntrk", svtxntrk);
  jetTree->SetBranchAddress("svtxdl", svtxdl);
  jetTree->SetBranchAddress("svtxdls", svtxdls);
  jetTree->SetBranchAddress("svtxdl2d", svtxdl2d);
  jetTree->SetBranchAddress("svtxdls2d", svtxdls2d);
  jetTree->SetBranchAddress("svtxm", svtxm);
  jetTree->SetBranchAddress("svtxpt", svtxpt);
  jetTree->SetBranchAddress("svtxeta", svtxeta);
  jetTree->SetBranchAddress("svtxphi", svtxphi);
  jetTree->SetBranchAddress("svtxmcorr", svtxmcorr);
  jetTree->SetBranchAddress("svtxTrkSumChi2", svtxTrkSumChi2);
  jetTree->SetBranchAddress("svtxTrkNetCharge", svtxTrkNetCharge);
  jetTree->SetBranchAddress("svtxNtrkInCone", svtxNtrkInCone);
  jetTree->SetBranchAddress("refparton_flavor", refparton_flavor);
  jetTree->SetBranchAddress("refparton_flavorForB", refparton_flavorForB);
  jetTree->SetBranchAddress("refpt",refpt);
  jetTree->SetBranchAddress("mupt", mupt);
  jetTree->SetBranchAddress("mue", mue);
  jetTree->SetBranchAddress("mudr", mudr);
  jetTree->SetBranchAddress("muptrel", muptrel);
  
  TTree *tkTree = (TTree *) fIn->Get("anaTrack/trackTree");
  Int_t nTrk;
  Float_t  trkPt[5000],trkPtError[5000],trkEta[5000],trkPhi[5000];
  UChar_t  trkNHit[5000],trkNlayer[5000],trkNdof[5000];
  Float_t  trkDxy1[5000],trkDz1[5000],trkChi2[5000];
  Float_t  trkDxyOverDxyError[5000],trkDzOverDzError[5000];
  Bool_t   highPurity[5000],tight[5000],loose[5000];
  tkTree->SetBranchAddress("nTrk", &nTrk);
  tkTree->SetBranchAddress("trkPt", trkPt);
  tkTree->SetBranchAddress("trkPtError", trkPtError);
  tkTree->SetBranchAddress("trkNHit", trkNHit);
  tkTree->SetBranchAddress("trkNlayer", trkNlayer);
  tkTree->SetBranchAddress("trkEta", trkEta);
  tkTree->SetBranchAddress("trkPhi", trkPhi);
  tkTree->SetBranchAddress("trkDxyOverDxyError", trkDxyOverDxyError);
  tkTree->SetBranchAddress("trkDzOverDzError", trkDzOverDzError);
  tkTree->SetBranchAddress("highPurity", highPurity);
  tkTree->SetBranchAddress("tight", tight);
  tkTree->SetBranchAddress("loose", loose);
  tkTree->SetBranchAddress("trkChi2", trkChi2);
  tkTree->SetBranchAddress("trkNdof", trkNdof);
  tkTree->SetBranchAddress("trkDxy1", trkDxy1);
  tkTree->SetBranchAddress("trkDz1", trkDz1);

  TTree *hiTree = (TTree *) fIn->Get("hiEvtAnalyzer/HiTree");
  int hiBin;
  float vz;
  hiTree->SetBranchStatus("*", 0);
  hiTree->SetBranchStatus("hiBin", 1);
  hiTree->SetBranchStatus("vz", 1);
  hiTree->SetBranchAddress("hiBin", &hiBin);
  hiTree->SetBranchAddress("vz", &vz);
  
  //
  // prepare the output tree
  //
  TFile* outFile_p = new TFile(outFileName, "RECREATE");
  TTree *skimTree = new TTree("jets", "jets");
  skimTree->SetDirectory(outFile_p);
  std::map<TString,float> skimVars;
  TString varNames[]={"jet_pt","jet_eta","jet_phi","jet_m","jet_rawpt",
		      "gen_pt",
		      "jet_tau1","jet_tau2","jet_tau3",
		      "jet_csvV1","jet_csvV2",
		      "jet_svntrk","jet_svdl","jet_svdls","jet_svdl2d","jet_svdls2d",
		      "jet_svpt","jet_svm","jet_svmcorr", "jet_sve2e","jet_svptrel",
		      "jet_mue2jete",  "jet_mupt", "jet_muptrel", "jet_mudr",
		      "isB", "isC", "isUDSG", "isUnmatched"};
  for(auto s : varNames)
    {
      skimVars[s]=0.;
      skimTree->Branch(s,&(skimVars[s]),s+"/F");
    } 
  
  int n_Cpfcand;
  float Cpfcan_pt[5000],Cpfcan_dpt2pt[5000],Cpfcan_erel[5000],Cpfcan_phirel[5000],Cpfcan_etarel[5000],Cpfcan_deltaR[5000],Cpfcan_dxy[5000],Cpfcan_dxysig[5000],Cpfcan_dz[5000],Cpfcan_dzsig[5000],Cpfcan_chi2[5000],Cpfcan_nlayer[5000],Cpfcan_quality[5000],Cpfcan_ptrel[5000];
  skimTree->Branch("n_Cpfcand", &n_Cpfcand,"n_Cpfcand/i");
  skimTree->Branch("Cpfcan_pt", Cpfcan_pt,"Cpfcan_pt[n_Cpfcand]/f");
  skimTree->Branch("Cpfcan_dpt2pt", Cpfcan_dpt2pt,"Cpfcan_dpt2pt[n_Cpfcand]/f");
  skimTree->Branch("Cpfcan_erel", Cpfcan_erel,"Cpfcan_erel[n_Cpfcand]/f");
  skimTree->Branch("Cpfcan_phirel",Cpfcan_phirel,"Cpfcan_phirel[n_Cpfcand]/f");
  skimTree->Branch("Cpfcan_etarel",Cpfcan_etarel,"Cpfcan_etarel[n_Cpfcand]/f");
  skimTree->Branch("Cpfcan_deltaR",Cpfcan_deltaR,"Cpfcan_deltaR[n_Cpfcand]/f");
  skimTree->Branch("Cpfcan_dxy",Cpfcan_dxy,"Cpfcan_dxy[n_Cpfcand]/f");  
  skimTree->Branch("Cpfcan_dxysig",Cpfcan_dxysig,"Cpfcan_dxysig[n_Cpfcand]/f");
  skimTree->Branch("Cpfcan_dz",Cpfcan_dz,"Cpfcan_dz[n_Cpfcand]/f");
  skimTree->Branch("Cpfcan_dzsig",Cpfcan_dzsig,"Cpfcan_dzsig[n_Cpfcand]/f");
  skimTree->Branch("Cpfcan_chi2",Cpfcan_chi2,"Cpfcan_chi2[n_Cpfcand]/f");
  skimTree->Branch("Cpfcan_ptrel",Cpfcan_ptrel,"Cpfcan_ptrel[n_Cpfcand]/f");
  skimTree->Branch("Cpfcan_nlayer",Cpfcan_nlayer,"Cpfcan_nlayer[n_Cpfcand]/f");
  skimTree->Branch("Cpfcan_quality",Cpfcan_quality,"Cpfcan_quality[n_Cpfcand]/f");

  
  
  int nEntries = (int)jetTree->GetEntries();
  int entryDiv = ((int)(nEntries/20));
  for(int entry = 0; entry < nEntries; entry++)
    { 
      if(entry%entryDiv == 0 && nEntries >= 10000) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
      jetTree->GetEntry(entry);
      tkTree->GetEntry(entry);
      hiTree->GetEntry(entry);
      
      if(TMath::Abs(vz) > 15) continue;

      skimVars["hiBin"] = hiBin;
      
      for(int jetIter = 0; jetIter < nref; jetIter++){
        if(jtpt[jetIter]<jetPtCut) continue;
	if(fabs(jteta[jetIter])>jetEtaCut) continue;

	TLorentzVector jP4(0,0,0,0);
	jP4.SetPtEtaPhiM( jtpt[jetIter], jteta[jetIter], jtphi[jetIter], jtm[jetIter] );
	
	skimVars["jet_rawpt"]   = rawpt[jetIter];
	skimVars["jet_pt"]      = jtpt[jetIter];
	skimVars["jet_eta"]     = jteta[jetIter];
	skimVars["jet_phi"]     = jtphi[jetIter];
	skimVars["jet_m"]       = jtm[jetIter];
	skimVars["jet_tau1"]    = jttau1[jetIter];
	skimVars["jet_tau2"]    = jttau2[jetIter];
	skimVars["jet_tau3"]    = jttau3[jetIter];
	skimVars["jet_csvV1"]   = discr_csvV1[jetIter];
	skimVars["jet_csvV2"]   = discr_csvV2[jetIter];
	skimVars["jet_svntrk"]  = svtxntrk[jetIter];
	skimVars["jet_svdl"]    = svtxdl[jetIter];
	skimVars["jet_svdls"]   = svtxdls[jetIter];
	skimVars["jet_svdl2d"]  = svtxdl2d[jetIter];
	skimVars["jet_svdls2d"] = svtxdls2d[jetIter];
	skimVars["jet_svpt"]    = svtxpt[jetIter];	
	TLorentzVector svtxP4(0,0,0,0);
	svtxP4.SetPtEtaPhiM(svtxpt[jetIter],svtxeta[jetIter],svtxphi[jetIter],svtxm[jetIter]);
	skimVars["jet_sve2e"] = svtxP4.E()/jP4.E();
	skimVars["jet_svptrel"] = ptRel(jP4,svtxP4);
	skimVars["jet_svm"]     = svtxm[jetIter];
	skimVars["jet_svmcorr"] = svtxmcorr[jetIter];
	skimVars["jet_svtxTrkNetCharge"] = svtxTrkNetCharge[jetIter];
	skimVars["jet_svtrkSumChi2"] = svtxTrkSumChi2[jetIter];
	skimVars["jet_svtrkInCone"] = svtxNtrkInCone[jetIter];
	skimVars["jet_mue2jete"]    = mue[jetIter]/jP4.E();
	skimVars["jet_mupt"]        = mupt[jetIter] ;
	skimVars["jet_muptrel"]     = muptrel[jetIter];
	skimVars["jet_mudr"]        = mudr[jetIter];
	int absid( abs(refparton_flavorForB[jetIter]) );
	skimVars["isB"]             = absid==5 ? 1. : 0.;
	skimVars["isC"]             = absid==4 ? 1. : 0.;
	skimVars["isUnmatched"]     = absid==0 ? 1. : 0.;
	skimVars["isUDSG"]          = (absid!=5 && absid!=4 && absid!=0);
	skimVars["gen_pt"]          = refpt[jetIter];

	float etasign = jteta[jetIter] <0 ? -1. : 1.;
	
	n_Cpfcand=0;
	for(int ipf=0; ipf<nTrk; ipf++)
	  {
	    float etarel=etasign*fabs(trkEta[ipf]-jteta[jetIter]);
	    float phirel=TVector2::Phi_mpi_pi(trkPhi[ipf]-jtphi[jetIter]);
	    float dR=sqrt(pow(etarel,2)+pow(phirel,2));
	    if(dR>0.3) continue;
	    

	    TLorentzVector tkP4(0,0,0,0);
	    tkP4.SetPtEtaPhiM( trkPt[ipf], trkEta[ipf], trkPhi[ipf], 0.139);
	    Cpfcan_pt[n_Cpfcand]=tkP4.Pt();
	    Cpfcan_dpt2pt[n_Cpfcand]=trkPtError[ipf]/tkP4.Pt();
	    Cpfcan_erel[n_Cpfcand]=tkP4.E()/jP4.E();
	    Cpfcan_phirel[n_Cpfcand]=phirel;
	    Cpfcan_etarel[n_Cpfcand]=etarel;
	    Cpfcan_ptrel[n_Cpfcand]=ptRel(jP4,tkP4);
	    Cpfcan_deltaR[n_Cpfcand]=dR;
	    Cpfcan_dxy[n_Cpfcand]=trkDxy1[ipf];
	    Cpfcan_dxysig[n_Cpfcand]=trkDxyOverDxyError[ipf];
	    Cpfcan_dz[n_Cpfcand]=trkDz1[ipf];					\
	    Cpfcan_dzsig[n_Cpfcand]=trkDzOverDzError[ipf];
	    Cpfcan_nlayer[n_Cpfcand]=trkNlayer[ipf];
	    Cpfcan_chi2[n_Cpfcand]=trkChi2[ipf]/trkNdof[ipf];
	    Cpfcan_quality[n_Cpfcand]=float(loose[ipf]+tight[ipf]*10+highPurity[ipf]*100);

	    n_Cpfcand++;
	  }

	skimTree->Fill();
      }
    }

  fIn->Close();
  outFile_p->cd();
  skimTree->SetDirectory(outFile_p);
  skimTree->Write();
  outFile_p->Close();

  cout << "[makeJetsSkim] all done" << endl;
}

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
float catchInfs(const float in,const float replace_value){
  if(in==in){
    if(std::isinf(in))
      return replace_value;
    else if(in < -1e32 || in > 1e32)
      return replace_value;
    return in;
  }
  return replace_value;
}

//
float catchInfsAndBound(const float in,const float replace_value, const float lowerbound, const float upperbound){
  float withoutinfs=catchInfs(in,replace_value);
  if(withoutinfs<lowerbound) return lowerbound;
  if(withoutinfs>upperbound) return upperbound;
  return withoutinfs;
}


//
class Jet
{
public:
  Jet() {}
  void reset() 
  {
    hibin=0;
    pt=0; eta=0;phi=0;m=0;rawpt=0;genpt=0;tau1=0;tau2=0;tau3=0;
    csvV1=0; csvV2=0;trackSumJetDeltaR=0;trackSip2dSigAboveCharm=0;trackSip3dSigAboveCharm=0;trackSip2dValAboveCharm=0;trackSip3dValAboveCharm=0;
    svntrk=0;svdl=0;svdls=0;svdl2d=0;svdls2d=0;svpt=0;svm=0;svmcorr=0;sve2e=0;svptrel=0;svdr2jet=0;
    mupt2jet=0;mupt=0;muptrel=0;mudr=0;
    ntk=0;
    isB=0;isC=0;isUDSG=0;isUnmatched=0;
  }

  float ptRel(TLorentzVector &part)
  {  
    //jet kinematics
    TLorentzVector jet(0,0,0,0);
    jet.SetPtEtaPhiM(pt,eta,phi,m);
    float lj_x = jet.Px();
    float lj_y = jet.Py();
    float lj_z = jet.Pz();
    
    // absolute values squared
    float lj2  = jet.Vect().Mag2();
    float lep2 = part.Vect().Mag2();
    
    // projection vec(mu) to lepjet axis
    float lepXlj = jet.Vect().Dot(part.Vect());
    
    // absolute value squared and normalized
    float pLrel2 = lepXlj*lepXlj/lj2;
    
    // lep2 = pTrel2 + pLrel2
    float pTrel2 = lep2-pLrel2;
    
    return (pTrel2 > 0) ? std::sqrt(pTrel2) : 0.0;
  }
  
  void instantiateTree(TTree *t)
  {
    t->Branch("hibin",&hibin,"hibin/I"); 
    t->Branch("event",&event,"event/l"); 
    t->Branch("run",&run,"run/i"); 
    t->Branch("lumi",&lumi,"lumi/i"); 
    t->Branch("jetIdx",&jetIdx,"jetIdx/I"); 
    t->Branch("pt",&pt,"pt/F"); 
    t->Branch("eta",&eta,"eta/F");
    t->Branch("phi",&phi,"phi/F");
    t->Branch("m",&m,"m/F");
    t->Branch("rawpt",&rawpt,"rawpt/F");
    t->Branch("genpt",&genpt,"genpt/F");    
    t->Branch("tau1",&tau1,"tau1/F");
    t->Branch("tau2",&tau2,"tau2/F");
    t->Branch("tau3",&tau3,"tau3/F");
    t->Branch("csvV1",&csvV1,"csvV1/F"); 
    t->Branch("csvV2",&csvV2,"csvV2/F");
    t->Branch("trackSumJetDeltaR",&trackSumJetDeltaR,"trackSumJetDeltaR/F");
    t->Branch("trackSip2dSigAboveCharm",&trackSip2dSigAboveCharm,"trackSip2dSigAboveCharm/F");
    t->Branch("trackSip3dSigAboveCharm",&trackSip3dSigAboveCharm,"trackSip3dSigAboveCharm/F");
    t->Branch("trackSip2dValAboveCharm",&trackSip2dValAboveCharm,"trackSip2dValAboveCharm/F");
    t->Branch("trackSip3dValAboveCharm",&trackSip3dValAboveCharm,"trackSip3dValAboveCharm/F");
    t->Branch("svntrk",&svntrk,"svntrk/I");
    t->Branch("svdl",&svdl,"svdl/F");
    t->Branch("svdls",&svdls,"svdls/F");
    t->Branch("svdl2d",&svdl2d,"svdl2d/F");
    t->Branch("svdls2d",&svdls2d,"svdls2d/F");
    t->Branch("svpt",&svpt,"svpt/F");
    t->Branch("svm",&svm,"svm/F");
    t->Branch("svmcorr",&svmcorr,"svmcorr/F");
    t->Branch("sve2e",&sve2e,"sve2e/F");
    t->Branch("svdr2jet",&svdr2jet,"svdr2jet/F");
    t->Branch("svptrel",&svptrel,"svptrel/F");
    t->Branch("svtksumchi2",&svtksumchi2,"svtksumchi2/F");
    t->Branch("svcharge",&svcharge,"svcharge/I");
    t->Branch("svtkincone",&svtkincone,"svtkincone/I");
    t->Branch("mupt2jet",&mupt2jet,"mupt2jet/F");
    t->Branch("mupt",&mupt,"mupt/F");
    t->Branch("muptrel",&muptrel,"muptrel/F");
    t->Branch("mudr",&mudr,"mudr/F");
    t->Branch("isB",&isB,"isB/I");
    t->Branch("isC",&isC,"isC/I");
    t->Branch("isUDSG",&isUDSG,"isUDSG/I");
    t->Branch("isUnmatched",&isUnmatched,"isUnmatched/I");
    t->Branch("ntk",&ntk,"ntk/I");
    t->Branch("tkptrel",tkptrel,"tkptrel/F");
    t->Branch("tkdr",tkdr,"tkdr/F");
    t->Branch("tkprob0",tkprob0,"tkprob0/F");
    t->Branch("tkprob1",tkprob1,"tkprob1/F");
    t->Branch("tkip2d",tkip2d,"tkip2d/F");
    t->Branch("tkip2dsig",tkip2dsig,"tkip2dsig/F");
    t->Branch("tkip3d",tkip3d,"tkip3d/F");
    t->Branch("tkip3dsig",tkip3dsig,"tkip3dsig/F");
    t->Branch("tkipdist2j",tkipdist2j,"tkipdist2j/F");
    t->Branch("tkipclosest2j",tkipclosest2j,"tkipclosest2j/F");
    t->Branch("tkptratio",tkptratio,"tkptratio/F");
    t->Branch("tkpparratio",tkpparratio,"tkpparratio/F");
  }

  //variables
  unsigned int run,lumi;
  unsigned long long event;
  int hibin,jetIdx;
  float pt, eta,phi,m,rawpt,genpt,tau1,tau2,tau3;
  float csvV1, csvV2;
  float trackSumJetDeltaR,trackSip2dSigAboveCharm,trackSip3dSigAboveCharm,trackSip2dValAboveCharm,trackSip3dValAboveCharm;
  float svdl,svdls,svdl2d,svdls2d,svpt,svm,svmcorr,sve2e,svptrel,svtksumchi2,svdr2jet;
  int svntrk,svcharge,svtkincone;
  float mupt2jet,mupt,muptrel,mudr;
  int ntk;
  float tkptrel[200],tkdr[200],tkdxy[200],tkdz[200],tkchi2[200],tkprob0[200],tkprob1[200],tkip2d[200],tkip2dsig[200],tkip3d[200],tkip3dsig[200],tkipdist2j[200],tkipclosest2j[200],tkptratio[200],tkpparratio[200];
  int tknhit[200],tknhitpixel[200],tknhitstrip[200];
  int isB,isC,isUDSG,isUnmatched;  
};



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
  float         rawpt[100],jtpt[100],jteta[100],jtphi[100],jtm[100],jtarea[100];
  float         jttau1[100],jttau2[100],jttau3[100],discr_csvV1[100],discr_csvV2[100];
  int           svtxntrk[100];
  float         svtxdl[100],svtxdls[100],svtxdl2d[100],svtxdls2d[100],svtxm[100],svtxpt[100],svtxeta[100],svtxphi[100],svtxmcorr[100];
  float svtxTrkSumChi2[100],svJetDeltaR[100];
  int svtxTrkNetCharge[100],svtxNtrkInCone[100];
  Float_t         trackSip2dSigAboveCharm[100],trackSip3dSigAboveCharm[100],trackSip2dValAboveCharm[100],trackSip3dValAboveCharm[100],trackSumJetDeltaR[100];
  int          refparton_flavor[100],refparton_flavorForB[100];
  Int_t           nIP;
  Int_t           ipJetIndex[4000];
  Float_t ipPt[4000],ipProb0[4000],ipProb1[4000],ip2d[4000],ip2dSig[4000],ip3d[4000],ip3dSig[4000],ipDist2Jet[4000],ipDist2JetSig[4000],ipClosest2Jet[4000],trackPtRel[4000],trackPPar[4000],trackPParRatio[4000],trackDeltaR[4000],trackPtRatio[4000];
  float  refpt[100];
  Float_t mupt[100],mudr[100],muptrel[100];
  
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
  jetTree->SetBranchAddress("trackSip2dSigAboveCharm", trackSip2dSigAboveCharm);
  jetTree->SetBranchAddress("trackSip3dSigAboveCharm", trackSip3dSigAboveCharm);
  jetTree->SetBranchAddress("trackSip2dValAboveCharm", trackSip2dValAboveCharm);
  jetTree->SetBranchAddress("trackSip3dValAboveCharm", trackSip3dValAboveCharm);
  jetTree->SetBranchAddress("trackSumJetDeltaR",trackSumJetDeltaR);
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
  jetTree->SetBranchAddress("svJetDeltaR", svJetDeltaR);
  jetTree->SetBranchAddress("refparton_flavor", refparton_flavor);
  jetTree->SetBranchAddress("refparton_flavorForB", refparton_flavorForB);
  jetTree->SetBranchAddress("refpt",refpt);
  jetTree->SetBranchAddress("mupt", mupt);
  jetTree->SetBranchAddress("mudr", mudr);
  jetTree->SetBranchAddress("muptrel", muptrel);
  jetTree->SetBranchAddress("nIP", &nIP);
  jetTree->SetBranchAddress("ipJetIndex", ipJetIndex);
  jetTree->SetBranchAddress("ipPt", ipPt);
  jetTree->SetBranchAddress("ipProb0", ipProb0);
  jetTree->SetBranchAddress("ipProb1", ipProb1);
  jetTree->SetBranchAddress("ip2d", ip2d);
  jetTree->SetBranchAddress("ip2dSig", ip2dSig);
  jetTree->SetBranchAddress("ip3d", ip3d);
  jetTree->SetBranchAddress("ip3dSig", ip3dSig);
  jetTree->SetBranchAddress("ipDist2Jet", ipDist2Jet);
  jetTree->SetBranchAddress("ipDist2JetSig", ipDist2JetSig);
  jetTree->SetBranchAddress("ipClosest2Jet", ipClosest2Jet);
  jetTree->SetBranchAddress("trackPtRel", trackPtRel);
  jetTree->SetBranchAddress("trackPPar", trackPPar);
  jetTree->SetBranchAddress("trackPParRatio", trackPParRatio);
  jetTree->SetBranchAddress("trackDeltaR", trackDeltaR);
  jetTree->SetBranchAddress("trackPtRatio", trackPtRatio);
  
  TTree *hiTree = (TTree *) fIn->Get("hiEvtAnalyzer/HiTree");
  int hiBin;
  float vz;
  unsigned int run,lumi;
  unsigned long long event;
  hiTree->SetBranchStatus("*", 0);
  hiTree->SetBranchStatus("hiBin", 1);
  hiTree->SetBranchStatus("vz", 1);
  hiTree->SetBranchStatus("run", 1);
  hiTree->SetBranchStatus("lumi", 1);
  hiTree->SetBranchStatus("evt", 1);
  hiTree->SetBranchAddress("hiBin", &hiBin);
  hiTree->SetBranchAddress("vz", &vz);
  hiTree->SetBranchAddress("run", &run);
  hiTree->SetBranchAddress("lumi", &lumi);
  hiTree->SetBranchAddress("evt", &event);

  //
  // prepare the output tree
  //
  TFile* outFile_p = new TFile(outFileName, "RECREATE");
  TTree *skimTree = new TTree("jets", "jets");
  skimTree->SetDirectory(outFile_p);
  
  Jet j;
  j.instantiateTree(skimTree);
  
  //loop over the events
  int nEntries = (int)jetTree->GetEntries();
  int entryDiv = ((int)(nEntries/20));
  for(int entry = 0; entry < nEntries; entry++)
    { 
      if(entry%entryDiv == 0 && nEntries >= 10000) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
      jetTree->GetEntry(entry);
      hiTree->GetEntry(entry);
      
      if(TMath::Abs(vz) > 15) continue;
      
      for(int jetIter = 0; jetIter < nref; jetIter++){
        if(jtpt[jetIter]<jetPtCut) continue;
	if(fabs(jteta[jetIter])>jetEtaCut) continue;
	
	TLorentzVector jP4(0,0,0,0);
	jP4.SetPtEtaPhiM( jtpt[jetIter], jteta[jetIter], jtphi[jetIter], jtm[jetIter] );


	j.event=event;
	j.run=run;
	j.lumi=lumi;
	j.hibin=hiBin;

	j.jetIdx=jetIter;
	j.pt=jtpt[jetIter];
	j.eta=jteta[jetIter];
	j.phi=jtphi[jetIter];
	j.m=jtm[jetIter];
	j.genpt=catchInfsAndBound(refpt[jetIter],-1,-1,999999.);
	j.rawpt=rawpt[jetIter];
	j.tau1=jttau1[jetIter];
	j.tau2=jttau2[jetIter];
	j.tau3=jttau3[jetIter];

	int absid( abs(refparton_flavorForB[jetIter]) );
	j.isB         = absid==5 ? 1. : 0.;
	j.isC         = absid==4 ? 1. : 0.;
	j.isUnmatched = absid==0 ? 1. : 0.;
	j.isUDSG      = (absid!=5 && absid!=4 && absid!=0);
 
	j.csvV1=catchInfsAndBound(discr_csvV1[jetIter],-2,-2,1);
	j.csvV2=catchInfsAndBound(discr_csvV2[jetIter],-2,-2,1);
	j.trackSumJetDeltaR=trackSumJetDeltaR[jetIter];
	j.trackSip2dSigAboveCharm=catchInfsAndBound(trackSip2dSigAboveCharm[jetIter],-50,-40,40);
	j.trackSip3dSigAboveCharm=catchInfsAndBound(trackSip3dSigAboveCharm[jetIter],-50,-40,40);
	j.trackSip2dValAboveCharm=catchInfsAndBound(trackSip2dValAboveCharm[jetIter],-0.2,-0.1,0.2);
	j.trackSip3dValAboveCharm=catchInfsAndBound(trackSip3dValAboveCharm[jetIter],-0.2,-0.1,0.2);

	TLorentzVector svtxP4(0,0,0,0);
	svtxP4.SetPtEtaPhiM(svtxpt[jetIter],svtxeta[jetIter],svtxphi[jetIter],svtxm[jetIter]);	
	j.svntrk=svtxntrk[jetIter];
	j.svdl=catchInfsAndBound(svtxdl[jetIter],-1,0,6);
	j.svdls=catchInfsAndBound(svtxdls[jetIter],-1,0,100);
	j.svdl2d=catchInfsAndBound(svtxdl2d[jetIter],-1,0,6);
	j.svdls2d=catchInfsAndBound(svtxdls2d[jetIter],-1,0,100);
	j.svpt=svtxpt[jetIter];
	j.svm=svtxm[jetIter];
	j.svmcorr=svtxmcorr[jetIter];
	j.sve2e=catchInfsAndBound(svtxP4.E()/jP4.E(),-1,0,2);
	j.svptrel=catchInfsAndBound(j.ptRel(svtxP4),-1,-1,8);
	j.svcharge=catchInfsAndBound(svtxTrkNetCharge[jetIter],0,-15,15);
	j.svtksumchi2=catchInfsAndBound(svtxTrkSumChi2[jetIter],-1,0,200);
	j.svtkincone=svtxNtrkInCone[jetIter];
	j.svdr2jet=svJetDeltaR[jetIter];

	j.mupt=catchInfsAndBound(mupt[jetIter],-1,0.,99999.);
	j.mupt2jet=j.mupt>0 ? catchInfsAndBound(mupt[jetIter]/jP4.Pt(),-1,0.,2.) : 0.;
	j.muptrel=j.mupt>0 ? catchInfsAndBound(muptrel[jetIter],-1,0.,8.) : 0.;
	j.mudr=j.mupt>0 ? catchInfsAndBound(mudr[jetIter],-1,0.,1.0) : 0.;
  
	j.ntk=0;
	for(int ip=0; ip<nIP; ip++)
	  {
	    if(ipJetIndex[ip]!=jetIter) continue;
	    j.tkptrel[j.ntk]       = catchInfsAndBound(trackPtRel[ip],-1,0,8);
	    j.tkdr[j.ntk]          = catchInfsAndBound(trackDeltaR[ip],-1,0,1);	    
	    j.tkprob0[j.ntk]       = ipProb0[ip];
	    j.tkprob1[j.ntk]       = ipProb1[ip];
	    j.tkip2d[j.ntk]        = catchInfsAndBound(ip2d[ip],0,-0.15,0.15);
	    j.tkip2dsig[j.ntk]     = catchInfsAndBound(ip2dSig[ip],0,-20,20);
	    j.tkip3d[j.ntk]        = catchInfsAndBound(ip3d[ip],0,-0.15,0.15);
	    j.tkip3dsig[j.ntk]     = catchInfsAndBound(ip3dSig[ip],0,-20,20);
	    j.tkipdist2j[j.ntk]    = catchInfsAndBound(ipDist2Jet[ip],0,-0.2,0.2);
	    j.tkipclosest2j[j.ntk] = catchInfsAndBound(ipClosest2Jet[ip],0,0,5);
	    j.tkptratio[j.ntk]     = catchInfsAndBound(trackPtRatio[ip],0,0,1);
	    j.tkpparratio[j.ntk]   = catchInfsAndBound(trackPParRatio[ip],0,0,1);
	    j.ntk++;
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

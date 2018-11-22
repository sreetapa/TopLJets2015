#ifndef electronidcriteria_h
#define electronidcriteria_h

//PbPb: https://twiki.cern.ch/twiki/bin/view/CMS/ElectronPbPb5TeV#3_Selection_for_different_centra
//pp:   https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns


const int nBarrelEndcap=2;
const int nCentEleId=4;
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

void setEleIdCuts() {
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
}

//For PbPb electron ID is centrality dependent
//centrality bins 0-20%, 20-40%, 40-70%, 70-100%
double centMinEleId[4] = {0.,20.,40.,70.};
double centMaxEleId[4] = {20.,40.,70.,100.};

int getCentBinEleId(double cent) {

  int centBin = -1;
  for(unsigned int i = 0; i<sizeof(centMinEleId)/sizeof(double); ++i) {
    if(cent>=centMinEleId[i] && cent<centMaxEleId[i]) centBin = i;
  }
  return centBin;
}


#endif

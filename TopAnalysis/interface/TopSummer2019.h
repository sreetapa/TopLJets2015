#ifndef _TopSummer2019_h_
#define _TopSummer2019_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

void RunTopSummer2019(const TString filename,
                      TString outname,
                      TH1F *normH, 
                      TH1F *genPU,
                      TString era,
                      Bool_t debug=false);

#endif

#/bin/bash
#root -b -q simpleTT2L.C+\(\"Data2018_EM_SS.root\",\"/eos/cms/store/cmst3/user/psilva/PbPb2018/Data2018PbPb_EM.root\",false,false,true\)
#root -b -q simpleTT2L.C+\(\"Data2018_ZEE_SS.root\",\"/eos/cms/store/cmst3/user/psilva/PbPb2018/Data2018PbPb_ZEE.root\",false,false,true\)
#root -b -q simpleTT2L.C+\(\"Data2018_EM.root\",\"/eos/cms/store/cmst3/user/psilva/PbPb2018/Data2018PbPb_EM.root\",false,false\)
#root -b -q simpleTT2L.C+\(\"Data2018_ZEE.root\",\"/eos/cms/store/cmst3/user/psilva/PbPb2018/Data2018PbPb_ZEE.root\",false,false\)
#root -b -q simpleTT2L.C+\(\"MCpp5p2TeV_DY.root\",\"/eos/cms/store/group/phys_top/gkrintir/TopHI/RunIIpp5Spring18DR/DYJetsToLL_MLL-50_TuneCP5_5020GeV-amcatnloFXFX-pythia8/merged/HiForestAODSIM_94.root\",true,true\)
root -b -q runTTto2Lselection.C+\(\"MCpp5p2TeV_TT.root\",\"/eos/cms/store/group/phys_top/gkrintir/TopHI/RunIIpp5Spring18DR/TT_TuneCP5_5p02TeV-powheg-pythia8/merged/HiForestAODSIM_94.root\",true,true\)

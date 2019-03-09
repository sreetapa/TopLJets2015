#!/bin/bash

#determine CMSSW config
SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`
ARCH=${SCRIPTPATH##/*/}
WORKDIR=${SCRIPTPATH}/..

#configure environment
cd $WORKDIR
export SCRAM_ARCH=$ARCH
eval `scram r -sh`

input=${1}
output=${2}
isMC=${3}
isPP=${4}
doSameSign=${5}
input=/eos/cms/store/group/phys_top/PbPbTTbar_2018/SkimElectrons_PromptRecov1/Chunk_55_ext0.root
output=plots
isMC=false
isPP=false
doSameSign=false

root -b -l << EOF
gSystem->Load("/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/fastjet/3.3.0-omkpbe/lib/libfastjet.so");
gSystem->Load("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/test/runTTto2Lselection_C.so");
runTTto2Lselection("${input}","${output}",${isMC},${isPP},${doSameSign});
.q
EOF

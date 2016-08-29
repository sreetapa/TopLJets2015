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

cone=${1}
isPbPb=${2}
evtsPerJob=${3}
startEvt=${4}
outFile=${5}
vecList=${6}
root -b -q .L makeMuJetsSkim.C+\(\"test\",${cone},${isPbPb},${evtsPerJob},${startEvt},false,false,\"${outFile}\",\"${vecList}\"\);
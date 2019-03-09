#!/bin/bash

HOME=`pwd`

#determine CMSSW config
SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`
ARCH=${SCRIPTPATH##/*/}
WORKDIR=${SCRIPTPATH}/..

#configure environment
cd $WORKDIR
export SCRAM_ARCH=$ARCH
eval `scram r -sh`

cd $HOME
input=${1}
output=${2}

#FIXME 
extraOpts=""
isMC=${3}
isPP=${4}
doSameSign=${5}

runTTto2Lselection --in ${input} --out tto2l.root ${extraOpts}
cp -v tto2l.root ${output}


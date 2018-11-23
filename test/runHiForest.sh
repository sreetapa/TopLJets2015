#!/bin/bash
cd /afs/cern.ch/user/p/psilva/work/PbPb/CMSSW_10_3_1/src/HeavyIonsAnalysis/
eval `scram r -sh`
cd -
cfg=${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/test/runForestAOD_pponAA_DATA_103X_PR.py
cmsRun ${cfg} inputFiles=${1} outputFile=${2} maxEvents=-1
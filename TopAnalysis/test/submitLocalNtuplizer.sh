#!/bin/bash

outdir=/store/cmst3/group/top/RunIIReReco/2017/newproton_calib

#low pileup runs
json=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_306896-307082_13TeV_PromptReco_#Collisions17_JSON_LowPU.txt
commonOpts="--addParent --proxy --lumiMask ${json} --output ${outdir}"
commonOpts="${commonOpts} --extraOpts runOnData=True,era=era2017"
jobList=(
    "--jobTag Data13TeV_2017H_SingleMuon_v2 --dataset /SingleMuon/Run2017H-17Nov2017-v2/MINIAOD"
    #"--jobTag Data13TeV_2017H_DoubleMuon    --dataset /HighEGJet/Run2017H-17Nov2017-v1/MINIAOD"  #only in T0/T1
    #"--jobTag Data13TeV_2017H_HighEGJet     --dataset /HighEGJet/Run2017H-17Nov2017-v1/MINIAOD"  #only in T0/T1
)

for ((i = 0; i < ${#jobList[@]}; i++))
do
    echo ${jobList[$i]}
    python scripts/submitLocalNtuplizer.py ${commonOpts} ${jobList[$i]}; 
done


#high pileup runs
#json=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt
#commonOpts="--addParent --proxy --lumiMask ${json} --output ${outdir}"
#commonOpts="${commonOpts} --extraOpts runOnData=True,era=era2017"
#jobList=(
#    "--jobTag Data13TeV_2017B_SingleMuon --dataset /SingleMuon/Run2017B-31Mar2018-v1/MINIAOD"
#    "--jobTag Data13TeV_2017C_SingleMuon --dataset /SingleMuon/Run2017C-31Mar2018-v1/MINIAOD"
#    "--jobTag Data13TeV_2017D_SingleMuon --dataset /SingleMuon/Run2017D-31Mar2018-v1/MINIAOD"
#    "--jobTag Data13TeV_2017E_SingleMuon --dataset /SingleMuon/Run2017E-31Mar2018-v1/MINIAOD"
#    "--jobTag Data13TeV_2017F_SingleMuon --dataset /SingleMuon/Run2017F-31Mar2018-v1/MINIAOD"
#    "--jobTag Data13TeV_2017B_DoubleMuon --dataset /DoubleMuon/Run2017B-31Mar2018-v1/MINIAOD"
#    "--jobTag Data13TeV_2017C_DoubleMuon --dataset /DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD"
#    "--jobTag Data13TeV_2017D_DoubleMuon --dataset /DoubleMuon/Run2017D-31Mar2018-v1/MINIAOD"
#    "--jobTag Data13TeV_2017E_DoubleMuon --dataset /DoubleMuon/Run2017E-31Mar2018-v1/MINIAOD"
#    "--jobTag Data13TeV_2017F_DoubleMuon --dataset /DoubleMuon/Run2017F-31Mar2018-v1/MINIAOD"
#)
#for ((i = 0; i < ${#jobList[@]}; i++))
#do
#    echo ${jobList[$i]}
#    python scripts/submitLocalNtuplizer.py ${commonOpts} ${jobList[$i]}; 
#done


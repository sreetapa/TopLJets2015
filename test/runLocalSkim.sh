WORKDIR=`pwd`

# example
# sh /afs/cern.ch/user/p/psilva/scratch0/CMSSW_7_5_8_patch3/src/topskim/test/runLocalSkim.sh /store/himc/HINPbPbWinter16DR/Pythia6_bJet170_pp502_Hydjet_MB/AODSIM/75X_mcRun2_HeavyIon_v13-v1/50000/00A14813-B70E-E611-8C4E-02163E014702.root hiftest.root /tmp/psilva

# TODOs 
# start by setting the proxy env
# move to storage area
# make mother script which loops over input files and sends jobs to batch
#

echo "Working directory is ${WORKDIR}"
cd /afs/cern.ch/user/p/psilva/scratch0/CMSSW_7_5_8_patch3/src/topskim
eval `scram r -sh`
cd ${WORKDIR}
input=${1}
outfile=${2}
outdir=${3}
echo "Forestizing ${input}"
cmsRun ${CMSSW_BASE}/src/topskim/test/runOverPbPb_MIX_75X.py inputFile=${input} outFilename=hif_${outfile}
echo "Making jet tree out of HIForest"
root -l ${CMSSW_BASE}/src/topskim/makeJetsSkim.C+\(\"${outfile}\",\"hif_${outfile}\"\) && .q;
echo "Moving result to ${outdir}"
rm hif_${outfile}
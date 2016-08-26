import ROOT
import os,sys
from subprocess import Popen, PIPE

if len(sys.argv)<3:
    print "python scripts/launchSkimFor.py directory evtsPerJob tag outDir"
    exit(-1)

baseDir    = 'eos/cms'+sys.argv[1]
evtsPerJob = int(sys.argv[2])
tag        = 'test' if len(sys.argv)<4 else sys.argv[3]
outDir     = 'results' if len(sys.argv)<5 else sys.argv[4]
cone       = 0.2
isPbPb     = 1
queue      = '8nh'
cwd        = os.getcwd()

#mount locally EOS                                                                                                                                                                                             
eos_cmd = '/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'
Popen([eos_cmd, ' -b fuse mount', 'eos'],stdout=PIPE).communicate()

#get
vecList=''
chain=ROOT.TChain('hiEvtAnalyzer/HiTree')
for url in os.listdir(baseDir):
    if '.root' in url : 
        chain.AddFile(baseDir+'/'+url)
        vecList+='root://eoscms//'+baseDir+'/'+url+','
nEvts=chain.GetEntries()
vecList=vecList[:-1]

#unmount locally EOS                                                                                                                                                                                           
Popen([eos_cmd, ' -b fuse umount', 'eos'],stdout=PIPE).communicate()

nJobs=nEvts/evtsPerJob
if nJobs*nEvts<evtsPerJob : nJobs=nJobs+1
print '#events=',nEvts,'will be analysed in',nEvts/evtsPerJob,'jobs'

os.system('mkdir -p '+outDir)
for i in xrange(0,nJobs):
    startEvt=i*evtsPerJob
    outFile='%s/%s/%s_%d.root'%(cwd,outDir,tag,startEvt)
    localRun='root -b -q .L makeMuJetsSkim.C+(\\\"test\\\",%f,%d,%d,%d,false,false,\\\"%s\\\",\\\"%s\\\")'%(cone,isPbPb,evtsPerJob,startEvt,outFile,vecList)
    cmd='bsub -q %s %s/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (queue,cwd,localRun)
    os.system(cmd)
    raw_input()

print 'The output of the jobs will be named as %s_i.root and can be found in %s'%(tag,outDir)

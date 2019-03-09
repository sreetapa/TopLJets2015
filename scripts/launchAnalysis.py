import ROOT
import os,sys
from subprocess import Popen, PIPE

if len(sys.argv)<3:
    print "python scripts/launchAnalysis.py input output [queue=local/condor_queue]"
    exit(-1)

input=sys.argv[1]
output=sys.argv[2]
queue=sys.argv[2] if len(sys.argv)>3 else 'local'
queue='tomorrow'

#check if library is there
if not os.path.isfile('runTTto2Lselection_C.so'):
    print 'Compiling analysis as shared library'
    ROOT.gSystem.Load("/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/fastjet/3.3.0-omkpbe/lib/libfastjet.so")
    ROOT.gSystem.CompileMacro('test/runTTto2Lselection.C','fk')


#create condor jobs
if queue != 'local':
    print 'Condor submission script is condor.sub'
    with open('condor.sub','w') as c:
        c.write('executable  = %s/HeavyIonsAnalysis/topskim/scripts/wrapAnalysis.sh\n'%os.environ['CMSSW_BASE'])
        c.write('output      = condor.out\n')
        c.write('error       = condor.err\n')
        c.write('log         = condor.log\n')
        c.write('arguments   = $(chunk) %s\n'%output)
        c.write('queue chunk matching (%s/*.root)\n'%input)



#nJobs=nEvts/evtsPerJob
#if nJobs*nEvts<evtsPerJob : nJobs=nJobs+1
#print '#events=',nEvts,'will be analysed in',nEvts/evtsPerJob,'jobs'
#
#os.system('mkdir -p '+outDir)
#for i in xrange(0,nJobs):
#    startEvt=i*evtsPerJob
#    outFile='%s/%s/%s_%d.root'%(cwd,outDir,tag,startEvt)    
#    cmd='bsub -q %s %s/scripts/wrapLocalAnalysisRun.sh %f %d %d %d %s %s' % (queue,cwd,cone,isPbPb,evtsPerJob,startEvt,outFile,vecList)
#    #cmd='%s/scripts/wrapLocalAnalysisRun.sh %f %d %d %d %s %s' % (cwd,cone,isPbPb,evtsPerJob,startEvt,outFile,vecList)
#    os.system(cmd)
#    
#print 'The output of the jobs will be named as %s_i.root and can be found in %s'%(tag,outDir)

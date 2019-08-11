import os
import sys
import subprocess
import optparse

def getDatasetComponents(opt):

    """ query components of the dataset """

    fList=[]

    #get files in the dataset
    print 'Querying',opt.dataset
    p = subprocess.Popen(['dasgoclient --query=\"file dataset=%s\"'%opt.dataset], 
                         stdout=subprocess.PIPE, 
                         stderr=subprocess.PIPE,
                         shell=True)
    out, err = p.communicate()

    dataset_files=out.split()
    nfiles=len(dataset_files)
    print 'I have %d jobs to submit...adding extra info'%nfiles

    #get parents, if required
    for i in range(nfiles):
        x=dataset_files[i]

        sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(nfiles))))
        sys.stdout.flush()

        if opt.addParent:
            p = subprocess.Popen(['dasgoclient --query=\"parent file=%s\"'%x],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True)
            out, err = p.communicate()
            fList.append( (x,out.split()) )
        else:
            fList.append( (x,[]) )
    
    return fList


def buildCondorFile(opt,fList,FarmDirectory):

    """ builds the condor file to submit the ntuplizer """

    jobTag=opt.jobTag

    cmssw=os.environ['CMSSW_BASE']

    secFileType=None

    #condor submission file
    condorFile='%s/condor_%s.sub'%(FarmDirectory,opt.jobTag)
    with open (condorFile,'w') as condor:
        condor.write('executable = {0}/worker_{1}.sh\n'.format(FarmDirectory,opt.jobTag))
        condor.write('output     = {0}/output_{1}.out\n'.format(FarmDirectory,opt.jobTag))
        condor.write('error      = {0}/output_{1}.err\n'.format(FarmDirectory,opt.jobTag))
        condor.write('log        = {0}/output_{1}.log\n'.format(FarmDirectory,opt.jobTag))
        condor.write('+JobFlavour = "nextweek"\n')
        OpSysAndVer = str(os.system('cat /etc/redhat-release'))
        if 'SLC' in OpSysAndVer:
            OpSysAndVer = "SLCern6"
        else:
            OpSysAndVer = "CentOS7"
        condor.write('requirements = (OpSysAndVer =?= "{0}")\n'.format(OpSysAndVer)) 

        for i in range(len(fList)):
            secFileList=','.join(fList[i][1])
            if '/RAW' in secFileList: secFileType='RAW'
            if '/AOD' in secFileList: secFileType='AOD'
            condor.write('arguments=%d %s %s\n'%(i,fList[i][0],secFileList))
            condor.write('queue 1\n')
                                
    #local worker script
    print 'Secondary file type',secFileType
    if secFileType == 'RAW':
        opt.extraOpts += ' redoProtonRecoFromRAW=True'
    if secFileType == 'AOD':
        opt.extraOpts += ' runWithAOD=True'
    print 'Additional extra opts will be set to',opt.extraOpts

    workerFile='%s/worker_%s.sh'%(FarmDirectory,opt.jobTag)
    with open(workerFile,'w') as worker:
        worker.write('#!/bin/bash\n')
        worker.write('WORKDIR=`pwd`\n')
        worker.write('echo "Working directory is ${WORKDIR}"\n')
        worker.write('cd %s\n'%cmssw)
        worker.write('eval `scram r -sh`\n')
        worker.write('cd ${WORKDIR}\n')
        worker.write('opts="lumiJson=%s inputFile=${2}"\n'%opt.lumiMask)
        worker.write('if [ ! -z "${3}" ]; then\n')
        worker.write('  opts="${opts} secInputFile=${3}"\n')
        worker.write('fi\n')
        if opt.proxy:
            worker.write('export X509_USER_PROXY=%s/myproxy509\n'%FarmDirectory)
        else:
            worker.write('echo "No proxy has been configured"\n')
        worker.write('opts="${opts} %s"\n'%opt.extraOpts)
        worker.write('echo "Running cmsRun with the following options: ${opts}"\n')
        worker.write('cmsRun %s/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py ${opts}\n'%cmssw)
        worker.write('xrdcp --force MiniEvents.root root://eoscms//%s/%s/MiniEvents_${1}.root\n'%(opt.output,opt.jobTag))
        worker.write('rm MiniEvents.root\n')

    return condorFile


def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--jobTag',      
                      dest='jobTag',      
                      help='job tag [%default]',  
                      default='MiniEvents',   
                      type='string')
    parser.add_option('--proxy',
                      dest='proxy',      
                      help='start proxy [%default]',  
                      default=False,
                      action='store_true')
    parser.add_option('--dryRun',
                      dest='dryRun',      
                      help='dry run (do not submit jobs to condor) [%default]',  
                      default=False,
                      action='store_true')
    parser.add_option('--output',    
                      dest='output',   
                      help='output to store [%default]',      
                      default='/store/cmst3/group/top/RunIIReReco/',
                      type='string')
    parser.add_option('--dataset',  
                      dest='dataset', 
                      help='dataset to process [%default]',   
                      default=None,    
                      type='string')
    parser.add_option('--addParent',   
                      dest='addParent',  
                      help='add parent [%default]',   
                      default=False,   
                      action='store_true')
    parser.add_option('--extraOpts',    
                      dest='extraOpts',
                      help='extra options to use in the ntuplizer cfg [%default]',
                      default='runOnData=True,era=era2017,runWithAOD=True',
                      type='string')
    parser.add_option('--lumiMask',        
                      dest='lumiMask',    
                      help='json with list of good lumis', 
                      default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt')
    (opt, args) = parser.parse_args()

    opt.extraOpts=' '.join(opt.extraOpts.split(','))

    #prepare directory with scripts
    cmssw=os.environ['CMSSW_BASE']
    FarmDirectory='%s/FarmLocalNtuple'%cmssw
    os.system('mkdir -p '+FarmDirectory)

    #start proxy
    if opt.proxy:
        os.system('voms-proxy-init --voms cms --out %s/myproxy509'%FarmDirectory)
        os.environ['X509_USER_PROXY']='%s/myproxy509'%FarmDirectory

    #get file list
    fList=getDatasetComponents(opt)

    #build condor submission script and launch jobs
    condor=buildCondorFile(opt,fList,FarmDirectory)
    if not opt.dryRun:
        os.system('condor_submit %s'%condor)

if __name__ == "__main__":
    sys.exit(main())

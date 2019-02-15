#!/usr/bin/env python

import os

def createJob(request,pset,dataset,lumiMask,lfnDirBase):
    
    """outputs a crag cfg file"""

    cfg_file='grid/%s_cfg.py'%request
    cfg=open(cfg_file,'w')
    cfg.write('from WMCore.Configuration import Configuration\n')
    cfg.write('import os\n')
    cfg.write('config = Configuration()\n')
    cfg.write('\n')
    cfg.write('config.section_("General")\n')
    cfg.write('config.General.requestName = "%s"\n' % request)
    cfg.write('config.General.workArea = "grid"\n')
    cfg.write('config.General.transferOutputs=True\n')
    cfg.write('\n')
    cfg.write('config.section_("JobType")\n')
    cfg.write('config.JobType.pluginName = "Analysis"\n')
    cfg.write('config.JobType.psetName = "%s"\n'%pset)
    cfg.write('config.JobType.disableAutomaticOutputCollection = False\n')
    cfg.write('config.JobType.pyCfgParams = [\'isPP=False\',\'maxEvents=-1\',\'outputFile=HiForest.root\']\n')
    cfg.write('config.JobType.outputFiles = [\'HiForest.root\']\n')
    cfg.write('\n')
    cfg.write('config.section_("Data")\n')
    cfg.write('config.Data.inputDataset = "%s"\n' % dataset)
    cfg.write('config.Data.inputDBS = "global"\n')
    cfg.write('config.Data.splitting = "Automatic"\n')
    cfg.write('config.Data.lumiMask = \'%s\'\n' %lumiMask)
    cfg.write('config.Data.ignoreLocality = False\n')    
    cfg.write('config.Data.publication = False\n')
    cfg.write('config.Data.outLFNDirBase = \"%s\"\n' % lfnDirBase)
    cfg.write('\n')
    cfg.write('config.section_("Site")\n')
    cfg.write('config.Site.storageSite = "T2_CH_CERN"\n')
    cfg.close()
    
    return cfg_file


#prepare output
os.system('mkdir -p grid')
pset="/afs/cern.ch/user/p/psilva/work/PbPb/CMSSW_10_3_1/src/HeavyIonsAnalysis/topskim/test/runForestAOD_pponAA_DATA_103X_PR.py"
lumiMask="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/DCSOnly/json_DCSONLY_HI.txt"
lfnDirBase="/store/group/cmst3/group/top/psilva/2018PbPb/"
submit=True
if submit:
    print 'Will submit crab jobs. Make sure you have sourced crab.sh'
    print 'source /cvmfs/cms.cern.ch/crab3/crab.sh'


for r,d in [
    ('DataPbPb_2018A_SingleMu_ZMM',   '/HISingleMuon/HIRun2018A-PbPbZMM-PromptReco-v2/RAW-RECO'),
    ('DataPbPb_2018A_DoubleMuon_ZMM', '/HIDoubleMuon/HIRun2018A-PbPbZMM-PromptReco-v2/RAW-RECO'),
    ('DataPbPb_2018A_HardProbes_ZEE', '/HIHardProbes/HIRun2018A-PbPbZEE-PromptReco-v2/RAW-RECO'),
    ('DataPbPb_2018A_HardProbes_EMu', '/HIHardProbes/HIRun2018A-PbPbEMu-PromptReco-v2/RAW-RECO'),
    ]:

    cfg=createJob(request=r,pset=pset,dataset=d,lumiMask=lumiMask,lfnDirBase=lfnDirBase)
    if submit : os.system('alias crab=\'/cvmfs/cms.cern.ch/crab3/crab-env-bootstrap.sh\' && crab submit -c %s' % cfg)

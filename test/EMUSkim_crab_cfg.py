from WMCore.Configuration import Configuration
import os
config = Configuration()

config.section_("General")
config.General.requestName = "Data2018_PbPbEmu"
config.General.workArea = "grid"
config.General.transferOutputs=True

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "/afs/cern.ch/user/p/psilva/work/PbPb/CMSSW_10_3_1/src/HeavyIonsAnalysis/topskim/test/runForestAOD_pponAA_DATA_103X_PR.py"
config.JobType.disableAutomaticOutputCollection = False
config.JobType.pyCfgParams = ['isPP=False','maxEvents=-1']

config.section_("Data")
config.Data.inputDataset = "/HIHardProbes/HIRun2018A-PbPbEMu-PromptReco-v1/RAW-RECO"
config.Data.inputDBS = "global"
config.Data.useParent = False
config.Data.splitting = "Automatic"
config.Data.ignoreLocality = False
config.Data.publication = False
config.Data.outLFNDirBase = "/store/group/cmst3/group/top/psilva/2018PbPb/"
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/DCSOnly/json_DCSONLY_HI.txt'



config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"

import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnData', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run this on real data"
                 )
options.register('runProtonFastSim', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Run proton fastsim for this angle"
                 )
options.register('runWithAOD', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "run with AOD"
                 )
options.register('noParticleLevel', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Do not run the particleLevel sequence"
                 )
options.register('era', 'era2017',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "era to use (configurations defined in python/EraConfig.py)"
                 )
options.register('outFilename', 'MiniEvents.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name"
                 )
options.register('baseJetCollection','slimmedJets',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Base jet collection"
                 )
options.register('inputFile', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "input file to process"
                 )
options.register('secInputFile', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "secondary input file to process"
                 )
options.register('lumiJson', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 'apply this lumi json file'
                 )
options.register('saveTree', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "save summary tree"
                 )
options.register('savePF', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'save PF candidates'
                 )
options.register('applyFilt', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'save PF candidates'
                 )
options.parseArguments()

#start process
from Configuration.StandardSequences.Eras import eras
process = cms.Process("MiniAnalysis", eras.ctpps_2016)      

#get the configuration to apply
from TopLJets2015.TopAnalysis.EraConfig import getEraConfiguration
globalTag, jecTag, jecDB, jerTag, jerDB = getEraConfiguration(era=options.era,isData=options.runOnData)

# Load the standard set of configuration modules
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#EGM customization
from TopLJets2015.TopAnalysis.customizeEGM_cff import customizeEGM
customizeEGM(process=process,era=options.era)

# global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag)
print 'Global tag is',globalTag

process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")
process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")
process.load("GeneratorInterface.RivetInterface.particleLevel_cfi") 

# particle level definitions
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
                                            inputPruned = cms.InputTag("prunedGenParticles"),
                                            inputPacked = cms.InputTag("packedGenParticles"),
                                            )
process.load('GeneratorInterface.RivetInterface.genParticles2HepMC_cfi')
process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")

#jet energy corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from TopLJets2015.TopAnalysis.customizeJetTools_cff import *
customizeJetTools(process=process,
                  jecDB=jecDB,
                  jecTag=jecTag,
                  jerDB=jerDB,
                  jerTag=jerTag,
                  baseJetCollection=options.baseJetCollection,
                  runOnData=options.runOnData)

#message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# set input to process
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/FE7ABAEB-4A42-E811-87A3-0CC47AD98D26.root'),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            )
if '2016' in options.era:
      process.source.fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv3/ST_t-channel_antitop_4f_mtop1715_inclusiveDecays_13TeV-powhegV2-madspin-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/120000/16CEB785-3FE6-E811-AAE8-FA163E9D74F8.root')

if options.runOnData:
      if '2017' in options.era:
            #process.source.fileNames = cms.untracked.vstring('/store/data/Run2017F/ZeroBias/MINIAOD/31Mar2018-v1/30000/8604984F-DF37-E811-ACFE-008CFA197B74.root')
            process.source.fileNames = cms.untracked.vstring('/store/data/Run2017B/SingleMuon/MINIAOD/31Mar2018-v1/90000/46CE7E24-EA37-E811-95CB-0025905A6132.root')
            if options.runWithAOD:
                  print 'Adding secondary filenames'
                  process.source.secondaryFileNames = cms.untracked.vstring(['/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/70003/C4FDB602-77D8-E711-B975-02163E011E63.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/70003/089281B2-74D8-E711-96D9-02163E019CF1.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40001/74593CEA-11D9-E711-A49E-FA163E62ECD8.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40001/4CC8D127-EED8-E711-91CC-02163E016CCF.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40001/3AFF8583-F0D8-E711-ADB8-FA163EA84AFA.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40001/24A6AA06-21D9-E711-B5AD-FA163E3BC51E.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40001/1A02BCE9-FDD8-E711-BB87-FA163E408AED.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40001/0A90781C-EBD8-E711-9198-FA163E9875C8.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40001/0481E927-F2D8-E711-A11A-FA163EFB3309.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40000/FC54E9BB-D8D8-E711-8076-FA163EC00925.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40000/DA1366B2-D2D8-E711-ACEE-FA163E24EC20.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40000/8CCE5537-D3D8-E711-B316-FA163E653C95.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40000/8492E4C5-C2D8-E711-9631-FA163EA53599.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40000/82906F7E-C4D8-E711-8A90-FA163EDC0DF0.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40000/38AFFEDF-C7D8-E711-85A5-FA163EE5AC62.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40000/16C2CF2F-CAD8-E711-8630-FA163E9C5649.root'
                                                                             ])
      elif '2018' in options.era:
            process.source.fileNames = cms.untracked.vstring('/store/data/Run2018D/SingleMuon/MINIAOD/PromptReco-v2/000/325/001/00000/E4A5114B-87A1-5244-A82F-C810E06A6EFA.root')
            if options.runWithAOD:
                  print 'Adding secondary filenames'
                  process.source.secondaryFileNames = cms.untracked.vstring([
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/E4C865F4-CAF0-A047-A4C9-055A8731DBD5.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/00E04CC2-4EFE-814D-AA1B-A260246B2B61.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/DE0267EE-4B93-A445-B368-741AFFBC3790.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/E0F4348B-08AB-274C-8A5A-1D339A3C7C65.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/E8370901-46DE-564F-BB1D-7FBE1C48FBFD.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/F1F34366-D7F5-974E-9B23-17D8B8B81B22.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/DD605145-693E-1D41-9B8F-732AF82190A4.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/CF5AF7A9-C27E-5D42-AE56-D430F8CFDEA3.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/CC6892AB-7588-D64D-976F-0F113395F712.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/B17BA02E-1988-BE45-9CB2-2B550EEC5CA0.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/A29DA982-5015-B943-BE0E-54D97CB85262.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/9F010069-21D9-CE4F-A25E-3CC724D25537.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/819397B5-9EBE-7240-9879-2E16CE5E5DC4.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/77BF9FFC-4A63-804E-B7A5-82D1D62B0465.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/5CDFB33C-26D6-9849-8C3E-8A3ADCBE9DCA.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/4EFF4077-DEE4-3242-AB47-47863E82E586.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/4CECB062-26A4-A041-86B0-0E10EFECFC52.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/1BD78D97-19BD-9346-AE92-6CC2C5447290.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/508C6852-B14B-C146-AC2B-CBAA0572262E.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/22E112CB-975A-CA4B-A3AB-497058F2CCB3.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/1BCF511B-C81F-A04F-BD01-7A40E1EC5F60.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/129C4BEF-C3A6-BA45-AF85-9352A7EF5D43.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/0FAB9271-24D9-7041-93AC-0DAE36520AEF.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/0B0BB432-3587-274F-A5B2-9E3D501D6BA3.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/BAD48B74-91BA-EC45-88DF-D12E12079F19.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/BE535651-B675-E54E-A0F8-BD5FD6EE030F.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/08D02A96-3056-374C-BA74-87A0F3E77BB6.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/041D9B36-E41A-DD4A-995A-4642D1CCD8F7.root',
                        '/store/data/Run2018D/SingleMuon/RAW/v1/000/325/001/00000/E4B7AF86-6865-BF46-A9FA-F10D6618F992.root'
                  ])
      else:
            process.source.fileNames = cms.untracked.vstring('/store/data/Run2016B/SingleMuon/MINIAOD/17Jul2018_ver2-v1/00000/0219DD4A-7D8C-E811-B7EB-001E67248A43.root')
            if options.runWithAOD:
                  print 'Adding secondary filenames'
                  process.source.secondaryFileNames = cms.untracked.vstring(['/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/110000/B6479036-9486-E711-A6CF-003048FFD734.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/110001/0A8A7032-9086-E711-BD1B-0025905A60C6.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/110001/22FFA42A-9086-E711-9B49-0CC47A7C3472.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70000/767AE9A5-BA81-E711-9A8E-0CC47A745282.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70000/D67D18A5-BA81-E711-B6C2-0CC47A7C361E.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70000/E06164A7-BA81-E711-B612-0025905B858A.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70001/302E2840-C481-E711-8913-002618FDA259.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70001/687FA66D-BC81-E711-B41D-0CC47A745282.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70001/6890466C-BC81-E711-9C62-0CC47A4C8E1C.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70001/E64E346F-BC81-E711-93BC-0CC47A4C8F18.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70001/EC3632F7-C081-E711-B219-0CC47A7452DA.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70001/F6499DE0-C281-E711-9835-0CC47A78A33E.root'])

if options.inputFile:
      fileList=[]
      secFileList=[]
      if '.root' in options.inputFile :
            fileList=[options.inputFile]
            if options.secInputFile and '.root' in options.secInputFile:
                  secFileList=[options.secInputFile]
                  print 'Will run also on',secFileList
      else:
            import os
            inDir=options.inputFile
            fileList = ['file:'+os.path.join(inDir,f) for f in os.listdir(inDir)]
      process.source.fileNames = cms.untracked.vstring(fileList)     
      process.source.secondaryFileNames = cms.untracked.vstring(secFileList)
      
print  "Processing",process.source.fileNames

#apply lumi json, if passed in command line
if options.lumiJson:
    print 'Lumi sections will be selected with',options.lumiJson
    from FWCore.PythonUtilities.LumiList import LumiList
    myLumis = LumiList(filename = options.lumiJson).getCMSSWString().split(',')
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.source.lumisToProcess.extend(myLumis)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outFilename)
                                  )

#analysis
from TopLJets2015.TopAnalysis.miniAnalyzer_cfi import  ANALYSISJETIDS,ANALYSISTRIGGERLISTS
process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')
print 'MiniAnalyzer configuration is as follows:'
process.analysis.saveTree  = cms.bool(options.saveTree)
process.analysis.savePF    = cms.bool(options.savePF)
process.analysis.applyFilt = cms.bool(options.applyFilt)
print '\t save tree=',options.saveTree,' save PF=',options.savePF
if 'era2017' in options.era:
      process.analysis.jetIdToUse=ANALYSISJETIDS[2017]
      process.analysis.triggersToUse=ANALYSISTRIGGERLISTS[2017]
      print '\t Using 2017 triggers/jet ids'
else:
      process.analysis.jetIdToUse=ANALYSISJETIDS[2016]
      process.analysis.triggersToUse=ANALYSISTRIGGERLISTS[2016]
      print '\t Using 2016 triggers/jet ids'
if options.runOnData:
      process.analysis.metFilterBits = cms.InputTag("TriggerResults","","RECO")
      print '\t will save met filter bits'

#schedule execution
toSchedule=[]
if process.egammaPostReco:
      toSchedule.append( process.egammaPostReco )
if process.updatedPatJetsUpdatedJECBTag:
      process.custom_jec_seq=cms.Sequence(process.QGTagger * process.patJetCorrFactorsUpdatedJECBTag * process.updatedPatJetsUpdatedJECBTag)
      process.custom_jec=cms.Path(process.custom_jec_seq)
      toSchedule.append( process.custom_jec)
if process.fullPatMetSequenceModifiedMET:
      process.custom_met=cms.Path(process.fullPatMetSequenceModifiedMET)
      toSchedule.append(process.custom_met)
if not (options.runOnData or options.noParticleLevel):
      process.mctruth=cms.Path(process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel)
      toSchedule.append( process.mctruth )

if options.runOnData:
      from TopLJets2015.TopAnalysis.protonReco_cfg import ctppsCustom
      ctppsCustom(process)
      process.ppsReco=cms.Path(process.ppsSeq)
      toSchedule.append(process.ppsReco)

if options.runProtonFastSim:
      from TopLJets2015.TopAnalysis.protonReco_cfg import setupProtonReco
      setupProtonReco(process,options.runProtonFastSim)
      toSchedule.append(process.pps_fastsim)
      toSchedule.append(process.pps_simulation_step)
      toSchedule.append(process.pps_reco_step)

      process.analysis.tagRecoProtons = cms.InputTag('ctppsProtonReconstructionOFDB')


process.ana=cms.Path(process.analysis)
toSchedule.append( process.ana )
if options.runOnData:
      process.analysis.tagRecoProtons = cms.InputTag('ctppsProtons')      
                           
process.schedule=cms.Schedule( (p for p in toSchedule) )
print process.schedule

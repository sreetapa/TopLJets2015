import FWCore.ParameterSet.Config as cms
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
# EGM corrections :  https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2

def customizeEGM(process,era):

    egmEra='2017-Nov17ReReco'
    if '2016' in era: egmEra='2016-Legacy'
    if '2018' in era: egmEra="2018-Prompt"
    setupEgammaPostRecoSeq(process,isMiniAOD=True,runEnergyCorrections=False,era=egmEra)
    #setupEgammaPostRecoSeq(process,isMiniAOD=True,applyEnergyCorrections=True,applyVIDOnCorrectedEgamma=True,era=egmEra)
    #process.electronMVAValueMapProducer.src = cms.InputTag("gedGsfElectrons","","MINIAOD")
    #process.photonMVAValueMapProducer.src = cms.InputTag("gedPhotons","","MINIAOD")
    
    process.egammaPostReco=cms.Path(process.egammaPostRecoSeq)

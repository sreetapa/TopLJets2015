import ROOT
import os
from subprocess import Popen, PIPE

qcdSamples=['/store/cmst3/group/hintt/psilva/PP5TeV/MC/Pythia8_bJet30_pp502_TuneCUETP8M1/crab_MCpp502_bJet_pthat30_pythia8/160824_090828/0000',
            '/store/cmst3/group/hintt/psilva/PP5TeV/MC/Pythia8_bJet50_pp502_TuneCUETP8M1/crab_MCpp502_bJet_pthat50_pythia8/160824_142558/0000',
            '/store/cmst3/group/hintt/psilva/PP5TeV/MC/Pythia8_bJet80_pp502_TuneCUETP8M1/crab_MCpp502_bJet_pthat80_pythia8/160824_142546/0000',
            '/store/cmst3/group/hintt/psilva/PP5TeV/MC/Pythia8_bJet100_pp502_TuneCUETP8M1/crab_MCpp502_bJet_pthat100_pythia8/160824_142629/0000',
            '/store/cmst3/group/hintt/psilva/PP5TeV/MC/Pythia8_bJet120_pp502_TuneCUETP8M1/crab_MCpp502_bJet_pthat120_pythia8/160823_153738/0000'            
            ]

ptHat=[30,50,80,100,120]
ptHatEventCounters=[0.]*5
ptHatTheoryXsec=[2829690.,271417.,4785.,1762.,746]


#mount locally EOS                                                                                                                                                                                             
eos_cmd = '/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'
Popen([eos_cmd, ' -b fuse mount', 'eos'],stdout=PIPE).communicate()

#loop over samples and fill chain of all QCD trees
chain=ROOT.TChain('hiEvtAnalyzer/HiTree')
for directory in qcdSamples:
    baseDir='eos/cms/%s'%directory
    for url in os.listdir(baseDir):
        chain.AddFile(baseDir+'/'+url)

#loop over all events to determine total number of events in the different pthat bins
for i in xrange(0,chain.GetEntries()):
    chain.GetEntry(i)
    pthat=chain.pthat
    if pthat<50    : ptHatEventCounters[0]=ptHatEventCounters[0]+1.0
    elif pthat<80  : ptHatEventCounters[1]=ptHatEventCounters[1]+1.0
    elif pthat<100 : ptHatEventCounters[2]=ptHatEventCounters[2]+1.0
    elif pthat<120 : ptHatEventCounters[3]=ptHatEventCounters[3]+1.0
    else           : ptHatEventCounters[4]=ptHatEventCounters[4]+1.0

#print results for the weights
for j in xrange(0,len(ptHat)):

    ptmin=ptHat[j]
    ptmax=float('inf')
    if j<len(ptHat)-1 : ptmax=ptHat[j+1]

    thXsec=ptHatTheoryXsec[j]
    if j<len(ptHat)-1 : thXsec=thXsec-ptHatTheoryXsec[j+1]

    wgt=thXsec/ptHatEventCounters[j]

    print ptmin,'<pthat<',ptmax,'xsec=',thXsec,'pb','#events=',ptHatEventCounters[j],'weight=',wgt

#unmount locally EOS                                                                                                                                                                                           
Popen([eos_cmd, ' -b fuse umount', 'eos'],stdout=PIPE).communicate()

#!/usr/bin/env python

import ROOT
import os,sys
import random

OUTSTORE='/store/cmst3/user/psilva/Mix_Pythia6_bJet_pp502_Hydjet_MB'

"""
"""
def mixJets(jets,nevts,baseOutput):

    evList=range(0,jets.GetEntries())
    random.shuffle(evList)

    subsetCtr=0
    while len(evList)>nevts:

        mixedJets=jets.CloneTree(0)

        for i in xrange(0,nevts):
            jets.GetEntry( evList.pop() )
            mixedJets.Fill()
        subsetCtr+=1

        outName='%s_%d.root'%(baseOutput,subsetCtr)
        fOut=ROOT.TFile(outName,'RECREATE')
        mixedJets.Write()
        fOut.Close()

        baseName=os.path.basename(outName)
        os.system('xrdcp -f %s root://eoscms//eos/cms/%s/%s'%(outName,OUTSTORE,baseName))
        os.system('rm -v %s'%outName)
    
    print '[mixJets] produced %d output files'%subsetCtr


"""
"""
def main():

    if len(sys.argv)<2:
        print 'mixTreesForTraining.py eventsPerFile dir1 dir2 ...'
        print '(check hardcoded OUTSTORE in the file)'
        sys.exit(0)

    nevts    = int(sys.argv[1])
    fileList = [os.path.join(baseDir,f) for baseDir in sys.argv[2:] for f in os.listdir(baseDir)]
    random.shuffle(fileList)

    print 'Found %d files to mix'%len(fileList)

    mixCounter,fCounter=0,0
    jets=ROOT.TChain('jets')
    while len(fileList)>0:
        jets.AddFile( fileList.pop() )
        fCounter+=1
        if jets.GetEntries()>nevts: 
            print 'Calling mixJets with %d jets to mix in chunks of %d, from %d files'%(jets.GetEntries(),nevts,fCounter)
            mixJets(jets=jets,
                    nevts=nevts,
                    baseOutput='/tmp/psilva/MixedJets_%d'%mixCounter)
            mixCounter+=1
            fCounter=0
            jets.Reset()
            print 'Done, jet chain has now %d files'%jets.GetListOfFiles().GetEntriesFast()
        



if __name__ == "__main__":
    main()

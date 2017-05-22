#!/usr/bin/env python

import ROOT
import os,sys
import random
import optparse

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

    usage = 'usage: %prog [options] {dir,file}1 {dir,file}2'
    parser = optparse.OptionParser(usage)
    parser.add_option('-n', '--nevts',   dest='nevts',    help='number of events per file [%default]',  default=200000, type='int')
    parser.add_option('-c', '--counter', dest='counter',  help='counter tag [%default]',                default=0,      type='int')
    parser.add_option(      '--sub',     dest='sub',      help='submit [%default]',                     default=False,  action='store_true')
    (opt, args) = parser.parse_args()
    

    #if arguments contain already ROOT files, just mix them
    if '.root' in args[0]:
        jets=ROOT.TChain('jets')
        for f in args: jets.AddFile( f )
        print 'Will mix jets from %d files and %d events'%(len(args),jets.GetEntries())
        mixJets(jets=jets, nevts=opt.nevts,baseOutput='/tmp/MixedJets_%d'%opt.counter)

    #otherwise create the files list and submit to the batch
    else:
        fileList = [os.path.join(x,f) for x in args for f in os.listdir(x)]
        random.shuffle(fileList)
        print 'Found %d files to mix'%len(fileList)

        mixCounter=0
        jets=ROOT.TChain('jets')
        mixFileList=[]
        while len(fileList)>0:
            url=fileList.pop()
            mixFileList.append( url )
            jets.AddFile( url )

            if jets.GetEntries()>opt.nevts: 

                scriptName='mixscript_%d.sh'%mixCounter
                with open(scriptName,'w') as fout:
                    fout.write('#!/bin/bash\n')
                    fout.write('cd %s/src/topskim/test/\n'%os.environ['CMSSW_BASE'])
                    fout.write('eval `scram r -sh`\n')
                    fout.write('python ${CMSSW_BASE}/src/topskim/test/mixTreesForTraining.py --nevts %d --counter %d %s\n'%(opt.nevts,mixCounter,' '.join(mixFileList)))

                if opt.sub:
                    os.system('chmod u+x %s'%scriptName)
                    os.system('bsub -q 2nd %s/src/topskim/test/%s'%(os.environ['CMSSW_BASE'],scriptName))                    

                mixCounter+=1
                mixFileList=[]
                jets.Reset()

        



if __name__ == "__main__":
    main()

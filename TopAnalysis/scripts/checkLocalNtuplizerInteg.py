import ROOT
import sys
import os
from subprocess import check_output
from checkLocalAnalysisInteg import checkIntegrity

#parse condor file
condor=sys.argv[1]
f=open(condor,'r')
lines = [line.rstrip('\n') for line in f]
nfiles = len([ x for x in lines if 'arguments' in x ])
f.close()

#parse worker file
worker=condor.replace('condor','worker')
worker=worker.replace('sub','sh')
f=open(worker,'r')
for line in f:
    if not 'xrdcp' in line : continue
    outfileName=line.split()[3]
    outfileName=outfileName.replace('${1}','{0}')
f.close()

#check if file is present
reSub=[]
for i in range(nfiles):
    res=checkIntegrity(outfileName.format(i),'analysis/data')
    if res[0] is None : continue
    reSub.append(i+1)

print len(reSub),'/',nfiles,'need to be resubmitted'

newCondor=condor.replace('.sub','_resub.sub')
with open(newCondor,'w') as f:
    iarg=0

    skipNext=False
    for il in range(len(lines)):
        l=lines[il]

        l=l.replace('longlunch','tomorrow')

        if skipNext:
            skipNext=False
            continue

        writeLine=False
        if 'arguments=' in l:
            iarg+=1
            if iarg in reSub:
                writeLine=True
            else:
                skipNext=True
        else:
            writeLine=True

        if not writeLine: continue
        f.write(l+'\n')
print 'New condor submission file @',newCondor


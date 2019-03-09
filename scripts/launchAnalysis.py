import ROOT
import os,sys
from subprocess import Popen, PIPE

if len(sys.argv)<4:
    print "python scripts/launchAnalysis.py input output tag"
    exit(-1)

input=sys.argv[1]
output=sys.argv[2]
tag=sys.argv[3]


#prepare output
os.system('mkdir -p %s'%output)

#create condor jobs
print 'Condor submission script is condor_%s.sub'%tag
with open('condor_%s.sub'%tag,'w') as c:
    c.write('executable  = %s/src/HeavyIonsAnalysis/topskim/scripts/wrapAnalysis.sh\n'%os.environ['CMSSW_BASE'])
    c.write('output      = condor_%s.out\n'%tag)
    c.write('error       = condor_%s.err\n'%tag)
    c.write('log         = condor_%s.log\n'%tag)
    
    chunks=os.listdir(input)
    for i in range(len(chunks)):
        c.write('arguments   = {0}/{1} {2}/{3}_{4}.root\n'.format(input,chunks[i],output,tag,i))
        c.write('queue 1\n')

#submit to condor
os.system('condor_submit condor_%s.sub'%tag)

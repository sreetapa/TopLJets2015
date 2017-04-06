#!/usr/bin/env python

import das_client
import sys
import os
import json
import optparse

def getListOfFiles(dataset):
      das_query_string = 'das_client.py --query="file dataset=%s | grep file.name" --limit=0' % (dataset)
      das_output = os.popen(das_query_string).read() #query DAS via bash
      das_output = das_output.replace("\"","")
      return das_output.split('\n')

def main():
      dataset=sys.argv[1]

      usage = 'usage: %prog [options]'
      parser = optparse.OptionParser(usage)
      parser.add_option('-q', '--queue'       ,dest='queue'       ,help='batch queue'    ,default='workday')
      parser.add_option('-p', '--startProxy'  ,dest='startProxy'  ,help='startProxy'     ,default=False,  action='store_true')
      (opt, args) = parser.parse_args()

      if opt.startProxy:
            os.system('voms-proxy-init --voms cms')
            os.system('cp /tmp/x509up_u100724 ~/private/cur_proxy')

      #prepare output
      tag=dataset.split('/')[1]
      outDir='/store/cmst3/group/hintt/psilva/%s/'%tag
      os.system('eos mkdir %s'%outDir)

      FARMDIR='FARM/%s'%tag
      os.system('mkdir -p %s'%FARMDIR)
      fileList=getListOfFiles(dataset)
      script='/afs/cern.ch/user/p/psilva/scratch0/CMSSW_7_5_8_patch3/src/topskim/test/runLocalSkim.sh'
      for i in xrange(0,len(fileList)):
            if len(fileList[i])==0: continue
            scriptArgs='-i %s -o JetTree_%d.root -s %s'%(fileList[i],i,outDir)            
            condorScript='%s/job_%d.sh'%(FARMDIR,i)
            with open(condorScript,'w') as f:
                  f.write('executable            = %s\n'%script)
                  f.write('arguments             = %s\n'%scriptArgs)
                  f.write('output                = %s/job_%d.$(ClusterId).$(ProcId).out\n'%(FARMDIR,i))
                  f.write('error                 = %s/job_%d.$(ClusterId).$(ProcId).err\n'%(FARMDIR,i))
                  f.write('log                   = %s/job_%d.$(ClusterId).$(ProcId).log\n'%(FARMDIR,i))
                  f.write('+JobFlavour = "%s"\n'%opt.queue)
                  f.write('queue\n')

            os.system('condor_submit %s'%condorScript)


if __name__ == "__main__":
    main()

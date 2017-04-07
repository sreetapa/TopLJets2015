import ROOT
import commands
import sys

def getPtEta(t,c):
    pteta={}
    for tag,cut in [ ('b','isB==1'),
                     ('unm','isUnmatched==1'),
                     ('light','isB==0 && isUnmatched==0')]:
        t.Draw('fabs(eta):min(pt,200) >> pteta(50,30,200,50,0,2.1)',cut,'colz')
        pteta[tag]=ROOT.TGraph2D(c.GetPrimitive('pteta'))
        pteta[tag].SetName(tag)
    return pteta

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',400,400)

inDir=sys.argv[1]
fileList=commands.getstatusoutput('eos ls %s'%inDir)[1].split('\n')
t=ROOT.TChain('jets')
for f in fileList: t.AddFile('root://eoscms//eos/cms/%s/%s'%(inDir,f))

#build pt/eta histos
pteta=getPtEta(t,c)

#start control histos
histos={}
colors={'b':2,'light':1,'unm':ROOT.kGray}
for flav in pteta:
    histos[flav+'_csvV2']=ROOT.TH1F(flav+'_csvV2',flav+';CSVv2;Jets',50,0,1)
    histos[flav+'_csvV2'].SetLineColor(colors[flav])
    histos[flav+'_csvV2'].Sumw2()
    histos[flav+'_csvV2'].SetDirectory(0)

#fill histos
for i in xrange(0,t.GetEntries()):
    t.GetEntry(i)
    flav='light'
    if t.isB==1 : flav='b'
    elif t.isUnmatched: flav='unm'
    
    wgt=pteta[flav].Interpolate(min(t.pt,200),abs(t.eta))
    if wgt==0 : continue

    histos[flav+'_csvV2'].Fill(t.csvV2,1./wgt) 

#show histos
c.Clear()
histos['unm_csvV2'].Draw('hist')
histos['light_csvV2'].Draw('histsame')
histos['b_csvV2'].Draw('histsame')
c.BuildLegend()
c.Modified()
c.Update()
raw_input()

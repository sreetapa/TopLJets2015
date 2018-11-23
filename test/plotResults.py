import ROOT
from Plot import *

lumi=(1e-6)*60 #*(1e3)
APb=208
k5p5to5p02=75./69.

samples=[('MCpp5p2TeV_DY.root',  'Z/#gamma^{*}',  2010*(APb**2)*lumi*k5p5to5p02,  852),         
         ('Data2018_EM_SS.root', 'Fakes (data)',  1,                              17),
         ('MCpp5p2TeV_TT.root',  't#bar{t}',      69*(APb**2)*lumi*k5p5to5p02,    633),
         ('Data2018_EM.root',    'Data',          None,                           1),         
         ('Data2018_ZEE.root',   'Data',          None,                           1),

         ]
dists=['mll','ptll','dphill','njets','nbjets','jpt','jeta','jcsv']
plots=[ ch+'_'+d for ch in ['zee','em'] for d in dists ]


def getPlot(pname,inF,scale=None):
    fIn=ROOT.TFile.Open(inF)
    h=fIn.Get(pname)
    try:
        if not 'jets' in pname:
            h.Rebin()
        h.SetName(pname+'_'+inF)
        h.SetDirectory(0)
        if scale: h.Scale(scale)
    except:
        h=None
    fIn.Close()
    return h

ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

for d in plots:
    p=Plot(d,com='5.02 TeV')
    histos=[]
    for f,proc,xsec,ci in samples:
        print xsec
        histos.append(getPlot(d,f,xsec))
        if histos[-1]:
            fixExtremities(histos[-1],False,False)
            p.add(h=histos[-1],
                  title=proc,
                  color=ci,
                  isData=False if xsec else True,
                  spImpose=False,
                  isSyst=False)
    p.show(outDir='./',lumi=lumi)
    p.reset()


          

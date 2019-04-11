import ROOT
import sys
import os
import numpy as np
from HeavyIonsAnalysis.topskim.LeptonObject import *

rhoTitle={'lep_chrho':'#rho_{ch}(l)',
          'lep_phrho':'#rho_{#gamma}(l)',
          'lep_nhrho':'#rho_{nh}(l)',
          'chrho_restr':'#rho_{ch}(|#eta|)',
          'phorho_restr':'#rho_{#gamma}(|#eta|)',
          'nhrho_restr':'#rho_{nh}(|#eta|)',
          'chrho':'#rho_{ch}',
          'phorho':'#rho_{#gamma}',
          'nhrho':'#rho_{nh}',
          'rho':'#rho',
          None:'None'}
isoTitle={'iso':'Total isolation',
          'chiso':'Charged isolation',
          'nhiso':'Neutral hadron isolation',
          'phiso':'Photon isolation',
          'isonocomp'         : "I",
          'chisonocomp'       : "I_{ch}",
          'isoglobalcomp'     : "[I-UE(#rho)]", 
          'chisoglobalcomp'   : "[I_{ch}-UE(#rho_{ch})]", 
          'isoglobalpartcomp' : "#Sigma_{k}[I_{k}-UE(#rho_{k})]",
          'isoregpartcomp'    : "#Sigma_{k}[I_{k}-UE(#rho_{k}^{region})]",
          'isoleppartcomp'    : "#Sigma_{k}[I_{k}-UE(#rho_{k}^{lepton})]",
          }


def dileptonFinder(eleColl):

    """takes the first lepton pair with a mass close to the Z"""

    dil=[]
    isZ=None
    isSS=None
    for i in range(len(eleColl)):
        for j in range(i+1,len(eleColl)):
            ll=eleColl[i].p4+eleColl[j].p4
            isSS=True if eleColl[i].charge*eleColl[j].charge>0 else False
            isZ=True if abs(ll.M()-91)<15 and not isSS else False
            if not isZ and not isSS: continue            
            dil=[eleColl[i],eleColl[j]]
    return dil,isZ,isSS

def getROC(sig,bkg,name):
    gr=ROOT.TGraph()
    gr.SetName(name+'_roc')
    tot_sig=sig.Integral()
    tot_bkg=bkg.Integral()
    nbins=sig.GetNbinsX()
    best_effsig=0
    best_effbkg=0
    best_xbin=nbins+1
    for xbin in range(nbins+1):
        effsig=sig.Integral(0,xbin+1)/tot_sig
        effbkg=bkg.Integral(0,xbin+1)/tot_bkg
        gr.SetPoint(xbin,effsig,effbkg)
        if abs(effsig-0.89)>abs(best_effsig-0.89) : continue
        best_effbkg=effbkg
        best_effsig=effsig
        best_xbin=xbin
    print name,best_effsig,best_effbkg,'cut<',best_xbin,sig.GetXaxis().GetBinCenter(best_xbin+1)
    return gr

def header(extraTxt=[]):
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.04)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.12,0.97,'#bf{CMS} #it{preliminary}')
    txt.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
    txt.DrawLatex(0.95,0.97,'#scale[0.8]{1618 #mub^{-1} (#sqrt{s_{NN}}=5.02 TeV)}')

    txt.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
    txt.SetTextSize(0.035)
    for i in range(len(extraTxt)):
        txt.DrawLatex(0.15,0.9-i*0.04,extraTxt[i])


#open all electron files
url=sys.argv[1]
t=ROOT.TChain('tree')
for f in os.listdir(url):
    if not 'SkimElectrons' in f or not '.root' in f: continue
    t.Add(os.path.join(url,f))

#loop over events
leptons=[]
for iev in range(t.GetEntries()):
    t.GetEntry(iev)

    #get dilepton
    eleColl=getElectrons(t)
    dil,isZ,isSS=dileptonFinder(eleColl)
    if not isZ and not isSS: continue

    #centrality
    cenbin=t.cenbin

    leptons.append( [dil[0],isZ,isSS,cenbin] )
    leptons.append( [dil[1],isZ,isSS,cenbin] )

print 'Found',len(leptons),'electrons for optimization'

#centrality bins
zCentralities=[]
for _,isZ,_,cenbin in leptons:
    if not isZ: continue
    zCentralities.append(cenbin)
zcenq=np.percentile(zCentralities,[0,20,40,60,80,100],axis=0)


ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kLightTemperature)
c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.05)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.12)
c.SetBottomMargin(0.1)
from array import array

rhoProf=ROOT.TF1('rhoprof','[0]/([1]+x)',0,100)
isoProf=ROOT.TF1('isoprof','[0]*pow(x+[1],2)+[2]*(x+[1])',0,1000)
#cenProf=ROOT.TF1('cenprof','[0]/(1+exp([1]*x))',0,100)
cenProf=ROOT.TF1('cenprof','exp([0]*x)*([1]*x*x+[2]*x+[3])',0,100)
finalIsols={
    #'isonocomp'        :[('iso',None)],
    #'chisonocomp'      :[('chiso',None)],
    'isoglobalcomp'    :[('iso','rho')],    
    #'chisoglobalcomp'  :[('chiso','chrho')],
    #'isoglobalpartcomp':[('chiso','chrho'),('nhiso','nhrho'),('phiso','phorho')],
    #'isoregpartcomp'   :[('chiso','chrho_restr'),('nhiso','nhrho_restr'),('phiso','phorho_restr')],
    #'isoleppartcomp'   :[('chiso','lep_chrho'),('nhiso','lep_nhrho'),('phiso','lep_phrho')]
    }
rocs=[]
rocshcen=[]
for isoKey in finalIsols:
    isoComponents=finalIsols[isoKey]
    isoVals=[0]*len(leptons)
    for iso,rho in isoComponents:

        isoCompVals=[]
        for l,isZ,isSS,cenbin in leptons:
            if not isZ: continue
            iiso=getattr(l,iso)
            irho=l.rho[rho] if rho else 0.
            isoCompVals.append([iiso,irho,cenbin])

        q=np.percentile(isoCompVals,range(0,110,10),axis=0)

        #parameterize rho vs centrality
        h=ROOT.TH2F('%svscen'%rho,
                    ';Centrality;%s'%rhoTitle[rho],
                    5,array('d',zcenq),
                    10,array('d',q[:,1]))
        for _,irho,cenbin in isoCompVals:
            h.Fill(cenbin,irho)
        c.SetRightMargin(0.12)
        h.Draw('colz')
        h.SetTitle('scatter')
        px=h.ProfileX()
        px.SetMarkerStyle(20)
        px.SetFillStyle(0)
        px.Draw('e1same')
        px.SetTitle('profile')
        px.Fit(rhoProf,'MRQ+')
        extraTxt=['#rho=%3.3f/(%3.3f+c)'%(rhoProf.GetParameter(0),rhoProf.GetParameter(1))]
        header(extraTxt)
        c.Modified()
        c.Update()
        c.SaveAs(h.GetName()+'.png')

        #parameterize average profile of component versus rho
        q=np.percentile(isoCompVals,range(0,110,10),axis=0)
        h=ROOT.TH2F('%svs%s'%(rho,iso),
                    ';%s;%s'%(rhoTitle[rho],isoTitle[iso]),
                    10,array('d',q[:,1]),
                    10,array('d',q[:,0]))
        for iiso,irho,_ in isoCompVals:
            h.Fill(irho,iiso)
        h.Draw('colz')
        h.SetTitle('scatter')
        px=h.ProfileX()
        px.SetMarkerStyle(20)
        px.SetFillStyle(0)
        px.Draw('e1same')
        px.SetTitle('profile')
        isoProf.SetParLimits(0,0,1e3)
        px.Fit(isoProf,'MRQ+')
        extraTxt=['UE=%3.4f#rho^{2}+%3.4f#rho+%3.2f'%(isoProf.GetParameter(0),isoProf.GetParameter(1),isoProf.GetParameter(2))]
        header(extraTxt)
        c.Modified()
        c.Update()
        c.SaveAs(h.GetName()+'.png')
        for i in range(len(leptons)):
            iiso=getattr(leptons[i][0],iso)
            if rho:
                iiso-=isoProf.Eval(leptons[i][0].rho[rho])
            isoVals[i]+=iiso #/leptons[i][0].p4.Pt()


    q=np.percentile(isoVals,range(0,110,10),axis=0)
    #h=ROOT.TH2F('cenvs%s'%isoKey,       ';Centrality;%s'%isoTitle[isoKey], 5,array('d',zcenq),50,-3,3)
    #hbkg=ROOT.TH2F('bkgcenvs%s'%isoKey, ';Centrality;%s'%isoTitle[isoKey], 5,array('d',zcenq),50,-3,3)

    h=ROOT.TH2F('cenvs%s'%isoKey,       ';Centrality;%s [GeV]'%isoTitle[isoKey], 5,array('d',zcenq), 10,array('d',q))
    hbkg=ROOT.TH2F('bkgcenvs%s'%isoKey, ';Centrality;%s [GeV]'%isoTitle[isoKey], 5,array('d',zcenq), 10,array('d',q))
    for i in range(len(leptons)):
        l,isZ,isSS,cenbin=leptons[i]
        if isZ:
            h.Fill(cenbin,isoVals[i])
        if isSS:
            hbkg.Fill(cenbin,isoVals[i])

    px=h.ProfileX()
    px.Fit(cenProf)

    huncor=ROOT.TH2F('cenvsuncor%s'%isoKey,       ';Centrality;%s/p_{T}'%isoTitle[isoKey], 5,array('d',zcenq), 30,-5,5)
    huncorbkg=ROOT.TH2F('bkgcenvsuncor%s'%isoKey, ';Centrality;%s/p_{T}'%isoTitle[isoKey], 5,array('d',zcenq), 30,-5,5)
    hcor=ROOT.TH2F('cenvscor%s'%isoKey,           ';Centrality;%s/p_{T}'%isoTitle[isoKey], 5,array('d',zcenq), 30,-5,5)
    hcorbkg=ROOT.TH2F('bkgcenvscor%s'%isoKey,     ';Centrality;%s/p_{T}'%isoTitle[isoKey], 5,array('d',zcenq), 30,-5,5)
    for i in range(len(leptons)):
        l,isZ,isSS,cenbin=leptons[i]
        pt=l.p4.Pt()
        iso=isoVals[i]
        corIso=iso-cenProf.Eval(cenbin)
        if isZ:
            huncor.Fill(cenbin,iso/pt)
            hcor.Fill(cenbin,corIso/pt)
        if isSS:
            huncorbkg.Fill(cenbin,iso/pt)
            hcorbkg.Fill(cenbin,corIso/pt)

    for h1,h2,f in [(h,hbkg,cenProf),
                    (huncor,huncorbkg,None),
                    (hcor,hcorbkg,None)]:

        c.SetRightMargin(0.03)
        h1.SetBarWidth(0.4);
        h1.SetBarOffset(-0.25);
        h2.SetBarWidth(0.4);
        h2.SetBarOffset(0.25);
        h2.SetLineColor(ROOT.kGreen+2);
        h1.Draw("candle1");
        h2.Draw("candle1 same");
        if f : f.Draw('same')
        leg=ROOT.TLegend(0.73,0.94,0.93,0.8)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        leg.SetFillStyle(0)
        h1.SetTitle('Z#rightarrowee')
        leg.AddEntry(h1,h1.GetTitle(),'lp')
        h2.SetTitle('same-sign ee')
        leg.AddEntry(h2,h2.GetTitle(),'lp')
        leg.Draw()
        extraTxt=[]
        if f: extraTxt=['I=exp(%3.3fc)[%3.3fc^{2}+%3.3fc+%3.3f]'%(f.GetParameter(0),f.GetParameter(1),f.GetParameter(2),f.GetParameter(3))]
        header(extraTxt=extraTxt)
        c.Modified()
        c.Update()
        c.RedrawAxis()
        c.SaveAs(h1.GetName()+'.png')

        if f: continue

        tag='' if 'uncor' in h1.GetName() else '*'

        rocs.append( getROC(h1.ProjectionY('sig'),h2.ProjectionY('bkg'),isoKey) )
        ci=len(rocs)%5+1
        ls=len(rocs)/4+1
        rocs[-1].SetLineColor(ci)
        rocs[-1].SetMarkerColor(ci)
        rocs[-1].SetLineWidth(2)
        rocs[-1].SetLineStyle(ls)
        rocs[-1].SetTitle(isoTitle[isoKey]+tag)
        rocshcen.append( getROC(h1.ProjectionY('sig',1,2),h2.ProjectionY('bkg',1,2),isoKey+'hcen') )
        rocshcen[-1].SetLineColor(ci)
        rocshcen[-1].SetMarkerColor(ci)
        rocshcen[-1].SetLineWidth(2)
        rocshcen[-1].SetLineStyle(ls)
        rocshcen[-1].SetTitle(isoTitle[isoKey]+tag)
    

import numpy as np
rocratio=[]
for i in range(len(rocs)):
    rocratio.append(rocs[i].Clone())
    rocratio[-1].SetName('rocratio%d'%i)
    rocratio[-1].Set(0)
    for eff in np.arange(0,1,0.02):
        den=rocs[i].Eval(eff)
        if den<=0 : continue
        num=rocshcen[i].Eval(eff)
        if num<=0: continue
        rocratio[-1].SetPoint(rocratio[-1].GetN(),eff,num/den)
c.SetGridy()
c.SetRightMargin(0.03)
for rocColl, rocName in [(rocs,'rocs'),
                         (rocshcen,'rocshcen'),
                         (rocratio,'rocsratio')]:

    mg=ROOT.TMultiGraph()
    for g in rocColl: mg.Add(g,'l')
    mg.Draw('al')
    mg.GetXaxis().SetRangeUser(0,1)
    mg.GetYaxis().SetRangeUser(0,1)
    mg.GetXaxis().SetTitle('Signal efficiency')
    mg.GetYaxis().SetTitle('Background efficiency')
    mg.GetXaxis().SetRangeUser(0,1)
    mg.GetYaxis().SetRangeUser(0,1)
    if 'ratio' in rocName:
        mg.GetYaxis().SetTitle('Central / Inclusive')
        mg.GetYaxis().SetRangeUser(0.5,2.5)
        leg=c.BuildLegend(0.65,0.94,0.92,0.6)
        l=ROOT.TLine(0.9,0.5,0.9,2.5)
        l.SetLineColor(ROOT.kGray)
        l.Draw()
    else:
        leg=c.BuildLegend(0.15,0.94,0.4,0.6)
        l=ROOT.TLine(0.9,0,0.9,1)
        l.SetLineColor(ROOT.kGray)
        l.Draw()
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.SetFillStyle(0)
   
    header()
    c.Modified()
    c.Update()
    c.SaveAs("iso%s.png"%rocName)

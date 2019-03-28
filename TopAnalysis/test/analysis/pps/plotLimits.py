import ROOT
import pickle
import sys

url=sys.argv[1]
id=int(sys.argv[2])
with open(url,'r') as cache:
    limits=[x for x in pickle.load(cache) if x[1]==id]
lumi=37.5 if id==23 else 2.64
searchTitle='pp#rightarrowppZX' if id==23 else 'pp#rightarrowpp#gammaX'


r68=ROOT.TGraphAsymmErrors()
r68.SetFillColor(8)
r68.SetFillStyle(1001)
r68.SetMarkerColor(8)
r68.SetLineColor(1)
r68.SetLineWidth(2)
r68.SetLineStyle(7)
r68.SetName('r68')
r68.SetTitle('Expected (68%)')
r95=ROOT.TGraphAsymmErrors()
r95.SetFillColor(5)
r95.SetFillStyle(1001)
r95.SetMarkerColor(5)
r95.SetLineColor(1)
r95.SetLineWidth(2)
r95.SetLineStyle(7)
r95.SetName('r95')
r95.SetTitle('Expected (95%)')
rmed=ROOT.TGraph()
rmed.SetFillColor(0)
rmed.SetFillStyle(0)
rmed.SetLineColor(1)
rmed.SetLineWidth(2)
rmed.SetLineStyle(7)
sig=rmed.Clone('sig')

for vals in limits:

    m=vals[0]
    np=rmed.GetN()
    rmed.SetPoint(np,m,vals[4])
    r68.SetPoint(np,m,vals[4])
    r68.SetPointError(np,10,10,vals[4]-vals[3],vals[5]-vals[3])
    r95.SetPoint(np,m,vals[4])
    r95.SetPointError(np,10,10,vals[4]-vals[2],vals[6]-vals[4])
    sig.SetPoint(np,m,vals[7])

rmed.Sort()
r68.Sort()
r95.Sort()
sig.Sort()

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.03)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.1)
mg=ROOT.TMultiGraph()
mg.Add(r95,'l3')
mg.Add(r68,'l3')
mg.Add(rmed,'l')
frame=ROOT.TH1F('frame',';m_{X} [GeV];95% CL limits on #sigma/#sigma_{fid}',1,700,1700)
frame.GetYaxis().SetRangeUser(0,2)
frame.SetBinContent(1,1)
frame.SetLineWidth(2)
frame.SetLineColor(ROOT.kRed)
frame.Draw('hist')
mg.Draw('l3')
leg=ROOT.TLegend(0.15,0.92,0.4,0.72)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.AddEntry(r68,r68.GetTitle(),'lf')
leg.AddEntry(r95,r95.GetTitle(),'lf')
leg.Draw()
tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.04)
tex.SetNDC()
tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{preliminary}')
tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
tex.DrawLatex(0.97,0.975,'#scale[0.9]{%s fb^{-1} (13 TeV)}'%lumi)
tex.DrawLatex(0.95,0.9,searchTitle)
c.Modified()
c.Update()   
c.RedrawAxis()
for ext in ['png','pdf']:
    c.SaveAs('limits_%d_hpur.%s'%(id,ext))

c.Clear()
frame.Draw()
frame.GetYaxis().SetTitle('Signal significance (asymptotic)')
frame.SetBinContent(1,5)
frame.GetYaxis().SetRangeUser(0,10)
sig.Draw('l')
tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.04)
tex.SetNDC()
tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{preliminary}')
tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
tex.DrawLatex(0.97,0.975,'#scale[0.9]{%s fb^{-1} (13 TeV)}'%lumi)
tex.DrawLatex(0.95,0.9,searchTitle)
c.Modified()
c.Update()   
c.RedrawAxis()
for ext in ['png','pdf']:
    c.SaveAs('significance_%d_hpur.%s'%(id,ext))
import ROOT as r

r.gROOT.ProcessLine(".L functions.cc+" )

def deactivateBranches(tree):
    tree.SetBranchStatus('*', 0)
    tree.SetBranchStatus('branchname'     , 1)
    return tree

r.TMVA.Tools.Instance()
 
# note that it seems to be mandatory to have an
# output file, just passing None to TMVA::Factory(..)
# does not work. Make sure you don't overwrite an
# existing file.

output_fn  = 'training.root'
output_f   = r.TFile(output_fn,'RECREATE')
 
factory = r.TMVA.Factory('TMVAClassification', output_f,
                         ':'.join([
                         '!V',
                         '!Silent',
                         'Color',
                         'DrawProgressBar',
                         'Transformations=I;D;P;G,D',
                         'AnalysisType=Classification'])  )

dataloader = r.TMVA.DataLoader('trainingV2_sixBestNew')

dataloader.AddVariable('lep_pt[0]'  , 'p_{T}^{lep1}'          , 'GeV' , 'F')
dataloader.AddVariable('apt'        , 'A_{pt}'                , ''    , 'F')
dataloader.AddVariable('llpt'       , 'p_{T}^{ll}'            , 'GeV' , 'F')
dataloader.AddVariable('abs(lleta)' , '|#eta^{ll}|'           , ''    , 'F')
dataloader.AddVariable('dphi'       , '|#Delta #phi|'         , 'rad' , 'F')
dataloader.AddVariable('dphilll2'   , '|#Delta #phi_{ll,l2}|' , 'rad' , 'F')

##  dataloader.AddSpectator('llm'                                         , 'spec m_{ll}'         , 'GeV' )#, 'F')

#dataloader.AddVariable('abs(dphi_2(lep_pt[0],lep_eta[0],lep_phi[0],lep_pt[1],lep_eta[1],lep_phi[1],1))', '|#Delta #phi_{ll,l1}|'    , 'rad'    , 'F')
#dataloader.AddVariable('abs(sumeta)'                          , '|#sum #eta_{i}|', ''    , 'F')
#dataloader.AddVariable('abs(lep_eta[0])+abs(lep_eta[1])'                          , '#sum |#eta_{i}|', ''    , 'F')
#dataloader.AddVariable('max(abs(lep_eta[0]),abs(lep_eta[1]))'                          , 'max(|#eta_{i}|)', ''    , 'F')
#dataloader.AddVariable('deta'                                        , '#Delta #eta'    , ''    , 'F')
#dataloader.AddVariable('lep_pt[1]'                                   , 'p_{T}^{lep2}'   , 'GeV' , 'F')
 
## add weights for signal and background
dataloader.SetBackgroundWeightExpression('weight/abs(weight)')
dataloader.SetSignalWeightExpression('weight/abs(weight)')

## get background tree and friends etc
bkg_tfile = r.TFile('../plots/DYJetsToLL_MLL-50_TuneCP5_5020GeV-amcatnloFXFX-pythia8.root')
bkg_tree = bkg_tfile.Get('tree')

## get signal tree and friends etc
sig_tfile = r.TFile('../plots/TT_TuneCP5_5p02TeV-powheg-pythia8.root')
sig_tree = sig_tfile.Get('tree')

dataloader.AddSignalTree    (sig_tree)
dataloader.AddBackgroundTree(bkg_tree)

# cuts defining the signal and background sample
sig_cutstring = 'llm > 106. || llm < 76.'
bkg_cutstring = 'llm > 106. || llm < 76.'
sigCut = r.TCut(sig_cutstring)
bgCut  = r.TCut(bkg_cutstring)

## just some histograms if you want to. not really necessary
h_ahisto = r.TH1F('h_ahisto', 'h_ahisto', 100, 0., 3.1416)
h_bhisto = r.TH1F('h_bhisto', 'h_bhisto', 100, 0., 3.1416)

sig_tree.Draw('dphi>>h_ahisto', sig_cutstring+'*weight/abs(weight)')
bkg_tree.Draw('dphi>>h_bhisto', bkg_cutstring+'*weight/abs(weight)')
h_ahisto.Scale(1./h_ahisto.Integral())
h_bhisto.Scale(1./h_bhisto.Integral())

#h_dphiROC = r.TH2F('h_dphiROC', 'ROC of dPhi', 100, 0., 1., 100, 0., 1.)
h_dphiROC = r.TH1F('h_dphiROC', 'ROC of dPhi', 100, 0., 1.)

for ib in range(1,h_ahisto.GetNbinsX()+1):
    sigeff = h_ahisto.Integral(1,ib)/h_ahisto.Integral()
    bkgeff = h_bhisto.Integral(1,ib)/h_bhisto.Integral()
    h_dphiROC.SetBinContent(h_dphiROC.GetXaxis().FindBin(sigeff), 1.-bkgeff)
    h_dphiROC.SetBinError  (h_dphiROC.GetXaxis().FindBin(sigeff), 0.001)

r.gStyle.SetOptStat(0)
c1 = r.TCanvas()
h_dphiROC.SetMarkerSize(1.2)
h_dphiROC.SetLineColor(r.kPink)
h_dphiROC.GetYaxis().SetRangeUser(0.2, 1.)
h_dphiROC.GetYaxis().SetTitle('backround rejection')
h_dphiROC.GetXaxis().SetTitle('signal efficiency')
h_dphiROC.Draw('l')

 
dataloader.PrepareTrainingAndTestTree(sigCut,   # signal events
                                   bgCut,    # background events
                                   ':'.join([
                                   'nTrain_Signal=0',
                                   'nTrain_Background=0',
                                   'SplitMode=Random',
                                   'NormMode=NumEvents',
                                   '!V' ]))

## define your methods: BDT, FISHER, LIKELIHOOD. along with the options

### some sort of bdt
bdt = factory.BookMethod(dataloader, r.TMVA.Types.kBDT, 'BDT', 
                         ':'.join([ '!H', 
                                    '!V', 
                                    'NTrees=850', 
                                    'CreateMVAPdfs',
                                    'MinNodeSize=0.05', 
                                    'MaxDepth=3', 
                                    'BoostType=AdaBoost', 
                                    'AdaBoostBeta=0.5', 
                                    'SeparationType=GiniIndex', 
                                    'nCuts=20', 
                                    'PruneMethod=NoPruning', ]))

bdtg = factory.BookMethod(dataloader, r.TMVA.Types.kBDT, 'BDTG', "!H:!V:NTrees=1000:CreateMVAPdfs:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

# Fisher discriminant (same as LD)
fisher = factory.BookMethod(dataloader, r.TMVA.Types.kFisher, 'Fisher', 
                            ':'.join(['H',
                                      '!V',
                                      'Fisher',
                                      'CreateMVAPdfs',
                                      'PDFInterpolMVAPdf=Spline2',
                                      'NbinsMVAPdf=50',
                                      'NsmoothMVAPdf=10']) )

## likelihood
lh = factory.BookMethod( dataloader, r.TMVA.Types.kLikelihood, 'LikelihoodD', 
                        ':'.join([ '!H', 
                                   '!V', 
                                   'CreateMVAPdfs',
                                   '!TRansformOutput', 
                                   'PDFInterpol=Spline2', 
                                   'NSmoothSig[0]=20', 
                                   'NSmooth=5', 
                                   'NAvEvtPerBin=50', 
                                   'VarTransform=Decorrelate' ]))

## ## cuts
## cuts = factory.BookMethod( dataloader, r.TMVA.Types.kCuts, 'CutsD',
##                            ':'.join(['!H', 
##                                      '!V', 
##                                      'FitMethod=MC', 
##                                      'EffSel', 
##                                      'SampleSize=200000',
##                                      'VarProp=FSmart',
##                                      'VarTransform=Decorrelate']) );
 
## do the training
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

output_f.Close()

r.TMVA.TMVAGui(output_fn)


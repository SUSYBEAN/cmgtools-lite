import ROOT

from CMGTools.XZZ2l2nu.plotting.TreePlotter import TreePlotter
from CMGTools.XZZ2l2nu.plotting.MergedPlotter import MergedPlotter
from CMGTools.XZZ2l2nu.plotting.StackPlotter import StackPlotter

ROOT.gROOT.ProcessLine('.x tdrstyle.C') 
tag='mc_'
cuts='(nllnunu)'
lumi=2.1
sepSig=True

vvPlotters=[]
vvSamples = ['WW','WZ','ZZ']

for sample in vvSamples:
    vvPlotters.append(TreePlotter(sample+'.root','tree'))
    vvPlotters[-1].setupFromFile(sample+'.pck')
    vvPlotters[-1].addCorrectionFactor('xsec','tree')
    vvPlotters[-1].addCorrectionFactor('genWeight','tree')
    vvPlotters[-1].addCorrectionFactor('puWeight','tree')

VV = MergedPlotter(vvPlotters)
VV.setFillProperties(1001,ROOT.kAzure-9)

zjetsPlotters=[]
zjetsSamples = ['DYJetsToLL_M50_HT100to200','DYJetsToLL_M50_HT200to400','DYJetsToLL_M50_HT400to600','DYJetsToLL_M50_HT600toInf']
zjetsSampleNames = ['ZJets1','ZJets2','ZJets3','ZJets4']

for sample in zjetsSamples:
    zjetsPlotters.append(TreePlotter(sample+'.root','tree'))
    zjetsPlotters[-1].setupFromFile(sample+'.pck')
    zjetsPlotters[-1].addCorrectionFactor('xsec','tree')
    zjetsPlotters[-1].addCorrectionFactor('genWeight','tree')
    zjetsPlotters[-1].addCorrectionFactor('puWeight','tree')

ZJets = MergedPlotter(zjetsPlotters)
ZJets.setFillProperties(1001,ROOT.kGreen+2)


sigPlotters=[]
sigSamples = ['RSGravToZZToZZinv_narrow_M-1000', 'RSGravToZZToZZinv_narrow_M-1400', 'RSGravToZZToZZinv_narrow_M-2000']
sigSampleNames = ['RS1000','RS1400','RS2000']

for sample in sigSamples:
    sigPlotters.append(TreePlotter(sample+'.root','tree'))
    sigPlotters[-1].setupFromFile(sample+'.pck')
    sigPlotters[-1].addCorrectionFactor('xsec','tree')
    sigPlotters[-1].addCorrectionFactor('genWeight','tree')
    sigPlotters[-1].addCorrectionFactor('puWeight','tree')
    sigPlotters[-1].setFillProperties(0,ROOT.kWhite)



mcStack = StackPlotter()
mcStack.addPlotter(VV, "WWWZZZ","WW/WZ/ZZ", "background")
mcStack.addPlotter(ZJets, "ZJets","Z+Jets", "background")

for i in range(len(sigSamples)):
  sigPlotters[i].setLineProperties(2,ROOT.kRed+i,2)
  mcStack.addPlotter(sigPlotters[i],sigSampleNames[i],sigSampleNames[i],'signal')  
  
 
mcStack.setLog(True)

mcStack.drawStack('llnunu_mt', cuts, str(lumi*1000), 80, 0.0, 2400.0, titlex = "M_{T}", units = "GeV",output=tag+'out_mt_stack.eps',separateSignal=sepSig)
mcStack.drawStack('llnunu_l2_pt', cuts, str(lumi*1000), 50, 0.0, 1500.0, titlex = "MET", units = "GeV",output=tag+'out_met_stack.eps',separateSignal=sepSig)
mcStack.drawStack('llnunu_l1_pt', cuts, str(lumi*1000), 50, 0.0, 1500.0, titlex = "P_{T}(Z)", units = "GeV",output=tag+'out_zpt_stack.eps',separateSignal=sepSig)


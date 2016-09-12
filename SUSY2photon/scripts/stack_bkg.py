import ROOT

from CMGTools.XZZ2l2nu.plotting.TreePlotter import TreePlotter
from CMGTools.XZZ2l2nu.plotting.MergedPlotter import MergedPlotter
from CMGTools.XZZ2l2nu.plotting.StackPlotter import StackPlotter

ROOT.gROOT.ProcessLine('.x tdrstyle.C') 
vvPlotters=[]
vvSamples = ['WW','WZ','ZZ']
tag='bkg_'
cuts='(nllnunu)'
lumi=2.318278305

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


bkgStack = StackPlotter()
bkgStack.addPlotter(VV, "WWWZZZ","WW/WZ/ZZ", "background")
bkgStack.addPlotter(ZJets, "ZJets","Z+Jets", "background")
 
bkgStack.setLog(True)

bkgStack.drawStack('llnunu_mt', cuts, str(lumi*1000), 100, 0.0, 2000.0, titlex = "M_{T}", units = "GeV",output=tag+'out_mt_stack.eps')
bkgStack.drawStack('llnunu_l2_pt', cuts, str(lumi*1000), 100, 0.0, 2000.0, titlex = "MET", units = "GeV",output=tag+'out_met_stack.eps')
bkgStack.drawStack('llnunu_l1_pt', cuts, str(lumi*1000), 100, 0.0, 2000.0, titlex = "P_{T}(Z)", units = "GeV",output=tag+'out_zpt_stack.eps')


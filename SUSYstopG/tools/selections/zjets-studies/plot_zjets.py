#!/usr/bin/env python

import ROOT
import os
from math import *
from python.TreePlotter import TreePlotter
from python.MergedPlotter import MergedPlotter
from python.StackPlotter import StackPlotter
from python.SetCuts import SetCuts

Channel=raw_input("Please choose a channel (el or mu): \n")
tag0='ZJstudy'
outdir='test'
indir="../AnalysisRegion"
lumi=2.169126704526
LogY=False

if not os.path.exists(outdir): os.system('mkdir '+outdir)

tag = tag0+'_'+'test'
if LogY: tag = tag+'_log'
outTag=outdir+'/'+tag

#  A | C
# ------- dPhi(Z,MET) = 2.5
#  B | D
# dR(ll)=0.8 / 1.0
mycuts=SetCuts()
tex_dic=mycuts.Tex_dic
cuts=mycuts.abcdCuts(Channel)
print cuts

### ----- Initialize (samples):

zjetsPlotters=[]
zjetsSamples = ['DYJetsToLL_M50'] # M50_BIG = M50 + M50_Ext, 150M evts

for sample in zjetsSamples:
    zjetsPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    zjetsPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
    zjetsPlotters[-1].addCorrectionFactor('xsec','tree')
    zjetsPlotters[-1].addCorrectionFactor('genWeight','tree')
    zjetsPlotters[-1].addCorrectionFactor('puWeight','tree')

ZJets = MergedPlotter(zjetsPlotters)
ZJets.setFillProperties(1001,ROOT.kGreen+2)

dataPlotters=[]
dataSamples = ['SingleElectron_Run2015C_25ns_16Dec',
               'SingleElectron_Run2015D_16Dec',
               'SingleMuon_Run2015C_25ns_16Dec',
               'SingleMuon_Run2015D_16Dec']
for sample in dataSamples:
    dataPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    
Data = MergedPlotter(dataPlotters)

otherBGPlotters=[]
otherBGSamples=['WZTo2L2Q','WZTo3LNu',
                 'ZZTo2L2Nu','ZZTo2L2Q','ZZTo4L',
                 'TTTo2L2Nu',
                 'WWTo2L2Nu','WWToLNuQQ','WZTo1L1Nu2Q']
for sample in otherBGSamples:
    otherBGPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    otherBGPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
    otherBGPlotters[-1].addCorrectionFactor('xsec','tree')
    otherBGPlotters[-1].addCorrectionFactor('genWeight','tree')
    otherBGPlotters[-1].addCorrectionFactor('puWeight','tree')
    

### ----- Execute (plotting):

ROOT.gROOT.ProcessLine('.x tdrstyle.C')
ROOT.gStyle.SetPadBottomMargin(0.2)
ROOT.gStyle.SetPadLeftMargin(0.15)

#ROOT.TH1.AddDirectory(ROOT.kFALSE) #in this way you could close the TFile after you registe the histograms
yields={}
err={}
nan=0
fout = ROOT.TFile(outTag+'.root','recreate')

for key in tex_dic:
    canvas = ROOT.TCanvas('c1', 'c1', 600,630)
    canvas.Print(outTag+'.ps[')

    # MET:
    h_met=ZJets.drawTH1('met_pt',cuts[key],str(lumi*1000),50,0,5,titlex='#Delta R(#mu,#mu)',units='',drawStyle="HIST")
    err[key]=ROOT.Double(0.0)
    yields[key]=h_met.IntegralAndError(0,1+h_met.GetNbinsX(),err[key])
    
    # deltaR(mu,mu)
    hdRZmm=ZJets.drawTH1('llnunu_l1_deltaR',cuts[key],str(lumi*1000),50,0,5,titlex='#Delta R(#mu,#mu)',units='',drawStyle="HIST")
    hdRZmm.SetName('hdRZmm')
    hdRZmm.GetYaxis().SetTitle("Events")
    
    canvas.Clear()
    hdRZmm.Draw()
    canvas.Print(outTag+'.ps')
    canvas.Clear()
    
    hdRZmm.Write()
    
    #deltaPhi(Zmumu,MET)
    hdPhiZmmMet=ZJets.drawTH1('llnunu_deltaPhi',cuts[key],str(lumi*1000),50,0,5,titlex='#Delta#Phi(Z_{#mu,#mu},MET)',units='',drawStyle="HIST")
    hdPhiZmmMet.SetName('hdPhiZmmMet')
    hdPhiZmmMet.GetYaxis().SetTitle("Events")
    
    canvas.Clear()
    hdPhiZmmMet.Draw()
    canvas.Print(outTag+'.ps]')
    #canvas.Clear()
    
    hdPhiZmmMet.Write()
    os.system('ps2pdf '+outTag+'.ps '+outTag+'.pdf')    

print yields
print err

### ----- Finalizing:

fout.Close()



#!/usr/bin/env python

import ROOT
import os
from math import *
from CMGTools.XZZ2l2nu.plotting.TreePlotter import TreePlotter
from CMGTools.XZZ2l2nu.plotting.MergedPlotter import MergedPlotter
from CMGTools.XZZ2l2nu.plotting.StackPlotter import StackPlotter

cutChain='zjetscut'
#cutChain='zjetscutmet50'

tag0='ZJstudy'
outdir='plots'
indir="/afs/cern.ch/work/h/heli/public/XZZ/76X"
lumi=2.169126704526
LogY=False

if not os.path.exists(outdir): os.system('mkdir '+outdir)

tag = tag0+'_'+cutChain
if LogY: tag = tag+'_log'

outTag=outdir+'/'+tag

cuts_zjets='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(abs(llnunu_l1_mass-91.18)<20))'
cuts_zjets_met50='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(abs(llnunu_l1_mass-91.18)<20&&met_pt>50))'

if cutChain=='zjetscut': cuts=cuts_zjets
elif cutChain=='zjetscutmet50': cuts=cuts_zjets_met50
else : cuts=cuts_loose


isEl='(abs(llnunu_l1_l1_pdgId)==11)'
isMu='(abs(llnunu_l1_l1_pdgId)==13)'

#  A | B
# ------- dPhi(Z,MET) = 2.5
#  C | D
# dR(ll)=0.8

RgA ='(llnunu_l1_deltaR<0.8&&llnunu_deltaPhi>2.5)'
RgB ='(llnunu_l1_deltaR>=0.8&&llnunu_deltaPhi>2.5)'
RgC ='(llnunu_l1_deltaR<0.8&&llnunu_deltaPhi<=2.5)'
RgD ='(llnunu_l1_deltaR>=0.8&&llnunu_deltaPhi<=2.5)'

ROOT.gROOT.ProcessLine('.x tdrstyle.C') 

ROOT.gStyle.SetPadBottomMargin(0.2)
ROOT.gStyle.SetPadLeftMargin(0.15)

zjetsPlotters=[]
zjetsSamples = ['DYJetsToLL_M50_BIG'] # M50_BIG = M50 + M50_Ext
#zjetsSamples = ['DYJetsToLL_M50'] # M50_BIG = M50 + M50_Ext


for sample in zjetsSamples:
    zjetsPlotters.append(TreePlotter(indir+'/'+sample+'/vvTreeProducer/tree.root','tree'))
    zjetsPlotters[-1].setupFromFile(indir+'/'+sample+'/skimAnalyzerCount/SkimReport.pck')
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
    dataPlotters.append(TreePlotter(indir+'/'+sample+'/vvTreeProducer/tree.root','tree'))

Data = MergedPlotter(dataPlotters)


fout = ROOT.TFile(outTag+'.root','recreate')
canvas = ROOT.TCanvas('c1', 'c1', 600,630)

canvas.Print(outTag+'.ps[')

# 1d variables

#deltaR(mu,mu)
hdRZmm=ZJets.drawTH1('llnunu_l1_deltaR',cuts+'&&(abs(llnunu_l1_l1_pdgId)==13)',str(lumi*1000),50,0,5,titlex='#Delta R(#mu,#mu)',units='',drawStyle="HIST")
hdRZmm.SetName('hdRZmm')
hdRZmm.GetYaxis().SetTitle("Events")

canvas.Clear()
hdRZmm.Draw()
canvas.Print(outTag+'.ps')
canvas.Clear()

hdRZmm.Write()


#deltaPhi(Zmumu,MET)
hdPhiZmmMet=ZJets.drawTH1('llnunu_deltaPhi',cuts+'&&(abs(llnunu_l1_l1_pdgId)==13)',str(lumi*1000),50,0,5,titlex='#Delta#Phi(Z_{#mu,#mu},MET)',units='',drawStyle="HIST")
hdPhiZmmMet.SetName('hdPhiZmmMet')
hdPhiZmmMet.GetYaxis().SetTitle("Events")

canvas.Clear()
hdPhiZmmMet.Draw()
canvas.Print(outTag+'.ps')
canvas.Clear()

hdPhiZmmMet.Write()


#deltaR(e,e)
hdRZee=ZJets.drawTH1('llnunu_l1_deltaR',cuts+'&&(abs(llnunu_l1_l1_pdgId)==11)',str(lumi*1000),50,0,5,titlex='#Delta R(e,e)',units='',drawStyle="HIST")
hdRZee.SetName('hdRZmm')
hdRZee.GetYaxis().SetTitle("Events")

canvas.Clear()
hdRZee.Draw()
canvas.Print(outTag+'.ps')
canvas.Clear()

hdRZee.Write()


#deltaPhi(Zee,MET)
hdPhiZeeMet=ZJets.drawTH1('llnunu_deltaPhi',cuts+'&&(abs(llnunu_l1_l1_pdgId)==11)',str(lumi*1000),50,0,5,titlex='#Delta#Phi(Z_{e,e},MET)',units='',drawStyle="HIST")
hdPhiZeeMet.SetName('hdPhiZeeMet')
hdPhiZeeMet.GetYaxis().SetTitle("Events")

canvas.Clear()
hdPhiZeeMet.Draw()
canvas.Print(outTag+'.ps')
canvas.Clear()

hdPhiZeeMet.Write()

# 4 METs

# mu
hMETAmu=ZJets.drawTH1('met_pt',cuts+'&&'+isMu+'&&'+RgA,str(lumi*1000),50,0,1000,titlex='MET',units='GeV',drawStyle="HIST")
hMETAmu.SetName("hMETAmu")
hMETBmu=ZJets.drawTH1('met_pt',cuts+'&&'+isMu+'&&'+RgB,str(lumi*1000),50,0,1000,titlex='MET',units='GeV',drawStyle="HIST")
hMETBmu.SetName("hMETBmu")
hMETCmu=ZJets.drawTH1('met_pt',cuts+'&&'+isMu+'&&'+RgC,str(lumi*1000),50,0,1000,titlex='MET',units='GeV',drawStyle="HIST")
hMETCmu.SetName("hMETCmu")
hMETDmu=ZJets.drawTH1('met_pt',cuts+'&&'+isMu+'&&'+RgD,str(lumi*1000),50,0,1000,titlex='MET',units='GeV',drawStyle="HIST")
hMETDmu.SetName("hMETDmu")

hMETAmu.GetYaxis().SetTitle("Events")
hMETBmu.GetYaxis().SetTitle("Events")
hMETCmu.GetYaxis().SetTitle("Events")
hMETDmu.GetYaxis().SetTitle("Events")

# integral and error
nMETAmuErr=ROOT.Double(0.0)
nMETBmuErr=ROOT.Double(0.0)
nMETCmuErr=ROOT.Double(0.0)
nMETDmuErr=ROOT.Double(0.0)
nMETAmu = hMETAmu.IntegralAndError(1,hMETAmu.GetNbinsX(),nMETAmuErr)
nMETBmu = hMETBmu.IntegralAndError(1,hMETBmu.GetNbinsX(),nMETBmuErr)
nMETCmu = hMETCmu.IntegralAndError(1,hMETCmu.GetNbinsX(),nMETCmuErr)
nMETDmu = hMETDmu.IntegralAndError(1,hMETDmu.GetNbinsX(),nMETDmuErr)

rMETmuAovC = nMETAmu/nMETCmu
rMETmuBovD = nMETBmu/nMETDmu
rMETmuAovCerr = sqrt((nMETCmu**2*nMETAmuErr**2+nMETAmu**2*nMETCmuErr**2)/nMETCmu**4)
rMETmuBovDerr = sqrt((nMETDmu**2*nMETBmuErr**2+nMETBmu**2*nMETDmuErr**2)/nMETDmu**4)

print '#####  Muon ######'
print ' A/C  | '+str(rMETmuAovC)+' +- '+str(rMETmuAovCerr)
print ' B/D  | '+str(rMETmuBovD)+' +- '+str(rMETmuBovDerr)
print '##################'

# predicted A
hMETAmuPred = hMETCmu.Clone("hMETAmuPred")
hMETAmuPred.Scale(rMETmuBovD)
hMETAmuPred.SetLineColor(2)
hMETAmuPred.SetMarkerColor(2)

#legend
lgMETAmu=ROOT.TLegend(0.7,0.7,0.85,0.86)
lgMETBmu=ROOT.TLegend(0.7,0.7,0.85,0.86)
lgMETCmu=ROOT.TLegend(0.7,0.7,0.85,0.86)
lgMETDmu=ROOT.TLegend(0.7,0.7,0.85,0.86)
lgMETAmu.SetName('lgMETAmu')
lgMETBmu.SetName('lgMETBmu')
lgMETCmu.SetName('lgMETCmu')
lgMETDmu.SetName('lgMETDmu')

lgMETAmu.AddEntry(hMETAmu, 'A (signal)', 'pl')
lgMETAmu.AddEntry(hMETAmuPred, 'Predicted', 'pl')
lgMETBmu.AddEntry(hMETAmu, 'B ', 'pl')
lgMETCmu.AddEntry(hMETBmu, 'C ', 'pl')
lgMETDmu.AddEntry(hMETDmu, 'D', 'pl')

# write
hMETAmu.Write()
hMETBmu.Write()
hMETCmu.Write()
hMETDmu.Write()
hMETAmuPred.Write()
lgMETAmu.Write()
lgMETBmu.Write()
lgMETCmu.Write()
lgMETDmu.Write()

# plot
canvas.Clear()
canvas.Divide(2,2)
canvas.cd(1)
hMETAmu.Draw()
hMETAmuPred.Draw('same')
lgMETAmu.Draw()
canvas.cd(2)
hMETBmu.Draw()
lgMETBmu.Draw()
canvas.cd(3)
hMETCmu.Draw()
lgMETCmu.Draw()
canvas.cd(4)
hMETDmu.Draw()
lgMETDmu.Draw()
canvas.Print(outTag+'.ps')
canvas.Clear()

# one by one
canvas.Clear()
hMETAmu.Draw()
hMETAmuPred.Draw('same')
lgMETAmu.Draw()
canvas.Print(outTag+'.ps')
canvas.Clear()

hMETBmu.Draw()
lgMETBmu.Draw()
canvas.Print(outTag+'.ps')
canvas.Clear()

hMETCmu.Draw()
lgMETCmu.Draw()
canvas.Print(outTag+'.ps')
canvas.Clear()

hMETDmu.Draw()
lgMETDmu.Draw()
canvas.Print(outTag+'.ps')
canvas.Clear()


# log plots
# plot
canvas.Clear()
canvas.SetLogy()
canvas.Divide(2,2)
canvas.cd(1)
hMETAmu.Draw()
hMETAmuPred.Draw('same')
lgMETAmu.Draw()
canvas.cd(2)
hMETBmu.Draw()
lgMETBmu.Draw()
canvas.cd(3)
hMETCmu.Draw()
lgMETCmu.Draw()
canvas.cd(4)
hMETDmu.Draw()
lgMETDmu.Draw()
canvas.Print(outTag+'.ps')
canvas.Clear()

# one by one
canvas.Clear()
hMETAmu.Draw()
hMETAmuPred.Draw('same')
lgMETAmu.Draw()
canvas.Print(outTag+'.ps')
canvas.Clear()

hMETBmu.Draw()
lgMETBmu.Draw()
canvas.Print(outTag+'.ps')
canvas.Clear()

hMETCmu.Draw()
lgMETCmu.Draw()
canvas.Print(outTag+'.ps')
canvas.Clear()

hMETDmu.Draw()
lgMETDmu.Draw()
canvas.Print(outTag+'.ps')
canvas.Clear()
canvas.SetLogy(False)


canvas.Print(outTag+'.ps]')
os.system('ps2pdf '+outTag+'.ps '+outTag+'.pdf')
fout.Close()



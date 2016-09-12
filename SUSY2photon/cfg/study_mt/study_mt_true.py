#!/usr/bin/env python

import ROOT
import os
from CMGTools.XZZ2l2nu.plotting.TreePlotter import TreePlotter
from CMGTools.XZZ2l2nu.plotting.MergedPlotter import MergedPlotter
from CMGTools.XZZ2l2nu.plotting.StackPlotter import StackPlotter

cutChain='loosecut'
#cutChain='tightzpt100'
#cutChain='tightzpt100met100'
#cutChain='tightmz15zpt225met165'
#cutChain='tightzpt100met100dphi'
#cutChain='tightmet250dphi'
#cutChain='zjetscut'
#cutChain='zjetscutmet50'

channel='both' # can be el or mu or both

sigXsecOne=True

outdir='plots'
indir="./"
lumi=2.169126704526
sepSig=True
LogY=True
DrawLeptons=False
doRatio=True
test=True
Blind=False
FakeData=False

elChannel='(abs(llnunu_l1_l1_pdgId)==11)'
muChannel='(abs(llnunu_l1_l1_pdgId)==13)'

if not os.path.exists(outdir): os.system('mkdir '+outdir)

tag = cutChain+'_'
tag = tag+channel+'_'
if LogY: tag = tag+'log_'
if sigXsecOne: tag=tag+'sig1_'

cuts_loose='(nllnunu)'
cuts_loose_zll='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>100&&llnunu_l1_mass>70&&llnunu_l1_mass<110))'
cuts_loose_zll_met100='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>100&&llnunu_l1_mass>70&&llnunu_l1_mass<110&&met_pt>100))'
cuts_loose_zll_mz15_zpt225_met165='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>225&&llnunu_l1_mass>76.2&&llnunu_l1_mass<106.2&&met_pt>165))'
cuts_loose_zll_met100_dphi='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>100&&llnunu_l1_mass>70&&llnunu_l1_mass<110&&met_pt>100&&cos(llnunu_deltaPhi)<-0.2))'
cuts_loose_zll_met250_dphi='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>100&&llnunu_l1_mass>70&&llnunu_l1_mass<110&&met_pt>250&&cos(llnunu_deltaPhi)<-0.2))'

cuts_zjets='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(abs(llnunu_l1_mass-91.18)<20))'
cuts_zjets_met50='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(abs(llnunu_l1_mass-91.18)<20&&met_pt>50))'

if cutChain=='loosecut': cuts=cuts_loose
elif cutChain=='tightzpt100': cuts=cuts_loose_zll
elif cutChain=='tightzpt100met100': cuts=cuts_loose_zll_met100
elif cutChain=='tightmz15zpt225met165': cuts=cuts_loose_zll_mz15_zpt225_met165
elif cutChain=='tightzpt100met100dphi': cuts=cuts_loose_zll_met100_dphi
elif cutChain=='tightmet250dphi': cuts=cuts_loose_zll_met250_dphi
elif cutChain=='zjetscut': cuts=cuts_zjets
elif cutChain=='zjetscutmet50': cuts=cuts_zjets_met50
else : cuts=cuts_loose


if channel=='el': cuts = cuts+'&&'+elChannel
elif channel=='mu': cuts = cuts+'&&'+muChannel


ROOT.gROOT.ProcessLine('.x tdrstyle.C') 

sigPlotters=[]
sigSamples = [
'BulkGravToZZToZlepZinv_narrow_600', 
#'BulkGravToZZToZlepZinv_narrow_800', 
#'BulkGravToZZToZlepZinv_narrow_1000', 
#'BulkGravToZZToZlepZinv_narrow_1200', 
#'BulkGravToZZToZlepZinv_narrow_1400', 
#'BulkGravToZZToZlepZinv_narrow_1600', 
#'BulkGravToZZToZlepZinv_narrow_1800', 
#'BulkGravToZZToZlepZinv_narrow_2000', 
#'BulkGravToZZToZlepZinv_narrow_2500', 
#'BulkGravToZZToZlepZinv_narrow_3000', 
#'BulkGravToZZToZlepZinv_narrow_3500', 
#'BulkGravToZZToZlepZinv_narrow_4000', 
#'BulkGravToZZToZlepZinv_narrow_4500', 
]


sigXsec = {
'BulkGravToZZToZlepZinv_narrow_600' : 1,
'BulkGravToZZToZlepZinv_narrow_800' : 1,
'BulkGravToZZToZlepZinv_narrow_1000' : 1,
'BulkGravToZZToZlepZinv_narrow_1200' : 1,
'BulkGravToZZToZlepZinv_narrow_1400' : 1,
'BulkGravToZZToZlepZinv_narrow_1600' : 1,
'BulkGravToZZToZlepZinv_narrow_1800' : 1,
'BulkGravToZZToZlepZinv_narrow_2000' : 1,
'BulkGravToZZToZlepZinv_narrow_2500' : 1,
'BulkGravToZZToZlepZinv_narrow_3000' : 1,
'BulkGravToZZToZlepZinv_narrow_3500' : 1,
'BulkGravToZZToZlepZinv_narrow_4000' : 1,
'BulkGravToZZToZlepZinv_narrow_4500' : 1,
}
sigMassWindow = {
'BulkGravToZZToZlepZinv_narrow_600' : [60,560,620],
'BulkGravToZZToZlepZinv_narrow_800' : [60,400,1000],
'BulkGravToZZToZlepZinv_narrow_1000' : [80,400,1200],
'BulkGravToZZToZlepZinv_narrow_1200' : [100,400,1400],
'BulkGravToZZToZlepZinv_narrow_1400' : [120,400,1600],
'BulkGravToZZToZlepZinv_narrow_1600' : [140,400,1800],
'BulkGravToZZToZlepZinv_narrow_1800' : [160,400,2000],
'BulkGravToZZToZlepZinv_narrow_2000' : [180,400,2200],
'BulkGravToZZToZlepZinv_narrow_2500' : [200,400,2800],
'BulkGravToZZToZlepZinv_narrow_3000' : [200,400,3500],
'BulkGravToZZToZlepZinv_narrow_3500' : [200,400,4000],
'BulkGravToZZToZlepZinv_narrow_4000' : [200,400,4500],
'BulkGravToZZToZlepZinv_narrow_4500' : [200,400,5000],
}
for sample in sigSamples:
    sigPlotters.append(TreePlotter(indir+'/'+sample+'/vvTreeProducer/tree.root','tree'))
    sigPlotters[-1].setupFromFile(indir+'/'+sample+'/skimAnalyzerCount/SkimReport.pck')
    sigPlotters[-1].addCorrectionFactor(str(sigXsec[sample]),'tree')
    sigPlotters[-1].addCorrectionFactor('genWeight','tree')
    sigPlotters[-1].addCorrectionFactor('puWeight','tree')
    sigPlotters[-1].setFillProperties(0,ROOT.kWhite)



fout = ROOT.TFile("study_mt_true.root", 'recreate')

c1 = ROOT.TCanvas("c1","c1")
c1.SetBottomMargin(0.15)
c1.SetLeftMargin(0.15)
c1.Print("study_mt_true.ps[")

h1 = []
h2 = []
h3 = []
h4 = []
lg = []

for idx,sample in enumerate(sigSamples):
    nbins = sigMassWindow[sample][0]
    xmin = sigMassWindow[sample][1]
    xmax = sigMassWindow[sample][2]
    h1.append(sigPlotters[idx].drawTH1(sample+'_hmt','genX_mt',cuts,str(lumi*1000),nbins,xmin,xmax,titlex = "True M_{T}",units = "GeV"))
    h2.append(sigPlotters[idx].drawTH1(sample+'_hmta','genX_mta',cuts,str(lumi*1000),nbins,xmin,xmax,titlex = "True M_{T}",units = "GeV"))
    h3.append(sigPlotters[idx].drawTH1(sample+'_hmtb','genX_mtb',cuts,str(lumi*1000),nbins,xmin,xmax,titlex = "True M_{T}",units = "GeV"))
    h4.append(sigPlotters[idx].drawTH1(sample+'_hmtc','genX_mtc',cuts,str(lumi*1000),nbins,xmin,xmax,titlex = "True M_{T}",units = "GeV"))

    h1[idx].SetMarkerSize(0.5)
    h2[idx].SetMarkerSize(0.5)
    h3[idx].SetMarkerSize(0.5)
    h4[idx].SetMarkerSize(0.5)

    h1[idx].SetLineColor(2)
    h1[idx].SetMarkerColor(2)
    h2[idx].SetLineColor(4)
    h2[idx].SetMarkerColor(4)
    h3[idx].SetLineColor(6)
    h3[idx].SetMarkerColor(6)
    h4[idx].SetLineColor(8)
    h4[idx].SetMarkerColor(8)

    lg.append(ROOT.TLegend(0.2,0.5,0.65,0.8))
    #lg[idx].SetName(sample+'_lg')
    #lg[idx].AddEntry(h1[idx],"MT","pl")
    #lg[idx].AddEntry(h2[idx],"MTa","pl")
    ##lg[idx].AddEntry(h3[idx],"MTb","pl")
    #lg[idx].AddEntry(h4[idx],"MTc","pl")

    lg[idx].AddEntry(h1[idx],"M(#nu#nu) = True M(#nu#nu)","pl")
    lg[idx].AddEntry(h2[idx],"M(#nu#nu) = 91.188 GeV","pl")
    lg[idx].AddEntry(h4[idx],"M(#nu#nu) = True M(ll)","pl")

    c1.Clear()
    h1[idx].Draw('l hist')
    h2[idx].Draw('l hist same')
    #h3[idx].Draw('l hist same')
    h4[idx].Draw('l hist same')
    lg[idx].Draw()
    c1.Print("study_mt_true.ps")
    c1.Clear()

    fout.cd()
    h1[idx].Write()
    h2[idx].Write()
    h3[idx].Write()
    h4[idx].Write()
    lg[idx].Write()

c1.Print("study_mt_true.ps]")

os.system('ps2pdf study_mt_true.ps')
fout.Close()



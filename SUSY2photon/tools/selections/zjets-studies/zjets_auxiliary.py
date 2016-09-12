#!/usr/bin/env python

import ROOT
import os
from python.TreePlotter import TreePlotter
from python.MergedPlotter import MergedPlotter
#from python.myStackPlotter import StackPlotter
from python.SimplePlot import *

outtxt = open('CorrFactor.txt', 'a')
if os.stat("CorrFactor.txt").st_size == 0:
    outtxt.write("GetCorrelationFactor(llnunu_l1_deltaR:llnunu_deltaPhi) \n" )
    sout="{0:>6}, {1:>6}, {2:>20}, {3:>20}, {4:>20}\n"
    outtxt.write(sout.format("ZpT", "MET", "inclusive", "mu", "electron"))
else: pass

Channel=raw_input("Please choose a channel (inclusive, el or mu): \n")
pdgID={'el':'11','mu':'13'}
if Channel=="": print "[info] all 3 channels will be plotted."
kincut=raw_input("Please choose ZpT (preselected at 100GeV) > a and met (preselected at 0) > b cuts (no cut applied if you skip): \n a,b = ")
ZpTCut='&&llnunu_l1_pt>'+kincut.split(',')[0]  if kincut.split(',')[0] else ''
MetCut='&&met_pt>'+kincut.split(',')[1]  if kincut.split(',')[1] else ''
#outtxt.write('[NEW] %s %s %s' % (ZpTCut,MetCut,'*'*20))# use for debug

outdir='./plots/aux'
indir="../AnalysisRegion"
#lumi=2.169126704526
if not os.path.exists(outdir): os.system('mkdir '+outdir)

#######----------- Prepare samples to plot:
zjetsPlotters=[]
zjetsSamples = ['DYJetsToLL_M50_BIG'] # M50_BIG = M50 + M50_Ext

for sample in zjetsSamples:
    zjetsPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    #zjetsPlotters[-1].addCorrectionFactor('1./SumWeights','tree') # refers to sum of genWeight
    # zjetsPlotters[-1].addCorrectionFactor('xsec','tree')
    zjetsPlotters[-1].addCorrectionFactor('(genWeight/abs(genWeight))','tree') # to keep the sign but not the value 
    zjetsPlotters[-1].addCorrectionFactor('puWeight','tree')
    zjetsPlotters[-1].addCorrectionFactor('triggersf','tree')
    zjetsPlotters[-1].addCorrectionFactor('llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf','tree')
        
ZJets = MergedPlotter(zjetsPlotters)
ZJets.setFillProperties(1001,ROOT.kGreen+2)
    
#######----------- Start Plotting:

lsChannel=[]
if Channel=="": lsChannel=['inclusive','el','mu']
else: lsChannel.append(Channel)

res=dict()
for Channel in lsChannel:
    print "[info] ", Channel, '*'*20
    ROOT.gROOT.ProcessLine('.x tdrstyle.C')

    if Channel=='inclusive':
        factor_cuts='(abs(llnunu_l1_mass-91.1876)<20.0'+ZpTCut+MetCut+')'
    else:
        factor_cuts='(abs(llnunu_l1_mass-91.1876)<20.0&&abs(llnunu_l1_l1_pdgId)=='+pdgID[Channel]+ZpTCut+MetCut+')'
        
        
    print '[info] cuts used here: ', factor_cuts
    
    h2_dR_dphi=ZJets.drawTH2("llnunu_l1_deltaR:llnunu_deltaPhi",factor_cuts,str(1),#lumi*1000
                             32,0,3.2, 32,0,3.2,
                             titlex = "#Delta#Phi_{Z,MET}",unitsx = "",
                             titley = "#Delta R_{l,l}", unitsy = "", drawStyle = "COLZ")

    ROOT.gStyle.SetPadBottomMargin(0.15);
    ROOT.gStyle.SetPadLeftMargin(0.15);
    ROOT.gStyle.SetPadRightMargin(0.12);
    ROOT.gStyle.SetTitleXOffset(0.5);
    ROOT.gStyle.SetTitleYOffset(0.5);
    c1=ROOT.TCanvas("c1","c1",1)

    h2_dR_dphi.Draw("COLZ")
    print h2_dR_dphi.GetEntries(), h2_dR_dphi.GetSumOfWeights()
    
    print "llnunu_l1_deltaR:llnunu_deltaPhi, correlation factor r = %.5f +- %.5f" %( h2_dR_dphi.GetCorrelationFactor(), GetCorrelationFactorError(h2_dR_dphi.GetCorrelationFactor(), h2_dR_dphi.GetSumOfWeights()))

    res[Channel]=("%.5f +- %.5f" % ( h2_dR_dphi.GetCorrelationFactor(), GetCorrelationFactorError(h2_dR_dphi.GetCorrelationFactor(), h2_dR_dphi.GetSumOfWeights())))
    
    c1.SaveAs(outdir+"/h2_dRll_dPhiZMet_"+Channel+'_'+kincut.split(',')[0]+'_'+kincut.split(',')[1]+".eps")

print res
sout2="{inclusive:>20}, {mu:>20}, {el:>20}\n"
outtxt.write("{0:>6}, {1:>6}, ".format(*kincut.split(','))+sout2.format(**res))


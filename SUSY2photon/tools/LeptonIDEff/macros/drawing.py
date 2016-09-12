#! /usr/bin/env python
from ROOT import *
import sys

if len(sys.argv)<2:
    print "help: comparing.py foldertag [title]"
    sys.exit()

gROOT.SetBatch()
gStyle.SetOptStat(000000)
gStyle.SetLegendBorderSize(0)
feff=TFile('/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/tools/LeptonIDEff/rootfiles/{0}/effratio_npnm_mu_set.root'.format(sys.argv[1]))
prstr="/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/tools/LeptonIDEff/plots/{0}/effsf.pdf".format(sys.argv[1])
c=TCanvas('c','c',1600,1600)
names=[' mc',' data',' mc gen']
effps=[feff.Get("hmc_eff"),feff.Get("hdt_eff"),feff.Get("hmcgen_eff"),feff.Get("hratio_dt_mc")]
#leg=TLegend(.75,.91,.9,1)    
leg=TLegend(.6,.7,.9,.89)    
for i in effps:
    for n in range(i.GetNbinsX()):
        i.SetBinError(n+1,1.12*i.GetBinError(n+1))
for i in range(3):
    effps[i].SetLineColor(i+1)
    effps[i].SetMarkerColor(i+1)
    effps[i].SetMarkerStyle(i+20)
    leg.AddEntry(effps[i],names[i])
    if not i:
        effps[i].SetTitle(sys.argv[1] if len(sys.argv)==2 else sys.argv[2])
        mi=min([effps[x].GetMinimum() for x in range(3)])
        ma=max([effps[x].GetMaximum() for x in range(3)])
        effps[i].SetMinimum(2*mi-ma)
        effps[i].SetMaximum(3*ma-2*mi)
        effps[i].Draw()
    else:
        effps[i].Draw('same')
leg.Draw('same')
c.Print(prstr+'(')
effps[3].SetTitle(sys.argv[1] if len(sys.argv)==2 else sys.argv[2])
mi=effps[3].GetMinimum()
ma=effps[3].GetMaximum()
effps[3].SetMinimum(2*mi-ma)
effps[3].SetMaximum(2*ma-mi)
effps[3].Draw()
c.Print(prstr+')')


#! /usr/bin/env python
from ROOT import *
import sys

if len(sys.argv)<2:
    print "help: comparing.py foldertag [title] [minimum z]"
    sys.exit()
title=sys.argv[2] if len(sys.argv)>2 else sys.argv[1]
miniz=float(sys.argv[3]) if len(sys.argv)>3 else .96
gROOT.SetBatch()
gStyle.SetOptStat(000000)
gStyle.SetLegendBorderSize(0)
feff=TFile('/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/tools/LeptonIDEff/rootfiles/{0}/effratio_npnm_mu_set.root'.format(sys.argv[1]))
c=TCanvas('c','c',1600,1600)
c.SetRightMargin(0.2)
names=[' mc',' mc gen',' data',' SF']
effps=[feff.Get("hmc_eff"),feff.Get("hmcgen_eff"),feff.Get("hdt_eff"),feff.Get("hratio_dt_mc")]
for i in range(len(effps)):
    prstr="/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/tools/LeptonIDEff/plots/{0}/effsf.pdf".format(sys.argv[1])
    effps[i].SetTitle('{0}{1};muon_Pt(GeV);muon_eta'.format(title,names[i]))
    effps[i].SetMinimum(miniz)
    effps[i].SetMaximum(1)
    effps[i].GetZaxis().SetTitleOffset(2)
    effps[i].Draw("colztext45")
    if not i:prstr+="("
    if i==len(effps)-1:prstr+=")"
    c.Print(prstr)

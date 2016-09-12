#! /usr/bin/env python
from ROOT import *
import sys

if len(sys.argv)<3:
    print "help: comparing.py foldertag roottag [minimum z]"
    print "roottag: SoftID/LooseID/MediumID/TightIDandIPCut"
    sys.exit()

miniz=float(sys.argv[3]) if len(sys.argv)>3 else .96
gROOT.SetBatch()
gStyle.SetOptStat(000000)
gStyle.SetLegendBorderSize(0)
fpog=TFile("/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/tools/LeptonIDEff/rootfiles/MuonID_Z_RunCD_Reco76X_Feb15.root")
feff=TFile('/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/tools/LeptonIDEff/rootfiles/{0}/effratio_npnm_mu_set.root'.format(sys.argv[1]))
c=TCanvas('c','c',1600,1600)
c.SetRightMargin(0.2)
pogdir='MC_NUM_{0}_DEN_genTracks_PAR_pt_spliteta_bin1'.format(sys.argv[2])
names=[' mc POG',' mc Wmass',' mc gen',' data POG',' data Wmass',' SF POG',' SF Wmass']
effps=[fpog.Get(pogdir+'/efficienciesMC/pt_abseta_MC'),feff.Get("hmc_eff"),feff.Get("hmcgen_eff"),fpog.Get(pogdir+'/efficienciesDATA/pt_abseta_DATA'),feff.Get("hdt_eff"),fpog.Get(pogdir+'/pt_abseta_ratio'),feff.Get("hratio_dt_mc")]
for i in range(len(effps)):
    prstr="/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/tools/LeptonIDEff/plots/{0}/sfcomp.pdf".format(sys.argv[1])
    effps[i].SetTitle(sys.argv[2].split('and')[0]+names[i]+';muon_Pt(GeV);muon_eta;Efficiency of {0}::above'.format(sys.argv[2].split('and')[0]))
    effps[i].SetMinimum(miniz)
    effps[i].SetMaximum(1)
    effps[i].GetZaxis().SetTitleOffset(2)
    effps[i].Draw("colztext45")
    if not i:prstr+="("
    if i==len(effps)-1:prstr+=")"
    c.Print(prstr)

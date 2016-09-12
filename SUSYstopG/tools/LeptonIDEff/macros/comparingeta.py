#! /usr/bin/env python
from ROOT import *
import sys

if len(sys.argv)!=3:
    print "help: comparing.py foldertag roottag"
    print "roottag: SoftID/LooseID/MediumID/TightIDandIDCut"
    sys.exit()

gROOT.SetBatch()
gStyle.SetOptStat(000000)
gStyle.SetLegendBorderSize(0)
fpog=TFile("/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/tools/LeptonIDEff/rootfiles/MuonID_Z_RunCD_Reco76X_Feb15.root")
feff=TFile('/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/tools/LeptonIDEff/rootfiles/{0}/effratio_npnm_mu_set.root'.format(sys.argv[1]))
c=TCanvas('c','c',1600,1600)
pogdir='MC_NUM_{0}_DEN_genTracks_PAR_eta'.format(sys.argv[2])
names=[' mc',' data',' SF']
pogps=[fpog.Get(pogdir+'/efficienciesMC/histo_eta_MC'),fpog.Get(pogdir+'/efficienciesDATA/histo_eta_DATA'),fpog.Get(pogdir+'/eta_ratio')]
effps=[feff.Get("hmc_eff"),feff.Get("hdt_eff"),feff.Get("hratio_dt_mc"),feff.Get("hmcgen_eff")]
for i in range(len(pogps)):
    prstr="/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/tools/LeptonIDEff/plots/{0}/sfcomp.pdf".format(sys.argv[1])
    pogps[i].SetLineColor(2)
    pogps[i].SetMarkerColor(2)
    pogps[i].SetMarkerStyle(20)
    pogps[i].SetTitle(sys.argv[2].split('and')[0]+names[i]+';muon_eta;')
    effps[i].SetLineColor(4)
    effps[i].SetMarkerColor(4)
    effps[i].SetMarkerStyle(22)
    for n in range(effps[i].GetNbinsX()):
        effps[i].SetBinError(n+1,1.12*effps[i].GetBinError(n+1))
    mi=min(pogps[i].GetMinimum(),effps[i].GetMinimum())
    ma=max(pogps[i].GetMaximum(),effps[i].GetMaximum())
    pogps[i].SetMinimum(2*mi-ma)
    pogps[i].SetMaximum(3*ma-2*mi)
#    leg=TLegend(.75,.91,.9,1)    
    leg=TLegend(.6,.7,.9,.89)    
    leg.AddEntry(pogps[i],'POG'+names[i])
    leg.AddEntry(effps[i],'Wmass'+names[i])
    pogps[i].Draw()
    effps[i].Draw('same')
    if not i:
        prstr+="("
        effps[3].SetLineColor(3)
        effps[3].SetMarkerColor(3)
        effps[3].SetMarkerStyle(23)
        effps[3].Draw("same")
        leg.AddEntry(effps[3],'Gen Efficiency')
    leg.Draw("same")
    if i==len(pogps)-1:prstr+=")"
    c.Print(prstr)

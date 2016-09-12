#!/usr/bin/env python

from ROOT import *
from array import array

#pu_list=[500,530,570,600,630,670,675,680,685,691,698]
#pu_list=[620, 625, 630, 635, 640, 645, 650, 655, 660, 665, 670, 675, 680, 685, 691, 698]
#pu_list=[660, 665, 670, 675, 680, 685, 690, 695, 700, 705, 710, 713, 715, 720, 725, 730, 735, 740, 745, 750, 755, 760, 765]
#pu_list=[575, 580, 585, 590, 595, 600, 605, 610, 615, 620, 625, 630, 635, 640, 645, 650, 655, 660, 665, 670, 675, 680, 685, 690, 695, 700, 705, 710, 713, 715, 720, 725, 730, 735, 740, 745, 750, 755, 760, 765]
pu_list=[691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705]


inputdir='/data/XZZ/80X_Ntuple/80X_20160618_Skim'


outtag='scan_pu'

DataHLT=True

if not DataHLT:
  outtag += '_nohlt'

dtfileName='SingleEMU_Run2016B_PromptReco_v1v2.root'
mcfileName='DYJetsToLL_M50_PUScan.root'

lumiTag = 'CMS 13 TeV 2016B L=2.6 fb^{-1}'

cuts_loose_z='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_mass>70&&llnunu_l1_mass<110))'

cuts=cuts_loose_z


gROOT.ProcessLine('.x tdrstyle.C')

dtfile = TFile(inputdir+'/'+dtfileName)
mcfile = TFile(inputdir+'/'+mcfileName)
dttree=dtfile.Get('tree')
mctree=mcfile.Get('tree')



fout = TFile(outtag+'.root','recreate')

plots = TCanvas('plots','plots',600,800)
plots.SetLeftMargin(0.15)
plots.SetBottomMargin(0.15)


lumipt = TPaveText(0.45,0.8,0.9,0.9,"brNDC")
lumipt.SetBorderSize(0)
lumipt.SetTextAlign(12)
lumipt.SetFillStyle(0)
lumipt.SetTextFont(42)
lumipt.SetTextSize(0.045)
lumipt.AddText(0.15,0.15, lumiTag)

plots.Print(outtag+'.ps[')


ress = array('d', 50*[0.0])

hs_dt = []
hs_mc = []
hs_res = []
lg_dtmc = []
chi2ndf_pu = []
pvalue_pu = []


for pu in pu_list:
    hdt = TH1D("hdt","hdt_"+str(pu),50,0,50)    
    hmc = TH1D("hmc","hmc_"+str(pu),50,0,50)    
    hdt.SetName('hdt_'+str(pu))
    hmc.SetName('hmc_'+str(pu))
    hdt.SetTitle('nVert: MB xsec= {pu:n}00 ub {lumitag:s}'.format(pu=pu,lumitag=lumiTag))
    hdt.Sumw2()
    hmc.Sumw2()
    if not DataHLT: dttree.Draw('nVert>>hdt_'+str(pu), '('+cuts+')', 'e')
    else: dttree.Draw('nVert>>hdt_'+str(pu), '('+cuts+'&&((HLT_MU&&!HLT_ELE)||(HLT_ELE)))', 'e')
    #mctree.Draw('nVert>>hmc_'+str(pu), '(('+cuts+')*(genWeight*ZPtWeight*puWeight'+str(pu)+'))', 'e')
    mctree.Draw('nVert>>hmc_'+str(pu), '(('+cuts+')*(genWeight*puWeight'+str(pu)+'))', 'e')
    hmc.Scale(hdt.Integral()/hmc.Integral())
    hdt.GetXaxis().SetTitle('N vertices')
    hmc.GetXaxis().SetTitle('N vertices')
    hdt.GetYaxis().SetTitle('Events')
    hmc.GetYaxis().SetTitle('Events')
    hdt.SetMarkerStyle(20)
    hmc.SetLineColor(4)
    hmc.SetMarkerColor(4)
    ymax = max(hdt.GetMaximum(),hmc.GetMaximum())*1.2
    hdt.GetYaxis().SetRangeUser(0,ymax)
    hmc.GetYaxis().SetRangeUser(0,ymax)

    pv = hdt.Chi2Test(hmc, "UW P", ress)
    #print 'p-value pu '+str(pu)+' : '+str(pv)
    pvalue_pu.append(pv)
    chi2ndf = hdt.Chi2Test(hmc, "UW P CHI2/NDF", ress)
    chi2ndf_pu.append(chi2ndf)
    hres = TH1D("hres_"+str(pu),"hres_"+str(pu), 50,0, 50)
    hres.Sumw2()
    for ib in range(hres.GetNbinsX()):
        hres.SetBinContent(ib+1, ress[ib])
        hres.SetBinError(ib+1, 1.0)
    hres.SetMarkerStyle(20)  
    hres.SetMarkerColor(2)
    hres.GetXaxis().SetTitle('N vertices')
    hres.GetYaxis().SetTitle('Normalized Residual (Pull)')
    hres.SetTitle('Pull: MB xsec= {pu:n}00 ub , Chi2/Ndf= {chi2ndf:.2f}, P-value= {pv:.3E}'.format(pu=pu,chi2ndf=chi2ndf,pv=pv))

    plots.Clear()
    plots.Divide(1,2)
    pad1=plots.cd(1)
    pad1.SetTopMargin(0.10)
    pad1.SetBottomMargin(0.15)
    pad1.SetLeftMargin(0.15)
    hdt.Draw()
    hmc.Draw('hist same')
    lumipt.Draw()
    pad2=plots.cd(2)
    pad2.SetTopMargin(0.10)
    pad2.SetBottomMargin(0.15)
    pad2.SetLeftMargin(0.15)
    hres.Draw()
    lumipt.Draw()
    plots.Print(outtag+'.ps')
    plots.Clear()

    hs_dt.append(hdt)
    hs_mc.append(hmc)
    hs_res.append(hres)
  
    fout.cd()
    hdt.Write()
    hmc.Write()
    hres.Write() 
    plots.Write('plots_'+str(pu))


gr_pvalue_pu = TGraph()
for i in range(len(pu_list)):
  gr_pvalue_pu.SetPoint(i,pu_list[i]*0.1,pvalue_pu[i])
gr_pvalue_pu.SetName('gr_pvalue_pu')
gr_pvalue_pu.SetTitle('Pileup Data/MC matching p-value vs. MB xsec')
gr_pvalue_pu.SetMarkerStyle(20)
gr_pvalue_pu.GetXaxis().SetTitle('MB xsec (mb)')
gr_pvalue_pu.GetYaxis().SetTitle('p-value')
gr_pvalue_pu.GetXaxis().SetRangeUser(min(pu_list)*0.1*0.9,max(pu_list)*0.1*1.1)
gr_pvalue_pu.GetYaxis().SetRangeUser(min(pvalue_pu)*1e-10,max(pvalue_pu)*1e20)
gr_pvalue_pu.Sort()

#fc_pvalue_pu = TF1("fc_pvalue_pu", "pol2", min(pu_list)*0.1*0.9, max(pu_list)*0.1*1.1)
#gr_pvalue_pu.Fit(fc_pvalue_pu)
#x_maxpvalue=fc_pvalue_pu.GetMaximumX()
#gr_pvalue_pu.SetTitle('Pileup Data/MC matching p-value vs. MB xsec, max at {xm:.3f} mb'.format(xm=x_maxpvalue))


gr_chi2ndf_pu = TGraph()
for i in range(len(pu_list)):
  gr_chi2ndf_pu.SetPoint(i,pu_list[i]*0.1,chi2ndf_pu[i])
gr_chi2ndf_pu.SetName('gr_chi2ndf_pu')
gr_chi2ndf_pu.SetTitle('Pileup Data/MC matching Chi2/Ndf vs. MB xsec')
gr_chi2ndf_pu.SetMarkerStyle(20)
gr_chi2ndf_pu.GetXaxis().SetTitle('MB xsec (mb)')
gr_chi2ndf_pu.GetYaxis().SetTitle('Chi2/Ndf')
gr_chi2ndf_pu.GetXaxis().SetRangeUser(min(pu_list)*0.1*0.9,max(pu_list)*0.1*1.1)
gr_chi2ndf_pu.GetYaxis().SetRangeUser(min(chi2ndf_pu)*0.98,max(chi2ndf_pu)*1.1)
gr_chi2ndf_pu.Sort()

fc_chi2ndf_pu = TF1("fc_chi2ndf_pu", "pol2", min(pu_list)*0.1*0.9, max(pu_list)*0.1*1.1)
gr_chi2ndf_pu.Fit(fc_chi2ndf_pu)
x_minchi2ndf=fc_chi2ndf_pu.GetMinimumX()
gr_chi2ndf_pu.SetTitle('Pileup Data/MC matching Chi2/Ndf vs. MB xsec, min at {xm:.3f} mb'.format(xm=x_minchi2ndf))

print '#########################################################################'
print ' Pileup Mini-Bias cross-section scan results: best MB xsec = {xm:.3f} mb'.format(xm=x_minchi2ndf)
print '#########################################################################'

plots.Clear()
plots.Divide(1,2)
pad1=plots.cd(1)
pad1.SetTopMargin(0.10)
pad1.SetBottomMargin(0.15)
pad1.SetLeftMargin(0.15)
pad1.SetLogy(1)
gr_pvalue_pu.Draw("apl")
lumipt.Draw()
pad2=plots.cd(2)
pad2.SetTopMargin(0.10)
pad2.SetBottomMargin(0.15)
pad2.SetLeftMargin(0.15)
gr_chi2ndf_pu.Draw("apl")
lumipt.Draw()
plots.Print(outtag+'.ps')
pad1.SetLogy(0)
plots.Clear()

fout.cd()
gr_pvalue_pu.Write()
gr_chi2ndf_pu.Write()
plots.Write('plots_gr_pvalue_chi2ndf_pu')


plots.Print(outtag+'.ps]')

gROOT.ProcessLine('.! ps2pdf '+outtag+'.ps')

fout.Close()

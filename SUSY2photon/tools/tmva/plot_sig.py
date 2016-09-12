#!/usr/bin/env python


from ROOT import *
import math

#tags = ['m1800']
tags = ['m600', 'm800', 'm1000', 'm1200', 'm1400', 'm1600', 'm1800', 'm2000', 'm2500', 'm3000', 'm3500', 'm4000']
xranges = {'m600':[0,2000], 'm800':[0,2000], 'm1000':[0,2000], 'm1200':[0,2000], 'm1400':[0,2000], 'm1600':[0,2000], 'm1800':[0,2500], 'm2000':[0,2500], 'm2500':[0,3000], 'm3000':[0,3500], 'm3500':[0,3700], 'm4000':[0,4200]} 

gROOT.ProcessLine('.x tdrstyle.C')

def hist_signif(name,Nsig,Nbkg,Hsig,Hbkg,NbinsX):

    Xmin = Hsig.GetXaxis().GetXmin()
    Xmax = Hsig.GetXaxis().GetXmax()
    Hsignif = TH1F(name, name, NbinsX, Xmin, Xmax)
    Hsignif.SetTitle("Significance vs. MVA cuts")

    for i in range(1, NbinsX+1):
        i_nsig = 0
        i_nbkg = 0
        i_bin_min = Hsignif.GetXaxis().GetBinLowEdge(i)
        i_bin_max = Hsignif.GetXaxis().GetBinUpEdge(i)
        for j in range(1,Hsig.GetNbinsX()+1):
            j_bin_center = Hsig.GetXaxis().GetBinCenter(j)
            if i_bin_min<j_bin_center<i_bin_max:
                i_nsig += Nsig*Hsig.GetBinContent(j)
                i_nbkg += Nbkg*Hbkg.GetBinContent(j)
        #print "i_nsig = "+str(i_nsig)+"; i_nbkg = "+str(i_nbkg)
        
        sb = i_nsig+(i_nbkg if i_nbkg>0.0 else 0.0)
        i_signif = i_nsig/math.sqrt(sb) 
        Hsignif.SetBinContent(i, i_signif)

    return Hsignif

def plot_signif(tag):
    infile_name = 'TMVA_'+tag+'.root'
    fin = TFile(infile_name)

    c1 = TCanvas("c1", "c1")
    testTree = fin.Get("TestTree")
    htest=TH1F("htest", "htest", 100,-100000,100000)
    testTree.Draw("BDT>>htest", "(classID==0)*(weight)")
    nsig=htest.Integral()*2.0
    testTree.Draw("BDT>>htest", "(classID==1)*(weight)")
    nbkg=htest.Integral()*2.0

    hsig = fin.Get("Method_BDT/BDT/MVA_BDT_effS")
    hbkg = fin.Get("Method_BDT/BDT/MVA_BDT_effB")
    hsigif = hist_signif(name="hsignif", Nsig=nsig, Nbkg=nbkg, Hsig=hsig, Hbkg=hbkg, NbinsX=100)

    BDT_cut="BDT>"+str(hsigif.GetBinCenter(hsigif.GetMaximumBin()))
    print "#### "+tag+" : "+BDT_cut

    c1.Clear()
    hsigif.GetXaxis().SetTitle("MVA cuts")
    hsigif.GetYaxis().SetTitle("significance")
    hsigif.Draw()
    c1.SaveAs("plot_signif_"+tag+".pdf")

    c1.Clear()
    c1.SetLogy()
    hsig.SetLineColor(2)
    hbkg.SetLineColor(4)
    hsig.GetYaxis().SetRangeUser(0.001,1.5)
    hbkg.GetYaxis().SetRangeUser(0.001,1.5)
    hsig.Draw()
    hbkg.Draw("same")
    c1.SaveAs("plot_eff_"+tag+".pdf") 

    fin.Close()

    return BDT_cut


def drawTH1(tag, name, var, xtitle, mvaCut, sigScale, nbins, xmin, xmax, logView):
    infile_name = 'TMVA_'+tag+'.root'
    fin = TFile(infile_name)
    testTree = fin.Get("TestTree")
    hsig = TH1F(name+"_sig", name, nbins, xmin, xmax) 
    hbkg = TH1F(name+"_bkg", name, nbins, xmin, xmax) 
    hsig.Sumw2()
    hbkg.Sumw2()
    testTree.Draw(var+'>>'+hsig.GetName(),'('+mvaCut+'&&classID==0)*(weight*'+str(sigScale)+')')
    testTree.Draw(var+'>>'+hbkg.GetName(),'('+mvaCut+'&&classID==1)*(weight)')

    hsig.SetLineColor(2)
    hsig.SetMarkerColor(2)
    hbkg.SetLineColor(4)
    hbkg.SetMarkerColor(4)
    
    hsig.GetXaxis().SetTitle(xtitle)
    hbkg.GetXaxis().SetTitle(xtitle)
   
    hsig.GetYaxis().SetTitle('Nevts')
    hbkg.GetYaxis().SetTitle('Nevts')

    lg = TLegend(0.6,0.7,0.89,0.89)
    lg.SetName(name+"_lg")
    lg.AddEntry(hsig, "Signal", "l")
    lg.AddEntry(hbkg, "Background", "l")
  
    ymax = max(hsig.GetMaximum(),hbkg.GetMaximum())*1.1
    ymin = min(hsig.GetMaximum(),hbkg.GetMaximum())*0.001

    hsig.GetYaxis().SetRangeUser(0.0001, ymax)
    hbkg.GetYaxis().SetRangeUser(0.0001, ymax)

    c1 = TCanvas(name+"_c1", "c1")
    if logView:
        c1.SetLogy()
        hsig.GetYaxis().SetRangeUser(ymin, ymax*10)
        hbkg.GetYaxis().SetRangeUser(ymin, ymax*10)
        
    hsig.Draw("hist")
    hbkg.Draw("hist same")
    lg.Draw("same")

    c1.SaveAs("plot_signif_"+name+"_"+tag+".pdf")

    fin.Close()


if __name__ == '__main__':
    for tag in tags:    
        BDT_cut = plot_signif(tag)
        drawTH1(tag, "mt_bdtcut", "llnunu_mt", "M_{T} (GeV)", BDT_cut, 10, 60, xranges[tag][0], xranges[tag][1], True)
        drawTH1(tag, "mt_nocut", "llnunu_mt", "M_{T} (GeV)", "BDT>-100000", 10, 60, xranges[tag][0], xranges[tag][1], True)


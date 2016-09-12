#!/usr/bin/env python

from ROOT import *
gROOT.ProcessLine('.x tdrstyle.C')

lumi=76.6
#pileup pdf from:  https://github.com/cms-sw/cmssw/blob/CMSSW_7_6_3/SimGeneral/MixingModule/python/mix_2015_25ns_FallMC_matchData_PoissonOOTPU_cfi.py
pdf = [0.000108643,
                0.000388957,
                0.000332882,
                0.00038397,
                0.000549167,
                0.00105412,
                0.00459007,
                0.0210314,
                0.0573688,
                0.103986,
                0.142369,
                0.157729,
                0.147685,
                0.121027,
                0.08855,
                0.0582866,
                0.0348526,
                0.019457,
                0.0107907,
                0.00654313,
                0.00463195,
                0.00370927,
                0.0031137,
                0.00261141,
                0.00215499,
                0.00174491,
                0.00138268,
                0.00106731,
                0.000798828,
                0.00057785,
                0.00040336,
                0.00027161,
                0.000176535,
                0.00011092,
                6.75502e-05,
                4.00323e-05,
                2.32123e-05,
                1.32585e-05,
                7.51611e-06,
                4.25902e-06,
                2.42513e-06,
                1.39077e-06,
                8.02452e-07,
                4.64159e-07,
                2.67845e-07,
                1.5344e-07,
                8.68966e-08,
                4.84931e-08,
                2.6606e-08,
                1.433e-08, 
                0.0, 0.0, 0.0]

fdt=TFile('pileup_DATA_76x_50ns.root')
pileup_dt = fdt.Get('pileup')
pileup_dt.Scale(1.0/pileup_dt.Integral())
# still save 25ns mc pileup, but data/mc ratio take the 50ns data
fmc=TFile('pileup_MC_76x_50ns.root','recreate')
pileup_mc=TH1F('pileup','pileup',52,0,52)
pileup_mc.Sumw2()
for i in range(52):    pileup_mc.SetBinContent(i+1,pdf[i])
pileup_mc.Scale(1.0/pileup_mc.Integral())
gStyle.SetStatStyle(0)
pileup_dt.SetStats(0)
pileup_mc.SetStats(0)

c1 = TCanvas("c1","c1",600,600)
c1.SetLeftMargin(0.15)
c1.SetBottomMargin(0.15)

pileup_dt.SetLineColor(kRed)
pileup_dt.SetMarkerColor(kRed)
pileup_dt.SetMarkerStyle(20)
pileup_dt.GetXaxis().SetTitle('N pileup')
pileup_dt.GetYaxis().SetTitle('Norm.')
pileup_mc.SetLineColor(kBlue)
pileup_mc.SetMarkerColor(kBlue)
pileup_mc.SetMarkerStyle(20)
pileup_mc.SetLineWidth(2)
pileup_mc.GetXaxis().SetTitle('N pileup')
pileup_mc.GetYaxis().SetTitle('Norm.')

pileup_mc.GetYaxis().SetTitleOffset(1.5)
pileup_dt.GetYaxis().SetTitleOffset(1.5)
lg = TLegend(0.6,0.7,0.85,0.85)
lg.AddEntry(pileup_dt, "Data 50ns", "pl")
lg.AddEntry(pileup_mc, "MC 25ns", "pl")

pt = TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC")
pt.SetBorderSize(0)
pt.SetTextAlign(12)
pt.SetFillStyle(0)
pt.SetTextFont(42)
pt.SetTextSize(0.03)
text = pt.AddText(0.15,0.3,"CMS Preliminary")
text = pt.AddText(0.55,0.3,"#sqrt{s} = 13 TeV, L = "+"{:.3}".format(float(lumi))+" pb^{-1}")

pileup_dt.SetMaximum(pileup_mc.GetBinContent(pileup_mc.GetMaximumBin())*1.2)

pileup_dt.Draw()
pileup_mc.Draw("HIST,SAME")
lg.Draw("same")
pt.Draw()
c1.SaveAs("pileup76x50ns.eps")
c1.SaveAs("pileup76x50ns.pdf")


# ratio
puweight = pileup_dt.Clone("puweight")
puweight.Divide(pileup_mc)
puweight.GetYaxis().SetTitle("PU Weight")
lg1 = TLegend(0.56,0.7,0.85,0.85)
lg1.AddEntry(puweight, "Data/MC ", "pl")
lg1.AddEntry(puweight, "2015 B+C 50ns", "")

c2 = TCanvas("c2","c2",600,600)
c2.SetLeftMargin(0.15)
c2.SetBottomMargin(0.15)
puweight.Draw()
lg1.Draw()
pt.Draw()
c2.SaveAs("puweight76x50ns.eps")
c2.SaveAs("puweight76x50ns.pdf")

# print ratio 
puwts = [puweight.GetBinContent(i+1) for i in range(52)]
print puwts

pileup_mc.Write()
puweight.Write()
fmc.Close()

#!/usr/bin/env python
from ROOT import *
import math, os, sys, copy

def GetCorrelationFactorError(corr, num):
    se=math.sqrt((1-corr**2)/(num-2))
    return se
    
def GetError(A, B, a=0., b=0.):
    # q=A/B, sigma(A)=a, sigma(B)=b, sigma(q)=error:
    if a*b==0:
        a=math.sqrt(A)
        b=math.sqrt(B)
    error=math.sqrt((a/B)**2+(b*A/B**2)**2)
    return error

def GetRatio_TH1(h1, h2, h2_isStack=False):
    ''' h1/h2 '''
    hratio = h1.Clone("hRatio")
    if h2_isStack:
        h2_new = h2.GetHistogram()
        h2_new.SetName('hstackmerge')
        for hist in h2.GetHists():
            h2_new.Add(hist)
    else:
        h2_new=h2.Clone("single h2 as denominator")
            
    for i in xrange(h1.GetXaxis().GetNbins()):
        N1 = h1.GetBinContent(i)
        N2 = h2_new.GetBinContent(i)
        E1 = h1.GetBinError(i)
        E2 = h2_new.GetBinError(i)
        RR = N1/N2 if N2>0 else 0
        EE = 0 if N2<=0 else GetError(N1,N2,E1,E2) #math.sqrt(N2*N2*E1*E1+N1*N1*E2*E2)/N2/N2

        hratio.SetBinContent(i, RR)
        hratio.SetBinError(i, EE)
        # if blinding and h1.GetBinCenter(i)>blindingCut: 
        #     hratio.SetBinContent(i, 0)
        #     hratio.SetBinError(i, 0)

    hratio.SetMarkerStyle(20)
    hratio.SetLineWidth(1)
    hratio.SetMarkerSize(1.)
    hratio.SetMarkerColor(kBlack)
    hratio.GetXaxis().SetTitle('')
    hratio.GetYaxis().SetTitle('Data/MC')
    hratio.GetYaxis().SetRangeUser(0.0,2.0)
    hratio.GetXaxis().SetLabelFont(42)
    hratio.GetXaxis().SetLabelOffset(0.007)
    hratio.GetXaxis().SetLabelSize(0.1)
    hratio.GetXaxis().SetTitleSize(0.05)
    hratio.GetXaxis().SetTitleOffset(1.15)
    hratio.GetXaxis().SetTitleFont(42)
    hratio.GetYaxis().SetLabelFont(42)
    hratio.GetYaxis().SetLabelOffset(0.01)
    hratio.GetYaxis().SetLabelSize(0.06)
    hratio.GetYaxis().SetTitleSize(0.12)
    hratio.GetYaxis().SetTitleOffset(0.5)
    hratio.GetYaxis().SetTitleFont(42)
    hratio.GetZaxis().SetLabelFont(42)
    hratio.GetZaxis().SetLabelOffset(0.007)
    hratio.GetZaxis().SetLabelSize(0.045)
    hratio.GetZaxis().SetTitleSize(0.05)
    hratio.GetZaxis().SetTitleFont(42)

    return hratio

def GetLegend(h1,label1,opt1,h2,label2,opt2):
    legend=TLegend(0.45,0.75,0.8,0.90,"","brNDC");
    legend.SetFillStyle(0); #set legend box transparent
    legend.SetBorderSize(0);
    legend.SetTextSize(0.05);
    legend.SetTextFont(42);
#    legend.SetNColumns(2);
    legend.SetLineColor(0);
    legend.AddEntry(h1,label1,opt1);
    legend.AddEntry(h2,label2,opt2);
    return legend

def drawStack_simple(frame, hstack, hdata, hratio, legend,
                     hstack_opt="nostack",
                     outDir="./", output="output", channel="",
                     xmin=50., xmax=500., xtitle="" ,units="",
                     lumi=2.169, notes=""):
    
    fout = TFile(outDir+'/'+output+'_'+channel.Data()+'.root', 'recreate')
        
    c1 = TCanvas(output+'_'+"c1", "c1", 600, 750); c1.Draw()
    c1.SetWindowSize(600 + (600 - c1.GetWw()), (750 + (750 - c1.GetWh())))
    p1 = TPad(output+'_'+"pad1","pad1",0,0.25,1,0.99)
    p1.SetBottomMargin(0.15)
    p1.SetLeftMargin(0.15)
    p1.Draw()
    p2 = TPad(output+'_'+"pad2","pad2",0,0,1,0.25)
    p2.SetTopMargin(0.03)
    p2.SetBottomMargin(0.3)
    p2.SetLeftMargin(0.15)
    p2.SetFillStyle(0)
    p2.Draw()
    
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)
    
    p1.cd()
    frame.Draw()
    hstack.Draw(hstack_opt+", same")
    hdata.Draw("Psame")
    legend.Draw("same")
    hstack.SetMinimum(0.01)
    hstack.GetHistogram().GetXaxis().SetRangeUser(xmin,xmax)

    maxi=hstack.GetHistogram().GetXaxis().GetXmax()
    mini=hstack.GetHistogram().GetXaxis().GetXmin()
    bins=hstack.GetHistogram().GetNbinsX()
    if len(units)>0:
        hstack.GetHistogram().GetXaxis().SetTitle(xtitle + " (" +units+")")
        hstack.GetHistogram().GetYaxis().SetTitle("Events / "+str((maxi-mini)/bins)+ " "+units)
    else:
        hstack.GetHistogram().GetXaxis().SetTitle(xtitle)
        hstack.GetHistogram().GetYaxis().SetTitle("Events")

    hratio.SetName(output+'_'+'hratio')
    hline = copy.deepcopy(hratio)
    for ii in range(hline.GetNbinsX()+1):
        hline.SetBinContent(ii,1.0)
        hline.SetBinError(ii,0.0)
    hline.SetLineColor(kRed)
    hline.SetFillStyle(0)
    p2.cd()
    hratio.Draw('AXIS')
    hline.Draw('HIST,SAME')
    hratio.Draw('P,SAME')
    hratio.GetXaxis().SetRangeUser(xmin,xmax)
    hline.GetXaxis().SetRangeUser(xmin,xmax)
    
    p1.Update()
    p2.Update()
    c1.Update()

    #-------- draw pave tex
    pt =TLatex() #(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC")
    pt.SetNDC()
    pt.SetTextAlign(12)
    pt.SetTextFont(42)
    pt.SetTextSize(0.03)

    p1.cd()
    pt.DrawLatex(0.20,0.97,"CMS Preliminary")
    pt.DrawLatex(0.55,0.97,"#sqrt{s} = 13 TeV, #intLdt = "+"{:.3}".format(float(lumi))+" fb^{-1}")
    pt.SetTextSize(0.04)
    pt.DrawLatex(0.25,0.86, notes)
    p1.SetLogy()
    p1.RedrawAxis()
    p1.Update()

    if channel!="":
        channel_tex="ee" if channel.Contains("el") else "#mu#mu"
        pt.DrawLatex(0.25,0.82, channel_tex+" channel")
            
    fout.cd()
    c1.Write()
    hratio.Write()
    hstack.Write()
    fout.Close()

    c1.Print(outDir+'/'+output+'_'+channel.Data()+'.eps')
    
    return


def drawCompare(hstack, hratio, legend,
                hstack_opt="nostack",outdir="./",tag="test",
                xmin=50., xmax=500., xtitle="" ,units="",
                lumi=2.169, notes=""):

    fout = TFile(outdir+'/'+tag+'.root', 'recreate')
    
    c1 = TCanvas(tag+"_c1", "c1", 600, 750); c1.Draw()
    c1.SetWindowSize(600 + (600 - c1.GetWw()), (750 + (750 - c1.GetWh())))
    p1 = TPad(tag+"_pad1","pad1",0,0.25,1,0.99)
    p1.SetBottomMargin(0.15)
    p1.SetLeftMargin(0.15)
    p1.Draw()
    p2 = TPad(tag+"_pad2","pad2",0,0,1,0.25)
    p2.SetTopMargin(0.03)
    p2.SetBottomMargin(0.3)
    p2.SetLeftMargin(0.15)
    p2.SetFillStyle(0)
    p2.Draw()
    
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)
    
    p1.cd()
    hstack.Draw(hstack_opt)
    legend.Draw("same")
        
    hstack.SetMinimum(0.01)
    hstack.GetHistogram().GetXaxis().SetRangeUser(xmin,xmax)
    
    maxi=hstack.GetHistogram().GetXaxis().GetXmax()
    mini=hstack.GetHistogram().GetXaxis().GetXmin()
    bins=hstack.GetHistogram().GetNbinsX()
    if len(units)>0:
        hstack.GetHistogram().GetXaxis().SetTitle(xtitle + " (" +units+")")
        hstack.GetHistogram().GetYaxis().SetTitle("Events / "+str((maxi-mini)/bins)+ " "+units)
    else:
        hstack.GetHistogram().GetXaxis().SetTitle(xtitle)
        hstack.GetHistogram().GetYaxis().SetTitle("Events")

    hratio.SetName(tag+'_'+'hratio')
    hline = hratio.Clone(tag+'_'+"hline")
    for ii in range(hline.GetNbinsX()+1):
        hline.SetBinContent(ii,1.0)
        hline.SetBinError(ii,0.0)
    hline.SetLineColor(kRed)
    hline.SetFillStyle(0)
    p2.cd()
    hratio.Draw('AXIS')
    hline.Draw('HIST,SAME')
    hratio.Draw('P,SAME')
    hratio.GetXaxis().SetRangeUser(xmin,xmax)

    p1.SetLogy()
    p1.RedrawAxis()
    p1.Update()
    p2.Update()
    c1.Update()
        
    #-------- draw pave tex
    pt =TLatex() #(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC")
    pt.SetNDC()
    pt.SetTextAlign(12)
    pt.SetTextFont(42)
    pt.SetTextSize(0.03)
    
    p1.cd()
    pt.DrawLatex(0.19,0.97,"CMS Preliminary")
    pt.DrawLatex(0.55,0.97,"#sqrt{s} = 13 TeV, #intLdt = "+"{:.3}".format(float(lumi))+" fb^{-1}")
    pt.SetTextSize(0.04)
    pt.DrawLatex(0.45,0.6, notes)
       
    fout.cd()
    c1.Write()
    hratio.Write()
    hstack.Write()
    fout.Close()
    
    c1.Print(outdir+'/'+tag+'.eps')
    
    return

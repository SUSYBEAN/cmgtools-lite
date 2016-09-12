#!/bin/env python

#import optparse
import os
import ROOT 
from Setup import *


#parser = optparse.OptionParser()
#parser.add_option("-m","--masses",dest="masses",default='limitPlot',help="Limit plot")


tag='xzz2l2nu'

indir='results'
outdir='plots'

OptCut='mtcut'

log=True

masses = [
         600,
         800,
          1000,
         1200,
          1400,
         1600,1800,
          2000,
         2500,
          3000,
         3500,
         4000,
#         4500
         ]

cuts=mtcuts
if OptCut=='mtcut': cuts=mtcuts
elif OptCut=='zptcut': cuts=zptcuts
elif OptCut=='metcut': cuts=metcuts



def GetLimitsInFile(filename,data,mass,cut):
    infile = ROOT.TFile(filename)
    tree = infile.Get('limit')
    if not (mass in data.keys()):
        data[mass]={}
    if not (cut in data[mass].keys()):
        data[mass][cut]={}
    for event in tree:
        if event.quantileExpected<0:
            data[mass][cut]['obs']=event.limit
        if event.quantileExpected>0.02 and event.quantileExpected<0.03:
            data[mass][cut]['-2sigma']=event.limit
        if event.quantileExpected>0.15 and event.quantileExpected<0.17:
            data[mass][cut]['-1sigma']=event.limit
        if event.quantileExpected>0.49 and event.quantileExpected<0.51:
            data[mass][cut]['exp']=event.limit
        if event.quantileExpected>0.83 and event.quantileExpected<0.85:
            data[mass][cut]['+1sigma']=event.limit
        if event.quantileExpected>0.974 and event.quantileExpected<0.976:
            data[mass][cut]['+2sigma']=event.limit
    infile.Close()


def PlotLimitsOfMassPoint(data,mass,label,xtitle,ytitle,xmin,xmax,ymin,ymax,log,plot):
    band68=ROOT.TGraphAsymmErrors()
    band68.SetName("band68")
    band95=ROOT.TGraphAsymmErrors()
    band95.SetName("band95")
    line_plus1=ROOT.TGraph()
    line_plus1.SetName("line_plus1")
    line_plus2=ROOT.TGraph()
    line_plus2.SetName("line_plus2")
    line_minus1=ROOT.TGraph()
    line_minus1.SetName("line_minus1")
    line_minus2=ROOT.TGraph()
    line_minus2.SetName("line_minus2")

    N=0
    for cut,info in data[mass].iteritems():
        band68.SetPoint(N,cut,info['exp'])
        band95.SetPoint(N,cut,info['exp'])
        line_plus1.SetPoint(N,cut,info['+1sigma'])
        line_plus2.SetPoint(N,cut,info['+2sigma'])
        line_minus1.SetPoint(N,cut,info['-1sigma'])
        line_minus2.SetPoint(N,cut,info['-2sigma'])
        band68.SetPointError(N,0.0,0.0,info['exp']-info['-1sigma'],info['+1sigma']-info['exp'])
        band95.SetPointError(N,0.0,0.0,info['exp']-info['-2sigma'],info['+2sigma']-info['exp'])
        N=N+1

    band68.Sort()
    band95.Sort()
    line_plus1.Sort()    
    line_plus2.Sort()    
    line_minus1.Sort()    
    line_minus2.Sort()    

    #plotting information
    c=ROOT.TCanvas("c","c",600,600)
    #c.SetWindowSize(600 + (600 - c.GetWw()), (750 + (750 - c.GetWh())))
    c.SetGrid()
    c.SetBottomMargin(0.15)
    c.SetLeftMargin(0.15)
    frame=c.DrawFrame(xmin,ymin,xmax,ymax)
    frame.GetXaxis().SetTitle(xtitle)
    frame.GetYaxis().SetTitle(ytitle)
    frame.GetXaxis().SetLabelFont(42)
    frame.GetXaxis().SetLabelOffset(0.007)
    frame.GetXaxis().SetLabelSize(0.03)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleOffset(1.15)
    frame.GetXaxis().SetTitleFont(42)
    frame.GetYaxis().SetLabelFont(42)
    frame.GetYaxis().SetLabelOffset(0.007)
    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.4)
    frame.GetYaxis().SetTitleFont(42)
    frame.GetZaxis().SetLabelFont(42)
    frame.GetZaxis().SetLabelOffset(0.007)
    frame.GetZaxis().SetLabelSize(0.045)
    frame.GetZaxis().SetTitleSize(0.05)
    frame.GetZaxis().SetTitleFont(42)

    band68.SetFillColor(ROOT.kGreen)
    band68.SetLineWidth(2)
    band68.SetLineColor(ROOT.kBlack)
    band68.SetLineStyle(2)
    band68.SetMarkerStyle(20)
 
    band95.SetFillColor(ROOT.kYellow)
    band95.SetLineStyle(2)

    line_plus1.SetLineWidth(1)
    line_plus1.SetLineColor(ROOT.kGreen-3)

    line_plus2.SetLineWidth(1)
    line_plus2.SetLineColor(ROOT.kYellow-6)

    line_minus1.SetLineWidth(1)
    line_minus1.SetLineColor(ROOT.kGreen-3)

    line_minus2.SetLineWidth(1)
    line_minus2.SetLineColor(ROOT.kYellow-6)

    c.cd()
    frame.Draw()
    band95.Draw("3same")
    band68.Draw("3same")
    band68.Draw("XLsame")
    line_plus1.Draw("Lsame")
    line_plus2.Draw("Lsame")
    line_minus1.Draw("Lsame")
    line_minus2.Draw("Lsame")

    c.SetLogy(log)
    c.Draw()

    #legend = ROOT.TLegend(0.62,0.3,0.92,0.60,"","brNDC")
    legend = ROOT.TLegend(0.32,0.6,0.62,0.80,"","brNDC")
    legend.SetName('legend')
    legend.SetBorderSize(0)
    legend.SetLineColor(1)
    legend.SetLineStyle(1)
    legend.SetLineWidth(1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetFillColor(ROOT.kWhite)

    legend.AddEntry(band68,'Expected', 'l')
    legend.AddEntry(band68,'Expected #pm 1 #sigma', 'f')
    legend.AddEntry(band95,'Expected #pm 2 #sigma', 'f')

    legend.Draw('same')

    pt =ROOT.TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC")
    pt.SetBorderSize(0)
    pt.SetTextAlign(12)
    pt.SetFillStyle(0)
    pt.SetTextFont(42)
    pt.SetTextSize(0.03)
    text = pt.AddText(0.01,0.3,label)
    pt.Draw()
    c.Draw()
    c.RedrawAxis()

    c.SaveAs(plot+'.pdf')
    c.SaveAs(plot+'.eps')
    c.SaveAs(plot+'.png')



  
data = {}

# for all mass points
for mass in masses:
    for cut in cuts[mass]:
        zptcut0=zptcuts[mass][0]
        metcut0=metcuts[mass][0]
        mtcut0=mtcuts[mass][0]
        if OptCut=='mtcut': mtcut0=cut
        elif OptCut=='zptcut': zptcut0=cut
        elif OptCut=='metcut': metcut0=cut
        # card tag
        #cardTag = '{tag}_m{mass:n}_zptcut{zptcut:n}_metcut{metcut:n}_mtcut{mtcut:n}'.format(tag=tag,mass=mass,zptcut=zptcut0,metcut=metcut0,mtcut=mtcut0)
        cardTag = '{tag}_m{mass:n}_zpt{zptcut:n}_metcut{metcut:n}_mtcut{mtcut:n}'.format(tag=tag,mass=mass,zptcut=zptcut0,metcut=metcut0,mtcut=mtcut0)
        # rootfile
        infilename = '{indir}/higgsCombine{cardTag}.Asymptotic.mH{mass:n}.root'.format(indir=indir,cardTag=cardTag,mass=mass)
        #
        GetLimitsInFile(infilename,data,mass,cut) 

    # end loop cuts
# end loop masses
      

#print data

for mass in masses:
    zptcut0=zptcuts[mass][0]
    metcut0=metcuts[mass][0]
    mtcut0=mtcuts[mass][0]
    xmin=cuts[mass][0]-10
    xmax=cuts[mass][-1]+10
    if log:
        ymin=min( x['exp'] for x in data[mass].values() )*0.1
        ymax=max( x['exp'] for x in data[mass].values() )*100
    else:
        ymin=min( x['exp'] for x in data[mass].values() )*0.2
        ymax=max( x['exp'] for x in data[mass].values() )*2.0
    cut='M_{T}'
    if OptCut=='mtcut': cut='M_{T}'
    elif OptCut=='zptcut': cut='P_{T}(Z)'
    elif OptCut=='metcut': cut='MET'
    label = cut+' cut optimization for M_{G}='+str(mass)+' GeV'
    xtitle = cut+' Cut (GeV)'
    ytitle='#sigma_{95%}/#sigma_{th.}'
    plot=outdir+'/'+'plotCutOpt_m'+str(mass)+'_zpt'
    plot='{outdir}/plotCutOpt_m{mass:n}_{optcut}_zpt{zptcut:n}_met{metcut:n}_mt{mtcut:n}'.format(outdir=outdir,mass=mass,zptcut=zptcut0,metcut=metcut0,mtcut=mtcut0,optcut=OptCut)
    PlotLimitsOfMassPoint(data,mass,label=label,xtitle=xtitle,ytitle=ytitle,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,log=log,plot=plot)






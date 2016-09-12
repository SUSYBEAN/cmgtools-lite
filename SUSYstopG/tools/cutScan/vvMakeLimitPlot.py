#!/usr/bin/env python

import ROOT
import optparse


ROOT.gROOT.ProcessLine('.x tdrstyle.C')

parser = optparse.OptionParser()
parser.add_option("-o","--output",dest="output",default='limitPlot',help="Limit plot")

parser.add_option("-x","--minX",dest="minX",type=float,help="minimum x",default=1000.0)
parser.add_option("-X","--maxX",dest="maxX",type=float,help="maximum x",default=5000.0)
parser.add_option("-y","--minY",dest="minY",type=float,help="minimum y",default=0.00001)
parser.add_option("-Y","--maxY",dest="maxY",type=float,help="maximum y",default=1000)
parser.add_option("-b","--blind",dest="blind",type=int,help="Not do observed ",default=1)
parser.add_option("-l","--log",dest="log",type=int,help="Log plot",default=1)

parser.add_option("-t","--titleX",dest="titleX",default='M_{Grav} (GeV)',help="title of x axis")
parser.add_option("-T","--titleY",dest="titleY",default='#sigma_{95\%}/#sigma_{th.} ',help="title of y axis")

parser.add_option("-p","--paveText",dest="label",default='CMS Preliminary , L = 2.17 fb^{-1}',help="label")

parser.add_option("-s","--scale",dest="scale",type=float,help="scale factor, such as 1/BR to get the xsec before decay",default=1)

parser.add_option("--drawXsec",dest="drawXsec", type=int,default=1,help="plot xsec")





(options,args) = parser.parse_args()
#define output dictionary


f=ROOT.TFile(args[0])
limit=f.Get("limit")
data={}

BulkGZZ2l2nuXsec = {
600:8.61578e-03,
800:1.57965e-03,
1000:4.21651e-04,
1200:1.39919e-04,
1400:5.32921e-05,
1600:2.24428e-05,
1800:1.01523e-05,
2000:4.86037e-06,
2500:9.08739e-07,
3000:1.98856e-07,
3500:4.87505e-08,
4000:1.25937e-08,
}

for event in limit:
    if not (event.mh in data.keys()):
        data[event.mh]={}

    print 'processing',event.mh,event.quantileExpected,event.limit     

    if event.quantileExpected<0 and options.blind==0:            
        data[event.mh]['obs']=event.limit
    if event.quantileExpected>0.02 and event.quantileExpected<0.03:            
        data[event.mh]['-2sigma']=event.limit
    if event.quantileExpected>0.15 and event.quantileExpected<0.17:            
        data[event.mh]['-1sigma']=event.limit
    if event.quantileExpected>0.49 and event.quantileExpected<0.51:            
        data[event.mh]['exp']=event.limit
    if event.quantileExpected>0.83 and event.quantileExpected<0.85:            
        data[event.mh]['+1sigma']=event.limit
    if event.quantileExpected>0.974 and event.quantileExpected<0.976:            
        data[event.mh]['+2sigma']=event.limit

# apply scale
for mass in data.keys():
    for key in data[mass].keys():
        data[mass][key] *= options.scale

# calculate xsec limit
if options.drawXsec:
    for mass in data.keys():
        for key in data[mass].keys():
            data[mass][key] *= BulkGZZ2l2nuXsec[mass]

print data

# signal xsec to be plotted
sig_xsec = BulkGZZ2l2nuXsec

# apply scale
for mass in sig_xsec.keys():
    sig_xsec[mass] *= options.scale

# if not drawXsec, plot the signal strength, which is 1 as expected.
if options.drawXsec==0:
    for mass in sig_xsec.keys():
        sig_xsec[mass] /= sig_xsec[mass]


band68=ROOT.TGraphAsymmErrors()
band68.SetName("band68")
band95=ROOT.TGraphAsymmErrors()
band95.SetName("band95")
bandObs=ROOT.TGraph()
bandObs.SetName("bandObs")

line_plus1=ROOT.TGraph()
line_plus1.SetName("line_plus1")

line_plus2=ROOT.TGraph()
line_plus2.SetName("line_plus2")

line_minus1=ROOT.TGraph()
line_minus1.SetName("line_minus1")

line_minus2=ROOT.TGraph()
line_minus2.SetName("line_minus2")

sigXsec=ROOT.TGraph()
sigXsec.SetName("sigXsec")

N=0
for mass,info in data.iteritems():
    band68.SetPoint(N,mass,info['exp'])
    band95.SetPoint(N,mass,info['exp'])
    line_plus1.SetPoint(N,mass,info['+1sigma'])
    line_plus2.SetPoint(N,mass,info['+2sigma'])
    line_minus1.SetPoint(N,mass,info['-1sigma'])
    line_minus2.SetPoint(N,mass,info['-2sigma'])

    if options.blind==0: bandObs.SetPoint(N,mass,info['obs'])
    band68.SetPointError(N,0.0,0.0,info['exp']-info['-1sigma'],info['+1sigma']-info['exp'])
    band95.SetPointError(N,0.0,0.0,info['exp']-info['-2sigma'],info['+2sigma']-info['exp'])

    if mass in BulkGZZ2l2nuXsec.keys():
        sigXsec.SetPoint(N,mass,sig_xsec[mass])

    N=N+1


band68.Sort()
band95.Sort()
bandObs.Sort()
line_plus1.Sort()    
line_plus2.Sort()    
line_minus1.Sort()    
line_minus2.Sort()    
sigXsec.Sort()


#plotting information

c=ROOT.TCanvas("c","c",600,600)
#c.SetWindowSize(600 + (600 - c.GetWw()), (750 + (750 - c.GetWh())))
c.SetGrid()
c.SetBottomMargin(0.15)
c.SetLeftMargin(0.15)
frame=c.DrawFrame(options.minX,options.minY,options.maxX,options.maxY)
frame.GetXaxis().SetTitle(options.titleX)
frame.GetYaxis().SetTitle(options.titleY)
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

bandObs.SetLineWidth(2)
#bandObs.SetLineColor(ROOT.kRed)

line_plus1.SetLineWidth(1)
line_plus1.SetLineColor(ROOT.kGreen-3)

line_plus2.SetLineWidth(1)
line_plus2.SetLineColor(ROOT.kYellow-6)

line_minus1.SetLineWidth(1)
line_minus1.SetLineColor(ROOT.kGreen-3)

line_minus2.SetLineWidth(1)
line_minus2.SetLineColor(ROOT.kYellow-6)

sigXsec.SetLineWidth(1)
sigXsec.SetLineColor(ROOT.kRed)

c.cd()
frame.Draw()
band95.Draw("3same")
band68.Draw("3same")
band68.Draw("XLsame")
line_plus1.Draw("Lsame")
line_plus2.Draw("Lsame")
line_minus1.Draw("Lsame")
line_minus2.Draw("Lsame")
sigXsec.Draw("Lsame")

c.SetLogy(options.log)
c.Draw()

if options.drawXsec==0:
    legend = ROOT.TLegend(0.62,0.3,0.92,0.60,"","brNDC")
else: 
    legend = ROOT.TLegend(0.62,0.5,0.92,0.8,"","brNDC")
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
if options.blind==0:
    legend.AddEntry(bandObs,'Observed', 'pl')
legend.AddEntry(sigXsec, "G_{Bulk}, #tilde{k}=0.5", 'l')

legend.Draw('same')

pt =ROOT.TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC")
pt.SetBorderSize(0)
pt.SetTextAlign(12)
pt.SetFillStyle(0)
pt.SetTextFont(42)
pt.SetTextSize(0.03)
text = pt.AddText(0.01,0.3,options.label)
pt.Draw()
c.Draw()
c.RedrawAxis()

if options.blind==0:
    bandObs.Draw("plsame")



fout=ROOT.TFile(options.output+'.root',"RECREATE")
fout.cd()

c.SaveAs(options.output+'.pdf')
c.Write()
band68.Write()
band95.Write()
bandObs.Write()
line_plus1.Write()    
line_plus2.Write()    
line_minus1.Write()    
line_minus2.Write()    

fout.Close()
f.Close()



#!/usr/bin/env python

from ROOT import *
import os
from Fitter import *

ROOT.gROOT.ProcessLine('.x tdrstyle.C') 

sigPlotters=[]
sigSamples = [
'BulkGravToZZToZlepZinv_narrow_600', 
#'BulkGravToZZToZlepZinv_narrow_800', 
#'BulkGravToZZToZlepZinv_narrow_1000', 
#'BulkGravToZZToZlepZinv_narrow_1200', 
#'BulkGravToZZToZlepZinv_narrow_1400', 
#'BulkGravToZZToZlepZinv_narrow_1600', 
#'BulkGravToZZToZlepZinv_narrow_1800', 
#'BulkGravToZZToZlepZinv_narrow_2000', 
#'BulkGravToZZToZlepZinv_narrow_2500', 
#'BulkGravToZZToZlepZinv_narrow_3000', 
#'BulkGravToZZToZlepZinv_narrow_3500', 
#'BulkGravToZZToZlepZinv_narrow_4000', 
#'BulkGravToZZToZlepZinv_narrow_4500', 
]


sigXsec = {
'BulkGravToZZToZlepZinv_narrow_600' : 1,
'BulkGravToZZToZlepZinv_narrow_800' : 1,
'BulkGravToZZToZlepZinv_narrow_1000' : 1,
'BulkGravToZZToZlepZinv_narrow_1200' : 1,
'BulkGravToZZToZlepZinv_narrow_1400' : 1,
'BulkGravToZZToZlepZinv_narrow_1600' : 1,
'BulkGravToZZToZlepZinv_narrow_1800' : 1,
'BulkGravToZZToZlepZinv_narrow_2000' : 1,
'BulkGravToZZToZlepZinv_narrow_2500' : 1,
'BulkGravToZZToZlepZinv_narrow_3000' : 1,
'BulkGravToZZToZlepZinv_narrow_3500' : 1,
'BulkGravToZZToZlepZinv_narrow_4000' : 1,
'BulkGravToZZToZlepZinv_narrow_4500' : 1,
}
sigMassWindow = {
'BulkGravToZZToZlepZinv_narrow_600' : [60,200,800],
'BulkGravToZZToZlepZinv_narrow_800' : [60,400,1000],
'BulkGravToZZToZlepZinv_narrow_1000' : [80,400,1200],
'BulkGravToZZToZlepZinv_narrow_1200' : [100,400,1400],
'BulkGravToZZToZlepZinv_narrow_1400' : [120,400,1600],
'BulkGravToZZToZlepZinv_narrow_1600' : [140,400,1800],
'BulkGravToZZToZlepZinv_narrow_1800' : [160,400,2000],
'BulkGravToZZToZlepZinv_narrow_2000' : [180,400,2200],
'BulkGravToZZToZlepZinv_narrow_2500' : [200,400,2800],
'BulkGravToZZToZlepZinv_narrow_3000' : [200,400,3500],
'BulkGravToZZToZlepZinv_narrow_3500' : [200,400,4000],
'BulkGravToZZToZlepZinv_narrow_4000' : [200,400,4500],
'BulkGravToZZToZlepZinv_narrow_4500' : [200,400,5000],
}


fin = ROOT.TFile('study_mt.root')
fout = ROOT.TFile("fit_mt.root", 'recreate')

c1 = ROOT.TCanvas("c1","c1")
c1.SetBottomMargin(0.15)
c1.SetLeftMargin(0.15)
c1.Print("fit_mt.ps[")

# initialize fitter
fitter = Fitter(poi = ['x'])

h1 = []
h2 = []
h3 = []
h4 = []
lg = []
frame = []
func1 = []
func2 = []
func3 = []
func4 = []
data1 = []
data2 = []
data3 = []
data4 = []

for idx,sample in enumerate(sigSamples):
    nbins = sigMassWindow[sample][0]
    xmin = sigMassWindow[sample][1]
    xmax = sigMassWindow[sample][2]
    h1.append(fin.Get(sample+'_hmt'))
    h2.append(fin.Get(sample+'_hmta'))
    h3.append(fin.Get(sample+'_hmtb'))
    h4.append(fin.Get(sample+'_hmtc'))

    lg.append(ROOT.TLegend(0.2,0.5,0.5,0.8))
    lg[idx].SetName(sample+'_lg')
    lg[idx].AddEntry(h1[idx],"MT","pl")
    lg[idx].AddEntry(h2[idx],"MTa","pl")
    lg[idx].AddEntry(h3[idx],"MTb","pl")
    lg[idx].AddEntry(h4[idx],"MTc","pl")

    frame.append(fitter.getFrame(name='frame_'+sample+'_mt',poi='x'))
    frame[idx].SetTitle('MT of '+sample) 
    frame[idx].SetAxisRange(xmin,xmax,'X')
    frame[idx].SetXTitle('M_{T} (GeV)')
    frame[idx].SetYTitle('Events')

    f_name_1 = 'f_doubleCB_'+sample+'_mt'
    h_name_1 = 'h_'+sample+'_mt'
    fitter.doubleCB(name=f_name_1)
    fitter.importBinnedData(h1[idx],name=h_name_1)
    fitter.fit(model=f_name_1,data=h_name_1)
    func1.append(fitter.getModel(f_name_1))
    data1.append(fitter.getData(h_name_1))
    data1[idx].plotOn(frame[idx],RooFit.DrawOption('hist'),RooFit.LineColor(2))
    func1[idx].plotOn(frame[idx],RooFit.DrawOption('hist'),RooFit.LineColor(2))

    f_name_2 = 'f_doubleCB_'+sample+'_mta'
    h_name_2 = 'h_'+sample+'_mta'
    fitter.doubleCB(name=f_name_2)
    fitter.importBinnedData(h2[idx],name=h_name_2)
    fitter.fit(model=f_name_2,data=h_name_2)
    func2.append(fitter.getModel(f_name_2))
    data2.append(fitter.getData(h_name_2))
    data2[idx].plotOn(frame[idx],RooFit.DrawOption('hist'),RooFit.LineColor(4))
    func2[idx].plotOn(frame[idx],RooFit.DrawOption('hist'),RooFit.LineColor(4))

    f_name_3 = 'f_doubleCB_'+sample+'_mtb'
    h_name_3 = 'h_'+sample+'_mtb'
    fitter.doubleCB(name=f_name_3)
    fitter.importBinnedData(h3[idx],name=h_name_3)
    fitter.fit(model=f_name_3,data=h_name_3)
    func3.append(fitter.getModel(f_name_3))
    data3.append(fitter.getData(h_name_3))
    data3[idx].plotOn(frame[idx],RooFit.DrawOption('hist'),RooFit.LineColor(6))
    func3[idx].plotOn(frame[idx],RooFit.DrawOption('hist'),RooFit.LineColor(6))

    f_name_4 = 'f_doubleCB_'+sample+'_mtc'
    h_name_4 = 'h_'+sample+'_mtc'
    fitter.doubleCB(name=f_name_4)
    fitter.importBinnedData(h4[idx],name=h_name_4)
    fitter.fit(model=f_name_4,data=h_name_4)
    func4.append(fitter.getModel(f_name_4))
    data4.append(fitter.getData(h_name_4))
    data4[idx].plotOn(frame[idx],RooFit.DrawOption('hist'),RooFit.LineColor(8))
    func4[idx].plotOn(frame[idx],RooFit.DrawOption('hist'),RooFit.LineColor(8))



    c1.Clear()
    frame[idx].Draw()
    lg[idx].Draw()
    c1.Print("fit_mt.ps")
    c1.Clear()

    fout.cd()
    data1[idx].Write()
    data2[idx].Write()
    data3[idx].Write()
    data4[idx].Write()
    func1[idx].Write()
    func2[idx].Write()
    func3[idx].Write()
    func4[idx].Write()
    frame[idx].Write()
    lg[idx].Write()

c1.Print("fit_mt.ps]")

os.system('ps2pdf fit_mt.ps')
fout.Close()



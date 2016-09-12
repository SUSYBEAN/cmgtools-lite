#! /usr/bin/env python
from ROOT import *
from array import array
import sys
if len(sys.argv)<3:
    print 'geteff.py rootfile histname [title]'
    sys.exit()

def geteff(a,b,n=4,nm=''):
    print 'getting efficiency',nm
    if a.GetSize()!=b.GetSize():
        print 'ERROR: Different bin number!'
        return TH1F()
    e=a.Clone(nm)
    for i in range(a.GetSize()):
        a0=a.GetBinContent(i)
        b0=b.GetBinContent(i)
        if b0:
            e.SetBinContent(i,100.*a0/b0)
            a1=b0-a0
            if a1<0: print i,b0,a0
            e.SetBinError(i,abs(a1*a0/b0**3)**.5*100)
#        print e.GetBinLowEdge(i),e.GetBinContent(i),a0,b0
    e.SetMarkerColor(n)
    e.SetLineColor(n)
    e.SetTitle('efficiency;Pt/GeV;%')
    return e

def getSF(a,b,n=4,nm=''):
    print 'getting scale factor',nm
    if a.GetSize()!=b.GetSize():
        print 'ERROR: Different bin number!'
        return TH1F()
    e=a.Clone(nm)
    for i in range(1,a.GetSize()-1):
        a0=a.GetBinContent(i)
        b0=b.GetBinContent(i)
        if a0 and b0:
            e.SetBinContent(i,a0/b0)
            e.SetBinError(i,((a.GetBinError(i)/a.GetBinContent(i))**2+(b.GetBinError(i)/b.GetBinContent(i))**2)**.5*e.GetBinContent(i))
    e.SetMarkerColor(n)
    e.SetLineColor(n)
    e.SetTitle('SF;Pt/GeV;')
    return e

if __name__=='__main__':
    gROOT.SetBatch()
    gStyle.SetOptStat(000000)
    gStyle.SetLegendBorderSize(0)
    t1=TFile(sys.argv[1])
    fo=TFile('output.root','update')
    hs1=t1.Get('hnp')
    hs2=t1.Get('hnm')
    ez1=geteff(hs1,hs2,nm=sys.argv[2])
    if len(sys.argv)>3:    ez1.SetTitle(sys.argv[3])
    ez1.Write(sys.argv[2],2)
    c=TCanvas()
    if 'TH2' in str(type(ez1)):ez1.Draw('colz')
    else:ez1.Draw()
    c.Print('plots/'+sys.argv[2]+'.png')
    fo.Close()
    t1.Close()


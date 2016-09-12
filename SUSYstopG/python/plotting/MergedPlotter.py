import ROOT
import sys
from array import array
import pickle
from CMGTools.XZZ2l2nu.plotting.PlotterBase import PlotterBase
class MergedPlotter(PlotterBase):

    def __init__(self,plotters):
        super(MergedPlotter,self).__init__()
        self.plotters=plotters
        self.corrFactors=plotters[0].corrFactors

    def applySmoothing(self):
        for plotter in self.plotters:
            plotter.applySmoothing()


    def drawTH1(self,name,var,cuts,lumi,bins,min,max,titlex = "",units = "",drawStyle = "HIST"):
        h=None
        for plotter in self.plotters:
            if h is None:
                h=plotter.drawTH1(name,var,cuts,lumi,bins,min,max,titlex,units,drawStyle)
            else:
                h.Add(plotter.drawTH1(name,var,cuts,lumi,bins,min,max,titlex,units,drawStyle))
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        h.GetXaxis().SetTitle(titlex)
        if units != '': h.GetXaxis().SetTitle(titlex+ " ("+units+")")
        return h

    def drawTH2(self,name,var,cuts,lumi,binsx,minx,maxx,binsy,miny,maxy,titlex = "",unitsx = "",titley = "",unitsy = "",drawStyle = "COLZ"):
        h=None
        for plotter in self.plotters:
            if h is None:
                h=plotter.drawTH2(name,var,cuts,lumi,binsx,minx,maxx,binsy,miny,maxy,titlex,unitsx,titley,unitsy,drawStyle)
            else:
                h.Add(plotter.drawTH2(name,var,cuts,lumi,binsx,minx,maxx,binsy,miny,maxy,titlex,unitsx,titley,unitsy,drawStyle))

#        h.SetLineStyle(self.linestyle)
#        h.SetLineColor(self.linecolor)
#        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        h.GetXaxis().SetTitle(titlex)
        if unitsx != '' : h.GetXaxis().SetTitle(titlex+ " ("+unitsx+")")
        h.GetYaxis().SetTitle(titley)
        if unitsy != '' : h.GetYaxis().SetTitle(titley+ " ("+unitsy+")")
        return h


    def drawProfile(self,name,var,cuts,lumi,binsx,minx,maxx,miny,maxy,titlex = "",unitsx = "",titley = "",unitsy = "",drawStyle = "COLZ"):
        h=None
        for plotter in self.plotters:
            if h is None:
                h=plotter.drawProfile(name,var,cuts,lumi,binsx,minx,maxx,miny,maxy,titlex,unitsx,titley,unitsy,drawStyle)
            else:
                h.Add(plotter.drawProfile(name,var,cuts,lumi,binsx,minx,maxx,miny,maxy,titlex,unitsx,titley,unitsy,drawStyle))

#        h.SetLineStyle(self.linestyle)
#        h.SetLineColor(self.linecolor)
#        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
        return h


    def drawTH3(self,name,var,cuts,lumi,binsx,minx,maxx,binsy,miny,maxy,binsz,minz,maxz,titlex = "",unitsx = "",titley = "",unitsy = "",drawStyle = "COLZ"):
        h=None
        for plotter in self.plotters:
            if h is None:
                h=plotter.drawTH3(name,var,cuts,lumi,binsx,minx,maxx,binsy,miny,maxy,binsz,minz,maxz,titlex,unitsx,titley,unitsy,drawStyle)
            else:
                h.Add(plotter.drawTH3(name,var,cuts,lumi,binsx,minx,maxx,binsy,miny,maxy,binsz,minz,maxz,titlex,unitsx,titley,unitsy,drawStyle))

        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
        return h


    def drawTH2Binned(self,name,var,cuts,lumi,binningx,binningy,titlex = "",unitsx = "",titley = "",unitsy = "",drawStyle = "COLZ"):
        h=None
        for plotter in self.plotters:
            if h is None:
                h=plotter.drawTH2Binned(name,var,cuts,lumi,binningx,binningy,titlex,unitsx,titley,unitsy,drawStyle)
            else:
                h.Add(plotter.drawTH2Binned(name,var,cuts,lumi,binningx,binningy,titlex,unitsx,titley,unitsy,drawStyle))

#        h.SetLineStyle(self.linestyle)
#        h.SetLineColor(self.linecolor)
#        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
        return h

    def drawTH1Binned(self,name,var,cuts,lumi,binningx,titlex = "",unitsx = "",drawStyle = "COLZ"):
        h=None
        for plotter in self.plotters:
            if h is None:
                h=plotter.drawTH1Binned(name,var,cuts,lumi,binningx,titlex,unitsx,drawStyle)
            else:
                h.Add(plotter.drawTH1Binned(name,var,cuts,lumi,binningx,titlex,unitsx,drawStyle))
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        return h


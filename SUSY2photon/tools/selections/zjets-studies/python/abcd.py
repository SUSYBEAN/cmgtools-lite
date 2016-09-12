#!/usr/bin/env python

import ROOT
import os,copy
from python.TreePlotter import TreePlotter
from python.MergedPlotter import MergedPlotter
from python.myStackPlotter import StackPlotter
from python.SimplePlot import *

from python.SetCuts import SetCuts


class abcdAnalyzer:
    def __init__(self, indir="../AnalysisRegion", outdir='plots',
                 lumi = 2.169126704526,  sepSig=True,
                 LogY=True,   doRatio=True):
        if not os.path.exists(outdir): os.system('mkdir '+outdir)

        self.mycuts= SetCuts()
        self.Channel=raw_input("Please choose a channel (el or mu): \n")
        self.outdir = outdir
        self.lumi = lumi
        self.sepSig = sepSig
        self.tex_dic = self.mycuts.Tex_dic
        self.whichregion=raw_input("Please choose a benchmarck Region (SR or VR): \n")
        
        self.cuts = self.mycuts.abcdCuts(self.Channel, self.whichregion)
        self.preCuts = self.mycuts.abcdCuts(self.Channel, self.whichregion, True)

        if self.whichregion=="VR":
            self.nbins, self.xMin, self.xMax = 10, 0, float(self.mycuts.met_pt)
        else:
            self.nbins, self.xMin, self.xMax = 10, 0, 500
                        
                        
        #######----------- Prepare samples to plot:
        zjetsPlotters=[]
        #zjetsSamples = ['DYJetsToLL_M50_HT100to200','DYJetsToLL_M50_HT200to400','DYJetsToLL_M50_HT400to600','DYJetsToLL_M50_HT600toInf']
        #zjetsSamples = ['DYJetsToLL_M50','DYJetsToLL_M50_Ext']
        zjetsSamples = ['DYJetsToLL_M50_BIG'] # M50_BIG = M50 + M50_Ext

        for sample in zjetsSamples:
            zjetsPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
            zjetsPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
            zjetsPlotters[-1].addCorrectionFactor('xsec','tree')
            zjetsPlotters[-1].addCorrectionFactor('genWeight','tree')
            zjetsPlotters[-1].addCorrectionFactor('puWeight','tree')
            zjetsPlotters[-1].addCorrectionFactor('triggersf','tree')
            zjetsPlotters[-1].addCorrectionFactor('llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf','tree')

        self.ZJets = MergedPlotter(zjetsPlotters)
        self.ZJets.setFillProperties(1001,ROOT.kGreen+2)


        wwPlotters=[]
        wwSamples = ['WWTo2L2Nu','WWToLNuQQ','WZTo1L1Nu2Q']
        
        for sample in wwSamples:
            wwPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
            wwPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
            wwPlotters[-1].addCorrectionFactor('xsec','tree')
            wwPlotters[-1].addCorrectionFactor('genWeight','tree')
            wwPlotters[-1].addCorrectionFactor('puWeight','tree')
            wwPlotters[-1].addCorrectionFactor('triggersf','tree')
            wwPlotters[-1].addCorrectionFactor('llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf','tree')
            
        self.WW = MergedPlotter(wwPlotters)
        self.WW.setFillProperties(1001,ROOT.kOrange)
            
            
        vvPlotters=[]
        vvSamples = ['WZTo2L2Q','WZTo3LNu',
                     'ZZTo2L2Nu',
                     'ZZTo2L2Q','ZZTo4L']
            
        for sample in vvSamples:
            vvPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
            vvPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
            vvPlotters[-1].addCorrectionFactor('xsec','tree')
            vvPlotters[-1].addCorrectionFactor('genWeight','tree')
            vvPlotters[-1].addCorrectionFactor('puWeight','tree')
            vvPlotters[-1].addCorrectionFactor('triggersf','tree')
            vvPlotters[-1].addCorrectionFactor('llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf','tree')
        
        self.VV = MergedPlotter(vvPlotters)
        self.VV.setFillProperties(1001,ROOT.kMagenta)
            
        wjetsPlotters=[]
        wjetsSamples = ['WJetsToLNu']
            
        for sample in wjetsSamples:
            wjetsPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
            wjetsPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
            wjetsPlotters[-1].addCorrectionFactor('xsec','tree')
            wjetsPlotters[-1].addCorrectionFactor('genWeight','tree')
            wjetsPlotters[-1].addCorrectionFactor('puWeight','tree')
            wjetsPlotters[-1].addCorrectionFactor('triggersf','tree')
            wjetsPlotters[-1].addCorrectionFactor('llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf','tree')

        self.WJets = MergedPlotter(wjetsPlotters)
        self.WJets.setFillProperties(1001,ROOT.kBlue-6)

        ttPlotters=[]
        ttSamples = ['TTTo2L2Nu']#,'TTZToLLNuNu','TTWJetsToLNu']

        for sample in ttSamples:
            ttPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
            ttPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
            ttPlotters[-1].addCorrectionFactor('xsec','tree')
            ttPlotters[-1].addCorrectionFactor('genWeight','tree')
            ttPlotters[-1].addCorrectionFactor('puWeight','tree')
            ttPlotters[-1].addCorrectionFactor('triggersf','tree')
            ttPlotters[-1].addCorrectionFactor('llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf','tree')
            
        self.TT = MergedPlotter(ttPlotters)
        self.TT.setFillProperties(1001,ROOT.kAzure-9)

        # --> define different background sets:
        nonZBGPlotters = []
        nonZBGSamples = wwSamples + vvSamples + wjetsSamples + ttSamples
        
        for sample in nonZBGSamples:
            nonZBGPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
            nonZBGPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
            nonZBGPlotters[-1].addCorrectionFactor('xsec','tree')
            nonZBGPlotters[-1].addCorrectionFactor('genWeight','tree')
            nonZBGPlotters[-1].addCorrectionFactor('puWeight','tree')
            nonZBGPlotters[-1].addCorrectionFactor('triggersf','tree')
            nonZBGPlotters[-1].addCorrectionFactor('llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf','tree')

        self.NonZBG = MergedPlotter(nonZBGPlotters)
        self.NonZBG.setFillProperties(1001,ROOT.kPink+6)

        nonResBGPlotters = []
        nonResBGSamples = wwSamples + wjetsSamples + ttSamples
        
        for sample in nonResBGSamples:
            nonResBGPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
            nonResBGPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
            nonResBGPlotters[-1].addCorrectionFactor('xsec','tree')
            nonResBGPlotters[-1].addCorrectionFactor('genWeight','tree')
            nonResBGPlotters[-1].addCorrectionFactor('puWeight','tree')
            nonResBGPlotters[-1].addCorrectionFactor('triggersf','tree')
            nonResBGPlotters[-1].addCorrectionFactor('llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf','tree')

        self.NonResBG = MergedPlotter(nonResBGPlotters)
        self.NonResBG.setFillProperties(1001,ROOT.kYellow)

        
        resBGPlotters = []
        resBGSamples = zjetsSamples + vvSamples

        for sample in resBGSamples:
            resBGPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
            resBGPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
            resBGPlotters[-1].addCorrectionFactor('xsec','tree')
            resBGPlotters[-1].addCorrectionFactor('genWeight','tree')
            resBGPlotters[-1].addCorrectionFactor('puWeight','tree')
            resBGPlotters[-1].addCorrectionFactor('triggersf','tree')
            resBGPlotters[-1].addCorrectionFactor('llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf','tree')

        self.ResBG = MergedPlotter(resBGPlotters)
        self.ResBG.setFillProperties(1001,ROOT.kRed)

        
        # --> Prepare the signal plotters:
        sigPlotters=[]
        sigSamples = [
            'BulkGravToZZToZlepZinv_narrow_800', 
            'BulkGravToZZToZlepZinv_narrow_1000', 
            'BulkGravToZZToZlepZinv_narrow_1200', 
        ]
        k=1000
        sigSampleNames = [
            str(k)+' x BulkG-800',
            str(k)+' x BulkG-1000',
            str(k)+' x BulkG-1200',
        ]
        sigXsec = {
            'BulkGravToZZToZlepZinv_narrow_800'  : 4.42472e-04*k,
            'BulkGravToZZToZlepZinv_narrow_1000' : 1.33926e-04*k,
            'BulkGravToZZToZlepZinv_narrow_1200' : 4.76544e-05*k,
        }

        for sample in sigSamples:
            sigPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
            sigPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
            sigPlotters[-1].addCorrectionFactor(str(sigXsec[sample]),'tree')
            sigPlotters[-1].addCorrectionFactor('genWeight','tree')
            sigPlotters[-1].addCorrectionFactor('puWeight','tree')
            sigPlotters[-1].addCorrectionFactor('triggersf','tree')
            sigPlotters[-1].addCorrectionFactor('llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf','tree')
            sigPlotters[-1].setFillProperties(0,ROOT.kWhite)

        # --> Prepare data plotters:    
        dataPlotters=[]
        dataSamples = ['SingleElectron_Run2015C_25ns_16Dec',
                       'SingleElectron_Run2015D_16Dec',
                       'SingleMuon_Run2015C_25ns_16Dec',
                       'SingleMuon_Run2015D_16Dec']
        for sample in dataSamples:
            dataPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
            
        self.Data = MergedPlotter(dataPlotters)
        self.Data.setFillProperties(1001,ROOT.kGreen+2)
        
        self.Stack = StackPlotter()
        self.Stack.addPlotter(self.Data, "data_obs", "Data", "data")
        self.Stack.addPlotter(self.WJets, "WJets","W+Jets", "background")
        self.Stack.addPlotter(self.WW, "WW","WW, WZ non-reson.", "background")
        self.Stack.addPlotter(self.TT, "TT","TT", "background")
        self.Stack.addPlotter(self.VV, "ZZ","ZZ, WZ reson.", "background")
        self.Stack.addPlotter(self.ZJets, "ZJets","Z+Jets", "background")
        
        for i in range(len(sigSamples)):
            sigPlotters[i].setLineProperties(2,ROOT.kRed+i,2)
            self.Stack.addPlotter(sigPlotters[i],sigSamples[i],sigSampleNames[i],'signal')  

 
        self.Stack.setLog(LogY)
        self.Stack.doRatio(doRatio)
        ROOT.gROOT.ProcessLine('.x tdrstyle.C')

    def GetStack(self):
        return self.Stack

#######----------- Start Plotting:
    def draw_preselection(self):

        tag='PreSelection_'+self.Channel+'_'
        #print self.Stack.log
        
        if self.Stack.log: tag=tag+'log_'
        
        self.Stack.drawStack('met_pt', self.preCuts, str(self.lumi*1000), self.nbins, self.xMin, self.xMax, titlex = "E_{T}^{miss}", units = "GeV",
                        output=tag+'met_low',outDir=self.outdir,separateSignal=self.sepSig,
                        drawtex=self.whichregion+" pre-selection", channel=self.Channel)
        self.Stack.drawStack('llnunu_l1_deltaR', self.preCuts, str(self.lumi*1000), 30,0,3, titlex = "#Delta R(l,l)", units = "",
                        output=tag+'dR_ll',outDir=self.outdir,separateSignal=self.sepSig,
                        drawtex=self.whichregion+" pre-selection", channel=self.Channel)
        self.Stack.drawStack('llnunu_deltaPhi', self.preCuts, str(self.lumi*1000), 35,0,3.5, titlex = "#Delta#Phi(Z,E_{T}^{miss})", units = "",
                        output=tag+'dPhi_llvv',outDir=self.outdir,separateSignal=self.sepSig,
                        drawtex=self.whichregion+" pre-selection", channel=self.Channel)
        self.Stack.drawStack('llnunu_l1_deltaPhi', self.preCuts, str(self.lumi*1000), 30,0,3, titlex = "#Delta#Phi(l,l)", units = "",
                        output=tag+'dPhi_ll',outDir=self.outdir,separateSignal=self.sepSig,
                        drawtex=self.whichregion+" pre-selection", channel=self.Channel)
        self.Stack.drawStack('llnunu_l1_mass', self.preCuts, str(self.lumi*1000), 50, 50, 150, titlex = "M_{ll}", units = "GeV",
                        output=tag+'zmass',outDir=self.outdir,separateSignal=self.sepSig,
                        drawtex=self.whichregion+" pre-selection", channel=self.Channel)
        self.Stack.drawStack('llnunu_mt', self.preCuts, str(self.lumi*1000), 30, 0.0, 600.0, titlex = "M_{T}^{ll#nu#nu}", units = "GeV",
                        output=tag+'mt_low',outDir=self.outdir,separateSignal=self.sepSig,
                        drawtex=self.whichregion+" pre-selection", channel=self.Channel)
        #self.Stack.drawStack('llnunu_mt', self.preCuts, str(self.lumi*1000), 60, 0.0, 1200.0, titlex = "M_{T}^{ll#nu#nu}", units = "GeV", output=tag+'mt',outDir=self.outdir,separateSignal=self.sepSig,drawtex=self.whichregion+" pre-selection", channel=self.Channel)
        #self.Stack.drawStack('llnunu_mt', self.preCuts, str(self.lumi*1000), 150, 0.0, 3000.0, titlex = "M_{T}^{ll#nu#nu}", units = "GeV", output=tag+'mt_high',outDir=self.outdir,separateSignal=self.sepSig,drawtex=self.whichregion+" pre-selection", channel=self.Channel)
    
        # merge all output plots into one pdf file
        os.system('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+self.outdir+'/'+tag+'all.pdf '+self.outdir+'/'+tag+'*.eps')
        print 'All plots merged in single pdf file '+tag+'.pdf .'
        # merge root file
        os.system('hadd -f '+self.outdir+'/'+tag+'all.root '+self.outdir+'/'+tag+'*.root')
        return

    def draw_BCD(self):
        #cuts=python.SetCuts.Cuts(self.Channel)

        for key in self.tex_dic:
            if ROOT.TString(key).Contains("SR"):
                continue
            else:
                tag = key+'_'+self.whichregion+'_'+self.Channel+'_'
                if self.Stack.log: tag=tag+'log_'
                else: pass
                self.Stack.drawStack('llnunu_mt', self.cuts[key], str(self.lumi*1000), 60, 0.0, 1200.0, titlex = "M_{T}", units = "GeV",
                                     output=tag+'mt',outDir=self.outdir,separateSignal=self.sepSig,
                                     drawtex=self.tex_dic[key],channel=self.Channel)
                self.Stack.drawStack('llnunu_mt', self.cuts[key], str(self.lumi*1000), 150, 0.0, 3000.0, titlex = "M_{T}", units = "GeV",
                                     output=tag+'mt_high',outDir=self.outdir,separateSignal=self.sepSig,
                                     drawtex=self.tex_dic[key],channel=self.Channel)
                self.Stack.drawStack('llnunu_l1_pt', self.cuts[key], str(self.lumi*1000), 100, 0, 1000, titlex = "p_{T}^{ll}", units = "GeV",
                                     output=tag+'pt',outDir=self.outdir,separateSignal=self.sepSig,
                                     drawtex=self.tex_dic[key],channel=self.Channel)
                self.Stack.drawStack('met_pt', self.cuts[key], str(self.lumi*1000), self.nbins, self.xMin, self.xMax, titlex = "E_{T}^{miss}", units = "GeV",
                                     output=tag+'met_low',outDir=self.outdir,separateSignal=self.sepSig,
                                     drawtex=self.tex_dic[key],channel=self.Channel)
                # merge all output plots into one pdf file
                os.system('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+self.outdir+'/'+tag+'all.pdf '+self.outdir+'/'+tag+'*.eps')
                print 'All plots merged in single pdf file '+tag+'.pdf .'
                # merge root file
                os.system('hadd -f '+self.outdir+'/'+tag+'all.root '+self.outdir+'/'+tag+'*.root')
        return

    def draw_A(self, useMC = False, sf_bd = 0.0):
        tag ='regionA_'+self.whichregion+'_'+self.Channel+'_'

        if not useMC:
            if sf_bd == 0.0: sf_bd=ROOT.Double(raw_input("Please give the scale factor derived from B,D regions: \n"))
            else: pass
    
            self.Stack.rmPlotter(self.ZJets, "ZJets","Z+Jets", "background")
            self.Stack.drawStack('met_pt', self.cuts['SR'], str(self.lumi*1000), self.nbins, self.xMin, self.xMax,
                                  titlex = "E_{T}^{miss}", units = "GeV",output=tag+'met_low',outDir=self.outdir,separateSignal=self.sepSig,
                                  drawtex=self.tex_dic['SR'],channel=self.Channel)
            
            file=ROOT.TFile(self.outdir+'/'+tag+'met_low.root')
    
            h_met_zjets_a=self.ZJets.drawTH1('met_pt',self.cuts['SR'],str(self.lumi*1000),self.nbins,self.xMin,self.xMax,titlex='E_{T}^{miss}',units='GeV',drawStyle="HIST")
            h_met_zjets_a.Scale(1./h_met_zjets_a.Integral())
            h_met_zjets_c=self.ZJets.drawTH1('met_pt',self.cuts['CRc'],str(self.lumi*1000),self.nbins,self.xMin,self.xMax,titlex='E_{T}^{miss}',units='GeV',drawStyle="HIST")
            h_met_zjets_c.Scale(1./h_met_zjets_c.Integral())
            h_met_zjets_a.Divide(h_met_zjets_c)
    
            h_met_zjets=self.Data.drawTH1('met_pt',self.cuts['CRc'],str('1'),self.nbins,self.xMin,self.xMax,titlex='E_{T}^{miss}',units='GeV',drawStyle="HIST")
            h_met_nonzjets=self.NonZBG.drawTH1('met_pt',self.cuts['CRc'],str(self.lumi*1000),self.nbins,self.xMin,self.xMax,titlex='E_{T}^{miss}',units='GeV',drawStyle="HIST")
            h_met_zjets.Add(h_met_nonzjets,-1) # subtract the non-z bkg MC from data-C
            h_met_zjets.Multiply(h_met_zjets_a)
            h_met_zjets.Scale(sf_bd)
            
            hframe=file.Get(tag+'met_low_frame')
            hs=file.Get(tag+'met_low_stack')
            hs.Add(h_met_zjets)
            hdata=file.Get(tag+'met_low_data')
            hratio=GetRatio_TH1(hdata,hs,True)
        
            legend=file.Get(tag+'met_low_legend')
            myentry=ROOT.TLegendEntry(h_met_zjets,"Z+jets(data-driven)","f")
            legend.GetListOfPrimitives().AddFirst(myentry)
    
            # Let's remove the signal entries in the legend 
            for i in legend.GetListOfPrimitives():
                if ROOT.TString(i.GetLabel()).Contains("Bulk"):
                    legend.GetListOfPrimitives().Remove(i)
                    
                    
            drawStack_simple(hframe, hs, hdata, hratio, legend,
                             hstack_opt="A, HIST",
                             outDir=self.outdir, output=self.whichregion+'_regionA_'+self.Channel+'_met_low',channel=ROOT.TString(self.Channel),
                             xmin=self.xMin, xmax=self.xMax, xtitle="E_{T}^{miss}" ,units="GeV",
                             lumi=self.lumi, notes="Region A ("+self.whichregion+")")
    
    
            self.Stack.addPlotter(self.ZJets, "ZJets","Z+Jets", "background")

        else:
            self.Stack.drawStack('met_pt', self.cuts['SR'], str(self.lumi*1000), self.nbins, self.xMin, self.xMax,
                                 titlex = "E_{T}^{miss}", units = "GeV",output=tag+'met_low',outDir=self.outdir,separateSignal=self.sepSig,
                                 drawtex=self.tex_dic['SR'],channel=self.Channel)
            
        return
        

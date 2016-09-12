#!/usr/bin/env python

import ROOT
import os
from TreePlotter import TreePlotter
from MergedPlotter import MergedPlotter
from StackPlotter import StackPlotter

#cutChain='loosecut'
cutChain='tightzpt100'
#cutChain='tightzpt100met100'
#cutChain='tightzpt100met100dphi'
#cutChain='tightmet250dphi'
#cutChain='zjetscut'
#cutChain='zjetscutmet50'

channel='both' # can be el or mu or both

outdir='plots'
indir="76X_JEC_Skim_BDT"
lumi=2.169126704526
sepSig=True
LogY=True
DrawLeptons=False
doRatio=True
test=True
Blind=True
FakeData=False

elChannel='(abs(llnunu_l1_l1_pdgId)==11)'
muChannel='(abs(llnunu_l1_l1_pdgId)==13)'

if not os.path.exists(outdir): os.system('mkdir '+outdir)

tag = cutChain+'_'
tag = tag+channel+'_'
if LogY: tag = tag+'log_'


cuts_loose='(nllnunu)'
cuts_loose_zll='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>100&&llnunu_l1_mass>70&&llnunu_l1_mass<110))'
cuts_loose_zll_met100='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>100&&llnunu_l1_mass>70&&llnunu_l1_mass<110&&met_pt>100))'
cuts_loose_zll_met100_dphi='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>100&&llnunu_l1_mass>70&&llnunu_l1_mass<110&&met_pt>100&&cos(llnunu_deltaPhi)<-0.2))'
cuts_loose_zll_met250_dphi='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>100&&llnunu_l1_mass>70&&llnunu_l1_mass<110&&met_pt>250&&cos(llnunu_deltaPhi)<-0.2))'

cuts_zjets='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(abs(llnunu_l1_mass-91.18)<20))'
cuts_zjets_met50='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(abs(llnunu_l1_mass-91.18)<20&&met_pt>50))'

if cutChain=='loosecut': cuts=cuts_loose
elif cutChain=='tightzpt100': cuts=cuts_loose_zll
elif cutChain=='tightzpt100met100': cuts=cuts_loose_zll_met100
elif cutChain=='tightzpt100met100dphi': cuts=cuts_loose_zll_met100_dphi
elif cutChain=='tightmet250dphi': cuts=cuts_loose_zll_met250_dphi
elif cutChain=='zjetscut': cuts=cuts_zjets
elif cutChain=='zjetscutmet50': cuts=cuts_zjets_met50
else : cuts=cuts_loose


if channel=='el': cuts = cuts+'&&'+elChannel
elif channel=='mu': cuts = cuts+'&&'+muChannel


ROOT.gROOT.ProcessLine('.x tdrstyle.C') 

NonResoPlotters=[]
NonResoSamples = ['NonReso']

for sample in NonResoSamples:
    NonResoPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    NonResoPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
    NonResoPlotters[-1].addCorrectionFactor('xsec','tree')
    NonResoPlotters[-1].addCorrectionFactor('genWeight','tree')
    NonResoPlotters[-1].addCorrectionFactor('puWeight','tree')

NonReso = MergedPlotter(NonResoPlotters)
NonReso.setFillProperties(1001,ROOT.kOrange)


ZResoPlotters=[]
ZResoSamples = ['ZReso']

for sample in ZResoSamples:
    ZResoPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    ZResoPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
    ZResoPlotters[-1].addCorrectionFactor('xsec','tree')
    ZResoPlotters[-1].addCorrectionFactor('genWeight','tree')
    ZResoPlotters[-1].addCorrectionFactor('puWeight','tree')

ZReso = MergedPlotter(ZResoPlotters)
ZReso.setFillProperties(1001,ROOT.kMagenta)

zjetsPlotters=[]
zjetsSamples = ['DYJetsToLL_M50_BIG'] 


for sample in zjetsSamples:
    zjetsPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    zjetsPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
    zjetsPlotters[-1].addCorrectionFactor('xsec','tree')
    zjetsPlotters[-1].addCorrectionFactor('genWeight','tree')
    zjetsPlotters[-1].addCorrectionFactor('puWeight','tree')

ZJets = MergedPlotter(zjetsPlotters)
ZJets.setFillProperties(1001,ROOT.kGreen+2)

sigPlotters=[]
sigSamples = [
'BulkGravToZZToZlepZinv_narrow_600',
#'BulkGravToZZToZlepZinv_narrow_800',
'BulkGravToZZToZlepZinv_narrow_1000',
#'BulkGravToZZToZlepZinv_narrow_1200',
#'BulkGravToZZToZlepZinv_narrow_1400',
#'BulkGravToZZToZlepZinv_narrow_1600', 
#'BulkGravToZZToZlepZinv_narrow_1800', 
'BulkGravToZZToZlepZinv_narrow_2000',
#'BulkGravToZZToZlepZinv_narrow_2500',
'BulkGravToZZToZlepZinv_narrow_3000',
#'BulkGravToZZToZlepZinv_narrow_3500', 
'BulkGravToZZToZlepZinv_narrow_4000', 
#'BulkGravToZZToZlepZinv_narrow_4500', 
]

k=1000
sigSampleNames = {
'BulkGravToZZToZlepZinv_narrow_600':str(k)+' x BulkG-600',
'BulkGravToZZToZlepZinv_narrow_800':str(k)+' x BulkG-800',
'BulkGravToZZToZlepZinv_narrow_1000':str(k)+' x BulkG-1000',
'BulkGravToZZToZlepZinv_narrow_1200':str(k)+' x BulkG-1200',
'BulkGravToZZToZlepZinv_narrow_1400':str(k)+' x BulkG-1400',
'BulkGravToZZToZlepZinv_narrow_1600':str(k)+' x BulkG-1600',
'BulkGravToZZToZlepZinv_narrow_1800':str(k)+' x BulkG-1800',
'BulkGravToZZToZlepZinv_narrow_2000':str(k)+' x BulkG-2000',
'BulkGravToZZToZlepZinv_narrow_2500':str(k)+' x BulkG-2500',
'BulkGravToZZToZlepZinv_narrow_3000':str(k)+' x BulkG-3000',
'BulkGravToZZToZlepZinv_narrow_3500':str(k)+' x BulkG-3500',
'BulkGravToZZToZlepZinv_narrow_4000':str(k)+' x BulkG-4000',
'BulkGravToZZToZlepZinv_narrow_4500':str(k)+' x BulkG-4500',
}
sigXsec = {
'BulkGravToZZToZlepZinv_narrow_600'  : 3.44631e-04*k,
'BulkGravToZZToZlepZinv_narrow_800'  : 6.31859e-05*k,
'BulkGravToZZToZlepZinv_narrow_1000' : 1.68661e-05*k,
'BulkGravToZZToZlepZinv_narrow_1200' : 5.59677e-06*k,
'BulkGravToZZToZlepZinv_narrow_1400' : 2.13168e-06*k,
'BulkGravToZZToZlepZinv_narrow_1600' : 8.97713e-07*k,
'BulkGravToZZToZlepZinv_narrow_1800' : 4.06090e-07*k,
'BulkGravToZZToZlepZinv_narrow_2000' : 1.94415e-07*k,
'BulkGravToZZToZlepZinv_narrow_2500' : 3.63496e-08*k,
'BulkGravToZZToZlepZinv_narrow_3000' : 7.95422e-09*k,
'BulkGravToZZToZlepZinv_narrow_3500' : 1.95002e-09*k,
'BulkGravToZZToZlepZinv_narrow_4000' : 5.03747e-10*k,
'BulkGravToZZToZlepZinv_narrow_4500' : 1*k,
}


for sample in sigSamples:
    sigPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    sigPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
    sigPlotters[-1].addCorrectionFactor(str(sigXsec[sample]),'tree')
    sigPlotters[-1].addCorrectionFactor('genWeight','tree')
    sigPlotters[-1].addCorrectionFactor('puWeight','tree')
    sigPlotters[-1].setFillProperties(0,ROOT.kWhite)


dataPlotters=[]
dataSamples = ['SingleElectron_Run2015C_25ns_16Dec',
'SingleElectron_Run2015D_16Dec',
'SingleMuon_Run2015C_25ns_16Dec',
'SingleMuon_Run2015D_16Dec']
for sample in dataSamples:
    dataPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))

Data = MergedPlotter(dataPlotters)




Stack = StackPlotter()
Stack.addPlotter(Data, "data_obs", "Data", "data")
Stack.addPlotter(NonReso, "NonReso","TT, WW, WZ non-reson.", "background")
Stack.addPlotter(ZReso, "ZReso","ZZ, WZ reson.", "background")
Stack.addPlotter(ZJets, "ZJets","Z+Jets", "background")


for i in range(len(sigSamples)):
  sigPlotters[i].setLineProperties(2,ROOT.kRed+i,2)
  Stack.addPlotter(sigPlotters[i],sigSamples[i],sigSampleNames[sigSamples[i]],'signal')  

 
Stack.setLog(LogY)
Stack.doRatio(doRatio)


Stack.drawStack('llnunu_mt', cuts, str(lumi*1000), 60, 0.0, 1200.0, titlex = "M_{T}", units = "GeV",output=tag+'mt',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=300)
#Stack.drawStack('llnunu_mt', cuts, str(lumi*1000), 150, 0.0, 3000.0, titlex = "M_{T}", units = "GeV",output=tag+'mt_high',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=300)
Stack.drawStack('met_pt', cuts, str(lumi*1000), 100, 0, 1000, titlex = "MET", units = "GeV",output=tag+'met',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=100)
#Stack.drawStack('met_pt', cuts, str(lumi*1000), 50, 0, 500, titlex = "MET", units = "GeV",output=tag+'met_low',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=100)
#Stack.drawStack('met_pt', cuts, str(lumi*1000), 150, 0, 1500, titlex = "MET", units = "GeV",output=tag+'met_high',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=100)
Stack.drawStack('bdt_m600', cuts, str(lumi*1000), 100, -1.0, 1.0, titlex = "BDT_m600", units = "",output=tag+'bdt_m600',outDir=outdir,separateSignal=sepSig)
Stack.drawStack('bdt_m800', cuts, str(lumi*1000), 100, -1.0, 1.0, titlex = "BDT_m800", units = "",output=tag+'bdt_m800',outDir=outdir,separateSignal=sepSig)
Stack.drawStack('bdt_m1000', cuts, str(lumi*1000), 100, -1.0, 1.0, titlex = "BDT_m1000", units = "",output=tag+'bdt_m1000',outDir=outdir,separateSignal=sepSig)
Stack.drawStack('bdt_m1200', cuts, str(lumi*1000), 100, -1.0, 1.0, titlex = "BDT_m1200", units = "",output=tag+'bdt_m1200',outDir=outdir,separateSignal=sepSig)
Stack.drawStack('bdt_m1400', cuts, str(lumi*1000), 100, -1.0, 1.0, titlex = "BDT_m1400", units = "",output=tag+'bdt_m1400',outDir=outdir,separateSignal=sepSig)
Stack.drawStack('bdt_m1600', cuts, str(lumi*1000), 100, -1.0, 1.0, titlex = "BDT_m1600", units = "",output=tag+'bdt_m1600',outDir=outdir,separateSignal=sepSig)
Stack.drawStack('bdt_m1800', cuts, str(lumi*1000), 100, -1.0, 1.0, titlex = "BDT_m1800", units = "",output=tag+'bdt_m1800',outDir=outdir,separateSignal=sepSig)
Stack.drawStack('bdt_m2000', cuts, str(lumi*1000), 100, -1.0, 1.0, titlex = "BDT_m2000", units = "",output=tag+'bdt_m2000',outDir=outdir,separateSignal=sepSig)
Stack.drawStack('bdt_m2500', cuts, str(lumi*1000), 100, -1.0, 1.0, titlex = "BDT_m2500", units = "",output=tag+'bdt_m2500',outDir=outdir,separateSignal=sepSig)
Stack.drawStack('bdt_m3000', cuts, str(lumi*1000), 100, -1.0, 1.0, titlex = "BDT_m3000", units = "",output=tag+'bdt_m3000',outDir=outdir,separateSignal=sepSig)
Stack.drawStack('bdt_m3500', cuts, str(lumi*1000), 100, -1.0, 1.0, titlex = "BDT_m3500", units = "",output=tag+'bdt_m3500',outDir=outdir,separateSignal=sepSig)
Stack.drawStack('bdt_m4000', cuts, str(lumi*1000), 100, -1.0, 1.0, titlex = "BDT_m4000", units = "",output=tag+'bdt_m4000',outDir=outdir,separateSignal=sepSig)



if not test :
    Stack.drawStack('llnunu_l1_mass', cuts, str(lumi*1000), 50, 50, 150, titlex = "M(Z)", units = "GeV",output=tag+'zmass',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_pt', cuts, str(lumi*1000), 75, 0.0, 1500.0, titlex = "P_{T}(Z)", units = "GeV",output=tag+'zpt',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('met_pt/sqrt(llnunu_l2_sumEt)', cuts, str(lumi*1000), 100, 0.0, 20.0, titlex = "MET/#sqrt{sumE_{T}}", units = "",output=tag+'metOvSqSET',outDir=outdir,separateSignal=sepSig)

    Stack.drawStack('llnunu_deltaPhi', cuts, str(lumi*1000), 50, 0, 5, titlex = "#Delta#phi(Z,MET)", units = "",output=tag+'dphiZMet',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('cos(llnunu_deltaPhi)', cuts, str(lumi*1000), 50, -1, 1, titlex = "cos(#Delta#phi(Z,MET))", units = "",output=tag+'CosdphiZMet',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('abs(llnunu_l1_pt+met_pt*cos(llnunu_deltaPhi))', cuts, str(lumi*1000), 100, 0, 1000, titlex = "#Delta P_{T}^{#parallel}(Z,MET)", units = "GeV",output=tag+'dPTPara',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('abs(met_pt*sin(llnunu_deltaPhi))', cuts, str(lumi*1000), 100, 0, 100, titlex = "#Delta P_{T}^{#perp}(Z,MET)", units = "GeV",output=tag+'dPTPerp',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('abs(llnunu_l1_pt+met_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts, str(lumi*1000), 100, 0, 5, titlex = "#Delta P_{T}^{#parallel}(Z,MET)/P_{T}(Z)", units = "",output=tag+'dPTParaRel',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('abs(met_pt*sin(llnunu_deltaPhi))/llnunu_l1_pt', cuts, str(lumi*1000), 100, 0, 5, titlex = "#Delta P_{T}^{#perp}(Z,MET)/P_{T}", units = "",output=tag+'dPTPerpRel',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('nVert', cuts, str(lumi*1000), 50, 0.0, 50.0, titlex = "N vertices", units = "",output=tag+'nVert',outDir=outdir,separateSignal=sepSig)
 
    Stack.drawStack('llnunu_l1_eta', cuts, str(lumi*1000), 100, -5.0, 5.0, titlex = "#eta(Z) ", units = "",output=tag+'zeta',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_phi', cuts, str(lumi*1000), 64, -3.2, 3.2, titlex = "#phi(Z)", units = "",output=tag+'zphi',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_mt', cuts, str(lumi*1000), 50, 0.0, 150.0, titlex = "M_{T}(Z)", units = "GeV",output=tag+'zmt',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_deltaPhi', cuts, str(lumi*1000), 50, 0.0, 5.0, titlex = "#Delta#phi", units = "",output=tag+'ZdeltaPhi',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('cos(llnunu_l1_deltaPhi)', cuts, str(lumi*1000), 50, -1, 1, titlex = "Cos(#Delta#phi)", units = "",output=tag+'CosZdeltaPhi',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_deltaR', cuts, str(lumi*1000), 50, 0.0, 5.0, titlex = "#Delta R", units = "",output=tag+'ZdeltaR',outDir=outdir,separateSignal=sepSig)
  
    Stack.drawStack('llnunu_l2_phi', cuts, str(lumi*1000), 100, -3.2, 3.2, titlex = "#phi(MET)", units = "",output=tag+'metPhi',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l2_sumEt', cuts, str(lumi*1000), 80, 0.0, 2400.0, titlex = "sumE_{T}", units = "GeV",output=tag+'metSumEt',outDir=outdir,separateSignal=sepSig)

if DrawLeptons:
    Stack.drawStack('llnunu_l1_l1_pt', cuts, str(lumi*1000), 100, 0.0, 1000.0, titlex = "P_{T}(l_{1})", units = "GeV",output=tag+'pTlep1',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l1_eta', cuts, str(lumi*1000), 60, -3.0, 3.0, titlex = "#eta(l_{1})", units = "",output=tag+'etalep1',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l1_phi', cuts, str(lumi*1000), 64, -3.2, 3.2, titlex = "#phi(l_{1})", units = "",output=tag+'philep1',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l1_miniRelIso', cuts, str(lumi*1000), 100, 0.0, 0.2, titlex = "miniISO_{rel}(l_{1})", units = "",output=tag+'miniISOlep1',outDir=outdir,separateSignal=sepSig)

    Stack.drawStack('llnunu_l1_l2_pt', cuts, str(lumi*1000), 100, 0.0, 1000.0, titlex = "P_{T}(l_{2})", units = "GeV",output=tag+'pTlep2',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l2_eta', cuts, str(lumi*1000), 60, -3.0, 3.0, titlex = "#eta(l_{2})", units = "",output=tag+'etalep2',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l2_phi', cuts, str(lumi*1000), 64, -3.2, 3.2, titlex = "#phi(l_{2})", units = "",output=tag+'philep2',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l2_miniRelIso', cuts, str(lumi*1000), 100, 0.0, 0.2, titlex = "miniISO_{rel}(l_{2})", units = "",output=tag+'miniISOlep2',outDir=outdir,separateSignal=sepSig)

    Stack.drawStack('llnunu_l1_l1_pt', cuts+"&&(abs(llnunu_l1_l1_pdgId)==13)", str(lumi*1000), 100, 0.0, 1000.0, titlex = "P_{T}(#mu_{1})", units = "GeV",output=tag+'pTlep1_mu',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l1_eta', cuts+"&&(abs(llnunu_l1_l1_pdgId)==13)", str(lumi*1000), 60, -3.0, 3.0, titlex = "#eta(#mu_{1})", units = "",output=tag+'etalep1_mu',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l1_phi', cuts+"&&(abs(llnunu_l1_l1_pdgId)==13)", str(lumi*1000), 64, -3.2, 3.2, titlex = "#phi(#mu_{1})", units = "",output=tag+'philep1_mu',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l1_miniRelIso', cuts+"&&(abs(llnunu_l1_l1_pdgId)==13)", str(lumi*1000), 100, 0.0, 0.2, titlex = "miniISO_{rel}(#mu_{1})", units = "",output=tag+'miniISOlep1_mu',outDir=outdir,separateSignal=sepSig)

    Stack.drawStack('llnunu_l1_l2_pt', cuts+"&&(abs(llnunu_l1_l1_pdgId)==13)", str(lumi*1000), 100, 0.0, 1000.0, titlex = "P_{T}(#mu_{2})", units = "GeV",output=tag+'pTlep2_mu',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l2_eta', cuts+"&&(abs(llnunu_l1_l1_pdgId)==13)", str(lumi*1000), 60, -3.0, 3.0, titlex = "#eta(#mu_{2})", units = "",output=tag+'etalep2_mu',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l2_phi', cuts+"&&(abs(llnunu_l1_l1_pdgId)==13)", str(lumi*1000), 64, -3.2, 3.2, titlex = "#phi(#mu_{2})", units = "",output=tag+'philep2_mu',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l2_miniRelIso', cuts+"&&(abs(llnunu_l1_l1_pdgId)==13)", str(lumi*1000), 100, 0.0, 0.2, titlex = "miniISO_{rel}(#mu_{2})", units = "",output=tag+'miniISOlep2_mu',outDir=outdir,separateSignal=sepSig)

    Stack.drawStack('llnunu_l1_l1_pt', cuts+"&&(abs(llnunu_l1_l1_pdgId)==11)", str(lumi*1000), 100, 0.0, 1000.0, titlex = "P_{T}(e_{1})", units = "GeV",output=tag+'pTlep1_el',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l1_eta', cuts+"&&(abs(llnunu_l1_l1_pdgId)==11)", str(lumi*1000), 60, -3.0, 3.0, titlex = "#eta(e_{1})", units = "",output=tag+'etalep1_el',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l1_phi', cuts+"&&(abs(llnunu_l1_l1_pdgId)==11)", str(lumi*1000), 64, -3.2, 3.2, titlex = "#phi(e_{1})", units = "",output=tag+'philep1_el',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l1_miniRelIso', cuts+"&&(abs(llnunu_l1_l1_pdgId)==11)", str(lumi*1000), 100, 0.0, 0.2, titlex = "miniISO_{rel}(e_{1})", units = "",output=tag+'miniISOlep1_el',outDir=outdir,separateSignal=sepSig)

    Stack.drawStack('llnunu_l1_l2_pt', cuts+"&&(abs(llnunu_l1_l1_pdgId)==11)", str(lumi*1000), 100, 0.0, 1000.0, titlex = "P_{T}(e_{2})", units = "GeV",output=tag+'pTlep2_el',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l2_eta', cuts+"&&(abs(llnunu_l1_l1_pdgId)==11)", str(lumi*1000), 60, -3.0, 3.0, titlex = "#eta(e_{2})", units = "",output=tag+'etalep2_el',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l2_phi', cuts+"&&(abs(llnunu_l1_l1_pdgId)==11)", str(lumi*1000), 64, -3.2, 3.2, titlex = "#phi(e_{2})", units = "",output=tag+'philep2_el',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l2_miniRelIso', cuts+"&&(abs(llnunu_l1_l1_pdgId)==11)", str(lumi*1000), 100, 0.0, 0.2, titlex = "miniISO_{rel}(e_{2})", units = "",output=tag+'miniISOlep2_el',outDir=outdir,separateSignal=sepSig)

    Stack.drawStack('llnunu_l1_l1_pt', cuts+"&&(abs(llnunu_l1_l1_pdgId)==13)", str(lumi*1000), 100, 45.0, 1000.0, titlex = "P_{T}(#mu_{1})", units = "GeV",output=tag+'pTlep1_mu_pthlt',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_l1_pt', cuts+"&&(abs(llnunu_l1_l1_pdgId)==11)", str(lumi*1000), 100, 105.0, 1000.0, titlex = "P_{T}(e_{1})", units = "GeV",output=tag+'pTlep1_el_pthlt',outDir=outdir,separateSignal=sepSig)


# merge all output plots into one pdf file
os.system('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+outdir+'/'+tag+'all.pdf '+outdir+'/'+tag+'*.eps')
print 'All plots merged in single pdf file '+tag+'.pdf .'
# merge root file
os.system('rm '+outdir+'/'+tag+'all.root ')
os.system('hadd -f '+outdir+'/'+tag+'all.root '+outdir+'/'+tag+'*.root')



#!/usr/bin/env python

import ROOT
import os, string, math, pickle
from CMGTools.XZZ2l2nu.plotting.TreePlotter import TreePlotter
from CMGTools.XZZ2l2nu.plotting.MergedPlotter import MergedPlotter
from CMGTools.XZZ2l2nu.plotting.StackPlotter import StackPlotter

tag='New80X_'

#cutChain='loosecut'
cutChain='tight'
#cutChain='tightzpt100'
#cutChain='tightzpt100met50'
#cutChain='tightzpt100met100'
#cutChain='tightzpt100met200'
#cutChain='tightzpt100met100dphi'
#cutChain='tightmet250dphi'
#cutChain='zjetscut'
#cutChain='zjetscutmet50'

useZPtWeight=True
DataPreHLT=False

channel='both' # can be el or mu or both
# samples tag
#stag='bkg'
#stag='zjets'
#stag='sig'
stag='all'


outdir='plots'
#indir='/data/XZZ/76X_Ntuple/76X_20160514_PUSkim2'
indir='/home/heli/work/XZZ/80X_Ntuple/80X_20160606_NoHLT_Skim'
#indir="/afs/cern.ch/work/m/mewu/public/76X_new"
#indir="/data/mewu/76X_new"
#lumi=2.169126704526
lumi=0.5893
sepSig=True
LogY=True
DrawLeptons=False
doRatio=True
test=True
Blind=False
FakeData=False
UseMETFilter=True
SignalAll1pb=False
puWeight='puWeight'
k=1 # signal scale

elChannel='(abs(llnunu_l1_l1_pdgId)==11)'
muChannel='(abs(llnunu_l1_l1_pdgId)==13)'

if not os.path.exists(outdir): os.system('mkdir '+outdir)

tag = tag+cutChain+'_'
tag = tag+puWeight+'_'


if UseMETFilter: tag = tag+'metfilter_'

if not Blind: tag = tag+'unblind_'

tag = tag+channel+'_'
if LogY: tag = tag+'log_'
if SignalAll1pb:
    tag += '1pb_'
else:
    tag += 'scale'+str(k)+'_'

metfilter='(Flag_EcalDeadCellTriggerPrimitiveFilter&&Flag_HBHENoiseIsoFilter&&Flag_goodVertices&&Flag_HBHENoiseFilter&&Flag_CSCTightHalo2015Filter&&Flag_eeBadScFilter)'

cuts_loose='(nllnunu)'
cuts_loose_z='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_mass>70&&llnunu_l1_mass<110))'
cuts_loose_zll='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>100&&llnunu_l1_mass>70&&llnunu_l1_mass<110))'
cuts_loose_zll_met50='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>100&&llnunu_l1_mass>70&&llnunu_l1_mass<110&&llnunu_l2_pt>50))'
cuts_loose_zll_met100='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>100&&llnunu_l1_mass>70&&llnunu_l1_mass<110&&llnunu_l2_pt>100))'
cuts_loose_zll_met200='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>100&&llnunu_l1_mass>70&&llnunu_l1_mass<110&&llnunu_l2_pt>200))'
cuts_loose_zll_met100_dphi='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>100&&llnunu_l1_mass>70&&llnunu_l1_mass<110&&llnunu_l2_pt>100&&cos(llnunu_deltaPhi)<-0.2))'
cuts_loose_zll_met250_dphi='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_pt>100&&llnunu_l1_mass>70&&llnunu_l1_mass<110&&llnunu_l2_pt>250&&cos(llnunu_deltaPhi)<-0.2))'

cuts_zjets='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(abs(llnunu_l1_mass-91.18)<20))'
cuts_zjets_met50='(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(abs(llnunu_l1_mass-91.18)<20&&llnunu_l2_pt>50))'

if cutChain=='loosecut': cuts=cuts_loose
elif cutChain=='tight': cuts=cuts_loose_z
elif cutChain=='tightzpt100': cuts=cuts_loose_zll
elif cutChain=='tightzpt100met50': cuts=cuts_loose_zll_met50
elif cutChain=='tightzpt100met100': cuts=cuts_loose_zll_met100
elif cutChain=='tightzpt100met200': cuts=cuts_loose_zll_met200
elif cutChain=='tightzpt100met100dphi': cuts=cuts_loose_zll_met100_dphi
elif cutChain=='tightmet250dphi': cuts=cuts_loose_zll_met250_dphi
elif cutChain=='zjetscut': cuts=cuts_zjets
elif cutChain=='zjetscutmet50': cuts=cuts_zjets_met50
else : cuts=cuts_loose


if channel=='el': cuts = cuts+'&&'+elChannel
elif channel=='mu': cuts = cuts+'&&'+muChannel

if UseMETFilter:
    cuts = '('+cuts+'&&'+metfilter+')'

ROOT.gROOT.ProcessLine('.x tdrstyle.C') 

#allSamples = ['WWTo2L2Nu','WWToLNuQQ','WZTo1L1Nu2Q',
#            'WZTo2L2Q','WZTo3LNu',
#'ZZTo2L2Nu',
#'ZZTo2L2Q','ZZTo4L',
#'DYJetsToLL_M50_BIG',
#'TTTo2L2Nu',
#'BulkGravToZZToZlepZinv_narrow_600',
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
#]

allPlotters = {}


wwPlotters=[]
wwSamples = ['WWTo2L2Nu','WWToLNuQQ','WZTo1L1Nu2Q']

for sample in wwSamples:
    wwPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    wwPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
    wwPlotters[-1].addCorrectionFactor('xsec','tree')
    wwPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    wwPlotters[-1].addCorrectionFactor(puWeight,'tree')
    wwPlotters[-1].addCorrectionFactor('(llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf)','tree')
    wwPlotters[-1].addCorrectionFactor('triggersf','tree')
    allPlotters[sample] = wwPlotters[-1]

WW = MergedPlotter(wwPlotters)
WW.setFillProperties(1001,ROOT.kOrange)


vvPlotters=[]
vvSamples = ['WZTo2L2Q','WZTo3LNu',
'ZZTo2L2Nu',
'ZZTo2L2Q','ZZTo4L']

for sample in vvSamples:
    vvPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    vvPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
    vvPlotters[-1].addCorrectionFactor('xsec','tree')
    vvPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    vvPlotters[-1].addCorrectionFactor(puWeight,'tree')
    vvPlotters[-1].addCorrectionFactor('(llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf)','tree')
    vvPlotters[-1].addCorrectionFactor('triggersf','tree')
    allPlotters[sample] = vvPlotters[-1]

VV = MergedPlotter(vvPlotters)
VV.setFillProperties(1001,ROOT.kMagenta)

wjetsPlotters=[]
wjetsSamples = ['WJetsToLNu']

for sample in wjetsSamples:
    wjetsPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    wjetsPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
    wjetsPlotters[-1].addCorrectionFactor('xsec','tree')
    wjetsPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    wjetsPlotters[-1].addCorrectionFactor(puWeight,'tree')
    wjetsPlotters[-1].addCorrectionFactor('(llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf)','tree')
    wjetsPlotters[-1].addCorrectionFactor('triggersf','tree')

WJets = MergedPlotter(wjetsPlotters)
WJets.setFillProperties(1001,ROOT.kBlue-6)


zjetsPlotters=[]
zjetsSamples = ['DYJetsToLL_M50_ZPt']
#zjetsSamples = ['DYJetsToLL_M50_MGMLM_Ext1']
#zjetsSamples = ['DYJetsToLL_M50']

#zjetsSamples = ['DYJetsToLL_M50_HT100to200','DYJetsToLL_M50_HT200to400','DYJetsToLL_M50_HT400to600','DYJetsToLL_M50_HT600toInf']
#zjetsSamples = ['DYJetsToLL_M50','DYJetsToLL_M50_Ext']
#zjetsSamples = ['DYJetsToLL_M50_BIG'] # M50_BIG = M50 + M50_Ext
#zjetsSamples = ['DYJetsToLL_M50'] # M50_BIG = M50 + M50_Ext


for sample in zjetsSamples:
    zjetsPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    zjetsPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
    if useZPtWeight:
        zjetsPlotters[-1].addCorrectionFactor('ZPtWeight','ZPtWeight')
    #zjetsPlotters[-1].addCorrectionFactor('PhiStarWeight','PhiStarWeight')
    #zjetsPlotters[-1].addCorrectionFactor('xsec','tree')
    zjetsPlotters[-1].addCorrectionFactor('(1921.8*3)','xsec') # FEWZ NNLO.results_z_m50_nnlo_inclusive_NNPDF30_nlo_as_0118
    #zjetsPlotters[-1].addCorrectionFactor('(1907.0*3)','xsec') # FEWZ NNLO.results_z_m50_nnlo_fsrOn_lowstat_inclusive_NNPDF30_nlo_as_0118
    #zjetsPlotters[-1].addCorrectionFactor('(1870.0*3)','xsec') # FEWZ NNLO.results_z_m50_nnlo_inclusive_NNPDF30_nlo_as_0118
    zjetsPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    zjetsPlotters[-1].addCorrectionFactor(puWeight,'puWeight')
    #zjetsPlotters[-1].addCorrectionFactor("puWeight655",'puWeight')
    #zjetsPlotters[-1].addCorrectionFactor("puWeight"+putag,'puWeight')
    #zjetsPlotters[-1].addCorrectionFactor('(llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf)','tree')
    #zjetsPlotters[-1].addCorrectionFactor('triggersf','tree')
    allPlotters[sample] = zjetsPlotters[-1]

ZJets = MergedPlotter(zjetsPlotters)
ZJets.setFillProperties(1001,ROOT.kGreen+2)

ttPlotters=[]
ttSamples = ['TTTo2L2Nu']#,'TTZToLLNuNu','TTWJetsToLNu']

for sample in ttSamples:
    ttPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    ttPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
    ttPlotters[-1].addCorrectionFactor('xsec','tree')
    ttPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    ttPlotters[-1].addCorrectionFactor(puWeight,'tree')
    ttPlotters[-1].addCorrectionFactor('(llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf)','tree')
    ttPlotters[-1].addCorrectionFactor('triggersf','tree')
    allPlotters[sample] = ttPlotters[-1]

TT = MergedPlotter(ttPlotters)
TT.setFillProperties(1001,ROOT.kAzure-9)

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
#'BulkGravToZZToZlepZinv_narrow_3000',
#'BulkGravToZZToZlepZinv_narrow_3500', 
#'BulkGravToZZToZlepZinv_narrow_4000', 
#'BulkGravToZZToZlepZinv_narrow_4500', 
]


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
4500:1.0,
}

sigXsec = {
'BulkGravToZZToZlepZinv_narrow_600'  : BulkGZZ2l2nuXsec[600]*k,
'BulkGravToZZToZlepZinv_narrow_800'  : BulkGZZ2l2nuXsec[800]*k,
'BulkGravToZZToZlepZinv_narrow_1000' : BulkGZZ2l2nuXsec[1000]*k,
'BulkGravToZZToZlepZinv_narrow_1200' : BulkGZZ2l2nuXsec[1200]*k,
'BulkGravToZZToZlepZinv_narrow_1400' : BulkGZZ2l2nuXsec[1400]*k,
'BulkGravToZZToZlepZinv_narrow_1600' : BulkGZZ2l2nuXsec[1600]*k,
'BulkGravToZZToZlepZinv_narrow_1800' : BulkGZZ2l2nuXsec[1800]*k,
'BulkGravToZZToZlepZinv_narrow_2000' : BulkGZZ2l2nuXsec[2000]*k,
'BulkGravToZZToZlepZinv_narrow_2500' : BulkGZZ2l2nuXsec[2500]*k,
'BulkGravToZZToZlepZinv_narrow_3000' : BulkGZZ2l2nuXsec[3000]*k,
'BulkGravToZZToZlepZinv_narrow_3500' : BulkGZZ2l2nuXsec[3500]*k,
'BulkGravToZZToZlepZinv_narrow_4000' : BulkGZZ2l2nuXsec[4000]*k,
'BulkGravToZZToZlepZinv_narrow_4500' : BulkGZZ2l2nuXsec[4500]*k,
}

if SignalAll1pb:
    for sig in sigSamples:
        sigXsec[sig] = 1.0
        sigSampleNames[sig] = string.replace(sigSampleNames[sig], str(k)+' x', '1 pb')

for sample in sigSamples:
    sigPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    sigPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
    sigPlotters[-1].addCorrectionFactor(str(sigXsec[sample]),'tree')
    sigPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    sigPlotters[-1].addCorrectionFactor(puWeight,'tree')
    sigPlotters[-1].addCorrectionFactor('(llnunu_l1_l1_lepsf*llnunu_l1_l2_lepsf)','tree')
    sigPlotters[-1].addCorrectionFactor('triggersf','tree')
    sigPlotters[-1].setFillProperties(0,ROOT.kWhite)
    allPlotters[sample] = sigPlotters[-1]


dataPlotters=[]
dataSamples = ['SingleMuon_Run2016B_PromptReco_v2', 'SingleElectron_Run2016B_PromptReco_v2']

if DataPreHLT:
    dataSamples = ['SingleMuon_Run2016B_PromptReco_v2_HLT_MU', 'SingleElectron_Run2016B_PromptReco_v2_HLT_ELE']

for sample in dataSamples:
    dataPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))

if not DataPreHLT:
    dataPlotters[0].addCorrectionFactor('(HLT_MU)','HLT')
    dataPlotters[1].addCorrectionFactor('(HLT_ELE&&!HLT_MU)','HLT')

Data = MergedPlotter(dataPlotters)


Stack = StackPlotter()
Stack.addPlotter(Data, "data_obs", "Data", "data")
#Stack.addPlotter(WJets, "WJets","W+Jets", "background")
Stack.addPlotter(WW, "VVNonReso","WW, WZ non-reson.", "background")
Stack.addPlotter(TT, "TT","TT", "background")
Stack.addPlotter(VV, "VVZReso","ZZ, WZ reson.", "background")
Stack.addPlotter(ZJets, "ZJets","Z+Jets", "background")


for i in range(len(sigSamples)):
  sigPlotters[i].setLineProperties(2,ROOT.kRed+i,2)
  Stack.addPlotter(sigPlotters[i],sigSamples[i],sigSampleNames[sigSamples[i]],'signal')  

 
Stack.setLog(LogY)
Stack.doRatio(doRatio)


###
fout = ROOT.TFile(tag+'.root', 'recreate')

# both
#h1 = Data.drawTH2('h_data_upara_zpt', 'llnunu_l1_pt:(llnunu_l1_pt+llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts, "(1)", 1500, -15, 15, 1000, 0, 1000,  titlex = "#Delta P_{T}^{#parallel}(Z,MET)/P_{T}(Z)",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")
#h2 = ZJets.drawTH2('h_zjets_upara_zpt', 'llnunu_l1_pt:(llnunu_l1_pt+llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts, str(lumi*1000), 1500, -15, 15, 1000, 0, 1000,  titlex = "#Delta P_{T}^{#parallel}(Z,MET)/P_{T}(Z)",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")
#h3 = Data.drawTH2('h_data_uperp_zpt_el', 'llnunu_l1_pt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts, "(1)", 1500, -15, 15, 1000, 0, 1000,  titlex = "#Delta P_{T}^{#perp}(Z,MET)/P_{T}",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")
#h4 = ZJets.drawTH2('h_zjets_uperp_zpt_el', 'llnunu_l1_pt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts, str(lumi*1000), 1500, -15, 15, 1000, 0, 1000,  titlex = "#Delta P_{T}^{#perp}(Z,MET)/P_{T}",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")


h1 = Data.drawTH2('h_data_upara_zpt', 'llnunu_l1_pt:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts, "(1)", 1500, -15, 15, 1000, 0, 1000,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")
h2 = ZJets.drawTH2('h_zjets_upara_zpt', 'llnunu_l1_pt:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts, str(lumi*1000), 1500, -15, 15, 1000, 0, 1000,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")


h3 = Data.drawTH2('h_data_uperp_zpt', 'llnunu_l1_pt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts, "(1)", 1500, -15, 15, 1000, 0, 1000,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")
h4 = ZJets.drawTH2('h_zjets_uperp_zpt', 'llnunu_l1_pt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts, str(lumi*1000), 1500, -15, 15, 1000, 0, 1000,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")

h1.Write()
h2.Write()
h3.Write()
h4.Write()

h1 = Data.drawTH2('h_data_upara_zpt_zeta', 'abs(llnunu_l1_eta):(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts, "(1)", 1500, -15, 15, 1000, 0, 10,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "#eta(Z)", unitsy = "",drawStyle = "COLZ")
h2 = ZJets.drawTH2('h_zjets_upara_zpt_zeta', 'abs(llnunu_l1_eta):(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts, str(lumi*1000), 1500, -15, 15, 1000, 0, 10,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "#eta(Z)", unitsy = "",drawStyle = "COLZ")


h3 = Data.drawTH2('h_data_uperp_zpt_zeta', 'abs(llnunu_l1_eta):(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts, "(1)", 1500, -15, 15, 1000, 0, 10,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "#eta(Z)", unitsy = "",drawStyle = "COLZ")
h4 = ZJets.drawTH2('h_zjets_uperp_zpt_zeta', 'abs(llnunu_l1_eta):(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts, str(lumi*1000), 1500, -15, 15, 1000, 0, 10,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "#eta(Z)", unitsy = "",drawStyle = "COLZ")

h1.Write()
h2.Write()
h3.Write()
h4.Write()

h1 = Data.drawTH2('h_data_upara_zpt_nvtx', 'nVert:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts, "(1)", 1500, -15, 15, 50, 0, 50,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "nVertex", unitsy = "",drawStyle = "COLZ")
h2 = ZJets.drawTH2('h_zjets_upara_zpt_nvtx', 'nVert:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts, str(lumi*1000), 1500, -15, 15, 50, 0, 50,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "N vertex", unitsy = "",drawStyle = "COLZ")


h3 = Data.drawTH2('h_data_uperp_zpt_nvtx', 'nVert:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts, "(1)", 1500, -15, 15, 50, 0, 50,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "N vertex", unitsy = "",drawStyle = "COLZ")
h4 = ZJets.drawTH2('h_zjets_uperp_zpt_nvtx', 'nVert:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts, str(lumi*1000), 1500, -15, 15, 50, 0, 50,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "N vertex", unitsy = "",drawStyle = "COLZ")

h1.Write()
h2.Write()
h3.Write()
h4.Write()

h1 = Data.drawTH2('h_data_upara_zpt_set', 'llnunu_l2_sumEt:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts, "(1)", 1500, -15, 15, 300, 0, 3000,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "sumEt", unitsy = "GeV",drawStyle = "COLZ")
h2 = ZJets.drawTH2('h_zjets_upara_zpt_set', 'llnunu_l2_sumEt:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts, str(lumi*1000), 1500, -15, 15, 300, 0, 3000,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "sumEt", unitsy = "GeV",drawStyle = "COLZ")


h3 = Data.drawTH2('h_data_uperp_zpt_set', 'llnunu_l2_sumEt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts, "(1)", 1500, -15, 15, 300, 0, 3000,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "sumEt", unitsy = "GeV",drawStyle = "COLZ")
h4 = ZJets.drawTH2('h_zjets_uperp_zpt_set', 'llnunu_l2_sumEt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts, str(lumi*1000), 1500, -15, 15, 300, 0, 3000,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "sumEt", unitsy = "GeV",drawStyle = "COLZ")

h1.Write()
h2.Write()
h3.Write()
h4.Write()


# electron
#h1 = Data.drawTH2('h_data_upara_zpt_el', 'llnunu_l1_pt:(llnunu_l1_pt+llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+elChannel, "(1)", 1500, -15, 15, 1000, 0, 1000,  titlex = "#Delta P_{T}^{#parallel}(Z,MET)/P_{T}(Z)",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")
#h2 = ZJets.drawTH2('h_zjets_upara_zpt_el', 'llnunu_l1_pt:(llnunu_l1_pt+llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+elChannel, str(lumi*1000), 1500, -15, 15, 1000, 0, 1000,  titlex = "#Delta P_{T}^{#parallel}(Z,MET)/P_{T}(Z)",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")


#h3 = Data.drawTH2('h_data_uperp_zpt_el', 'llnunu_l1_pt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+elChannel, "(1)", 1500, -15, 15, 1000, 0, 1000,  titlex = "#Delta P_{T}^{#perp}(Z,MET)/P_{T}",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")
#h4 = ZJets.drawTH2('h_zjets_uperp_zpt_el', 'llnunu_l1_pt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+elChannel, str(lumi*1000), 1500, -15, 15, 1000, 0, 1000,  titlex = "#Delta P_{T}^{#perp}(Z,MET)/P_{T}",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")

h1 = Data.drawTH2('h_data_upara_zpt_el', 'llnunu_l1_pt:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+elChannel, "(1)", 1500, -15, 15, 1000, 0, 1000,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")
h2 = ZJets.drawTH2('h_zjets_upara_zpt_el', 'llnunu_l1_pt:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+elChannel, str(lumi*1000), 1500, -15, 15, 1000, 0, 1000,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")


h3 = Data.drawTH2('h_data_uperp_zpt_el', 'llnunu_l1_pt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+elChannel, "(1)", 1500, -15, 15, 1000, 0, 1000,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")
h4 = ZJets.drawTH2('h_zjets_uperp_zpt_el', 'llnunu_l1_pt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+elChannel, str(lumi*1000), 1500, -15, 15, 1000, 0, 1000,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")

h1.Write()
h2.Write()
h3.Write()
h4.Write()

h1 = Data.drawTH2('h_data_upara_zpt_zeta_el', 'abs(llnunu_l1_eta):(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+elChannel, "(1)", 1500, -15, 15, 1000, 0, 10,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "#eta(Z)", unitsy = "",drawStyle = "COLZ")
h2 = ZJets.drawTH2('h_zjets_upara_zpt_zeta_el', 'abs(llnunu_l1_eta):(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+elChannel, str(lumi*1000), 1500, -15, 15, 1000, 0, 10,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "#eta(Z)", unitsy = "",drawStyle = "COLZ")


h3 = Data.drawTH2('h_data_uperp_zpt_zeta_el', 'abs(llnunu_l1_eta):(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+elChannel, "(1)", 1500, -15, 15, 1000, 0, 10,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "#eta(Z)", unitsy = "",drawStyle = "COLZ")
h4 = ZJets.drawTH2('h_zjets_uperp_zpt_zeta_el', 'abs(llnunu_l1_eta):(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+elChannel, str(lumi*1000), 1500, -15, 15, 1000, 0, 10,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "#eta(Z)", unitsy = "",drawStyle = "COLZ")

h1.Write()
h2.Write()
h3.Write()
h4.Write()


h1 = Data.drawTH2('h_data_upara_zpt_nvtx_el', 'nVert:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+elChannel, "(1)", 1500, -15, 15, 50, 0, 50,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "nVertex", unitsy = "",drawStyle = "COLZ")
h2 = ZJets.drawTH2('h_zjets_upara_zpt_nvtx_el', 'nVert:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+elChannel, str(lumi*1000), 1500, -15, 15, 50, 0, 50,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "N vertex", unitsy = "",drawStyle = "COLZ")


h3 = Data.drawTH2('h_data_uperp_zpt_nvtx_el', 'nVert:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+elChannel, "(1)", 1500, -15, 15, 50, 0, 50,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "N vertex", unitsy = "",drawStyle = "COLZ")
h4 = ZJets.drawTH2('h_zjets_uperp_zpt_nvtx_el', 'nVert:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+elChannel, str(lumi*1000), 1500, -15, 15, 50, 0, 50,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "N vertex", unitsy = "",drawStyle = "COLZ")

h1.Write()
h2.Write()
h3.Write()
h4.Write()

h1 = Data.drawTH2('h_data_upara_zpt_set_el', 'llnunu_l2_sumEt:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+elChannel, "(1)", 1500, -15, 15, 300, 0, 3000,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "sumEt", unitsy = "GeV",drawStyle = "COLZ")
h2 = ZJets.drawTH2('h_zjets_upara_zpt_set_el', 'llnunu_l2_sumEt:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+elChannel, str(lumi*1000), 1500, -15, 15, 300, 0, 3000,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "sumEt", unitsy = "GeV",drawStyle = "COLZ")


h3 = Data.drawTH2('h_data_uperp_zpt_set_el', 'llnunu_l2_sumEt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+elChannel, "(1)", 1500, -15, 15, 300, 0, 3000,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "sumEt", unitsy = "GeV",drawStyle = "COLZ")
h4 = ZJets.drawTH2('h_zjets_uperp_zpt_set_el', 'llnunu_l2_sumEt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+elChannel, str(lumi*1000), 1500, -15, 15, 300, 0, 3000,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "sumEt", unitsy = "GeV",drawStyle = "COLZ")

h1.Write()
h2.Write()
h3.Write()
h4.Write()

# muon
#h1 = Data.drawTH2('h_data_upara_zpt_mu', 'llnunu_l1_pt:(llnunu_l1_pt+llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+muChannel, "(1)", 1500, -15, 15, 1000, 0, 1000,  titlex = "#Delta P_{T}^{#parallel}(Z,MET)/P_{T}(Z)",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")
#h2 = ZJets.drawTH2('h_zjets_upara_zpt_mu', 'llnunu_l1_pt:(llnunu_l1_pt+llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+muChannel, str(lumi*1000), 1500, -15, 15, 1000, 0, 1000,  titlex = "#Delta P_{T}^{#parallel}(Z,MET)/P_{T}(Z)",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")


#h3 = Data.drawTH2('h_data_uperp_zpt_mu', 'llnunu_l1_pt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+muChannel, "(1)", 1500, -15, 15, 1000, 0, 1000,  titlex = "#Delta P_{T}^{#perp}(Z,MET)/P_{T}",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")
#h4 = ZJets.drawTH2('h_zjets_uperp_zpt_mu', 'llnunu_l1_pt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+muChannel, str(lumi*1000), 1500, -15, 15, 1000, 0, 1000,  titlex = "#Delta P_{T}^{#perp}(Z,MET)/P_{T}",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")

h1 = Data.drawTH2('h_data_upara_zpt_mu', 'llnunu_l1_pt:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+muChannel, "(1)", 1500, -15, 15, 1000, 0, 1000,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")
h2 = ZJets.drawTH2('h_zjets_upara_zpt_mu', 'llnunu_l1_pt:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+muChannel, str(lumi*1000), 1500, -15, 15, 1000, 0, 1000,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")


h3 = Data.drawTH2('h_data_uperp_zpt_mu', 'llnunu_l1_pt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+muChannel, "(1)", 1500, -15, 15, 1000, 0, 1000,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")
h4 = ZJets.drawTH2('h_zjets_uperp_zpt_mu', 'llnunu_l1_pt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+muChannel, str(lumi*1000), 1500, -15, 15, 1000, 0, 1000,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "P_{T}(Z)", unitsy = "",drawStyle = "COLZ")

h1.Write()
h2.Write()
h3.Write()
h4.Write()

h1 = Data.drawTH2('h_data_upara_zpt_zeta_mu', 'abs(llnunu_l1_eta):(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+muChannel, "(1)", 1500, -15, 15, 1000, 0, 10,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "#eta(Z)", unitsy = "",drawStyle = "COLZ")
h2 = ZJets.drawTH2('h_zjets_upara_zpt_zeta_mu', 'abs(llnunu_l1_eta):(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+muChannel, str(lumi*1000), 1500, -15, 15, 1000, 0, 10,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "#eta(Z)", unitsy = "",drawStyle = "COLZ")


h3 = Data.drawTH2('h_data_uperp_zpt_zeta_mu', 'abs(llnunu_l1_eta):(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+muChannel, "(1)", 1500, -15, 15, 1000, 0, 10,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "#eta(Z)", unitsy = "",drawStyle = "COLZ")
h4 = ZJets.drawTH2('h_zjets_uperp_zpt_zeta_mu', 'abs(llnunu_l1_eta):(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+muChannel, str(lumi*1000), 1500, -15, 15, 1000, 0, 10,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "#eta(Z)", unitsy = "",drawStyle = "COLZ")

h1.Write()
h2.Write()
h3.Write()
h4.Write()


h1 = Data.drawTH2('h_data_upara_zpt_nvtx_mu', 'nVert:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+muChannel, "(1)", 1500, -15, 15, 50, 0, 50,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "nVertex", unitsy = "",drawStyle = "COLZ")
h2 = ZJets.drawTH2('h_zjets_upara_zpt_nvtx_mu', 'nVert:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+muChannel, str(lumi*1000), 1500, -15, 15, 50, 0, 50,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "N vertex", unitsy = "",drawStyle = "COLZ")


h3 = Data.drawTH2('h_data_uperp_zpt_nvtx_mu', 'nVert:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+muChannel, "(1)", 1500, -15, 15, 50, 0, 50,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "N vertex", unitsy = "",drawStyle = "COLZ")
h4 = ZJets.drawTH2('h_zjets_uperp_zpt_nvtx_mu', 'nVert:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+muChannel, str(lumi*1000), 1500, -15, 15, 50, 0, 50,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "N vertex", unitsy = "",drawStyle = "COLZ")

h1.Write()
h2.Write()
h3.Write()
h4.Write()

h1 = Data.drawTH2('h_data_upara_zpt_set_mu', 'llnunu_l2_sumEt:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+muChannel, "(1)", 1500, -15, 15, 300, 0, 3000,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "sumEt", unitsy = "GeV",drawStyle = "COLZ")
h2 = ZJets.drawTH2('h_zjets_upara_zpt_set_mu', 'llnunu_l2_sumEt:(llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts+'&&'+muChannel, str(lumi*1000), 1500, -15, 15, 300, 0, 3000,  titlex = "MET^{#parallel}/P_{T}(Z)",unitsx = "", titley = "sumEt", unitsy = "GeV",drawStyle = "COLZ")


h3 = Data.drawTH2('h_data_uperp_zpt_set_mu', 'llnunu_l2_sumEt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+muChannel, "(1)", 1500, -15, 15, 300, 0, 3000,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "sumEt", unitsy = "GeV",drawStyle = "COLZ")
h4 = ZJets.drawTH2('h_zjets_uperp_zpt_set_mu', 'llnunu_l2_sumEt:(llnunu_l2_pt*sin(llnunu_l1_phi-llnunu_l2_phi))/llnunu_l1_pt', cuts+'&&'+muChannel, str(lumi*1000), 1500, -15, 15, 300, 0, 3000,  titlex = "MET^{#perp}/P_{T}",unitsx = "", titley = "sumEt", unitsy = "GeV",drawStyle = "COLZ")

h1.Write()
h2.Write()
h3.Write()
h4.Write()

'''
for (plotter,typeP,label,name) in zip(Stack.plotters,Stack.types,Stack.labels,Stack.names):

    lumi_str = str(lumi*1000)
    if typeP=='data': lumi_str = '(1)'
        
    h1 = plotter.drawTH1('h_'+name+'_el', 'llnunu_mta',cuts+'&&'+elChannel,lumi_str,500,0,5000,titlex = "M_{T}",units= "GeV", drawStyle = "HIST")
    h1.Write()

    h2 = plotter.drawTH1('h_'+name+'_mu', 'llnunu_mta',cuts+'&&'+muChannel,lumi_str,500,0,5000,titlex = "M_{T}",units= "GeV", drawStyle = "HIST")
    h2.Write()
'''    
'''
# electron
h_Data_el = Data.drawTH1('h_data_obs_el','llnunu_mta',cuts+'&&'+elChannel,'(1)',500,0,5000,titlex = "M_{T}",units= "GeV", drawStyle = "HIST")
h_Data_el.Write()

h_VVNonReso_el = WW.drawTH1('h_VVNonReso_el','llnunu_mta',cuts+'&&'+elChannel,str(lumi*1000),250,0,2500,300,0,3000,500,0,5000,titlex = "M_{T}",unitsx = "GeV",drawStyle = "HIST")
h_VVNonReso_el.Write()

h_VVZReso_el = VV.drawTH3('h3d_VVZReso_el','llnunu_mta:llnunu_l2_pt:llnunu_l1_pt',cuts+'&&'+elChannel,str(lumi*1000),250,0,2500,300,0,3000,500,0,5000,titlex = "P_{T}(Z)",unitsx = "GeV",titley = "MET",unitsy = "GeV",drawStyle = "COLZ")
h_VVZReso_el.Write()

h_ZJets_el = ZJets.drawTH3('h3d_ZJets_el','llnunu_mta:llnunu_l2_pt:llnunu_l1_pt',cuts+'&&'+elChannel,str(lumi*1000),250,0,2500,300,0,3000,500,0,5000,titlex = "P_{T}(Z)",unitsx = "GeV",titley = "MET",unitsy = "GeV",drawStyle = "COLZ")
h_ZJets_el.Write()

h_TT_el = TT.drawTH3('h3d_TT_el','llnunu_mta:llnunu_l2_pt:llnunu_l1_pt',cuts+'&&'+elChannel,str(lumi*1000),250,0,2500,300,0,3000,500,0,5000,titlex = "P_{T}(Z)",unitsx = "GeV",titley = "MET",unitsy = "GeV",drawStyle = "COLZ")
h_TT_el.Write()

# muon
h_Data_mu = Data.drawTH3('h3d_data_obs_mu','llnunu_mta:llnunu_l2_pt:llnunu_l1_pt',cuts+'&&'+muChannel,'(1)',250,0,2500,300,0,3000,500,0,5000,titlex = "P_{T}(Z)",unitsx = "GeV",titley = "MET",unitsy = "GeV",drawStyle = "COLZ")
h_Data_mu.Write()

h_VVNonReso_mu = WW.drawTH3('h3d_VVNonReso_mu','llnunu_mta:llnunu_l2_pt:llnunu_l1_pt',cuts+'&&'+muChannel,str(lumi*1000),250,0,2500,300,0,3000,500,0,5000,titlex = "P_{T}(Z)",unitsx = "GeV",titley = "MET",unitsy = "GeV",drawStyle = "COLZ")
h_VVNonReso_mu.Write()

h_VVZReso_mu = VV.drawTH3('h3d_VVZReso_mu','llnunu_mta:llnunu_l2_pt:llnunu_l1_pt',cuts+'&&'+muChannel,str(lumi*1000),250,0,2500,300,0,3000,500,0,5000,titlex = "P_{T}(Z)",unitsx = "GeV",titley = "MET",unitsy = "GeV",drawStyle = "COLZ")
h_VVZReso_mu.Write()

h_ZJets_mu = ZJets.drawTH3('h3d_ZJets_mu','llnunu_mta:llnunu_l2_pt:llnunu_l1_pt',cuts+'&&'+muChannel,str(lumi*1000),250,0,2500,300,0,3000,500,0,5000,titlex = "P_{T}(Z)",unitsx = "GeV",titley = "MET",unitsy = "GeV",drawStyle = "COLZ")
h_ZJets_mu.Write()

h_TT_mu = TT.drawTH3('h3d_TT_mu','llnunu_mta:llnunu_l2_pt:llnunu_l1_pt',cuts+'&&'+muChannel,str(lumi*1000),250,0,2500,300,0,3000,500,0,5000,titlex = "P_{T}(Z)",unitsx = "GeV",titley = "MET",unitsy = "GeV",drawStyle = "COLZ")
h_TT_mu.Write()

for i in range(len(sigSamples)):

    h_Sig_el = sigPlotters[i].drawTH3('h3d_'+sigSamples[i]+'_el','llnunu_mta:llnunu_l2_pt:llnunu_l1_pt',cuts+'&&'+elChannel,str(lumi*1000),250,0,2500,300,0,3000,500,0,5000,titlex = "P_{T}(Z)",unitsx = "GeV",titley = "MET",unitsy = "GeV",drawStyle = "COLZ")
    h_Sig_el.Write()

    h_Sig_mu = sigPlotters[i].drawTH3('h3d_'+sigSamples[i]+'_mu','llnunu_mta:llnunu_l2_pt:llnunu_l1_pt',cuts+'&&'+muChannel,str(lumi*1000),250,0,2500,300,0,3000,500,0,5000,titlex = "P_{T}(Z)",unitsx = "GeV",titley = "MET",unitsy = "GeV",drawStyle = "COLZ")
    h_Sig_mu.Write()
'''
fout.Close()



############
'''
# QCD uncertainty

# samples
#samples = wwSamples+vvSamples+zjetsSamples+ttSamples+sigSamples
 
if stag=='bkg': samples = wwSamples+vvSamples+ttSamples
elif stag=='zjets': samples = zjetsSamples
elif stag=='sig': samples = sigSamples
else: samples = wwSamples+vvSamples+zjetsSamples+ttSamples+sigSamples
 
#samples=['ZZTo2L2Nu']
 
# LHE weights input dir
lhedir='LHEWeights'


# qcd weights id
wtIds_qcd = [1,2,3,4,6,8] 

fout = ROOT.TFile("uncert_acc_qcd_"+stag+"_"+channel+".root", "recreate")



hunc_qcd = {}
all_var_qcd = {}
all_unc_qcd = {}
#for sample in allPlotters.keys():
#for sample in [sample]: 
for sample in samples: 

    print 'QCD for '+sample

    # open pickle file
    lhepkl= open(lhedir+'/'+sample+'/LHEWeightAnalyzer/LHEWeightReport.pck')

    data = pickle.load(lhepkl)

    hunc_qcd[sample] = ROOT.TH1F(sample+'_hunc_qcd', sample+'_hunc_qcd', 200,0,2)
 
    h = allPlotters[sample].drawTH1(sample+'_mh_qcd_wt_0','llnunu_mta',cuts,str(lumi*1000),50,0,5000,titlex = "M_{T}",units = "GeV",drawStyle = "HIST")
    n0 = h.Integral()
    print ' n0(pass) = '+str(n0)

    # if no event, no uncertainty
    if n0<=0.0 : 
        print ' Sample '+sample+' has Zero entry in singal region'
        continue 
    all_var_qcd[sample] = []

    for idx in wtIds_qcd:
        print ' - QCD '+str(idx)

        # lhe idx
        lheidx = idx+1001 # for background
        lheId0 = 1001
        if 'BulkGrav' in sample:
            lheidx = idx+1 # for signal
            lheId0 = 1
        # modify weight 
        allPlotters[sample].changeCorrectionFactor('(genWeight*LHEweight_wgt['+str(idx)+']/LHEweight_wgt[0])','genWeight')
        # get n = ni(pass)
        h = allPlotters[sample].drawTH1(sample+'_mh_qcd_wt_'+str(idx),'llnunu_mta',cuts,str(lumi*1000),50,0,5000,titlex = "M_{T}",units = "GeV",drawStyle = "HIST") 
        n = h.Integral()

        # frac = ni(pass) / n0(pass)
        frac = n/n0
        # get tfrac = ni(total) / n0(total)
        tfrac = float(data[str(lheidx)][1])/float(data[str(lheId0)][1])

        # unc = [ni(pass)/ni(total)] / [n0(pass)/n0(total)] 
        #     = [ni(pass)/n0(pass)] / [ni(total)/n0(total)]
        #     = frac / tfrac
        unc = frac/tfrac
        print '  ni(pass) = {nip}, ni(pass)/n0(pass) = {frac}, ni(total)/n0(total) = {tfrac}, unc = {unc}'.format(nip=n,frac=frac,tfrac=tfrac,unc=unc)


        hunc_qcd[sample].Fill(unc)
        all_var_qcd[sample].append(unc)

    qcd_up = max(all_var_qcd[sample])
    qcd_dn = min(all_var_qcd[sample])
    qcd_unc = abs(qcd_up-qcd_dn)/2.0
    print '  qcd_unc = '+str(qcd_unc)
    hunc_qcd[sample].Write()

    all_unc_qcd[sample] = qcd_unc


print 'unc_qcd =',all_unc_qcd
print 'var_qcd =',all_var_qcd


fout.Close()

# print table
print '\\begin{table}[htdp]'
print '\\caption{QCD uncertainties on acceptance for background.}'
print '\\begin{center}'
print '\\begin{footnotesize}'
print '\\begin{tabular}{c c  }'
print '\\hline'
print ' MC sample & unc.   \\\\'
print '\\hline'
for bkg in samples:
    if not (bkg in all_unc_qcd.keys()):
        all_unc_qcd[bkg] = 0.0
    print '{bkg} & {bpdf:.3f} \\\\'.format(bkg=bkg,bpdf=all_unc_qcd[bkg])
print '\\hline'
print '\\end{tabular}'
print '\\end{footnotesize}'
print '\\end{center}'
print '\\label{default}'
print '\\end{table}'
'''



#################



#Stack.drawStack('llnunu_mta', cuts, str(lumi*1000), 60, 0.0, 1200.0, titlex = "M_{T}", units = "GeV",output=tag+'mt',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=300)
#Stack.drawStack('llnunu_mt', cuts, str(lumi*1000), 150, 0.0, 3000.0, titlex = "M_{T}", units = "GeV",output=tag+'mt_high',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=300)
#Stack.drawStack('llnunu_l2_pt', cuts, str(lumi*1000), 100, 0, 1000, titlex = "MET", units = "GeV",output=tag+'met',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=200)
#Stack.drawStack('llnunu_l2_pt', cuts, str(lumi*1000), 50, 0, 500, titlex = "MET", units = "GeV",output=tag+'met_low',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=200)
#Stack.drawStack('llnunu_l2_pt', cuts, str(lumi*1000), 150, 0, 1500, titlex = "MET", units = "GeV",output=tag+'met_high',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=200)
#Stack.drawStack('met_JECJER_pt', cuts, str(lumi*1000), 50, 0, 500, titlex = "MET", units = "GeV",output=tag+'met_jecjer_low',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=200)
#Stack.drawStack('nVert', cuts, str(lumi*1000), 50, 0.0, 50.0, titlex = "N vertices", units = "",output=tag+'nVert',outDir=outdir,separateSignal=sepSig)
#Stack.drawStack('llnunu_l2_sumEt', cuts, str(lumi*1000), 100, 0.0, 3000.0, titlex = "sumE_{T}", units = "GeV",output=tag+'metSumEt',outDir=outdir,separateSignal=sepSig)
#Stack.drawStack('llnunu_l1_pt', cuts, str(lumi*1000), 75, 0.0, 1500.0, titlex = "P_{T}(Z)", units = "GeV",output=tag+'zpt',outDir=outdir,separateSignal=sepSig)


if not test :
    Stack.drawStack('llnunu_mta', cuts, str(lumi*1000), 60, 0.0, 1200.0, titlex = "M_{T}", units = "GeV",output=tag+'mt',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=300)
    Stack.drawStack('llnunu_l2_pt', cuts, str(lumi*1000), 50, 0, 500, titlex = "MET", units = "GeV",output=tag+'met_low',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=200)
    Stack.drawStack('llnunu_l1_pt', cuts, str(lumi*1000), 75, 0.0, 1500.0, titlex = "P_{T}(Z)", units = "GeV",output=tag+'zpt',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_mass', cuts, str(lumi*1000), 50, 50, 150, titlex = "M(Z)", units = "GeV",output=tag+'zmass',outDir=outdir,separateSignal=sepSig)
    #Stack.drawStack('llnunu_l1_pt', cuts, str(lumi*1000), 75, 0.0, 1500.0, titlex = "P_{T}(Z)", units = "GeV",output=tag+'zpt',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l2_pt/sqrt(llnunu_l2_sumEt)', cuts, str(lumi*1000), 100, 0.0, 20.0, titlex = "MET/#sqrt{sumE_{T}}", units = "",output=tag+'metOvSqSET',outDir=outdir,separateSignal=sepSig)

    Stack.drawStack('llnunu_deltaPhi', cuts, str(lumi*1000), 50, 0, 5, titlex = "#Delta#phi(Z,MET)", units = "",output=tag+'dphiZMet',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('cos(llnunu_deltaPhi)', cuts, str(lumi*1000), 50, -1, 1, titlex = "cos(#Delta#phi(Z,MET))", units = "",output=tag+'CosdphiZMet',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('abs(llnunu_l1_pt+llnunu_l2_pt*cos(llnunu_deltaPhi))', cuts, str(lumi*1000), 100, 0, 1000, titlex = "#Delta P_{T}^{#parallel}(Z,MET)", units = "GeV",output=tag+'dPTPara',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('abs(llnunu_l2_pt*sin(llnunu_deltaPhi))', cuts, str(lumi*1000), 100, 0, 100, titlex = "#Delta P_{T}^{#perp}(Z,MET)", units = "GeV",output=tag+'dPTPerp',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('abs(llnunu_l1_pt+llnunu_l2_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt', cuts, str(lumi*1000), 100, 0, 5, titlex = "#Delta P_{T}^{#parallel}(Z,MET)/P_{T}(Z)", units = "",output=tag+'dPTParaRel',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('abs(llnunu_l2_pt*sin(llnunu_deltaPhi))/llnunu_l1_pt', cuts, str(lumi*1000), 100, 0, 5, titlex = "#Delta P_{T}^{#perp}(Z,MET)/P_{T}", units = "",output=tag+'dPTPerpRel',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('nVert', cuts, str(lumi*1000), 50, 0.0, 50.0, titlex = "N vertices", units = "",output=tag+'nVert',outDir=outdir,separateSignal=sepSig)
 
    Stack.drawStack('llnunu_l1_eta', cuts, str(lumi*1000), 100, -5.0, 5.0, titlex = "#eta(Z) ", units = "",output=tag+'zeta',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_phi', cuts, str(lumi*1000), 64, -3.2, 3.2, titlex = "#phi(Z)", units = "",output=tag+'zphi',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_mt', cuts, str(lumi*1000), 50, 0.0, 150.0, titlex = "M_{T}(Z)", units = "GeV",output=tag+'zmt',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_deltaPhi', cuts, str(lumi*1000), 50, 0.0, 5.0, titlex = "#Delta#phi", units = "",output=tag+'ZdeltaPhi',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('cos(llnunu_l1_deltaPhi)', cuts, str(lumi*1000), 50, -1, 1, titlex = "Cos(#Delta#phi)", units = "",output=tag+'CosZdeltaPhi',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_deltaR', cuts, str(lumi*1000), 50, 0.0, 5.0, titlex = "#Delta R", units = "",output=tag+'ZdeltaR',outDir=outdir,separateSignal=sepSig)
  
    Stack.drawStack('llnunu_l2_phi', cuts, str(lumi*1000), 100, -3.2, 3.2, titlex = "#phi(MET)", units = "",output=tag+'metPhi',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l2_sumEt', cuts, str(lumi*1000), 80, 0.0, 3000.0, titlex = "sumE_{T}", units = "GeV",output=tag+'metSumEt',outDir=outdir,separateSignal=sepSig)

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
#os.system('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+outdir+'/'+tag+'all.pdf '+outdir+'/'+tag+'*.eps')
#print 'All plots merged in single pdf file '+tag+'.pdf .'
# merge root file
#os.system('rm '+outdir+'/'+tag+'all.root ')
#os.system('hadd -f '+outdir+'/'+tag+'all.root '+outdir+'/'+tag+'*.root')



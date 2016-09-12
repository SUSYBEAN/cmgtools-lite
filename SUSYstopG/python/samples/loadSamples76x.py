#####################
# load samples 
#####################

import PhysicsTools.HeppyCore.framework.config as cfg
import os

# Load backgrounds
from CMGTools.XZZ2l2nu.samples.samples_13TeV_RunIIFall15MiniAODv2 import *
# Load signals
from CMGTools.XZZ2l2nu.samples.samples_13TeV_signal76X import *
# Load Data 
from CMGTools.XZZ2l2nu.samples.samples_13TeV_DATA2015 import *
# Load triggers
from CMGTools.XZZ2l2nu.samples.triggers_13TeV_Spring15 import *
# Load reference sample for jet energy corrections and jet resolution
from CMGTools.XZZ2l2nu.samples.samples_JERCReference import * 

# mc reference sample for jet energy corrections and jet resolution
jercRefMC=[JERCRef_MC, JERCRef_MC_eos]
jercRefdt=[JERCRef_data]

# backgrounds
backgroundSamples = [
DYJetsToLL_M50,
DYJetsToLL_M50_Ext,
DYJetsToLL_M50_PtZ100,
WJetsToLNu,
TTTo2L2Nu,
ZZTo2L2Nu,
ZZTo4L,
ZZTo2L2Q,
WWTo2L2Nu,
WWToLNuQQ,
WZTo1L1Nu2Q,
WZTo2L2Q,
WZTo3LNu,
TTZToLLNuNu,
TTWJetsToLNu,
ggZZTo2e2nu,
ggZZTo2mu2nu,
]

extraBackgroundMC = [
# DYJetsToLL_M50_HT100to200,
# DYJetsToLL_M50_HT200to400,
# DYJetsToLL_M50_HT400to600,
# DYJetsToLL_M50_HT600toInf,
WW,
WZ,
ZZ,
TT,
]

# signals
signalSamples = [
BulkGravToZZToZlepZinv_narrow_600,
BulkGravToZZToZlepZinv_narrow_800,
BulkGravToZZToZlepZinv_narrow_1000,
BulkGravToZZToZlepZinv_narrow_1200,
BulkGravToZZToZlepZinv_narrow_1400,
BulkGravToZZToZlepZinv_narrow_1600,
BulkGravToZZToZlepZinv_narrow_1800,
BulkGravToZZToZlepZinv_narrow_2000,
BulkGravToZZToZlepZinv_narrow_2500,
BulkGravToZZToZlepZinv_narrow_3000,
BulkGravToZZToZlepZinv_narrow_3500,
BulkGravToZZToZlepZinv_narrow_4000,
BulkGravToZZToZlepZinv_narrow_4500,
]

BulkGravToZZToZlepZhad = [
BulkGravToZZToZlepZhad_narrow_600,
BulkGravToZZToZlepZhad_narrow_800,
BulkGravToZZToZlepZhad_narrow_1000,
BulkGravToZZToZlepZhad_narrow_1200,
BulkGravToZZToZlepZhad_narrow_1400,
BulkGravToZZToZlepZhad_narrow_1600,
BulkGravToZZToZlepZhad_narrow_1800,
BulkGravToZZToZlepZhad_narrow_2000,
BulkGravToZZToZlepZhad_narrow_2500,
BulkGravToZZToZlepZhad_narrow_3000,
BulkGravToZZToZlepZhad_narrow_3500,
BulkGravToZZToZlepZhad_narrow_4000,
BulkGravToZZToZlepZhad_narrow_4500,
]

# MC samples
otherMcSamples = BulkGravToZZToZlepZhad  
mcSamples = signalSamples + backgroundSamples + jercRefMC + extraBackgroundMC
#mcSamples = backgroundSamples

# data
SingleMuon=[SingleMuon_Run2015C_25ns_16Dec,SingleMuon_Run2015D_16Dec]
SingleElectron=[SingleElectron_Run2015C_25ns_16Dec,SingleElectron_Run2015D_16Dec]
SinglePhoton=[SinglePhoton_Run2015C_25ns_16Dec,SinglePhoton_Run2015D_16Dec]
MuEG=[MuonEG_Run2015C_25ns_16Dec,MuonEG_Run2015D_16Dec]

for s in SingleMuon:
    s.triggers = triggers_1mu_noniso
    s.vetoTriggers = []
for s in SingleElectron:
    s.triggers = triggers_1e_noniso
    s.vetoTriggers = triggers_1mu_noniso
for s in SinglePhoton:
    s.trigers = triggers_photon_unbias


dataSamples=SingleMuon+SingleElectron #+jercRefdt

# JSON
# https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2015Analysis#ReReco_at_25_ns
# note, but the golden one doesn't change wrt. the prompt reco one as expected. 
silverJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Silver.txt'
goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt'

jsonFile = goldenJson


from CMGTools.XZZ2l2nu.setup.Efficiencies import *
dataDir = "$CMSSW_BASE/src/CMGTools/XZZ2l2nu/data"

for comp in mcSamples+otherMcSamples:
    comp.isMC = True
    comp.isData = False
    comp.splitFactor = 250
    comp.puFileMC=dataDir+"/pileup_MC_76x.root"
    comp.puFileData=dataDir+"/pileup_DATA.root"
    comp.eSFinput=dataDir+"/CutBasedID_LooseWP_76X_18Feb.txt_SF2D.root"
    comp.efficiency = eff2012
    comp.triggers=triggers_1mu_noniso+triggers_1e_noniso
    comp.globalTag = "Summer15_25nsV6_MC"

for comp in dataSamples:
    comp.splitFactor = 250
    comp.isMC = False
    comp.isData = True
    comp.json = jsonFile
    comp.globalTag = "Summer15_25nsV6_DATA"


#####################
# load samples 
#####################

import PhysicsTools.HeppyCore.framework.config as cfg
import os

# Load backgrounds
from CMGTools.XZZ2l2nu.samples.samples_13TeV_RunIISpring15MiniAODv2 import *
# Load signals
from CMGTools.XZZ2l2nu.samples.samples_13TeV_signal import *
# Load Data 
from CMGTools.XZZ2l2nu.samples.samples_13TeV_DATA2015 import *
# Load triggers
from CMGTools.XZZ2l2nu.samples.triggers_13TeV_Spring15 import *


# backgrounds
backgroundSamples = [
DYJetsToLL_M50_HT100to200,
DYJetsToLL_M50_HT200to400,
DYJetsToLL_M50_HT400to600,
DYJetsToLL_M50_HT600toInf,
WW,
WZ,
ZZ,
DYJetsToLL_M50,
WJetsToLNu,
TT_pow,
TTLep_pow,
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
]


# signals
signalSamples = [
#RSGravToZZToZZinv_narrow_600,
RSGravToZZToZZinv_narrow_800,
RSGravToZZToZZinv_narrow_1000,
RSGravToZZToZZinv_narrow_1200,
RSGravToZZToZZinv_narrow_1400,
RSGravToZZToZZinv_narrow_1600,
RSGravToZZToZZinv_narrow_2000,
RSGravToZZToZZinv_narrow_2500,
RSGravToZZToZZinv_narrow_3000,
RSGravToZZToZZinv_narrow_3500,
RSGravToZZToZZinv_narrow_4500,
#BulkGravToZZ_narrow_600,
#BulkGravToZZ_narrow_800,
#BulkGravToZZ_narrow_1000,
#BulkGravToZZ_narrow_1200,
#BulkGravToZZ_narrow_1400,
#BulkGravToZZ_narrow_1600,
#BulkGravToZZ_narrow_1800,
#BulkGravToZZ_narrow_2000,
#BulkGravToZZ_narrow_2500,
#BulkGravToZZ_narrow_3000,
#BulkGravToZZ_narrow_3500,
#BulkGravToZZ_narrow_4000,
#BulkGravToZZ_narrow_4500,
]

# mc samples
mcSamples = signalSamples + backgroundSamples

# other mc samples
bulkJetsSamples = [
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

otherMcSamples = bulkJetsSamples

# data
SingleMuon=[SingleMuon_Run2015D_Promptv4,SingleMuon_Run2015D_05Oct]
SingleElectron=[SingleElectron_Run2015D_Promptv4,SingleElectron_Run2015D_05Oct]
SinglePhoton=[SinglePhoton_Run2015D_Promptv4,SinglePhoton_Run2015D_05Oct]# added by Mengqing

for s in SingleMuon:
    s.triggers = triggers_1mu_noniso
    s.vetoTriggers = []
for s in SingleElectron:
    s.triggers = triggers_1e_noniso
    s.vetoTriggers = triggers_1mu_noniso

dataSamples=SingleMuon+SingleElectron

# JSON
silverJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt'
goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'

jsonFile = goldenJson


from CMGTools.XZZ2l2nu.setup.Efficiencies import *
dataDir = "$CMSSW_BASE/src/CMGTools/XZZ2l2nu/data"

for comp in mcSamples+otherMcSamples:
    comp.isMC = True
    comp.isData = False
    comp.splitFactor = 250
    comp.puFileMC=dataDir+"/pileup_MC.root"
    comp.puFileData=dataDir+"/pileup_DATA.root"
    comp.efficiency = eff2012
    comp.triggers=triggers_1mu_noniso+triggers_1e_noniso
    comp.globalTag = "Summer15_25nsV6_MC"

for comp in dataSamples:
    comp.splitFactor = 250
    comp.isMC = False
    comp.isData = True
    comp.json = jsonFile
    comp.globalTag = "Summer15_25nsV6_DATA"


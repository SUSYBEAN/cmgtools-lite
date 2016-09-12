##########################################################
##      configuration for XZZ2l2nu 
##########################################################

import CMGTools.XZZ2l2nu.fwlite.Config as cfg
from CMGTools.XZZ2l2nu.fwlite.Config import printComps
from CMGTools.XZZ2l2nu.RootTools import *
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption


#Load all common analyzers
from CMGTools.XZZ2l2nu.analyzers.coreXZZ_cff import *

#-------- SAMPLES AND TRIGGERS -----------
from CMGTools.XZZ2l2nu.samples.loadSamples76x import *
selectedComponents = mcSamples+dataSamples

triggerFlagsAna.triggerBits ={
    "ISOMU":triggers_1mu_iso,
    "MU":triggers_1mu_noniso,
    "ISOELE":triggers_1e,
    "ELE":triggers_1e_noniso,
    "HT800":triggers_HT800,
    "HT900":triggers_HT900,
    "JJ":triggers_dijet_fat,
    "MET90":triggers_met90_mht90+triggers_metNoMu90_mhtNoMu90,
    "MET120":triggers_metNoMu120_mhtNoMu120
}

#-------- Analyzer
from CMGTools.XZZ2l2nu.analyzers.treeXZZ_cff import *

# make VV selection loosened as ElMu selection
leptonicVAna.selectVBoson = (lambda x: x.mass()>30.0 and x.mass()<180.0)
leptonicVAna.selectFakeBoson = (lambda x: x.mass()>30.0 and x.mass()<180.0)
multiStateAna.selectPairLLNuNu = (lambda x: x.leg1.mass()>30.0 and x.leg1.mass()<180.0)
multiStateAna.selectPairElMuNuNu = (lambda x: x.leg1.mass()>30.0 and x.leg1.mass()<180.0)

vvSkimmer.required = ['LLNuNu', 'ElMuNuNu']
leptonicVAna.doElMu = True

vvTreeProducer.globalVariables = [
         NTupleVariable("nLL",lambda ev: len(ev.LL) , int),
         NTupleVariable("nElMu",lambda ev: len(ev.ElMu) , int),       
         NTupleVariable("nLLNuNu",lambda ev: len(ev.LLNuNu) , int),
         NTupleVariable("nVert",  lambda ev: len(ev.goodVertices), int, help="Number of good vertices"),
         NTupleVariable("triggersf",  lambda x : getattr(x,'trgsf',1), help="singleelectron/muon trigger sf"),
         NTupleVariable("triggersfUp",  lambda x : getattr(x,'trgsfUp',1), help="singleelectron/muon trigger sf upper"),
         NTupleVariable("triggersfLo",  lambda x : getattr(x,'trgsfLo',1), help="singleelectron/muon trigger sf lower"),
     ]

vvTreeProducer.collections = {
         "genLeptons" : NTupleCollection("genLep", genParticleType, 10, help="Generated leptons (e/mu) from W/Z decays"),
         "LLNuNu"     : NTupleCollection("llnunu",LLNuNuType ,5, help="VV candidate with di-lepton and MET"),
         "ElMuNuNu"     : NTupleCollection("elmununu",LLNuNuType ,5, help="Fake VV candidate with el-mu and MET"),
     }
jetAna.smearJets=False

#-------- SEQUENCE
#sequence = cfg.Sequence(coreSequence+[vvSkimmer,vvTreeProducer])

coreSequence = [
    skimAnalyzer,
    genAna,
    #lheWeightAna,
    jsonAna,
    triggerAna,
    pileUpAna,
    vertexAna,
    lepAna,
    jetAna,
    metAna,
    leptonicVAna,
    multiStateAna,
    eventFlagsAna,
#    triggerFlagsAna,
]
    
#sequence = cfg.Sequence(coreSequence+[vvTreeProducer])
#sequence = cfg.Sequence(coreSequence+[vvSkimmer,vvTreeProducer])
sequence = cfg.Sequence(coreSequence+[vvSkimmer,fullTreeProducer])
 

#-------- HOW TO RUN
test = 1
if test==1:
    # test a single component, using a single thread.
    #selectedComponents = dataSamples
    #selectedComponents = mcSamples
    #selectedComponents = [SingleMuon_Run2015D_Promptv4,SingleElectron_Run2015D_Promptv4]
    #selectedComponents = [SingleElectron_Run2015D_Promptv4,SingleElectron_Run2015D_05Oct]
    #selectedComponents = [RSGravToZZToZZinv_narrow_800]
    #selectedComponents = [WWTo2L2Nu] #DYJetsToLL_M50
    selectedComponents = backgroundSamples # signalSamples
    #selectedComponents = [BulkGravToZZToZlepZinv_narrow_1400, BulkGravToZZToZlepZinv_narrow_1600]
    #selectedComponents = [BulkGravToZZToZlepZinv_narrow_800]
    #selectedComponents = [WWToLNuQQ, WZTo1L1Nu2Q, WZTo2L2Q, ZZTo2L2Nu]
    #selectedComponents = [ZZTo2L2Nu,DYJetsToLL_M50]
    #selectedComponents = [BulkGravToZZ_narrow_800]
    #selectedComponents = [BulkGravToZZToZlepZhad_narrow_800]
    for c in selectedComponents:
        #c.files = c.files[0]
        c.splitFactor = (len(c.files)/20 if len(c.files)>20 else 1)
        #c.splitFactor = 1
        #c.triggers=triggers_1mu_noniso
        #c.triggers=triggers_1e_noniso

## output histogram
outputService=[]
from PhysicsTools.HeppyCore.framework.services.tfile import TFileService
output_service = cfg.Service(
    TFileService,
    'outputfile',
    name="outputfile",
    fname='vvTreeProducer/tree.root',
    option='recreate'
    )
outputService.append(output_service)

from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
from CMGTools.TTHAnalysis.tools.EOSEventsWithDownload import EOSEventsWithDownload
event_class = EOSEventsWithDownload
event_class = Events
if getHeppyOption("nofetch"):
    event_class = Events
config = cfg.Config( components = selectedComponents,
                     sequence = sequence,
                     services = [],
                     events_class = event_class)




